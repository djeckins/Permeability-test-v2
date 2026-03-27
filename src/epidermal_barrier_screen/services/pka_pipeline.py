"""Per-pKa pipeline: collect all pKa values, compute per-pKa ionization & logD.

This module sits on top of the existing ionization.py cascade and produces
an exploded per-pKa table alongside the aggregate molecule-level results.
"""
from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from typing import Any

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Per-pKa entry dataclass
# ---------------------------------------------------------------------------
@dataclass
class PkaEntry:
    """One pKa value with full provenance and per-pKa calculations."""
    pka_value: float
    source: str  # "pubchem", "drugbank_live", "chembl", "qupkake", "input"
    source_record_id: str | None = None
    evidence_type: str = "unknown"  # "experimental", "computed", "ML_predicted", "user_provided"
    site_type: str = "unknown"  # "acid", "base", "unknown"
    site_label: str | None = None
    raw_text: str | None = None
    confidence_score: str = "moderate"

    # Per-pKa calculations (filled by compute_per_pka)
    ph: float | None = None
    neutral_fraction: float | None = None
    ionized_fraction: float | None = None
    logp: float | None = None
    logd: float | None = None
    calculation_status: str = "pending"
    warning: str | None = None


def _evidence_type_for_source(source: str) -> str:
    """Map source to evidence type."""
    if source in ("pubchem", "drugbank_live"):
        return "experimental"
    if source == "chembl":
        return "computed"
    if source == "qupkake":
        return "ML_predicted"
    if source == "input":
        return "user_provided"
    return "unknown"


# ---------------------------------------------------------------------------
# Per-pKa ionization & logD calculations
# ---------------------------------------------------------------------------
def _neutral_fraction_acid(pka: float, ph: float) -> float:
    """Fraction unionized for an acid: HA <-> H+ + A-."""
    try:
        ratio = 10 ** (ph - pka)
        return 1.0 / (1.0 + ratio)
    except OverflowError:
        return 0.0


def _neutral_fraction_base(pka: float, ph: float) -> float:
    """Fraction unionized for a base: BH+ <-> B + H+."""
    try:
        ratio = 10 ** (pka - ph)
        return 1.0 / (1.0 + ratio)
    except OverflowError:
        return 0.0


def _logd_acid(logp: float, pka: float, ph: float) -> float:
    """logD for acidic pKa: logD = logP - log10(1 + 10^(pH-pKa))."""
    try:
        return logp - math.log10(1.0 + 10 ** (ph - pka))
    except (OverflowError, ValueError):
        return logp


def _logd_base(logp: float, pka: float, ph: float) -> float:
    """logD for basic pKa: logD = logP - log10(1 + 10^(pKa-pH))."""
    try:
        return logp - math.log10(1.0 + 10 ** (pka - ph))
    except (OverflowError, ValueError):
        return logp


def compute_per_pka(entry: PkaEntry, ph: float, logp: float) -> PkaEntry:
    """Fill in per-pKa ionization and logD for a single entry."""
    entry.ph = ph
    entry.logp = logp

    if entry.site_type == "acid":
        entry.neutral_fraction = round(_neutral_fraction_acid(entry.pka_value, ph), 6)
        entry.ionized_fraction = round(1.0 - entry.neutral_fraction, 6)
        entry.logd = round(_logd_acid(logp, entry.pka_value, ph), 4)
        entry.calculation_status = "ok"
    elif entry.site_type == "base":
        entry.neutral_fraction = round(_neutral_fraction_base(entry.pka_value, ph), 6)
        entry.ionized_fraction = round(1.0 - entry.neutral_fraction, 6)
        entry.logd = round(_logd_base(logp, entry.pka_value, ph), 4)
        entry.calculation_status = "ok"
    else:
        entry.neutral_fraction = None
        entry.ionized_fraction = None
        entry.logd = None
        entry.calculation_status = "unknown_site_type"
        entry.warning = "site_type is unknown; cannot compute directional ionization"

    return entry


# ---------------------------------------------------------------------------
# Collect all pKa entries from an existing analyze_ionization result
# ---------------------------------------------------------------------------
def collect_pka_entries(
    ion_result: dict[str, Any],
    sites: list[Any] | None = None,
) -> list[PkaEntry]:
    """Extract all pKa values from an analyze_ionization result into PkaEntry list.

    This reads the existing acidic_pka_list, basic_pka_list, and database
    metadata fields already computed by ionization.py.
    """
    entries: list[PkaEntry] = []
    pka_source = ion_result.get("pka_source", "not_found")
    evidence = _evidence_type_for_source(pka_source)

    # Parse acidic pKa list
    acidic_str = ion_result.get("acidic_pka_list")
    if acidic_str:
        for val_str in str(acidic_str).split(";"):
            val_str = val_str.strip()
            if not val_str:
                continue
            try:
                v = float(val_str)
            except ValueError:
                continue
            entries.append(PkaEntry(
                pka_value=v,
                source=pka_source,
                source_record_id=_get_source_id(ion_result, pka_source),
                evidence_type=evidence,
                site_type="acid",
                confidence_score=ion_result.get("pka_confidence", "moderate"),
            ))

    # Parse basic pKa list
    basic_str = ion_result.get("basic_pka_list")
    if basic_str:
        for val_str in str(basic_str).split(";"):
            val_str = val_str.strip()
            if not val_str:
                continue
            try:
                v = float(val_str)
            except ValueError:
                continue
            entries.append(PkaEntry(
                pka_value=v,
                source=pka_source,
                source_record_id=_get_source_id(ion_result, pka_source),
                evidence_type=evidence,
                site_type="base",
                confidence_score=ion_result.get("pka_confidence", "moderate"),
            ))

    # Also collect any DB pKa values from other sources for provenance
    # (PubChem pKa values that weren't the primary source)
    pubchem_str = ion_result.get("pubchem_pka_values")
    if pubchem_str and pka_source != "pubchem":
        for val_str in str(pubchem_str).split(";"):
            val_str = val_str.strip()
            if not val_str:
                continue
            try:
                v = float(val_str)
            except ValueError:
                continue
            # Don't duplicate
            if not any(abs(e.pka_value - v) < 0.01 for e in entries):
                entries.append(PkaEntry(
                    pka_value=v,
                    source="pubchem",
                    source_record_id=str(ion_result.get("pubchem_cid")) if ion_result.get("pubchem_cid") else None,
                    evidence_type="experimental",
                    site_type=_guess_site_type_from_value(v, ion_result),
                    confidence_score="high",
                ))

    # DrugBank pKa values (if not primary source)
    if pka_source != "drugbank_live":
        for db_key, st in [("acidic_pka_drugbank", "acid"), ("basic_pka_drugbank", "base")]:
            val = ion_result.get(db_key)
            if val is not None:
                try:
                    v = float(val)
                except (ValueError, TypeError):
                    continue
                if not any(abs(e.pka_value - v) < 0.01 for e in entries):
                    entries.append(PkaEntry(
                        pka_value=v,
                        source="drugbank_live",
                        source_record_id=ion_result.get("drugbank_url"),
                        evidence_type="experimental",
                        site_type=st,
                        confidence_score="high" if ion_result.get("drugbank_match_status") == "exact" else "low",
                    ))

    # ChEMBL acd_pka (if not primary source)
    if pka_source != "chembl" and ion_result.get("chembl_acd_pka") is not None:
        v = float(ion_result["chembl_acd_pka"])
        if not any(abs(e.pka_value - v) < 0.01 for e in entries):
            entries.append(PkaEntry(
                pka_value=v,
                source="chembl",
                source_record_id=ion_result.get("chembl_id"),
                evidence_type="computed",
                site_type=_guess_site_type_from_value(v, ion_result),
                confidence_score="moderate",
            ))

    return entries


def _get_source_id(ion_result: dict[str, Any], source: str) -> str | None:
    if source == "pubchem":
        cid = ion_result.get("pubchem_cid")
        return str(cid) if cid else None
    if source == "drugbank_live":
        return ion_result.get("drugbank_url")
    if source == "chembl":
        return ion_result.get("chembl_id")
    return None


def _guess_site_type_from_value(pka_val: float, ion_result: dict[str, Any]) -> str:
    """Classify site_type using ionization_class from structural analysis."""
    ion_class = ion_result.get("ionization_class", "")
    if ion_class == "acid":
        return "acid"
    if ion_class == "base":
        return "base"
    if ion_class == "ampholyte":
        # For ampholytes, lower pKa values tend to be acidic
        return "acid" if pka_val < 7.0 else "base"
    return "unknown"


def build_pka_detail_rows(
    entries: list[PkaEntry],
    ph: float,
    logp: float,
    identity: dict[str, Any],
) -> list[dict[str, Any]]:
    """Build exploded per-pKa rows with identity columns + per-pKa calculations.

    Parameters
    ----------
    entries : list of PkaEntry (already collected)
    ph : user pH
    logp : logP value for logD calculation
    identity : dict with keys like name, canonical_smiles, inchikey, etc.

    Returns
    -------
    list of flat dicts, one per pKa entry
    """
    rows: list[dict[str, Any]] = []
    for entry in entries:
        compute_per_pka(entry, ph, logp)
        row = {
            "name": identity.get("name"),
            "canonical_smiles": identity.get("canonical_smiles"),
            "inchikey": identity.get("inchikey"),
            "pka_value": entry.pka_value,
            "pka_source": entry.source,
            "source_record_id": entry.source_record_id,
            "evidence_type": entry.evidence_type,
            "site_type": entry.site_type,
            "site_label": entry.site_label,
            "raw_text": entry.raw_text,
            "confidence_score": entry.confidence_score,
            "pH": entry.ph,
            "neutral_fraction": entry.neutral_fraction,
            "ionized_fraction": entry.ionized_fraction,
            "logP": entry.logp,
            "logD": entry.logd,
            "calculation_status": entry.calculation_status,
            "warning": entry.warning,
        }
        rows.append(row)
    return rows
