"""Per-pKa pipeline: preserve all values and compute row-level ionization/logD."""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from epidermal_barrier_screen.predictors.pka_predictor import PkaObservation


@dataclass
class PkaEntry:
    pka_value: float
    source: str
    source_record_id: str | None = None
    evidence_type: str = "unknown"
    site_type: str = "unknown"
    site_label: str | None = None
    raw_text: str | None = None
    confidence_score: str = "moderate"
    ph: float | None = None
    neutral_fraction: float | None = None
    ionized_fraction: float | None = None
    logp: float | None = None
    logd: float | None = None
    calculation_status: str = "pending"
    warning: str | None = None


def _neutral_fraction_acid(pka: float, ph: float) -> float:
    return 1.0 / (1.0 + 10 ** (ph - pka))


def _neutral_fraction_base(pka: float, ph: float) -> float:
    return 1.0 / (1.0 + 10 ** (pka - ph))


def _logd_acid(logp: float, pka: float, ph: float) -> float:
    return logp - math.log10(1.0 + 10 ** (ph - pka))


def _logd_base(logp: float, pka: float, ph: float) -> float:
    return logp - math.log10(1.0 + 10 ** (pka - ph))


def compute_per_pka(entry: PkaEntry, ph: float, logp: float) -> PkaEntry:
    entry.ph = ph
    entry.logp = logp
    if entry.site_type == "acid":
        entry.neutral_fraction = round(_neutral_fraction_acid(entry.pka_value, ph), 6)
        entry.ionized_fraction = round(1 - entry.neutral_fraction, 6)
        entry.logd = round(_logd_acid(logp, entry.pka_value, ph), 4)
        entry.calculation_status = "ok"
    elif entry.site_type == "base":
        entry.neutral_fraction = round(_neutral_fraction_base(entry.pka_value, ph), 6)
        entry.ionized_fraction = round(1 - entry.neutral_fraction, 6)
        entry.logd = round(_logd_base(logp, entry.pka_value, ph), 4)
        entry.calculation_status = "ok"
    else:
        entry.calculation_status = "unknown_site_type"
        entry.warning = "site_type unknown; logD left null"
    return entry


def collect_pka_entries(ion_result: dict[str, Any], sites: list[Any] | None = None) -> list[PkaEntry]:
    entries: list[PkaEntry] = []
    observations = ion_result.get("pka_observations") or []
    for obs in observations:
        if isinstance(obs, PkaObservation):
            entries.append(
                PkaEntry(
                    pka_value=obs.pka_value,
                    source=obs.source,
                    source_record_id=obs.source_record_id,
                    evidence_type=obs.evidence_type,
                    site_type=obs.site_type,
                    site_label=obs.site_label,
                    raw_text=obs.raw_text,
                    confidence_score=obs.confidence_score,
                )
            )
    if entries:
        return entries

    # backward compatibility fallback
    for source_key, stype in (("acidic_pka_list", "acid"), ("basic_pka_list", "base")):
        raw = ion_result.get(source_key)
        if not raw:
            continue
        for token in str(raw).split(";"):
            token = token.strip()
            if not token:
                continue
            try:
                val = float(token)
            except ValueError:
                continue
            entries.append(PkaEntry(pka_value=val, source=ion_result.get("pka_source", "unknown"), site_type=stype))
    return entries


def build_pka_detail_rows(entries: list[PkaEntry], ph: float, logp: float, identity: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for entry in entries:
        compute_per_pka(entry, ph, logp)
        rows.append(
            {
                "name": identity.get("name"),
                "input_name": identity.get("input_name"),
                "input_type": identity.get("input_type"),
                "matched_name": identity.get("matched_name"),
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
        )
    return rows
