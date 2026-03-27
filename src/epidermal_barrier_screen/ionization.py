from __future__ import annotations

"""Ionization and pKa orchestration (DB-first, QupKake fallback only)."""

from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Literal
import logging
import math

from rdkit import Chem
from rdkit.Chem import inchi

from epidermal_barrier_screen.adapters import chembl_adapter, drugbank_adapter, pubchem_adapter
from epidermal_barrier_screen.chemistry.ionization import (
    dominant_charge_class,
    expected_charge_acid,
    expected_charge_base,
    neutral_fraction_acid,
    neutral_fraction_base,
)
from epidermal_barrier_screen.predictors.pka_predictor import PkaObservation
from epidermal_barrier_screen.predictors.qupkake_adapter import predict_qupkake

log = logging.getLogger(__name__)
DEFAULT_PH: float = 5.5
IonType = Literal["acid", "base"]


@dataclass(frozen=True)
class IonizableSite:
    name: str
    ion_type: IonType
    atom_indices: tuple[int, ...]


_ACID_PATTERNS: list[tuple[str, str]] = [
    ("[CX3](=O)[OX2H1]", "carboxylic_acid"),
    ("[PX4](=O)([OX2H1])[OX2H1,OX1-]", "phosphate"),
    ("[SX4](=O)(=O)[OX2H1]", "sulfonic_acid"),
    ("[c][OX2H1]", "phenol"),
    ("[SX2H1]", "thiol"),
    ("[NX3][SX4](=O)(=O)[#6]", "sulfonamide"),
]
_BASE_PATTERNS: list[tuple[str, str]] = [
    ("[NX3H2;!$(NC=O);!$(NS(=O)(=O));!$(N[#6]=[!#6])][#6]", "primary_amine"),
    ("[NX3H1;!$(NC=O);!$(NS(=O)(=O));!$([nH]);!$(N[#6]=[!#6])]([#6])[#6]", "secondary_amine"),
    ("[NX3H0;!$(NC=O);!$(NS(=O)(=O));!$(N[#6]=[!#6])]([#6])([#6])[#6]", "tertiary_amine"),
    ("[nX2H0;r5,r6]", "pyridine_like"),
    ("[nH]1ccnc1", "imidazole_like"),
    ("[$([CX3](=[NX2])([NX3])[NX3])]", "guanidine"),
    ("[$([CX3](=[NX2])[NX3])]", "amidine"),
    ("[NX3H2]c1ccccc1", "aniline_like"),
]


@lru_cache(maxsize=1)
def _compiled_patterns() -> list[tuple[Chem.Mol, str, IonType]]:
    compiled: list[tuple[Chem.Mol, str, IonType]] = []
    for smarts, name in _ACID_PATTERNS:
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None:
            compiled.append((patt, name, "acid"))
    for smarts, name in _BASE_PATTERNS:
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None:
            compiled.append((patt, name, "base"))
    return compiled


def detect_ionizable_sites(mol: Chem.Mol) -> list[IonizableSite]:
    sites: list[IonizableSite] = []
    seen: set[tuple[str, tuple[int, ...]]] = set()
    for patt, name, ion_type in _compiled_patterns():
        for match in mol.GetSubstructMatches(patt):
            key = (name, tuple(sorted(match)))
            if key in seen:
                continue
            seen.add(key)
            sites.append(IonizableSite(name=name, ion_type=ion_type, atom_indices=tuple(match)))
    return sites


def classify_ionization(sites: list[IonizableSite]) -> str:
    if not sites:
        return "non_ionizable"
    has_acid = any(s.ion_type == "acid" for s in sites)
    has_base = any(s.ion_type == "base" for s in sites)
    if has_acid and has_base:
        return "ampholyte"
    return "acid" if has_acid else "base"


def _site_type_from_structure(default_type: str, ion_class: str, pka_value: float) -> str:
    if default_type in {"acid", "base"}:
        return default_type
    if ion_class in {"acid", "base"}:
        return ion_class
    if ion_class == "ampholyte":
        return "acid" if pka_value < 7.0 else "base"
    return "unknown"


def _round_list(vals: list[float]) -> str | None:
    return "; ".join(f"{v:.2f}" for v in vals) if vals else None


def _calculate_representative(acidic: list[float], basic: list[float]) -> tuple[float | None, str | None]:
    if len(acidic) + len(basic) == 1:
        val = acidic[0] if acidic else basic[0]
        return round(val, 2), None
    if len(acidic) + len(basic) > 1:
        return None, "multiprotic_or_no_single_representative_pKa"
    return None, None


def _site_pka_lists_from_source(
    *,
    sites: list[IonizableSite],
    canonical_smiles: str | None,
    name: str | None,
    inchikey: str | None,
    input_pka: float | None = None,
    input_pka_acidic: float | None = None,
    input_pka_basic: float | None = None,
) -> tuple[list[float], list[float], str, dict[str, Any], list[PkaObservation]]:
    """Return all DB pKa observations, fallback to QupKake only if DB empty."""

    db_meta: dict[str, Any] = {
        "chembl_id": None,
        "chembl_name": None,
        "chembl_acd_pka": None,
        "pubchem_cid": None,
        "pubchem_pka_values": None,
        "drugbank_match_status": "no_match",
        "drugbank_name": None,
        "drugbank_url": None,
        "acidic_pka_drugbank": None,
        "basic_pka_drugbank": None,
        "physiological_charge_drugbank": None,
    }

    observations: list[PkaObservation] = []
    ion_class = classify_ionization(sites)

    if input_pka_acidic is not None:
        observations.append(PkaObservation(input_pka_acidic, "input", evidence_type="user_provided", site_type="acid", confidence_score="high"))
    if input_pka_basic is not None:
        observations.append(PkaObservation(input_pka_basic, "input", evidence_type="user_provided", site_type="base", confidence_score="high"))
    if input_pka is not None and not observations:
        observations.append(PkaObservation(input_pka, "input", evidence_type="user_provided", site_type=ion_class if ion_class in {"acid", "base"} else "unknown", confidence_score="high"))
    if observations:
        acid = [o.pka_value for o in observations if o.site_type == "acid"]
        base = [o.pka_value for o in observations if o.site_type == "base"]
        unknown = [o.pka_value for o in observations if o.site_type == "unknown"]
        if unknown:
            if ion_class == "base":
                base.extend(unknown)
            else:
                acid.extend(unknown)
        return acid, base, "input", db_meta, observations

    pubchem_cid, pubchem_obs = pubchem_adapter.lookup_pka(canonical_smiles, name, inchikey)
    db_meta["pubchem_cid"] = pubchem_cid
    if pubchem_obs:
        db_meta["pubchem_pka_values"] = [obs.pka_value for obs in pubchem_obs]
        observations.extend(pubchem_obs)

    query = name or inchikey or canonical_smiles or ""
    if query:
        db_result = drugbank_adapter.lookup(query)
        db_meta["drugbank_match_status"] = db_result.match_status
        db_meta["drugbank_name"] = db_result.name
        db_meta["drugbank_url"] = db_result.url
        db_meta["acidic_pka_drugbank"] = db_result.acidic_pka
        db_meta["basic_pka_drugbank"] = db_result.basic_pka
        db_meta["physiological_charge_drugbank"] = db_result.physiological_charge
        if db_result.pka_observations:
            observations.extend(db_result.pka_observations)

    chembl = chembl_adapter.lookup(canonical_smiles, name, inchikey)
    db_meta["chembl_id"] = chembl.chembl_id
    db_meta["chembl_name"] = chembl.pref_name
    if chembl.pka_observations:
        db_meta["chembl_acd_pka"] = chembl.pka_observations[0].pka_value
        observations.extend(chembl.pka_observations)

    source = "not_found"
    db_observations = [o for o in observations if o.source in {"pubchem", "drugbank_live", "chembl"}]

    if db_observations:
        observations = db_observations
        source = "database"
    else:
        if canonical_smiles:
            try:
                qpk = predict_qupkake(canonical_smiles)
                observations = qpk
                source = "qupkake" if qpk else "not_found"
            except Exception as exc:
                log.warning("QupKake fallback failed for %s: %s", canonical_smiles, exc)
                source = "not_found"
                observations = []

    acidic: list[float] = []
    basic: list[float] = []
    for obs in observations:
        kind = _site_type_from_structure(obs.site_type, ion_class, obs.pka_value)
        obs.site_type = kind
        if kind == "acid":
            acidic.append(obs.pka_value)
        elif kind == "base":
            basic.append(obs.pka_value)

    if source == "database":
        pka_source = "pubchem" if any(o.source == "pubchem" for o in observations) else (
            "drugbank_live" if any(o.source == "drugbank_live" for o in observations) else "chembl"
        )
    else:
        pka_source = source

    return acidic, basic, pka_source, db_meta, observations


def analyze_ionization(
    *,
    mol: Chem.Mol,
    canonical_smiles: str,
    clogp: float,
    ph: float,
    name: str | None = None,
    input_pka: float | None = None,
    input_pka_acidic: float | None = None,
    input_pka_basic: float | None = None,
    input_logd_7_4: float | None = None,
) -> dict[str, Any]:
    try:
        inchikey_value = inchi.MolToInchiKey(mol)
    except Exception:
        inchikey_value = None

    sites = detect_ionizable_sites(mol)
    ion_class = classify_ionization(sites)
    acidic_pkas, basic_pkas, pka_source, db_meta, observations = _site_pka_lists_from_source(
        sites=sites,
        canonical_smiles=canonical_smiles,
        name=name,
        inchikey=inchikey_value,
        input_pka=input_pka,
        input_pka_acidic=input_pka_acidic,
        input_pka_basic=input_pka_basic,
    )

    pka_confidence = (
        "high" if pka_source in {"input", "pubchem", "drugbank_live"}
        else "moderate" if pka_source in {"chembl", "qupkake"}
        else "none" if ion_class != "non_ionizable"
        else "high"
    )

    predicted_pka, pka_note = _calculate_representative(acidic_pkas, basic_pkas)

    warnings_list: list[str] = []
    if not (0 <= float(ph) <= 14):
        warnings_list.append("invalid_pH_range")
    if pka_source == "not_found" and ion_class != "non_ionizable":
        warnings_list.append("no_pKa_from_databases_or_qupkake")

    if ion_class == "non_ionizable":
        fraction_unionized = 1.0
        fraction_ionized = 0.0
        expected_net_charge = 0.0
        logd = clogp
        logd_method = "neutral_fallback_equals_cLogP"
        ionization_status = "ok"
    elif not (acidic_pkas or basic_pkas):
        fraction_unionized = None
        fraction_ionized = None
        expected_net_charge = None
        logd = None
        logd_method = "unavailable_due_to_missing_pKa"
        ionization_status = "uncertain"
    else:
        unionized = 1.0
        expected_net_charge = 0.0
        for pka in acidic_pkas:
            unionized *= neutral_fraction_acid(pka, ph)
            expected_net_charge += expected_charge_acid(pka, ph)
        for pka in basic_pkas:
            unionized *= neutral_fraction_base(pka, ph)
            expected_net_charge += expected_charge_base(pka, ph)

        fraction_unionized = round(max(0.0, min(1.0, unionized)), 4)
        fraction_ionized = round(1.0 - fraction_unionized, 4)
        expected_net_charge = round(expected_net_charge, 4)
        logd = round(clogp + math.log10(max(fraction_unionized, 1e-12)), 4)
        logd_method = "pKa_corrected"
        ionization_status = "ok"

    if input_logd_7_4 is not None and abs(float(ph) - 7.4) <= 0.15:
        logd = round(float(input_logd_7_4), 4)
        logd_method = "input_experimental_logD_7_4"

    dominant = dominant_charge_class(expected_net_charge)

    result: dict[str, Any] = {
        "inchikey": inchikey_value,
        "ionization_class": ion_class,
        "ionizable_group_count": len(sites),
        "ionizable_groups": "; ".join(s.name for s in sites) if sites else None,
        "predicted_pka": predicted_pka,
        "acidic_pka_list": _round_list(acidic_pkas),
        "basic_pka_list": _round_list(basic_pkas),
        "pka_source": pka_source,
        "pka_prediction_method": pka_source,
        "pka_confidence": pka_confidence,
        "pka_note": pka_note,
        "lookup_match_name": db_meta.get("chembl_name") or db_meta.get("drugbank_name") or name,
        "chembl_id": db_meta.get("chembl_id"),
        "chembl_name": db_meta.get("chembl_name"),
        "chembl_acd_pka": db_meta.get("chembl_acd_pka"),
        "pubchem_cid": db_meta.get("pubchem_cid"),
        "pubchem_pka_values": _round_list(db_meta.get("pubchem_pka_values") or []),
        "drugbank_match_status": db_meta.get("drugbank_match_status"),
        "drugbank_name": db_meta.get("drugbank_name"),
        "drugbank_url": db_meta.get("drugbank_url"),
        "acidic_pka_drugbank": db_meta.get("acidic_pka_drugbank"),
        "basic_pka_drugbank": db_meta.get("basic_pka_drugbank"),
        "physiological_charge_drugbank": db_meta.get("physiological_charge_drugbank"),
        "protonation_state_method": None,
        "dominant_state_pH": None,
        "dominant_charge_class_pH": None,
        "expected_net_charge_pH": None,
        "fraction_unionized": fraction_unionized,
        "fraction_ionized": fraction_ionized,
        "expected_net_charge": expected_net_charge,
        "dominant_charge_class": dominant,
        "ionization_status": ionization_status,
        "logd": logd,
        "logd_method": logd_method,
        "used_qupkake_fallback": pka_source == "qupkake",
        "logp_value": clogp,
        "logp_source": "rdkit_crippen",
        "logp_evidence_type": "computed",
        "pka_values": _round_list(acidic_pkas + basic_pkas),
        "warnings": "; ".join(warnings_list) if warnings_list else None,
        "pka_observations": observations,
    }
    return result
