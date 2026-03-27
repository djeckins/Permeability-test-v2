from __future__ import annotations

import logging
from dataclasses import dataclass
from functools import lru_cache
from urllib.parse import quote_plus

from epidermal_barrier_screen.adapters.common import get_json
from epidermal_barrier_screen.predictors.pka_predictor import PkaObservation

log = logging.getLogger(__name__)
_API = "https://www.ebi.ac.uk/chembl/api/data"


@dataclass
class ChemblResult:
    chembl_id: str | None
    pref_name: str | None
    pka_observations: list[PkaObservation]


def _extract(mol_data: dict) -> ChemblResult:
    props = mol_data.get("molecule_properties") or {}
    pka_val = props.get("acd_pka")
    observations: list[PkaObservation] = []
    try:
        if pka_val is not None:
            observations.append(
                PkaObservation(
                    pka_value=float(pka_val),
                    source="chembl",
                    source_record_id=mol_data.get("molecule_chembl_id"),
                    evidence_type="computed",
                    site_type="unknown",
                    raw_text=f"acd_pka={pka_val}",
                    confidence_score="moderate",
                )
            )
    except (TypeError, ValueError):
        pass
    return ChemblResult(
        chembl_id=mol_data.get("molecule_chembl_id"),
        pref_name=mol_data.get("pref_name"),
        pka_observations=observations,
    )


@lru_cache(maxsize=512)
def lookup(canonical_smiles: str | None, name: str | None, inchikey: str | None) -> ChemblResult:
    if inchikey:
        data = get_json(f"{_API}/molecule/{inchikey}.json")
        if isinstance(data, dict):
            result = _extract(data)
            if result.pka_observations:
                return result

    if canonical_smiles:
        url = f"{_API}/molecule.json?molecule_structures__canonical_smiles__flexmatch={quote_plus(canonical_smiles)}&limit=1"
        data = get_json(url)
        if isinstance(data, dict) and data.get("molecules"):
            result = _extract(data["molecules"][0])
            if result.pka_observations:
                return result

    if name:
        data = get_json(f"{_API}/molecule/search.json?q={quote_plus(name)}&limit=3")
        if isinstance(data, dict):
            for mol in data.get("molecules", []):
                result = _extract(mol)
                if result.pka_observations:
                    return result

    return ChemblResult(chembl_id=None, pref_name=None, pka_observations=[])
