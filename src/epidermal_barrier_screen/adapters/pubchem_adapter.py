from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from functools import lru_cache
from urllib.parse import quote_plus

from epidermal_barrier_screen.adapters.common import get_json
from epidermal_barrier_screen.predictors.pka_predictor import PkaObservation

log = logging.getLogger(__name__)
_PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_PUGVIEW = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"


@dataclass
class PubChemIdentity:
    cid: int | None
    canonical_smiles: str | None
    inchikey: str | None
    molecular_formula: str | None
    iupac_name: str | None
    synonyms: list[str]


def _as_float(val: object) -> float | None:
    try:
        return float(val)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


@lru_cache(maxsize=512)
def resolve_name(name: str) -> PubChemIdentity | None:
    encoded = quote_plus(name.strip())
    cid_data = get_json(f"{_PUG}/compound/name/{encoded}/cids/JSON")
    if not isinstance(cid_data, dict):
        log.info("PubChem name lookup returned no JSON for %s", name)
        return None
    cids = cid_data.get("IdentifierList", {}).get("CID", [])
    if not cids:
        log.info("PubChem name lookup found no CID for %s", name)
        return None
    return fetch_identity_by_cid(int(cids[0]))


@lru_cache(maxsize=512)
def fetch_identity_by_cid(cid: int) -> PubChemIdentity | None:
    props_url = (
        f"{_PUG}/compound/cid/{cid}/property/"
        "CanonicalSMILES,InChIKey,MolecularFormula,IUPACName/JSON"
    )
    props_data = get_json(props_url)
    if not isinstance(props_data, dict):
        return None
    props = (props_data.get("PropertyTable", {}).get("Properties", []) or [None])[0]
    if not isinstance(props, dict):
        return None
    syn_data = get_json(f"{_PUG}/compound/cid/{cid}/synonyms/JSON")
    synonyms: list[str] = []
    if isinstance(syn_data, dict):
        infos = syn_data.get("InformationList", {}).get("Information", [])
        if infos:
            synonyms = infos[0].get("Synonym", [])[:20]
    return PubChemIdentity(
        cid=cid,
        canonical_smiles=props.get("CanonicalSMILES"),
        inchikey=props.get("InChIKey"),
        molecular_formula=props.get("MolecularFormula"),
        iupac_name=props.get("IUPACName"),
        synonyms=synonyms,
    )


@lru_cache(maxsize=512)
def resolve_inchikey(inchikey: str) -> PubChemIdentity | None:
    data = get_json(f"{_PUG}/compound/inchikey/{inchikey}/cids/JSON")
    if not isinstance(data, dict):
        return None
    cids = data.get("IdentifierList", {}).get("CID", [])
    if not cids:
        return None
    return fetch_identity_by_cid(int(cids[0]))


@lru_cache(maxsize=1024)
def cid_from_smiles(smiles: str) -> int | None:
    data = get_json(f"{_PUG}/compound/smiles/{quote_plus(smiles)}/cids/JSON")
    if not isinstance(data, dict):
        return None
    cids = data.get("IdentifierList", {}).get("CID", [])
    return int(cids[0]) if cids else None


@lru_cache(maxsize=1024)
def fetch_pka_by_cid(cid: int) -> list[PkaObservation]:
    data = get_json(f"{_PUGVIEW}/data/compound/{cid}/JSON?heading=Dissociation+Constants")
    if not isinstance(data, dict):
        log.info("PubChem PUG-View returned no pKa payload for CID %s", cid)
        return []

    values: list[PkaObservation] = []

    def walk(sections: list[dict]) -> None:
        for sec in sections:
            for sub in sec.get("Section", []):
                if isinstance(sub, dict):
                    walk([sub])
            for info in sec.get("Information", []):
                value = info.get("Value", {})
                for swm in value.get("StringWithMarkup", []):
                    txt = swm.get("String", "")
                    if not txt:
                        continue
                    for m in re.finditer(r"(?<!\d)(-?\d{1,2}(?:\.\d+)?)(?!\d)", txt):
                        parsed = _as_float(m.group(1))
                        if parsed is None or parsed < -2 or parsed > 16:
                            continue
                        values.append(
                            PkaObservation(
                                pka_value=parsed,
                                source="pubchem",
                                source_record_id=str(cid),
                                evidence_type="experimental",
                                site_type="unknown",
                                raw_text=txt,
                                confidence_score="high",
                            )
                        )

    walk(data.get("Record", {}).get("Section", []))

    dedup: dict[float, PkaObservation] = {}
    for item in values:
        key = round(item.pka_value, 3)
        dedup.setdefault(key, item)
    return list(dedup.values())


def lookup_pka(canonical_smiles: str | None, name: str | None, inchikey: str | None) -> tuple[int | None, list[PkaObservation]]:
    cid: int | None = None
    if inchikey:
        identity = resolve_inchikey(inchikey)
        cid = identity.cid if identity else None
    if cid is None and canonical_smiles:
        cid = cid_from_smiles(canonical_smiles)
    if cid is None and name:
        identity = resolve_name(name)
        cid = identity.cid if identity else None
    if cid is None:
        return None, []
    return cid, fetch_pka_by_cid(cid)
