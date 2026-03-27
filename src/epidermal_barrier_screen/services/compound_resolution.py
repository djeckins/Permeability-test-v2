"""Resolve compound names/synonyms to canonical SMILES via PubChem.

If the user enters a compound name (e.g. "caffeic acid"), this module
resolves it to a canonical SMILES string so the rest of the pipeline can
work from structure.  It also retrieves InChIKey, molecular formula, and
synonyms when available.
"""
from __future__ import annotations

import logging
import os
import re
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Any
from urllib.parse import quote_plus

import requests
from rdkit import Chem

log = logging.getLogger(__name__)

_HTTP_TIMEOUT = float(os.getenv("EPIDERMAL_HTTP_TIMEOUT", "8.0"))
_DISABLE_LIVE = os.getenv("EPIDERMAL_DISABLE_LIVE_LOOKUP", "0") == "1"

_PUBCHEM_PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

_HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (compatible; PermeabilityScreen/0.3; "
        "+https://github.com/djeckins/Permeability-test)"
    ),
}


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------
@dataclass
class ResolvedCompound:
    """Identity information for a resolved compound."""
    input_text: str
    input_type: str  # "smiles", "inchikey", "name"
    matched_name: str | None = None
    canonical_smiles: str | None = None
    inchikey: str | None = None
    molecular_formula: str | None = None
    pubchem_cid: int | None = None
    source_name: str | None = None  # e.g. "pubchem", "rdkit"
    source_identifier: str | None = None  # e.g. CID as string
    resolution_confidence: str = "none"  # "high", "moderate", "low", "none"
    resolution_notes: str | None = None
    synonyms: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Input type detection
# ---------------------------------------------------------------------------
_INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
# Rough: SMILES contains digits, brackets, =, #, (, ), c, C, etc. No spaces
_SMILES_CHARS = set("CNOPSFIBrcnopsfi[]()=#+-\\/@.%0123456789lae")


_SMILES_STRUCTURAL_RE = re.compile(r"[=\(\)#\[\]\\/@\d]")
_ONLY_ALPHA_RE = re.compile(r"^[a-zA-Z]+$")


def detect_input_type(text: str) -> str:
    """Classify input as 'smiles', 'inchikey', or 'name'.

    Uses fast heuristics first to avoid expensive RDKit parsing on names.
    """
    text = text.strip()
    if not text:
        return "name"
    if _INCHIKEY_RE.match(text):
        return "inchikey"
    # Space → almost certainly a name (compound names, not SMILES)
    if " " in text:
        return "name"
    # Pure alphabetic strings are almost always names, not SMILES
    # (valid SMILES like "CCO" are short; longer alpha-only = word/name)
    if _ONLY_ALPHA_RE.match(text) and len(text) > 4:
        return "name"
    # If it has structural SMILES characters (=, (, ), #, [, ], digits, etc.)
    # it's very likely SMILES
    if _SMILES_STRUCTURAL_RE.search(text):
        mol = Chem.MolFromSmiles(text)
        if mol is not None:
            return "smiles"
        return "name"
    # Short alpha-only (≤4 chars): could be short SMILES like "CCO", "CCCC"
    if all(ch in _SMILES_CHARS for ch in text):
        mol = Chem.MolFromSmiles(text)
        if mol is not None:
            return "smiles"
    return "name"


# ---------------------------------------------------------------------------
# PubChem resolution helpers
# ---------------------------------------------------------------------------
@lru_cache(maxsize=512)
def _pubchem_get_json(url: str) -> dict | None:
    """GET a PubChem URL and return parsed JSON."""
    if _DISABLE_LIVE:
        return None
    try:
        resp = requests.get(url, headers=_HEADERS, timeout=_HTTP_TIMEOUT)
        if not resp.ok:
            log.debug("PubChem HTTP %d for %s", resp.status_code, url)
            return None
        return resp.json()
    except Exception as exc:
        log.debug("PubChem request failed: %s", exc)
        return None


def _pubchem_resolve_name(name: str) -> dict[str, Any] | None:
    """Resolve a compound name to CID + canonical SMILES via PubChem PUG REST."""
    encoded = quote_plus(name.strip())
    # Get CID first
    url = f"{_PUBCHEM_PUG}/compound/name/{encoded}/cids/JSON"
    data = _pubchem_get_json(url)
    if not isinstance(data, dict):
        return None
    cids = data.get("IdentifierList", {}).get("CID", [])
    if not cids:
        return None
    cid = cids[0]

    # Get full compound properties
    prop_url = (
        f"{_PUBCHEM_PUG}/compound/cid/{cid}/property/"
        "CanonicalSMILES,IsomericSMILES,InChIKey,MolecularFormula,IUPACName/JSON"
    )
    props_data = _pubchem_get_json(prop_url)
    if not isinstance(props_data, dict):
        return None
    props_list = props_data.get("PropertyTable", {}).get("Properties", [])
    if not props_list:
        return None
    props = props_list[0]

    # Get synonyms (first few)
    syn_url = f"{_PUBCHEM_PUG}/compound/cid/{cid}/synonyms/JSON"
    syn_data = _pubchem_get_json(syn_url)
    synonyms: list[str] = []
    if isinstance(syn_data, dict):
        info_list = syn_data.get("InformationList", {}).get("Information", [])
        if info_list:
            synonyms = info_list[0].get("Synonym", [])[:10]

    return {
        "cid": cid,
        "canonical_smiles": props.get("CanonicalSMILES"),
        "isomeric_smiles": props.get("IsomericSMILES"),
        "inchikey": props.get("InChIKey"),
        "molecular_formula": props.get("MolecularFormula"),
        "iupac_name": props.get("IUPACName"),
        "synonyms": synonyms,
    }


def _pubchem_resolve_inchikey(inchikey: str) -> dict[str, Any] | None:
    """Resolve InChIKey to CID + canonical SMILES via PubChem."""
    url = f"{_PUBCHEM_PUG}/compound/inchikey/{inchikey}/cids/JSON"
    data = _pubchem_get_json(url)
    if not isinstance(data, dict):
        return None
    cids = data.get("IdentifierList", {}).get("CID", [])
    if not cids:
        return None
    cid = cids[0]

    prop_url = (
        f"{_PUBCHEM_PUG}/compound/cid/{cid}/property/"
        "CanonicalSMILES,IsomericSMILES,InChIKey,MolecularFormula,IUPACName/JSON"
    )
    props_data = _pubchem_get_json(prop_url)
    if not isinstance(props_data, dict):
        return None
    props_list = props_data.get("PropertyTable", {}).get("Properties", [])
    if not props_list:
        return None
    props = props_list[0]

    return {
        "cid": cid,
        "canonical_smiles": props.get("CanonicalSMILES"),
        "isomeric_smiles": props.get("IsomericSMILES"),
        "inchikey": props.get("InChIKey"),
        "molecular_formula": props.get("MolecularFormula"),
        "iupac_name": props.get("IUPACName"),
        "synonyms": [],
    }


# ---------------------------------------------------------------------------
# Main resolution function
# ---------------------------------------------------------------------------
def resolve_compound(text: str) -> ResolvedCompound:
    """Resolve input text (name, SMILES, or InChIKey) to a canonical structure.

    For SMILES input: canonicalize with RDKit, generate InChIKey.
    For name/InChIKey input: resolve via PubChem to get canonical SMILES.
    """
    text = text.strip()
    input_type = detect_input_type(text)
    result = ResolvedCompound(input_text=text, input_type=input_type)

    if input_type == "smiles":
        mol = Chem.MolFromSmiles(text)
        if mol is None:
            result.resolution_notes = "invalid_smiles"
            return result
        result.canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
        try:
            from rdkit.Chem import inchi as rdkit_inchi
            result.inchikey = rdkit_inchi.MolToInchiKey(mol)
        except Exception:
            pass
        result.source_name = "rdkit"
        result.resolution_confidence = "high"
        return result

    if input_type == "inchikey":
        result.inchikey = text
        pub = _pubchem_resolve_inchikey(text)
        if pub and pub.get("canonical_smiles"):
            result.canonical_smiles = pub["canonical_smiles"]
            result.pubchem_cid = pub["cid"]
            result.molecular_formula = pub.get("molecular_formula")
            result.matched_name = pub.get("iupac_name")
            result.source_name = "pubchem"
            result.source_identifier = str(pub["cid"])
            result.resolution_confidence = "high"
        else:
            result.resolution_notes = "inchikey_not_found_in_pubchem"
        return result

    # input_type == "name"
    pub = _pubchem_resolve_name(text)
    if pub and pub.get("canonical_smiles"):
        result.canonical_smiles = pub["canonical_smiles"]
        result.inchikey = pub.get("inchikey")
        result.pubchem_cid = pub["cid"]
        result.molecular_formula = pub.get("molecular_formula")
        result.matched_name = pub.get("iupac_name") or text
        result.source_name = "pubchem"
        result.source_identifier = str(pub["cid"])
        result.synonyms = pub.get("synonyms", [])

        # Confidence: check if input name appears in synonyms
        name_lower = text.lower().strip()
        syn_lower = [s.lower() for s in result.synonyms]
        if name_lower in syn_lower or name_lower == (result.matched_name or "").lower():
            result.resolution_confidence = "high"
        else:
            result.resolution_confidence = "moderate"
            result.resolution_notes = "name_matched_but_not_in_top_synonyms"
    else:
        result.resolution_notes = "name_not_found_in_pubchem"
        result.resolution_confidence = "none"

    return result
