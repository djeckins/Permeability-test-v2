from __future__ import annotations

"""Ionization analysis: ChEMBL/PubChem → DrugBank → QupKake → heuristic → Dimorphite-DL → ionization → logD.

Pipeline
--------
1. Canonicalize molecule (caller provides canonical SMILES).
2. Try ChEMBL API lookup (by InChIKey / canonical SMILES / name).
3. If no ChEMBL match → try PubChem PUG-View (Dissociation Constants).
4. If no PubChem match → try DrugBank web lookup.
5. If no DrugBank match → try QupKake ML-based site-aware pKa predictor.
6. If QupKake unavailable → site-aware heuristic pKa predictor (SMARTS fallback).
7. Run Dimorphite-DL around user pH for protonation-state enumeration.
8. Decide whether one representative pKa is chemically meaningful.
9. Compute pH-specific ionization (fraction unionized, net charge).
10. Compute pH-specific logD.
11. Return full provenance metadata.

Notes
-----
- pH is always passed explicitly by the caller; never hardcoded.
- ChEMBL and PubChem are queried first as structured database APIs.
- DrugBank matching is conservative: only exact or highly reliable matches are used.
- QupKake provides site-level micro-pKa predictions (acidic & basic).
- For polyphenols / multiprotic molecules, pKa column is left blank.
- When pKa is uncertain, ionization and logD are left blank rather than faked.
"""

from dataclasses import dataclass, field
from functools import lru_cache
from typing import Any, Literal
import json as _json
import logging
import math
import os
import re
from urllib.parse import quote_plus, urljoin

import requests
from bs4 import BeautifulSoup
from rdkit import Chem
from rdkit.Chem import inchi

log = logging.getLogger(__name__)

DEFAULT_PH: float = 5.5
IonType = Literal["acid", "base"]

# ---------------------------------------------------------------------------
# HTTP / API configuration
# ---------------------------------------------------------------------------
_HTTP_TIMEOUT = float(os.getenv("EPIDERMAL_HTTP_TIMEOUT", "8.0"))
_DISABLE_LIVE_LOOKUP = os.getenv("EPIDERMAL_DISABLE_LIVE_LOOKUP", "0") == "1"
# Legacy env-var alias keeps backward compat
if os.getenv("EPIDERMAL_DISABLE_DRUGBANK_LOOKUP", "0") == "1":
    _DISABLE_LIVE_LOOKUP = True

_HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (compatible; PermeabilityScreen/0.3; "
        "+https://github.com/djeckins/Permeability-test)"
    ),
    "Accept-Language": "en-US,en;q=0.9",
}

# ---------------------------------------------------------------------------
# ChEMBL API configuration
# ---------------------------------------------------------------------------
_CHEMBL_API_BASE = "https://www.ebi.ac.uk/chembl/api/data"
_DISABLE_CHEMBL = os.getenv("EPIDERMAL_DISABLE_CHEMBL", "0") == "1"

# ---------------------------------------------------------------------------
# PubChem API configuration
# ---------------------------------------------------------------------------
_PUBCHEM_PUG_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_PUBCHEM_PUGVIEW_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"
_DISABLE_PUBCHEM = os.getenv("EPIDERMAL_DISABLE_PUBCHEM", "0") == "1"

# ---------------------------------------------------------------------------
# DrugBank configuration
# ---------------------------------------------------------------------------
_DRUGBANK_SEARCH_URL = "https://go.drugbank.com/unearth/q?query={query}&searcher=drugs"
_DRUGBANK_BASE_URL = "https://go.drugbank.com"
_MAX_CANDIDATES = int(os.getenv("EPIDERMAL_DRUGBANK_MAX_CANDIDATES", "5"))
_DISABLE_DRUGBANK = os.getenv("EPIDERMAL_DISABLE_DRUGBANK", "0") == "1"

# Match quality thresholds for DrugBank
_DRUGBANK_EXACT_THRESHOLD = 90   # score >= 90 → exact
_DRUGBANK_ACCEPT_THRESHOLD = 60  # score >= 60 → accepted; < 60 → rejected


# ---------------------------------------------------------------------------
# Ionizable site dataclass
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class IonizableSite:
    name: str
    ion_type: IonType
    atom_indices: tuple[int, ...]
    heuristic_pka: float
    source: str = "heuristic"


# ---------------------------------------------------------------------------
# SMARTS patterns for site-aware heuristic pKa prediction
# ---------------------------------------------------------------------------
_ACID_PATTERNS: list[tuple[str, str, float]] = [
    ("[CX3](=O)[OX2H1]", "carboxylic_acid", 4.2),
    ("[PX4](=O)([OX2H1])[OX2H1,OX1-]", "phosphate", 2.1),
    ("[SX4](=O)(=O)[OX2H1]", "sulfonic_acid", -1.0),
    ("[c][OX2H1]", "phenol", 9.9),
    ("[SX2H1]", "thiol", 10.4),
    ("[NX3][SX4](=O)(=O)[#6]", "sulfonamide", 9.5),
]

_BASE_PATTERNS: list[tuple[str, str, float]] = [
    ("[NX3H2;!$(NC=O);!$(NS(=O)(=O));!$(N[#6]=[!#6])][#6]", "primary_amine", 10.6),
    ("[NX3H1;!$(NC=O);!$(NS(=O)(=O));!$([nH]);!$(N[#6]=[!#6])]([#6])[#6]", "secondary_amine", 10.8),
    ("[NX3H0;!$(NC=O);!$(NS(=O)(=O));!$(N[#6]=[!#6])]([#6])([#6])[#6]", "tertiary_amine", 9.8),
    ("[nX2H0;r5,r6]", "pyridine_like", 3.3),
    ("[nH]1ccnc1", "imidazole_like", 6.9),
    ("[$([CX3](=[NX2])([NX3])[NX3])]", "guanidine", 13.5),
    ("[$([CX3](=[NX2])[NX3])]", "amidine", 11.5),
    ("[NX3H2]c1ccccc1", "aniline_like", 5.0),
]


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------
def _safe_float(value: Any) -> float | None:
    if value in (None, "", "None", "nan"):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _parse_pka_list(value: Any) -> list[float]:
    if value in (None, "", "None", "nan"):
        return []
    if isinstance(value, (int, float)):
        return [float(value)]
    raw = str(value).replace("|", ";").replace(",", ";")
    vals: list[float] = []
    for part in raw.split(";"):
        part = part.strip()
        if not part:
            continue
        try:
            vals.append(float(part))
        except ValueError:
            continue
    return vals


def _normalize_name(value: str | None) -> str:
    return re.sub(r"[^a-z0-9]+", "", (value or "").lower())


def _normalize_text(value: str | None) -> str:
    return re.sub(r"\s+", " ", (value or "").strip().lower())


# ---------------------------------------------------------------------------
# Compiled SMARTS patterns (cached)
# ---------------------------------------------------------------------------
@lru_cache(maxsize=1)
def _compiled_patterns() -> list[tuple[Chem.Mol, str, IonType, float]]:
    compiled: list[tuple[Chem.Mol, str, IonType, float]] = []
    for smarts, name, pka in _ACID_PATTERNS:
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None:
            compiled.append((patt, name, "acid", pka))
    for smarts, name, pka in _BASE_PATTERNS:
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None:
            compiled.append((patt, name, "base", pka))
    return compiled


# ---------------------------------------------------------------------------
# Generic HTTP helper
# ---------------------------------------------------------------------------
@lru_cache(maxsize=512)
def _fetch_json(url: str) -> dict | list | None:
    """GET *url* and parse JSON. Returns None on any error."""
    if _DISABLE_LIVE_LOOKUP:
        return None
    try:
        resp = requests.get(url, headers=_HEADERS, timeout=_HTTP_TIMEOUT)
        if not resp.ok:
            return None
        return resp.json()
    except Exception:
        return None


# ---------------------------------------------------------------------------
# ChEMBL API lookup
# ---------------------------------------------------------------------------
@lru_cache(maxsize=256)
def _chembl_lookup_by_inchikey(inchikey: str) -> dict[str, Any] | None:
    """Query ChEMBL molecule endpoint by standard InChIKey."""
    if _DISABLE_CHEMBL or _DISABLE_LIVE_LOOKUP or not inchikey:
        return None
    url = f"{_CHEMBL_API_BASE}/molecule/{inchikey}.json"
    data = _fetch_json(url)
    if not isinstance(data, dict):
        return None
    return _extract_chembl_pka(data)


@lru_cache(maxsize=256)
def _chembl_search_by_smiles(canonical_smiles: str) -> dict[str, Any] | None:
    """Search ChEMBL by canonical SMILES (flexmatch)."""
    if _DISABLE_CHEMBL or _DISABLE_LIVE_LOOKUP or not canonical_smiles:
        return None
    encoded = quote_plus(canonical_smiles)
    url = (
        f"{_CHEMBL_API_BASE}/molecule.json"
        f"?molecule_structures__canonical_smiles__flexmatch={encoded}&limit=1"
    )
    data = _fetch_json(url)
    if not isinstance(data, dict):
        return None
    mols = data.get("molecules", [])
    if not mols:
        return None
    return _extract_chembl_pka(mols[0])


@lru_cache(maxsize=256)
def _chembl_search_by_name(name: str) -> dict[str, Any] | None:
    """Search ChEMBL molecule by name."""
    if _DISABLE_CHEMBL or _DISABLE_LIVE_LOOKUP or not name:
        return None
    encoded = quote_plus(name.strip())
    url = f"{_CHEMBL_API_BASE}/molecule/search.json?q={encoded}&limit=3"
    data = _fetch_json(url)
    if not isinstance(data, dict):
        return None
    mols = data.get("molecules", [])
    if not mols:
        return None
    # Pick first result whose pref_name matches closely
    name_norm = _normalize_name(name)
    for mol_entry in mols:
        pref = mol_entry.get("pref_name") or ""
        if _normalize_name(pref) == name_norm:
            result = _extract_chembl_pka(mol_entry)
            if result:
                return result
    # Fallback to first result if it has pKa
    return _extract_chembl_pka(mols[0])


def _extract_chembl_pka(mol_data: dict) -> dict[str, Any] | None:
    """Extract pKa fields from a ChEMBL molecule JSON object."""
    props = mol_data.get("molecule_properties") or {}
    acd_pka = _safe_float(props.get("acd_pka"))
    acd_logp = _safe_float(props.get("acd_logp"))

    chembl_id = mol_data.get("molecule_chembl_id")
    pref_name = mol_data.get("pref_name")

    if acd_pka is None:
        return None

    # ChEMBL acd_pka is the "most acidic pKa" from ACD/Labs or ChemAxon.
    # It can be acidic or basic depending on the molecule.
    # We store it and let the caller decide.
    return {
        "chembl_id": chembl_id,
        "chembl_name": pref_name,
        "chembl_acd_pka": acd_pka,
        "chembl_acd_logp": acd_logp,
    }


def _chembl_lookup(
    *,
    canonical_smiles: str | None,
    name: str | None,
    inchikey: str | None,
) -> dict[str, Any] | None:
    """Try ChEMBL lookup: InChIKey → SMILES → name."""
    if _DISABLE_CHEMBL or _DISABLE_LIVE_LOOKUP:
        return None

    # InChIKey is the most reliable identifier
    if inchikey:
        result = _chembl_lookup_by_inchikey(inchikey)
        if result:
            return result

    if canonical_smiles:
        result = _chembl_search_by_smiles(canonical_smiles)
        if result:
            return result

    if name and name.strip():
        result = _chembl_search_by_name(name)
        if result:
            return result

    return None


# ---------------------------------------------------------------------------
# PubChem API lookup
# ---------------------------------------------------------------------------
@lru_cache(maxsize=256)
def _pubchem_get_cid_by_inchikey(inchikey: str) -> int | None:
    """Resolve InChIKey → PubChem CID."""
    if not inchikey:
        return None
    url = f"{_PUBCHEM_PUG_BASE}/compound/inchikey/{inchikey}/cids/JSON"
    data = _fetch_json(url)
    if not isinstance(data, dict):
        return None
    cids = data.get("IdentifierList", {}).get("CID", [])
    return cids[0] if cids else None


@lru_cache(maxsize=256)
def _pubchem_get_cid_by_smiles(smiles: str) -> int | None:
    """Resolve canonical SMILES → PubChem CID."""
    if not smiles:
        return None
    encoded = quote_plus(smiles)
    url = f"{_PUBCHEM_PUG_BASE}/compound/smiles/{encoded}/cids/JSON"
    data = _fetch_json(url)
    if not isinstance(data, dict):
        return None
    cids = data.get("IdentifierList", {}).get("CID", [])
    return cids[0] if cids else None


@lru_cache(maxsize=256)
def _pubchem_get_cid_by_name(name: str) -> int | None:
    """Resolve compound name → PubChem CID."""
    if not name:
        return None
    encoded = quote_plus(name.strip())
    url = f"{_PUBCHEM_PUG_BASE}/compound/name/{encoded}/cids/JSON"
    data = _fetch_json(url)
    if not isinstance(data, dict):
        return None
    cids = data.get("IdentifierList", {}).get("CID", [])
    return cids[0] if cids else None


@lru_cache(maxsize=256)
def _pubchem_get_pka_from_pugview(cid: int) -> dict[str, Any] | None:
    """Fetch pKa (Dissociation Constants) from PubChem PUG-View annotations."""
    url = (
        f"{_PUBCHEM_PUGVIEW_BASE}/data/compound/{cid}/JSON"
        f"?heading=Dissociation+Constants"
    )
    data = _fetch_json(url)
    if not isinstance(data, dict):
        return None

    # Navigate the PUG-View nested structure to find pKa text values.
    # Structure: Record → Section[] → Section[] → Information[] → Value
    pka_values: list[float] = []
    pka_texts: list[str] = []

    def _walk_sections(sections: list) -> None:
        for sec in sections:
            if not isinstance(sec, dict):
                continue
            # Recurse into sub-sections
            if "Section" in sec:
                _walk_sections(sec["Section"])
            # Extract Information entries
            for info in sec.get("Information", []):
                val = info.get("Value", {})
                # Try StringWithMarkup (most common for pKa annotations)
                for swm in val.get("StringWithMarkup", []):
                    text = swm.get("String", "")
                    if text:
                        pka_texts.append(text)
                        # Extract numeric pKa values from text
                        for m in re.finditer(
                            r"(?:pKa|pka|PKA)\s*[=:≈~]?\s*(-?\d+(?:\.\d+)?)",
                            text, re.IGNORECASE,
                        ):
                            v = _safe_float(m.group(1))
                            if v is not None:
                                pka_values.append(v)
                        # Also match bare floats that look like pKa
                        for m in re.finditer(
                            r"(?<!\d)(-?\d{1,2}\.\d{1,2})(?!\d)",
                            text,
                        ):
                            v = _safe_float(m.group(1))
                            if v is not None and -2 <= v <= 16:
                                if v not in pka_values:
                                    pka_values.append(v)
                # Try Number directly
                num = val.get("Number")
                if isinstance(num, (list,)):
                    for n in num:
                        v = _safe_float(n)
                        if v is not None:
                            pka_values.append(v)
                elif num is not None:
                    v = _safe_float(num)
                    if v is not None:
                        pka_values.append(v)

    record = data.get("Record", {})
    _walk_sections(record.get("Section", []))

    if not pka_values:
        return None

    return {
        "pubchem_cid": cid,
        "pubchem_pka_values": pka_values,
        "pubchem_pka_texts": pka_texts,
    }


def _pubchem_lookup(
    *,
    canonical_smiles: str | None,
    name: str | None,
    inchikey: str | None,
) -> dict[str, Any] | None:
    """Try PubChem pKa lookup: InChIKey → SMILES → name → PUG-View."""
    if _DISABLE_PUBCHEM or _DISABLE_LIVE_LOOKUP:
        return None

    cid: int | None = None

    if inchikey:
        cid = _pubchem_get_cid_by_inchikey(inchikey)
    if cid is None and canonical_smiles:
        cid = _pubchem_get_cid_by_smiles(canonical_smiles)
    if cid is None and name and name.strip():
        cid = _pubchem_get_cid_by_name(name)
    if cid is None:
        return None

    return _pubchem_get_pka_from_pugview(cid)


# ---------------------------------------------------------------------------
# QupKake ML pKa predictor (optional)
# ---------------------------------------------------------------------------
_DISABLE_QUPKAKE = os.getenv("EPIDERMAL_DISABLE_QUPKAKE", "0") == "1"


@lru_cache(maxsize=1)
def _qupkake_available() -> bool:
    """Check whether QupKake is importable."""
    if _DISABLE_QUPKAKE:
        return False
    try:
        import qupkake  # noqa: F401
        return True
    except ImportError:
        return False


@lru_cache(maxsize=256)
def _qupkake_predict(canonical_smiles: str) -> dict[str, Any] | None:
    """Run QupKake on a SMILES and return site-level acidic/basic pKa lists.

    QupKake returns micro-pKa values for each ionizable site (atom-level).
    We collect them into acidic and basic lists.
    """
    if not _qupkake_available():
        return None

    try:
        from qupkake import Predictor

        predictor = Predictor()
        results = predictor.predict(canonical_smiles)

        if results is None or not hasattr(results, "__iter__"):
            return None

        acidic_pkas: list[float] = []
        basic_pkas: list[float] = []

        # QupKake returns per-site predictions; iterate and collect
        for site in results:
            pka_val = None
            site_type = None

            if isinstance(site, dict):
                pka_val = _safe_float(site.get("pka") or site.get("pKa"))
                site_type = site.get("type", "").lower()
            elif hasattr(site, "pka"):
                pka_val = _safe_float(getattr(site, "pka", None))
                site_type = getattr(site, "type", "").lower()

            if pka_val is None:
                continue

            if site_type in ("acid", "acidic", "deprotonation"):
                acidic_pkas.append(pka_val)
            elif site_type in ("base", "basic", "protonation"):
                basic_pkas.append(pka_val)
            else:
                # Heuristic: pKa < 7 → likely acidic, >= 7 → likely basic
                if pka_val < 7.0:
                    acidic_pkas.append(pka_val)
                else:
                    basic_pkas.append(pka_val)

        if not acidic_pkas and not basic_pkas:
            return None

        return {
            "acidic_pkas": sorted(acidic_pkas),
            "basic_pkas": sorted(basic_pkas, reverse=True),
            "method": "qupkake",
        }
    except Exception as exc:
        log.debug("QupKake failed for %s: %s", canonical_smiles, exc)
        return None


# ---------------------------------------------------------------------------
# DrugBank web-scraping lookup
# ---------------------------------------------------------------------------
@lru_cache(maxsize=256)
def _fetch_url(url: str) -> str | None:
    if _DISABLE_LIVE_LOOKUP:
        return None
    try:
        resp = requests.get(url, headers=_HEADERS, timeout=_HTTP_TIMEOUT)
        if not resp.ok:
            return None
        return resp.text
    except Exception:
        return None


def _extract_first_float(label: str, text: str) -> float | None:
    patterns = [
        rf"{re.escape(label)}\s*[,\:]?\s*(-?\d+(?:\.\d+)?)",
        rf"{re.escape(label)}\s+(-?\d+(?:\.\d+)?)",
    ]
    for pat in patterns:
        m = re.search(pat, text, flags=re.IGNORECASE)
        if m:
            return _safe_float(m.group(1))
    return None


def _extract_drugbank_detail(html: str, query: str) -> dict[str, Any] | None:
    soup = BeautifulSoup(html, "html.parser")
    page_name = None
    h1 = soup.find("h1")
    if h1:
        page_name = h1.get_text(" ", strip=True)
    if not page_name:
        og = soup.find("meta", attrs={"property": "og:title"})
        if og and og.get("content"):
            page_name = str(og.get("content")).strip()
    if not page_name:
        title = soup.title.get_text(" ", strip=True) if soup.title else ""
        page_name = title.split(": Uses")[0].strip() if title else None

    text = soup.get_text(" ", strip=True)
    acidic = _extract_first_float("pKa (Strongest Acidic)", text)
    basic = _extract_first_float("pKa (Strongest Basic)", text)
    phys_charge = _extract_first_float("Physiological Charge", text)

    if acidic is None and basic is None and phys_charge is None:
        return None

    query_norm = _normalize_name(query)
    page_name_norm = _normalize_name(page_name)
    page_text_norm = _normalize_text(text[:12000])

    score = 0
    if query_norm and page_name_norm == query_norm:
        score = 100
    elif query_norm and query_norm in page_name_norm:
        score = 90
    elif query_norm and page_name_norm in query_norm and page_name_norm:
        score = 80
    elif query_norm and query.lower() in page_text_norm:
        score = 60

    return {
        "name": page_name,
        "strongest_acidic_pka": acidic,
        "strongest_basic_pka": basic,
        "physiological_charge": phys_charge,
        "score": score,
    }


@lru_cache(maxsize=256)
def _drugbank_search_candidates(query: str) -> tuple[str, ...]:
    if _DISABLE_LIVE_LOOKUP:
        return tuple()
    query = (query or "").strip()
    if not query:
        return tuple()

    url = _DRUGBANK_SEARCH_URL.format(query=quote_plus(query))
    html = _fetch_url(url)
    if not html:
        return tuple()

    soup = BeautifulSoup(html, "html.parser")
    urls: list[str] = []
    seen: set[str] = set()

    for a in soup.find_all("a", href=True):
        href = str(a["href"]).strip()
        if not re.match(r"^/drugs/DB[0-9A-Z]+", href):
            continue
        abs_url = urljoin(_DRUGBANK_BASE_URL, href)
        if abs_url in seen:
            continue
        seen.add(abs_url)
        urls.append(abs_url)
        if len(urls) >= _MAX_CANDIDATES:
            break

    return tuple(urls)


@lru_cache(maxsize=256)
def _live_drugbank_lookup_by_query(query: str) -> dict[str, Any] | None:
    candidates = _drugbank_search_candidates(query)
    if not candidates:
        return None

    best: dict[str, Any] | None = None
    best_score = -1

    for url in candidates:
        html = _fetch_url(url)
        if not html:
            continue
        data = _extract_drugbank_detail(html, query=query)
        if not data:
            continue
        if data["score"] > best_score:
            best_score = int(data["score"])
            best = {**data, "source_url": url, "query": query}
            if best_score >= 100:
                break

    if best is None:
        return None
    if best_score < _DRUGBANK_ACCEPT_THRESHOLD:
        return None
    return best


def _live_drugbank_lookup(
    *,
    canonical_smiles: str | None,
    name: str | None,
    inchikey: str | None,
) -> dict[str, Any] | None:
    """Try DrugBank lookup with priority: name → InChIKey → canonical SMILES."""
    if _DISABLE_DRUGBANK or _DISABLE_LIVE_LOOKUP:
        return None

    queries: list[str] = []
    if name and str(name).strip():
        queries.append(str(name).strip())
    if inchikey and str(inchikey).strip():
        queries.append(str(inchikey).strip())
    if canonical_smiles and str(canonical_smiles).strip():
        queries.append(str(canonical_smiles).strip())

    for q in queries:
        data = _live_drugbank_lookup_by_query(q)
        if data is not None:
            return data
    return None


def _drugbank_match_status(entry: dict[str, Any] | None) -> str:
    """Classify DrugBank match quality."""
    if entry is None:
        return "no_match"
    score = entry.get("score", 0)
    if score >= _DRUGBANK_EXACT_THRESHOLD:
        return "exact"
    if score >= _DRUGBANK_ACCEPT_THRESHOLD:
        return "uncertain"
    return "no_match"


# ---------------------------------------------------------------------------
# Part 2 — Site-aware heuristic pKa prediction
# ---------------------------------------------------------------------------
def detect_ionizable_sites(mol: Chem.Mol) -> list[IonizableSite]:
    sites: list[IonizableSite] = []
    seen: set[tuple[str, tuple[int, ...]]] = set()
    for patt, name, ion_type, heuristic_pka in _compiled_patterns():
        for match in mol.GetSubstructMatches(patt):
            key = (name, tuple(sorted(match)))
            if key in seen:
                continue
            seen.add(key)
            sites.append(
                IonizableSite(
                    name=name,
                    ion_type=ion_type,
                    atom_indices=tuple(match),
                    heuristic_pka=heuristic_pka,
                )
            )
    return sites


def classify_ionization(sites: list[IonizableSite]) -> str:
    if not sites:
        return "non_ionizable"
    has_acid = any(s.ion_type == "acid" for s in sites)
    has_base = any(s.ion_type == "base" for s in sites)
    if has_acid and has_base:
        return "ampholyte"
    if has_acid:
        return "acid"
    return "base"


def _count_phenol_sites(sites: list[IonizableSite]) -> int:
    return sum(1 for s in sites if s.name == "phenol")


def _has_carboxylic_acid(sites: list[IonizableSite]) -> bool:
    return any(s.name == "carboxylic_acid" for s in sites)


# ---------------------------------------------------------------------------
# Part 3 — Single representative pKa decision
# ---------------------------------------------------------------------------
def _decide_representative_pka(
    *,
    acidic_pkas: list[float],
    basic_pkas: list[float],
    ion_class: str,
    sites: list[IonizableSite],
) -> tuple[float | None, str | None]:
    """Decide whether a single representative pKa is chemically meaningful.

    Returns
    -------
    (pka_value_or_None, pka_note_or_None)
    """
    if ion_class == "non_ionizable":
        return None, None

    n_acidic = len(acidic_pkas)
    n_basic = len(basic_pkas)
    n_phenol = _count_phenol_sites(sites)
    has_cooh = _has_carboxylic_acid(sites)

    # Case D — ampholyte: both acidic and basic sites
    if n_acidic > 0 and n_basic > 0:
        return None, "multiprotic_or_no_single_representative_pKa"

    # Case D — multiple acidic sites
    if n_acidic > 1:
        # Special case: carboxylic acid + phenols → COOH dominates first ionization
        if has_cooh and n_phenol >= 1 and n_acidic == (1 + n_phenol):
            # The carboxylic acid is the strongest (lowest pKa); use it
            cooh_pkas = [s.heuristic_pka for s in sites if s.name == "carboxylic_acid"]
            if cooh_pkas:
                return round(min(acidic_pkas), 2), None
        # Polyphenol or multi-acidic: no single pKa
        if n_phenol >= 2:
            return None, "multiprotic_or_no_single_representative_pKa"
        # Other multi-acidic: no single pKa
        return None, "multiprotic_or_no_single_representative_pKa"

    # Case D — multiple basic sites
    if n_basic > 1:
        return None, "multiprotic_or_no_single_representative_pKa"

    # Case B — simple monoprotic acid
    if n_acidic == 1 and n_basic == 0:
        return round(acidic_pkas[0], 2), None

    # Case C — simple monoprotic base
    if n_basic == 1 and n_acidic == 0:
        return round(basic_pkas[0], 2), None

    return None, None


# ---------------------------------------------------------------------------
# Part 4 — Dimorphite-DL protonation states
# ---------------------------------------------------------------------------
def _run_dimorphite(canonical_smiles: str, ph: float) -> dict[str, Any]:
    """Run Dimorphite-DL around user pH and extract dominant protonation info."""
    result: dict[str, Any] = {
        "protonation_state_method": None,
        "dominant_state_pH": None,
        "dominant_charge_class_pH": None,
        "expected_net_charge_pH": None,
        "dimorphite_states": None,
    }

    try:
        from dimorphite_dl import protonate_smiles
    except ImportError:
        result["protonation_state_method"] = "dimorphite_unavailable"
        return result

    try:
        states = protonate_smiles(
            canonical_smiles,
            ph_min=max(0.0, ph - 0.5),
            ph_max=min(14.0, ph + 0.5),
            precision=1.0,
        )
    except Exception as exc:
        log.debug("Dimorphite-DL failed for %s: %s", canonical_smiles, exc)
        result["protonation_state_method"] = "dimorphite_error"
        return result

    if not states:
        result["protonation_state_method"] = "dimorphite_no_output"
        return result

    result["protonation_state_method"] = "dimorphite_dl"
    result["dimorphite_states"] = states

    # Parse charge from each state
    charges: list[int] = []
    for smi in states:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
            charges.append(charge)

    if charges:
        # Most common charge = dominant state
        from collections import Counter
        charge_counts = Counter(charges)
        dominant_charge = charge_counts.most_common(1)[0][0]
        result["dominant_state_pH"] = round(ph, 2)
        result["expected_net_charge_pH"] = dominant_charge
        result["dominant_charge_class_pH"] = _dominant_charge_class(float(dominant_charge))

    return result


# ---------------------------------------------------------------------------
# pH-specific ionization math
# ---------------------------------------------------------------------------
def _neutral_fraction_acid(pka: float, ph: float) -> float:
    try:
        ratio = 10 ** (ph - pka)
        return 1.0 / (1.0 + ratio)
    except OverflowError:
        return 0.0


def _neutral_fraction_base(pka: float, ph: float) -> float:
    try:
        ratio = 10 ** (pka - ph)
        return 1.0 / (1.0 + ratio)
    except OverflowError:
        return 0.0


def _mean_charge_acid(pka: float, ph: float) -> float:
    try:
        ratio = 10 ** (ph - pka)
        return -ratio / (1.0 + ratio)
    except OverflowError:
        return -1.0


def _mean_charge_base(pka: float, ph: float) -> float:
    try:
        ratio = 10 ** (pka - ph)
        return ratio / (1.0 + ratio)
    except OverflowError:
        return 1.0


def _dominant_charge_class(expected_net_charge: float | None) -> str | None:
    if expected_net_charge is None:
        return None
    if abs(expected_net_charge) < 0.2:
        return "neutral"
    if 0.2 <= expected_net_charge < 1.2:
        return "+1"
    if -1.2 < expected_net_charge <= -0.2:
        return "-1"
    if expected_net_charge >= 1.2:
        return ">=+2"
    return "<=-2"


# ---------------------------------------------------------------------------
# pKa source resolution: input → ChEMBL → PubChem → DrugBank → QupKake → heuristic
# ---------------------------------------------------------------------------
def _site_pka_lists_from_source(
    *,
    sites: list[IonizableSite],
    canonical_smiles: str | None,
    name: str | None,
    inchikey: str | None,
    input_pka: float | None = None,
    input_pka_acidic: float | None = None,
    input_pka_basic: float | None = None,
) -> tuple[list[float], list[float], str, dict[str, Any]]:
    """Return acidic/basic pKa lists, source tag, and database metadata.

    Cascade: input → ChEMBL → PubChem → DrugBank → QupKake → heuristic.

    Returns
    -------
    acidic_list, basic_list, pka_source, db_meta
    """
    db_meta: dict[str, Any] = {
        # ChEMBL
        "chembl_id": None,
        "chembl_name": None,
        "chembl_acd_pka": None,
        # PubChem
        "pubchem_cid": None,
        "pubchem_pka_values": None,
        # DrugBank
        "drugbank_match_status": "not_attempted",
        "drugbank_name": None,
        "drugbank_url": None,
        "acidic_pka_drugbank": None,
        "basic_pka_drugbank": None,
        "physiological_charge_drugbank": None,
    }

    # ── Priority 0: user-provided pKa values from input file ──
    if input_pka_acidic is not None or input_pka_basic is not None:
        acidic = [input_pka_acidic] if input_pka_acidic is not None else []
        basic = [input_pka_basic] if input_pka_basic is not None else []
        if input_pka is not None and not acidic and not basic:
            ion_class = classify_ionization(sites)
            if ion_class == "base":
                basic = [input_pka]
            else:
                acidic = [input_pka]
        return acidic, basic, "input", db_meta

    if input_pka is not None:
        ion_class = classify_ionization(sites)
        if ion_class == "base":
            return [], [input_pka], "input", db_meta
        return [input_pka], [], "input", db_meta

    # Helper to classify a single pKa value as acidic or basic
    def _assign_single_pka(pka_val: float) -> tuple[list[float], list[float]]:
        ion_class = classify_ionization(sites)
        if ion_class == "base":
            return [], [pka_val]
        return [pka_val], []

    # ── Priority 1: ChEMBL API ──
    chembl = _chembl_lookup(
        canonical_smiles=canonical_smiles,
        name=name,
        inchikey=inchikey,
    )
    if chembl is not None:
        db_meta["chembl_id"] = chembl.get("chembl_id")
        db_meta["chembl_name"] = chembl.get("chembl_name")
        db_meta["chembl_acd_pka"] = chembl.get("chembl_acd_pka")

        pka_val = chembl["chembl_acd_pka"]
        if pka_val is not None:
            acidic, basic = _assign_single_pka(pka_val)
            return acidic, basic, "chembl", db_meta

    # ── Priority 2: PubChem PUG-View (Dissociation Constants) ──
    pubchem = _pubchem_lookup(
        canonical_smiles=canonical_smiles,
        name=name,
        inchikey=inchikey,
    )
    if pubchem is not None:
        db_meta["pubchem_cid"] = pubchem.get("pubchem_cid")
        pka_vals = pubchem.get("pubchem_pka_values", [])
        db_meta["pubchem_pka_values"] = pka_vals

        if pka_vals:
            # PubChem may return multiple pKa values; separate into acidic/basic
            # using site detection as a guide
            ion_class = classify_ionization(sites)
            if ion_class == "base":
                return [], pka_vals, "pubchem", db_meta
            elif ion_class == "ampholyte":
                # With multiple values, assign lowest as acidic, highest as basic
                sorted_vals = sorted(pka_vals)
                mid = len(sorted_vals) // 2
                return sorted_vals[:mid] or sorted_vals[:1], sorted_vals[mid:] or [], "pubchem", db_meta
            else:
                # acid or non_ionizable with pKa data → treat as acidic
                return pka_vals, [], "pubchem", db_meta

    # ── Priority 3: DrugBank web-scraping ──
    entry = _live_drugbank_lookup(
        canonical_smiles=canonical_smiles,
        name=name,
        inchikey=inchikey,
    )

    match_status = _drugbank_match_status(entry)
    db_meta["drugbank_match_status"] = match_status

    if entry is not None:
        db_meta["drugbank_name"] = entry.get("name")
        db_meta["drugbank_url"] = entry.get("source_url")
        db_meta["acidic_pka_drugbank"] = entry.get("strongest_acidic_pka")
        db_meta["basic_pka_drugbank"] = entry.get("strongest_basic_pka")
        db_meta["physiological_charge_drugbank"] = entry.get("physiological_charge")

    if match_status == "exact":
        acidic_db: list[float] = []
        basic_db: list[float] = []
        acid_pka = _safe_float(entry.get("strongest_acidic_pka"))  # type: ignore[union-attr]
        base_pka = _safe_float(entry.get("strongest_basic_pka"))  # type: ignore[union-attr]
        if acid_pka is not None:
            acidic_db.append(acid_pka)
        if base_pka is not None:
            basic_db.append(base_pka)
        if acidic_db or basic_db:
            return acidic_db, basic_db, "drugbank_live", db_meta

    # If DrugBank match is uncertain, do NOT use its pKa — fall through
    # but keep the metadata for provenance

    # ── Priority 4: QupKake ML-based site-aware predictor ──
    if canonical_smiles:
        qupkake_result = _qupkake_predict(canonical_smiles)
        if qupkake_result is not None:
            return (
                qupkake_result["acidic_pkas"],
                qupkake_result["basic_pkas"],
                "qupkake",
                db_meta,
            )

    # ── Priority 5: SMARTS-based heuristic pKa prediction (fallback) ──
    acidic = [s.heuristic_pka for s in sites if s.ion_type == "acid"]
    basic = [s.heuristic_pka for s in sites if s.ion_type == "base"]
    return acidic, basic, "heuristic_site_rules", db_meta


# ---------------------------------------------------------------------------
# Confidence assessment
# ---------------------------------------------------------------------------
def _assess_pka_confidence(
    *,
    pka_source: str,
    ion_class: str,
    n_acidic: int,
    n_basic: int,
    n_phenol: int,
) -> str:
    """Return a confidence label for the pKa prediction."""
    if pka_source == "input":
        return "high"
    if pka_source == "pubchem":
        return "high"  # experimental annotations
    if pka_source == "drugbank_live":
        return "high"
    if pka_source == "chembl":
        return "moderate"  # ChemAxon predicted values
    if pka_source == "qupkake":
        return "moderate"  # ML-predicted, better than heuristic

    # Heuristic site rules or non-ionizable
    if ion_class == "non_ionizable":
        return "high"

    # Simple monoprotic → moderate confidence for heuristic
    if (n_acidic == 1 and n_basic == 0) or (n_basic == 1 and n_acidic == 0):
        if n_phenol == 0:
            return "moderate"
        # Single phenol is less reliable
        return "low"

    # Polyphenols
    if n_phenol >= 2:
        return "low"

    # Multiprotic / complex
    if n_acidic + n_basic > 2:
        return "low"

    return "moderate"


# ---------------------------------------------------------------------------
# Main analysis function
# ---------------------------------------------------------------------------
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
    """Analyze ionization and pH-corrected lipophilicity for one molecule.

    Implements the full pipeline:
    ChEMBL → PubChem → DrugBank → QupKake → heuristic pKa → Dimorphite-DL → ionization → logD
    """
    try:
        inchikey = inchi.MolToInchiKey(mol)
    except Exception:
        inchikey = None

    sites = detect_ionizable_sites(mol)
    ion_class = classify_ionization(sites)

    # --- pKa resolution ---
    acidic_pkas, basic_pkas, pka_source, db_meta = _site_pka_lists_from_source(
        sites=sites,
        canonical_smiles=canonical_smiles,
        name=name,
        inchikey=inchikey,
        input_pka=input_pka,
        input_pka_acidic=input_pka_acidic,
        input_pka_basic=input_pka_basic,
    )

    n_phenol = _count_phenol_sites(sites)
    pka_confidence = _assess_pka_confidence(
        pka_source=pka_source,
        ion_class=ion_class,
        n_acidic=len(acidic_pkas),
        n_basic=len(basic_pkas),
        n_phenol=n_phenol,
    )

    # --- Dimorphite-DL protonation states ---
    dimorphite_result = _run_dimorphite(canonical_smiles, ph)

    # --- Representative pKa decision ---
    predicted_pka, pka_note = _decide_representative_pka(
        acidic_pkas=acidic_pkas,
        basic_pkas=basic_pkas,
        ion_class=ion_class,
        sites=sites,
    )

    # Merge DrugBank pka_note
    if db_meta.get("drugbank_match_status") == "uncertain":
        extra = "drugbank_match_uncertain_pKa_not_used"
        pka_note = f"{pka_note}; {extra}" if pka_note else extra

    # --- pH-specific ionization ---
    ionization_status: str = "ok"

    if ion_class == "non_ionizable":
        fraction_unionized: float | None = 1.0
        fraction_ionized: float | None = 0.0
        expected_net_charge: float | None = 0.0
        logd: float | None = clogp
        logd_method = "neutral_fallback_equals_cLogP"
    elif not acidic_pkas and not basic_pkas:
        # Ionizable by site detection but no pKa values available
        fraction_unionized = None
        fraction_ionized = None
        expected_net_charge = None
        logd = None
        logd_method = "unavailable_due_to_missing_pKa"
        ionization_status = "uncertain"
    elif pka_confidence == "low" and pka_source == "heuristic_site_rules":
        # Low confidence heuristic pKa → do NOT fake ionization or logD
        # Per spec: when pKa is uncertain, leave ionization and logD blank
        fraction_unionized = None
        fraction_ionized = None
        expected_net_charge = None
        logd = None
        logd_method = "unavailable_due_to_uncertain_ionization"
        ionization_status = "uncertain"
    else:
        # Normal computation with reasonable confidence
        unionized = 1.0
        for pka in acidic_pkas:
            unionized *= _neutral_fraction_acid(pka, ph)
        for pka in basic_pkas:
            unionized *= _neutral_fraction_base(pka, ph)
        fraction_unionized = max(0.0, min(1.0, unionized))
        fraction_ionized = 1.0 - fraction_unionized

        expected_net_charge = 0.0
        for pka in acidic_pkas:
            expected_net_charge += _mean_charge_acid(pka, ph)
        for pka in basic_pkas:
            expected_net_charge += _mean_charge_base(pka, ph)

        logd = round(clogp + math.log10(max(fraction_unionized, 1e-12)), 4)
        logd_method = "pKa_corrected"

    # Override with experimental logD if pH matches
    if input_logd_7_4 is not None and abs(float(ph) - 7.4) <= 0.15:
        logd = round(float(input_logd_7_4), 4)
        logd_method = "input_experimental_logD_7_4"

    # --- Round final values ---
    expected_net_charge = None if expected_net_charge is None else round(float(expected_net_charge), 4)
    fraction_unionized = None if fraction_unionized is None else round(float(fraction_unionized), 4)
    fraction_ionized = None if fraction_ionized is None else round(float(fraction_ionized), 4)
    dominant_charge_class = _dominant_charge_class(expected_net_charge)

    # --- Build result ---
    ionizable_groups = [s.name for s in sites]

    result: dict[str, Any] = {
        "inchikey": inchikey,
        "ionization_class": ion_class,
        "ionizable_group_count": len(ionizable_groups),
        "ionizable_groups": "; ".join(ionizable_groups) if ionizable_groups else None,
        # pKa columns
        "predicted_pka": predicted_pka,
        "acidic_pka_list": "; ".join(f"{v:.2f}" for v in acidic_pkas) if acidic_pkas else None,
        "basic_pka_list": "; ".join(f"{v:.2f}" for v in basic_pkas) if basic_pkas else None,
        "pka_source": pka_source,
        "pka_prediction_method": pka_source,
        "pka_confidence": pka_confidence,
        "pka_note": pka_note,
        "lookup_match_name": (
            db_meta.get("chembl_name")
            or db_meta.get("drugbank_name")
            or name
        ),
        # ChEMBL metadata
        "chembl_id": db_meta.get("chembl_id"),
        "chembl_name": db_meta.get("chembl_name"),
        "chembl_acd_pka": db_meta.get("chembl_acd_pka"),
        # PubChem metadata
        "pubchem_cid": db_meta.get("pubchem_cid"),
        "pubchem_pka_values": (
            "; ".join(f"{v:.2f}" for v in db_meta["pubchem_pka_values"])
            if db_meta.get("pubchem_pka_values")
            else None
        ),
        # DrugBank metadata
        "drugbank_match_status": db_meta["drugbank_match_status"],
        "drugbank_name": db_meta["drugbank_name"],
        "drugbank_url": db_meta["drugbank_url"],
        "acidic_pka_drugbank": db_meta["acidic_pka_drugbank"],
        "basic_pka_drugbank": db_meta["basic_pka_drugbank"],
        "physiological_charge_drugbank": db_meta["physiological_charge_drugbank"],
        # Dimorphite-DL
        "protonation_state_method": dimorphite_result["protonation_state_method"],
        "dominant_state_pH": dimorphite_result["dominant_state_pH"],
        "dominant_charge_class_pH": dimorphite_result["dominant_charge_class_pH"],
        "expected_net_charge_pH": dimorphite_result["expected_net_charge_pH"],
        # Ionization
        "fraction_unionized": fraction_unionized,
        "fraction_ionized": fraction_ionized,
        "expected_net_charge": expected_net_charge,
        "dominant_charge_class": dominant_charge_class,
        "ionization_status": ionization_status,
        # LogD
        "logd": logd,
        "logd_method": logd_method,
    }

    return result
