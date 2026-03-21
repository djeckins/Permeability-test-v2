from __future__ import annotations

"""Ionization analysis: DrugBank lookup → site-aware pKa heuristic → Dimorphite-DL → pH-specific ionization → logD.

Pipeline
--------
1. Canonicalize molecule (caller provides canonical SMILES).
2. Try exact DrugBank lookup (name → InChIKey → canonical SMILES).
3. If no reliable match → enhanced site-aware heuristic pKa predictor.
4. Run Dimorphite-DL around user pH for protonation-state enumeration.
5. Decide whether one representative pKa is chemically meaningful.
6. Compute pH-specific ionization (fraction unionized, net charge).
7. Compute pH-specific logD.
8. Return full provenance metadata.

Notes
-----
- pH is always passed explicitly by the caller; never hardcoded.
- DrugBank matching is conservative: only exact or highly reliable matches are used.
- For polyphenols / multiprotic molecules, pKa column is left blank.
"""

from dataclasses import dataclass, field
from functools import lru_cache
from typing import Any, Literal
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
# DrugBank configuration
# ---------------------------------------------------------------------------
_DRUGBANK_SEARCH_URL = "https://go.drugbank.com/unearth/q?query={query}&searcher=drugs"
_DRUGBANK_BASE_URL = "https://go.drugbank.com"
_HTTP_TIMEOUT = float(os.getenv("EPIDERMAL_DRUGBANK_TIMEOUT", "8.0"))
_MAX_CANDIDATES = int(os.getenv("EPIDERMAL_DRUGBANK_MAX_CANDIDATES", "5"))
_DISABLE_LIVE_LOOKUP = os.getenv("EPIDERMAL_DISABLE_DRUGBANK_LOOKUP", "0") == "1"
_HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (compatible; EpidermalBarrierScreen/0.3; "
        "+https://github.com/djeckins/Permeability-test)"
    ),
    "Accept-Language": "en-US,en;q=0.9",
}

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
# Part 1 — DrugBank lookup
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
    if _DISABLE_LIVE_LOOKUP:
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
# Part 1+2 — pKa source resolution (input → DrugBank → heuristic)
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
    """Return acidic/basic pKa lists, source tag, and DrugBank metadata.

    Returns
    -------
    acidic_list, basic_list, pka_source, drugbank_meta
    """
    drugbank_meta: dict[str, Any] = {
        "drugbank_match_status": "not_attempted",
        "drugbank_name": None,
        "drugbank_url": None,
        "acidic_pka_drugbank": None,
        "basic_pka_drugbank": None,
        "physiological_charge_drugbank": None,
    }

    # Priority 0: user-provided pKa values from input file
    if input_pka_acidic is not None or input_pka_basic is not None:
        acidic = [input_pka_acidic] if input_pka_acidic is not None else []
        basic = [input_pka_basic] if input_pka_basic is not None else []
        if input_pka is not None and not acidic and not basic:
            ion_class = classify_ionization(sites)
            if ion_class == "base":
                basic = [input_pka]
            else:
                acidic = [input_pka]
        return acidic, basic, "input", drugbank_meta

    if input_pka is not None:
        ion_class = classify_ionization(sites)
        if ion_class == "base":
            return [], [input_pka], "input", drugbank_meta
        return [input_pka], [], "input", drugbank_meta

    # Priority 1: Live DrugBank lookup
    entry = _live_drugbank_lookup(
        canonical_smiles=canonical_smiles,
        name=name,
        inchikey=inchikey,
    )

    match_status = _drugbank_match_status(entry)
    drugbank_meta["drugbank_match_status"] = match_status

    if entry is not None:
        drugbank_meta["drugbank_name"] = entry.get("name")
        drugbank_meta["drugbank_url"] = entry.get("source_url")
        drugbank_meta["acidic_pka_drugbank"] = entry.get("strongest_acidic_pka")
        drugbank_meta["basic_pka_drugbank"] = entry.get("strongest_basic_pka")
        drugbank_meta["physiological_charge_drugbank"] = entry.get("physiological_charge")

    if match_status == "exact":
        # Use DrugBank pKa values
        acidic: list[float] = []
        basic: list[float] = []
        acid_pka = _safe_float(entry.get("strongest_acidic_pka"))  # type: ignore[union-attr]
        base_pka = _safe_float(entry.get("strongest_basic_pka"))  # type: ignore[union-attr]
        if acid_pka is not None:
            acidic.append(acid_pka)
        if base_pka is not None:
            basic.append(base_pka)
        if acidic or basic:
            return acidic, basic, "drugbank_live", drugbank_meta

    # If DrugBank match is uncertain, do NOT use its pKa — fall through to heuristic
    # but keep the metadata for provenance

    # Priority 2: Site-aware heuristic pKa prediction
    acidic = [s.heuristic_pka for s in sites if s.ion_type == "acid"]
    basic = [s.heuristic_pka for s in sites if s.ion_type == "base"]
    return acidic, basic, "heuristic_site_rules", drugbank_meta


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
    if pka_source == "drugbank_live":
        return "high"

    # Heuristic site rules
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
    DrugBank → heuristic pKa → Dimorphite-DL → ionization → logD
    """
    try:
        inchikey = inchi.MolToInchiKey(mol)
    except Exception:
        inchikey = None

    sites = detect_ionizable_sites(mol)
    ion_class = classify_ionization(sites)

    # --- pKa resolution ---
    acidic_pkas, basic_pkas, pka_source, drugbank_meta = _site_pka_lists_from_source(
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
    if drugbank_meta.get("drugbank_match_status") == "uncertain":
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
        # Low confidence heuristic pKa → compute but mark as uncertain
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

        # Still compute logD but mark as uncertain
        logd = round(clogp + math.log10(max(fraction_unionized, 1e-12)), 4)
        logd_method = "pKa_corrected_low_confidence"
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
        "lookup_match_name": drugbank_meta.get("drugbank_name") or name,
        # DrugBank metadata
        "drugbank_match_status": drugbank_meta["drugbank_match_status"],
        "drugbank_name": drugbank_meta["drugbank_name"],
        "drugbank_url": drugbank_meta["drugbank_url"],
        "acidic_pka_drugbank": drugbank_meta["acidic_pka_drugbank"],
        "basic_pka_drugbank": drugbank_meta["basic_pka_drugbank"],
        "physiological_charge_drugbank": drugbank_meta["physiological_charge_drugbank"],
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
