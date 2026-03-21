from __future__ import annotations

"""Ionization analysis with live DrugBank lookup and site-aware fallback.

Priority order for pKa / charge information:
1. User-provided pKa values from input records
2. Live lookup on public DrugBank search/detail pages
3. Internal site-aware heuristic pKa assignment from SMARTS patterns

Notes
-----
- pH is always passed explicitly by the caller.
- Live DrugBank lookup is intentionally conservative: it queries DrugBank
  by the molecule name first, then falls back to InChIKey / canonical SMILES.
- If DrugBank does not return a usable match, the code falls back to the
  internal site-aware heuristic rules.
"""

from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Literal
import math
import os
import re
from urllib.parse import quote_plus, urljoin

import requests
from bs4 import BeautifulSoup
from rdkit import Chem
from rdkit.Chem import inchi

DEFAULT_PH: float = 5.5
IonType = Literal["acid", "base"]

_DRUGBANK_SEARCH_URL = "https://go.drugbank.com/unearth/q?query={query}&searcher=drugs"
_DRUGBANK_BASE_URL = "https://go.drugbank.com"
_HTTP_TIMEOUT = float(os.getenv("EPIDERMAL_DRUGBANK_TIMEOUT", "8.0"))
_MAX_CANDIDATES = int(os.getenv("EPIDERMAL_DRUGBANK_MAX_CANDIDATES", "5"))
_DISABLE_LIVE_LOOKUP = os.getenv("EPIDERMAL_DISABLE_DRUGBANK_LOOKUP", "0") == "1"
_HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (compatible; EpidermalBarrierScreen/0.2; "
        "+https://github.com/djeckins/Permeability-test)"
    ),
    "Accept-Language": "en-US,en;q=0.9",
}


@dataclass(frozen=True)
class IonizableSite:
    name: str
    ion_type: IonType
    atom_indices: tuple[int, ...]
    heuristic_pka: float
    source: str = "heuristic"


# Typical site pKa values used only as rough fallbacks.
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
    # Be conservative: accept only reasonable matches.
    if best_score < 60:
        return None
    return best


def _live_drugbank_lookup(
    *,
    canonical_smiles: str | None,
    name: str | None,
    inchikey: str | None,
) -> dict[str, Any] | None:
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


def _site_pka_lists_from_source(
    *,
    sites: list[IonizableSite],
    canonical_smiles: str | None,
    name: str | None,
    inchikey: str | None,
    input_pka: float | None = None,
    input_pka_acidic: float | None = None,
    input_pka_basic: float | None = None,
) -> tuple[list[float], list[float], str, str | None, str | None]:
    """Return acidic/basic pKa lists and metadata.

    Returns
    -------
    acidic_list, basic_list, pka_source, matched_entry_name, lookup_notes
    """
    if input_pka_acidic is not None or input_pka_basic is not None:
        acidic = [input_pka_acidic] if input_pka_acidic is not None else []
        basic = [input_pka_basic] if input_pka_basic is not None else []
        if input_pka is not None and not acidic and not basic:
            ion_class = classify_ionization(sites)
            if ion_class == "base":
                basic = [input_pka]
            else:
                acidic = [input_pka]
        return acidic, basic, "input", name, None

    if input_pka is not None:
        ion_class = classify_ionization(sites)
        if ion_class == "base":
            return [], [input_pka], "input", name, None
        return [input_pka], [], "input", name, None

    entry = _live_drugbank_lookup(
        canonical_smiles=canonical_smiles,
        name=name,
        inchikey=inchikey,
    )
    if entry is not None:
        acidic = _parse_pka_list(entry.get("acidic_pka_list"))
        basic = _parse_pka_list(entry.get("basic_pka_list"))
        acid_pka = _safe_float(entry.get("strongest_acidic_pka"))
        base_pka = _safe_float(entry.get("strongest_basic_pka"))
        if acid_pka is not None and not acidic:
            acidic.append(acid_pka)
        if base_pka is not None and not basic:
            basic.append(base_pka)
        note = entry.get("source_url")
        if acidic or basic or _safe_float(entry.get("physiological_charge")) is not None:
            return acidic, basic, "drugbank_live", entry.get("name") or name, note

    acidic = [s.heuristic_pka for s in sites if s.ion_type == "acid"]
    basic = [s.heuristic_pka for s in sites if s.ion_type == "base"]
    return acidic, basic, "heuristic_site_rules", name, None


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
    """Analyze ionization and pH-corrected lipophilicity for one molecule."""
    try:
        inchikey = inchi.MolToInchiKey(mol)
    except Exception:
        inchikey = None

    sites = detect_ionizable_sites(mol)
    ion_class = classify_ionization(sites)

    acidic_pkas, basic_pkas, pka_source, matched_name, lookup_notes = _site_pka_lists_from_source(
        sites=sites,
        canonical_smiles=canonical_smiles,
        name=name,
        inchikey=inchikey,
        input_pka=input_pka,
        input_pka_acidic=input_pka_acidic,
        input_pka_basic=input_pka_basic,
    )

    ionizable_groups = [s.name for s in sites]
    acidic_neutral = [_neutral_fraction_acid(pka, ph) for pka in acidic_pkas]
    basic_neutral = [_neutral_fraction_base(pka, ph) for pka in basic_pkas]

    if not acidic_pkas and not basic_pkas and ion_class == "non_ionizable":
        fraction_unionized = 1.0
        fraction_ionized = 0.0
        expected_net_charge = 0.0
        logd = clogp
        logd_method = "neutral_fallback_equals_cLogP"
        predicted_pka = None
        pka_note = None
    else:
        unionized = 1.0
        for fn in acidic_neutral:
            unionized *= fn
        for fn in basic_neutral:
            unionized *= fn
        fraction_unionized = max(0.0, min(1.0, unionized)) if (acidic_pkas or basic_pkas) else None
        fraction_ionized = None if fraction_unionized is None else 1.0 - fraction_unionized

        expected_net_charge = 0.0
        for pka in acidic_pkas:
            expected_net_charge += _mean_charge_acid(pka, ph)
        for pka in basic_pkas:
            expected_net_charge += _mean_charge_base(pka, ph)

        if pka_source == "drugbank_live" and not acidic_pkas and not basic_pkas:
            # kept for completeness; current parser fills strongest acidic/basic when available
            live = _live_drugbank_lookup(canonical_smiles=canonical_smiles, name=name, inchikey=inchikey)
            phys_charge = _safe_float(live.get("physiological_charge") if live else None)
            expected_net_charge = phys_charge
            fraction_unionized = None
            fraction_ionized = None
            logd = None
            logd_method = "unavailable_due_to_missing_pKa"
        else:
            if fraction_unionized is None:
                logd = None
                logd_method = "unavailable_due_to_missing_pKa"
            else:
                logd = round(clogp + math.log10(max(fraction_unionized, 1e-12)), 4)
                logd_method = "pKa_corrected"

        all_pkas = acidic_pkas + basic_pkas
        predicted_pka = round(min(all_pkas, key=lambda x: abs(x - ph)), 2) if all_pkas else None
        pka_note = lookup_notes

    if input_logd_7_4 is not None and abs(float(ph) - 7.4) <= 0.15:
        logd = round(float(input_logd_7_4), 4)
        logd_method = "input_experimental_logD_7_4"

    expected_net_charge = None if expected_net_charge is None else round(float(expected_net_charge), 4)
    fraction_unionized = None if fraction_unionized is None else round(float(fraction_unionized), 4)
    fraction_ionized = None if fraction_ionized is None else round(float(fraction_ionized), 4)
    dominant_charge_class = _dominant_charge_class(expected_net_charge)

    return {
        "inchikey": inchikey,
        "ionization_class": ion_class,
        "ionizable_group_count": len(ionizable_groups),
        "ionizable_groups": "; ".join(ionizable_groups) if ionizable_groups else None,
        "acidic_pka_list": "; ".join(f"{v:.2f}" for v in acidic_pkas) if acidic_pkas else None,
        "basic_pka_list": "; ".join(f"{v:.2f}" for v in basic_pkas) if basic_pkas else None,
        "predicted_pka": predicted_pka,
        "pka_source": pka_source,
        "lookup_match_name": matched_name,
        "pka_note": pka_note,
        "fraction_unionized": fraction_unionized,
        "fraction_ionized": fraction_ionized,
        "expected_net_charge": expected_net_charge,
        "dominant_charge_class": dominant_charge_class,
        "logd": logd,
        "logd_method": logd_method,
    }
