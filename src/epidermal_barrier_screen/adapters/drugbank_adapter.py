from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from functools import lru_cache
from urllib.parse import quote_plus, urljoin

from bs4 import BeautifulSoup

from epidermal_barrier_screen.adapters.common import get_text
from epidermal_barrier_screen.predictors.pka_predictor import PkaObservation

log = logging.getLogger(__name__)
_BASE = "https://go.drugbank.com"
_SEARCH = "https://go.drugbank.com/unearth/q?query={query}&searcher=drugs"


@dataclass
class DrugBankResult:
    match_status: str
    name: str | None
    url: str | None
    acidic_pka: float | None
    basic_pka: float | None
    physiological_charge: float | None
    pka_observations: list[PkaObservation]


def _extract_first_float(label: str, text: str) -> float | None:
    for pat in [
        rf"{re.escape(label)}\s*[,:]?\s*(-?\d+(?:\.\d+)?)",
        rf"{re.escape(label)}\s+(-?\d+(?:\.\d+)?)",
    ]:
        match = re.search(pat, text, flags=re.IGNORECASE)
        if match:
            try:
                return float(match.group(1))
            except ValueError:
                return None
    return None


@lru_cache(maxsize=512)
def _candidate_urls(query: str) -> tuple[str, ...]:
    html = get_text(_SEARCH.format(query=quote_plus(query)))
    if not html:
        return tuple()
    soup = BeautifulSoup(html, "html.parser")
    found: list[str] = []
    seen: set[str] = set()
    for a_tag in soup.find_all("a", href=True):
        href = str(a_tag["href"])
        if not re.match(r"^/drugs/DB[0-9A-Z]+", href):
            continue
        full = urljoin(_BASE, href)
        if full in seen:
            continue
        seen.add(full)
        found.append(full)
        if len(found) >= 5:
            break
    return tuple(found)


@lru_cache(maxsize=512)
def lookup(query: str) -> DrugBankResult:
    urls = _candidate_urls(query)
    if not urls:
        return DrugBankResult("no_match", None, None, None, None, None, [])

    best_score = -1
    best: DrugBankResult | None = None

    query_norm = re.sub(r"[^a-z0-9]+", "", query.lower())
    for url in urls:
        html = get_text(url)
        if not html:
            continue
        soup = BeautifulSoup(html, "html.parser")
        title = soup.find("h1")
        name = title.get_text(" ", strip=True) if title else None
        text = soup.get_text(" ", strip=True)

        acidic = _extract_first_float("pKa (Strongest Acidic)", text)
        basic = _extract_first_float("pKa (Strongest Basic)", text)
        charge = _extract_first_float("Physiological Charge", text)

        if acidic is None and basic is None and charge is None:
            continue

        name_norm = re.sub(r"[^a-z0-9]+", "", (name or "").lower())
        score = 60 if query.lower() in text.lower() else 0
        if name_norm == query_norm:
            score = 100
        elif query_norm and query_norm in name_norm:
            score = max(score, 90)

        observations: list[PkaObservation] = []
        if acidic is not None:
            observations.append(
                PkaObservation(
                    pka_value=acidic,
                    source="drugbank_live",
                    source_record_id=url,
                    evidence_type="experimental",
                    site_type="acid",
                    raw_text="pKa (Strongest Acidic)",
                    confidence_score="high" if score >= 90 else "low",
                )
            )
        if basic is not None:
            observations.append(
                PkaObservation(
                    pka_value=basic,
                    source="drugbank_live",
                    source_record_id=url,
                    evidence_type="experimental",
                    site_type="base",
                    raw_text="pKa (Strongest Basic)",
                    confidence_score="high" if score >= 90 else "low",
                )
            )

        match_status = "exact" if score >= 90 else ("uncertain" if score >= 60 else "no_match")
        candidate = DrugBankResult(match_status, name, url, acidic, basic, charge, observations)
        if score > best_score:
            best_score = score
            best = candidate

    if best is None:
        return DrugBankResult("no_match", None, None, None, None, None, [])
    if best_score < 60:
        return DrugBankResult("no_match", best.name, best.url, best.acidic_pka, best.basic_pka, best.physiological_charge, [])
    return best
