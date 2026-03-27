from __future__ import annotations

import logging
import os
from functools import lru_cache

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

log = logging.getLogger(__name__)
_HTTP_TIMEOUT = float(os.getenv("EPIDERMAL_HTTP_TIMEOUT", "8.0"))
_DISABLE_LIVE = os.getenv("EPIDERMAL_DISABLE_LIVE_LOOKUP", "0") == "1"

_HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (compatible; PermeabilityScreen/0.3; "
        "+https://github.com/djeckins/Permeability-test)"
    ),
    "Accept-Language": "en-US,en;q=0.9",
}


@lru_cache(maxsize=1)
def get_session() -> requests.Session:
    session = requests.Session()
    retries = Retry(
        total=2,
        backoff_factor=0.4,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET",),
    )
    adapter = HTTPAdapter(max_retries=retries)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def get_json(url: str) -> dict | list | None:
    if _DISABLE_LIVE:
        log.info("Live lookup disabled for URL: %s", url)
        return None
    try:
        resp = get_session().get(url, headers=_HEADERS, timeout=_HTTP_TIMEOUT)
        if not resp.ok:
            log.warning("HTTP %s for %s", resp.status_code, url)
            return None
        return resp.json()
    except Exception as exc:
        log.warning("JSON request failed for %s: %s", url, exc)
        return None


def get_text(url: str) -> str | None:
    if _DISABLE_LIVE:
        log.info("Live lookup disabled for URL: %s", url)
        return None
    try:
        resp = get_session().get(url, headers=_HEADERS, timeout=_HTTP_TIMEOUT)
        if not resp.ok:
            log.warning("HTTP %s for %s", resp.status_code, url)
            return None
        return resp.text
    except Exception as exc:
        log.warning("Text request failed for %s: %s", url, exc)
        return None
