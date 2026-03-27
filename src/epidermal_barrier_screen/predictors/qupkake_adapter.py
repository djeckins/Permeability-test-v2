from __future__ import annotations

import logging
import os
from functools import lru_cache

from epidermal_barrier_screen.predictors.pka_predictor import PkaObservation

log = logging.getLogger(__name__)
_DISABLE_QUPKAKE = os.getenv("EPIDERMAL_DISABLE_QUPKAKE", "0") == "1"


@lru_cache(maxsize=1)
def qupkake_available() -> bool:
    if _DISABLE_QUPKAKE:
        return False
    try:
        import qupkake  # noqa: F401

        return True
    except ImportError:
        return False


def predict_qupkake(canonical_smiles: str) -> list[PkaObservation]:
    """Predict pKa with QupKake only.

    Raises RuntimeError when QupKake is unavailable so callers can log transparently.
    """

    if not canonical_smiles:
        return []
    if not qupkake_available():
        raise RuntimeError("QupKake is not available in this environment")

    from qupkake import Predictor

    predictor = Predictor()
    raw_results = predictor.predict(canonical_smiles)
    if raw_results is None:
        return []

    observations: list[PkaObservation] = []
    for site in raw_results:
        pka_val = None
        site_type = "unknown"
        site_label = None
        if isinstance(site, dict):
            pka_val = site.get("pka", site.get("pKa"))
            site_type = str(site.get("type", "unknown")).lower()
            site_label = site.get("label")
        else:
            pka_val = getattr(site, "pka", None)
            site_type = str(getattr(site, "type", "unknown")).lower()
            site_label = getattr(site, "label", None)

        try:
            numeric = float(pka_val)
        except (TypeError, ValueError):
            continue

        normalized_type = "unknown"
        if site_type in {"acid", "acidic", "deprotonation"}:
            normalized_type = "acid"
        elif site_type in {"base", "basic", "protonation"}:
            normalized_type = "base"

        observations.append(
            PkaObservation(
                pka_value=numeric,
                source="qupkake",
                source_record_id=canonical_smiles,
                evidence_type="ML_predicted",
                site_type=normalized_type,
                site_label=site_label,
                raw_text=str(site),
                confidence_score="moderate",
            )
        )
    return observations
