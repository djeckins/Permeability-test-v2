from __future__ import annotations

import math


def logd_for_acid(logp: float, pka: float, ph: float) -> float:
    return logp - math.log10(1.0 + (10 ** (ph - pka)))


def logd_for_base(logp: float, pka: float, ph: float) -> float:
    return logp - math.log10(1.0 + (10 ** (pka - ph)))
