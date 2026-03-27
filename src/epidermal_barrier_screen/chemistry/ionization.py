from __future__ import annotations

import math


def neutral_fraction_acid(pka: float, ph: float) -> float:
    try:
        return 1.0 / (1.0 + (10 ** (ph - pka)))
    except OverflowError:
        return 0.0


def neutral_fraction_base(pka: float, ph: float) -> float:
    try:
        return 1.0 / (1.0 + (10 ** (pka - ph)))
    except OverflowError:
        return 0.0


def ionized_fraction_from_neutral(neutral_fraction: float) -> float:
    return 1.0 - neutral_fraction


def expected_charge_acid(pka: float, ph: float) -> float:
    try:
        ratio = 10 ** (ph - pka)
        return -ratio / (1.0 + ratio)
    except OverflowError:
        return -1.0


def expected_charge_base(pka: float, ph: float) -> float:
    try:
        ratio = 10 ** (pka - ph)
        return ratio / (1.0 + ratio)
    except OverflowError:
        return 1.0


def dominant_charge_class(expected_net_charge: float | None) -> str | None:
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
