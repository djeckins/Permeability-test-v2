from __future__ import annotations

from epidermal_barrier_screen.services.pka_pipeline import PkaEntry, compute_per_pka


def test_acid_pka4_ph4_neutral_half() -> None:
    entry = PkaEntry(pka_value=4.0, source="test", site_type="acid")
    compute_per_pka(entry, ph=4.0, logp=1.0)
    assert entry.neutral_fraction is not None
    assert abs(entry.neutral_fraction - 0.5) < 0.01


def test_base_pka8_ph8_neutral_half() -> None:
    entry = PkaEntry(pka_value=8.0, source="test", site_type="base")
    compute_per_pka(entry, ph=8.0, logp=1.0)
    assert entry.neutral_fraction is not None
    assert abs(entry.neutral_fraction - 0.5) < 0.01
