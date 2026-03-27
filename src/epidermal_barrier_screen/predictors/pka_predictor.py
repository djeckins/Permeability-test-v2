from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol


@dataclass
class PkaObservation:
    """One pKa observation with provenance metadata."""

    pka_value: float
    source: str
    source_record_id: str | None = None
    evidence_type: str = "unknown"
    site_type: str = "unknown"
    site_label: str | None = None
    raw_text: str | None = None
    confidence_score: str = "moderate"


class PkaPredictor(Protocol):
    """Protocol for pKa predictors."""

    def predict(self, canonical_smiles: str) -> list[PkaObservation]:
        """Return predicted pKa observations for a canonical SMILES."""
