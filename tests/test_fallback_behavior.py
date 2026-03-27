from __future__ import annotations

from epidermal_barrier_screen.io import parse_input
from epidermal_barrier_screen.screen import screen_records


def test_db_empty_triggers_qupkake(monkeypatch) -> None:
    from epidermal_barrier_screen import ionization

    monkeypatch.setattr(ionization.pubchem_adapter, "lookup_pka", lambda *args, **kwargs: (None, []))
    monkeypatch.setattr(
        ionization.drugbank_adapter,
        "lookup",
        lambda *args, **kwargs: ionization.drugbank_adapter.DrugBankResult("no_match", None, None, None, None, None, []),
    )
    monkeypatch.setattr(
        ionization.chembl_adapter,
        "lookup",
        lambda *args, **kwargs: ionization.chembl_adapter.ChemblResult(None, None, []),
    )
    monkeypatch.setattr(
        ionization,
        "predict_qupkake",
        lambda _smiles: [
            ionization.PkaObservation(
                pka_value=4.2,
                source="qupkake",
                evidence_type="ML_predicted",
                site_type="acid",
                confidence_score="moderate",
            )
        ],
    )

    records = parse_input("smiles", "OC(=O)c1ccccc1")
    df = screen_records(records, ph=7.4)
    assert df.loc[0, "pka_source"] == "qupkake"
    assert bool(df.loc[0, "used_qupkake_fallback"]) is True
    assert "heuristic" not in str(df.loc[0, "pka_prediction_method"]).lower()
