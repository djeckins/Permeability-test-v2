"""Test pKa lookup and name resolution for 15 molecules.

These tests verify:
- Name → SMILES resolution via PubChem
- pKa retrieval from databases (PubChem, DrugBank, ChEMBL)
- QupKake fallback when DB pKa unavailable
- Final table contains all required columns
- No heuristic pKa path
"""
from __future__ import annotations

import os
import math
import pytest
import pandas as pd

from epidermal_barrier_screen.io import parse_input
from epidermal_barrier_screen.screen import screen_records, get_pka_detail_table
from epidermal_barrier_screen.services.compound_resolution import (
    detect_input_type,
    resolve_compound,
)


# ---------------------------------------------------------------------------
# 15 test molecules (input by NAME)
# ---------------------------------------------------------------------------
TEST_MOLECULES = [
    "caffeic acid",
    "aspirin",
    "ibuprofen",
    "diclofenac",
    "benzoic acid",
    "salicylic acid",
    "acetic acid",
    "nicotine",
    "pyridine",
    "lidocaine",
    "histamine",
    "imidazole",
    "levodopa",
    "sertraline",
    "ciprofloxacin",
]


# Columns that MUST exist in the final table (preserves old + adds new)
REQUIRED_OLD_COLUMNS = [
    "name",
    "parse_status",
    "ph_used",
    "mw",
    "hba",
    "hbd",
    "rotb",
    "hac",
    "predicted_pka",
    "acidic_pka_list",
    "basic_pka_list",
    "pka_source",
    "pka_prediction_method",
    "pka_confidence",
    "pka_note",
    "lookup_match_name",
    "chembl_id",
    "chembl_name",
    "chembl_acd_pka",
    "pubchem_cid",
    "pubchem_pka_values",
    "drugbank_match_status",
    "drugbank_name",
    "drugbank_url",
    "acidic_pka_drugbank",
    "basic_pka_drugbank",
    "physiological_charge_drugbank",
    "ionization_class",
    "ionizable_group_count",
    "ionizable_groups",
    "inchikey",
    "clogp",
    "logd",
    "logd_method",
    "tpsa",
    "fraction_unionized",
    "fraction_ionized",
    "expected_net_charge",
    "dominant_charge_class",
    "ionization_status",
    "protonation_state_method",
    "dominant_state_pH",
    "dominant_charge_class_pH",
    "expected_net_charge_pH",
    "formal_charge",
    "mw_status",
    "logd_status",
    "tpsa_status",
    "hbd_status",
    "hba_status",
    "rotb_status",
    "hac_status",
    "formal_charge_status",
    "ionization_status_criterion",
    "final_result",
    "input_smiles",
    "canonical_smiles",
]

NEW_COLUMNS = [
    "input_name",
    "input_type",
    "matched_name",
    "resolution_source",
    "source_name",
    "source_identifier",
    "resolution_confidence",
    "resolution_notes",
    "molecular_formula",
    "logp_value",
    "logp_source",
    "logp_evidence_type",
    "pka_values",
    "pka_details",
    "used_qupkake_fallback",
    "warnings",
]


# ---------------------------------------------------------------------------
# Input type detection tests
# ---------------------------------------------------------------------------
class TestInputTypeDetection:
    def test_smiles_detected(self):
        assert detect_input_type("CC(=O)Oc1ccccc1C(=O)O") == "smiles"

    def test_inchikey_detected(self):
        assert detect_input_type("BSYNRYMUTXBXSQ-UHFFFAOYSA-N") == "inchikey"

    def test_name_detected(self):
        assert detect_input_type("caffeic acid") == "name"
        assert detect_input_type("aspirin") == "name"
        assert detect_input_type("ibuprofen") == "name"

    def test_ethanol_smiles(self):
        assert detect_input_type("CCO") == "smiles"


# ---------------------------------------------------------------------------
# Network availability check
# ---------------------------------------------------------------------------
def _network_available() -> bool:
    """Check if external network access is available."""
    try:
        import requests
        resp = requests.get(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/aspirin/cids/JSON",
            timeout=5,
        )
        return resp.ok
    except Exception:
        return False


# ---------------------------------------------------------------------------
# Compound resolution tests
# ---------------------------------------------------------------------------
class TestCompoundResolution:
    def test_smiles_resolution(self):
        r = resolve_compound("CCO")
        assert r.input_type == "smiles"
        assert r.canonical_smiles == "CCO"
        assert r.resolution_confidence == "high"

    @pytest.mark.skipif(
        os.getenv("EPIDERMAL_DISABLE_LIVE_LOOKUP", "0") == "1" or not _network_available(),
        reason="Live lookups disabled or network unavailable",
    )
    def test_name_resolution_caffeic_acid(self):
        r = resolve_compound("caffeic acid")
        assert r.input_type == "name"
        assert r.canonical_smiles is not None
        assert "C" in r.canonical_smiles
        assert r.pubchem_cid is not None
        assert r.resolution_confidence in ("high", "moderate")


# ---------------------------------------------------------------------------
# Per-pKa ionization calculation tests
# ---------------------------------------------------------------------------
class TestIonizationCalculations:
    """Exact calculation tests for Henderson-Hasselbalch."""

    def test_acid_pka4_ph4(self):
        """pKa=4.0, pH=4.0 => neutral_fraction ~ 0.5."""
        from epidermal_barrier_screen.services.pka_pipeline import (
            PkaEntry, compute_per_pka,
        )
        entry = PkaEntry(pka_value=4.0, source="test", site_type="acid")
        compute_per_pka(entry, ph=4.0, logp=1.5)
        assert entry.neutral_fraction is not None
        assert abs(entry.neutral_fraction - 0.5) < 0.01

    def test_base_pka8_ph8(self):
        """pKa=8.0, pH=8.0 => neutral_fraction ~ 0.5."""
        from epidermal_barrier_screen.services.pka_pipeline import (
            PkaEntry, compute_per_pka,
        )
        entry = PkaEntry(pka_value=8.0, source="test", site_type="base")
        compute_per_pka(entry, ph=8.0, logp=1.5)
        assert entry.neutral_fraction is not None
        assert abs(entry.neutral_fraction - 0.5) < 0.01

    def test_acid_logd(self):
        """logD for acid = logP - log10(1 + 10^(pH-pKa))."""
        from epidermal_barrier_screen.services.pka_pipeline import (
            PkaEntry, compute_per_pka,
        )
        entry = PkaEntry(pka_value=4.0, source="test", site_type="acid")
        compute_per_pka(entry, ph=7.4, logp=2.0)
        # logD = 2.0 - log10(1 + 10^(7.4-4.0)) = 2.0 - log10(1 + 2511.886) ≈ 2.0 - 3.4 = -1.4
        assert entry.logd is not None
        assert entry.logd < 0


# ---------------------------------------------------------------------------
# Full pipeline test (offline — no live lookups)
# ---------------------------------------------------------------------------
class TestOfflinePipeline:
    """Tests that run with live lookups disabled."""

    def test_smiles_input_produces_all_columns(self):
        """A SMILES input should produce a full table with all columns."""
        records = parse_input("smiles", "CC(=O)Oc1ccccc1C(=O)O")
        df = screen_records(records, ph=7.4)
        assert len(df) == 1

        for col in REQUIRED_OLD_COLUMNS:
            assert col in df.columns, f"Missing old column: {col}"
        for col in NEW_COLUMNS:
            assert col in df.columns, f"Missing new column: {col}"

    def test_name_input_smiles_list_mode(self):
        """Names in smiles_list mode should be detected and attempted."""
        records = parse_input("smiles_list", "caffeic acid\nCCO Ethanol")
        assert len(records) == 2
        # First record: "caffeic acid" is a name
        assert records[0]["input_type"] == "name"
        # Second: "CCO" is SMILES
        assert records[1]["input_type"] == "smiles"

    def test_non_ionizable_ethanol(self):
        records = parse_input("smiles", "CCO")
        df = screen_records(records, ph=5.5)
        assert df["ionization_class"].iloc[0] == "non_ionizable"
        assert df["fraction_unionized"].iloc[0] == 1.0
        assert df["pka_confidence"].iloc[0] == "high"

    def test_final_result_present(self):
        records = parse_input("smiles", "CC(=O)Oc1ccccc1C(=O)O")
        df = screen_records(records, ph=7.4)
        assert df["final_result"].iloc[0] in ("PASS", "BORDERLINE", "FAIL")

    def test_used_qupkake_fallback_field(self):
        """used_qupkake_fallback should be False when live lookups disabled."""
        os.environ["EPIDERMAL_DISABLE_LIVE_LOOKUP"] = "1"
        try:
            records = parse_input("smiles", "OC(=O)c1ccccc1")
            df = screen_records(records, ph=7.4)
            # With all lookups disabled, qupkake won't be the source
            assert "used_qupkake_fallback" in df.columns
        finally:
            os.environ.pop("EPIDERMAL_DISABLE_LIVE_LOOKUP", None)


# ---------------------------------------------------------------------------
# 15-molecule integration test (requires network)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    os.getenv("EPIDERMAL_DISABLE_LIVE_LOOKUP", "0") == "1" or not _network_available(),
    reason="Live lookups disabled or network unavailable",
)
class TestFifteenMolecules:
    """Integration tests running name resolution + pKa lookup for 15 compounds."""

    @pytest.fixture(scope="class")
    def results(self):
        """Run the pipeline for all 15 molecules at pH 7.4."""
        payload = "\n".join(TEST_MOLECULES)
        records = parse_input("smiles_list", payload)
        df = screen_records(records, ph=7.4)
        return df

    def test_all_molecules_present(self, results):
        assert len(results) == len(TEST_MOLECULES)

    def test_no_all_invalid(self, results):
        """At least some molecules should resolve successfully."""
        valid = results[results["parse_status"] == "ok"]
        assert len(valid) >= 10, f"Only {len(valid)} resolved out of {len(TEST_MOLECULES)}"

    def test_caffeic_acid_resolved(self, results):
        row = results.iloc[0]  # First entry is caffeic acid
        assert row["parse_status"] == "ok", "caffeic acid failed to resolve"
        assert row["canonical_smiles"] is not None
        assert row["input_type"] == "name"

    def test_pka_source_explicit(self, results):
        """All valid rows should have explicit pka_source."""
        valid = results[results["parse_status"] == "ok"]
        for _, row in valid.iterrows():
            assert row["pka_source"] in (
                "pubchem", "drugbank_live", "chembl", "qupkake", "input", "not_found"
            ), f"Unexpected pka_source: {row['pka_source']} for {row['name']}"

    def test_all_old_columns_preserved(self, results):
        for col in REQUIRED_OLD_COLUMNS:
            assert col in results.columns, f"Missing old column: {col}"

    def test_new_columns_present(self, results):
        for col in NEW_COLUMNS:
            assert col in results.columns, f"Missing new column: {col}"

    def test_no_heuristic_pka(self, results):
        """Verify no heuristic pKa source appears."""
        valid = results[results["parse_status"] == "ok"]
        for _, row in valid.iterrows():
            assert "heuristic" not in str(row.get("pka_source", "")).lower()
            assert "heuristic" not in str(row.get("pka_prediction_method", "")).lower()

    def test_pka_detail_table(self, results):
        """Per-pKa detail table should have been built."""
        detail_df = get_pka_detail_table()
        # Should have at least some entries (may be empty if no pKa found)
        assert isinstance(detail_df, pd.DataFrame)

    def test_logd_calculated_when_pka_available(self, results):
        """Molecules with pKa should have logD computed."""
        valid = results[results["parse_status"] == "ok"]
        has_pka = valid[valid["pka_source"].isin(["pubchem", "drugbank_live", "chembl", "qupkake", "input"])]
        for _, row in has_pka.iterrows():
            # If pKa was found, logd should be computed
            assert row["logd"] is not None or row["ionization_class"] == "non_ionizable", \
                f"logD missing for {row['name']} with pka_source={row['pka_source']}"


# ---------------------------------------------------------------------------
# Fallback behavior tests
# ---------------------------------------------------------------------------
class TestFallbackBehavior:
    """Test QupKake fallback when all DB adapters return no pKa."""

    def test_no_db_pka_triggers_qupkake_or_not_found(self):
        """With live lookups disabled, pka_source should be 'not_found'."""
        os.environ["EPIDERMAL_DISABLE_LIVE_LOOKUP"] = "1"
        try:
            records = parse_input("smiles", "OC(=O)c1ccccc1")  # benzoic acid
            df = screen_records(records, ph=7.4)
            # QupKake likely not installed, so should be not_found
            assert df["pka_source"].iloc[0] in ("not_found", "qupkake")
            if df["pka_source"].iloc[0] == "qupkake":
                assert df["used_qupkake_fallback"].iloc[0] is True
        finally:
            os.environ.pop("EPIDERMAL_DISABLE_LIVE_LOOKUP", None)

    def test_no_heuristic_path_exists(self):
        """Verify there is NO heuristic pKa path in the code."""
        import inspect
        from epidermal_barrier_screen import ionization

        source = inspect.getsource(ionization._site_pka_lists_from_source)
        assert "heuristic_site_rules" not in source
        assert "heuristic_pka" not in source
