"""Screen molecule records against epidermal barrier passage criteria."""
from __future__ import annotations

from typing import Any

import pandas as pd

from epidermal_barrier_screen.descriptors import calculate
from epidermal_barrier_screen.ionization import DEFAULT_PH, analyze_ionization


# ---------------------------------------------------------------------------
# Per-criterion status functions
# Each returns one of: "optimal", "suboptimal", "poor"
# ---------------------------------------------------------------------------


def _mw_status(v: float) -> str:
    if v < 300:
        return "optimal"
    if v <= 500:
        return "suboptimal"
    return "poor"


def _logd_status(v: float | None) -> str:
    if v is None:
        return "poor"
    if 1.0 <= v <= 3.0:
        return "optimal"
    if 0.5 <= v < 1.0 or 3.0 < v <= 5.0:
        return "suboptimal"
    return "poor"


def _tpsa_status(v: float) -> str:
    if v < 60:
        return "optimal"
    if v <= 130:
        return "suboptimal"
    return "poor"


def _hbd_status(v: int) -> str:
    if 0 <= v <= 3:
        return "optimal"
    if 4 <= v <= 5:
        return "suboptimal"
    return "poor"


def _hba_status(v: int) -> str:
    if 2 <= v <= 8:
        return "optimal"
    if 8 < v <= 10:
        return "suboptimal"
    return "poor"


def _rotb_status(v: int) -> str:
    if v < 10:
        return "optimal"
    if v <= 15:
        return "suboptimal"
    return "poor"


def _hac_status(v: int) -> str:
    if v < 30:
        return "optimal"
    if v <= 50:
        return "suboptimal"
    return "poor"


def _charge_status_from_expected(expected_net_charge: float | None, fallback_formal_charge: int) -> str:
    if expected_net_charge is None:
        v = fallback_formal_charge
    else:
        abs_q = abs(expected_net_charge)
        if abs_q < 0.2:
            v = 0
        elif abs_q <= 1.2:
            v = 1
        else:
            v = 2
    if v == 0:
        return "optimal"
    if v == 1:
        return "suboptimal"
    return "poor"


def _ionization_status(fraction_unionized: float | None) -> str:
    if fraction_unionized is None:
        return "poor"
    if fraction_unionized >= 0.8:
        return "optimal"
    if fraction_unionized >= 0.5:
        return "suboptimal"
    return "poor"


# ---------------------------------------------------------------------------
# Overall result
# ---------------------------------------------------------------------------


def _final_result(statuses: list[str]) -> str:
    """Derive overall PASS / BORDERLINE / FAIL from per-criterion statuses.

    PASS       – at most 1 non-optimal criterion (suboptimal or poor)
    BORDERLINE – more than 1 non-optimal criterion and at most 2 poor criteria
    FAIL       – 3 or more poor criteria
    """
    n_poor = statuses.count("poor")
    n_suboptimal = statuses.count("suboptimal")
    n_non_optimal = n_poor + n_suboptimal

    if n_non_optimal <= 1:
        return "PASS"
    if n_poor <= 2:
        return "BORDERLINE"
    return "FAIL"


_STATUS_COLUMNS = [
    "mw_status",
    "logd_status",
    "tpsa_status",
    "hbd_status",
    "hba_status",
    "rotb_status",
    "hac_status",
    "formal_charge_status",
    "ionization_status",
]


def screen_records(records: list[dict[str, Any]], ph: float = DEFAULT_PH) -> pd.DataFrame:
    rows = []
    for rec in records:
        row: dict[str, Any] = {
            "name": rec.get("name"),
            "input_smiles": rec.get("input_smiles"),
            "canonical_smiles": rec.get("canonical_smiles"),
            "parse_status": rec.get("parse_status", "invalid"),
        }

        if rec.get("parse_status") != "ok" or rec.get("mol") is None:
            row["final_result"] = "invalid_input"
            rows.append(row)
            continue

        mol = rec["mol"]
        desc = calculate(mol)
        row.update(desc)

        ion = analyze_ionization(
            mol=mol,
            canonical_smiles=rec.get("canonical_smiles") or "",
            clogp=desc["clogp"],
            ph=float(ph),
            name=rec.get("name"),
            input_pka=rec.get("input_pka"),
            input_pka_acidic=rec.get("input_pka_acidic"),
            input_pka_basic=rec.get("input_pka_basic"),
            input_logd_7_4=rec.get("input_logd_7_4"),
        )
        row.update(ion)
        row["ph_used"] = round(float(ph), 2)

        row["mw_status"] = _mw_status(desc["mw"])
        row["logd_status"] = _logd_status(row.get("logd"))
        row["tpsa_status"] = _tpsa_status(desc["tpsa"])
        row["hbd_status"] = _hbd_status(desc["hbd"])
        row["hba_status"] = _hba_status(desc["hba"])
        row["rotb_status"] = _rotb_status(desc["rotb"])
        row["hac_status"] = _hac_status(desc["hac"])
        row["formal_charge_status"] = _charge_status_from_expected(
            row.get("expected_net_charge"), desc["formal_charge"]
        )
        row["ionization_status"] = _ionization_status(row.get("fraction_unionized"))

        statuses = [row[c] for c in _STATUS_COLUMNS]
        row["final_result"] = _final_result(statuses)
        rows.append(row)

    col_order = [
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
        "lookup_match_name",
        "pka_note",
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
        "formal_charge",
        "mw_status",
        "logd_status",
        "tpsa_status",
        "hbd_status",
        "hba_status",
        "rotb_status",
        "hac_status",
        "formal_charge_status",
        "ionization_status",
        "final_result",
        "input_smiles",
        "canonical_smiles",
    ]

    df = pd.DataFrame(rows)
    for col in col_order:
        if col not in df.columns:
            df[col] = None
    return df[col_order]
