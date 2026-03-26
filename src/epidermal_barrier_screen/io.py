"""Parse various input formats into a list of molecule records.

Supports SMILES, compound names (auto-resolved via PubChem), SDF, and ZIP.
"""
from __future__ import annotations

import io
import logging
import zipfile
from typing import Any

from rdkit import Chem

from epidermal_barrier_screen.services.compound_resolution import (
    detect_input_type,
    resolve_compound,
)

log = logging.getLogger(__name__)

_PKA_KEYS = ("pKa", "PKA", "predicted_pKa", "input_pka")
_PKA_ACID_KEYS = (
    "pKa_strongest_acidic",
    "PKA_STRONGEST_ACIDIC",
    "acidic_pKa",
    "input_pka_acidic",
)
_PKA_BASE_KEYS = (
    "pKa_strongest_basic",
    "PKA_STRONGEST_BASIC",
    "basic_pKa",
    "input_pka_basic",
)
_LOGD_KEYS = ("LogD", "LOGD", "logD_7_4", "input_logd_7_4")
_NAME_KEYS = ("_Name", "Name", "TITLE", "ID", "Compound")


def _try_float(value: str | None) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (ValueError, TypeError):
        return None


def _sdf_prop(mol: Chem.Mol, keys: tuple[str, ...]) -> str | None:
    for key in keys:
        if mol.HasProp(key):
            return mol.GetProp(key)
    return None


def _record_from_mol(mol: Chem.Mol | None, input_smiles: str | None = None) -> dict[str, Any]:
    if mol is None:
        return {
            "name": None,
            "input_smiles": input_smiles,
            "canonical_smiles": None,
            "mol": None,
            "parse_status": "invalid",
            "input_pka": None,
            "input_pka_acidic": None,
            "input_pka_basic": None,
            "input_logd_7_4": None,
            # Resolution fields
            "input_name": None,
            "input_type": "smiles",
            "matched_name": None,
            "resolution_source": None,
            "resolution_confidence": "none",
            "resolution_notes": "invalid_smiles",
        }

    name = _sdf_prop(mol, _NAME_KEYS) or input_smiles
    canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
    return {
        "name": name,
        "input_smiles": input_smiles or canonical_smiles,
        "canonical_smiles": canonical_smiles,
        "mol": mol,
        "parse_status": "ok",
        "input_pka": _try_float(_sdf_prop(mol, _PKA_KEYS)),
        "input_pka_acidic": _try_float(_sdf_prop(mol, _PKA_ACID_KEYS)),
        "input_pka_basic": _try_float(_sdf_prop(mol, _PKA_BASE_KEYS)),
        "input_logd_7_4": _try_float(_sdf_prop(mol, _LOGD_KEYS)),
        # Resolution fields
        "input_name": None,
        "input_type": "smiles",
        "matched_name": None,
        "resolution_source": "rdkit",
        "resolution_confidence": "high",
        "resolution_notes": None,
    }


def _parse_smiles_or_name(text: str) -> dict[str, Any]:
    """Parse input text, auto-detecting whether it's SMILES or a compound name."""
    text = text.strip()
    if not text:
        return _record_from_mol(None, input_smiles=None)

    input_type = detect_input_type(text)

    if input_type == "smiles":
        mol = Chem.MolFromSmiles(text)
        record = _record_from_mol(mol, input_smiles=text)
        record["input_type"] = "smiles"
        return record

    # It's a name or InChIKey — resolve via PubChem
    resolved = resolve_compound(text)

    if resolved.canonical_smiles:
        mol = Chem.MolFromSmiles(resolved.canonical_smiles)
    else:
        mol = None

    if mol is None:
        record = _record_from_mol(None, input_smiles=text)
        record["input_name"] = text
        record["input_type"] = input_type
        record["resolution_notes"] = resolved.resolution_notes or "resolution_failed"
        return record

    canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
    record = {
        "name": resolved.matched_name or text,
        "input_smiles": text,
        "canonical_smiles": canonical_smiles,
        "mol": mol,
        "parse_status": "ok",
        "input_pka": None,
        "input_pka_acidic": None,
        "input_pka_basic": None,
        "input_logd_7_4": None,
        # Resolution fields
        "input_name": text,
        "input_type": input_type,
        "matched_name": resolved.matched_name,
        "resolution_source": resolved.source_name,
        "resolution_confidence": resolved.resolution_confidence,
        "resolution_notes": resolved.resolution_notes,
    }

    # Store PubChem CID if we got it during resolution
    if resolved.pubchem_cid:
        record["resolved_pubchem_cid"] = resolved.pubchem_cid
    if resolved.inchikey:
        record["resolved_inchikey"] = resolved.inchikey
    if resolved.molecular_formula:
        record["resolved_molecular_formula"] = resolved.molecular_formula

    return record


def _parse_smiles(smiles: str) -> dict[str, Any]:
    """Legacy: parse a known SMILES string."""
    smiles = smiles.strip()
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    return _record_from_mol(mol, input_smiles=smiles or None)


def _records_from_sdf_supplier(supplier: Chem.ForwardSDMolSupplier) -> list[dict[str, Any]]:
    return [_record_from_mol(mol) for mol in supplier]


def parse_input(
    mode: str,
    payload: str | bytes,
    *,
    filename: str | None = None,
) -> list[dict[str, Any]]:
    if mode == "smiles":
        return [_parse_smiles_or_name(str(payload))]

    if mode == "smiles_list":
        records = []
        for raw_line in str(payload).splitlines():
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            # Strategy: first try to parse the first token as SMILES.
            # If the first token looks like a name, treat the WHOLE line as a
            # name (compound names may contain spaces, e.g. "caffeic acid").
            name_hint: str | None = None

            if "," in line and not line.startswith("["):
                # Comma-separated: first field is SMILES or name, rest is hint
                first, rest = line.split(",", 1)
                text = first.strip()
                name_hint = rest.strip() or None
            else:
                parts = line.split(None, 1)
                first_token = parts[0]
                first_type = detect_input_type(first_token)

                if first_type == "smiles":
                    # First token is valid SMILES; rest is name hint
                    text = first_token
                    name_hint = parts[1].strip() if len(parts) > 1 else None
                else:
                    # First token is not SMILES — treat the entire line as a name
                    text = line
                    name_hint = None

            record = _parse_smiles_or_name(text)

            if name_hint and (record["name"] == text or record["name"] is None):
                record["name"] = name_hint
            records.append(record)
        return records

    if mode == "sdf":
        if not isinstance(payload, (bytes, bytearray)):
            raise TypeError("SDF mode expects bytes payload")
        supplier = Chem.ForwardSDMolSupplier(io.BytesIO(payload), removeHs=False)
        return _records_from_sdf_supplier(supplier)

    if mode == "sdf_zip":
        if not isinstance(payload, (bytes, bytearray)):
            raise TypeError("sdf_zip mode expects bytes payload")
        records: list[dict[str, Any]] = []
        with zipfile.ZipFile(io.BytesIO(payload)) as zf:
            sdf_names = [n for n in zf.namelist() if n.lower().endswith(".sdf")]
            if not sdf_names:
                raise ValueError("ZIP archive contains no .sdf files")
            for sdf_name in sorted(sdf_names):
                sdf_bytes = zf.read(sdf_name)
                supplier = Chem.ForwardSDMolSupplier(io.BytesIO(sdf_bytes), removeHs=False)
                records.extend(_records_from_sdf_supplier(supplier))
        return records

    raise ValueError(f"Unknown input mode: {mode!r}")
