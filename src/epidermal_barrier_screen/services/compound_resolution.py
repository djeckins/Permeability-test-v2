"""Resolve names/synonyms/InChIKey to canonical SMILES (database-first)."""
from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field

from rdkit import Chem

from epidermal_barrier_screen.adapters import pubchem_adapter

log = logging.getLogger(__name__)


@dataclass
class ResolvedCompound:
    input_text: str
    input_type: str
    matched_name: str | None = None
    canonical_smiles: str | None = None
    inchikey: str | None = None
    molecular_formula: str | None = None
    pubchem_cid: int | None = None
    source_name: str | None = None
    source_identifier: str | None = None
    resolution_confidence: str = "none"
    resolution_notes: str | None = None
    synonyms: list[str] = field(default_factory=list)


_INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
_SMILES_STRUCTURAL_RE = re.compile(r"[=\(\)#\[\]\\/@\d]")
_ONLY_ALPHA_RE = re.compile(r"^[a-zA-Z]+$")
_SMILES_CHARS = set("CNOPSFIBrcnopsfi[]()=#+-\\/@.%0123456789lae")


def detect_input_type(text: str) -> str:
    text = text.strip()
    if not text:
        return "name"
    if _INCHIKEY_RE.match(text):
        return "inchikey"
    if " " in text:
        return "name"
    if _ONLY_ALPHA_RE.match(text) and len(text) > 4:
        return "name"
    if _SMILES_STRUCTURAL_RE.search(text):
        return "smiles" if Chem.MolFromSmiles(text) else "name"
    if all(ch in _SMILES_CHARS for ch in text):
        return "smiles" if Chem.MolFromSmiles(text) else "name"
    return "name"


def _name_variants(name: str) -> list[str]:
    stripped = name.strip()
    lowered = stripped.lower()
    compact = re.sub(r"[^a-z0-9]+", " ", lowered).strip()
    variants = [stripped]
    for candidate in (lowered, compact, compact.replace(" ", "")):
        if candidate and candidate not in variants:
            variants.append(candidate)
    return variants


def resolve_compound(text: str) -> ResolvedCompound:
    text = text.strip()
    input_type = detect_input_type(text)
    result = ResolvedCompound(input_text=text, input_type=input_type)

    if input_type == "smiles":
        mol = Chem.MolFromSmiles(text)
        if mol is None:
            result.resolution_notes = "invalid_smiles"
            return result
        result.canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
        try:
            from rdkit.Chem import inchi as rdkit_inchi

            result.inchikey = rdkit_inchi.MolToInchiKey(mol)
        except Exception as exc:
            log.warning("Failed to generate InChIKey from SMILES: %s", exc)
        result.source_name = "rdkit"
        result.source_identifier = "rdkit_canonicalization"
        result.resolution_confidence = "high"
        return result

    if input_type == "inchikey":
        result.inchikey = text
        identity = pubchem_adapter.resolve_inchikey(text)
        if identity and identity.canonical_smiles:
            result.canonical_smiles = identity.canonical_smiles
            result.pubchem_cid = identity.cid
            result.molecular_formula = identity.molecular_formula
            result.matched_name = identity.iupac_name
            result.source_name = "pubchem"
            result.source_identifier = str(identity.cid)
            result.resolution_confidence = "high"
            return result
        result.resolution_notes = "inchikey_not_found_in_pubchem"
        return result

    # name / synonym path
    for variant in _name_variants(text):
        identity = pubchem_adapter.resolve_name(variant)
        if not identity or not identity.canonical_smiles:
            continue

        result.canonical_smiles = identity.canonical_smiles
        result.inchikey = identity.inchikey
        result.pubchem_cid = identity.cid
        result.molecular_formula = identity.molecular_formula
        result.matched_name = identity.iupac_name or variant
        result.synonyms = identity.synonyms
        result.source_name = "pubchem"
        result.source_identifier = str(identity.cid) if identity.cid else None

        synonyms_lower = {s.lower() for s in identity.synonyms}
        if text.lower() in synonyms_lower or text.lower() == (result.matched_name or "").lower():
            result.resolution_confidence = "high"
            result.resolution_notes = "exact_or_synonym_match"
        elif variant != text:
            result.resolution_confidence = "moderate"
            result.resolution_notes = f"resolved_via_normalized_variant:{variant}"
        else:
            result.resolution_confidence = "moderate"
            result.resolution_notes = "resolved_name_without_explicit_synonym_match"
        return result

    result.resolution_notes = "name_not_found_in_pubchem"
    return result
