"""Calculate physicochemical descriptors using RDKit."""
from __future__ import annotations

from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.rdchem import Mol


def calculate(mol: Mol) -> dict:
    """Return a dict of physicochemical descriptors for *mol*.

    Parameters
    ----------
    mol:
        A valid RDKit molecule.

    Returns
    -------
    dict with keys: mw, clogp, tpsa, hbd, hba, rotb, hac, formal_charge
    """
    return {
        "mw": round(Descriptors.ExactMolWt(mol), 4),
        "clogp": round(Descriptors.MolLogP(mol), 4),
        "tpsa": round(Descriptors.TPSA(mol), 4),
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "rotb": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "hac": mol.GetNumHeavyAtoms(),
        "formal_charge": sum(atom.GetFormalCharge() for atom in mol.GetAtoms()),
    }
