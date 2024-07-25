"""Provides a utility function that converts a list of SMILES strings to RDKit
Mol objects."""
import numpy as np
import rdkit
from rdkit import Chem

def get_mol_list(smiles_strings, names=None):
	"""
	Takes in a series of SMILES strings and converts them into RDKit Mol
	objects, also providing error messages for any invalid smiles strings. If
	also provided with a series of titles to be associated with the molecules,
	then the titles for all valid SMILES string will also be returned.
	Args:
		smiles_string: list of str, list of SMILES strings.
		names: optional, list of str, any titles that are to be associated with
			the valid molecules.
	Returns:
		mols: list of Mol objects, a list of molecules from all valid SMILES
			strings.
		err: list of str, a list of error messages associated with each SMILES
			string.
		valid_names: list of str, a list of titles associated with all valid
			molecules.
	"""
	mols = []
	errs = []
	valid_names = []
	for i, smiles in enumerate(smiles_strings):
		mol = Chem.MolFromSmiles(smiles, sanitize=False)
		if mol is None:
			errs.append("ERROR: %s - Invalid smiles string." % smiles)
			continue
		sanitize_result = Chem.SanitizeMol(mol, catchErrors=True)
		if sanitize_result == 0:
			mols.append(mol)
			errs.append(None)
			if names is not None:
				valid_names.append(names[i])
		else:
			errs.append("ERROR: %s - SanitizeMol failed. Flag: %s" % (smiles, str(rdkit.Chem.rdmolops.SanitizeFlags.values[sanitize_result])))
	if names is None:
		return mols, errs
	else:
		return mols, errs, valid_names
