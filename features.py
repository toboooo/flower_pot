import numpy as np
import rdkit
from rdkit import Chem
from rdkit import AllChem
from rdkit.Chem import Descriptors, rdFreeSASA
from rdkit.Chem.rdMolDescriptors import CalcMQNs, CalcAUTOCORR2D, BCUT2D
from rdkit.Chem.rdMolDescriptors import CalcAUTOCORR3D, CalcRDF, CalcMORSE, CalcWHIM, CalcGETAWAY
from rdkit.Chem.Descriptors3D import *

def get_mol_list(smiles_strings):
	mols = []
	for smiles in smiles_strings:
		mol = Chem.MolFromSmiles(smiles, sanitize=False)
		if mol is None:
			print("ERROR: %s - Invalid smiles string.")
			continue
		sanitize_result = Chem.SanitizeMol(mol, catchErrors=True)
		if sanitize_result == 0:
			mols.append(mol)
		else:
			print("ERROR: %s - SanitizeMol failed. Flag: %s" % (smiles, str(rdkit.Chem.rdmolops.SanitizeFlags.values[sanitize_result])))
	return mols

def get_descriptors(mols):
	functions = [rdFreeSASA.CalcSASA, CalcMQNs, CalcAUTOCORR2D, BCUT2D]
	descriptors = [Descriptors.CalcMolDescriptors(mol) for mol in mols]
	descriptor_array = np.array([[mol_descrs[descr] for descr in mol_descrs] for mol_descrs in descriptors])
	for function in functions:
		additional_descrs = np.array([function(mol) for mol in mols])
		descriptor_array = np.hstack((descriptor_array, additional_descrs.reshape(additional_descrs.shape[0], 1)))
	descr_names = [name for name in descriptors[0].keys()]
	to_keep = ~np.all(descriptor_array == descriptor_array[0,:], axis=0)
	ipc_index = descr_names.index("Ipc")
	to_keep[ipc_index] = False
	return descriptor_array[:,to_keep]

def get_fingerprints(mols, radius=2, fp_length=2048):
	fingerprints_array = np.zeros((len(mols), fp_length))
	fp_array = np.zeros(fp_length)
	for m, mol in enumerate(mols):
		fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=fp_length)
		Chem.DataStructs.ConvertToNumpyArray(fingerprint, fp_array)
		fingerprints_array[m] = fp_array
	return fingerprints_array

def get_3d_descriptors(mols, mmff_optimise=False):
	descriptors = []
	functions = [Asphericity, Eccentricity, InertialShapeFactor, NPR1, NPR2, PMI1, PMI2, PMI3, RadiusOfGyration, SpherocityIndex, CalcAUTOCORR3D, CalcRDF, CalcMORSE, CalcWHIM, CalcGETAWAY]
	for mol in mols:
		mol3d = Chem.AddHs(mol)
		AllChem.EmbedMolecule(mol3d)
		if mmff_optimise:
			AllChem.MMFFOptimizeMolecule(mol3d)
		descriptors.append([function(mol3d) for function in functions])
	descriptor_array = np.array(descriptors)
	to_keep = ~np.all(descriptor_array == descriptor_array[0,:], axis=0)
	return descriptor_array[:,to_keep]
