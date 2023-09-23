import rdkit
from rdkit import Chem
from rdkit import AllChem
from rdkit.Chem import Descriptors, rdFreeSASA, Descriptors3D

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

def get_descriptors(mols, sasa=True):
	descriptors = [Descriptors.CalcMolDescriptors(mol) for mol in mols]
	descriptor_array = np.array([[mol_descrs[descr] for descr in mol_descrs] for mol_descrs in descriptors])
	if sasa:
		freesasas = np.array([rdFreeSASA.CalcSASA(mol) for mol in mols])
		descriptor_array = np.hstack((descriptor_array, freesasas.reshape(freesasas.shape[0], 1)))
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
