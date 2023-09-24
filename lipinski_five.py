from rdkit import Chem
from rdkit.Chem import Descriptors

def lipinski_five_test(mols):
	results = []
	for mol in mols:
		MW = Descriptors.MolWt(mol)
		HBA = Descriptors.NOCount(mol)
		HBD = Descriptors.NHOHCount(mol)
		LogP = Descriptors.MolLogP(mol)
		conditions = [MW <= 500, HBA <= 10, HBD <= 5, LogP <= 5]
		results.append([conditions.count(True) >= 3] + conditions)
	return results
