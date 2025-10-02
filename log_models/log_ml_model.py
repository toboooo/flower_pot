import numpy as np
from rdkit.Chem.rdMolDescriptors import *
from rdkit.Chem import rdFingerprintGenerator

class MLPWrapper():
	def __init__(self, weights, biases, x_mean, x_std):
		self.weights = list(weights[f"arr_{i}"] for i in range(len(weights)))
		self.biases = list(biases[f"arr_{i}"] for i in range(len(biases)))
		self.x_mean = x_mean
		self.x_std = x_std
		self.n_layers = len(weights)

	@staticmethod
	def fingerprint_molecules(mols):
		fingerprints = []
		mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=1024)
		for mol in mols:
			fingerprints.append(mfpgen.GetCountFingerprintAsNumPy(mol))
		return np.array(fingerprints)

	@staticmethod
	def featurize_molecules(mols):
		values = []
		for mol in mols:
			mol_features = []
			crippen_descriptors = CalcCrippenDescriptors(mol)
			mol_features.append(crippen_descriptors[0])
			mol_features.append(crippen_descriptors[1])
			mol_features.append(CalcTPSA(mol))
			mol_features.append(CalcLabuteASA(mol))
			mol_features.append(CalcFractionCSP3(mol))
			mol_features.append(CalcNumAliphaticCarbocycles(mol))
			mol_features.append(CalcNumAliphaticHeterocycles(mol))
			mol_features.append(CalcNumAliphaticRings(mol))
			mol_features.append(CalcNumAromaticCarbocycles(mol))
			mol_features.append(CalcNumAromaticHeterocycles(mol))
			mol_features.append(CalcNumAromaticRings(mol))
			mol_features.append(CalcNumHBA(mol))
			mol_features.append(CalcNumHBD(mol))
			values.append(np.array(mol_features))
		values = np.array(values)
		fingerprints = MLPWrapper.fingerprint_molecules(mols)
		values = np.concatenate((values, fingerprints), axis=1)
		return values, (1, 2, 3)

	def relu(self, z):
		return np.maximum(z, 0, out=z)

	def predict(self, mols):
		if self.x_mean is None or self.x_std is None:
			return None
		x_values, scale_indices = MLPWrapper.featurize_molecules(mols)
		x_values[:,scale_indices] = (x_values[:,scale_indices] - self.x_mean) / self.x_std
		activation = x_values
		for i in range(self.n_layers):
			activation = np.dot(activation, self.weights[i])
			activation += self.biases[i]
			if i != self.n_layers - 1:
				self.relu(activation)
		return activation.ravel()


if __name__ == "__main__":
	import re
	import csv
	import matplotlib.pyplot as plt
	from sklearn.neural_network import MLPRegressor
	from rdkit import Chem

	METALS = ("Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv")

	def filter_training_data(filename, prop_key):
		metal_regexes = tuple(re.compile(f"\[{metal}.*\]") for metal in METALS)
		molecules = []
		property_values = []
		with open(filename, "r") as file:
			reader = csv.DictReader(file)
			for line in reader:
				no_metal = True
				smiles_string = line["smiles"]
				for metal_regex in metal_regexes:
					metal_match = metal_regex.search(smiles_string)
					if metal_match:
						no_metal = False
						break
				if not no_metal:
					continue
				mol = Chem.MolFromSmiles(smiles_string)
				if mol is None:
					continue
				if CalcNumHeavyAtoms(mol) >= 3:
					molecules.append(mol)
					property_values.append(float(line[prop_key]))
		return molecules, property_values

	np.random.seed(5)
	logs_network = MLPRegressor(hidden_layer_sizes=(512,256), activation="relu", alpha=0.001, solver="sgd", early_stopping=True, max_iter=2000)
	sol_mols, solubilities = filter_training_data("solubility.csv", "logs")
	shuffled_indices = np.arange(len(sol_mols))
	np.random.shuffle(shuffled_indices)
	shuffled_sol_mols = []
	for index in shuffled_indices:
		shuffled_sol_mols.append(sol_mols[index])
	logs = np.array(solubilities)
	logs = logs[shuffled_indices]
	train_thresh = int(0.8 * len(sol_mols))
	mols_train, mols_test = shuffled_sol_mols[:train_thresh], shuffled_sol_mols[train_thresh:]
	logs_train, logs_test = logs[:train_thresh], logs[train_thresh:]
	x_train, scale_indices = MLPWrapper.featurize_molecules(mols_train)
	x_mean, x_std = x_train[:,scale_indices].mean(axis=0), x_train[:,scale_indices].std(axis=0)
	x_train[:,scale_indices] = (x_train[:,scale_indices] - x_mean) / x_std
	x_test, scale_indices = MLPWrapper.featurize_molecules(mols_test)
	x_test[:,scale_indices] = (x_test[:,scale_indices] - x_mean) / x_std
	logs_network.fit(x_train, logs_train)
	train_preds = logs_network.predict(x_train)
	preds = logs_network.predict(x_test)
	print("LogS Train MAE = %.5f log units" % np.mean(np.abs(train_preds - logs_train)))
	print("LogS Test MAE = %.5f log units" % np.mean(np.abs(preds - logs_test)))
	plt.plot(train_preds, logs_train, "o", label="Train mols")
	plt.plot(preds, logs_test, "o", label="Test mols")
	plt.xlabel("Predicted LogS")
	plt.ylabel("Actual LogS")
	plt.legend()
	plt.savefig("../results/logs_network_plot.png")
	plt.clf()
	logs = np.array(solubilities)
	features, scale_indices = MLPWrapper.featurize_molecules(sol_mols)
	mean, std = features[:,scale_indices].mean(axis=0), features[:,scale_indices].std(axis=0)
	features[:,scale_indices] = (features[:,scale_indices] - mean) / std
	print("Fitting full neural network for solubility...")
	logs_network.fit(features, logs)
	print("Done! Saving LogS network weights.")
	np.savez("logs_network_weights.npz", *logs_network.coefs_)
	np.savez("logs_network_biases.npz", *logs_network.intercepts_)
	print("LogS feature means: %s" % mean)
	print("LogS feature stds: %s" % std)

	logd_network = MLPRegressor(hidden_layer_sizes=(512,256), activation="relu", alpha=0.001, solver="sgd", early_stopping=True, max_iter=2000)
	lipo_mols, lipophilicities = filter_training_data("lipophilicity.csv", "logd")
	shuffled_indices = np.arange(len(lipo_mols))
	np.random.shuffle(shuffled_indices)
	shuffled_lipo_mols = []
	for index in shuffled_indices:
		shuffled_lipo_mols.append(lipo_mols[index])
	logd = np.array(lipophilicities)
	logd = logd[shuffled_indices]
	train_thresh = int(0.8 * len(lipo_mols))
	mols_train, mols_test = shuffled_lipo_mols[:train_thresh], shuffled_lipo_mols[train_thresh:]
	logd_train, logd_test = logd[:train_thresh], logd[train_thresh:]
	x_train, scale_indices = MLPWrapper.featurize_molecules(mols_train)
	x_mean, x_std = x_train[:,scale_indices].mean(axis=0), x_train[:,scale_indices].std(axis=0)
	x_train[:,scale_indices] = (x_train[:,scale_indices] - x_mean) / x_std
	x_test, scale_indices = MLPWrapper.featurize_molecules(mols_test)
	x_test[:,scale_indices] = (x_test[:,scale_indices] - x_mean) / x_std
	logd_network.fit(x_train, logd_train)
	train_preds = logd_network.predict(x_train)
	preds = logd_network.predict(x_test)
	print("LogD Train MAE = %.5f log units" % np.mean(np.abs(train_preds - logd_train)))
	print("LogD Test MAE = %.5f log units" % np.mean(np.abs(preds - logd_test)))
	plt.plot(train_preds, logd_train, "o", label="Train mols")
	plt.plot(preds, logd_test, "o", label="Test mols")
	plt.xlabel("Predicted LogD")
	plt.ylabel("Actual LogD")
	plt.legend()
	plt.savefig("../results/logd_network_plot.png")
	plt.clf()
	logd = np.array(lipophilicities)
	features, scale_indices = MLPWrapper.featurize_molecules(lipo_mols)
	mean, std = features[:,scale_indices].mean(axis=0), features[:,scale_indices].std(axis=0)
	features[:,scale_indices] = (features[:,scale_indices] - mean) / std
	print("Fitting full neural network for lipophilicity...")
	logd_network.fit(features, logd)
	print("Done! Saving LogD network weights.")
	np.savez("logd_network_weights.npz", *logd_network.coefs_)
	np.savez("logd_network_biases.npz", *logd_network.intercepts_)
	print("LogD feature means: %s" % mean)
	print("LogD feature stds: %s" % std)
