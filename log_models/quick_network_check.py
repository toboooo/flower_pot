import numpy as np
from rdkit.Chem.rdMolDescriptors import *
from rdkit.Chem import rdFingerprintGenerator

class MLPWrapper():
	def __init__(self, weights, biases, x_mean, x_std, activation="relu"):
		self.weights = list(weights[f"arr_{i}"] for i in range(len(weights)))
		self.biases = list(biases[f"arr_{i}"] for i in range(len(biases)))
		self.x_mean = x_mean
		self.x_std = x_std
		self.n_layers = len(weights)
		if activation == "relu":
			self.activation_func = self.relu
		elif activation == "tanh":
			self.activation_func = self.tanh

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
			mol_features.append(CalcCrippenDescriptors(mol)[0])
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
		return values, (1, 2)

	def relu(self, z):
		return np.maximum(z, 0, out=z)

	def tanh(self, z):
		return np.tanh(z, out=z)

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
				self.activation_func(activation)
		return activation.ravel()


if __name__ == "__main__":
	import sys
	import csv
	import matplotlib.pyplot as plt
	from sklearn.neural_network import MLPRegressor
	from sklearn.model_selection import GridSearchCV
	sys.path.append("..")
	from utils.mol_list import get_mol_list
	mlp_model = MLPRegressor(solver="sgd", hidden_layer_sizes=(512,256), activation="relu", early_stopping=True, max_iter=5000, alpha=0.001, verbose=True)
	np.random.seed(5)
	sol_smiles = []
	solubilities = []
	with open("solubility.csv", "r") as file:
		reader = csv.DictReader(file)
		for line in reader:
			sol_smiles.append(line["smiles"])
			solubilities.append(float(line["logs"]))
	shuffled_indices = np.arange(len(sol_smiles))
	np.random.shuffle(shuffled_indices)
	shuffled_sol_smiles = []
	for index in shuffled_indices:
		shuffled_sol_smiles.append(sol_smiles[index])
	sol_mols, errors = get_mol_list(shuffled_sol_smiles)
	logs = np.array(solubilities)
	logs = logs[shuffled_indices]
	train_thresh = int(0.8 * len(sol_mols))
	mols_train, mols_test = sol_mols[:train_thresh], sol_mols[train_thresh:]
	logs_train, logs_test = logs[:train_thresh], logs[train_thresh:]
	x_train, scale_indices = MLPWrapper.featurize_molecules(mols_train)
	x_mean, x_std = x_train[:,scale_indices].mean(axis=0), x_train[:,scale_indices].std(axis=0)
	x_train[:,scale_indices] = (x_train[:,scale_indices] - x_mean) / x_std
	x_test, scale_indices = MLPWrapper.featurize_molecules(mols_test)
	x_test[:,scale_indices] = (x_test[:,scale_indices] - x_mean) / x_std
	mlp_model.fit(x_train, logs_train)
	train_preds = mlp_model.predict(x_train)
	preds = mlp_model.predict(x_test)
	print(np.min(preds))
	strange_index = np.where(preds == np.min(preds))[0][0]
	print(strange_index)
	print(shuffled_sol_smiles[train_thresh+strange_index])
	print("LogS Train MAE = %.5f log units" % np.mean(np.abs(train_preds - logs_train)))
	print("LogS Test MAE = %.5f log units" % np.mean(np.abs(preds - logs_test)))
	exit()

	mlp_model = MLPRegressor(solver="sgd", hidden_layer_sizes=(1000,), activation="relu", early_stopping=True, max_iter=5000, alpha=0.001)
	lipo_smiles = []
	lipophilicities = []
	with open("lipophilicity.csv", "r") as file:
		reader = csv.DictReader(file)
		for line in reader:
			lipo_smiles.append(line["smiles"])
			lipophilicities.append(float(line["logd"]))
	shuffled_indices = np.arange(len(lipo_smiles))
	np.random.shuffle(shuffled_indices)
	shuffled_lipo_smiles = []
	for index in shuffled_indices:
		shuffled_lipo_smiles.append(lipo_smiles[index])
	lipo_mols, errors = get_mol_list(shuffled_lipo_smiles)
	logd = np.array(lipophilicities)
	logd = logd[shuffled_indices]
	train_thresh = int(0.8 * len(lipo_mols))
	mols_train, mols_test = lipo_mols[:train_thresh], lipo_mols[train_thresh:]
	logd_train, logd_test = logd[:train_thresh], logd[train_thresh:]
	x_train, scale_indices = MLPWrapper.featurize_molecules(mols_train)
	x_mean, x_std = x_train[:,scale_indices].mean(axis=0), x_train[:,scale_indices].std(axis=0)
	x_train[:,scale_indices] = (x_train[:,scale_indices] - x_mean) / x_std
	x_test, scale_indices = MLPWrapper.featurize_molecules(mols_test)
	x_test[:,scale_indices] = (x_test[:,scale_indices] - x_mean) / x_std
	mlp_model.fit(x_train, logd_train)
	train_preds = mlp_model.predict(x_train)
	preds = mlp_model.predict(x_test)
	print("LogD Train MAE = %.5f log units" % np.mean(np.abs(train_preds - logd_train)))
	print("LogD Test MAE = %.5f log units" % np.mean(np.abs(preds - logd_test)))
