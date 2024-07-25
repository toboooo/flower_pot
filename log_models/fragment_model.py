"""Contains the implementation of the Crippen-style atomic contribution model to
predict LogD and LogS, as well as code to optimise it for the lipophilicity and
solubility datasets."""
import os
import numpy as np
from rdkit import Chem

class FragmentLogModel():
	"""Implements the atomic contribution model to predict LogS or LogD."""
	def __init__(self, fragments_file, coefs=None):
		"""
		Class constructor for the atomic contribution model.
		Args:
			fragments_file: str, file name ('Crippen.txt') that contains the
				SMARTS patterns for each atom type in the model.
			coefs: optional, numpy array that contains the atomic contribution
				coeficients.
		"""
		if os.path.exists(fragments_file):
			with open(fragments_file, "r") as file:
				self.patterns = []
				for line in file:
					if line.startswith("#"):
						continue
					smarts = line.split("\t")[1]
					if smarts != "":
						pattern = Chem.MolFromSmarts(smarts)
						self.patterns.append(pattern)
		else:
			print("ERROR: Could not find file %s." % fragments_file)
			self.patterns = None
		self.coefs = coefs

	def load_coef_file(self, coefs_filename):
		"""
		Allows coefficients stored in a binary numpy .npy file to be loaded into
		the model.
		Args:
			coefs_filename: str, the file name for the binary file containing
				atomic contribution coefficients.
		"""
		if self.patterns is not None:
			coefs = np.load(coefs_filename)
			if coefs.shape[0] == len(self.patterns):
				self.coefs = coefs
			else:
				print("ERROR: Different number of fragments and coefficients.")
		else:
			print("ERROR: Fragment SMARTS have not been initialised.")

	def get_fragment_counts(self, mols):
		"""
		For a series of molecules, determines the counts of each type of atom
		known to the model in each molecule.
		Args:
			mols: list of Mol objects, a list of molecules for which predictions
				are to be made.
		Returns:
			x: 2d numpy array where each row corresponds to a molecule in mols
				and each column corresponds to the count of each atom type.
		"""
		if self.patterns is not None:
			x = np.zeros((len(mols), len(self.patterns)))
			for i, mol in enumerate(mols):
				for j, pattern in enumerate(self.patterns):
					matches = mol.GetSubstructMatches(pattern)
					x[i,j] = len(matches)
			return x
		else:
			print("ERROR: Fragment SMARTS have not been initialised.")
			return None

	def fit(self, mols, y):
		"""
		Optimises the atomic contribution coeffients using the least-squares
		solver from numpy to obtain the best fit to the experimental data.
		Args:
			mols: list of Mol objects, list of molecules from which the model is
				to be trained.
			y: numpy array of experimental target values for each molecule.
		"""
		x = self.get_fragment_counts(mols)
		if x is not None:
			lstsq_result = np.linalg.lstsq(x, y, rcond=None)
			self.coefs = lstsq_result[0]
		else:
			print("ERROR: Failed to get fragment counts.")

	def predict(self, mols):
		"""
		Predicts LogS or LogD values for the input molecules using the fitted
		model coefficients.
		Args:
			mols: list of Mol objects, a list of molecules for which predictions
				are to be made.
		Returns:
			Predicted LogS or LogD values for each molecule.
		"""
		if self.coefs is not None:
			x = self.get_fragment_counts(mols)
			if x is not None:
				return np.dot(x, self.coefs)
			else:
				print("ERROR: Failed to get fragment counts.")
		else:
			print("ERROR: Fragment model has not yet been fit.")


if __name__ == "__main__":
	import sys
	import csv
	import matplotlib.pyplot as plt
	sys.path.append("..")
	from mol_list import get_mol_list

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
	logs_model = FragmentLogModel("Crippen.txt")
	train_thresh = int(0.8 * len(sol_mols))
	mols_train, mols_test = sol_mols[:train_thresh], sol_mols[train_thresh:]
	logs_train, logs_test = logs[:train_thresh], logs[train_thresh:]
	logs_model.fit(mols_train, logs_train)
	train_predictions = logs_model.predict(mols_train)
	predictions = logs_model.predict(mols_test)
	print("Test MAE = %.5f" % np.mean(np.abs(logs_test - predictions)))
	plt.plot(train_predictions, logs_train, "o", label="Train mols")
	plt.plot(predictions, logs_test, "o", label="Test mols")
	plt.title("LogS")
	plt.xlabel("Predicted LogS")
	plt.ylabel("Actual LogS")
	plt.legend()
	plt.show()
	logs_model.fit(sol_mols, logs)
	with open("sol_coefs.npy", "wb") as file:
		np.save(file, logs_model.coefs)

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
	logd_model = FragmentLogModel("Crippen.txt")
	train_thresh = int(0.8 * len(lipo_mols))
	mols_train, mols_test = lipo_mols[:train_thresh], lipo_mols[train_thresh:]
	logd_train, logd_test = logd[:train_thresh], logd[train_thresh:]
	logd_model.fit(mols_train, logd_train)
	train_predictions = logd_model.predict(mols_train)
	predictions = logd_model.predict(mols_test)
	print("Test MAE = %.5f" % np.mean(np.abs(logd_test - predictions)))
	plt.plot(train_predictions, logd_train, "o", label="Train mols")
	plt.plot(predictions, logd_test, "o", label="Test mols")
	plt.title("LogD")
	plt.xlabel("Predicted LogD")
	plt.ylabel("Actual LogD")
	plt.legend()
	plt.show()
	logd_model.fit(lipo_mols, logd)
	with open("lipo_coefs.npy", "wb") as file:
		np.save(file, logd_model.coefs)
