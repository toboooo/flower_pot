# Ellipse perimeter approximation: https://www.mathsisfun.com/geometry/ellipse-perimeter.html
import csv
import numpy as np
from sklearn.metrics import matthews_corrcoef
from rdkit.Chem import Descriptors
from features import get_mol_list

def get_start_vectors(grid_size=9, x_bounds=(0,175), y_bounds=(-1,8), d_factors=(1,1.5,2,3)):
	x = np.linspace(x_bounds[0], x_bounds[1], grid_size)
	x2, x1 = np.meshgrid(x, x)
	x1, x2 = x1[np.triu_indices_from(x1)].flatten(), x2[np.triu_indices_from(x2)].flatten()
	y = np.linspace(y_bounds[0], y_bounds[1], grid_size)
	y2, y1 = np.meshgrid(y, y)
	y1, y2 = y1[np.tril_indices_from(y1)].flatten(), y2[np.tril_indices_from(y2)].flatten()
	x_indices, y_indices = np.meshgrid(np.arange(x1.shape[0]), np.arange(x1.shape[0]))
	x_indices, y_indices = x_indices.flatten(), y_indices.flatten()
	x1, y1, x2, y2 = x1[x_indices], y1[y_indices], x2[x_indices], y2[y_indices]
	xy = np.column_stack((x1, y1, x2, y2))
	c = np.linalg.norm(xy[:,(0,1)] - xy[:,(2,3)], axis=1).reshape(xy.shape[0], 1)
	d = np.array(d_factors)
	dc = (c * d).flatten()
	xy_indices = np.repeat(np.arange(xy.shape[0]), d.shape[0])
	return np.column_stack((xy[xy_indices], dc))

def get_egg_features(mols):
	features = []
	for mol in mols:
		features.append([Descriptors.TPSA(mol), Descriptors.MolLogP(mol)])
	return np.array(features)

def read_from_file(filename):
	file = open(filename, "r")
	csv_reader = csv.DictReader(file)
	smiles = []
	labels = []
	for line in csv_reader:
		smiles.append(line["SMILES"])
		if "hia" in filename:
			prop = "HIA"
		elif "bbb" in filename:
			prop = "BBB"
		if line[prop] == "-":
			labels.append(0)
		elif line[prop] == "+":
			labels.append(1)
	return smiles, labels

def normalise_features(features):
	stddev = np.std(features, axis=0)
	return features / stddev, stddev

def unnormalise_features(features, stddev):
	return features * stddev
