import multiprocessing
import csv
import numpy as np
from rdkit.Chem import Descriptors
from features import get_mol_list

def get_start_vectors(grid_size, x_bounds, y_bounds, d_factors):
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
	dc = (c * d_factors).flatten()
	xy_indices = np.repeat(np.arange(xy.shape[0]), d_factors.shape[0])
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
	return smiles, np.array(labels)

def normalise_features(features):
	stddev = np.std(features, axis=0)
	return features / stddev, stddev

def ellipse_classify(mol_features, ellipse_vectors):
	distances = np.linalg.norm(ellipse_vectors[:,np.newaxis,(0,1)] - mol_features[np.newaxis,:,:], axis=2)
	distances += np.linalg.norm(ellipse_vectors[:,np.newaxis,(2,3)] - mol_features[np.newaxis,:,:], axis=2)
	return (distances < ellipse_vectors[:,4].reshape(ellipse_vectors.shape[0], 1)).astype(int)

def calc_mcc(classifications, labels):
	tp = np.sum(np.logical_and(classifications == 1, labels == 1), axis=1)
	tn = np.sum(np.logical_and(classifications == 0, labels == 0), axis=1)
	fp = np.sum(np.logical_and(classifications == 1, labels == 0), axis=1)
	fn = np.sum(np.logical_and(classifications == 0, labels == 1), axis=1)
	mcc = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
	mcc = np.where(mcc == 0.0, 1.0, mcc)
	mcc = ((tp * tn) - (fp * fn)) / mcc
	return mcc

def calc_area(ellipse_vectors):
	c2 = np.sum(np.square(ellipse_vectors[:,(0,1)] - ellipse_vectors[:,(2,3)]), axis=1) / 4
	a = ellipse_vectors[:,4] / 2
	b = np.sqrt(np.round(a**2 - c2, decimals=14))
	return np.pi * a * b

def optimise_ellipse(mol_features, ellipse_vectors, labels, area_penalty=0.05, change_size=0.5, random_mult=0.025):
	classifications = ellipse_classify(mol_features, ellipse_vectors)
	mcc = calc_mcc(classifications, labels)
	area = calc_area(ellipse_vectors)
	best_score = mcc - area_penalty * area
	best_mcc = np.max(mcc)
	best_vector = ellipse_vectors[np.argmax(mcc)]
	for i in range(100000):
		change = np.random.uniform(-change_size, change_size, size=ellipse_vectors.shape)
		tentative_vectors = ellipse_vectors + change
		classifications = ellipse_classify(mol_features, tentative_vectors)
		mcc = calc_mcc(classifications, labels)
		area = calc_area(tentative_vectors)
		score = mcc - area_penalty * area
		tentative_randoms = random_mult * np.random.uniform(0.0, 1.0, size=score.shape)
		accept_indices = np.logical_and((score + tentative_randoms) > best_score, tentative_vectors[:,4] > 0.0)
		ellipse_vectors[accept_indices] = tentative_vectors[accept_indices]
		best_score[accept_indices] = score[accept_indices]
		if np.max(mcc) > best_mcc:
			best_mcc = np.max(mcc)
			best_vector = tentative_vectors[np.argmax(mcc)]

		print("iter: %d, best score = %.7f, best MCC = %.7f, vector: %s" % (i, np.max(best_score), best_mcc, str(best_vector)))

	return best_vector, best_mcc

def run_egg_optimisation(mol_features, ellipse_vectors, labels, stddev, seed, type, area_penalty=0.05, change_size=0.5, random_mult=0.025):
	np.random.seed(seed)
	best_vector, best_mcc = optimise_ellipse(mol_features, ellipse_vectors, labels, area_penalty=area_penalty, change_size=change_size, random_mult=random_mult)
	denormed_vector = best_vector.copy()
	denormed_vector[4] /= np.linalg.norm(denormed_vector[(0,1)] - denormed_vector[(2,3)])
	denormed_vector[(0,2)] *= stddev[0]
	denormed_vector[(1,3)] *= stddev[1]
	print(type, "MCC:", best_mcc, best_vector, denormed_vector)


if __name__ == "__main__":
	hia_smiles, hia_labels = read_from_file("hia_mols.csv")
	bbb_smiles, bbb_labels = read_from_file("bbb_mols.csv")
	hia_mols = get_mol_list(hia_smiles)
	bbb_mols = get_mol_list(bbb_smiles)
	hia_features = get_egg_features(hia_mols)
	bbb_features = get_egg_features(bbb_mols)
	hia_features, hia_stddev = normalise_features(hia_features)
	bbb_features, bbb_stddev = normalise_features(bbb_features)
	grid_size = 9
	hia_x_bounds = np.array([0, 175]) / hia_stddev[0]
	hia_y_bounds = np.array([-1, 8]) / hia_stddev[1]
	bbb_x_bounds = np.array([0, 120]) / bbb_stddev[0]
	bbb_y_bounds = np.array([-1, 8]) / bbb_stddev[1]
	d_factors = np.array([1.0, 1.5, 2.0, 3.0])
	hia_start_vectors = get_start_vectors(grid_size, hia_x_bounds, hia_y_bounds, d_factors)
	bbb_start_vectors = get_start_vectors(grid_size, bbb_x_bounds, bbb_y_bounds, d_factors)
	for i in range(10):
		p = multiprocessing.Process(target=run_egg_optimisation, args=(hia_features, hia_start_vectors, hia_labels, hia_stddev, i, "HIA"))
		p.start()
	for i in range(10):
		p = multiprocessing.Process(target=run_egg_optimisation, args=(bbb_features, bbb_start_vectors, bbb_labels, bbb_stddev, i, "BBB"))
		p.start()
