"""Performs Monte Carlo optimisation of the HARDBOILED-EGG model starting with
the ellipse parameters from the original BOILED-EGG publication
(doi.org/10.1002/cmdc.201600182) using the TPSA and LogP values from RDKit."""
import sys
import csv
import numpy as np
from rdkit.Chem import Descriptors
sys.path.append("..")
from mol_list import get_mol_list

def get_start_vectors(major_length, minor_length, angle, centre):
	"""
	Generates an array containing 10 repeats of the ellipse parameter vector
	that will be optimised for the HARDBOILED-EGG model with RDKit features.
	Specifically, the ellipse parameters are the x- and y-coordinates of the two
	foci and the length across the major axis of the ellipse.
	Args:
		major_length: float, the length across the major axis.
		minor_length: float, the length across the minor axis.
		angle: float, the tilt angle of the major axis above the x-axis in
			radians.
		centre: tuple of floats, the coordinates of the centre of the ellipse.
	Returns:
		2d array containing 10 repeats of the ellipse parameter vectors that are
		to be optimised.
	"""
	a = major_length / 2
	b = minor_length / 2
	c = np.sqrt(a**2 - b**2)
	x1 = centre[0] - c * np.cos(angle)
	x2 = centre[0] + c * np.cos(angle)
	y1 = centre[1] - c * np.sin(angle)
	y2 = centre[1] + c * np.sin(angle)
	ellipse_vector = np.array([[x1, y1, x2, y2, major_length]])
	return np.repeat(ellipse_vector, 10, axis=0)

def get_egg_features(mols):
	"""
	Calculates TPSA and LogP values for all provided molecules to be used as
	input features for the HARDBOILED-EGG model.
	Args:
		mols: list of Mol objects, a list of molecules.
	Returns:
		2d numpy array containing TSPA values for all molecules in the first
		column and LogP values in the second column.
	"""
	features = []
	for mol in mols:
		features.append([Descriptors.TPSA(mol), Descriptors.MolLogP(mol)])
	return np.array(features)

def read_from_file(filename):
	"""
	Reads through HIA and BBB data set files and returns all SMILES strings and
	permeation classification label.
	Args:
		filename: str, name of HIA or BBB data set file.
	Returns:
		smiles: list of str, list of all molecule SMILES strings in the file.
		labels: np.array, permeation classification of each molecule (0
			indicates non-permeating, 1 indicates permeating).
	"""
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

def ellipse_classify(mol_features, ellipse_vectors):
	"""
	Determines whether a set of data points (TPSA and LogP values) lies inside
	an ellipse described by its ellipse vector.
	Args:
		mol_features: 2d numpy array containing TSPA values for all molecules in
			the first column and LogP values in the second column.
		ellipse_vectors: numpy array containing the ellipse parameters.
			Specifically, the ellipse parameters are the x- and y-coordinates of
			the two foci and the length across the major axis of the ellipse.
	Returns:
		A numpy array in which 0 elements correspond to molecules that lie
		outside the ellipse described by the ellipse vectors and 1 for the
		molecules that are inside the ellipse.
	"""
	distances = np.linalg.norm(ellipse_vectors[:,np.newaxis,(0,1)] - mol_features[np.newaxis,:,:], axis=2)
	distances += np.linalg.norm(ellipse_vectors[:,np.newaxis,(2,3)] - mol_features[np.newaxis,:,:], axis=2)
	return (distances < ellipse_vectors[:,4].reshape(ellipse_vectors.shape[0], 1)).astype(int)

def calc_mcc(classifications, labels):
	"""
	Calculates the MCC score from the permeation classifications from an ellipse
	and the experimental permeation labels from the data set.
	Args:
		classifications: numpy array containing the ellipse classifications.
		labels: numpy array containing the experimental permeation labels.
	Returns:
		mcc: float, the MCC score.
	"""
	tp = np.sum(np.logical_and(classifications == 1, labels == 1), axis=1)
	tn = np.sum(np.logical_and(classifications == 0, labels == 0), axis=1)
	fp = np.sum(np.logical_and(classifications == 1, labels == 0), axis=1)
	fn = np.sum(np.logical_and(classifications == 0, labels == 1), axis=1)
	mcc = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
	mcc = np.where(mcc == 0.0, 1.0, mcc)
	mcc = ((tp * tn) - (fp * fn)) / mcc
	return mcc

def calc_area(ellipse_vectors):
	"""
	Calculates the area of the ellipse described by its parameters, used to
	penalise overly large ellipses.
	Args:
		ellipse_vectors: numpy array containing the ellipse parameters. 
			Specifically, the ellipse parameters are the x- and y-coordinates of
			the two foci and the length across the major axis of the ellipse.
	Returns:
		float, the area of the ellipse.
	"""
	c2 = np.sum(np.square(ellipse_vectors[:,(0,1)] - ellipse_vectors[:,(2,3)]), axis=1) / 4
	a = ellipse_vectors[:,4] / 2
	b = np.sqrt(np.abs(a**2 - c2))
	return np.pi * a * b

def optimise_ellipse(mol_features, ellipse_vectors, labels, area_penalty=0.05, change_size=0.5, random_mult=0.025):
	"""
	Performs 10 independent runs of 100,000 iteration Monte Carlo optimisation
	on the initial ellipse parameters from the original BOILED-EGG model, now
	with the RDKit values of TPSA and LogP. The ellipse parameter vector that
	is encountered at any point during the optimisation that yields the best MCC
	score is returned.
	Args:
		mol_features: 2d numpy array containing TSPA values for all molecules in
			the first column and LogP values in the second column.
		ellipse_vectors: numpy array containing the ellipse parameters. 
			Specifically, the ellipse parameters are the x- and y-coordinates of
			the two foci and the length across the major axis of the ellipse.
		labels: numpy array containing the experimental permeation labels.
		area_penalty: float, multiplicative factor that is applied to the
			ellipse area to penalise large ellipses.
		change_size: float, proposed updates to ellipse parameter values will be
			within the range -change_size to +change_size.
		random_mult: float, multiplicative factor that is applied to a randomly
			sampled number that is added to the current ellipse parameters'
			score to allow some potentially lower scoring parameter values to be
			accepted.
	Returns:
		best_vector: numpy array containing the ellipse parameters that achieved
			the highest MCC score at any point in the optimisation.
		best_mcc: float, the MCC score of the best classifying ellipse
			encountered at any point in the optimisation.
	"""
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
	return best_vector, best_mcc

def run_egg_optimisation(mol_features, ellipse_vectors, labels, seed, egg_type, area_penalty=0.05, change_size=0.5, random_mult=0.025):
	"""
	Runs the ellipse optimisation routine and prints the results.
	Args:
		mol_features: 2d numpy array containing TSPA values for all molecules in
			the first column and LogP values in the second column.
		ellipse_vectors: numpy array containing the ellipse parameters. 
			Specifically, the ellipse parameters are the x- and y-coordinates of
			the two foci and the length across the major axis of the ellipse.
		labels: numpy array containing the experimental permeation labels.
		seed: int, a random seed for the Monte Carlo optimisation.
		egg_type: str, indicates whether optimising the HIA or BBB ellipse.
		area_penalty: float, multiplicative factor that is applied to the
			ellipse area to penalise large ellipses.
		change_size: float, proposed updates to ellipse parameter values will be
			within the range -change_size to +change_size.
		random_mult: float, multiplicative factor that is applied to a randomly
			sampled number that is added to the current ellipse parameters'
			score to allow some potentially lower scoring parameter values to be
			accepted.
	"""
	np.random.seed(seed)
	best_vector, best_mcc = optimise_ellipse(mol_features, ellipse_vectors, labels, area_penalty=area_penalty, change_size=change_size, random_mult=random_mult)
	print(egg_type, "MCC:", best_mcc, best_vector)


if __name__ == "__main__":
	hia_smiles, hia_labels = read_from_file("hia_mols.csv")
	bbb_smiles, bbb_labels = read_from_file("bbb_mols.csv")
	hia_mols, hia_err_msgs = get_mol_list(hia_smiles)
	bbb_mols, bbb_err_msgs = get_mol_list(bbb_smiles)
	hia_features = get_egg_features(hia_mols)
	bbb_features = get_egg_features(bbb_mols)
	hia_start_vectors = get_start_vectors(142.081, 8.74, np.radians(-1.031325), (71.051, 2.292))
	bbb_start_vectors = get_start_vectors(82.061, 5.557, np.radians(-0.171887), (38.117, 3.177))
	run_egg_optimisation(hia_features, hia_start_vectors, hia_labels, 5, "HIA", change_size=1.0)
	run_egg_optimisation(bbb_features, bbb_start_vectors, bbb_labels, 5, "BBB", change_size=0.5)
