"""Provides the functions to implement the substructure filtering routines."""
import os
import sys
import csv
import pathlib
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import Draw

def resource_path(relative_path):
	"""
	Returns a path to a specified file, with the base path modified to allow
	pyinstaller to be able to locate the file.
	Args:
		relative_path: str, the relative path for the required file.
	Returns:
		Path to requested file, modified for pyinstaller if needed.
	"""
	try:
		base_path = sys._MEIPASS
	except Exception:
		base_path = os.path.abspath(".")
	return os.path.join(base_path, relative_path)

def get_filters(filter_name):
	"""
	Read in the SMARTS strings for the substructures for the specified filter.
	Args:
		filter_name: str, the name of the desired substructure filter.
	Returns:
		patterns: list of Mol objects, a list of the substructures that may be
			pattern matched against a complete structure.
		pattern_names: list of str, a list of the names of each substructure in
			the filter.
	"""
	smarts_strings = []
	pattern_names = []
	file = open(resource_path("substructure/" + filter_name.lower() + ".csv"), "r")
	csv_reader = csv.DictReader(file)
	for line in csv_reader:
		pattern_names.append(line["name"])
		smarts_strings.append(line["smarts"])
	file.close()
	patterns = [Chem.MolFromSmarts(smarts) for smarts in smarts_strings]
	return patterns, pattern_names

def check_filters(mols, file_names, filter, print_images=False):
	"""
	Applies a substructure filter to a series of molecules, returns counts of
	substructure hits and optionally draws images with hits highlighted.
	Args:
		mols: list of Mol objects, a list of molecules for which filters are to
			be applied.
		file_names: list of str, a list of the names associated with each
			molecule.
		filter: str, the name of the particular filter to be applied.
		print_images: optional, bool, whether or not images will be drawn.
	Returns:
		image_names: list of str, the names of the image files that were drawn.
		hit_names: list of str, the names of the substructures that were
			detected in each molecule.
		hit_counts: list of int, the numbers of substructures in the filter that
			each molecule contained.
	"""
	pathlib.Path("images").mkdir(parents=True, exist_ok=True)
	image_names = []
	hit_names = []
	hit_counts = []
	patterns, pattern_names = get_filters(filter)
	for mol, file_name in zip(mols, file_names):
		matched_substructs = []
		substruct_names = []
		hit_count = 0
		for pattern, name in zip(patterns, pattern_names):
			if mol.HasSubstructMatch(pattern):
				substruct_names.append(name)
				hit_count += 1
				if print_images:
					matched_substructs.append(mol.GetSubstructMatch(pattern))
		hit_names.append(substruct_names)
		hit_counts.append(hit_count)
		if print_images and len(matched_substructs) > 0:
			img = Draw.MolsToGridImage([mol for i in range(len(matched_substructs))], molsPerRow=len(matched_substructs), subImgSize=(200,200), highlightAtomLists=matched_substructs, legends=[filter + ": " + substruct_name for substruct_name in substruct_names])
			timecode = datetime.now().strftime("%d-%m-%H%M")
			image_name = os.path.join("images", timecode + "_" + file_name + "_" + filter + ".png")
			image_names.append(image_name)
			img.save(image_name)
		else:
			image_names.append(None)
	return image_names, hit_names, hit_counts
