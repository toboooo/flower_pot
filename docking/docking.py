"""Provides the functions that perform molecular docking calculations using
Autodock Vina and GOLD from the command line."""
import os
import shutil
import platform
from subprocess import call
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
from meeko import PDBQTMolecule
from meeko import RDKitMolCreate

def perform_vina_docking(mols, file_names, protein):
	"""
	Embeds each molecule in 3D coordinates, optimizes its geometry with MMFF
	and runs a molecular docking simulation with Autodock Vina, if it is
	available on the system PATH.
	Args:
		mols: list of Mol objects, a list of molecules for which docking
			calculations are to be performed.
		file_names: list of str, a list of the names associated with each
			molecule.
	Returns:
		docking_scores: list of float, the binding affinities in kcal/mol of
			each input molecule, as calculated with the Vina forcefield.
	"""
	# Check that Vina is available
	if shutil.which("vina") is None:
		print("ERROR: 'vina' command was not found on PATH.")
		return None
	# Check required docking files are present
	receptor_file = os.path.join("docking", protein, protein + ".pdbqt")
	config_file = os.path.join("docking", protein, protein + "_config.txt")
	if not os.path.exists(receptor_file):
		print("ERROR: Required docking data file %s not found." % receptor_file)
		return None
	if not os.path.exists(config_file):
		print("ERROR: Required docking data file %s not found." % config_file)
		return None
	# Add hydrogens and optimise 3D geometries of the molecules before saving
	# to PDBQT file.
	ligand_filenames = list()
	os.makedirs("ligand_best_poses", exist_ok=True)
	for i in range(len(mols)):
		mols[i] = Chem.AddHs(mols[i])
		AllChem.EmbedMolecule(mols[i])
		ff_result = AllChem.MMFFOptimizeMolecule(mols[i])
		# Try optimising again with more iterations if needed
		if ff_result == 1:
			print("WARNING: Force field geometry optimisation did not converge for molecule %s, trying again with more steps..." % Chem.MolToSmiles(mols[i]))
			AllChem.MMFFOptimizeMolecule(mols[i], maxIters=2000)
		elif ff_result == -1:
			print("WARNING: Force field could not be set up for molecule %s. Unable to optimise 3D geometry." % Chem.MolToSmiles(mols[i]))
		preparator = MoleculePreparation(merge_these_atom_types=[])
		mol_setups = preparator.prepare(mols[i])
		for setup in mol_setups:
			pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
			if is_ok:
				pdbqt_filename = os.path.join("ligand_best_poses", file_names[i] + ".pdbqt")
				pdbqt_file = open(pdbqt_filename, "w")
				pdbqt_file.write(pdbqt_string)
				pdbqt_file.close()
				ligand_filenames.append(pdbqt_filename)
			else:
				print("ERROR: Could not create PDBQT file for molecule %s. Error message: %s." % (Chem.MolToSmiles(mols[i]), error_msg))
				ligand_filenames.append(None)
	# Perform the docking calculations
	docking_scores = list()
	for ligand_filename in ligand_filenames:
		if ligand_filename is None:
			docking_scores.append(None)
			continue
		output_filename = ligand_filename.replace(".pdbqt", "_out.pdbqt")
		print("Docking %s with %s..." % (ligand_filename, protein))
		call(("vina", "--receptor", receptor_file, "--ligand", ligand_filename, "--config", config_file, "--out", output_filename))
		if os.path.exists(output_filename):
			found_score = False
			with open(output_filename, "r") as output_file:
				for line in output_file:
					if "VINA RESULT" in line:
						line = line.split()
						docking_scores.append(float(line[3]))
						found_score = True
						break
			if not found_score:
				print("WARNING: Could not find VINA RESULT in output file: %s. Could not get docking score." % output_filename)
				docking_scores.append(None)
			else:
				pdbqt_mol = PDBQTMolecule.from_file(output_filename, skip_typing=True)
				rdkit_mol = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)[0]
				sdf_savename = output_filename.replace("_out.pdbqt", "_%s_vina.sdf" % protein)
				with Chem.SDWriter(sdf_savename) as writer:
					writer.write(rdkit_mol)
			try:
				os.remove(output_filename)
			except FileNotFoundError:
				continue
		else:
			print("ERROR: Could not perform docking with vina")
			docking_scores.append(None)
		# Clean up ligand files
		try:
			os.remove(ligand_filename)
		except FileNotFoundError:
			continue
	return docking_scores

def perform_gold_docking(mols, file_names, protein, gold_dir=""):
	"""
	Embeds each molecule in 3D coordinates, optimizes its geometry with MMFF
	and runs a molecular docking simulation with GOLD, if the installation
	directory can be found.
	Args:
		mols: list of Mol objects, a list of molecules for which docking
			calculations are to be performed.
		file_names: list of str, a list of the names associated with each
			molecule.
	Returns:
		docking_scores: list of float, the GOLD docking scores for each input
			molecule.
	"""
	# Check that GOLD is available
	if gold_dir != "" and os.path.exists(os.path.join(gold_dir, "GOLD")):
		gold_dir = os.path.join(gold_dir, "GOLD")
		print("Found GOLD installation at: %s" % gold_dir)
	elif gold_dir != "" and not os.path.exists(os.path.join(gold_dir, "GOLD")):
		print("WARNING: Could not find GOLD installation from specified path: %s." % gold_dir)
		gold_dir = None
	else:
		gold_dir = None
	if gold_dir is None and os.path.exists("GOLD_INSTALLDIR.txt"):
		path_file = open("GOLD_INSTALLDIR.txt")
		gold_dir = next(path_file).rstrip("\n")
		path_file.close()
		if os.path.exists(os.path.join(gold_dir, "GOLD")):
			gold_dir = os.path.join(gold_dir, "GOLD")
			print("Found GOLD installation at: %s, from GOLD_INSTALLDIR.txt" % gold_dir)
		else:
			print("WARNING: Could not find GOLD installation at: %s, from file GOLD_INSTALLDIR.txt" % gold_dir)
			gold_dir = None
	if platform.system() == "Windows":
		guess_dir = "C:\\Program Files\\CCDC\\ccdc-software\\gold"
		bin_path = "gold\\d_win32\\bin\\gold_win32.exe"
	elif platform.system() == "Darwin":
		guess_dir = "/Applications/CCDC"
		bin_path = "bin/gold_auto"
	elif platform.system() == "Linux":
		guess_dir = os.path.join("/home", os.getlogin(), "CCDC")
		bin_path = "bin/gold_auto"
	if gold_dir is None:
		print("Searching for GOLD installation in guess directory: %s..." % guess_dir)
		if not os.path.exists(guess_dir):
			print("ERROR: Could not find GOLD installation directory.")
			return None
		else:
			for folder in os.listdir(guess_dir):
				print(folder)
				if folder == "GOLD":
					gold_dir = os.path.join(guess_dir, folder)
					print("Found GOLD installation directory at: %s!" % gold_dir)
					break
		if gold_dir is None:
			print("ERROR: Could not find GOLD installation directory.")
			return None
		gold_command = os.path.join(gold_dir, bin_path)
	else:
		gold_command = os.path.join(gold_dir, bin_path)
	if not os.path.exists(gold_command):
		print("ERROR: GOLD installation directory found, but executable %s was not." % gold_command)
		return None
	# Check required data files are present
	if not os.path.exists(os.path.join("docking", protein, protein + ".mol2")):
		print("ERROR: Required docking data file %s not found." % os.path.join("docking", protein, protein + ".mol2"))
		return None
	if not os.path.exists(os.path.join("docking", protein, "gold.conf")):
		print("ERROR: Required docking data file %s not found." % os.path.join("docking", protein, "gold.conf"))
		return None
	# Add hydrogens and optimise 3D geometries of the molecules before saving
	docking_scores = list()
	for i in range(len(mols)):
		mols[i] = Chem.AddHs(mols[i])
		AllChem.EmbedMolecule(mols[i])
		ff_result = AllChem.MMFFOptimizeMolecule(mols[i])
		# Try optimising again with more iterations if needed
		if ff_result == 1:
			AllChem.MMFFOptimizeMolecule(mols[i], maxIters=2000)
		elif ff_result == -1:
			print("WARNING: Force field could not be set up for molecule %s. Unable to optimise 3D geometry." % Chem.MolToSmiles(mols[i]))
		with Chem.SDWriter("ligand.sdf") as sdf_writer:
			sdf_writer.write(mols[i])
		print("Docking molecule: %s with target protein: %s" % (Chem.MolToSmiles(mols[i]), protein))
		call((gold_command, os.path.join("docking", protein, "gold.conf")))
		if not os.path.exists("bestranking.lst"):
			print("ERROR: GOLD did not write file 'bestranking.lst' for molecule: %s." % Chem.MolToSmiles(mols[i]))
			docking_scores.append(None)
		else:
			found_score = False
			with open("bestranking.lst", "r") as ranking_file:
				for line in ranking_file:
					if line[0] == "#" or line == "\n":
						continue
					line = line.rstrip("\n").split()
					try:
						docking_scores.append(float(line[0]))
						found_score = True
						break
					except ValueError:
						print("WARNING: Invalid docking score found in 'bestranking.lst'.")
			if not found_score:
				print("ERROR: GOLD was unable to perform docking for molecule: %s with target protein: %s." % (Chem.MolToSmiles(mols[i]), protein))
				docking_scores.append(None)
	try:
		os.remove("ligand.sdf")
		os.remove("gold.log")
		os.remove("bestranking.lst")
	except FileNotFoundError:
		pass
	return docking_scores

def estimate_ic50(docking_scores, protein, program="vina"):
	"""
	Provides a pedagogical estimation of IC50 of a drug molecule
	by calculating the difference in binding affinities between
	the current molecule and a reference, and then passing this
	difference into an exponentially decaying curve, such that
	initial improvements in a molecule's binding affinity are
	rewarded with large decreases in IC50, whereas over-
	optimisations beyond the reference structure only yield
	small returns.
	Args:
		docking_scores: list of floats, the binding affinities
			of the molecules for which IC50 values are
			to be returned.
		program: str, the program used to calculate the docking scores.
	Returns:
		list of float, a list of the pedagogical IC50 values.
	"""
	prefactor = 0.5
	scale = 10
	ic50_estimations = list()
	if program == "vina":
		reference_values = {"26S_proteasome": -6.8, "BTHalpha": -12.1, "CF-IIbeta": -9.2, "CHK1": -10.4, "CYP17a": -11.6}
		reference = reference_values[protein]
		for score in docking_scores:
			if score is not None:
				difference = reference - score
				ic50 = prefactor * 2**(-difference / scale)
				ic50_estimations.append(ic50)
			else:
				ic50_estimations.append(None)
	elif program == "gold":
		reference_values = {"26S_proteasome": 71.78, "BTHalpha": 118.88, "CF-IIbeta": 96.62, "CHK1": 39.72, "CYP17a": 23.18}
		reference = reference_values[protein]
		for score in docking_scores:
			if score is not None:
				difference = score - reference
				ic50 = prefactor * 2**(-difference / scale)
				ic50_estimations.append(ic50)
			else:
				ic50_estimations.append(None)
	return ic50_estimations
