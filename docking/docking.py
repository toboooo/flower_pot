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

OVERWRITE_INSTRUCTIONS = {"SAVE OPTIONS": "save_score_in_file = 1\nsave_protein_torsions = 1\nclean_up_option delete_empty_directories\nclean_up_option delete_redundant_log_files\nclean_up_option save_top_n_solutions 1\nclean_up_option delete_all_initialised_ligands\n\n",
"WRITE OPTIONS": "write_options = NO_LOG_FILES NO_LINK_FILES NO_RNK_FILES NO_GOLD_SOLN_LIGAND_MOL2_FILES NO_GOLD_PROTEIN_MOL2_FILE NO_LGFNAME_FILE NO_PLP_MOL2_FILES NO_PID_FILE NO_SEED_LOG_FILE NO_GOLD_ERR_FILE NO_FIT_PTS_FILES NO_ASP_MOL2_FILES NO_GOLD_LOG_FILE\n\n",
"DATA FILES": "ligand_data_file ligand.sdf 10\nparam_file = DEFAULT\nset_ligand_atom_types = 1\nset_protein_atom_types = 0\ndirectory = .\ntordist_file = DEFAULT\nmake_subdirs = 0\nsave_lone_pairs = 1\nread_fitpts = 0\n\n"}

def make_vina_ligands(mols, file_names):
	"""
	Writes pdbqt files for all provided ligands for later docking with Vina.
	Args:
		mols: list of Mol objects, a list of molecules for which docking
			calculations are to be performed.
		file_names: list of str, a list of the names associated with each
			of the molecules.
	Return:
		ligand_filenames: list of str, a list of pdbqt files for each provided
			molecule.
	"""
	ligand_filenames = list()
	for i in range(len(mols)):
		# Add hydrogens and optimise 3D geometries of the molecules before saving to PDBQT file.
		mols[i] = Chem.AddHs(mols[i])
		AllChem.EmbedMolecule(mols[i])
		ff_result = AllChem.MMFFOptimizeMolecule(mols[i])
		# Try optimising again with more iterations if needed
		if ff_result == 1:
			print("WARNING: Force field geometry optimisation did not converge for molecule %s, trying again with more steps..." % file_names[i])
			AllChem.MMFFOptimizeMolecule(mols[i], maxIters=2000)
		elif ff_result == -1:
			print("WARNING: Force field could not be set up for molecule %s. Unable to optimise 3D geometry." % file_names[i])
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
				print("ERROR: Could not create PDBQT file for molecule %s. Error message: %s." % (file_names[i], error_msg))
				ligand_filenames.append(None)
	return ligand_filenames

def run_vina_docking(ligand_filenames, receptor_file, config_file, protein):
	"""
	Runs Vina docking calculations for each prepared ligand for the specified
	receptor and configuration.
	Args:
		ligand_filenames: list of str, paths for each of the prepared pdbqt
			files for each molecule to be docked.
		receptor_file: str, path for the protein pdbqt file on which docking is
			to be performed.
		config_file: str, path for the configuration file for the Vina docking.
		protein: str, name of the protein against which docking calculations are
			to be performed.
	Returns:
		docking_scores: list of float, the calculated binding affinities for
			each molecule.
	"""
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
				with Chem.SDWriter(sdf_savename) as sdf_writer:
					sdf_writer.write(rdkit_mol)
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
		protein: str, name of the protein against which docking calculations are
			to be performed.
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
	os.makedirs("ligand_best_poses", exist_ok=True)
	# Prepare ligand files for all molecules
	ligand_filenames = make_vina_ligands(mols, file_names)
	# Perform the docking calculations
	docking_scores = run_vina_docking(ligand_filenames, receptor_file, config_file, protein)
	return docking_scores

def custom_vina_docking(mols, file_names, protein_filename, config_file):
	"""
	Perfoms Vina docking calculations using a user specified protein and
	configuration, rather than one of five built-in proteins.
	Args:
		mols: list of Mol objects, a list of molecules for which docking
			calculations are to be performed.
		file_names: list of str, a list of the names associated with each
			molecule.
		protein_filename: str, path to a pdbqt file of the protein for which
			docking calculations are to be performed.
		config_file: str, path to a Vina configuration file for the custom
			calculation.
	Returns:
		docking_scores: list of float, the binding affinities in kcal/mol of
			each input molecule, using the user provided protein and
			configuration.
	"""
	# Check that Vina is available, note protein and config files should already have been checked for existence
	if shutil.which("vina") is None:
		print("ERROR: 'vina' command was not found on PATH.")
		return None
	os.makedirs("ligand_best_poses", exist_ok=True)
	# Prepare ligand files for all molecules
	ligand_filenames = make_vina_ligands(mols, file_names)
	# Perform the docking calculations
	docking_scores = run_vina_docking(ligand_filenames, protein_filename, config_file, "custom")
	return docking_scores

def check_gold_avail(gold_dir):
	"""
	Check if GOLD is installed and return a path to the GOLD command if found.
	Args:
		gold_dir: str, optionally supplied directory for GOLD installation.
	Returns:
		gold_command: str, full path to GOLD executable, if found.
	"""
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
		print("WARNING: GOLD usage has only been tested on Windows. Please note that using GOLD through Flower Pot may not work on macOS.")
		guess_dir = "/Applications/CCDC"
		bin_path = "bin/gold_auto"
	elif platform.system() == "Linux":
		print("WARNING: GOLD usage has only been tested on Windows. Please note that using GOLD through Flower Pot may not work on Linux.")
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
	return gold_command

def make_gold_mol(mol, file_name):
	"""
	Generates 3D coordinates for a molecule using MMFF94 and writes the
	coordinates to an sdf file for GOLD.
	Args:
		mol: Mol object, the current molecule for which docking is to be
			performed.
		file_name: str, name to which the current molecule is associated.
	"""
	mol = Chem.AddHs(mol)
	AllChem.EmbedMolecule(mol)
	ff_result = AllChem.MMFFOptimizeMolecule(mol)
	# Try optimising again with more iterations if needed
	if ff_result == 1:
		AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
	elif ff_result == -1:
		print("WARNING: Force field could not be set up for molecule %s. Unable to optimise 3D geometry." % file_name)
	with Chem.SDWriter("ligand.sdf") as sdf_writer:
		sdf_writer.write(mol)

def process_gold_results(docking_scores, protein, file_name):
	"""
	Gather results from GOLD output files and record docking scores if
	successful.
	Args:
		docking_scores: list of floats, a list of currently collected docking
			scores.
		protein: str, name of protein.
		file_name: str, the name associated with the current molecule to be
			docked.
	"""
	if not os.path.exists("bestranking.lst"):
		print("ERROR: GOLD did not write file 'bestranking.lst' for molecule: %s." % file_name)
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
					if len(line) >= 10:
						output_file = line[9].strip("'")
					else:
						print("WARNING: No GOLD output file was found for molecule: %s." % file_name)
						break
					if os.path.exists(output_file):
						shutil.move(output_file, "ligand_best_poses/%s" % (file_name + "_" + protein + "_gold.sdf"))
					else:
						print("WARNING: No GOLD output file was found for molecule: %s." % file_name)
					break
				except ValueError:
					print("WARNING: Invalid docking score found in 'bestranking.lst'.")
		if not found_score:
			print("ERROR: GOLD was unable to perform docking for molecule: %s with target protein: %s." % (file_name, protein))
			docking_scores.append(None)

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
		gold_dir: str, optionally supplied directory for GOLD installation.
	Returns:
		docking_scores: list of float, the GOLD docking scores for each input
			molecule.
	"""
	# Check that GOLD is available
	gold_command = check_gold_avail(gold_dir)
	if gold_command is None:
		return None
	elif not os.path.exists(gold_command):
		print("ERROR: GOLD installation directory found, but executable %s was not." % gold_command)
		return None
	# Check required data files are present
	if not os.path.exists(os.path.join("docking", protein, protein + ".mol2")):
		print("ERROR: Required docking data file %s not found." % os.path.join("docking", protein, protein + ".mol2"))
		return None
	if not os.path.exists(os.path.join("docking", protein, "gold.conf")):
		print("ERROR: Required docking data file %s not found." % os.path.join("docking", protein, "gold.conf"))
		return None
	os.makedirs("ligand_best_poses", exist_ok=True)
	docking_scores = list()
	for i in range(len(mols)):
		make_gold_mol(mols[i], file_names[i])
		print("Docking molecule: %s with target protein: %s" % (file_names[i], protein))
		call((gold_command, os.path.join("docking", protein, "gold.conf")))
		process_gold_results(docking_scores, protein, file_names[i])
	try:
		os.remove("ligand.sdf")
		os.remove("bestranking.lst")
	except FileNotFoundError:
		pass
	return docking_scores

def process_gold_conf(gold_conf_file):
	"""
	Overwrites the 'save options', 'write options' and 'data file' fields of a
	gold.conf file for custom docking jobs such that GOLD cleans up after itself
	and this program can process the output. The user MUST make sure that the
	gold.conf file is valid and otherwise suited to their requirements.
	Args:
		gold_conf_file: str, the file name of the gold.conf file for a custom
			job. This function assumes that the file exists and is valid.
	"""
	overwrite_fields = {"SAVE OPTIONS": False, "WRITE OPTIONS": False,  "DATA FILES": False}
	conf_file = open(gold_conf_file, "r")
	lines = conf_file.readlines()
	conf_file.close()
	shutil.move(gold_conf_file, gold_conf_file + ".old")
	conf_file = open(gold_conf_file, "w")
	lines = iter(lines)
	while True:
		try:
			line = next(lines)
			possible_field = line.strip()
			if possible_field in overwrite_fields.keys() and not overwrite_fields[possible_field]:
				conf_file.write(line)
				conf_file.write(OVERWRITE_INSTRUCTIONS[possible_field])
				overwrite_fields[possible_field] = True
				line = next(lines).strip()
				while line != "" and line[0].isalpha():
					line = next(lines).strip()
			else:
				conf_file.write(line)
		except StopIteration:
			break
	for field in overwrite_fields:
		if not overwrite_fields[field]:
			conf_file.write("  " + field + "\n")
			conf_file.write(OVERWRITE_INSTRUCTIONS[field])
	conf_file.close()

def custom_gold_docking(mols, file_names, conf_filename, gold_dir):
	"""
	Performs GOLD docking score calculations using user specified protein files
	and configuration, rather than one of five in-built proteins.
	Args:
		mols: list of Mol objects, a list of molecules for which docking
			calculations are to be performed.
		file_names: list of str, a list of the names associated with each
			molecule.
		conf_filename: str, path to a gold.conf file for the custom calculation.
		gold_dir: str, optionally supplied directory for GOLD installation.
	Returns:
		docking_scores: list of float, the docking scores for each molecule
			calculated for the custom protein.
	"""
	# Check that GOLD is available, note that protein mol2 and gold.conf files should have already been checked for existence
	gold_command = check_gold_avail(gold_dir)
	if gold_command is None:
		return None
	elif not os.path.exists(gold_command):
		print("ERROR: GOLD installation directory found, but executable %s was not." % gold_command)
		return None
	# Edit gold.conf file
	process_gold_conf(conf_filename)
	os.makedirs("ligand_best_poses", exist_ok=True)
	docking_scores = list()
	for i in range(len(mols)):
		make_gold_mol(mols[i], file_names[i])
		print("Docking molecule: %s with configuration file: %s" % (file_names[i], conf_filename))
		call((gold_command, conf_filename))
		process_gold_results(docking_scores, "custom", file_names[i])
	try:
		os.remove("ligand.sdf")
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
	vina_scale = 0.2
	gold_scale = 7.5
	ic50_estimations = list()
	if program == "vina":
		reference_values = {"26S_proteasome": -6.8, "BTHalpha": -10.7, "CF-IIbeta": -9.1, "CHK1": -9.8, "CYP17a": -10.3}
		reference = reference_values[protein]
		for score in docking_scores:
			if score is not None:
				difference = reference - score
				ic50 = prefactor * 2**(-difference / vina_scale)
				ic50_estimations.append(ic50)
			else:
				ic50_estimations.append(None)
	elif program == "gold":
		reference_values = {"26S_proteasome": 71.78, "BTHalpha": 118.88, "CF-IIbeta": 96.62, "CHK1": 39.72, "CYP17a": 23.18}
		reference = reference_values[protein]
		for score in docking_scores:
			if score is not None:
				difference = score - reference
				ic50 = prefactor * 2**(-difference / gold_scale)
				ic50_estimations.append(ic50)
			else:
				ic50_estimations.append(None)
	return ic50_estimations
