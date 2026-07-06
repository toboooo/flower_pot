"""Functions for the assay simulation pop-up window. Heavily based on code by
Patrick Taylor (github.com/PatrickJTaylor) for the publication:
https://pubs.acs.org/doi/10.1021/acs.jchemed.3c00066"""
import os
import copy
import math
import tkinter as tk
import tkinter.ttk as ttk
from datetime import datetime
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageTk
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
NavigationToolbar2Tk)
from matplotlib.figure import Figure

PLOT_N_POINTS = 100

REFERENCE_MOLS = {"5IF3 Vina score": Chem.MolFromSmiles(
		"CC(C)C[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)c1cnccn1)C(O)O"
	),
	"5IF3 GOLD score": Chem.MolFromSmiles(
		"CC(C)C[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)c1cnccn1)B(O)O"
	),
	"1IEP Vina score": Chem.MolFromSmiles(
		"Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"
	),
	"1IEP GOLD score": Chem.MolFromSmiles(
		"Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
	),
	"2W26 Vina score": Chem.MolFromSmiles(
		"O=C(NC[C@H]1CN(c2ccc(N3CCOCC3=O)cc2)C(=O)O1)c1ccc(Cl)s1"
	),
	"2W26 GOLD score": Chem.MolFromSmiles(
		"O=C(NC[C@H]1CN(c2ccc(N3CCOCC3=O)cc2)C(=O)O1)c1ccc(Cl)s1"
	),
	"1ZYS Vina score": Chem.MolFromSmiles(
		"CN1CCN(c2ccc(-c3cnc4[nH]cc(NC(=O)c5cccnc5)c4c3)cc2)CC1"
	),
	"1ZYS GOLD score": Chem.MolFromSmiles(
		"CN1CCN(c2ccc(-c3cnc4[nH]cc(NC(=O)c5cccnc5)c4c3)cc2)CC1"
	),
	"3RUK Vina score": Chem.MolFromSmiles(
		"C[C@]12CC[C@H](O)CC1=CC[C@@H]1[C@@H]2CC[C@]2(C)C(c3cccnc3)=CC[C@@H]12"
	),
	"3RUK GOLD score": Chem.MolFromSmiles(
		"C[C@]12CC[C@H](O)CC1=CC[C@@H]1[C@@H]2CC[C@]2(C)C(c3cccnc3)=CC[C@@H]12"
	)
}

REFERENCE_SCORES = {"5IF3 Vina score": -6.8, "1IEP Vina score": -10.7,
	"2W26 Vina score": -9.1, "1ZYS Vina score": -9.8, "3RUK Vina score": -10.3,
	"5IF3 GOLD score": 71.78, "1IEP GOLD score": 118.88,
	"2W26 GOLD score": 96.62, "1ZYS GOLD score": 39.72,
	"3RUK GOLD score": 23.18}

def round_sigfigs(x, sigfigs=3):
	return round(x, sigfigs - int(math.floor(math.log10(abs(x)))) - 1)

def validate_input_fields():
	"""
	Validates that the input variables for the simulator are correct before they
	may be used to update the assay plot.
	Returns:
		bool
	"""
	global ic50_box_value
	try:
		input_ic50 = ic50_box_value.get()
	except tk._tkinter.TclError: 
		print("ERROR: Invalid IC50 value - not a valid number.")
		return False
	if input_ic50 <= 0.0:
		print("ERROR: Invalid IC50 value - must be greater than zero.")
		return False
	global start_conc_box
	try:
		start_conc = float(start_conc_box.get())
	except ValueError:
		print("ERROR: Invalid starting concentration - not a valid number.")
		return False
	if start_conc <= 0.0:
		print("ERROR: Invalid start concentration - must be greater than zero.")
		return False
	global dilution_factor_box
	try:
		dilution_factor = float(dilution_factor_box.get())
	except ValueError:
		print("ERROR: Invalid dilution factor - not a valid number.")
		return False
	if dilution_factor <= 0:
		print("ERROR: Invalid dilution factor - must be greater than zero.")
		return False
	elif dilution_factor >= 1.0:
		print("ERROR: Invalid dilution factor - must be less than one.")
		return False
	global n_dilutions_box
	try:
		n_dilutions = int(n_dilutions_box.get())
	except ValueError:
		print("ERROR: Invalid number of dilutions - not a valid number.")
		return False
	if n_dilutions <= 0:
		print("ERROR: Invalid number of dilutions - must be greater than zero.")
		return False
	global min_response_box
	try:
		min_response = float(min_response_box.get())
	except ValueError:
		print("ERROR: Invalid minimum response value - not a valid number.")
		return False
	global max_response_box
	try:
		max_response = float(max_response_box.get())
	except ValueError:
		print("ERROR: Invalid maximum response value - not a valid number.")
		return False
	if max_response <= min_response:
		print("ERROR: Invalid maximum response value - must be greater than "\
			"minimum response value.")
		return False
	global num_repeats_box
	try:
		num_repeats = int(num_repeats_box.get())
	except ValueError:
		print("ERROR: Invalid number of repeats - not a valid number.")
		return False
	if num_repeats <= 0:
		print("ERROR: Invalid number of repeats - must be greater than zero.")
		return False
	global noise_sd_box
	try:
		noise_sd = float(noise_sd_box.get())
	except ValueError:
		print("ERROR: Invalid noise standard deviation - not a valid number.")
		return False
	if noise_sd <= 0.0:
		print("ERROR: Invalid noise standard deviation - must be greater than "\
			"zero.")
		return False
	return True

def line_function(x, c):
	"""
	Calculates a horizontal line function given an intercept value.
	Args:
		x: numpy array of float, the input variable values.
		c: float, the intercept value of the line.
	Returns:
		numpy array of float.
	"""
	return np.repeat(c, x.shape)

def sigmoid_function(x, m, bl, bu):
	"""
	Calculates a (log-domain) sigmoid function given a midpoint, and a lower and
	upper bound.
	Args:
		x: numpy array of float, the input variabl values.
		m: float, the midpoint of the curve.
		bl: float, the lower bound of the curve.
		bu: float, the upper bound of the curve.
	Returns:
		numpy array of float.
	"""
	return bl + (bu - bl) / (1.0 + (x / m))

def update_simulation_plot(docking_score_keys):
	"""
	Handles the updates to the simulation plot when its variables are changed.
	Args:
		docking_score_keys: tuple of str, the keys for the properties dictionary
			corresponding to the calculated docking scores.
	"""
	global ic50_box_value
	global ic50_guesses
	global protein_index, mol_index
	global start_conc_box, dilution_factor_box, n_dilutions_box
	global min_response_box, max_response_box, num_repeats_box, noise_sd_box
	global function_type
	if not validate_input_fields():
		return
	ic50_value = ic50_box_value.get()
	ic50_guesses[docking_score_keys[protein_index]][mol_index] = ic50_value
	start_conc = float(start_conc_box.get())
	dilution_factor = float(dilution_factor_box.get())
	n_dilutions = int(n_dilutions_box.get())
	concs = np.array([start_conc * dilution_factor**i \
		for i in range(n_dilutions)])
	if function_type == "sigmoid":
		min_response = float(min_response_box.get())
		max_response = float(max_response_box.get())
		response = sigmoid_function(concs, ic50_value, min_response,
			max_response)
	elif function_type == "line":
		response = line_function(concs, ic50_value)
	num_repeats = int(num_repeats_box.get())
	noise_sd = float(noise_sd_box.get())
	clean_responses = np.tile(response, (num_repeats, 1))
	noisy_responses = clean_responses + np.random.normal(loc=0.0,
		scale=noise_sd, size=clean_responses.shape)
	mean_noisy_responses = np.mean(noisy_responses, axis=0)
	std_noisy_responses = np.std(noisy_responses, axis=0)
	curve_concs = np.logspace(np.log10(np.min(concs)), np.log10(np.max(concs)),
		PLOT_N_POINTS)
	if function_type == "sigmoid":
		curve_response = sigmoid_function(curve_concs, ic50_value, min_response,
			max_response)
	elif function_type == "line":
		curve_response = line_function(curve_concs, ic50_value)
	global assay_plot_canvas, assay_plot_ax
	assay_plot_ax.clear()
	assay_plot_ax.set_xlabel("Concentration / M")
	assay_plot_ax.set_ylabel("Assay Response")
	assay_plot_ax.set_xscale("log")
	assay_plot_ax.plot(curve_concs, curve_response, color="C0")
	assay_plot_ax.errorbar(concs, mean_noisy_responses, std_noisy_responses,
		fmt='o', color="C0")
	assay_plot_canvas.draw()

def update_protein_fields(file_names, properties, docking_score_keys):
	"""
	Updates the data displayed information corresponding to the now current
	protein after the 'Prev Protein' or 'Next protein' buttons are clicked.
	Args:
		file_names: list of str, the names of the files for each molecule.
		properties: dict, the computed properties of all of the molecules.
		docking_score_keys: tuple of str, the keys for the properties dictionary
			corresponding to the calculated docking scores.
	"""
	global mol_index, protein_index
	global protein_name_label
	global molecule_image_label, molecule_score_label
	global reference_image, reference_image_label, reference_score_label
	docking_score_key = docking_score_keys[protein_index]
	protein_name_label.configure(text=docking_score_key.rstrip(" score"))
	molecule_score_label.configure(text=file_names[mol_index] + '\n' + \
		docking_score_key + ":\n" + \
		str(properties[docking_score_key][mol_index]))
	if REFERENCE_MOLS.get(docking_score_key) is not None:
		reference_image = ImageTk.PhotoImage(Draw.MolToImage(
			REFERENCE_MOLS[docking_score_key], size=(100,100)))
		reference_image_label.configure(image=reference_image, text=None)
		reference_score_label.configure(text="Reference\n" + \
			docking_score_key + ":\n" + \
			str(REFERENCE_SCORES[docking_score_key]))
	else:
		reference_image = ImageTk.PhotoImage(Image.new("RGBA", (100, 100),
			(0, 0, 0, 0)))
		reference_image_label.configure(image=reference_image, text=None)
		reference_score_label.configure(text="No reference\ncompound for\n" + \
			docking_score_key.rstrip(" score"))
	global ic50_box_value, ic50_guesses
	ic50_box_value.set(ic50_guesses[docking_score_key][mol_index])
	update_simulation_plot(docking_score_keys)

def set_prev_protein(file_names, properties, docking_score_keys):
	"""
	Handles the functionality of the 'Prev Protein' button to change to the
	data for the previous protein.
	Args:
		file_names: list of str, the names of the files for each molecule.
		properties: dict, the computed properties of all of the molecules.
		docking_score_keys: tuple of str, the keys for the properties dictionary
			corresponding to the calculated docking scores.
	"""
	global protein_index
	if len(docking_score_keys) > 1:
		protein_index = (protein_index - 1) % len(docking_score_keys)
		update_protein_fields(file_names, properties, docking_score_keys)

def set_next_protein(file_names, properties, docking_score_keys):
	"""
	Handles the functionality of the 'Next Protein' button to change to the
	data for the next protein.
	Args:
		file_names: list of str, the names of the files for each molecule.
		properties: dict, the computed properties of all of the molecules.
		docking_score_keys: tuple of str, the keys for the properties dictionary
			corresponding to the calculated docking scores.
	"""
	global protein_index
	if len(docking_score_keys) > 1:
		protein_index = (protein_index + 1) % len(docking_score_keys)
		update_protein_fields(file_names, properties, docking_score_keys)

def update_molecule_fields(mols, file_names, properties, docking_score_keys):
	"""
	Updates the data displayed information corresponding to the now current
	protein after the 'Prev Protein' or 'Next protein' buttons are clicked.
	Args:
		mols: list of RDKit Mol objects, the molecules from which the properties
			were calculated.
		file_names: list of str, the names of the files for each molecule.
		properties: dict, the computed properties of all of the molecules.
		docking_score_keys: tuple of str, the keys for the properties dictionary
			corresponding to the calculated docking scores.
	"""
	global mol_index, protein_index
	global molecule_name_label
	global molecule_image, molecule_image_label, molecule_score_label
	docking_score_key = docking_score_keys[protein_index]
	molecule_name_label.configure(text=file_names[mol_index])
	flat_mol = Chem.MolFromSmiles(Chem.MolToSmiles(mols[mol_index]))
	molecule_image = ImageTk.PhotoImage(Draw.MolToImage(flat_mol,
		size=(100,100)))
	molecule_image_label.configure(image=molecule_image, text=None)
	molecule_score_label.configure(text=file_names[mol_index] + '\n' + \
		docking_score_key + ":\n" + \
		str(properties[docking_score_key][mol_index]))
	global ic50_box_value, ic50_guesses
	ic50_box_value.set(ic50_guesses[docking_score_key][mol_index])
	update_simulation_plot(docking_score_keys)

def set_prev_molecule(mols, file_names, properties, docking_score_keys):
	"""
	Handles the functionality of the 'Prev Molecule' button to change to the
	data for the next molecule.
	Args:
		mols: list of RDKit Mol objects, the molecules from which the properties
			were calculated.
		file_names: list of str, the names of the files for each molecule.
		properties: dict, the computed properties of all of the molecules.
		docking_score_keys: tuple of str, the keys for the properties dictionary
			corresponding to the calculated docking scores.
	"""
	global mol_index
	if len(mols) > 1:
		mol_index = (mol_index - 1) % len(mols)
		update_molecule_fields(mols, file_names, properties, docking_score_keys)

def set_next_molecule(mols, file_names, properties, docking_score_keys):
	"""
	Handles the functionality of the 'Next Molecule' button to change to the
	data for the next molecule.
	Args:
		mols: list of RDKit Mol objects, the molecules from which the properties
			were calculated.
		file_names: list of str, the names of the files for each molecule.
		properties: dict, the computed properties of all of the molecules.
		docking_score_keys: tuple of str, the keys for the properties dictionary
			corresponding to the calculated docking scores.
	"""
	global mol_index
	if len(mols) > 1:
		mol_index = (mol_index + 1) % len(mols)
		update_molecule_fields(mols, file_names, properties, docking_score_keys)

def reset_ic50_values(docking_score_keys):
	"""
	Resets all of the IC50 values to their original guesses.
	Args:
		docking_score_keys: tuple of str, the keys for the properties dictionary
			corresponding to the calculated docking scores.
	"""
	global ic50_box_value
	global ic50_guesses, original_ic50_guesses
	global protein_index, mol_index
	ic50_guesses = copy.deepcopy(original_ic50_guesses)
	ic50_box_value.set(
		ic50_guesses[docking_score_keys[protein_index]][mol_index])
	update_simulation_plot(docking_score_keys)

def reset_all_variables(docking_score_keys):
	"""
	Resets all of the variables controlling the assay plot.
	Args:
		docking_score_keys: tuple of str, the keys for the properties dictionary
			corresponding to the calculated docking scores.
	"""
	global start_conc_box, dilution_factor_box, n_dilutions_box
	global min_response_box, max_response_box, num_repeats_box, noise_sd_box
	start_conc_box.delete(0, "end")
	start_conc_box.insert(0, 0.005)
	dilution_factor_box.delete(0, "end")
	dilution_factor_box.insert(0, 0.2)
	n_dilutions_box.delete(0, "end")
	n_dilutions_box.insert(0, 8)
	min_response_box.delete(0, "end")
	min_response_box.insert(0, 50)
	max_response_box.delete(0, "end")
	max_response_box.insert(0, 200)
	num_repeats_box.delete(0, "end")
	num_repeats_box.insert(0, 3)
	noise_sd_box.delete(0, "end")
	noise_sd_box.insert(0, 15)
	update_simulation_plot(docking_score_keys)

def set_function_type(type_str, docking_score_keys):
	"""
	Sets the global function type for the assay simulation.
	Args:
		type_str: str, string representing the type of the assay function,
			either 'sigmoid' or 'line'.
		docking_score_keys: tuple of str, the keys for the properties dictionary
			corresponding to the calculated docking scores.
	"""
	global function_type
	function_type = type_str
	update_simulation_plot(docking_score_keys)

def write_assay_data(file_names, docking_score_keys):
	"""
	Writes the simulated assay data to the output file and textbox.
	Args:
		file_names: list of str, the names of the files for each molecule.
		docking_score_keys: tuple of str, the keys for the properties dictionary
			corresponding to the calculated docking scores.
	"""
	global protein_index, mol_index
	global ic50_guesses
	global ic50_box_value, start_conc_box, dilution_factor_box, n_dilutions_box
	global num_repeats_box, min_response_box, max_response_box, noise_sd_box
	global write_all_button_value
	global function_type
	ic50_value = ic50_box_value.get()
	docking_score_key = docking_score_keys[protein_index]
	ic50_guesses[docking_score_key][mol_index] = ic50_value
	start_conc = float(start_conc_box.get())
	dilution_factor = float(dilution_factor_box.get())
	n_dilutions = int(n_dilutions_box.get())
	concs = np.array([start_conc * dilution_factor**i \
		for i in range(n_dilutions)])
	min_response = float(min_response_box.get())
	max_response = float(max_response_box.get())
	num_repeats = int(num_repeats_box.get())
	noise_sd = float(noise_sd_box.get())
	os.makedirs("assay_data", exist_ok=True)
	if not write_all_button_value.get():
		mol_names = (file_names[mol_index],)
		ic50s = (ic50_value,)
	else:
		mol_names = file_names
		ic50s = ic50_guesses[docking_score_key]
	for mol_name, ic50_value in zip(mol_names, ic50s):
		if function_type == "sigmoid":
			response = sigmoid_function(concs, ic50_value, min_response,
				max_response)
		elif function_type == "line":
			response = line_function(concs, ic50_value)
		clean_responses = np.tile(response, (num_repeats, 1))
		noisy_responses = clean_responses + np.random.normal(loc=0.0,
			scale=noise_sd, size=clean_responses.shape)
		noisy_responses = np.transpose(noisy_responses)
		protein_name = docking_score_key.rstrip(" score").replace(' ', '_')
		output_filename = mol_name + '_' + protein_name + "_assay.csv"
		timestr = datetime.now().strftime("%d-%m-%H%M") + "_"
		output_filename = timestr + output_filename
		output_filename = os.path.join("assay_data", output_filename)
		print("Writing assay data to file %s..." % output_filename)
		output_file = open(output_filename, "w")
		output_line = "conc," + ",".join("repeat_" + str(i + 1) \
			for i in range(num_repeats)) + "\n"
		output_file.write(output_line)
		for i in range(n_dilutions):
			output_line = str(round_sigfigs(concs[i], 3)) + "," + \
				",".join(str(round_sigfigs(noisy_responses[i,j], 3)) \
				for j in range(num_repeats)) + "\n"
			output_file.write(output_line)
		output_file.close()

def flash_write_all_warning():
	"""
	Shows the warning about using potentially guessed IC50 values when writing
	the assay data for all compounds.
	"""
	global write_all_button_value, write_all_warn_label
	if write_all_button_value.get():
		write_all_warn_label.configure(text="Warning: The current (maybe "\
			"guessed)\nIC50s will be used for all compounds", fg="red")
	else:
		write_all_warn_label.configure(text="\n")

def build_simulator_window(main_window, mols, file_names, properties,
docking_score_keys):
	"""
	Builds the assay simulation window 
	Args:
		main_window: tk.Tk window object, the main window from which this
			simulator will be spawned.
		mols: list of RDKit Mol objects, the molecules from which the properties
			were calculated.
		file_names: list of str, the names of the files for each molecule.
		properties: dict, the computed properties of all of the molecules.
		docking_score_keys: tuple of str, the keys for the properties dictionary
			corresponding to the calculated docking scores.
	"""
	global simulator_frame
	simulator_window = tk.Toplevel(main_window)
	simulator_window.resizable(False, False)
	simulator_window.title("Assay Simulator")
	simulator_frame = tk.Frame(simulator_window)
	simulator_frame.grid(row=0, column=0)

	global mol_index
	global protein_index
	mol_index = 0
	protein_index = 0
	docking_score_key = docking_score_keys[protein_index]
	protein_program = docking_score_key.rstrip(" score")

	global protein_name_label, molecule_name_label
	tk.Label(simulator_frame, text="Protein & Program:").grid(row=0, column=0,
		sticky="nsew")
	protein_name_label = tk.Label(simulator_frame, text=protein_program)
	protein_name_label.grid(row=0, column=1, sticky="nsew")
	tk.Button(simulator_frame, text="Prev Protein",
		command=lambda: set_prev_protein(file_names, properties,
		docking_score_keys)).grid(row=1, column=0, sticky="nsew")
	tk.Button(simulator_frame, text="Next Protein",
		command=lambda: set_next_protein(file_names, properties,
		docking_score_keys)).grid(row=1, column=1, sticky="nsew")
	tk.Label(simulator_frame, text="Current Molecule:").grid(row=2, column=0,
		sticky="nsew")
	molecule_name_label = tk.Label(simulator_frame, text=file_names[mol_index])
	molecule_name_label.grid(row=2, column=1, sticky="nsew")
	tk.Button(simulator_frame, text="Prev Molecule",
		command=lambda: set_prev_molecule(mols, file_names, properties,
		docking_score_keys)).grid(row=3, column=0, sticky="nsew")
	tk.Button(simulator_frame, text="Next Molecule",
		command=lambda: set_next_molecule(mols, file_names, properties,
		docking_score_keys)).grid(row=3, column=1, sticky="nsew")

	global molecule_image, molecule_image_label, molecule_score_label
	global reference_image, reference_image_label, reference_score_label
	flat_mol = Chem.MolFromSmiles(Chem.MolToSmiles(mols[mol_index]))
	molecule_image = ImageTk.PhotoImage(Draw.MolToImage(flat_mol,
		size=(100,100)))
	molecule_image_label = tk.Label(simulator_frame, image=molecule_image,
		text=None)
	molecule_image_label.grid(row=4, column=0, sticky="nsew")
	molecule_score_label = tk.Label(simulator_frame,
		text=file_names[mol_index] + '\n' + docking_score_key + ":\n" + \
		str(properties[docking_score_key][mol_index]))
	molecule_score_label.grid(row=4, column=1, sticky="nsew")
	if REFERENCE_MOLS.get(docking_score_key) is not None:
		reference_image = ImageTk.PhotoImage(Draw.MolToImage(
			REFERENCE_MOLS[docking_score_key], size=(100,100)))
		reference_image_label = tk.Label(simulator_frame, image=reference_image,
			text=None)
		reference_image_label.grid(row=5, column=0, sticky="nsew")
		reference_score_label = tk.Label(simulator_frame,
			text="Reference\n" + docking_score_key + ":\n" + \
			str(REFERENCE_SCORES[docking_score_key]))
	else:
		reference_image = ImageTk.PhotoImage(Image.new("RGBA", (100, 100),
			(0, 0, 0, 0)))
		reference_image_label.configure(image=reference_image, text=None)
		reference_score_label.configure(text="No reference\ncompound for\n" + \
			docking_score_key.rstrip(" score"))
	reference_score_label.grid(row=5, column=1, sticky="nsew")

	tk.Label(simulator_frame, text="IC50 / M:").grid(row=6, column=0,
		sticky="nsew")
	tk.Label(simulator_frame, text="Start conc / M:").grid(row=7, column=0,
		sticky="nsew")
	tk.Label(simulator_frame, text="Dilution factor:").grid(row=8, column=0,
		sticky="nsew")
	tk.Label(simulator_frame, text="Dilutions:").grid(row=9, column=0,
		sticky="nsew")
	tk.Label(simulator_frame, text="Min response:").grid(row=10, column=0,
		sticky="nsew")
	tk.Label(simulator_frame, text="Max response:").grid(row=11, column=0,
		sticky="nsew")
	tk.Label(simulator_frame, text="Num repeats:").grid(row=12, column=0,
		sticky="nsew")
	tk.Label(simulator_frame, text="Noise SD:").grid(row=13, column=0,
		sticky="nsew")

	global ic50_guesses, original_ic50_guesses
	ic50_guesses = dict()
	for docking_key in docking_score_keys:
		ic50_guesses[docking_key] = list()
		if REFERENCE_SCORES.get(docking_key) is not None:
			reference_score = REFERENCE_SCORES[docking_key]
		else:
			reference_score = None
		for name, docking_score in zip(file_names, properties[docking_key]):
			# With an unavailable docking score, assume the score is one-fifth
			# that of the reference
			if docking_score is None:
				print("WARNING: No docking score available for compound '%s'. "\
					"Assuming its docking score is one-fifth that of the "\
					"reference." % name)
				docking_score = 0.2 * reference_score
			if reference_score is not None:
				ic50_guess = 5e-7**(docking_score / reference_score)
				ic50_guess = round_sigfigs(ic50_guess, 3)
				ic50_guesses[docking_key].append(ic50_guess)
			else:
				ic50_guesses[docking_key].append(5e-7)
	original_ic50_guesses = copy.deepcopy(ic50_guesses)

	global ic50_box, start_conc_box, dilution_factor_box, n_dilutions_box
	global min_response_box, max_response_box, num_repeats_box, noise_sd_box
	global ic50_box_value
	ic50_frame = tk.Frame(simulator_frame)
	ic50_frame.grid(row=6, column=1, sticky="nsew")
	ic50_box_value = tk.DoubleVar()
	ic50_box = tk.Entry(ic50_frame, textvariable=ic50_box_value, width=10)
	ic50_box_value.set(ic50_guesses[docking_score_key][0])
	ic50_box.bind("<Return>",
		lambda event: update_simulation_plot(docking_score_keys))
	ic50_box.bind("<KP_Enter>",
		lambda event: update_simulation_plot(docking_score_keys))
	ic50_box.grid(row=0, column=0, sticky="nsew")
	tk.Button(ic50_frame, text=">", height=1, width=1,
		command=lambda: update_simulation_plot(docking_score_keys)).grid(row=0,
		column=1, sticky="nsew")
	start_conc_frame = tk.Frame(simulator_frame)
	start_conc_frame.grid(row=7, column=1, sticky="nsew")
	start_conc_box = tk.Spinbox(start_conc_frame, from_=0.001, to=0.1,
		increment=0.001, repeatdelay=200, repeatinterval=50, width=10)
	start_conc_box.bind("<Return>",
		lambda event: update_simulation_plot(docking_score_keys))
	start_conc_box.bind("<KP_Enter>",
		lambda event: update_simulation_plot(docking_score_keys))
	start_conc_box.delete(0, "end")
	start_conc_box.insert(0, 0.005)
	start_conc_box.grid(row=0, column=0, sticky="nsew")
	tk.Button(start_conc_frame, text=">", height=1, width=1,
		command=lambda: update_simulation_plot(docking_score_keys)).grid(row=0,
		column=1, sticky="nsew")
	dilution_factor_frame = tk.Frame(simulator_frame)
	dilution_factor_frame.grid(row=8, column=1, sticky="nsew")
	dilution_factor_box = tk.Spinbox(dilution_factor_frame, from_=0.1, to=0.9,
		increment=0.1, repeatdelay=200, repeatinterval=50, width=10)
	dilution_factor_box.bind("<Return>",
		lambda event: update_simulation_plot(docking_score_keys))
	dilution_factor_box.bind("<KP_Enter>",
		lambda event: update_simulation_plot(docking_score_keys))
	dilution_factor_box.delete(0, "end")
	dilution_factor_box.insert(0, 0.2)
	dilution_factor_box.grid(row=0, column=0, sticky="nsew")
	tk.Button(dilution_factor_frame, text=">", height=1, width=1,
		command=lambda: update_simulation_plot(docking_score_keys)).grid(row=0,
		column=1, sticky="nsew")
	n_dilutions_frame = tk.Frame(simulator_frame)
	n_dilutions_frame.grid(row=9, column=1, sticky="nsew")
	n_dilutions_box = tk.Spinbox(n_dilutions_frame, from_=1, to=20,
		increment=1, repeatdelay=200, repeatinterval=50, width=10)
	n_dilutions_box.bind("<Return>",
		lambda event: update_simulation_plot(docking_score_keys))
	n_dilutions_box.bind("<KP_Enter>",
		lambda event: update_simulation_plot(docking_score_keys))
	n_dilutions_box.delete(0, "end")
	n_dilutions_box.insert(0, 8)
	n_dilutions_box.grid(row=0, column=0, sticky="nsew")
	tk.Button(n_dilutions_frame, text=">", height=1, width=1,
		command=lambda: update_simulation_plot(docking_score_keys)).grid(row=0,
		column=1, sticky="nsew")
	min_response_frame = tk.Frame(simulator_frame)
	min_response_frame.grid(row=10, column=1, sticky="nsew")
	min_response_box = tk.Spinbox(min_response_frame, from_=0, to=200,
		increment=10, repeatdelay=200, repeatinterval=50, width=10)
	min_response_box.bind("<Return>",
		lambda event: update_simulation_plot(docking_score_keys))
	min_response_box.bind("<KP_Enter>",
		lambda event: update_simulation_plot(docking_score_keys))
	min_response_box.delete(0, "end")
	min_response_box.insert(0, 50)
	min_response_box.grid(row=0, column=0, sticky="nsew")
	tk.Button(min_response_frame, text=">", height=1, width=1,
		command=lambda: update_simulation_plot(docking_score_keys)).grid(row=0,
		column=1, sticky="nsew")
	max_response_frame = tk.Frame(simulator_frame)
	max_response_frame.grid(row=11, column=1, sticky="nsew")
	max_response_box = tk.Spinbox(max_response_frame, from_=10, to=300,
		increment=10, repeatdelay=200, repeatinterval=50, width=10)
	max_response_box.bind("<Return>",
		lambda event: update_simulation_plot(docking_score_keys))
	max_response_box.bind("<KP_Enter>",
		lambda event: update_simulation_plot(docking_score_keys))
	max_response_box.delete(0, "end")
	max_response_box.insert(0, 200)
	max_response_box.grid(row=0, column=0, sticky="nsew")
	tk.Button(max_response_frame, text=">", height=1, width=1,
		command=lambda: update_simulation_plot(docking_score_keys)).grid(row=0,
		column=1, sticky="nsew")
	num_repeats_frame = tk.Frame(simulator_frame)
	num_repeats_frame.grid(row=12, column=1, sticky="nsew")
	num_repeats_box = tk.Spinbox(num_repeats_frame, from_=1, to=20,
		increment=1, repeatdelay=200, repeatinterval=50, width=10)
	num_repeats_box.bind("<Return>",
		lambda event: update_simulation_plot(docking_score_keys))
	num_repeats_box.bind("<KP_Enter>",
		lambda event: update_simulation_plot(docking_score_keys))
	num_repeats_box.delete(0, "end")
	num_repeats_box.insert(0, 3)
	num_repeats_box.grid(row=0, column=0, sticky="nsew")
	tk.Button(num_repeats_frame, text=">", height=1, width=1,
		command=lambda: update_simulation_plot(docking_score_keys)).grid(row=0,
		column=1, sticky="nsew")
	noise_sd_frame = tk.Frame(simulator_frame)
	noise_sd_frame.grid(row=13, column=1, sticky="nsew")
	noise_sd_box = tk.Spinbox(noise_sd_frame, from_=1, to=49,
		increment=2, repeatdelay=200, repeatinterval=50, width=10)
	noise_sd_box.bind("<Return>",
		lambda event: update_simulation_plot(docking_score_keys))
	noise_sd_box.bind("<KP_Enter>",
		lambda event: update_simulation_plot(docking_score_keys))
	noise_sd_box.delete(0, "end")
	noise_sd_box.insert(0, 15)
	noise_sd_box.grid(row=0, column=0, sticky="nsew")
	tk.Button(noise_sd_frame, text=">", height=1, width=1,
		command=lambda: update_simulation_plot(docking_score_keys)).grid(row=0,
		column=1, sticky="nsew")

	reset_ic50_button = tk.Button(simulator_frame, text="Reset IC50s",
		command=lambda: reset_ic50_values(docking_score_keys)).grid(row=14,
		column=0, sticky="nsew")
	reset_variables_button = tk.Button(simulator_frame, text="Reset Variables",
		command=lambda: reset_all_variables(docking_score_keys)).grid(row=14,
		column=1, sticky="nsew")

	global function_type
	function_type = "sigmoid"
	tk.Label(simulator_frame, text="Mode:").grid(row=0, column=2, sticky="w")
	sigmoid_button = tk.Button(simulator_frame, text="Sigmoid",
		command=lambda: set_function_type("sigmoid", docking_score_keys))
	sigmoid_button.grid(row=0, column=3, sticky="w")
	line_button = tk.Button(simulator_frame, text="Line",
		command=lambda: set_function_type("line", docking_score_keys))
	line_button.grid(row=0, column=4, sticky="w")

	global assay_plot_canvas, assay_plot_ax
	fig = Figure(figsize=(5, 4), dpi=100)
	assay_plot_ax = fig.add_subplot()
	assay_plot_canvas = FigureCanvasTkAgg(fig, master=simulator_frame)
	update_simulation_plot(docking_score_keys)
	toolbar = NavigationToolbar2Tk(assay_plot_canvas, simulator_frame,
		pack_toolbar=False)
	toolbar.update()
	assay_plot_canvas.get_tk_widget().grid(row=1, column=2, rowspan=12,
		columnspan=20, sticky="nsew")
	toolbar.grid(row=13, column=2, columnspan=20, sticky="nsew")

	global write_all_button_value, write_all_warn_label
	write_button = tk.Button(simulator_frame, text="Write Data",
		command=lambda: write_assay_data(file_names, docking_score_keys))
	write_button.grid(row=14, column=2, sticky="nsew")
	write_all_button_value = tk.BooleanVar(value=False)
	write_all_button = ttk.Checkbutton(simulator_frame, text="All compounds",
		variable=write_all_button_value, command=flash_write_all_warning)
	write_all_button.grid(row=14, column=3, sticky="nsew")
	write_all_warn_label = tk.Label(simulator_frame, text="\n")
	write_all_warn_label.grid(row=14, column=4, columnspan=18, sticky="nsew")


if __name__ == "__main__":
	window = tk.Tk()
	window.title("Main")
	button = tk.Button(text="Click",
		command=lambda: build_simulator_window(window, 
			[Chem.MolFromSmiles("ClC=CC=O"),
				Chem.MolFromSmiles("c1ccccc1C(O)=O"),
				Chem.MolFromSmiles("c1cccc2ccccc12")],
			["molecule1", "molecule2", "naphalene"],
			{"5IF3 Vina score": [-0.45, -0.71, -7.7],
				"1IEP GOLD score": [50.22, 1.36, 117.0],
				"Custom GOLD score": [44.02, 12.22, 75.73]},
			("5IF3 Vina score", "1IEP GOLD score", "Custom GOLD score")))
	button.grid(row=0, column=0)
	window.mainloop()
