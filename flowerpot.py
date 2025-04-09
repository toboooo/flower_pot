"""The main script that implements the tkinter GUI and handles the calculations
of the physical properties the user is interested in."""
import os
import re
import csv
import tkinter as tk
import tkinter.ttk as ttk
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import askdirectory
import numpy as np
from PIL import Image, ImageTk
from rdkit import Chem
from rdkit.Chem import Descriptors
from openpyxl import load_workbook, Workbook
from utils.mol_list import get_mol_list
from boiled_egg.print_egg import print_eggs
from substructure.filters import check_filters
from log_models.fragment_model import FragmentLogModel
from log_models.log_ml_model import MLPWrapper
from docking.docking import perform_vina_docking, perform_gold_docking, estimate_ic50, custom_gold_docking, custom_vina_docking

FRAME_X_PADDING = 1
FRAME_Y_PADDING = 1
MIN_SIZE_X = 1100
MIN_SIZE_Y = 530
TEXTBOX_WIDTH = 40
TEXTBOX_HEIGHT = 30
IMAGE_FRAME_X = 200
IMAGE_FRAME_Y = 220

LOGS_MEAN = np.array([59.85272991, 102.0204244])
LOGS_STD = np.array([57.07686877, 61.51980086])
LOGD_MEAN = np.array([78.00117738, 158.2493462])
LOGD_STD = np.array([33.83845809, 47.38073137])

def get_filename(filename_var, filename_entry):
	"""Implements the 'Browse' button to select an input file."""
	filename = askopenfilename()
	filename_var.set(filename)
	filename_entry.xview_moveto(1)

def get_directory(directory_var, directory_entry):
	"""Implements the 'Browse' button to select a directory."""
	dir_name = askdirectory()
	directory_var.set(dir_name)
	directory_entry.xview_moveto(1)

def activate_print_egg_button():
	"""Ensures that if the 'print egg' option is selected, the permeation
	calculation button will also be selected."""
	if print_egg_value.get():
		egg_value.set(True)
		print_egg_warning.config(text="Warning: Drawing lots of EGG images can take a long time.\nUse small numbers of molecules for better performance.", fg="orange")
	else:
		print_egg_warning.config(text="\n")

def activate_egg_button():
	"""Ensures that if the 'Permeation' option is deselected, the 'print egg'
	option will also be deselected."""
	if not egg_value.get():
		print_egg_value.set(False)
		print_egg_warning.config(text="\n")

def activate_ic50_button():
	"""Ensures that if the 'IC50' option is selected, the docking calculation
	button will also be selected."""
	if ic50_value.get():
		docking_value.set(True)
		docking_warning.config(text="Warning: Molecular docking is computationally intensive.\nOnly use for a small number of molecules at a time.", fg="red")

def activate_docking_button():
	"""Ensures that if the docking calculation option is deselected, the IC50
	estimation button will also be deselected, and also highlights a warning
	message to the user about the computational cost of docking."""
	if docking_value.get():
		docking_warning.config(text="Warning: Molecular docking is computationally intensive.\nOnly use for a small number of molecules at a time.", fg="red")
	else:
		docking_warning.config(text="\n")
		ic50_value.set(False)

def backward():
	"""Implements the functioning of the previous button for the pop-up image window."""
	global mol_no
	global image_list
	global file_names
	global smiles_heading
	global image_labels
	mol_no = (mol_no - 1) % len(file_names)
	smiles_heading.configure(text=file_names[mol_no])
	if len(image_list[mol_no]) > 0:
		picture_frame.grid_propagate(True)
		for i in range(len(image_labels)):
			if i < len(image_list[mol_no]):
				image_labels[i].configure(image=image_list[mol_no][i], text=None)
			else:
				image_labels[i].configure(image="", text=None)
		picture_frame.configure(height=picture_frame["height"], width=picture_frame["width"])
	else:
		picture_frame.grid_propagate(False)
		image_labels[0].configure(image="", text="No images for this molecule.")
		picture_frame.configure(width=IMAGE_FRAME_X, height=IMAGE_FRAME_Y)
		for label in image_labels[1:]:
			label.configure(image="", text=None)

def forward():
	"""Implements the functioning of the next button for the pop-up image window."""
	global mol_no
	global image_list
	global file_names
	global smiles_heading
	global image_labels
	mol_no = (mol_no + 1) % len(file_names)
	smiles_heading.configure(text=file_names[mol_no])
	if len(image_list[mol_no]) > 0:
		picture_frame.grid_propagate(True)
		for i in range(len(image_labels)):
			if i < len(image_list[mol_no]):
				image_labels[i].configure(image=image_list[mol_no][i], text=None)
			else:
				image_labels[i].configure(image="", text=None)
		picture_frame.configure(height=picture_frame["height"], width=picture_frame["width"])
	else:
		picture_frame.grid_propagate(False)
		image_labels[0].configure(image="", text="No images for this molecule.")
		picture_frame.configure(width=IMAGE_FRAME_X, height=IMAGE_FRAME_Y)
		for label in image_labels[1:]:
			label.configure(image="", text=None)

def calc_properties():
	"""
	Handles all inputs and outputs to the program, as well as the calculations
	of all physical properties that selected using the checkbutto options.

	This function's first step is to read in any SMILES strings in the
	program's input text box, provided the button that states this input
	box is to be ignored is not checked.

	If provided with an existing spreadsheet file name, the function will then
	attempt to find a column titled 'smiles' from which to read the molecule
	SMILES strings. Failing to find a column heading named 'smiles' (case-
	insensitive), the function will then look for the first column in the first
	row of the spreadsheet that contains a valid SMILES string and will use this
	column to read in the SMILES strings.

	Then, the function will perform all calculations that were requested by the
	user using the check boxes. First LogD, LogP, LogS and TPSA will be
	calculated for all of the input molecules that had valid SMILES strings. The
	user requested model type (Crippen linear model or neural network) will be
	used to calculate LogD and LogS.

	Next, the estimation of intestinal absorption and blood-brain barrier
	permeation using the HARDBOILED-EGG model will be performed, the
	substructure filters will be applied to all of the input molecules and then
	a pop-up window will be created that will display the visualisations of
	these results.

	The last optional property calculations that will be performed are 
	molecular docking simulations using either of the programs Autodock Vina or
	GOLD, if the executable binaries can be located. The docking scores may also
	be used in the pedagogical calculations of the IC50 value of a molecule,
	based on the difference between it's binding affinity and that of a
	reference structure. The user also has the option to provide customised
	Vina and GOLD job files, which will be run after the jobs specifying an
	in-built protein target.

	The function will write all of the results to the output text box on the
	right hand side of the GUI. Unless it is requested by using a check button
	that errors from incorrect SMILES strings are ignored, then error messages
	for these strings will also be written to the output text box.

	Finally, the function will write all of the calculated results to an output
	file, if a file name is provided. If requested as a .xlsx file, the results
	will be written to an Excel spreadsheet on a worksheet named 'Sheet1'. If
	the file extension is not one of '.csv' or '.xlsx', the .csv file extension
	will be added to the provided name and results will be written in the csv
	format.
	"""
	# Delete output text to show the job has started
	output_text.configure(state="normal")
	output_text.delete("1.0", "end")
	output_text.update_idletasks()

	smiles_strings = []
	names = []
	# Get input from textbox
	if not ignore_textbox_value.get():
		smiles_input = smiles_text.get("1.0", "end-1c")
		for smiles in smiles_input.splitlines():
			if smiles != "":
				smiles_strings.append(smiles)
				names.append(None)
	# Read input spreadsheet
	input_filename = csv_filename.get()
	if os.path.exists(input_filename):
		print("Attempting to read input from input file '%s'." % input_filename)
		try:
			if input_filename.endswith(".csv"):
				file = open(input_filename, "r")
				csv_reader = csv.reader(file)
				first_line = next(csv_reader)
				smiles_pos = None
				name_pos = None
				# Check if the file has headings and find the smiles and name positions
				lower_first_line = [item.lower() for item in first_line]
				if "smiles" in lower_first_line:
					smiles_pos = lower_first_line.index("smiles")
					if "name" in lower_first_line:
						name_pos = lower_first_line.index("name")
					elif "names" in lower_first_line:
						name_pos = lower_first_line.index("names")
				# Try looking for a valid smiles string instead
				else:
					print("Could not find 'smiles' column heading, looking for valid SMILES string in first row instead...")
					for i, item in enumerate(first_line):
						mol = Chem.MolFromSmiles(item, sanitize=False)
						if mol is not None and item != "":
							smiles_pos = i
							smiles_strings.append(item)
							names.append(None)
							break
					if smiles_pos is None:
						print("ERROR: Failed to find valid SMILES string in first line of file %s." % input_filename)
				if smiles_pos is not None:
					for line in csv_reader:
						if smiles_pos < len(line):
							smiles_strings.append(line[smiles_pos])
						else:
							print("WARNING: Line was shorter than expected position of SMILES string.")
							continue
						if name_pos is not None:
							if name_pos < len(line):
								names.append(line[name_pos])
							else:
								print("WARNING: Line was shorter than expected position of molecule name.")
								names.append(None)
						else:
							names.append(None)
				file.close()
			elif input_filename.endswith(".xlsx"):
				workbook = load_workbook(input_filename, read_only=True)
				worksheet = workbook[workbook.sheetnames[0]]
				smiles_pos = None
				name_pos = None
				rows = list(worksheet.rows)
				# Get first line of spreadsheet, record all cell values as strings
				first_line = [str(cell.value) for cell in rows[0]]
				lower_first_line = [item.lower() for item in first_line]
				# Check if the file has the smiles and names headings and note
				# their positions
				if "smiles" in lower_first_line:
					smiles_pos = lower_first_line.index("smiles")
					if "name" in lower_first_line:
						name_pos = lower_first_line.index("name")
					elif "names" in lower_first_line:
						name_pos = lower_first_line.index("names")
				# Try looking for a valid smiles string in the first line instead
				else:
					print("Could not find 'smiles' column heading, looking for valid SMILES string in first row instead...")
					for i, item in enumerate(first_line):
						mol = Chem.MolFromSmiles(item, sanitize=False)
						if mol is not None:
							smiles_pos = i
							smiles_strings.append(item)
							names.append(None)
							break
					if smiles_pos is None:
						print("ERROR: Failed to find valid SMILES string in first line of file %s." % input_filename)
				if smiles_pos is not None:
					for row in rows[1:]:
						if smiles_pos < len(row):
							smiles_strings.append(row[smiles_pos].value)
						else:
							print("WARNING: Line was shorter than expected position of SMILES string.")
							continue
						if name_pos is not None:
							if name_pos < len(row):
								names.append(row[name_pos].value)
							else:
								print("WARNING: Line was shorter than expected position of molecule name.")
								names.append(None)
						else:
							names.append(None)
				workbook.close()
			else:
				print("ERROR: Unsupported input file type, please use either a .csv file or a .xlsx spreadsheet.")
		except PermissionError:
			print("ERROR: Input file '%s' could not be read: Permission denied." % input_filename)
	elif not os.path.exists(input_filename) and input_filename != "":
		print("ERROR: Input file '%s' not found." % input_filename)

	# Process input smiles strings and associate each molecule with a name
	mols, errs, names = get_mol_list(smiles_strings, names=names)
	if len(mols) > 0:
		# Get list of strings for filenames
		global file_names
		file_names = []
		i = 0
		for name in names:
			if name is None or name == "" or name.isspace():
				file_names.append("molecule" + str(i+1))
				i += 1
			else:
				# Make sure there are no operating system invalid characters in name
				file_names.append(re.sub("[\\:/]", "", name))
		# Ensure that the same name does not appear twice
		for i, name in enumerate(file_names):
			n = 1
			for j, other_name in enumerate(file_names[i+1:], start=i+1):
				if other_name == name:
					other_name += "-" + str(n)
					file_names[j] = other_name
					n += 1

		# Initialise storage of properties
		properties = dict()
		heading = "name,smiles,"
		key_list = []

		# Numerical properties
		if any(wanted.get() for prop, wanted in numerical_properties.items()):
			if log_predictor_value.get() == 0:
				logs_model = FragmentLogModel("log_models/Crippen.txt")
				logs_model.load_coef_file("log_models/sol_coefs.npy")
				logd_model = FragmentLogModel("log_models/Crippen.txt")
				logd_model.load_coef_file("log_models/lipo_coefs.npy")
			elif log_predictor_value.get() == 1:
				logs_weights = np.load("log_models/logs_network_weights.npz")
				logs_biases = np.load("log_models/logs_network_biases.npz")
				logd_weights = np.load("log_models/logd_network_weights.npz")
				logd_biases = np.load("log_models/logd_network_biases.npz")
				logs_model = MLPWrapper(logs_weights, logs_biases, LOGS_MEAN, LOGS_STD, activation="relu")
				logd_model = MLPWrapper(logd_weights, logd_biases, LOGD_MEAN, LOGD_STD, activation="relu")
			for prop in sorted(numerical_properties.keys()):
				if numerical_properties[prop].get():
					heading += prop + ","
					key_list.append(prop)
					properties[prop] = []
			for prop in sorted(properties.keys()):
				if prop == "LogD":
					predictions = logd_model.predict(mols)
					properties[prop] = list(float(pred) for pred in predictions)
				elif prop == "LogS":
					predictions = logs_model.predict(mols)
					properties[prop] = list(float(pred) for pred in predictions)
				elif prop == "LogP":
					for mol in mols:
						properties[prop].append(Descriptors.MolLogP(mol))
				elif prop == "tPSA":
					for mol in mols:
						properties[prop].append(Descriptors.TPSA(mol))

		# HARDBOILED-EGG
		if egg_value.get():
			egg_image_files, hia_perms, bbb_perms = print_eggs(mols, file_names, print_egg_value.get())
			heading += "Gastrointestinal permeation,Blood-brain barrier permeation,"
			key_list.append("HIA")
			key_list.append("BBB")
			properties["HIA"] = hia_perms
			properties["BBB"] = bbb_perms

		# Substructure filters
		filter_image_files = []
		filter_names = []
		for filt in sorted(filter_values.keys()):
			if filter_values[filt].get():
				image_files, hit_names, hit_counts = check_filters(mols, file_names, filt, print_images=True)
				filter_image_files.append(image_files)
				filter_names.append(hit_names)
				heading += filt + " hit count,"
				key_list.append(filt)
				properties[filt] = hit_counts

		# Docking score and IC50 estimation
		if docking_value.get():
			# Dock with in-built proteins first
			for protein in sorted(protein_selection.keys()):
				if protein_selection[protein].get():
					if program_value.get() == 0:
						docking_scores = perform_gold_docking(mols, file_names, protein, gold_install_path.get())
					elif program_value.get() == 1:
						docking_scores = perform_vina_docking(mols, file_names, protein)
					if docking_scores is not None:
						heading += protein + " docking score,"
						key_list.append(protein + "_docking_score")
						properties[protein + "_docking_score"] = docking_scores
						if ic50_value.get():
							estimated_ic50s = estimate_ic50(docking_scores, protein, "gold" if program_value.get() == 0 else "vina")
							heading += protein + " IC50,"
							key_list.append(protein + "_ic50")
							properties[protein + "_ic50"] = estimated_ic50s
					else:
						print("ERROR: Could not perform docking score calculations for protein: %s." % protein)
			# Custom GOLD docking
			if select_custom_gold.get():
				if os.path.exists(conf_filename.get()):
					custom_gold_scores = custom_gold_docking(mols, file_names, conf_filename.get(), gold_install_path.get())
					if custom_gold_scores is not None:
						heading += "Custom GOLD score,"
						key_list.append("gold_docking_score")
						properties["gold_docking_score"] = custom_gold_scores
					else:
						print("ERROR: Could not perform docking score calculations with GOLD.")
				else:
					print("WARNING: gold.conf file '%s' was not found. Cannot perform custom GOLD docking." % conf_filename.get())
			# Custom Vina docking
			if select_custom_vina.get():
				if os.path.exists(pdbqt_filename.get()) and os.path.exists(config_filename.get()):
					custom_vina_scores = custom_vina_docking(mols, file_names, pdbqt_filename.get(), config_filename.get())
					if custom_vina_scores is not None:
						heading += "Custom Vina score,"
						key_list.append("vina_docking_score")
						properties["vina_docking_score"] = custom_vina_scores
					else:
						print("ERROR: Could not perform docking score calculations with Autodock Vina.")
				else:
					if not os.path.exists(pdbqt_filename.get()):
						print("WARNING: PBDQT file '%s' was not found. Cannot perform custom Vina docking." % pdbqt_filename.get())
					if not os.path.exists(config_filename.get()):
						print("WARNING: Vina config file '%s' was not found. Cannot perform custom Vina docking." % config_filename.get())

		# Pop-up image window
		if print_egg_value.get() or any(filter_values[filt].get() for filt in filter_values.keys()):
			global mol_no
			mol_no = 0
			global image_list
			image_list = []
			for i in range(len(file_names)):
				mol_image_list = []
				if print_egg_value.get():
					mol_image_list.append(ImageTk.PhotoImage(Image.open(egg_image_files[i])))
				if any(filter_values[filt].get() for filt in filter_values.keys()):
					for filter_image_list in filter_image_files:
						if filter_image_list[i] != None:
							mol_image_list.append(ImageTk.PhotoImage(Image.open(filter_image_list[i])))
				image_list.append(mol_image_list)

			image_window = tk.Toplevel(window)
			image_window.resizable(False, False)
			image_window.title("EGG and Filter Images")
			global picture_frame
			picture_frame = tk.Frame(image_window)
			picture_frame.grid(row=0, column=0, columnspan=2, rowspan=2, sticky="nsew")
			global smiles_heading
			smiles_heading = tk.Label(picture_frame, text=file_names[mol_no])
			smiles_heading.grid(row=0, column=0, columnspan=2, sticky="nsew")
			global image_labels
			image_labels = []
			for i in range(len(filter_values.keys()) + 1):
				label = tk.Label(picture_frame)
				label.grid(row=1, column=i, sticky="nsew")
				image_labels.append(label)
			if len(image_list[mol_no]) > 0:
				picture_frame.grid_propagate(True)
				for i in range(len(image_labels)):
					if i < len(image_list[mol_no]):
						image_labels[i].configure(image=image_list[mol_no][i], text=None)
					else:
						image_labels[i].configure(image="", text=None)
				picture_frame.configure(height=picture_frame["height"], width=picture_frame["width"])
			else:
				picture_frame.grid_propagate(False)
				image_labels[0].configure(image="", text="No images for this molecule.")
				picture_frame.configure(width=IMAGE_FRAME_X, height=IMAGE_FRAME_Y)
				for label in image_labels[1:]:
					label.configure(image="", text=None)
			button_frame = tk.Frame(image_window)
			button_frame.grid(row=2, column=0, columnspan=2, sticky="nws")
			backward_button = tk.Button(button_frame, text="Previous", command=backward)
			backward_button.grid(row=2, column=0, sticky="ew")
			forward_button = tk.Button(button_frame, text="Next", command=forward)
			forward_button.grid(row=2, column=1, sticky="ew")

		# Write output text
		output_text.insert(tk.END, heading + "\n")
		write_to_file = True
		is_excel = False
		out_foldername = output_foldername.get()
		out_filename = output_filename.get()
		try:
			if out_filename.endswith(".xlsx"):
				is_excel = True
			elif out_filename == "":
				out_filename = "output.csv"
				print("Writing data to output.csv...")
			elif not out_filename.endswith(".csv"):
				out_filename += ".csv"
				print("Output file name does not end with either '.xlsx' or '.csv', defaulting to a CSV file '%s'." % out_filename)
			out_filename = os.path.join(out_foldername, out_filename)
			if is_excel:
				output_workbook = Workbook(write_only=True)
				# output_file is actually a spreadsheet sheet if writing
				# to an excel file
				output_file = output_workbook.create_sheet("Sheet1")
			else:
				output_file = open(out_filename, "w")
		except FileNotFoundError:
			print("ERROR: Could not open file %s for output. Perhaps check the directory path is spelt correctly.\nContinuing without writing to file..." % out_filename)
			write_to_file = False
		except PermissionError:
			print("ERROR: File '%s' could not be opened for writing output: Permission denied.\nContinuing without writing to file..." % out_filename)
			write_to_file = False
		if write_to_file:
			if is_excel:
				output_file.append(heading.rstrip(",").split(","))
			else:
				output_file.write(heading + "\n")
		m = 0
		for smiles, err in zip(smiles_strings, errs):
			if err is None:
				line = [file_names[m], smiles]
				for key in key_list:
					prop = properties[key][m]
					if type(prop) == float:
						line.append("%.2f" % prop)
					else:
						line.append(str(prop))
				output_text.insert(tk.END, ",".join(line) + "\n")
				output_text.see(tk.END)
				if write_to_file:
					if is_excel:
						output_file.append(line)
					else:
						output_file.write(",".join(line) + "\n")
				m += 1
			elif type(err) == str and not quiet_errors_value.get():
				output_text.insert(tk.END, err + "\n")
				output_text.see(tk.END)
				if write_to_file:
					if is_excel:
						output_file.append([err])
					else:
						output_file.write(err + "\n")
		if write_to_file:
			if is_excel:
				try:
					output_workbook.save(out_filename)
				except FileNotFoundError:
					print("ERROR: File '%s' could not be saved: File not found." % out_filename)
				except PermissionError:
					print("ERROR: Input file '%s' could not be read: Permission denied." % input_filename)
				print("Writing output to file '%s'." % out_filename)
				output_workbook.close()
			else:
				print("Writing output to file '%s'." % out_filename)
				output_file.close()
		output_text.configure(state="disabled")

	else:
		print("ERROR: No valid SMILES strings were found.")


if __name__ == "__main__":
	window = tk.Tk()
	window.title("Flower Pot")
	window.minsize(MIN_SIZE_X, MIN_SIZE_Y)

	# SMILES input textbox on left of GUI
	smiles_frame = tk.Frame()
	smiles_frame.grid(row=0, column=0, sticky="nsew", padx=FRAME_X_PADDING, pady=FRAME_Y_PADDING)
	smiles_label = ttk.Label(smiles_frame, text="SMILES Input:")
	smiles_label.grid(row=0, column=0, sticky="w")
	smiles_text = tk.Text(smiles_frame, width=TEXTBOX_WIDTH, height=TEXTBOX_HEIGHT)
	smiles_text.grid(row=1, column=0, sticky="nsew")
	smiles_frame.grid_columnconfigure(0, weight=1)
	smiles_frame.grid_rowconfigure(1, weight=1)

	# Configuration options in middle of GUI
	option_frame = tk.Frame()
	option_frame.grid(row=0, column=1, sticky="nw", padx=FRAME_X_PADDING, pady=FRAME_Y_PADDING)
	# Input and output file options
	input_label = ttk.Label(option_frame, text="Input File Path:")
	input_label.grid(row=0, column=0, sticky="nw")
	csv_filename = tk.StringVar()
	csv_entry = ttk.Entry(option_frame, textvariable=csv_filename, state="disabled")
	csv_entry.grid(row=1, column=0, sticky="nw")
	browse_button = ttk.Button(option_frame, text="Browse", command=lambda: get_filename(csv_filename, csv_entry))
	browse_button.grid(row=1, column=1, sticky="nw")
	ignore_textbox_value = tk.BooleanVar()
	ignore_textbox_button = ttk.Checkbutton(option_frame, variable=ignore_textbox_value, text="Only use input from file")
	ignore_textbox_button.grid(row=2, column=0, sticky="nw")
	quiet_errors_value = tk.BooleanVar()
	quiet_errors_button = ttk.Checkbutton(option_frame, variable=quiet_errors_value, text="Quieten invalid SMILES errors")
	quiet_errors_button.grid(row=3, column=0, sticky="nw")
	output_folder_label = ttk.Label(option_frame, text="Output Folder:")
	output_folder_label.grid(row=4, column=0, sticky="nw")
	output_file_label = ttk.Label(option_frame, text="Output File Name:")
	output_file_label.grid(row=4, column=2, sticky="nw")
	output_foldername = tk.StringVar()
	output_folder_entry = ttk.Entry(option_frame, textvariable=output_foldername, state="disabled")
	output_folder_entry.grid(row=5, column=0, sticky="nw")
	output_browse_button = ttk.Button(option_frame, text="Select Folder", command=lambda: get_directory(output_foldername, output_folder_entry))
	output_browse_button.grid(row=5, column=1, sticky="nw")
	output_filename = tk.StringVar()
	output_entry = ttk.Entry(option_frame, textvariable=output_filename)
	output_entry.grid(row=5, column=2, sticky="nw")

	# Numerical properties section
	property_label = ttk.Label(option_frame, text="Numerical Properties:")
	property_label.grid(row=6, column=0, sticky="nw")
	numerical_properties = {"LogP": tk.BooleanVar(), "LogD": tk.BooleanVar(), "LogS": tk.BooleanVar(), "tPSA": tk.BooleanVar()}
	for i, num_prop in enumerate(("LogD", "LogS", "LogP", "tPSA")):
		checkbox = ttk.Checkbutton(option_frame, variable=numerical_properties[num_prop], text=num_prop)
		checkbox.grid(row=7, column=i, sticky="w")
	log_predictor_label = ttk.Label(option_frame, text="LogD/LogS Predictor:")
	log_predictor_label.grid(row=8, column=0, sticky="nw")
	log_predictor_value = tk.IntVar()
	crippen_button = ttk.Checkbutton(option_frame, variable=log_predictor_value,onvalue=0, offvalue=0, text="Crippen Linear Model")
	ml_button = ttk.Checkbutton(option_frame, variable=log_predictor_value, onvalue=1, offvalue=1, text="Machine Learning")
	crippen_button.grid(row=9, column=0, sticky="nw")
	ml_button.grid(row=9, column=1, sticky="nw")

	# HARDBOILED-EGG section
	egg_label = ttk.Label(option_frame, text="HARDBOILED-EGG:")
	egg_label.grid(row=10, column=0, sticky="w")
	egg_value = tk.BooleanVar()
	print_egg_value = tk.BooleanVar()
	egg_button = ttk.Checkbutton(option_frame, variable=egg_value, text="Permeation", command=activate_egg_button)
	egg_button.grid(row=11, column=0, sticky="w")
	print_egg_button = ttk.Checkbutton(option_frame, variable=print_egg_value, text="Print EGG", command=activate_print_egg_button)
	print_egg_button.grid(row=11, column=1, sticky="w")
	print_egg_warning = tk.Label(option_frame, text="\n")
	print_egg_warning.grid(row=12, column=0, columnspan=3, sticky="nsew")

	# Substructure filters section
	filter_label = ttk.Label(option_frame, text="Substructure Filters:")
	filter_label.grid(row=13, column=0, sticky="w")
	filter_values = {"PAINS_A": tk.BooleanVar(), "PAINS_B": tk.BooleanVar(), "PAINS_C": tk.BooleanVar(), "BRENK": tk.BooleanVar(), "NIH": tk.BooleanVar()}
	for i, filt in enumerate(sorted(list(filter_values.keys()))):
		checkbox = ttk.Checkbutton(option_frame, variable=filter_values[filt], text=filt)
		checkbox.grid(row=14, column=i, sticky="w")

	# Docking section
	docking_label = ttk.Label(option_frame, text="Molecular Docking:")
	docking_label.grid(row=15, column=0, sticky="w")
	docking_value = tk.BooleanVar()
	docking_button = ttk.Checkbutton(option_frame, variable=docking_value, text="Calculate Docking Score", command=activate_docking_button)
	docking_button.grid(row=16, column=0, sticky="w")
	ic50_value = tk.BooleanVar()
	ic50_button = ttk.Checkbutton(option_frame, variable=ic50_value, text="Estimate IC50", command=activate_ic50_button)
	ic50_button.grid(row=16, column=1, sticky="w")
	docking_warning = tk.Label(option_frame, text="\n")
	docking_warning.grid(row=17, column=0, columnspan=3, sticky="nsew")
	program_label = ttk.Label(option_frame, text="Docking Program:")
	program_label.grid(row=18, column=0, sticky="w")
	program_value = tk.IntVar()
	gold_button = ttk.Checkbutton(option_frame, variable=program_value, onvalue=0, offvalue=0, text="GOLD")
	vina_button = ttk.Checkbutton(option_frame, variable=program_value, onvalue=1, offvalue=1, text="AutoDock Vina")
	gold_button.grid(row=19, column=0, sticky="w")
	vina_button.grid(row=19, column=1, sticky="w")
	gold_box_label = tk.Label(option_frame, text="GOLD Installation Folder:")
	gold_box_label.grid(row=20, column=0, sticky="w")
	gold_install_path = tk.StringVar()
	gold_install_box = ttk.Entry(option_frame, textvariable=gold_install_path)
	gold_install_box.grid(row=20, column=1, sticky="w")
	gold_install_button = ttk.Button(option_frame, text="Browse", command=lambda: get_directory(gold_install_path, gold_install_box))
	gold_install_button.grid(row=20, column=2, sticky="w")
	proteins_label = tk.Label(option_frame, text="Target Proteins:")
	proteins_label.grid(row=21, column=0, sticky="w")
	protein_selection = {"26S_proteasome": tk.BooleanVar(), "BTHalpha": tk.BooleanVar(), "CF-IIbeta": tk.BooleanVar(), "CHK1": tk.BooleanVar(), "CYP17a": tk.BooleanVar()}
	for i, protein in enumerate(sorted(list(protein_selection.keys()))):
		protein_button = ttk.Checkbutton(option_frame, variable=protein_selection[protein], text=protein)
		protein_button.grid(row=22, column=i, sticky="w")

	# Custom docking section
	select_custom_gold = tk.BooleanVar()
	custom_gold_button = ttk.Checkbutton(option_frame, text="Custom GOLD Docking", variable=select_custom_gold)
	custom_gold_button.grid(row=23, column=0, sticky="nw")
	conf_file_label = tk.Label(option_frame, text="gold.conf File:")
	conf_file_label.grid(row=24, column=0, sticky="nw")
	conf_filename = tk.StringVar()
	conf_file_box = ttk.Entry(option_frame, textvariable=conf_filename)
	conf_file_box.grid(row=24, column=1, sticky="nw")
	conf_file_button = ttk.Button(option_frame, text="Browse", command=lambda: get_filename(conf_filename, conf_file_box))
	conf_file_button.grid(row=24, column=2, sticky="nw")
	select_custom_vina = tk.BooleanVar()
	custom_vina_button = ttk.Checkbutton(option_frame, text="Custom Vina Docking", variable=select_custom_vina)
	custom_vina_button.grid(row=25, column=0, sticky="nw")
	pdbqt_file_label = tk.Label(option_frame, text="Protein PDBQT File:")
	pdbqt_file_label.grid(row=26, column=0, sticky="nw")
	pdbqt_filename = tk.StringVar()
	pdbqt_file_box = ttk.Entry(option_frame, textvariable=pdbqt_filename)
	pdbqt_file_box.grid(row=26, column=1, sticky="nw")
	pdbqt_file_button = ttk.Button(option_frame, text="Browse", command=lambda: get_filename(pdbqt_filename, pdbqt_file_box))
	pdbqt_file_button.grid(row=26, column=2, sticky="nw")
	config_file_label = tk.Label(option_frame, text="Vina Config File:")
	config_file_label.grid(row=27, column=0, sticky="nw")
	config_filename = tk.StringVar()
	config_file_box = ttk.Entry(option_frame, textvariable=config_filename)
	config_file_box.grid(row=27, column=1, sticky="nw")
	config_file_button = ttk.Button(option_frame, text="Browse", command=lambda: get_filename(config_filename, config_file_box))
	config_file_button.grid(row=27, column=2, sticky="nw")

	go_button = ttk.Button(option_frame, text="Go", command=calc_properties)
	go_button.grid(row=28, column=0, sticky="nw")

	# Output textbox on right of GUI
	output_frame = tk.Frame()
	output_frame.grid(row=0, column=2, sticky="nsew", padx=FRAME_X_PADDING, pady=FRAME_Y_PADDING)
	output_label = tk.Label(output_frame, text="Calculation Output:")
	output_label.grid(row=0, column=0, sticky="nw")
	output_text = tk.Text(output_frame, state="disabled", width=TEXTBOX_WIDTH, height=TEXTBOX_HEIGHT)
	output_text.grid(row=1, column=0, sticky="nsew")
	output_frame.grid_columnconfigure(0, weight=1)
	output_frame.grid_rowconfigure(1, weight=1)

	window.columnconfigure((0, 1), weight=1)
	window.columnconfigure(2, weight=3)
	window.rowconfigure(0, weight=1)

	window.mainloop()
