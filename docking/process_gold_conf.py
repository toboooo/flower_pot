import sys
import shutil

overwrite_fields = {"SAVE OPTIONS": False, "WRITE OPTIONS": False,  "DATA FILES": False}
overwrite_instructions = {"SAVE OPTIONS": "save_score_in_file = 1\nsave_protein_torsions = 1\nclean_up_option delete_empty_directories\nclean_up_option delete_redundant_log_files\nclean_up_option save_top_n_solutions 1\nclean_up_option delete_all_initialised_ligands\n\n",
"WRITE OPTIONS": "write_options = NO_LOG_FILES NO_LINK_FILES NO_RNK_FILES NO_GOLD_SOLN_LIGAND_MOL2_FILES NO_GOLD_PROTEIN_MOL2_FILE NO_LGFNAME_FILE NO_PLP_MOL2_FILES NO_PID_FILE NO_SEED_LOG_FILE NO_GOLD_ERR_FILE NO_FIT_PTS_FILES NO_ASP_MOL2_FILES NO_GOLD_LOG_FILE\n\n",
"DATA FILES": "ligand_data_file ligand.sdf 10\nparam_file = DEFAULT\nset_ligand_atom_types = 1\nset_protein_atom_types = 0\ndirectory = .\ntordist_file = DEFAULT\nmake_subdirs = 0\nsave_lone_pairs = 1\nread_fitpts = 0\n\n"}

def process_gold_conf(gold_conf_file):
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
				conf_file.write(overwrite_instructions[possible_field])
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
			conf_file.write(overwrite_instructions[field])
	conf_file.close()

if __name__ == "__main__":
	process_gold_conf(sys.argv[1])
