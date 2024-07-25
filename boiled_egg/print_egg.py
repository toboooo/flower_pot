"""Provides the functions that calculate the HARDBOILED-EGG model predictions of
intestinal and blood-brain barrier permeation, as well as drawing the
corresponding plots if requested."""
import math
import pathlib
import multiprocessing
from datetime import datetime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from rdkit.Chem import Draw
from boiled_egg.reboil_egg import get_egg_features

CORE_COUNT_DIV = 2
MIN_MULTI_IMAGES = 50

def make_egg_plots(mols, egg_features, image_names, hia_params, bbb_params):
	"""
	Makes a plot of each molecule's position within the HARDBOILED-EGG model
	along with a drawing of that molecule in the top right corner of the plot.
	Args:
		mols: list of Mol objects, a list of the RDKit molecule objects.
		egg_features: 2d np.array, array containing the TPSA values of each of
			the molecules in the first column and the LogP values in the second
			column.
		image_names: list of str, a list of the file names to which the egg
			images will be saved.
		hia_params: tuple of float, ellipse parameters for the HIA classifying
			ellipse. In order these are: the x- and y-coordinates of the centre,
			the lengths of the major and minor axes and the tilt angle of the
			major axis above the x-axis.
		bbb_params: tuple of float, ellipse parameters for the BBB classifying
			ellipse. In order these are: the x- and y-coordinates of the centre,
			the lengths of the major and minor axes and the tilt angle of the
			major axis above the x-axis.
	"""
	fig, ax = plt.subplots(figsize=(5,4))
	newax = fig.add_axes([0.6, 0.6, 0.4, 0.4], anchor="NE", zorder=1)
	for mol, feat_vec, image_name in zip(mols, egg_features, image_names):
		ax.set_facecolor("grey")
		ax.set_xlabel("TPSA")
		ax.set_ylabel("LogP")
		ax.set_xlim(-10, 200)
		ax.set_ylim(-4, 8)
		hia_ellipse = Ellipse((hia_params[0], hia_params[1]), 2*hia_params[2], 2*hia_params[3], angle=math.degrees(hia_params[4]), fill=True, color="white")
		ax.add_patch(hia_ellipse)
		bbb_ellipse = Ellipse((bbb_params[0], bbb_params[1]), 2*bbb_params[2], 2*bbb_params[3], angle=math.degrees(bbb_params[4]), fill=True, color="yellow")
		ax.add_patch(bbb_ellipse)
		ax.scatter(feat_vec[0], feat_vec[1], s=30.0, c="black", marker="x")
		mol_img = Draw.MolToImage(mol)
		newax.imshow(mol_img)
		newax.axis("off")
		plt.savefig(image_name)
		print("Finished %s!" % image_name)
		ax.clear()
		newax.clear()

def print_eggs(mols, file_names, make_plots=False):
	"""
	Calculates the permeation predictions of the HARDBOILED-EGG model and
	handles the drawing of the egg images. If there are more images to draw than
	the value stored in MIN_MULTI_IMAGES, then half of the available CPU cores
	will be used by the multiprocessing library to speed up the image drawing.
	Args:
		mols: list of Mol objects, a list of RDKit molecule objects.
		file_names: list of str, a list of the names associated with the
		molecules.
		make_plots: optional, bool, whether to draw the model plots.
	Return:
		image_names: list of str, a list of all the names of the egg images.
		hia_perms: list of int, a list of the HIA permeation predictions for
			each molecule (0 indicates non-permeating, 1 indicates permeating).
		bbb_perms: list of int, a list of the BBB permeation predictions for
			each molecule (0 indicates non-permeating, 1 indicates permeating).
	"""
	pathlib.Path("images").mkdir(parents=True, exist_ok=True)
	hia_ev = [-2.64497678, 3.62300812, 141.27928864, 1.11967626, 144.22556707]
	hia_a = hia_ev[4] / 2
	hia_c = math.sqrt((hia_ev[0] - hia_ev[2])**2 + (hia_ev[1] - hia_ev[3])**2) / 2
	hia_b = math.sqrt(hia_a**2 - hia_c**2)
	hia_x = (hia_ev[0] + hia_ev[2]) / 2
	hia_y = (hia_ev[1] + hia_ev[3]) / 2
	hia_theta = math.atan((hia_ev[1] - hia_ev[3]) / (hia_ev[0] - hia_ev[2]))
	bbb_ev = [-2.48156312, 3.20928199, 79.25676677, 4.3657658, 81.98390281]
	bbb_a = bbb_ev[4] / 2
	bbb_c = math.sqrt((bbb_ev[0] - bbb_ev[2])**2 + (bbb_ev[1] - bbb_ev[3])**2) / 2
	bbb_b = math.sqrt(bbb_a**2 - bbb_c**2)
	bbb_x = (bbb_ev[0] + bbb_ev[2]) / 2
	bbb_y = (bbb_ev[1] + bbb_ev[3]) / 2
	bbb_theta = math.atan((bbb_ev[1] - bbb_ev[3]) / (bbb_ev[0] - bbb_ev[2]))
	egg_features = get_egg_features(mols)
	hia_perms = []
	bbb_perms = []
	image_names = []
	for file_name in file_names:
		timecode = datetime.now().strftime("%d-%m-%H%M")
		image_names.append("images/" + timecode + "_" + file_name + "_egg.png")
	for mol, feat_vec in zip(mols, egg_features):
		hia_focus_dists = math.sqrt((hia_ev[0] - feat_vec[0])**2 + (hia_ev[1] - feat_vec[1])**2) + math.sqrt((hia_ev[2] - feat_vec[0])**2 + (hia_ev[3] - feat_vec[1])**2)
		if hia_focus_dists < hia_ev[4]:
			hia_perms.append(1)
		else:
			hia_perms.append(0)
		bbb_focus_dists = math.sqrt((bbb_ev[0] - feat_vec[0])**2 + (bbb_ev[1] - feat_vec[1])**2) + math.sqrt((bbb_ev[2] - feat_vec[0])**2 + (bbb_ev[3] - feat_vec[1])**2)
		if bbb_focus_dists < bbb_ev[4]:
			bbb_perms.append(1)
		else:
			bbb_perms.append(0)
	if make_plots:
		if len(image_names) < MIN_MULTI_IMAGES:
			make_egg_plots(mols, egg_features, image_names, (hia_x, hia_y, hia_a, hia_b, hia_theta), (bbb_x, bbb_y, bbb_a, bbb_b, bbb_theta))
		else:
			n_available_cores = multiprocessing.cpu_count()
			if n_available_cores > 1:
				n_used_cores = n_available_cores // CORE_COUNT_DIV
			else:
				n_used_cores = n_available_cores
			seq_len = math.ceil(len(image_names) / n_used_cores)
			pool = multiprocessing.Pool(n_used_cores)
			processes = [pool.apply_async(make_egg_plots, args=(mols[seq_len*i:seq_len*i+seq_len], egg_features[seq_len*i:seq_len*i+seq_len], image_names[seq_len*i:seq_len*i+seq_len], (hia_x, hia_y, hia_a, hia_b, hia_theta), (bbb_x, bbb_y, bbb_a, bbb_b, bbb_theta))) for i in range(n_used_cores)]
			for p in processes:
				p.get()
			pool.close()
			pool.join()
	return image_names, hia_perms, bbb_perms
