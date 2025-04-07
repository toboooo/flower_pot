"""A simple script that plots the optimised HARDBOILED-EGG ellipses with the HIA
and BBB data sets superimposed on top of them."""
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
sys.path.append("..")
from utils.mol_list import get_mol_list
from reboil_egg import read_from_file, get_egg_features

hia_smiles, hia_labels = read_from_file("hia_mols.csv")
bbb_smiles, bbb_labels = read_from_file("bbb_mols.csv")
hia_mols, hia_err_msgs = get_mol_list(hia_smiles)
bbb_mols, bbb_err_msgs = get_mol_list(bbb_smiles)
hia_features = get_egg_features(hia_mols)
bbb_features = get_egg_features(bbb_mols)
hia_colors = ["blue", "green"]
bbb_colors = ["magenta", "orange"]

ev_hia = np.array([-2.64497678, 3.62300812, 141.27928864, 1.11967626, 144.22556707])

c = np.sqrt((ev_hia[0] - ev_hia[2])**2 + (ev_hia[1] - ev_hia[3])**2) / 2
a = ev_hia[4] / 2
b = np.sqrt(a**2 - c**2)
x = (ev_hia[0] + ev_hia[2]) / 2
y = (ev_hia[1] + ev_hia[3]) / 2
theta = np.arctan((ev_hia[1] - ev_hia[3]) / (ev_hia[0] - ev_hia[2]))
ellipse_hia = Ellipse((x, y), 2*a, 2*b, angle=np.degrees(theta), fill=False, color="black")
filled_ellipse_hia = Ellipse((x, y), 2*a, 2*b, angle=np.degrees(theta), fill=True, color="white")

ev_bbb = np.array([-2.48156312, 3.20928199, 79.25676677, 4.3657658, 81.98390281])

c = np.sqrt((ev_bbb[0] - ev_bbb[2])**2 + (ev_bbb[1] - ev_bbb[3])**2) / 2
a = ev_bbb[4] / 2
b = np.sqrt(a**2 - c**2)
x = (ev_bbb[0] + ev_bbb[2]) / 2
y = (ev_bbb[1] + ev_bbb[3]) / 2
theta = np.arctan((ev_bbb[1] - ev_bbb[3]) / (ev_bbb[0] - ev_bbb[2]))
ellipse_bbb = Ellipse((x, y), 2*a, 2*b, angle=np.degrees(theta), fill=False, color="black")
filled_ellipse_bbb = Ellipse((x, y), 2*a, 2*b, angle=np.degrees(theta), fill=True, color="yellow")

plt.grid()
plt.scatter(hia_features[:,0], hia_features[:,1], c=hia_labels, cmap=matplotlib.colors.ListedColormap(hia_colors), alpha=0.5)
plt.gca().add_patch(ellipse_hia)
plt.xlabel("TPSA")
plt.ylabel("LogP")
plt.show()

plt.grid()
plt.scatter(bbb_features[:,0], bbb_features[:,1], c=bbb_labels, cmap=matplotlib.colors.ListedColormap(bbb_colors), alpha=0.5)
plt.gca().add_patch(ellipse_bbb)
plt.xlabel("TPSA")
plt.ylabel("LogP")
plt.show()

plt.gca().add_patch(filled_ellipse_hia)
plt.gca().add_patch(filled_ellipse_bbb)
plt.scatter(hia_features[:,0], hia_features[:,1], c=hia_labels, cmap=matplotlib.colors.ListedColormap(hia_colors), alpha=0.5)
plt.scatter(bbb_features[:,0], bbb_features[:,1], c=bbb_labels, cmap=matplotlib.colors.ListedColormap(bbb_colors), alpha=0.8)
plt.xlim((-10, 210))
plt.ylim((-3.5, 8.0))
plt.xlabel("TPSA")
plt.ylabel("LogP")
plt.gca().set_facecolor("grey")
plt.show()
