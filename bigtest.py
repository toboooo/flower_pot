import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.model_selection import train_test_split
from sklearn.metrics import matthews_corrcoef, mean_absolute_error
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.svm import SVC, SVR
from mlpdropout import MLPDropoutClassifier, MLPDropoutRegressor
from features import get_mol_list, get_descriptors, get_fingerprints, get_3d_descriptors

# LogD
def read_lipophilicity(filename):
	smiles = []
	logd = []
	with open(filename, "r", newline="") as file:
		csv_reader = csv.DictReader(file, delimiter=",")
		for line in csv_reader:
			smiles.append(line["smiles"])
			logd.append(float(line["exp"]))
	return smiles, np.array(logd)

smiles, logds = read_lipophilicity("Lipophilicity.csv")
mols = get_mol_list(smiles)
#print("Generating 2D descriptors...")
#des_2d = get_descriptors(mols)
#des = get_descriptors(mols)
#print("Generating 3D descriptors...")
#des_3d = get_3d_descriptors(mols)
#des = np.hstack((des_2d, des_3d))
fps = get_fingerprints(mols)

#X_train, X_test, y_train, y_test = train_test_split(des, logds, test_size=0.2, random_state=5)
X_train, X_test, y_train, y_test = train_test_split(fps, logds, test_size=0.2, random_state=5)

#scaler_X = StandardScaler().fit(X_train)
scaler_y = StandardScaler().fit(y_train.reshape(-1,1))
#X_train = scaler_X.transform(X_train)
#X_test = scaler_X.transform(X_test)
y_train = scaler_y.transform(y_train.reshape(-1,1)).ravel()

#rfr = RandomForestRegressor(n_estimators=100, max_depth=15, min_samples_leaf=10, max_features=0.5)
#svr = SVR(kernel="rbf", C=1.0, gamma=1.0, epsilon=0.1)
mlp = MLPDropoutRegressor(hidden_layer_sizes=(512,512,512), alpha=1e-4, dropout=None)

#print("Training RFR...")
#rfr.fit(X_train, y_train)
#pred = rfr.predict(X_test)
#pred = scaler_y.inverse_transform(pred.reshape(-1,1)).ravel()
#print("RFR MAE:", mean_absolute_error(pred, y_test))
#plt.plot(pred, y_test, "o")
#plt.show()

#print("Training SVR...")
#svr.fit(X_train, y_train)
#pred = svr.predict(X_test)
#pred = scaler_y.inverse_transform(pred.reshape(-1,1)).ravel()
#print("SVR MAE:", mean_absolute_error(pred, y_test))
#plt.plot(pred, y_test, "o")
#plt.show()

print("Training MLP...")
mlp.fit(X_train, y_train)
pred = mlp.predict(X_test)
pred = scaler_y.inverse_transform(pred.reshape(-1,1)).ravel()
print("MLP MAE:", mean_absolute_error(pred, y_test))
print("MLP R2:", pearsonr(pred, y_test)[0]**2)
plt.plot(pred, y_test, "o")
plt.show()
