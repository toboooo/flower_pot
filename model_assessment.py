import numpy as np
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.metrics import balanced_accuracy_score, mean_absolute_error
from sklearn.preprocessing import StandardScaler
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.svm import SVC, SVR
from mlpdropout import MLPDropoutClassifier, MLPDropoutRegressor
from features import get_descriptors, get_fingerprints, get_3d_descriptors

def hyperparameter_search(model, X_train, y_train, parameter_ranges, classification=False, n_jobs=1):
	if classification:
		grid_search = GridSearchCV(model, parameter_ranges, scoring="balanced_accuracy", cv=5, refit=True, n_jobs=n_jobs)
	else:
		grid_search = GridSearchCV(model, parameter_ranges, scoring="neg_mean_absolute_error", cv=5, refit=True, n_jobs=n_jobs)
	grid_search.fit(X_train, y_train)
	return grid_search

def model_test(model, X_train, X_test, y_train, y_test, parameter_ranges, classification=False):
	search_result = hyperparameter_search(model, X_train, y_train, parameter_ranges, classication=classification)
	predictions = search_result.best_estimator_.predict(X_test)
	if classification:
		score = balanaced_accuracy_score(predictions, y_test)
	else:
		score = -mean_absolute_error(predictions, y_test)
	return score, search_result.best_params_

def select_model(X_train, X_test, y_train, y_test, classification=False):
	if classification:
		models = [RandomForestClassifier, SVC, MLPDropoutClassifier]
		parameter_ranges_list = [{"n_estimators": [50, 100, 200], "max_depth": [5, 10, 20], "min_samples_leaf": [5, 20, 100], "max_features": [0.2, 0.5, 0.8]}, {"kernel": ["linear", "rbf", "sigmoid"], "C": [0.001, 0.1, 1.0, 10.0], "gamma": [0.001, 0.1, 1.0, 10.0]}, {"hidden_layer_sizes": [(256, 256), (512, 512), (1024, 1024), (2048, 2048), (256, 256, 256, 256), (512, 512, 512, 512)], "alpha": [0.0, 1e-6, 1e-5, 1e-4], "learning_rate": [0.01, 0.05, 0.1], "dropout_rate": [None, 0.2, 0.5]}]
	else:
		models = [RandomForestRegressor, SVR, MLPDropoutRegressor]
		parameter_ranges_list = [{"n_estimators": [50, 100, 200], "max_depth": [5, 10, 20], "min_samples_leaf": [5, 20, 100], "max_features": [0.2, 0.5, 0.8]}, {"kernel": ["linear", "rbf", "sigmoid"], "C": [0.001, 0.1, 1.0, 10.0], "gamma": [0.001, 0.1, 1.0, 10.0], "epsilon": [0.001, 0.1, 1.0]}, {"hidden_layer_sizes": [(256, 256), (512, 512), (1024, 1024), (2048, 2048), (256, 256, 256, 256), (512, 512, 512, 512)], "alpha": [0.0, 1e-6, 1e-5, 1e-4], "learning_rate": [0.01, 0.05, 0.1], "dropout_rate": [None, 0.2, 0.5]}]
	model_names = ["RF", "SVM", "DNN"]
	best_score = 0 if classification else -float("inf")
	best_parameters = None
	best_model = None
	for i, model in enumerate(models):
		test_score, parameters = model_test(model, X_train, X_test, y_train, y_test, parameter_ranges_list[i], classification=classification)
		if test_score > best_score:
			best_model = model_names[i]
			best_parameters = parameters
	return best_model, best_score, best_parameters

def select_features(mols, y, parameter_ranges, classification=False):
	descriptors = get_descriptors(mols)
	fingerprints = get_fingerprints(mols)
	descriptors_3d = get_3d_descriptors(mols)
	descriptor_sets = [descriptors, fingerprints, descriptors_3d]
	for i, descriptor_set in enumerate(descriptor_sets):
		if i == 1:
			continue
		descriptor_transformer = ColumnTransformer([("standard_scaler", StandardScaler(), np.arange(descriptor_set.shape[1]))], remainder="passthrough")
		descriptor_sets[i] = descriptor_transformer.fit_transform(descriptor_set)
	descriptor_sets += [np.hstack((descriptor_sets[0], descriptor_sets[1])), np.hstack((descriptor_sets[0], descriptor_sets[2])), np.hstack((descriptor_sets[1], descriptor_sets[2])), np.hstack((descriptor_sets[0], descriptor_sets[1], descriptor_sets[2]))]
	descriptor_names = ["descriptors", "fingerprints", "3d", "descriptors+fingerprints", "descriptors+3d", "fingerprints+3d", "descriptors+fingerprints+3d"]
	for X, name in zip(descriptor_sets, descriptor_names):
		X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1, shuffle=False)
		best_model, best_score, best_hyperparameters = select_model(X_train, X_test, y_train, y_test, classification=classification)
		print(names)
		print("Best model: %s" % best_model)
		print("Test score: %f" % best_score)
		print("Hyperparameters:")
		print(best_parameters)
		print("\n")
