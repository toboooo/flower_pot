from sklearn.datasets import make_classification, make_regression
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier, MLPRegressor
from mlpdropout import MLPDropoutClassifier, MLPDropoutRegressor

X, y = make_classification(n_samples=100, random_state=1)
X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, random_state=1)
mlp = MLPClassifier(random_state=1, max_iter=300).fit(X_train, y_train)
print(mlp.score(X_test, y_test))
mlp = MLPDropoutClassifier(dropout=0.2, random_state=1, max_iter=300).fit(X_train, y_train)
print(mlp.score(X_test, y_test))

X, y = make_regression(n_samples=200, random_state=1)
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=1)
mlp = MLPRegressor(random_state=1, max_iter=10000).fit(X_train, y_train)
print(mlp.score(X_test, y_test))
mlp = MLPDropoutRegressor(dropout=0.05, random_state=1, max_iter=10000).fit(X_train, y_train)
print(mlp.score(X_test, y_test))
