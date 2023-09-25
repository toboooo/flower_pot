#https://datascience.stackexchange.com/questions/117082/how-can-i-implement-dropout-in-scikit-learn
import warnings
from sklearn.neural_network import MLPClassifier, MLPRegressor
from sklearn.neural_network._stochastic_optimizers import AdamOptimizer
from sklearn.neural_network._base import ACTIVATIONS, DERIVATIVES, LOSS_FUNCTIONS
from sklearn.utils import shuffle, gen_batches, check_random_state, _safe_indexing
from sklearn.utils.extmath import safe_sparse_dot
from sklearn.exceptions import ConvergenceWarning

class MLPDropoutClassifier(MLPClassifier):
	def __init__(self, hidden_layer_sizes=(100,), activation="relu", *,
		solver="adam", alpha=0.0001, batch_size="auto",
		learning_rate="constant", learning_rate_init=0.001, power_t=0.5,
		max_iter=200, shuffle=True, random_state=None, tol=1e-4, verbose=False,
		warm_start=False, momentum=0.9, nesterovs_momentum=True,
		early_stopping=False, validation_fraction=0.1, beta_1=0.9, beta_2=0.999,
		epsilon=1e-8, n_iter_no_change=10, max_fun=15000, dropout=None):
		super().__init__(hidden_layer_sizes=hidden_layer_sizes,
			activation=activation, solver=solver, alpha=alpha, 
			batch_size=batch_size, learning_rate=learning_rate,
			learning_rate_init=learning_rate_init, power_t=power_t,
			max_iter=max_iter, shuffle=shuffle, random_state=random_state,
			tol=tol, verbose=verbose, warm_start=warm_start, momentum=momentum,
			nesterovs_momentum=nesterovs_momentum,
			early_stopping=early_stopping,
			validation_fraction=validation_fraction, beta_1=beta_1,
			beta_2=beta_2, epsilon=epsilon, n_iter_no_change=n_iter_no_change,
			max_fun=max_fun)

	def _fit_stochastic(self, X, y, activations, deltas, coef_grads, intercept_grads, layer_units, incremental):
		params = self.coefs_ + self.intercepts_
		if not incremental or not hasattr(self, "_optimizer"):
			if self.solver == "sgd":
				self._optimizer = SGDOptimizer(params, self.learning_rate_init, self.learning_rate, self.momentum, self.nestervos_momentum, self.power_t)
			elif self.solver == "adam":
				self._optimizer = AdamOptimizer(params, self.learning_rate_init, self.beta_1, self.beta_2, self.epsilon)
		early_stopping = self.early_stopping and not incremental
		if early_stopping:
			should_stratify = is_classifier(self) and self.n_outputs_ == 1
			stratify = y if should_stratify else None
			X, X_val, y, y_val = train_test_split(X, y, random_state=self._random_state, test_size=self.validation_function, stratify=stratify)
			if is_classifier(self):
				y_val = self._label_binarizer.inverse_transform(y_val)
		else:
			X_val = None
			y_val = None
		n_samples = X.shape[0]
		sample_idx = np.arange(n_samples, dtype=int)
		if self.batch_size = "auto":
			batch_size = min(200, n_samples)
		else:
			if self.batch_size < 1 or self.batch_size > n_samples:
				warnings.warn("Got batch_size less 1 or larger than sample size.")
				batch_size = np.clip(self.batch_size, 1, n_samples)
		try:
			for it in range(self.max_iter):
				if self.shuffle:
					sample_idx = shuffle(sample_idx, random_state=self._random_state)
				accumulated_loss = 0.0
				for batch_slice in gen_batches(n_samples, batch_size):
					if self.shuffle:
						X_batch = _safe_indexing(X, sample_idx[batch_slice])
						y_batch = y[sample_idx[batch_size]]
					else:
						X_batch = X[batch_slice]
						y_batch = y[batch_slice]
				activations[0] = X_batch
				batch_loss, coef_grads, intercept_grads = self._backprop(X_batch, y_batch, activations, layer_units, deltas, coef_grads, intercept_grads)
				accumulated_loss += batch_loss * (batch_slice.stop - batch_slice.start)
				grads = coef_grads + intercept_grads
				self._optimizer.update_params(params, grads)
			self.n_iter_ += 1
			self.loss_ = accumulated_loss / X.shape[0]
			self.t_ += n_samples
			self.loss_curve_.append(self.loss_)
			if self.verbose:
				print("Iteration %d, loss = %.8f" % (self.n_iter_, self.loss_))
			self._update_no_improvement_count(early_stopping, X_val, y_val)
			self._optimizer.iteration_ends(self.t_)
			if self._no_improvement_count > self.n_iter_no_change:
				if early_stopping:
					msg = "Validation score did not improve more than tol=%f for %d consecutive epochs." % (self.tol, self.n_iter_no_change)
				else:
					msg = "Training loss did not improve more than tol=%f for %d consecutive epochs." % (self.tol, self.n_iter_no_change)
				is_stopping = self._optimizer.trigger_stopping(msg, self.verbose)
				if is_stopping:
					break
				else:
					self._no_improvemetn_count = 0
			if incremental:
				break
			if self.n_iter_ == self.max_iter:
				warnings.warn("Stochastic Optimizer: Maximum iterations (%d) reached and the optimization hasn't yet converged." % self.max_iter, ConvergenceWarning)
		except KeyboardInterrupt:
			warnings.warn("Training interrupted by keyboard.")
		if early_stopping:
			self.coefs_ = self._best_coefs
			self.intercepts_ = self._best_intercepts

	def _backprop(self, X, y, activations, layer_units, deltas, coef_grads, intercept_grads):
		n_samples = X.shape[0]
		dropout_mask = None
		if self.dropout != None:
			if 0 < self.dropout < 1:
				keep_probability = 1 - self.dropout
				dropout_masks = [np.ones(layer_units[0])]
				for units in layer_units[1:-1]:
					if self.random_state != None:
						layer_mask = (self._random_state.random(units) < keep_probability).astype(int) / keep_probability
					else:
						layer_mask = (np.random.rand(units) < keep_probability).astype(int) / keep_probability
					dropout_masks.append(layer_mask)
			else:
				raise ValueError("Dropout must be between zero and one.")
		activations = self._forward_pass(activations, dropout_masks)
		loss_func_name = self.loss
		if loss_func_name == "log_loss" and self.out_activation_ == "logistic":
			loss_func_name = "binary_log_loss"
		loss = LOSS_FUNCTIONS[loss_func_name](y, activations[-1])
		values = 0
		for s in self.coefs_:
			s = s.ravel()
			values += np.dot(s, s)
		loss += (0.5 * self.alpha) * values / n_samples
		last = self.n_layers_ - 2
		deltas[last] = activations[-1] - y
		self._compute_loss_grad(last, n_samples, activations, deltas, coef_grads, intercept_grads)
		inplace_derivative = DERIVATIVES[self.activation]
		for i in range(self.n_layers_ - 2, 0, -1):
			deltas[i-1] = safe_sparse_dot(delta[i], self.coefs_[i].T)
			inplace_derivative(activations[i], deltas[i-1])
			self._compute_loss_grad(i - 1, n_samples, activations, deltas, coef_grads, intercept_grads)
		if dropout_masks != None:
			for layer in range(len(coef_grads) - 1):
				mask = (~(dropout_masks[layer+1] == 0)).astype(int)
				coef_grads[layer] = coef_grads[layer] * mask[None,:]
				coef_grads[layer+1] = coef_grads[layer+1] * mask.reshape(-1,1)
				intercept_grads[layer] = intercept_grads[layer] * mask
		return loss, coef_grads, intercept_grads

	def _forward_pass(self, activations, dropout_masks=None):
		hidden_activation = ACTIVATIONS[self.activation]
		for i in range(self.n_layers_ - 1):
			activations[i+1] = safe_sparse_dot(activations[i], self.coefs_[i])
			activations[i+1] += self.intercepts_[i]
			if i + 1 != self.n_layers_ - 1:
				hidden_activation(activations[i+1])
			if i + 1 != self.n_layers_ - 1 and dropout_masks != None:
				check1 = activations[i].copy()
				activatons[i+1] = activations[i+1] * dropout_masks[i+1][None,:]
		output_activation = ACTIVATIONS[self.out_activation_]
		output_activation(activation[i+1])
		return activations
