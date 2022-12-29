%% Locally Weighted Regression
%  Performs locally weighted regression learning on the input data and returns the learned weights.
%  N denotes the number of data points.
%  M denotes the dimensionality of the input data.
%  K denotes the number of kernels.
%  The kernels are allocated in the data point space and provide as output the activation value of the
%  kernel for a specific data point.
%  The data point can be anything, e.g. a 1-D data point representing time (or a substitute variable for time),
%  or a 3-D data point representing Cartesian position, or Cartesian force, or it can be any other M-D datapoint
%  representing anything you like.
%  @param[in] Psi: K X N matrix, where the k-th row contains the k-th kernel function values for all data points.
%  @param[in] X: M X N matrix, where the j-th column corresponds to the j-th data point.
%  @param[in] Fd: 1 x N row vector, where the j-th value corresponds to the desired output for the j-th data point.
%  @param[in] zero_tol: Tollerance value to avoid divisions by zero.
%  @param[out] w: K x N matrix with the learned weights.
%
function w = LWR(Psi, X, Fd, zero_tol)
  
  if (nargin < 4), zero_tol = 0.0; end

  N_kernels = size(Psi,1);
  w_dim = size(X,1);
  
  w = zeros(N_kernels, w_dim);

  for k=1:N_kernels
      X_Psi = X .* repmat(Psi(k,:), size(X,1), 1);
      w(k,:) = ((X_Psi*X' + zero_tol) \ X_Psi*Fd')';
      
%       w(k,:) = (pinv(X_Psi*X' + zero_tol) * X_Psi*Fd')';
  end

end


