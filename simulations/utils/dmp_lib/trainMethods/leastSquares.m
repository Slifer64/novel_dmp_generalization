%% Performs kernel based LS (Least Squares)
%  Performs least squares learning on the input data and returns the learned weights.
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
function w = leastSquares(Psi, X, Fd, zero_tol)
  
if (nargin < 4), zero_tol = 0.0; end

N_kernels = size(Psi,1);
w_dim = size(X,1);
n_data = size(Psi,2);

% w = zeros(N_kernels, w_dim);

% normalize Psi
Psi = Psi ./ (repmat(sum(Psi,1),size(Psi,1),1) + zero_tol);

H = zeros(N_kernels*w_dim, n_data);
k = 1;
for i=1:w_dim
    H(k:k+N_kernels-1,:) = Psi.*repmat(X(i,:),N_kernels,1);
    k = k+N_kernels;
end

w = reshape((Fd/H)', N_kernels, w_dim);

end
