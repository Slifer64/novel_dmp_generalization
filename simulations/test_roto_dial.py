import numpy as np

def fnAR(x):
    '''
    Accelerated rotation of vector x to the direction of axis x_0.
    '''
    x = np.ndarray.flatten(x)
    n = np.size(x)
    R = np.eye(n)
    step = 1
    while (step < n):
        A = np.eye(n)
        it = 0
        while(it < n - step):
            r2 = x[it] * x[it] + x[it + step] * x[it + step]
            if (r2 > 0):
                r = np.sqrt(r2)
                pcos = x[it] / r
                psin = - x[it + step] / r

                # Base 2-dimensional rotation
                A[it, it] = pcos
                A[it, it + step] = - psin
                A[it + step, it] = psin
                A[it + step, it + step] = pcos
            it = it + 2 * step
        step = step * 2
        x = np.dot(A, x)
        R = np.dot(A, R)
    return R

def roto_dilatation(x0, x1):
    '''
    Return the roto-dilatation matrix which maps x0 to x1.
    '''

    # Write the inputs as 1D array
    x0 = np.ndarray.flatten(x0)
    x1 = np.ndarray.flatten(x1)
    if (np.size(x0) != np.size(x1)):
        raise ValueError(
            'The two vectors must have the same number of components')

    # Extract the norms
    norm_0 = np.linalg.norm(x0, 2)
    norm_1 = np.linalg.norm(x1, 2)

    # Normalize the two vectors
    x0_norm = x0 / norm_0
    x1_norm = x1 / norm_1
    Mx0 = fnAR(x0_norm)
    Mx1 = fnAR(x1_norm)

    print('Mx0 = ', Mx0)
    print('Mx1 = ', Mx1)

    M = np.dot(np.transpose(Mx1), Mx0)
    M *= norm_1 / norm_0
    return M


if __name__ == '__main__':

    x0 = np.array([0.5, 0.6, 0.2])
    x1 = np.array([0.2, 0.8, 0.6])

    M = roto_dilatation(x0, x1)

    x1_hat = np.dot(M, x0)

    print(x1)
    print(x1_hat)
