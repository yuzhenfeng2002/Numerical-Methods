import numpy as np


def lspoly(X, Y, d):
    """
    Least squares fit
    :param X: a 1*N numpy array contains abscissas
    :param Y: a 1*N numpy array contains ordinates
    :param d: degree of the polynomial
    :return: the coefficients of least-squares polynomial
    """
    n = len(X)
    F = np.zeros((n, d + 1))
    F[:, 0:1] = np.ones((n, 1))
    F[:, 1] = X.T
    for i in range(2, d + 1):
        F[:, i] = F[:, i-1] * X.T
    C = np.linalg.inv(F.T @ F) @ F.T @ Y.T
    return C


def csfit(X, Y, dx0, dxn):
    """
    Clamped cubic spline
    :param X: a 1*N numpy array contains abscissas
    :param Y: a 1*N numpy array contains ordinates
    :param dx0: S'(a)
    :param dxn: S'(b)
    :return: the coefficient list for N-1 cubic spline polynomials
    """
    n = len(X)
    A = np.zeros((n-2, n-2))
    H = np.diff(X)
    D = np.diff(Y) / H
    V = 6 * np.diff(D)
    V[0] = V[0] - 3 * (D[0] - dx0)
    V[-1] = V[-1] - 3 * (dxn - D[-1])
    V = V.T
    A[0][0] = 3 * H[0] / 2 + 2 * H[1]
    A[0][1] = H[1]
    A[-1][-2] = H[-2]
    A[-1][-1] = 2 * H[-2] + 3 * H[-1] / 2
    for k in range(1, n-3):
        A[k][k-1] = H[k-1]
        A[k][k] = 2 * (H[k-1] + H[k])
        A[k][k+1] = H[k]
    M = np.zeros(n)
    M[1:-1] = np.linalg.inv(A) @ V.T
    M[0] = 3 * (D[0] - dx0) / H[0] - M[1] / 2
    M[-1] = 3 * (dxn - D[-1]) / H[-1] - M[-2] / 2
    S = np.zeros((n-1, 4))
    for k in range(n-1):
        S[k][0] = Y[k]
        S[k][1] = D[k] - H[k] * (2 * M[k] + M[k+1]) / 6
        S[k][2] = M[k] / 2
        S[k][3] = (M[k+1] - M[k]) / 6 / H[k]
    return S


def tpcoeff(X, Y, M):
    """
    Fourier series
    :param X: a 1*N numpy array contains abscissas
    :param Y: a 1*N numpy array contains abscissas
    :param M: degree of the trigonometric polynomial
    :return:
    """
    n = len(X) - 1
    A = np.zeros(M + 1)
    B = np.zeros(M + 1)
    Yends = (Y[0] + Y[-1]) / 2
    Y[0] = Y[-1] = Yends
    A[0] = np.sum(Y)
    for i in range(M):
        A[i + 1] = np.cos((i + 1) * X) @ (Y.T)
        B[i + 1] = np.sin((i + 1) * X) @ (Y.T)
    A = A * 2 / n
    B = B * 2 / n
    A[0] = A[0] / 2
    return A, B
