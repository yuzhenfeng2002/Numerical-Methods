from numpy.polynomial import Polynomial as P
import numpy as np


def lagran(X, Y):
    """
    Lagrange approximation
    :param X: a 1*N numpy array contains abscissas
    :param Y: a 1*N numpy array contains ordinates
    :return: (the coefficients of Lagrange polynomial, the coefficients of L_Nk)
    """
    n = len(X)
    L = np.zeros((n, n))
    for k in range(n):
        V = P([1])
        for j in range(n):
            if k != j:
                V = V * P([-X[j], 1]) / (X[k] - X[j])
        V = V.coef
        L[k, :] = V
    C = np.diag(Y) @ L
    C = np.sum(C, axis=0)
    return C, L


def newpoly(X, Y):
    """
    Newton approximation
    :param X: a 1*N numpy array contains abscissas
    :param Y: a 1*N numpy array contains ordinates
    :return: (the coefficients of Newton polynomial, the divided-difference table)
    """
    n = len(X)
    D = np.zeros((n, n))
    D[:, 0] = Y.T
    for j in range(1, n):
        for k in range(j, n):
            D[k, j] = (D[k, j-1] - D[k-1, j-1]) / (X[k] - X[k-j])
    C = P([D[n-1, n-1]])
    for k in range(n-2, -1, -1):
        C = C * P([-X[k], 1])
        m = len(C.coef)
        C = C + P([D[k, k]])
    C = C.coef
    return C, D
