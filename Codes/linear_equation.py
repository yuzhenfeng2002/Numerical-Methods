def backsub(U: list, y: list):
    """
    To solve linear equation like Ux = y
    :param U: an N*N matrix
    :param y: an N*1 matrix
    :return: x, an N*1 matrix
    """
    n = len(y)
    x = [[y[n - 1][0] / U[n - 1][n - 1]]]
    for k in range(n - 2, -1, -1):
        ux = 0
        for i in range(k + 1, n):
            ux += x[n - 1 - i][0] * U[k][i]
        x.append([(y[k][0] - ux) / U[k][k]])
    x.reverse()
    return x


def forsub(L: list, y: list):
    """
    To solve linear equation like Lx = y
    :param L: an N*N matrix
    :param y: an N*1 matrix
    :return: x, an N*1 matrix
    """
    n = len(y)
    x = [[y[0][0] / L[0][0]]]
    for k in range(1, n):
        lx = 0
        for i in range(k):
            lx += x[i][0] * L[k][i]
        x.append([(y[k][0] - lx) / L[k][k]])
    return x


def uptrbk(A: list, b: list):
    """
    To solve linear equation like Ax = b
    :param A: an N*N matrix
    :param b: an N*1 matrix
    :return: x, an N*1 matrix
    """
    n = len(b)
    for p in range(n - 1):
        max_pivot = abs(A[p][p]); max_index = p
        for i in range(p + 1, n):
            if abs(A[i][p]) > max_pivot:
                max_pivot = abs(A[i][p])
                max_index = i
        if max_pivot == 0:
            print("A was singular!")
            return
        if max_index != p:
            A[max_index], A[p] = A[p], A[max_index]
            b[max_index], b[p] = b[p], b[max_index]
        for i in range(p + 1, n):
            factor = A[i][p] / A[p][p]
            A[i][p] = 0
            b[i][0] = b[i][0] - factor * b[p][0]
            for j in range(p + 1, n):
                A[i][j] = A[i][j] - factor * A[p][j]
    return backsub(A, b)


def lufact(A, b):
    """
    To solve linear equation like Ax = b
    :param A: an N*N matrix
    :param b: an N*1 matrix
    :return: x, an N*1 matrix
    """
    n = len(b)

    # LU decomposition
    R = list(range(n))
    for p in range(n - 1):
        max_pivot = abs(A[p][p]); max_index = p
        for i in range(p + 1, n):
            if abs(A[i][p]) > max_pivot:
                max_pivot = abs(A[i][p])
                max_index = i
        if max_pivot == 0:
            print("A was singular!")
            return
        if max_index != p:
            A[max_index], A[p] = A[p], A[max_index]
            R[max_index], R[p] = R[p], R[max_index]
        for i in range(p + 1, n):
            factor = A[i][p] / A[p][p]
            A[i][p] = factor
            for j in range(p + 1, n):
                A[i][j] = A[i][j] - factor * A[p][j]

    # Solve L(Ux) = b
    y = [[b[R[0]][0]]]
    for k in range(1, n):
        ly = 0
        for i in range(k):
            ly += y[i][0] * A[k][i]
        y.append([(b[R[k]][0] - ly)])

    # Solve Ux = y
    x = [[y[n - 1][0] / A[n - 1][n - 1]]]
    for k in range(n - 2, -1, -1):
        ux = 0
        for i in range(k + 1, n):
            ux += x[n - 1 - i][0] * A[k][i]
        x.append([(y[k][0] - ux) / A[k][k]])
    x.reverse()
    return x


def jacobi(A, b, P, delta, max1):
    """
    To solve linear equation like Ax = b
    :param A: an N*N matrix
    :param b: an N*1 matrix
    :param P: an N*1 matrix
    :param delta: the tolerance
    :param max1: the maximum number of the iterations
    :return: (the number of iterations, the approximation to the solution, the absolute error, the relative error)
    """
    n = len(b)
    for k in range(max1):
        X = []
        err = 0
        norm_X = 0
        for j in range(n):
            ap = 0
            for i in range(n):
                if i != j:
                    ap += A[j][i] * P[i][0]
            x = (b[j][0] - ap) / A[j][j]
            norm_X += x
            err += abs(x - P[j][0])
            X.append([x])
        try:
            relerr = err / (norm_X)
        except ZeroDivisionError:
            relerr = delta + 1
        P = X
        if err < delta or relerr < delta:
            break
    return k + 1, X, err, relerr


def gseid(A, b, P, delta, max1):
    """
    To solve linear equation like Ax = b
    :param A: an N*N matrix
    :param b: an N*1 matrix
    :param P: an N*1 matrix
    :param delta: the tolerance
    :param max1: the maximum number of the iterations
    :return: (the number of iterations, the approximation to the solution, the absolute error, the relative error)
    """
    n = len(b)
    for k in range(max1):
        X = []
        err = 0
        norm_X = 0
        P0 = P
        for j in range(n):
            ap = 0
            for i in range(n):
                if i != j:
                    ap += A[j][i] * P0[i][0]
            x = (b[j][0] - ap) / A[j][j]
            norm_X += x
            err += abs(x - P[j][0])
            X.append([x])
            P0[j][0] = x
        try:
            relerr = err / (norm_X)
        except ZeroDivisionError:
            relerr = delta + 1
        P = X
        if err < delta or relerr < delta:
            break
    return k + 1, X, err, relerr
