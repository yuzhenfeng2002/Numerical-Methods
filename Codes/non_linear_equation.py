def fixpt(g, p0, tol, max1):
    """
    To solve fixed point problem
    :param g: the iteration function
    :param p0: the initial point
    :param tol: the tolerance of the error
    :param max1: the maximum number of the iterations
    :return: (the number of iterations, the approximation to the fixed point, the sequence pn, the absolute error, the relative error)
    """
    P = [p0]
    n = 0
    for k in range(max1):
        pk = g(P[-1])
        P.append(pk)
        n += 1
        err = abs(P[-1] - P[-2])
        try:
            relerr = abs(P[-1] - P[-2]) / abs(P[-1])
        except ZeroDivisionError:
            relerr = tol + 1
        if err < tol or relerr < tol:
            break
    if k + 1 == max1:
        print("Maximum number of iterations exceeded!")
    return k + 1, P[-1], P, err, relerr


def bisect(f, a, b, delta):
    """
    To find one root of an equation
    :param f: the function
    :param a: the left point
    :param b: the right point
    :param delta: the tolerance
    :return: (the number of iterations, the approximation to the solution x0, the absolute error, the value f(x0))
    """
    if f(a)*f(b) > 0:
        print("Note: f(a)*f(b) > 0")
        return
    # max1 = round((math.log(b - a) - math.log(delta)) / math.log(2))
    k = 0
    while True:
        k += 1
        c = (a + b) / 2
        if f(c) == 0:
            a = c
            b = c
        elif f(a)*f(c) > 0:
            a = c
        else:
            b = c
        if b - a < delta:
            break
    c = (a + b) / 2
    err = abs(b - a)
    return k, c, err, f(c)


def regula(f, a, b, epsilon, max1):
    """
    To find one root of an equation
    :param f: the function
    :param a: the left point
    :param b: the right point
    :param epsilon: the tolerance
    :param max1: the maximum number of the iterations
    :return: (the number of iterations, the approximation to the solution x0, the value f(x0))
    """
    if f(a)*f(b) > 0:
        print("Note: f(a)*f(b) > 0")
        return
    for k in range(max1):
        c = b - (f(b) * (b - a)) / (f(b) - f(a))
        if f(c) == 0:
            a = c
            b = c
        elif f(a)*f(c) > 0:
            a = c
        else:
            b = c
        if abs(f(c)) < epsilon:
            break
    c = (a + b) / 2
    err = abs(b - a)
    return k + 1, c, f(c)


def approot(f, X, epsilon):
    """
    To get the approximate roots of an equation
    :param f: the function
    :param X: the vector of abscissas
    :param epsilon: the tolerance
    :return: the vector of approximate root(s)
    """
    Y = list(map(f, X))
    yrange = max(Y) - min(Y)
    epsilon2 = yrange * epsilon
    n = len(X)
    R = []
    X.append(X[-1]); Y.append(Y[-1])
    for i in range(1, len(X)-1):
        if Y[i]*Y[i-1] <= 0:
            R.append((X[i] + X[i-1]) / 2)
        s = (Y[i] - Y[i-1]) * (Y[i+1] - Y[i])
        if abs(Y[i]) < epsilon2 and s <= 0:
            R.append(X[i])
    return R


def newton(f, df, p0, delta, epsilon, max1):
    """
    To find one root of an equation
    :param f: the function
    :param df: the derivative of f
    :param p0: the initial approximate root
    :param delta: the tolerance for p0
    :param epsilon: the tolerance for f(p0)
    :param max1: the maximum number of the iterations
    :return: (the number of iterations, the approximation to the solution x0, the absolute error, the value f(x0))
    """
    for k in range(max1):
        p1 = p0 - f(p0) / df(p0)
        err = abs(p0 - p1)
        try:
            relerr = err / abs(p1)
        except ZeroDivisionError:
            relerr = delta + 1
        p0 = p1
        if err < delta or relerr < delta or abs(f(p0)) < epsilon:
            break
    return k + 1, p0, err, f(p0)


def secant(f, p0, p1, delta, epsilon, max1):
    """
    To find one root of an equation
    :param f: the function
    :param p0: the initial approximate root
    :param p1: the initial approximate root
    :param delta: the tolerance for p0
    :param epsilon: the tolerance for f(p0)
    :param max1: the maximum number of the iterations
    :return: (the number of iterations, the approximation to the solution x0, the absolute error, the value f(x0))
    """
    for k in range(max1):
        p2 = p1 - (f(p1) * (p0 - p1)) / (f(p0) - f(p1))
        err = abs(p2 - p1)
        try:
            relerr = err / abs(p2)
        except ZeroDivisionError:
            relerr = delta + 1
        p0 = p1
        p1 = p2
        if err < delta or relerr < delta or abs(f(p1)) < epsilon:
            break
    return k + 1, p1, err, f(p1)
