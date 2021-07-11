def fixpt(g, p0, tol, max1):
    P = [p0]
    n = 0
    for k in range(max1):
        pk = g(P[-1])
        P.append(pk)
        n += 1
        err = abs(P[-1] - P[-2])
        relerr = abs(P[-1] - P[-2]) / abs(P[-1])
        if  err < tol or relerr < tol:
            break
    if n == max1:
        print("Maximum number of iterations exceeded!")
    return P
