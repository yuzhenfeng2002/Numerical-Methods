# 求解$AX=B$的高斯消元和$LU$分解法

Theorem 1 若$\det(U) = a_{11}a_{22}\cdots a_{NN} ≠ 0$，则$UX=Y$存在唯一解
$$
x_N = \frac{a_{NN}}{b_N}\\
x_k = \frac{b_k - \sum \limits_{j=k+1}^N{a_{kj}x_j}}{a_{kk}}, k < N
$$

## 1 高斯消元法

### 1.1 高斯消元

Theorem 2 若$A$是非奇异矩阵，则存在线性方程组$UX=Y$与线性方程组$AX=B$等价，且$u_{kk} ≠ 0$. 于是，可以采用高斯消元法得到$UX=Y$，即可用Theorem 1中的方法求解$X$. 

### 1.2 主元选择

值得注意的是，在高斯消元的时候，会遇到待选取的主元$a_{pp}^{(p)}=0$的情况，这个时候可以采取平凡选主元策略，寻找第$p$行下满足$a_{pp}^{(p)}≠0$的第一行并交换. 

为了**减少误差的传播**，可以使用**偏序选主元策略**或按**比例偏序选主元策略**（平衡策略）. 

偏序选主元策略即寻找$a_{kp}$作为主元，$k$满足：
$$
|a_{kp}| = \max_{p\leq i\leq N}{|a_{ip}|}
$$
按比例偏序选主元策略即寻找$a_{kp}$作为主元，$k$满足：
$$
\frac{|a_{kp}|}{s_k} = \max_{p\leq i\leq N}{\frac{|a_{ip}|}{s_i}}
$$
其中$s_i$为在高斯消去的过程中每一行绝对值最大的元素. 

### 1.3 病态情况

若存在矩阵$B$，当$B$或$A$中的系数的 微小变化使得$X$变化很大时，则称矩阵$A$为病态矩阵，方程组$AX=B$称为病态方程组. 

当$A$的行列式接近于零的时候，可能发生病态情况. 

## 2 $LU$分解法

$$
AX = B,A = LU\\
\Rightarrow LUX = L(UX) = B\\
\Rightarrow UX = L^{-1}B\\
\Rightarrow X = U^{-1}L^{-1}B
$$

一般情况下，通过无行交换变换的高斯消元法可以进行$LU$分解. 但有的时候，需要进行行变换：
$$
AX = B \Rightarrow PAX = PB, PA = LU\\
\Rightarrow LUX = PB
$$

## 3 复杂度分析

对于高斯消元：

```python
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
```

需要的计算为：
$$
\sum \limits _{p=0}^{N-2}{((N-p-1) + (N-p-1)^2)} = O(n^3)
$$
对于回代法：

```python
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
```

需要的计算为：
$$
\sum \limits _{k=1}^{n-1}{(k-1)} = O(n^2)
$$
因此，计算量主要是由高斯消元造成的. 因此，对于固定的$A$，在$b$变化的情况下，只需要进行一次高斯消元，接着根据不同的$b$进行回代；也就是在实际计算中，往往采用$LU$分解法求解线性方程组. MATLAB中的`inv(A)`与`det(A)`也利用$LU$分解法. 