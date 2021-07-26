# 求解$AX=B$的迭代法

雅可比迭代：
$$
x^{(k+1)}_j = \frac{b_j - \sum \limits_{i\in[1,N]-\{j\}}{a_{ji}x_i^{(k)}}}{a_{jj}}
$$
高斯-赛德尔迭代：
$$
x^{(k+1)}_j = \frac{b_j - \sum \limits_{i=1}^{j-1}{a_{ji}x_i^{(k+1)}}-\sum \limits_{i=j+1}^{N}{a_{ji}x_i^{(k)}}}{a_{jj}}
$$
Definition 1 对于$N\times N$的矩阵$A$，如果$|a_{kk} > \sum \limits_{j=1,j≠k}^{N}{|a_{kj}|}|$，则称矩阵$A$具有严格对角优势. 

Theorem 1 矩阵$A$具有严格对角优势，则$AX=B$有唯一解$X=P$，且对于任意初始向量$P_0$，雅可比迭代序列都将收敛到$P$，高斯-赛德尔迭代也会收敛. 

判断序列是否收敛需要比较向量之间的距离，计算欧几里得距离需要较大的计算量，因此引入另一种范数$||X||_1=\sum \limits_{j=1}^N{|x_j|}$.

这样的迭代法可以扩展到非线性方程组，以下将以三维非线性方程组为例：

不动点迭代：
$$
p_{k+1} = g_1(p_k, q_k, r_k)\\
q_{k+1} = g_2(p_k, q_k, r_k)\\
r_{k+1} = g_3(p_k, q_k, r_k)\\
$$
赛德尔迭代：
$$
p_{k+1} = g_1(p_k, q_k, r_k)\\
q_{k+1} = g_2(p_{k+1}, q_k, r_k)\\
r_{k+1} = g_3(p_{k+1}, q_{k+1}, r_k)\\
$$
非线性方程组的牛顿法：

设有
$$
u = f_1(x, y)\\
v = f_2(x, y)
$$
于是
$$
\mathrm{d}u = \frac{\partial f_1}{\partial x}\mathrm{d}x + \frac{\partial f_1}{\partial y}\mathrm{d}y\\
\mathrm{d}v = \frac{\partial f_2}{\partial x}\mathrm{d}x + \frac{\partial f_2}{\partial y}\mathrm{d}y
$$
故
$$
\mathrm{d}X = J^{-1}\mathrm{d}F\\
X_k = X_{k-1} + J^{-1}(X_{k-1})(O-F(X_{k-1})) = X_{k-1} - J^{-1}(X_{k-1})F(X_{k-1})
$$
这就是一元牛顿法的一般化. 

