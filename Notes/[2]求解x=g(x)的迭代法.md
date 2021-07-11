# 求解$x=g(x)$的迭代法

Definition 1 函数$g(x)$**不动点**是指一个实数$P$ s.t. $P = g(P)$.

Definition 2 迭代$p_{n+1} = g(p_{n}), n = 0,1,\cdots$称为不动点迭代。

Theorem 1 $g$是一个连续函数，$p_n$是由不动点迭代生成的数列。如果$\lim_{n\rightarrow +\infin} {p_n}= P$，则$P$是$g(x)$的不动点。

Proof:
$$
g(P) = g(\lim_{n\rightarrow +\infin} {p_n}) = \lim_{n\rightarrow +\infin}g(p_n) = \lim_{n\rightarrow +\infin}p_{n+1} = \lim_{n\rightarrow +\infin}p_{n} = P
$$
Theorem 2 函数$g\in C[a,b]$. 若$x\in [a, b]$时，$g(x)\in [a, b]$，则$g$在$[a, b]$内存在不动点；且若$g\in D[a,b]$且存在$0<K<1$，使$|g\prime(x)| \leq K < 1$，则$g$在$[a, b]$内有唯一的不动点。

Proof:

Of theorem 2.1
$$
g(a)-a \geq 0\\
g(b)-b \leq 0\\
[g(a)-a][g(b)-b] \leq 0
$$
Of theorem 2.2

Assume that $g(a)=a, g(b) = b$, there exists $x$ s.t. $g\prime(x) = \frac{g(b)-g(a)}{b-a} = 1$.

Theorem 3 

设有（i）$g,g\prime \in C[a, b]$，（ii）$K>0$，（iii）$p_0 \in (a, b)$，（iv）$g(x)\in [a, b]$ s.t. $x\in [a, b]$.

若$|g\prime(x)| \leq K < 1$，则迭代$p_{n+1} = g(p_{n})$将收敛到唯一的不动点$P$，$P$被称为吸引不动点；若$|g\prime(x)| > 1$，则迭代$p_{n+1} = g(p_{n})$将不会收敛到$P$，$P$被称为排斥不动点。

Proof:
$$
|p_{n} - P| = |g(p_{n-1})\cdot g(P)| = |g\prime (x_0)\cdot (p_{n-1} - P)| = |g\prime (x_0)|\cdot |p_{n-1}-P| \leq K|p_{n-1}-P|\\
0 \leq |p_{n} - P| \leq K^n|p_{0}-P|\\
\lim_{n \rightarrow +\infin}|p_{n} - P| = 0
$$
