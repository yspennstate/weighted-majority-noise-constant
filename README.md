# The sharp small-noise constant for weighted majority

Yitzchak Shmalo

Flip each coordinate of a uniform point of {-1,1}^n independently with probability eps. Peres proved that a linear threshold function then changes its value with probability at most 2 sqrt(eps), with asymptotic constant sqrt(2/pi), and asked whether the constant can be lowered to 2/pi, the value for simple majority. We prove that every linear threshold function changes with probability at most (2/pi) sqrt(eps) + C eps^(3/5) for an absolute constant C. This proves the weak form of the Majority Is Least Stable conjecture recorded in O'Donnell's list of open problems. The sequences of threshold functions that attain the constant are exactly those in which the threshold b tends to zero and, for every fixed c > 0, the weights of size at least c sqrt(eps) carry a vanishing share of the squared norm. The bound of Peres does not extend from symmetric to merely centered inputs.

- `main.tex`, `main.pdf`: the paper.
- `check.py`: computes both sides of the identity of Lemma 2.1 in exact arithmetic for six linear threshold functions, the inequality of Lemma 3.3 for 60 random weight vectors, the noise sensitivity of simple majority against the main term of Proposition 3.2, and the bound of Proposition 5.1 at n = 256 (Python 3, standard library only; under a minute). Its output is `check_output.txt`.
