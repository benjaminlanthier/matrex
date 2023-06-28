# MatRexAlgs
`MatRexAlgs` stands for Matrix Reordering Algorithms.

The implemented algorithm right now is the Modified Sloan Row Ordering algorithm, or MSRO. By switching the rows (or the columns) order in the initial matrix wisely, it minimizes the mean/max front size of this matrix. This algorithm was initially implemented in Fortran in the HSL library and it is called the [mc62](https://www.hsl.rl.ac.uk/catalogue/mc62.html).

## Installation

You can install `MatRexAlgs` using `pip`:

```bash
pip install matrex
```

# The MSRO algorithm
The logic and the structure of the MSRO algorithm are described in this [paper](https://onlinelibrary.wiley.com/doi/abs/10.1002/(SICI)1099-1506(199904/05)6:3%3C189::AID-NLA160%3E3.0.CO;2-C), by Jennifer A. Scott. The difference between this and the implemented algorithm is that in this implementation, the priority function is the one described in this [paper](https://www.hsl.rl.ac.uk/specs/mc62.pdf), which is found in the HSL library website for the mc62 algorithm. To make it short, the priority function in the first paper is the following :

$$
\begin{equation}
    P_i = -W_1 \text{rcgain}_i + W_2 d(i, e),
\end{equation}
$$

and the one implemented here is the following :

$$
\begin{equation}
    P_i = -W_1 \text{rcgain}_i + W_2 d(i, e) - W_3 \text{nold}_i.
\end{equation}
$$

In both those equations, we have that :
 * $\text{rcgain}_i$ is "the increases to the row and column front sizes resulting from assembling row i next" [1]
 * $d(i, e)$ is the distance between the row $i$ and the row $e$, where $e$ is the target row, found by using the pseudodiameter of the row graph of the input matrix.
 * $\text{nold}_i$ is "the number of variables in row $i$ that are candidates for elimination and have already been brought into the front" [2].

# Dependencies
Packages needed to run this algorithm :
 * `numpy`
 * `networkx`
 * _Optional_ : `matplotlib`

# References
[1] Scott, Jennifer A. ‘A New Row Ordering Strategy for Frontal Solvers’. Numerical Linear Algebra with Applications, vol. 6, no. 3, Apr. 1999, pp. 189–211. DOI.org (Crossref), https://doi.org/10.1002/(SICI)1099-1506(199904/05)6:3<189::AID-NLA160>3.0.CO;2-C.

[2] HSL, a collection of Fortran codes for large-scale scientific computation. See http://www.hsl.rl.ac.uk/
