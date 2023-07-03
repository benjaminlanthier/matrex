# MatRexAlgs
`MatRexAlgs` stands for Matrix Reordering Algorithms.

The library implements the Modified Sloan Row Ordering (MSRO) algorithm. By permuting the rows of a matrix, the algorithm minimizes the mean/max front size of this matrix. The front size minimization can be viewed in the following example :

![Alt text](/front_size.png)
![image](https://github.com/benjaminlanthier/MatRexAlgs/assets/91567620/5a367873-34aa-403e-b053-3c666c78a1b5)


**(JC: It would be nice to have a picture here or some explanation of the front size so that this term makes sense.)** The original Fortran implementation is called [mc62](https://www.hsl.rl.ac.uk/catalogue/mc62.html) and is part of the HSL library.

## Installation

You can install `MatRexAlgs` using `pip`:

```bash
pip install matrex
```

# The MSRO algorithm
This [paper](https://onlinelibrary.wiley.com/doi/abs/10.1002/(SICI)1099-1506(199904/05)6:3%3C189::AID-NLA160%3E3.0.CO;2-C) by Jennifer A. Scott describes the logic and structure of the MSRO algorithm. It iteratively selects rows to be part of the new permutation according to a priority function. In the paper, the authors use the following priority function:

$$
\begin{equation}
    P_i = -W_1 \text{rcgain}_i + W_2 d(i, e),
\end{equation}
$$

We implement a slightly different priority function described in the [manual for mc62](https://www.hsl.rl.ac.uk/specs/mc62.pdf):

$$
\begin{equation}
    P_i = -W_1 \text{rcgain}_i + W_2 d(i, e) - W_3 \text{nold}_i.
\end{equation}
$$

In both equations, we have that:
 * $\text{rcgain}_i$ is "the increases to the row and column front sizes resulting from assembling row $i$ next" [1]
 * $d(i, e)$ is the distance between the row $i$ and the row $e$, where $e$ is the target row, found by using the pseudodiameter of the row graph of the input matrix.
 * $\text{nold}_i$ is "the number of variables in row $i$ that are candidates for elimination and have already been brought into the front" [2].

# Dependencies
Packages needed to run this algorithm :
 * `numpy`
 * `networkx`
 * _Optional_ : `matplotlib` (for visualizing the example)

# References
[1] Scott, Jennifer A. ‘A New Row Ordering Strategy for Frontal Solvers’. Numerical Linear Algebra with Applications, vol. 6, no. 3, Apr. 1999, pp. 189–211. [DOI](https://doi.org/10.1002/(SICI)1099-1506(199904/05)6:3<189::AID-NLA160>3.0.CO;2-C).

[2] HSL, a collection of Fortran codes for large-scale scientific computation. See their [site](https://www.hsl.rl.ac.uk/catalogue/mc62.html).
