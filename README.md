# MatRexAlgs
`MatRexAlgs` stands for Matrix Reordering Algorithms.

The implemented algorithm right now is called [MSRO](https://onlinelibrary.wiley.com/doi/abs/10.1002/(SICI)1099-1506(199904/05)6:3%3C189::AID-NLA160%3E3.0.CO;2-C). By switching the columns (or the rows) order in the initial matrix, it minimizes the mean/max front size of this matrix. This algorithm was initially implemented in Fortran in the [HSL](https://www.hsl.rl.ac.uk/catalogue/mc62.html) library.

Packages needed to run this algorithm :
 * `numpy`
 * `networkx`
 * `matplotlib`