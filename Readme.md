# QuasiconvexLSE
This is an R package to compute the multivariate quasiconvex (and decreasing) nonparametric LSE as described in "Least Squares Estimation of a Monotone Quasiconvex Regression Function" by Somabha Mukherjee, Rohit K. Patra, Andrew L. Johnson, and Hiroshi Morita, which can be found at the following link:

https://arxiv.org/abs/2003.04433

To download the R package use the following in R:


```
library(devtools)
install_github(repo = "somabhamukherjee/QuasiconvexLSE")
library(QuasiconvexLSE)
```

Note that the above package requires CPLEX and the R-package RCPLEX. 
Files in the folder titled "ReplicationCode" replicate Figure 5 of https://arxiv.org/abs/2003.04433. 

**References**

Somabha Mukherjee, Rohit K. Patra, Andrew L. Johnson, and Hiroshi Morita. **Least Squares Estimation of a Monotone Quasiconvex Regression Function**. 2020. arXiv:2003.04433

