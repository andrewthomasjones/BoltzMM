
<!-- README.md is generated from README.Rmd. Please edit that file -->
BoltzMM
=======

The BoltzMM package allows for computation of probability mass functions of fully-visible Boltzmann machines via and . Random data can be generated using . Maximum pseudolikelihood estimation of parameters via the MM algorithm can be conducted using . Computation of partial derivatives and Hessians can be performed via and . Covariance estimation and normal standard errors can be computed using and .

Installation
------------

You can install BoltzMM from github with:

``` r
# install.packages("devtools")
devtools::install_github("andrewthomasjones/BoltzMM")
```

Examples
--------

Compute the probability of every length n=3 binary spin vector under bvec and Mmat:

``` r
bvec <- c(0,0.5,0.25)
Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
allpfvbm(bvec,Mmat)
```

Generate num=1000 random strings of n=3 binary spin variables under bvec and Mmat.

``` r
num <- 1000
bvec <- c(0,0.5,0.25)
Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
data <- rfvbm(num,bvec,Mmat)
```

Fit a fully visible Boltzmann machine to data, starting from parameters bvec and Mmat.

``` r
fitfvbm(data,bvec,Mmat)
```

For more examples see individual help files.
