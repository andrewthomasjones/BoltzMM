
<img src="http://www.r-pkg.org/badges/version-last-release/BoltzMM"></img></a> [![Downloads from the RStudio CRAN mirror](http://cranlogs.r-pkg.org/badges/BoltzMM)](https://CRAN.R-project.org/package=BoltzMM) [![Build Status](https://travis-ci.org/andrewthomasjones/BoltzMM.svg?branch=master)](https://travis-ci.org/andrewthomasjones/BoltzMM)

<!-- README.md is generated from README.Rmd. Please edit that file -->
BoltzMM
=======

The BoltzMM package allows for computation of probability mass functions of fully-visible Boltzmann machines via `pfvbm` and `allpfvbm`. Random data can be generated using `rfvbm`. Maximum pseudolikelihood estimation of parameters via the MM algorithm can be conducted using `fitfvbm`. Computation of partial derivatives and Hessians can be performed via `fvbmpartiald` and `fvbmHessian`. Covariance estimation and normal standard errors can be computed using `fvbmcov` and `fvbmstderr`.

Installation
------------

Installation
------------

If `devtools` has already been installed, then the most current build of `BoltzMM` can be obtained via the command:

``` r
devtools::install_github('andrewthomasjones/BoltzMM',build_vignettes = T)
```

The latest stable build of `BoltzMM` can be obtain from CRAN via the command:

``` r
install.packages("BoltzMM", repos='http://cran.us.r-project.org')
```

An archival build of `BoltzMM` is available at <https://zenodo.org/record/1317784>. Manual installation instructions can be found within the *R* installation and administration manual <https://cran.r-project.org/doc/manuals/r-release/R-admin.html>.

Examples
--------

Compute the probability of every length n=3 binary spin vector under bvec and Mmat:

``` r
library(BoltzMM)
set.seed(1)

bvec <- c(0,0.5,0.25)
Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
allpfvbm(bvec,Mmat)
#>           [,1]       [,2]      [,3]      [,4]       [,5]       [,6]
#> [1,] 0.0666189 0.04465599 0.1213876 0.1213876 0.07362527 0.07362527
#>           [,7]      [,8]
#> [1,] 0.2001342 0.2985652
```

Generate num=1000 random strings of n=3 binary spin variables under bvec and Mmat.

``` r
library(BoltzMM)
set.seed(1)

num <- 1000
bvec <- c(0,0.5,0.25)
Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
data <- rfvbm(num,bvec,Mmat)

head(data)
#>      [,1] [,2] [,3]
#> [1,]    1    1   -1
#> [2,]   -1   -1    1
#> [3,]   -1    1    1
#> [4,]    1    1    1
#> [5,]   -1    1   -1
#> [6,]    1    1    1
```

Fit a fully visible Boltzmann machine to data, starting from parameters bvec and Mmat.

``` r
library(BoltzMM)
set.seed(1)

bvec <- c(0,0.5,0.25)
Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
data <- rfvbm(num,bvec,Mmat)

fitfvbm(data,bvec,Mmat)
#> $pll
#> [1] -1892.661
#> 
#> $bvec
#> [1] 0.02607382 0.46484595 0.27640931
#> 
#> $Mmat
#>           [,1]      [,2]      [,3]
#> [1,] 0.0000000 0.1179001 0.1444486
#> [2,] 0.1179001 0.0000000 0.0351134
#> [3,] 0.1444486 0.0351134 0.0000000
#> 
#> $itt
#> [1] 5
```

For more examples see individual help files.

Unit testing
------------

Using the package `testthat`, we have conducted the following unit test for the GitHub build, on the date: 10 November, 2018. The testing files are contained in the [tests](https://github.com/andrewthomasjones/BoltzMM/tree/master/tests) folder of the respository.

``` r

## Load 'BoltzMM' library.
library(BoltzMM)

## Load 'testthat' library.
library(testthat)

## Test 'BoltzMM'.
#test_package("BoltzMM")
```

Bug reporting and contributions
-------------------------------

Thank you for your interest in `BoltzMM`. If you happen to find any bugs in the program, then please report them on the Issues page (<https://github.com/andrewthomasjones/BoltzMM/issues>). Support can also be sought on this page. Furthermore, if you would like to make a contribution to the software, then please forward a pull request to the owner of the repository.
