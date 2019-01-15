---
title: 'BoltzMM: an R package for maximum pseudolikelihood estimation of fully-visible
  Boltzmann machines'
authors:
- affiliation: 1
  name: Andrew T. Jones
- affiliation: 2
  name: Jessica J. Bagnall
- affiliation: 2
  name: Hien D. Nguyen
  orcid: 0000-0002-9958-432X
date: "11 January, 2019"
bibliography: paper.bib
tags:
- artificial neural network
- graphical model
- maximum pseudolikelihood estimation
- multivariate binary data
- probability mass function
- R
affiliations:
- index: 1
  name: School of Mathematics and Physics, University of Queensland, St. Lucia 4072,
    Queensland Australia
- index: 2
  name: Department of Mathematics and Statistics, La Trobe University, Bundoora 3086,
    Victoria Australia
---

# Summary

The *Boltzmann machine* (BM), first considered by @Ackley1985, is a rich parametric probabilistic artificial neural network that is able densely represent any *probability mass function* (PMF) over the support of *spin-binary strings* $\mathbb{X}=\{-1,+1\}^d$ [see, e.g., @Le-Roux:2008aa]. Unfortunately, due to its complex form, and latent variable construction, it is often difficult to efficiently and meaningfully conduct parameter estimation when $d$ is large and sample size $n$. Furthermore, it is often not necessary to consider such elaborate forms for practical modeling.

In @Hyvarinen2006, the *fully-visible BM* (FVBM) was first considered as a tractable simplification of the BM. Let $\boldsymbol{X}\in\mathbb{X}$ be such that the PMF can be written as
$$f\left(\boldsymbol{x};\boldsymbol{\theta}\right)=\mathbb{P}\left(\boldsymbol{X}=\boldsymbol{x}\right)=\exp\left(\frac{1}{2}\boldsymbol{x}^{\top}\mathbf{M}\boldsymbol{x}+\mathbf{b}^{\top}\boldsymbol{x}\right)/z\left(\boldsymbol{\theta}\right),$$
where
$$z\left(\boldsymbol{\theta}\right)=\sum_{\boldsymbol{\xi}\in\mathbb{X}}\exp\left(\frac{1}{2}\boldsymbol{\xi}^{\top}\mathbf{M}\boldsymbol{\xi}+\mathbf{b}^{\top}\boldsymbol{\xi}\right),$$
is a normalizing constant, $\mathbf{b}\in\mathbb{R}^d$, and $\mathbf{M}\in\mathbb{R}^{d\times d}$ is a symmetric matrix with zero diagonal. We put the unique elements of the *bias vector* $\mathbf{b}$ and the *interaction matrix* $\mathbf{M}$ into the *parameter vector* $\boldsymbol{\theta}$. Interpretations of the parameter vector elements can be found in Section 2 of @Bagnall:2018aa.

The FVBM model can be viewed as a multivariate extension of the classical Bernoulli distribution, on spin-binary random variable strings (i.e, random variables $\boldsymbol{X}\in\mathbb{X}$, rather $\boldsymbol{X}\in\{0,1\}^d$), and can be demonstrated to be equivalent to the logistic multivariate binary model that was proposed by @Cox:1972aa, and can be considered as a fully-connected binary graphical model [see @Bagnall:2018aa for more details].

Let $\boldsymbol{X}_{1},\dots,\boldsymbol{X}_{n}\in\mathbb{X}$ be a sample of $n$ *independent and identically distributed* (IID) observations from an FVBM with unknown parameter vector $\boldsymbol{\theta}_{0}$. In @Hyvarinen2006, it was proved that $\boldsymbol{\theta}_{0}$ could be consistently estimated from $\boldsymbol{X}_{1},\dots,\boldsymbol{X}_{n}\in\mathbb{X}$ via the so-called *maximum pseudolikelihood estimator* (MPLE) of the type described in @Lindsay1988 and @Arnold1991. Furthermore, in @NguyenWood2016, an efficient and globally convergent *minorization-maximization* [MM; see @Hunter2004; @Nguyen:2017aa] was proposed for the computation of the MPLE. The asymptotic normality of the MPLE was proved in @Nguyen2016, which allows for the use of the MPLE of the FVBM as both a point estimator and a method for constructing useful hypothesis tests regarding relationships between binary random variables.

The `BoltzMM` package (Version 0.1.3; https://CRAN.R-project.org/package=BoltzMM) for the `R` statistical programming environment [@R-Core-Team:2018aa] provides a complete suite of functions for application and estimation of FVBM models. The most fundamental of the functions are `pvfbm` and `rfvbm`, which allow for the computation of the probability of any binary string, given valid input bias vector and interaction matrix, and the random generation of binary strings, given input parameter elements, respectively. Since the support $\mathbb{X}$ is finite, it is often desirable to compute the probabilities of all possible outcomes, which can be achieved using the function `allpfvbm`. It is also often desirable to compute the marginal probability of any particular string index
$$f_{i}\left(x_{i};\boldsymbol{\theta}\right)=\mathbb{P}\left(X_{i}=x_{i}\right),$$
where $i=1,\dots,d$, $X_i\in\{-1,1\}$, and $\boldsymbol{X}^{\top}=\left(X_{1},\dots X_{d}\right)$, for any parameter vector. This can be achieved using the function `marginpfvbm`.

The function `fitfvbm` implements the MM algorithm of @NguyenWood2016 in order to compute the MPLE of an unknown FVBM parameter vector $\boldsymbol{\theta}_0$, from an $n$ observation realization of the sample $\boldsymbol{X}_{1},\dots,\boldsymbol{X}_{n}$ (i.e., $\boldsymbol{x}_{1},\dots,\boldsymbol{x}_{n}$). That is, `fitfvbm` takes in the data as a $n\times d$ dimensional matrix, where each row is a spin-binary string of length $d$. Using the `fitfvbm` function output, and the asymptotic normality results of @Nguyen2016, standard errors can be computed using the function `fvbmstderr`, with the functions `fvbmcov`, `fvbmHess`, and `fvbmpartiald` providing intermediate calculations required by `fvbmstderr`.

As described in @Nguyen2016 and demonstrated in @Bagnall:2018aa, hypothesis testing can be conducted using Wald-type hypothesis tests. These Wald-type tests can be constructed for testing of hypotheses of the form
$$\text{H}_{0}\text{: }\theta_{0k}=\theta_{k}\text{, versus }\text{H}_{1}\text{: }\theta_{0k}\ne\theta_{k}\text{,}$$
where $\theta_{0k}$ is the $k$th element of the unknown parameter vector $\boldsymbol{\theta}_0$, and $\theta_k$ is some hypothesised value that it may take. The $p$-values and $z$-scores for these hypothesis tests can be computed using the outputs from `fitfvbm` and `fvbmstderr`, via the function `fvbmtests`.

All of the functions that have been comprehensively documented, and the documentation and example applications of each function can be accessed via the `help` function in `R`. Furthermore, the analysis of data arising from the voting patterns of the Senate of the 45th Australian Parliament in 2016, using the `BoltzMM` package, can be found in @Bagnall:2018aa. The report clearly describes the manner in which the FVBM can be used to conduct meaningful inference regarding multivariate binary data. The data that were analysed in the report are included in the package and can be accessed via the command `data(senate)`.

The `BoltzMM` package is programmed in native `R`, with particularly computationally intensive subroutines programmed in `C` and integrated via the `Rcpp` and `RcppArmadillo` packages of @Eddelbuettel2013. The package is therefore sufficiently suitable for analysis of moderately large dimensional data sets. On https://cran.r-project.org, the only other package that provides facilities for estimation of BMs, of any kind, is the `deepnet` package of @Rong:2014aa. Here, facilities for estimation of the so-called *restricted BM* (RBM) are included, via the functions `rbm.down`, `rbm.train`, and `rbm.up`. The RBM is an alternative simplification to the BM, which preserves the denseness properties BM but can be parameterized by fewer parameter elements. However, the estimation of RBM models is as intractable as for BM models, and due to the latent variable construction preventing the establishment of distributional results for the estimators of the RBM parameter elements, its use for statistical inference remains limited. Thus we view the functions from `deepnet` as being complimentary but not overlapping with `BoltzMM`.

Users can obtain the latest build of `BoltzMM` on GitHub (https://github.com/andrewthomasjones/BoltzMM). The latest stable build can be obtained from CRAN (https://CRAN.R-project.org/package=BoltzMM), and an archival build can be obtained from Zenodo (http://doi.org/10.5281/zenodo.2538256). Algorithm derivations and theoretical results regarding the implented functions can be found in @Nguyen2016 and @NguyenWood2016. An example application of the `BoltzMM` package for the analysis of real data is comprehensively described in @Bagnall:2018aa.   Thorough descriptions of the package functions appear in the manual, which can be accessed at https://cran.r-project.org/web/packages/BoltzMM/BoltzMM.pdf. Bug reports and other feedback can be directed to the GitHub issues page (https://github.com/andrewthomasjones/BoltzMM/issues).

# Acknowledgements
Jessica Bagnall is funded by a Research Training Program Stipend scholarship from La Trobe University. Hien Nguyen is funded under Australian Research Council (ARC) grant numbers DE170101134 and DP180101192.

# References
