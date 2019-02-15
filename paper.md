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

The `BoltzMM` package (Version 0.1.3; https://CRAN.R-project.org/package=BoltzMM) for the `R` statistical programming environment [@R-Core-Team:2018aa] provides a complete suite of functions for application and estimation of *fully-visible Boltzmann Machine* (FVBM) models. The FVBM was first considered in @Hyvarinen2006, and can be described as follows.

Let $\boldsymbol{X}\in\mathbb{X}$ (where $\mathbb{X}=\{-1,+1\}^d$ is the set of *spin-binary strings* of dimension $d$) be a random variable with *probability mass function* (PMF) of the form
$$f\left(\boldsymbol{x};\boldsymbol{\theta}\right)=\mathbb{P}\left(\boldsymbol{X}=\boldsymbol{x}\right)=\exp\left(\frac{1}{2}\boldsymbol{x}^{\top}\mathbf{M}\boldsymbol{x}+\mathbf{b}^{\top}\boldsymbol{x}\right)/z\left(\boldsymbol{\theta}\right),$$
where
$$z\left(\boldsymbol{\theta}\right)=\sum_{\boldsymbol{\xi}\in\mathbb{X}}\exp\left(\frac{1}{2}\boldsymbol{\xi}^{\top}\mathbf{M}\boldsymbol{\xi}+\mathbf{b}^{\top}\boldsymbol{\xi}\right),$$
is a normalizing constant, $\mathbf{b}\in\mathbb{R}^d$, and $\mathbf{M}\in\mathbb{R}^{d\times d}$ is a symmetric matrix with zero diagonal. We put the unique elements of the *bias vector* $\mathbf{b}$ and the *interaction matrix* $\mathbf{M}$ into the *parameter vector* $\boldsymbol{\theta}$. Interpretations of the parameter vector elements of an FVBM can be found in Section 2 of @Bagnall:2018aa.

The FVBM model can be viewed as a multivariate extension of the classical Bernoulli distribution, on spin-binary random variable strings (i.e, random variables $\boldsymbol{X}\in\mathbb{X}$, rather $\boldsymbol{X}\in\{0,1\}^d$), and can be demonstrated to be equivalent to the logistic multivariate binary model that was proposed by @Cox:1972aa. It can also be considered as a fully-connected binary graphical model [see @Bagnall:2018aa for more details]. Simply put, the FVBM is a model that allows for the characterization of interactions between multiple binary random variables, simultaneously.

Let $\boldsymbol{X}_{1},\dots,\boldsymbol{X}_{n}\in\mathbb{X}$ be a sample of $n$ *independent and identically distributed* (IID) observations from an FVBM with unknown parameter vector $\boldsymbol{\theta}_{0}$. In @Hyvarinen2006, it was proved that $\boldsymbol{\theta}_{0}$ could be consistently estimated from $\boldsymbol{X}_{1},\dots,\boldsymbol{X}_{n}\in\mathbb{X}$ via the so-called *maximum pseudolikelihood estimator* (MPLE) of the type described in @Lindsay1988. Furthermore, in @NguyenWood2016, a *minorization-maximization* [MM; see @Hunter2004; @Nguyen:2017aa] was proposed for the computation of the MPLE. The asymptotic normality of the MPLE was proved in @Nguyen2016, which allows for the use of the MPLE of the FVBM as both a point estimator and a method for constructing useful hypothesis tests regarding relationships between binary random variables.

The `BoltzMM` package makes the technical results above accessible. The most fundamental of the functions of `BoltzMM` are `pvfbm` and `rfvbm`, which allow for the computation of the probability of any binary string, given valid input bias vector and interaction matrix, and the random generation of binary strings, given input parameter elements, respectively. Since the support $\mathbb{X}$ is finite, it is often desirable to compute the probabilities of all possible outcomes, which can be achieved using the function `allpfvbm`. It is also often desirable to compute the marginal probability of any particular string index
$$f_{i}\left(x_{i};\boldsymbol{\theta}\right)=\mathbb{P}\left(X_{i}=x_{i}\right),$$
where $i=1,\dots,d$, $X_i\in\{-1,1\}$, and $\boldsymbol{X}^{\top}=\left(X_{1},\dots X_{d}\right)$, for any parameter vector. This can be achieved using the function `marginpfvbm`.

The function `fitfvbm` implements the MM algorithm of @NguyenWood2016 in order to compute the MPLE of an unknown FVBM parameter vector $\boldsymbol{\theta}_0$, from an $n$ observation realization of the sample $\boldsymbol{X}_{1},\dots,\boldsymbol{X}_{n}$ (i.e., $\boldsymbol{x}_{1},\dots,\boldsymbol{x}_{n}$). That is, `fitfvbm` takes in the data as a $n\times d$ dimensional matrix, where each row is a spin-binary string of length $d$. Using the `fitfvbm` function output, and the asymptotic normality results of @Nguyen2016, standard errors can be computed using the function `fvbmstderr`, with the functions `fvbmcov`, `fvbmHess`, and `fvbmpartiald` providing intermediate calculations required by `fvbmstderr`.

As described in @Nguyen2016 and demonstrated in @Bagnall:2018aa, hypothesis testing can be conducted using Wald-type hypothesis tests. These Wald-type tests can be constructed for testing of hypotheses of the form
$$\text{H}_{0}\text{: }\theta_{0k}=\theta_{k}\text{, versus }\text{H}_{1}\text{: }\theta_{0k}\ne\theta_{k}\text{,}$$
where $\theta_{0k}$ is the $k$th element of the unknown parameter vector $\boldsymbol{\theta}_0$, and $\theta_k$ is some hypothesized value that it may take. The $p$-values and $z$-scores for these hypothesis tests can be computed using the outputs from `fitfvbm` and `fvbmstderr`, via the function `fvbmtests`.

All of the functions have been documented, and the documentation and examples for each function can be accessed via the `help` function in `R`. The analysis of data arising from the voting patterns of the Senate of the 45th Australian Parliament in 2016, using the `BoltzMM` package, can be found in @Bagnall:2018aa. The report describes the manner in which the FVBM can be used to conduct inference regarding multivariate binary data. The data that were analysed in the report are included in the package and can be accessed via the command `data(senate)`.

The `BoltzMM` package is programmed in native `R`, with computationally intensive subroutines programmed in `C` and integrated via the `Rcpp` and `RcppArmadillo` packages of @Eddelbuettel2013. The package is therefore suitable for analysis of large data sets. On https://cran.r-project.org, the only other package that provides facilities for estimation of *Boltzmann machine* (BM) models, of any kind, is the `deepnet` package of @Rong:2014aa. Here, facilities for estimation of the so-called *restricted Boltzmann machine* (RBM) are included, via the functions `rbm.down`, `rbm.train`, and `rbm.up`. 

The RBM is an alternative simplification, to the FVBM, of the BM models first considered by @Ackley1985. The BM is a parametric probabilistic model that can closely approximate any PMF on $\mathbb{X}$ [see, e.g., @Le-Roux:2008aa]. Unfortunately, due to the intractability of its requirement for hidden (or unobserved) variables, it is often difficult to efficiently and meaningfully conduct parameter estimation when both dimension $d$ and sample size $n$ are large. Like the BM models, estimation of RBM models is also intractable, due to the requirement for hidden random variables. These hidden variables also prevent the establishment of asymptotics for the RBM parameter estimators, which makes the statistical inference of such models limited. Thus we view the functions from `deepnet` as being complementary but not overlapping with `BoltzMM`.

Users can obtain the latest build of `BoltzMM` on GitHub. The latest stable build can be obtained from CRAN (https://CRAN.R-project.org/package=BoltzMM), and an archival build can be obtained from Zenodo. Algorithm derivations and theoretical results regarding the implemented functions can be found in @Nguyen2016 and @NguyenWood2016. An example application of the `BoltzMM` package for the analysis of real data is described in @Bagnall:2018aa. Thorough descriptions of the package functions appear in the manual, which can be accessed at https://cran.r-project.org/web/packages/BoltzMM/BoltzMM.pdf. Bug reports and other feedback can be directed to the GitHub issues page.

# Acknowledgements
Jessica Bagnall is funded by a Research Training Program Stipend scholarship from La Trobe University. Hien Nguyen is funded under Australian Research Council (ARC) grant numbers DE170101134 and DP180101192.

# References
