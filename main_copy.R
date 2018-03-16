#'Standard errors for the parameter elements of a fitted fully-visible Boltzmann machine.
#'@description Computes the normal approximation standard errors from the sandwich estimator of the covariance matrix for a maximum pseudolikelihood estimated fully-visible Boltzmann machine.
#'@param data An N by n matrix, where each of the N rows contains a length n string of spin variables  (i.e. each element is -1 or 1).
#'@param covarmat A covariance matrix generated from \code{fvbmcov}.
#'@return A list containing 2 objects: a vector containing the standard errors corresponding to the bias parameters \code{bvec_se}, and a matrix containing the standard errors corresponding to the interaction parameters \code{Mmat_se}.
#'@references H.D. Nguyen and I.A. Wood (2016), Asymptotic normality of the maximum pseudolikelihood estimator for fully-visible Boltzmann machines, IEEE Transactions on Neural Networks and Learning Systems, vol. 27, pp. 897-902.
#'@examples # Generate num=1000 random strings of n=3 binary spin variables under bvec and Mmat.
#'@author Andrew T. Jones and Hien D. Nguyen
#'num <- 1000
#'bvec <- c(0,0.5,0.25)
#'Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
#'data <- rfvbm(num,bvec,Mmat)
#'# Fit a fully visible Boltzmann machine to data, starting from parameters bvec and Mmat.
#'model <- fitfvbm(data,bvec,Mmat)
#'# Compute the sandwich covariance matrix using the data and the model.
#'covarmat <- fvbmcov(data,model,fvbmHess)
#'# Compute the standard errors of the parameter elements according to a normal approximation.
#'fvbmstderr(data,covarmat)
#'@export
fvbmstderr <- function(data,covarmat) {
  N <- dim(data)[1]
  D <- dim(data)[2]

  stderr <- sqrt(diag(covarmat))/sqrt(N)

  bvec <- stderr[c(1:D)]
  Mmat <- matrix(0,D,D)
  Mmat[lower.tri(Mmat)] <- stderr[-c(1:D)]
  Mmat <- Mmat + t(Mmat)
  return(list(bvec_se = bvec, Mmat_se = Mmat))
}


#'Hessian of the log-pseudolikelihood function for a fitted fully-visible Boltzmann machine.
#'@description Computes the Hessian with respect to all unique parameter elements of the bias vector and interaction matrix of a fully-visible Boltzmann machine, for some random length n string of spin variables (i.e. each element is -1 or 1) and some fitted parameter values.
#'@param data An N by n matrix, where each of the N rows contains a length n string of spin variables  (i.e. each element is -1 or 1).
#'@param model List generated from \code{fitfvbm}.
#'@return The n+choose(n,2) by n+choose(n,2) Hessian matrix, summed over the N rows of \code{data} and evaluated at the fitted parameter values provided in \code{model}. Each row (column) is a unique element of the bias vector and interaction matrix. The rows are arranged in lexicographical order with the bias elements first, followed by the interaction elements. For example, if n=3, the order would be bias[1], bias[2] bias[3], interaction[1,2], interaction[1,3], and interaction[2,3].
#'@references H.D. Nguyen and I.A. Wood (2016), Asymptotic normality of the maximum pseudolikelihood estimator for fully-visible Boltzmann machines, IEEE Transactions on Neural Networks and Learning Systems, vol. 27, pp. 897-902.
#'@author Andrew T. Jones and Hien D. Nguyen
#'@examples # Generate num=1000 random strings of n=3 binary spin variables under bvec and Mmat.
#'num <- 1000
#'bvec <- c(0,0.5,0.25)
#'Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
#'data <- rfvbm(num,bvec,Mmat)
#'# Fit a fully visible Boltzmann machine to data, starting from parameters bvec and Mmat.
#'model <- fitfvbm(data,bvec,Mmat)
#'# Compute the Hessian matrix summed over all num rows of data.
#'fvbmHess(data,model)
#'@export
fvbmHess <- function(data, model) {
  bvec <- model[[2]]
  Mmat <- model[[3]]

  N <- dim(data)[1]
  D <- length(bvec)

  HessComps <- list()
  HessComps[[1]] <- matrix(0,D+1,D+1)
  for (jj in 1:D) {
    HessComps[[jj]] <- matrix(0,D+1,D+1)
    for (ii in 1:N) {
      x_bar <- as.matrix(c(1,data[ii,]),D+1,1)
      HessComps[[jj]] <- HessComps[[jj]] - x_bar%*%t(x_bar)/
        cosh(sum(Mmat[jj,]*data[ii,])+bvec[jj])^2
    }
  }

  Index <- matrix(0,D,D)
  Index[lower.tri(Index)] <- 1:(D*(D-1)/2)
  Index <- Index + t(Index)

  BigHess <- matrix(0,D+D*(D-1)/2,D+D*(D-1)/2)
  for (jj in 1:D) {
    WHICH <- which(Index[lower.tri(Index)]%in%Index[jj,])
    #Index[lower.tri(Index)] is 1:(D*(D-1)/2)
    #Index[jj,] is row j of Index


    NonZero <- HessComps[[jj]][-c(jj+1),]
    NonZero <- NonZero[,-c(jj+1)]
    BigHess[c(jj,D+WHICH),c(jj,D+WHICH)] <- BigHess[c(jj,D+WHICH),c(jj,D+WHICH)] +
      NonZero
  }

  return(BigHess)
}
