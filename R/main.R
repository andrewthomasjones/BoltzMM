# Get Standard errors
# takes an input covariance matrix from the fvbmcov function and sample data "data"
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


### Computes the Hessian of an fvbm's model parameters
# Takes input data (a data matrix) and model (an object generated from fitfvbm)
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
