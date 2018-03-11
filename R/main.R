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
