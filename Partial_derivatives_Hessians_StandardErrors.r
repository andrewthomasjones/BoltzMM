### Computes the partial derivatives of an fvbm's model parameters
# Takes input data (a data matrix) and model (an object generated from fitfvbm)
#'@export
fvbmpartiald_R <- function(data, model)
{
  bvec <- model[[2]]
  Mmat <- model[[3]]

  N <- dim(data)[1]
  D <- length(bvec)


  partiald <- list()
  partiald[[1]] <- c()

  for (jj in 1:D) {
    partiald[[jj]] <- rep(0,1+D)
    for (ii in 1:N) {
      Holder <- data[ii,]
      inter <- Holder[jj]-tanh(sum(Holder*Mmat[,jj])+bvec[jj])
      partiald[[jj]] <- partiald[[jj]] + c(inter,Holder*inter)
    }
  }
  bvecpartial <- rep(0,D)
  Mmatpartial <- matrix(0,D,D)
  for (jj in 1:D) {
    bvecpartial[jj] <- partiald[[jj]][1]
    Mmatpartial[jj,] <- partiald[[jj]][2:(D+1)]
  }
  Mmatpartial <- Mmatpartial + t(Mmatpartial)
  diag(Mmatpartial) <- rep(0,D)
  return(list(bvecpartial=bvecpartial,Mmatpartial=Mmatpartial))
}






# Covariance Matrix
# Takes input data (a data matrix) and model (an object generated from fitfvbm)
#'@export
fvbmcov_R <- function(data,model) {
  N <- dim(data)[1]
  D <- dim(data)[2]

  I_1 <- -(1/N)*fvbmHess(data,model)
  I_2 <- matrix(0,D+D*(D-1)/2,D+D*(D-1)/2)
  for (ii in 1:N) {
    Single <- matrix(data[ii,],1,D)
    Partial_res <- fvbmpartiald_R(Single,model)
    Extract <- matrix(c(Partial_res$bvecpartial,
                                     Partial_res$Mmatpartial[lower.tri(Partial_res$Mmatpartial)]),
                      D+D*(D-1)/2,1)
    I_2 <- I_2 + Extract%*%t(Extract)
  }
  I_2 <- (1/N)*I_2
  Covar <- pracma::pinv(I_1)%*%I_2%*%pracma::pinv(I_1)
  return(Covar)
}

