# Covariance Matrix
# Takes input data (a data matrix) and model (an object generated from fitfvbm)
#'@export
fvbmcov_p <- function(data,model) {
  N <- dim(data)[1]
  D <- dim(data)[2]

  I_1 <- -(1/N)*fvbmHess(data,model)
  I_2 <- matrix(0,D+D*(D-1)/2,D+D*(D-1)/2)
  for (ii in 1:N) {
    Single <- matrix(data[ii,],1,D)
    Partial_res <- fvbmpartiald(Single,model)
    Extract <- matrix(c(Partial_res$bvec,
                        Partial_res$Mmat[lower.tri(Partial_res$Mmat)]),
                      D+D*(D-1)/2,1)
    I_2 <- I_2 + Extract%*%t(Extract)
  }
  I_2 <- (1/N)*I_2
  Covar <- pracma::pinv(I_1)%*%I_2%*%pracma::pinv(I_1)
  return(Covar)
}

# Covariance Matrix
# Takes input data (a data matrix) and model (an object generated from fitfvbm)
#'@export
fvbmcov_i <- function(data,model) {
  N <- dim(data)[1]
  D <- dim(data)[2]

  I_1 <- -(1/N)*fvbmHess(data,model)
  I_2 <- matrix(0,D+D*(D-1)/2,D+D*(D-1)/2)
  for (ii in 1:N) {
    Single <- matrix(data[ii,],1,D)
    Partial_res <- fvbmpartiald(Single,model)
    Extract <- matrix(c(Partial_res$bvec,
                        Partial_res$Mmat[lower.tri(Partial_res$Mmat)]),
                      D+D*(D-1)/2,1)
    I_2 <- I_2 + Extract%*%t(Extract)
  }
  I_2 <- (1/N)*I_2
  Covar <- solve(I_1)%*%I_2%*%solve(I_1)
  return(Covar)
}


#bad
DATA <- rfvbm(1000,c(2,4,5),matrix(0.5,3,3)+diag(-0.5,3))
FIT <- fitfvbm(DATA,c(2,4,5),matrix(0.5,3,3)+diag(-0.5,3))
COV1<-fvbmcov_p(DATA,FIT)
COV2<-fvbmcov_i(DATA,FIT)
COV3<-fvbmcov(DATA,FIT,fvbmHess)

#good
DATA <- rfvbm(1000,c(0,1,0),matrix(c(0,1,1,1,0,1,1,1,0),3,3))
FIT <- fitfvbm(DATA,c(0,1,0),matrix(c(0,1,1,1,0,1,1,1,0),3,3))
COV1<-fvbmcov_p(DATA,FIT)
COV2<-fvbmcov_i(DATA,FIT)
COV3<-fvbmcov(DATA,FIT,fvbmHess)


