
# To do: 1. pfvbm, rfvbm
# 2. Fit fvbm
# 3. prbm, rrbm
# 4. Fit rbm
# 5. pbm, rbm
# 6. Fit bm

# pfvbm -- Generate probability of string xval occuring with a fvbm model
# with bias bvec and relationship matrix Mmat.
# Note that xval is a string of {-1,1}^n and Mmat is a symmetric matrix with zeros on the diag.
# function with inputs: bvec, Mmat, xval
pfvbm <- function(xval,bvec,Mmat) {
  nn <- length(bvec)
  zeta <- eval(parse(text=paste('expand.grid(',paste(rep('c(-1,1)',nn),collapse = ','),
                                ',stringsAsFactors = F)',collapse ='')))
  zeta <- as.matrix(zeta)
  norm <- 0
  for (ii in 1:(2^nn)) {
    norm <- norm + exp(0.5*zeta[ii,]%*%Mmat%*%zeta[ii,]+sum(bvec*zeta[ii,]))
  }
  xval <- matrix(xval,1,nn)
  prob <- exp(0.5*xval%*%Mmat%*%t(xval)+sum(bvec*xval))
  prob <- prob/norm
  return(prob)
}

# allpfvbm -- Return the probabilities for every string in {-1,1}^n from a fvbm model
# with bias bvec and relationship matrix Mmat.
# Order is as given by expand.grid(c(-1,1),c(-1,1),...).
# function with inputs: bvec, Mmat
allpfvbm <- function(bvec,Mmat) {
  nn <- length(bvec)
  zeta <- eval(parse(text=paste('expand.grid(',paste(rep('c(-1,1)',nn),collapse = ','),
                                ',stringsAsFactors = F)',collapse ='')))
  zeta <- as.matrix(zeta)
  probvec <- c()
  norm <- 0
  for (ii in 1:(2^nn)) {
    prob <- exp(0.5*zeta[ii,]%*%Mmat%*%zeta[ii,]+sum(bvec*zeta[ii,]))
    probvec[ii] <- prob
    norm <- norm + prob
  }
  return(probvec/norm)
}

#rfvbm -- Generate random data from fvbm model with bias bvec and relationship matrix Mmat.
# num is the number of observations to be drawn.
# function with inputs bvec, Mmat, num
rfvbm <- function(num,bvec,Mmat) {
  nn <- length(bvec)
  zeta <- eval(parse(text=paste('expand.grid(',paste(rep('c(-1,1)',nn),collapse = ','),
                                ',stringsAsFactors = F)',collapse ='')))
  zeta <- as.matrix(zeta)
  cumprob <- cumsum(allpfvbm(bvec,Mmat))
  returnmat <- matrix(NA,num,nn)
  for (ii in 1:num) {
    returnmat[ii,] <- zeta[which(cumprob>=stats::runif(1))[1],]
  }
  return(returnmat)
}

#fitfvbm -- Fit fvbm with starting parameters bvec and Mmat, using the MM algorithm.
# Function returns estimated parameters and final pseudo-log-likelihood.
# Takes input bvec, Mmat, data, delta_crit.
# delta_crit a termination criterion based on the relative error of the distance between parameter iterates.
fitfvbm <- function(data,bvec,Mmat,delta_crit=0.001) {

  # New parameters transfer into old parameters
  N <- dim(data)[1]
  D <- length(bvec)
  MM <- Mmat
  BB <- matrix(bvec,D,1)
  delta <- Inf
  old_par <- par <- c(BB,MM)

  while (delta > delta_crit)
  {
    old_par <- par

    for (jj in 1:D)
    {
      DERIV <- 0
      for (ii in 1:N)
      {
        DERIV <- DERIV + DATA[ii,jj] - tanh(MM[,jj]%*%DATA[ii,]+BB[jj,1])
      }
      BB[jj,1] <- BB[jj,1] + DERIV/N
    }


    for (jj in 1:D)
    {
      for (kk in jj:D)
      {
        if (jj != kk)
        {
          DERIV <- 0
          for (ii in 1:N)
          {
            DERIV <- DERIV + 2*DATA[ii,jj]*DATA[ii,kk] - DATA[ii,kk]*tanh(MM[,jj]%*%DATA[ii,]+BB[jj,1]) - DATA[ii,jj]*tanh(MM[,kk]%*%DATA[ii,]+BB[kk,1])
          }
          MM[jj,kk] <- MM[kk,jj] <- DERIV/(2*N) + MM[jj,kk]
        }
      }
    }
    par <- c(BB,MM)

    delta <- sqrt(sum((par-old_par)^2))/max(sqrt(sum(old_par)^2),1)
  }

  LIKE <- 0
  for (ii in 1:N)
  {
    for (jj in 1:D)
    {
      LIKE <- LIKE + DATA[ii,jj]*MM[,jj]%*%DATA[ii,] + BB[jj,1]*DATA[ii,jj] - log(cosh(MM[,jj]%*%DATA[ii,]+BB[jj,1])) - log(2)
    }
  }

  return(list(pll=LIKE,bvec=BB,Mmat=MM))
}
