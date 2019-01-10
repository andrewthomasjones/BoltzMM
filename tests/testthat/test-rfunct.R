context("Check native R Functions")
library(BoltzMM)

test_that("Check that fvbmHess calculates Hesssian Correctly.",{
  model<-list()
  model$pll<-NA
  model$bvec<-c(0,0.5,0.25)
  model$Mmat<- matrix(0.1,3,3) - diag(0.1,3,3)
  model$itt<-NA
  set.seed(1)
  data <- matrix(sample(c(-1,1),300, replace =TRUE), ncol=3)
  HessResult<-fvbmHess(data, model)

  tmp1 <- "./fvbmHess"

  # The first run always succeeds, but warns
  expect_known_output(HessResult, tmp1, print = TRUE)

})

test_that("Check that fvbmstderr calculates stderr Correctly.",{
  set.seed(1)
  data <- matrix(sample(c(-1,1),300, replace =TRUE), ncol=3)
  bvec <- c(0,0.5,0.25)
  Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
  covarmat <- diag(length(bvec)[1]*2)
  stderrResult<-fvbmstderr(data, covarmat)
  tmp2 <- "./fvbmstderr"

  # The first run always succeeds, but warns
  expect_known_output(stderrResult, tmp2, print = TRUE)

})

test_that("Check that marginpfvbm calculates marginal probailities Correctly.",{
  bvec <- c(0,0.5,0.25)
  Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
  marginResult<-marginpfvbm(bvec, Mmat)
  tmp2 <- "./marginpfvbm"

  # The first run always succeeds, but warns
  expect_known_output(marginResult, tmp2, print = TRUE)

})

test_that("Check that fvbmtests calculates scores and p-values correctly.",{
  set.seed(1)
  bvec <- c(0,0.5,0.25)
  Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
  num <- 1000
  data <- rfvbm(num,bvec,Mmat)
  model <- fitfvbm(data,bvec,Mmat)
  nullmodel <- list(bvec = c(0,0,0), Mmat = matrix(0,3,3))
  testResult<-fvbmtests(data,model,nullmodel)
  tmp2 <- "./fvbmtests"
  # The first run always succeeds, but warns
  expect_known_output(testResult, tmp2, print = TRUE)

})
