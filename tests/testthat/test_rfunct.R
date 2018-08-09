context("Check native R Functions")

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
