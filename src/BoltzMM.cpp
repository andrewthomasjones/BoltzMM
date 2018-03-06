// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

//' @importFrom Rcpp sourceCpp
//' @useDynLib BoltzMM
#include "RcppArmadillo.h"

Rcpp::NumericVector export_vec(arma::vec y)
{
  Rcpp::NumericVector tmp = Rcpp::wrap(y);
  tmp.attr("dim") = R_NilValue;
  return tmp;
}

// [[Rcpp::export]]
double mahalanobis_HD(arma::rowvec y, arma::rowvec mu, arma::mat sigma)
{
  arma::rowvec d = y-mu;
  arma::mat sigma1 =  sigma;
  sigma1.diag() =  1/sigma.diag();
  double delta = (as_scalar(d*sigma1*d.t()));
  return delta;
}
