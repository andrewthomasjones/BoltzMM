// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]


//' @importFrom Rcpp sourceCpp
//' @useDynLib BoltzMM
//'
#include "RcppArmadillo.h"
#include <bitset>
#include <boost/dynamic_bitset.hpp>


Rcpp::NumericVector export_vec(arma::vec y)
{
  Rcpp::NumericVector tmp = Rcpp::wrap(y);
  tmp.attr("dim") = R_NilValue;
  return tmp;
}


// [[Rcpp::export]]
arma::vec bin_vec(int y,  int n)
{
  const boost::dynamic_bitset<> b(n, y);
  arma::vec x = arma::zeros(n);
  for(int i=0; i<n; i++){
    x(i) = 2*(b[i]-0.5) ;
  }
  return(x);
}

// # pfvbm -- Generate probability of string xval occuring with a fvbm model
// # with bias bvec and relationship matrix Mmat.
// # Note that xval is a string of {-1,1}^n and Mmat is a symmetric matrix with zeros on the diag.
// # function with inputs: bvec, Mmat, xval
//'@export
// [[Rcpp::export]]
double pfvbm(arma::vec xval, arma::vec bvec, arma::mat Mmat) {
    int n = bvec.n_elem;
    // add checks that dimensions on xval, bvec and Mmat all match up
    // and are in range etc

    double prob = 0.0;
    double norm = 0.0;

    int count = std::pow(2,n);

    for(int i=0; i<count; i++){
        arma::vec zeta_i = bin_vec(i,n);
        norm += as_scalar(arma::exp(0.5*zeta_i.t()*Mmat*zeta_i+arma::dot(bvec,zeta_i)));
    }

    prob = as_scalar(arma::exp(0.5*xval.t()*Mmat*xval+ arma::dot(bvec,xval)));
    prob /= norm;
    return prob;
}

// # allpfvbm -- Return the probabilities for every string in {-1,1}^n from a fvbm model
// # with bias bvec and relationship matrix Mmat.
// # Order is as given by expand.grid(c(-1,1),c(-1,1),...).
// # function with inputs: bvec, Mmat
//'@export
// [[Rcpp::export]]
arma::rowvec allpfvbm( arma::vec bvec, arma::mat Mmat) {
  int n = bvec.n_elem;
  // add checks that dimensions on xval, bvec and Mmat all match up
  // and are in range etc

  double norm = 0.0;

  int count = std::pow(2,n);
  arma::rowvec probvec = arma::zeros(count).t();

  for(int i=0; i<count; i++){
    arma::vec zeta_i = bin_vec(i,n);
    double prob = as_scalar(arma::exp(0.5*zeta_i.t()*Mmat*zeta_i+arma::dot(bvec, zeta_i)));
    probvec(i) = prob;
    norm +=  prob;
  }

  probvec /= norm;
  return(probvec);
}

// #rfvbm -- Generate random data from fvbm model with bias bvec and relationship matrix Mmat.
// # num is the number of observations to be drawn.
// # function with inputs bvec, Mmat, num
//'@export
// [[Rcpp::export]]
arma::mat rfvbm(int num,arma::vec bvec, arma::mat Mmat) {
  // add checks that dimensions on xval, bvec and Mmat all match up
  // and are in range etc

  //symetric matrix speed ups??

  int n = bvec.n_elem;
  arma::mat returnmat   = arma::zeros(num,n);
  arma::rowvec cumprob = arma::cumsum(allpfvbm(bvec,Mmat));

  arma::vec random_nums = arma::randu(num);
  for(int i=0; i<num; i++){
    int j = as_scalar(find(cumprob > random_nums(i), 1, "first"));
    arma::vec zeta_j = bin_vec(j,n);
    returnmat.row(i) = zeta_j.t();
  }

  return(returnmat);
}


