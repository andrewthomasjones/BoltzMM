// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]


//' @importFrom Rcpp sourceCpp
//' @useDynLib BoltzMM
//'
#include "RcppArmadillo.h"
#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include <cmath>

//for returning a vector within a list
Rcpp::NumericVector export_vec(arma::vec y)
{
  Rcpp::NumericVector tmp = Rcpp::wrap(y);
  tmp.attr("dim") = R_NilValue;
  return tmp;
}


//creating binary vectors
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
    double prob = 0.0;
    double norm = 0.0;
    int count = std::pow(2,n);

    if(xval.n_elem!= n | Mmat.n_rows!=n | Mmat.n_cols!=n | Mmat.n_rows!=Mmat.n_cols ){
      Rcpp::Rcerr << "Input variable dimensions do not match";
    }else{

        for(int i=0; i<count; i++){
            arma::vec zeta_i = bin_vec(i,n);
            norm += as_scalar(arma::exp(0.5*zeta_i.t()*Mmat*zeta_i+arma::dot(bvec,zeta_i)));
        }

        prob = as_scalar(arma::exp(0.5*xval.t()*Mmat*xval+ arma::dot(bvec,xval)));
        prob /= norm;
    }
    return prob;
}

// # allpfvbm -- Return the probabilities for every string in {-1,1}^n from a fvbm model
// # with bias bvec and relationship matrix Mmat.
// # Order is as given by expand.grid(c(-1,1),c(-1,1),...).
// # function with inputs: bvec, Mmat
//'@export
// [[Rcpp::export]]
arma::rowvec allpfvbm(arma::vec bvec, arma::mat Mmat) {
    int n = bvec.n_elem;
    double norm = 0.0;
    int count = std::pow(2,n);
    arma::rowvec probvec = arma::zeros(count).t();

    if(Mmat.n_rows!=n | Mmat.n_cols!=n | Mmat.n_rows!=Mmat.n_cols ){
      Rcpp::Rcerr << "Input variable dimensions do not match";
    }else{
        for(int i=0; i<count; i++){
          arma::vec zeta_i = bin_vec(i,n);
          double prob = as_scalar(arma::exp(0.5*zeta_i.t()*Mmat*zeta_i+arma::dot(bvec, zeta_i)));
          probvec(i) = prob;
          norm +=  prob;
        }

        probvec /= norm;
    }
    return(probvec);
}

// #rfvbm -- Generate random data from fvbm model with bias bvec and relationship matrix Mmat.
// # num is the number of observations to be drawn.
// # function with inputs bvec, Mmat, num
//'@export
// [[Rcpp::export]]
arma::mat rfvbm(int num, arma::vec bvec, arma::mat Mmat) {
  int n = bvec.n_elem;
  arma::mat returnmat   = arma::zeros(num,n);

  if(Mmat.n_rows!=n | Mmat.n_cols!=n | Mmat.n_rows!=Mmat.n_cols ){
    Rcpp::Rcerr << "Input variable dimensions do not match";
  }else{

      arma::rowvec cumprob = arma::cumsum(allpfvbm(bvec,Mmat));
      //need to fix this, cant call system RNG if CRAN
      arma::vec random_nums = arma::randu(num);
      for(int i=0; i<num; i++){
        int j = arma::as_scalar(find(cumprob > random_nums(i), 1, "first"));
        arma::vec zeta_j = bin_vec(j,n);
        returnmat.row(i) = zeta_j.t();
      }
  }
  return(returnmat);
}


// #fitfvbm -- Fit fvbm with starting parameters bvec and Mmat, using the MM algorithm.
// # Function returns estimated parameters and final pseudo-log-likelihood.
// # Takes input bvec, Mmat, data, delta_crit.
// # delta_crit a termination criterion based on the relative error of the distance between parameter iterates.
//'@export
// [[Rcpp::export]]
Rcpp::List fitfvbm(arma::mat data, arma::vec bvec, arma::mat Mmat, double delta_crit = 0.001, int max_it = 1000){
    //New parameters transfer into old parameters
    int N = data.n_rows;
    int D = bvec.n_elem;
    int itt = 0;

    if(Mmat.n_rows!=D | Mmat.n_cols!=D | data.n_cols !=D | data.n_cols!=Mmat.n_cols | data.n_cols!=Mmat.n_rows |Mmat.n_rows!=Mmat.n_cols){
        Rcpp::Rcerr << "Input variable dimensions do not match";
        return(Rcpp::List::create());
    }

    arma::mat MM = Mmat;
    arma::mat dataj = data;
    arma::mat datak = data;
    arma::vec temp = arma::zeros(N);
    double delta = delta_crit+10.0;
    arma::mat par  = arma::join_rows(bvec,MM);
    arma::mat old_par = par;

    double DERIV = 0.0;
    double LIKE = 0.0;
    double sumj = 0.0;
    double sumk = 0.0;

    while (delta > delta_crit & itt<max_it)
    {
        itt++;
        old_par = par;

        for(int j=0; j<D; j++)
        {
            DERIV = arma::sum(data.col(j));
            dataj = data.each_row()%MM.col(j).t();
            DERIV -= arma::as_scalar(arma::sum(arma::tanh(sum(dataj,1) + bvec(j))));
            bvec(j) +=  DERIV/N;
        }



        for(int j=0; j<D; j++)
        {
            for(int k=(j+1); k<D; k++)
            {
                DERIV = 2.0*arma::dot(data.col(j),data.col(k));

                dataj = data.each_row()%MM.col(j).t();
                datak = data.each_row()%MM.col(k).t();

                sumk = arma::as_scalar(arma::dot(data.col(k),(arma::tanh(sum(dataj,1) + bvec(j)))));
                sumj = arma::as_scalar(arma::dot(data.col(j),(arma::tanh(sum(datak,1) + bvec(k)))));

                DERIV-=(sumk+sumj);

                MM(j,k) += DERIV/(2.0*N);
                MM(k,j) = MM(j,k);
            }
        }

        par = arma::join_rows(bvec,MM);
        delta = std::sqrt(arma::accu(arma::pow((par-old_par),2.0)))/std::max(std::sqrt(std::pow(arma::accu(old_par),2.0)),1.0);
    }

    LIKE = 0.0;
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<D; j++)
        {
            LIKE += data(i,j)*arma::dot(MM.col(j),data.row(i)) + bvec(j)*data(i,j) - std::log(std::cosh(arma::dot(MM.col(j),data.row(i))+bvec(j))) - std::log(2.0);
        }
    }

    Rcpp::List retList = Rcpp::List::create(
        Rcpp::Named("pll")= LIKE,
        Rcpp::Named("bvec")= export_vec(bvec),
        Rcpp::Named("Mmat")= MM,
        Rcpp::Named("itt")= itt
    );

    return(retList);
}

