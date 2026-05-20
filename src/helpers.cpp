#include <Rcpp.h>
#include "helpers.h"

int sample_j(int n){
  return (int) ( Rcpp::runif(1)(0) * n);
}