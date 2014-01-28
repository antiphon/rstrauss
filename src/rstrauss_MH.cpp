#include <Rcpp.h>
#include "potentials.h"
using namespace Rcpp;

// [[Rcpp::export]]
List rstrauss_MH(int dim, int n, double gamma, double R, NumericVector win, int iter, int dbg) {
  RNGScope scope;
  //initial pattern
  NumericVector x, y, z;
  x  = runif(n, win[0], win[1]);
  y  = runif(n, win[2], win[3]);
  if(dim==3) z = runif(n, win[2], win[3]);
  // then we loop
  int acc = 0;
  double E_old, E_new;
  double xnew, ynew, znew;
  double xold, yold, zold;
  double alpha;
  int j;
  for(int i=0; i < iter; i++) {
    xnew = runif(1, win[0], win[1])(0);
    ynew = runif(1, win[2], win[3])(0);
    j = i%n;
    E_old = potential(x, y, z, gamma, R, j, n, dim);
    xold = x(j);
    x(j) = xnew;
    yold = y(j);
    y(j) = ynew;
    if(dim==3) {
      zold = z(j);
      znew = runif(1, win[4], win[5])(0);
      z(j) = znew;
    }
    E_new = potential(x, y, z, gamma, R, j, n, dim);
    if(E_old == 0 & E_new > 0) {alpha = 1;}
    else if(E_new == 0) { alpha = 0;}
    else {alpha = E_new/E_old; }
    if(runif(1)(0) < alpha) {
      acc += 1;
    }
    else {
      x(j) = xold;
      y(j) = yold;
      if(dim==3) z(j) = zold;
    }
    if(dbg) printf("\r %i/%i", i+1, iter);
  }
  if(dbg) printf("\n");
  // and we are done
  List  xyz =   List::create(x, y, z);
  return(xyz);
}

