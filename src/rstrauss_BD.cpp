#include <Rcpp.h>
#include "potentials.h"
#include "helpers.h"
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
List rstrauss_BD(int dim, double beta, double gamma, double R, NumericVector win, 
                 int iter, int dbg) {
  RNGScope scope;
  int acc = 0;
  double alpha, Delta;
  double xnew, ynew, znew;
  double xold, yold, zold;
  double pdeath = 0.5 ; // constant
  int j;
  //initial pattern
  std::vector<double> x, y, z;
  x.push_back( runif(1, win[0], win[1])(0) );
  y.push_back( runif(1, win[2], win[3])(0) );
  if(dim==3) z.push_back( runif(1, win[4], win[5])(0) );
  int n = 1;
  
  // volume of the window
  double Volume = (win[1]-win[0])*(win[3]-win[2]);
  if(dim==3) Volume *= win[5]-win[4];
  //
  // main loop
  for(int i=0; i < iter; i++) {
    // birth or death
    if(runif(1)(0) < pdeath) { // death
      //printf("killing\n");
      j = sample_j(n);
      //
      Delta = potential(x, y, z, gamma, R, j, n, dim);
      alpha = (n/Volume)/beta * Delta;
      if(runif(1)(0) < alpha) { // occurs
        //printf("killed\n");
        x.erase(x.begin()+j);
        y.erase(y.begin()+j);
        if(dim==3) z.erase(z.begin()+j);
        n -= 1;
      }
    }
    else{ // birth, oh joy
      //printf("living\n");
      xnew = runif(1, win[0], win[1])(0);
      ynew = runif(1, win[2], win[3])(0);
      if(dim==3) znew = runif(1, win[4], win[5])(0);
      x.push_back(xnew);
      y.push_back(ynew);
      if(dim==3) z.push_back(znew);
      j = x.size()-1;
      Delta = potential(x, y, z, gamma, R, j, n+1, dim);
      alpha = (Volume/(n+1.0))*beta * Delta;
      if(runif(1)(0) < alpha) { //occurs
        //printf("birth\n");
        n += 1;
      }
      else { 
        x.pop_back();
        y.pop_back();
        if(dim==3) z.pop_back();
      }
    }
    
    
    if(dbg) printf("\r %i/%i [n=%i]", i+1, iter, n);
  }
  if(dbg) printf("\n");
  // and we are done
  List xyz =   List::create(wrap(x), wrap(y), wrap(z));
  return(xyz);
}

