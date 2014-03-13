#include <Rcpp.h>
#include <vector>

#include "potentials.h"
#include "helpers.h"
#include "Pp.h"

using namespace Rcpp;

// [[Rcpp::export]]
List rstrauss_BD(double beta, double gamma, double R, NumericVector win, 
                 int toroidal, int iter, int dbg, double blocking) {
  RNGScope scope;
  double alpha, Delta;
  double xnew, ynew, znew;
  double xold, yold, zold;
  double pdeath = 0.5 ; // constant
  int i,j;
  
  int dim = 2;
  if(win.size()>4) dim = 3;
  //initial pattern
  std::vector<double> window;
  for(j=0; j < win.size(); j++) window.push_back(win(j));
  Pp X(window, toroidal);
  
  if(blocking > 0) X.start_blocking(blocking);
  
  xnew = runif(1, win[0], win[1])(0) ;
  ynew = runif(1, win[2], win[3])(0) ;
  if(dim==3)  znew = runif(1, win[4], win[5])(0) ;
  
  int new_id = X.push_back(xnew, ynew, znew);
  
  // volume of the window
  double Volume = (win[1]-win[0])*(win[3]-win[2]);
  if(dim==3) Volume *= win[5]-win[4];
  //
  int n=X.size();
  // main loop
  for(i=0; i < iter; i++) {
     n = X.size();
    // birth or death
    if(runif(1)(0) < pdeath) { // death?
//      printf("killing\n");
      j = sample_j(n);
      //
      Delta = potential(X, gamma, R, j);
      alpha = (n/Volume)/beta * Delta;
      if(runif(1)(0) < alpha) { // death occurs
//        printf("killed\n");
        X.remove(&j);
      }
    }
    else{ // birth, oh joy
//      printf("not killing\n");
      xnew = runif(1, win[0], win[1])(0);
      ynew = runif(1, win[2], win[3])(0);
      if(dim==3) znew = runif(1, win[4], win[5])(0);
      new_id = X.push_back(xnew, ynew, znew);
      j = n;
      Delta = potential(X, gamma, R, j);
      alpha = (Volume/(n+1.0)) * beta * Delta;
      if(runif(1)(0) < alpha) { //occurs
//        printf("born\n");
      }
      else {
        X.pop_back();
      }
    }
    
    
    if(dbg) printf("\r %i/%i [n=%i]", i+1, iter, X.size());
  }
  if(dbg) printf("\n");
  // and we are done. Compile results:
  NumericVector x(X.size()), y(X.size()), z;
  
  if(dim==3) z = rep(0, X.size());
  
  for(i=0; X.size()>i; i++) {
    x(i)=X.points.at(i).getX(); 
    y(i)=X.points.at(i).getY(); 
    if(dim==3) z(i)= X.points.at(i).getZ(); 
  }
  
  List xyz =   List::create(x, y, z);
  //
  return(xyz);
}

