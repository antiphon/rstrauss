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
  double xnew, ynew, znew=0;
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
  
  Rprintf("%f %f, %f %f\n", win(0), win(1), win(2), win(3));
  
  xnew = R::runif(win(0), win(1));
  ynew = R::runif(win(2), win(3));
  if(dim==3)  znew = R::runif(win[4], win[5]) ;
  
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
    if(R::runif(0,1) < pdeath) { // death?
//      printf("killing\n");
      j = sample_j(n);
      //
      Delta = potential(X, gamma, R, j);
      alpha = (n/Volume)/beta * Delta;
      if(R::runif(0,1) < alpha) { // death occurs
//        printf("killed\n");
        X.remove(&j);
      }
    }
    else{ // birth, oh joy
//      printf("not killing\n");
      xnew = R::runif(win(0), win(1));
      ynew = R::runif(win(2), win(3));
      
      //Rprintf("%f %f\n", xnew, ynew);
      if(dim==3) {
        znew = R::runif(win(4), win(5));
        new_id = X.push_back(xnew, ynew, znew);
      }
      else new_id = X.push_back(xnew, ynew);
      j = n;
      Delta = potential(X, gamma, R, j);
      alpha = (Volume/(n+1.0)) * beta * Delta;
      if(R::runif(0,1) < alpha) { //occurs
        //Rprintf("b:(%f, %f)\n", xnew, ynew);
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
  NumericVector x(X.size()), y(X.size()), z(0);
  
  if(dim==3) z = rep(0, X.size());
  
  for(i=0; i < X.size(); i++) {
    x(i) = X.points.at(i).getX(); 
    y(i) = X.points.at(i).getY(); 
    if(dim==3) z(i)= X.points.at(i).getZ(); 
  }
  List xyz =   List::create(x, y, z);
  //
  return(xyz);
}

