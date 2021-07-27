#include <Rcpp.h>
#include <vector>

#include "potentials.h"
#include "helpers.h"
#include "Pp.h"

using namespace Rcpp;

// [[Rcpp::export]]
List rstrauss_BD_delta(double beta, double gamma, double delta, 
                       double R, NumericVector win, 
                       int toroidal, int iter, int dbg, double blocking) {
  
  RNGScope scope;
  double alpha, Delta, pen;
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
  
  xnew = R::runif(win[0], win[1]) ;
  ynew = R::runif(win[2], win[3]) ;
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
    pen = n > 1 ? n*(n-2) : 0;
    pen = exp(delta*pen);
    // birth or death
    if(R::runif(0,1) < pdeath) { // death?
//      printf("killing\n");
      j = sample_j(n);
      //
      Delta =  potential(X, gamma, R, j)/pen;
      alpha = (n/Volume)/beta / Delta;
      if(R::runif(0,1) < alpha) { // death occurs
//        printf("killed\n");
        X.remove(&j);
      }
    }
    else{ // birth, oh joy
//      printf("not killing\n");
      xnew = R::runif(win[0], win[1]);
      ynew = R::runif(win[2], win[3]);
      if(dim==3) {
        znew = R::runif(win[4], win[5]);
        new_id = X.push_back(xnew, ynew, znew);
      }
      else new_id = X.push_back(xnew, ynew);
      j = n;
      Delta = pen * potential(X, gamma, R, j);
      alpha = (Volume/(n+1.0)) * beta * Delta;
      if(R::runif(0,1) < alpha) { //occurs
//        printf("born\n");
      }
      else {
        X.pop_back();
      }
    }
    
    
    if(dbg) Rprintf("\r %i/%i [n=%i]", i+1, iter, X.size());
  }
  if(dbg) Rprintf("\n");
  // and we are done. Compile results:
  NumericVector x(X.size()), y(X.size()), z(0);
  
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

