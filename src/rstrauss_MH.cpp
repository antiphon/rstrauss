#include <Rcpp.h>
#include <vector>
#include "potentials.h"
#include "helpers.h"
#include "Pp.h"

using namespace Rcpp;

// [[Rcpp::export]]
List rstrauss_MH(int n, double gamma, double R, NumericVector win, 
                 int toroidal, int iter, int dbg, double blocking, NumericMatrix start) {
  RNGScope scope;
  int i, j;
  
  int dim = 2;
  if(win.size()>4) dim = 3;
  
  
  
  //initial pattern
  std::vector<double> window;
  for(j=0; j < win.size(); j++) window.push_back(win(j));
  Pp X(window, toroidal);
  
  if(blocking > 0) X.start_blocking(blocking);
  
  double xnew, ynew, znew;
  
//  if(start_is_empty){ 
//    printf("Got null start.");
//    for(i=0; i < n; i++){
//      xnew = runif(1, win[0], win[1])(0) ;
//      ynew = runif(1, win[2], win[3])(0) ;
//      if(dim==3)  znew = runif(1, win[4], win[5])(0) ;
//      int new_id = X.push_back(xnew, ynew, znew);
//    }
//  }
//  else{
  for(i=0; i < n; i++){
    xnew = start(i,0);
    ynew = start(i,1);
    if(dim==3)  znew = start(i,2);
    int new_id = X.push_back(xnew, ynew, znew);
  }
//  }
  
  // then we loop
  int acc = 0;
  double E_old, E_new;
  double xold, yold, zold;
  double alpha;
  
  for(i=0; i < iter; i++) {
    j = sample_j(n); //i%n;
//    printf("potential\n");
    E_old = potential(X, gamma, R, j);
    
    xnew = runif(1, win[0], win[1])(0);
    ynew = runif(1, win[2], win[3])(0);
    xold = X.getX(&j);
    yold = X.getY(&j);
    if(dim==3) {
      znew = runif(1, win[4], win[5])(0);
      zold = X.getZ(&j);
    }
//    printf("moving\n");
    X.move_cache(&j, xnew, ynew, znew);
//    printf("potential\n");
    E_new = potential(X, gamma, R, j);
    if(E_old == 0 & E_new > 0) {alpha = 1;}
    else if(E_new == 0) { alpha = 0;}
    else {alpha = E_new/E_old; }
    if(runif(1)(0) < alpha) {
      acc += 1;
    }
    else {
//      printf("moving back\n");
//      X.move(&j, xold, yold, zold);
      X.move_back();
    }
    if(dbg) printf("\r %i/%i", i+1, iter);
  }
  if(dbg) printf(" MH done.\n");
  
  // and we are done. Compile results:
  NumericVector x(X.size()), y(X.size()), z;
  if(dim==3) z = rep(0, X.size());
  for(i=0; i < X.size(); i++) {
    x(i)=X.points.at(i).getX(); 
    y(i)=X.points.at(i).getY(); 
    if(dim==3) z(i)= X.points.at(i).getZ(); 
  }
  List xyz =   List::create(x, y, z);
  return(xyz);
}

