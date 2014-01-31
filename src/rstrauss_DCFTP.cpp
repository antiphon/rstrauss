#include <Rcpp.h>
#include <vector>

#include "Pplite.h"
#include "potentials.h"
#include "helpers.h"

using namespace Rcpp;

// [[Rcpp::export]]
List rstrauss_DCFTP(double beta, double gamma, double R, 
                    NumericVector win, 
                    int toroidal,int T0, int dbg, int maxtry, double blocking) {
  RNGScope scope;
  int acc = 0;
  double alpha, alpha2, Delta;
  double xnew, ynew, znew;
  double pdeath = 0.5 ; // constant
  int i, j, n=2;
  
  int dim = 2;
  if(win.size()>4) dim = 3;
  std::vector<double> window;
  for(j=0; j < win.size(); j++) window.push_back(win(j));
  
  // volume of the window
  double Volume = (win[1]-win[0])*(win[3]-win[2]);
  if(dim==3) Volume *= win[5]-win[4];
  
  // storages
  std::vector<int> D, L, U;
  std::vector<bool> it_was_death;
  std::vector<double> M;
  std::vector<int> which_id;
  Pplite xyz(window, toroidal);
  std::vector<int> id;
  
  // TODO if(blocking > 0) xyz.start_blocking(blocking);
  
  //initial pattern D = D_0
  int nD = rpois(1, beta*Volume)(0);
  
  for(i=0; i < nD; i++) {
    xnew = runif(1, win[0], win[1])(0);
    ynew = runif(1, win[2], win[3])(0);
    if(dim==3) znew = runif(1, win[4], win[5])(0);
    xyz.push_back(xnew, ynew, znew);
    D.push_back(i);
  }
  
  int k = 1, l;
  
  // main loop
  bool loop = true;
  while(loop) {
    int time = T0*k;
    
    //
    // extend dominating process
    for(i=which_id.size(); i < time; i++ ) {
      // birth or death
      pdeath = 1 - beta/(beta + D.size());
      if(runif(1)(0) <  pdeath){ // Dominating process backwards death
        j = sample_j(D.size());
        which_id.push_back(D.at(j));
        D.erase(D.begin()+j);
        M.push_back(runif(1)(0));
        it_was_death.push_back(true);
      }
      else { // Dominating process backwards birth
        xnew = runif(1, win[0], win[1])(0);
        ynew = runif(1, win[2], win[3])(0);
        if(dim==3) znew=runif(1, win[4], win[5])(0);
        xyz.push_back(xnew, ynew, znew);
        j = xyz.size()-1;
        which_id.push_back(j);
        D.push_back(j);
        M.push_back(0);
        it_was_death.push_back(false);
      }
      if(dbg) printf("      \r[T=%i] [D: %i]", time, i);
    }
     
    //
    // Run the forward chains
    U = D;
    L.clear();
    
    for(i=time-1; i>=0; i--) {
      //
      // what happened
      if(it_was_death.at(i)) { // bw death = fw birth. This is the hardest part.
        l = which_id.at(i);
        // do we get a birth in upper chain:
        alpha =  potential(xyz, L, gamma, R, l);
        // do we get a birth in the lower chain:
        alpha2 = potential(xyz, U, gamma, R, l);
        if(M.at(i) < alpha) {
          U.push_back(l);
        }
        if(M.at(i) < alpha2) {
          L.push_back(l);
        }
      }
      else{ // bw birth = fw death
        l = which_id.at(i);
        for(j=0; j<U.size();j++) if(U.at(j)==l) {U.erase(U.begin()+j); break;} // remove from U
        for(j=0; j<L.size();j++) if(L.at(j)==l) {L.erase(L.begin()+j); break;} // remove from L
      }
      //
      if(dbg) printf("        \r[T=%i] [D: %i] [UL: %i]", time, time, i);
    }
//    
    
    
    
    // check coalescence
    if(L.size() == U.size()) {
      bool fail=false;
      int in_both = 0;
      for(i=0; i< L.size(); i++) {
        for(j=0; j < U.size(); j++) {
          if(L.at(i) == U.at(j)) {
            in_both ++;
            break;
          }
        }
      }
      fail = in_both == L.size();
      loop = !fail;
    }
    if(loop & dbg) 
        printf("  No coalescence. n(L)=%i,  n(U)=%i . \n", (int) L.size(), (int) U.size());
    if(!loop &dbg) printf("Coalescence! n=%i \n", (int) L.size());
    if(time > maxtry) {printf("Warning: did not converge.\n"); loop = false;}
    k *= 2;    
  }
  
  if(dbg) printf("\n");
  // and we are done. build the result:
  std::vector<double> rx(L.size()), ry(L.size()), rz;
  for(i=0; i < L.size(); i++) {
    j = L.at(i);
    rx.at(i) = xyz.getX(&j);
    ry.at(i) = xyz.getY(&j);
  }
  if(dim==3) {
    rz.resize(L.size());
    for(i=0; i < L.size(); i++) {
      j = L.at(i);
      rz.at(i)= xyz.getZ(&j);
    }
  }
//  std::vector<double> rx(2), ry(2), rz;
  List res =   List::create(wrap(rx), wrap(ry), wrap(rz));
  return(res);
}

