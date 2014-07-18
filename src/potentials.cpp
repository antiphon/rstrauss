#include <math.h>
#include "potentials.h"

//double potential(std::vector<double> x, 
//                std::vector<double> y, std::vector<double> z, 
//                double gamma, double R, int i, int d, int toroidal) {
//  int j;
//  double S=0, R2 = R*R, dij;
//  for(j=0; j < x.size(); j++) {
//    if(i!=j) {
//      dij = pow(x.at(i) - x.at(j), 2) + pow(y.at(i) - y.at(j), 2);
//      if(d == 3) dij += pow(z.at(i) - z.at(j), 2);
//      if(dij < R2) S += 1;
//    }
//  }
//  return pow(gamma, S);
//}
//
//
//double potential(std::vector<double> x, 
//                std::vector<double> y, 
//                std::vector<double> z, 
//                std::vector<int> these,
//                double gamma, double R, int i, int d, int toroidal) {
//  int j,k;
//  double S=0, R2 = R*R, dij;
//  for(k=0; k < these.size(); k++) {
//    j = these.at(k);
//    if(i!=j) {
//      dij = pow(x.at(i) - x.at(j), 2) + pow(y.at(i) - y.at(j), 2);
//      if(d == 3) dij += pow(z.at(i) - z.at(j), 2);
//      if(dij < R2) S += 1;
//    }
//  }
//  return pow(gamma, S);
//}

double potential(Pp X, double gamma, double R, int i) {
  int j;
  double S=0, R2 = R*R, dij;
  
  return pow(gamma, X.Rneighbours(i, R));
  
  if(X.is_blocked()){
    //printf("* getting neighbours\n");
    std::vector<int> neighs = X.block_neighbours(&i);
    //printf("* computing S\n");
    for(j=0; j < neighs.size(); j++) {
      if(neighs.at(j) != i) {
        dij = X.getDist(&i, &neighs.at(j));
          if(dij < R) S += 1;
      }
    }
  }
  else{
    for(j=0; j < X.size(); j++) {
      if(j != i) {
        dij = X.getDist(&i, &j);
          if(dij < R) S += 1;
      }
    }
  }
  return pow(gamma, S);
}

double potential(Pplite X, std::vector<int> these, double gamma, double R, int i) {
  int j, k;
  double S=0, dij;
  for(k=0; k < these.size(); k++) {
    j = these.at(k);
    if(i!=j) {
      dij = X.getDist(&i, &j);
        if(dij < R) S += 1;
      }
  }
  
  return pow(gamma, S);
}


