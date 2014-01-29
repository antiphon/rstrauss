#include <Math.h>
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
  for(j=0; j < X.size(); j++) {
  if(i!=j) {
    dij = X.getDist(&i, &j);
      if(dij < R) S += 1;
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


