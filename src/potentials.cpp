#include "potentials.h"

double potential(NumericVector x, NumericVector y, NumericVector z, double gamma, double R, int i, int n, int d) {
  int j;
  double S=0, R2 = R*R, dij;
  for(j=0; j < n; j++) {
    if(i!=j) {
      dij = pow(x(i) - x(j), 2) + pow(y(i) - y(j), 2);
      if(d == 3) dij += pow(z(i) - z(j), 2);
      if(dij < R2) S += 1;
    }
  }
  return pow(gamma, S);
}

double potential(std::vector<double> x, 
                std::vector<double> y, std::vector<double> z, 
                double gamma, double R, int i, int n, int d) {
  int j;
  double S=0, R2 = R*R, dij;
  for(j=0; j < n; j++) {
    if(i!=j) {
      dij = pow(x.at(i) - x.at(j), 2) + pow(y.at(i) - y.at(j), 2);
      if(d == 3) dij += pow(z.at(i) - z.at(j), 2);
      if(dij < R2) S += 1;
    }
  }
  return pow(gamma, S);
}


double potential(std::vector<double> x, 
                std::vector<double> y, 
                std::vector<double> z, 
                std::vector<int> these,
                double gamma, double R, int i, int d) {
  int j,k;
  double S=0, R2 = R*R, dij;
  for(k=0; k < these.size(); k++) {
    j = these.at(k);
    if(i!=j) {
      dij = pow(x.at(i) - x.at(j), 2) + pow(y.at(i) - y.at(j), 2);
      if(d == 3) dij += pow(z.at(i) - z.at(j), 2);
      if(dij < R2) S += 1;
    }
  }
  return pow(gamma, S);
}


