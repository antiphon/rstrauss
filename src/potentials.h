#include <Rcpp.h>
#include <vector>
using namespace Rcpp;



double potential(NumericVector x, NumericVector y, NumericVector z, 
                 double gamma, double R, int i, int n, int d);
                 
double potential(std::vector<double> x, std::vector<double> y,
                 std::vector<double> z, 
                 double gamma, double R, int i, int n, int d);
                 
double potential(List x, List y,
                 List z, 
                 double gamma, double R, int i, int n, int d);                 

double potential(std::vector<double> x, std::vector<double> y,
                 std::vector<double> z, std::vector<int> these,
                 double gamma, double R, int i, int d);
