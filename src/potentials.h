#include <vector>
#include "Pp.h"
#include "Pplite.h"

double potential(std::vector<double> x, std::vector<double> y,
                 std::vector<double> z, 
                 double gamma, double R, int i, int d, int toroidal);
                 
double potential(std::vector<double> x, std::vector<double> y,
                 std::vector<double> z, std::vector<int> these,
                 double gamma, double R, int i, int d, int toroidal);


double potential(Pp X, double gamma, double R, int i);


double potential(Pplite , std::vector<int> , double gamma, double R, int i);
