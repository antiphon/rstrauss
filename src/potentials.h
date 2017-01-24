#include <vector>
#include "Pp.h"
#include "Pplite.h"
#include "Pplited.h"

double potential(std::vector<double> x, std::vector<double> y,
                 std::vector<double> z, 
                 double gamma, double R, int i, int d, int toroidal);
                 
double potential(std::vector<double> x, std::vector<double> y,
                 std::vector<double> z, std::vector<int> these,
                 double gamma, double R, int i, int d, int toroidal);


double potential(Pp X, double gamma, double R, int i);

double potential_lj(Pp X, double sigma, double epsilon, int i);

double potential(Pplite , std::vector<int> , double gamma, double R, int i);

double potential(Pplited , double gamma, double R, int i);
 