#include <math.h>
#include <vector>
#include "Pplited.h"

/********************************************************************************************/
Pplited::~Pplited() 
{
}
/********************************************************************************************/
Pplited::Pplited(std::vector<double> win, int toroidal) {
  N=0;
  tor=toroidal;
  dimension = (int) win.size()/2;
  cache = new double[dimension + 1];
  window = win;
  if(tor) dist = &Pplited::distToroidal;
    else dist = &Pplited::distEuclidian;
  x.resize(dimension);
}
/********************************************************************************************/
int Pplited::size() {
  return N;
}
double Pplited::getCoord(int *i, int *l) { return x.at(*l).at(*i);} ;

/********************************************************************************************/
int Pplited::push_back(double *x) {
  for(int l = 0; l < dimension; l++) this->x.at(l).push_back(x[l]);
  N = this->x.at(0).size();
  return (N);
}
/********************************************************************************************/
void Pplited::pop_back() {
  for(int l=0; l < dimension; l++) x.at(l).pop_back();
  N = x.at(0).size();
}
  
/********************************************************************************************/
void Pplited::move(int *i, double *x){
  if(*i < N){
    for(int l=0; l < dimension; l++)this->x.at(l).at(*i) = x[l];
  };
}

void Pplited::move_cache(int *i, double *x){
  
  cache[dimension] = (double)*i;
  for(int l=0; l< dimension; l++) cache[l] = getCoord(i, &l);
  move(i, x);
}

void Pplited::move_back(){
  int i = (int) cache[dimension];
  move(&i, cache);
}
/********************************************************************************************/
void Pplited::remove(int *i) {
  for(int l=0; l < dimension; l++) x.at(l).erase(x.at(l).begin()+*i);
  N = x.at(0).size();
}
/********************************************************************************************/
double Pplited::getDist(int *i, int *j) {
  return (this->*dist)(i, j);
}
/********************************************************************************************/
int Pplited::Rneighbours(int i, double R) {
  int S=0;
  for(int j=0; j < N; j++) {
    if(j != i) {
      if( (this->*dist)(&i, &j) < R ) S += 1;
    }
  }
  return S;
}
/********************************************************************************************/
double Pplited::distEuclidian(int *i, int *j)
{
  if(*i==*j) return 0.0;
  double d=0;
  for(int l=0; l < dimension; l++) d += pow(x.at(l).at(*i)-x.at(l).at(*j), 2);
  return sqrt(d);
}

double Pplited::distToroidal(int *i, int *j) {
  if(*i==*j) return 0.0;
  
  double dx;
  double d=0;
  for(int l=0; l < dimension; l++) {
    dx = fabs(x.at(l).at(*i)-x.at(l).at(*j));
    dx = fmin( window[2*l+1]-window[2*l]-dx, dx ); 
    d += dx * dx;
  }
  return sqrt(d);
}
/********************************************************************************************/
