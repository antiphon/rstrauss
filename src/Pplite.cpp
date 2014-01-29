#include <Math.h>
#include <vector>
#include "Pplite.h"

/********************************************************************************************/
Pplite::~Pplite()
{
}
/********************************************************************************************/
Pplite::Pplite(std::vector<double> win, int toroidal) {
  N=0;
//  printf("Initializing Pplite");
  tor=toroidal;
  dimension = (int) win.size()/2;
  window = win;
//  for(int i=0; i < win.size();i++) window[i] = win.at(i);
  if(dimension==3) {
    if(tor) dist = &Pplite::distToroidal3D;
    else dist = &Pplite::distEuclidian3D;
  }
  else{
    if(tor) dist = &Pplite::distToroidal2D;
    else dist = &Pplite::distEuclidian2D;
  }
}
/********************************************************************************************/
int Pplite::size() {
  return N;
}
double Pplite::getX(int *i) { return x.at(*i);} ;
double Pplite::getY(int *i) { return y.at(*i);} ;
double Pplite::getZ(int *i) { return z.at(*i);} ;

/********************************************************************************************/
int Pplite::push_back(double x, double y) {
  this->x.push_back(x);
  this->y.push_back(y);
  N = this->x.size();
  return (N);
}

/********************************************************************************************/
int Pplite::push_back(double x, double y, double z) {
//  printf("Adding point");  
  this->push_back(x, y);
  if(dimension == 3) {
    this->z.push_back(z);
  };
  N = this->x.size();
//  printf(" (ok, size=%i)", N);
  return (N);
}
  
/********************************************************************************************/
void Pplite::pop_back() {
  x.pop_back();
  y.pop_back();
  if(dimension==3) z.pop_back();
  N = x.size();
}
  
/********************************************************************************************/
void Pplite::move(int *i, double x, double y){
  if(*i < N){
    this->x.at(*i) = x;
    this->y.at(*i) = y;
  };
}
/********************************************************************************************/
void Pplite::move(int *i, double x, double y, double z){
  if(*i < N){
    this->x.at(*i) = x;
    this->y.at(*i) = y;
    if(dimension==3) this->z.at(*i) = z;
  };
}

/********************************************************************************************/
void Pplite::remove(int *i) {
  x.erase(x.begin()+*i);
  y.erase(y.begin()+*i);
  if(dimension==3) z.erase(z.begin()+*i);
  N = x.size();
}
/********************************************************************************************/
double Pplite::getDist(int *i, int *j) {
  return (this->*dist)(i, j);
}
/********************************************************************************************/
double Pplite::distEuclidian2D(int *i, int *j)
{
  if(*i==*j) return 0.0;
  return sqrt(
  			pow( x.at(*i)- x.at(*j), 2.0) +
				pow( y.at(*i)- y.at(*j), 2.0) );
}

double Pplite::distToroidal2D(int *i, int *j) {
  //return Pplite::distEuclidian2D(i,j);
  if(*i==*j) return 0.0;
  double dx = fabs(x.at(*i)-x.at(*j));
  double dy = fabs(y.at(*i)-y.at(*j));
  dx = fmin( window[1]-window[0]-dx, dx ); 
  dy = fmin( window[3]-window[2]-dy, dy ); 
  	return	sqrt(
					pow( dx, 2.0 ) +
  				pow( dy, 2.0 )
          //pow( fminf( window[1]-window[0]-fabs(x.at(*i)-x.at(*j)) , fabs(x.at(*i)-x.at(*j)) ) ,2.0F) +
					//pow( fminf( window[3]-window[2]-fabs(y.at(*i)-y.at(*j)) , fabs(y.at(*i)-y.at(*j)) ) ,2.0F)
          );
  
}

/********************************************************************************************/
double Pplite::distEuclidian3D(int *i, int *j)
{
  if(*i==*j) return 0.0;
	return sqrt(
      				pow( x.at(*i)- x.at(*j), 2.0) +
      				pow( y.at(*i)- y.at(*j), 2.0) +
      				pow( z.at(*i)- z.at(*j), 2.0)   
            );
}

double Pplite::distToroidal3D(int *i, int *j) {
  if(*i==*j) return 0.0;
  	return	sqrt(
    					pow( fminf( window[1]-window[0]-fabs(x.at(*i)-x.at(*j)) , fabs(x.at(*i)-x.at(*j)) ) ,2.0) +
    					pow( fminf( window[3]-window[2]-fabs(y.at(*i)-y.at(*j)) , fabs(y.at(*i)-y.at(*j)) ) ,2.0) +
    					pow( fminf( window[5]-window[4]-fabs(z.at(*i)-z.at(*j)) , fabs(z.at(*i)-z.at(*j)) ) ,2.0)   
            );
  
}
