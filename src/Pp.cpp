#include <Math.h>
#include <vector>
#include "Pp.h"
#include "Point.h"
/********************************************************************************************/
  Pp::~Pp()
{
  }
/********************************************************************************************/
Pp::Pp(std::vector<double> win, int toroidal) {
  N=0;
//  printf("Initializing Pp");
  tor=toroidal;
  dimension = (int) win.size()/2;
	window = win;
  
  dist = &Pp::distEuclidian2D;
  
  if(dimension==3) {
    if(tor) dist = &Pp::distToroidal3D;
    else dist = &Pp::distEuclidian3D;
  }
  else{
    if(tor) dist = &Pp::distToroidal2D;
    else dist = &Pp::distEuclidian2D;
  }
}
/********************************************************************************************/
int Pp::size() {
  return N;
}
double Pp::getX(int *i) { return points.at(*i).getX();} ;
double Pp::getY(int *i) { return points.at(*i).getY();} ;
double Pp::getZ(int *i) { return points.at(*i).getZ();} ;

/********************************************************************************************/
int Pp::push_back(double x, double y) {
  return Pp::push_back(x, y, 0);
  
}

/********************************************************************************************/
int Pp::push_back(double x, double y, double z) {
//  printf("Adding point");
  Point *p;
  if(dimension == 3) p = new Point(x, y, z);
  else p = new Point(x, y);
  points.push_back(*p);
  N = points.size();
//  printf(" (ok, size=%i)", N);
  return (N-1);
}
  
/********************************************************************************************/
void Pp::pop_back() {
//  printf("Removing last point");
  points.pop_back();
  N = points.size();
//  printf(" (ok, new size=%i)\n", N);
}
  
/********************************************************************************************/
void Pp::move(int *i, double x, double y){
  points.at(*i).move(&x, &y);
}
/********************************************************************************************/
void Pp::move(int *i, double x, double y, double z){
  if(dimension==3) points.at(*i).move(&x, &y, &z);
  else points.at(*i).move(&x, &y);
}

/********************************************************************************************/
void Pp::remove(int i) {
//  printf("Removing point at location %i", i);
  points.erase(points.begin()+i);
  N = points.size();
//  printf(" (ok, new size=%i)\n", N);
}
/********************************************************************************************/
double Pp::getDist(int *i, int *j) {
  return (this->*dist)(i,j);
}
/********************************************************************************************/
double Pp::distEuclidian2D(int *i, int *j)
{
  if(*i==*j) return 0.0;
  return sqrt(
				pow( points.at(*i).getX()- points.at(*j).getX(), 2.0) +
				pow( points.at(*i).getY()- points.at(*j).getY(), 2.0) );
}

double Pp::distToroidal2D(int *i, int *j) {
  if(*i==*j) return 0.0;
  	return	sqrt(
					pow( fminf( window[1]-window[0]-fabs(points.at(*i).getX()-points.at(*j).getX()) , fabs(points.at(*i).getX()-points.at(*j).getX()) ) ,2.0F) +
					pow( fminf( window[3]-window[2]-fabs(points.at(*i).getY()-points.at(*j).getY()) , fabs(points.at(*i).getY()-points.at(*j).getY()) ) ,2.0F)
          );
  
}

/********************************************************************************************/
double Pp::distEuclidian3D(int *i, int *j)
{
  if(*i==*j) return 0.0;
	return sqrt(
				pow( points.at(*i).getX()- points.at(*j).getX(), 2.0) +
				pow( points.at(*i).getY()- points.at(*j).getY(), 2.0) +
				pow( points.at(*i).getZ()- points.at(*j).getZ(), 2.0)   );
}

double Pp::distToroidal3D(int *i, int *j) {
  if(*i==*j) return 0.0;
  	return	sqrt(
					pow( fminf( window[1]-window[0]-fabs(points.at(*i).getX()-points.at(*j).getX()) , fabs(points.at(*i).getX()-points.at(*j).getX()) ) ,2.0F) +
					pow( fminf( window[3]-window[2]-fabs(points.at(*i).getY()-points.at(*j).getY()) , fabs(points.at(*i).getY()-points.at(*j).getY()) ) ,2.0F) +
					pow( fminf( window[5]-window[4]-fabs(points.at(*i).getZ()-points.at(*j).getZ()) , fabs(points.at(*i).getZ()-points.at(*j).getZ()) ) ,2.0F)   );
  
}
