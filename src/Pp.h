/*
difference to Pp.h: use vectors of x,y,z instead of Point-vector

*/

#include <vector>
#include "Point.h"

#ifndef PP_H_
#define PP_H_

class Pp {
  int N;
	int tor;
	std::vector<double> window;
  int dimension;
	double (Pp::*dist)(int*, int*);
  double distEuclidian2D(int*, int*);
  double distToroidal2D (int*, int*);
	double distEuclidian3D(int*, int*);
  double distToroidal3D (int*, int*);
public:
  std::vector<Point> points;
  Pp(std::vector<double> win, int toroidal);
  virtual ~Pp();
  
  int size();
  
	int push_back(double x, double y);
  int push_back(double x, double y, double z);
  
  double getX(int *);
  double getY(int *);
  double getZ(int *);
  
  void pop_back();
  
  void remove(int i);
  
  double getDist(int *, int *);
  void move(int *, double, double);
  void move(int *, double, double, double);
};

#endif