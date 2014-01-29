#include <vector>
#include "Point.h"

#ifndef PPL_H_
#define PPL_H_

class Pplite {
  int N;
	int tor;
	std::vector<double> window;

  int dimension;
	double (Pplite::*dist)(int*, int*);
  double distEuclidian2D(int*, int*);
  double distToroidal2D (int*, int*);
	double distEuclidian3D(int*, int*);
  double distToroidal3D (int*, int*);

  std::vector<double> x, y, z;
  
public:
  Pplite(std::vector<double> win, int toroidal);
  
  virtual ~Pplite();
  
  int size();
  
	int push_back(double x, double y);
  int push_back(double x, double y, double z);
  
  double getX(int *);
  double getY(int *);
  double getZ(int *);
  
  void pop_back();
  
  void remove(int *i);
  
  double getDist(int *, int *);
  void move(int *, double, double);
  void move(int *, double, double, double);
};

#endif