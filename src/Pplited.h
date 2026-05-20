// arbitrary dimension
#include <vector>

#ifndef PPLD_H_
#define PPLD_H_

class Pplited {
  int N;
	int tor;
	std::vector<double> window;
	double *cache;
	
  int dimension;
  
	double (Pplited::*dist)(int*, int*);
  double distEuclidian(int*, int*);
  double distToroidal (int*, int*);
	
  std::vector<std::vector <double> > x;
  
public:
  Pplited(std::vector<double> win, int toroidal);
  
  virtual ~Pplited();
  
  int size();
  
	int push_back(double *x);
  
  double getCoord(int *, int *);
  
  void pop_back();
  
  void remove(int *i);
  
  double getDist(int *, int *);
  void move(int *, double *);
  void move_cache(int *, double *);
  void move_back();
  int Rneighbours(int i, double R);
    
  
};

#endif