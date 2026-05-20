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
	
  bool use_blocking;
  std::vector<int> blocks;  
  std::vector<double> blocksides;
  
  std::vector<std::vector<int> > blockmembers;
  
  std::vector<double> cache;
  
  std::vector<double> window;
  int dimension;
	double (Pp::*dist)(int*, int*);
  bool (Pp::*distLT)(int *, int*, double*);
  double distEuclidian2D(int*, int*);
  double distToroidal2D (int*, int*);
	double distEuclidian3D(int*, int*);
  double distToroidal3D (int*, int*);
  bool distLessThan2D(int *, int *, double *);
  bool distLessThanToroidal2D(int *, int *, double *);
  bool distLessThan3D(int *, int *, double *);
  bool distLessThanToroidal3D(int *, int *, double *);

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
  
  void start_blocking(double );
  bool block_neighbours(int *, int *);
  void set_block(int *);
  std::vector<int> get_block_vec(int *);
  int get_block_index(int *);
  
  int block_ijk_to_index(int *i, int *j, int *k);
  int block_ijk_to_index(int i, int j, int k);
  
  std::vector<int> block_index_to_ijk(int *index);
  
  std::vector<int> block_neighbours(int *);
  
  bool is_blocked();
  
  void pop_back();
  
  void remove(int *i);
  void remove_blockmembership(int *i);
  void shift_blocking_indices_by_one(int *i);
  
  
  double getDist(int *, int *);
  void move(int *, double, double);
  void move(int *, double, double, double);
  void move_cache(int *, double, double, double);
  void move_back();
  
  int Rneighbours(int i, double);
};

int min(int, int);
int max(int, int);

#endif  