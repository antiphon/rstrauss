#include <vector>

#ifndef POINT_H_
#define POINT_H_

class Point
{
  double x;
	double y;
	double z;
	double mark;
	int id;
  std::vector<int> neighbours;

public:
	Point();
	Point(double, double);
	Point(double, double, double);
	virtual ~Point();

	double  getX();
	double  getY();
	double  getZ();
	double  getMark();
	void    setMark(double *);

	int 	getId();
	void 	setId(int *);

  int	 	nsize();
	int   getNeighbour(int *);
	void  addNeighbour(int *i);
	void  removeNeighbour(int *i);
	void  clearNeighbourhood();

	void  move(double*, double*);
	void	move(double*, double*, double*);
};

#endif /*POINT_H_*/
