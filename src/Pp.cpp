#include <Math.h>
#include <stdlib.h>
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
  use_blocking = false;
  
  dist = &Pp::distEuclidian2D;
  cache.resize(5);
  
  // set distance metrics. Avoid if's as these are called alot.
  if(dimension==3) {
    if(tor) {
      dist = &Pp::distToroidal3D;
      distLT = &Pp::distLessThanToroidal3D;
    }
    else {
      dist = &Pp::distEuclidian3D;
      distLT = &Pp::distLessThan3D;
    }
  }
  else{
    if(tor) {
      dist = &Pp::distToroidal2D;
      distLT = &Pp::distLessThanToroidal2D;
    }
    else {
      dist = &Pp::distEuclidian2D;
      distLT = &Pp::distLessThan2D;
    }
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
  int idx=N-1;
  set_block(&idx);
  return idx;
}
  
/********************************************************************************************/
void Pp::pop_back() {
  int m = points.size()-1;
  remove(&m);
}
  
/********************************************************************************************/
void Pp::move(int *i, double x, double y){
  points.at(*i).move(&x, &y);
  remove_blockmembership(i);
  set_block(i);
}
/********************************************************************************************/
void Pp::move(int *i, double x, double y, double z){
  if(dimension==3) points.at(*i).move(&x, &y, &z);
  else points.at(*i).move(&x, &y);
  remove_blockmembership(i);
  set_block(i);
}

void Pp::move_cache(int *i, double x, double y, double z){
  cache.at(0) = (double)*i;
  cache.at(1) = points.at(*i).getX();
  cache.at(2) = points.at(*i).getY();
  cache.at(3) = points.at(*i).getZ();
  if(use_blocking) cache.at(4) = points.at(*i).getId();
  move(i, x, y, z);
}

/********************************************************************************************/
void Pp::move_back(){
  double x=cache.at(1), y=cache.at(2), z=cache.at(3);
  int i = cache.at(0);
  if(dimension==3) points.at(i).move(&x, &y, &z);
  else points.at(i).move(&x, &y);
  remove_blockmembership(&i);
  if(use_blocking){
    int k = cache.at(4);
    points.at(i).setId(&k);
    // book keeping
    blockmembers.at(k).push_back(i);
  }
}
/********************************************************************************************/
void Pp::remove_blockmembership(int *i){
  if(use_blocking){
    int j, found=0, k = points.at(*i).getId();
//    printf("* removing bm: %i , %i ", *i, k);
    for(j=0; j < blockmembers.at(k).size(); j++) if(blockmembers.at(k).at(j)==*i){found=1; break;}
    if(found) {
      blockmembers.at(k).erase( blockmembers.at(k).begin()+j );
//      printf(" (ok)\n");
    }
    else{ printf("** %i not in right block ?\n", *i);}
  }
}
/********************************************************************************************/
void Pp::shift_blocking_indices_by_one(int *i) {
      if(*i < N-1){
        // must update all numbers from i->n! @TODO find a better way
        int j,k;        
        for(j=0; j < blockmembers.size(); j++)
        for(k=0; k < blockmembers.at(j).size(); k++) 
          if(blockmembers.at(j).at(k) > *i) blockmembers.at(j).at(k)--;
      }
}

/********************************************************************************************/
void Pp::remove(int *i) {
//  printf("Removing point at location %i", i);
  remove_blockmembership(i);
  shift_blocking_indices_by_one(i);
  points.erase(points.begin()+*i);
  N = points.size();
//  printf(" (ok, new size=%i)\n", N);
}


/********************************************************************************************/
void Pp::start_blocking(double max_range) {
  
  // compute grid: in each dimension, what is the min grid size > max_range
  int nx, ny, nz, nblocks=0;
  if(max_range > window.at(1) - window.at(0)) nx = 1;
  else  nx = (int)  (window.at(1)-window.at(0))/max_range;
  double ax = (window.at(1)-window.at(0))/nx;
  if(max_range > window.at(3) - window.at(2)) ny = 1;
  else ny = (int)  (window.at(3)-window.at(2))/max_range;
  double ay = (window.at(3)-window.at(2))/ny;
  blocksides.push_back(ax);
  blocksides.push_back(ay);
  blocks.push_back(nx);
  blocks.push_back(ny);
  nblocks = nx*ny;
  if(dimension==3) {
    if(max_range > window.at(5) - window.at(4)) nz = 1;
    nz = (int)  (window.at(5)-window.at(4))/max_range ;
    double az = (window.at(5)-window.at(4))/nz;
    blocksides.push_back(az);
    blocks.push_back(nz);
    nblocks *= nz;
  }
  blockmembers.resize(nblocks);
  //printf("Blocking ready(%i, %i)\n", blockmembers.size(), blockmembers.at(0).size());
  use_blocking = true;
};
/********************************************************************************************/
std::vector<int> Pp::get_block_vec(int *i) {
  std::vector<int> ijk(3);
  if(!use_blocking) return ijk;// if not blocking.
  // set the block id of point i.
  // figure out which block is the i in:
  // in x-dimension:
  ijk.at(0) = (int) ceil( blocks.at(0) * ( points.at(*i).getX() - window.at(0) )/(window.at(1)-window.at(0)) ) - 1;
  ijk.at(1) = (int) ceil( blocks.at(1) * ( points.at(*i).getY() - window.at(2) ) / (window.at(3)-window.at(2)) ) - 1;
  if(dimension==3) ijk.at(2) = (int) ceil( blocks.at(2) * ( points.at(*i).getZ() - window.at(4) )/ ( window.at(5)-window.at(4) ) ) - 1;
  return ijk;
}
/********************************************************************************************/
int Pp::get_block_index(int *i) {
  if(!use_blocking) return 0;// if not blocking.
  std::vector<int> ijk = get_block_vec(i);
  int k = Pp::block_ijk_to_index(&ijk[0], &ijk[1], &ijk[2]);
  return k;
}
/********************************************************************************************/
int Pp::block_ijk_to_index(int *i, int *j, int *k){
  // convert (i,j,k) to block index in grid blocking
  // (0,0,0), (nx-1, 0, 0), (0, ny-1, 0) ...  are the corners
  return (*i) + blocks.at(0)*(*j) + blocks.at(0)*blocks.at(1)*(*k);
}
int Pp::block_ijk_to_index(int i, int j, int k){
  return block_ijk_to_index(&i, &j, &k);
}
/********************************************************************************************/
std::vector<int> Pp::block_index_to_ijk(int *index){
  std::vector<int> ijk(3);
  int a, ab, idx= *index;
  a = blocks.at(0);
  ab = blocks.at(0)*blocks.at(1);
  int k=0;
  if(dimension==3) k = idx / ab;
  int j = (idx-k*ab) / a;
  int i = idx - k*ab - j*a;
  ijk.at(0)=i;
  ijk.at(1)=j;
  ijk.at(2)=k;
  return ijk;
}
/********************************************************************************************/
void Pp::set_block(int *i) {
  if(!use_blocking) return;// if not blocking.
  // set the block id of point i.
  // figure out which block is the i in:
  int k = Pp::get_block_index(i);
  // use the id variables of Point.
  points.at(*i).setId(&k);
  // book keeping
  blockmembers.at(k).push_back(*i);
//  printf("p: %i -> b: %i/%i\n ", *i, k, (int) blockmembers.at(k).size());
}

/********************************************************************************************/
// return the point indices that are in the neighbouring blocks
std::vector<int> Pp::block_neighbours(int *ii) {
  std::vector<int> np;
  int i,j,k;
  // just in case, should be taken care elsewhere as this is slow
  if(!use_blocking) {
    np.resize(points.size());
    for(i=0; i < points.size(); i++) np.at(i) = i;
    return np;
  }
  int idx = points.at(*ii).getId();
  // the neighbouring blocks:
  std::vector<int> nb, ijk;
  
  //printf("** idx: %i\n ", idx);
  //printf("** pushed to np: %i\n ", (int)nb.size());
  ijk = block_index_to_ijk(&idx);
  //printf("** ijk: %i-%i\n", ijk.at(0), ijk.at(1));
  int ni, nj, nk=0;
  for(i=-1; i < 2; i++){
    ni = ijk.at(0)+i;
    if(tor) {
      if(ni<0) ni = blocks.at(0)+ni;
      else if(ni >= blocks.at(0)) ni = ni-blocks.at(0);
    }
    if( (ni >-1) & (ni < blocks.at(0))){
      for(j=-1; j < 2; j++){
        nj = ijk.at(1)+j;
        if(tor) {
          if(nj<0) nj = blocks.at(1)+nj;
          else if(nj >= blocks.at(1)) nj = nj-blocks.at(1);
        }
        if( (nj >-1) & (nj < blocks.at(1))){
          if(dimension==3){
            for(k=-1; k < 2; k++){  
              nk = ijk.at(2)+k;
              if(tor) {
                if(nk<0) nk = blocks.at(2)+nk;
                else if(nk >= blocks.at(2))nk = nj-blocks.at(2);
              }
              if( (nk >-1) & (nk < blocks.at(2)) ) {
                nb.push_back( block_ijk_to_index(ni, nj, nk) );
              }
            }
          }
          else{
            nb.push_back( block_ijk_to_index(ni, nj, nk) );            
          }
        }
      }
    }
  }
  //printf("** collection, nb: %i\n", (int) nb.size());
  int b;
  for(i=0; i < nb.size(); i++){
    b = nb.at(i);
    //printf("** b: %i , n=%i\n",  b, (int) blockmembers.at(b).size());
    for(j=0; j< blockmembers.at(b).size(); j++)
      np.push_back(blockmembers.at(b).at(j));
  }
  //printf("** got np:%i\n", (int) np.size());
  return np;
}
/********************************************************************************************/
bool Pp::is_blocked() {
  return use_blocking;
}
/********************************************************************************************/
bool Pp::block_neighbours(int *i, int *j) {
  if(!use_blocking) return true; // if not blocking ~ one block.
  // figure out whether point i is a block-neighbour of point j,
  // i.e. if the block i is in is adjacent or equal to the block where j is.
  int idx1, idx2;
  std::vector<int> a1, a2;
  idx1 = points.at(*i).getId();
  idx2 = points.at(*j).getId();
  if(idx1 == idx2) return true;
  int a, ab, dz=0, k1=0, k2=0;
  a = blocks.at(0);
  ab = blocks.at(0)*blocks.at(1);
  if(dimension==3){
    k1 = idx1 / ab;
    k2 = idx2 / ab;
    dz = abs(k1-k2);
    if(tor) dz = min( abs(blocks.at(2)-dz), dz);
    if(dz > 1) return false;
  }
  int j1 = (idx1-k1*ab) / a;
  int j2 = (idx2-k2*ab) / a;
  int dy = abs(j1-j2);
  if(tor) dy = min( abs(blocks.at(1)-dy), dy);
  if( dy > 1) return false;
  
  int i1 = idx1 - k1*ab - j1*a;
  int i2 = idx2 - k2*ab - j2*a;
  int dx = abs(i1-i2);
  if(tor) dx = min( abs(blocks.at(0)-dx), dx);
  if( dx > 1) return false;
  return true;
};
  

int min(int a, int b) {
  if(a < b) return a;
  return b;
}

int max(int a, int b) {
  if(a > b) return a;
  return b;
}


/********************************************************************************************/
double Pp::getDist(int *i, int *j) {
  return (this->*dist)(i,j);
}
/********************************************************************************************/
int Pp::Rneighbours(int i, double R) {
  int j, S=0;
  if(use_blocking){
    std::vector<int> neighs = block_neighbours(&i);
    //printf("* computing S\n");
    for(j=0; j < neighs.size(); j++) {
      if(neighs.at(j) != i) {
          if( (this->*distLT)(&i, &neighs.at(j), &R) ) S += 1;
      }
    }
  }
  else{
    for(j=0; j < points.size(); j++) {
      if(j != i) {
          if( (this->*distLT)(&i, &j, &R) ) S += 1;
      }
    }
  }
  return S;
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
/********************************************************************************************/
bool Pp::distLessThan2D(int *i, int *j, double *R) {
   double R2=(*R)*(*R);
   double dx = pow(points.at(*i).getX() - points.at(*j).getX(), 2);
   if(dx > R2) return false;
   double dy = pow(points.at(*i).getY() - points.at(*j).getY(), 2);
   if(dy > R2) return false;
   if( dx + dy > R2 ) return false;
   return true;
}

bool Pp::distLessThan3D(int *i, int *j, double *R) {
   double R2=(*R)*(*R);
   double dx = pow(points.at(*i).getX() - points.at(*j).getX(), 2);
   if(dx > R2) return false;
   double dy = pow(points.at(*i).getY() - points.at(*j).getY(), 2);
   if(dy > R2) return false;
   double dz = pow(points.at(*i).getZ() - points.at(*j).getZ(), 2);
   if(dz > R2) return false;
   if( dx + dy + dz > R2 ) return false;
   return true;
}
bool Pp::distLessThanToroidal2D(int *i, int *j, double *R) {
   double R2=(*R)*(*R);
   double dx = pow( fminf( window[1]-window[0]-fabs(points.at(*i).getX()-points.at(*j).getX()) , fabs(points.at(*i).getX()-points.at(*j).getX()) ) ,2);
   if(dx > R2) return false;
   double dy = pow( fminf( window[3]-window[2]-fabs(points.at(*i).getY()-points.at(*j).getY()) , fabs(points.at(*i).getY()-points.at(*j).getY()) ) ,2);
   if(dy > R2) return false;
   if( dx + dy > R2 ) return false;
   return true;
}
bool Pp::distLessThanToroidal3D(int *i, int *j, double *R) {
   double R2=(*R)*(*R);
   double dx = pow( fminf( window[1]-window[0]-fabs(points.at(*i).getX()-points.at(*j).getX()) , fabs(points.at(*i).getX()-points.at(*j).getX()) ) ,2);
   if(dx > R2) return false;
   double dy = pow( fminf( window[3]-window[2]-fabs(points.at(*i).getY()-points.at(*j).getY()) , fabs(points.at(*i).getY()-points.at(*j).getY()) ) ,2);
   if(dy > R2) return false;
   double dz = pow( fminf( window[5]-window[4]-fabs(points.at(*i).getZ()-points.at(*j).getZ()) , fabs(points.at(*i).getZ()-points.at(*j).getZ()) ) ,2);
   if(dz > R2) return false;
   if( dx + dy + dz > R2 ) return false;
   return true;
}


