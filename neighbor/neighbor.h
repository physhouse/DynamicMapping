#ifndef _NEIGHBOR_H_
#define _NEIGHBOR_H_

// Full Neighbor List for fast evaluation of C, M, N matrices
#include "pointers.h"

class Neighbor : protected Pointers
{
public:
  Neighbor(class Mapping* map);
  ~Neighbor();
  void     init();
  void     cleanup();
  void     buildCellList();
  void     buildNeighborList();
  int      atom2Cell(double* r);
  bool     inRange(double* R, double* r);

  int 	   cg_num;
  double   rcut;
  int*     numNeighbors;
  int**    list;	// Neighbor List
  int*	   numFgNeighbors;
  int**    fgList;
  int*     cellHead;	// Linked Cell List
  int*	   cellList;	// Linked Cell List
  int	   maxNeighbors;
  int	   dimCell;
  int	   numCells;
  double   cellLength;
};

#endif
