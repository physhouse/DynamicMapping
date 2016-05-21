#include "neighbor.h"
#include "fg_atoms.h"
#include "cg_sites.h"
#include "mapping.h"
#include "pointers.h"
#include <cstdio>
#include <cstdlib>

Neighbor::Neighbor(Mapping* map) : Pointers(map) {}
Neighbor::~Neighbor() { cleanup(); }

void Neighbor::init()
{
   cg_num = cg_sites->cg_num;
   maxNeighbors = fg_atoms->fg_num / 4;

   numNeighbors = new int[cg_num];
   memset(numNeighbors, 0, sizeof(int) * cg_num);
   list = new int*[cg_num];
   for (int i=0; i<cg_num; i++)
   {
      list[i] = new int[maxNeighbors];
      memset(list[i], 0, sizeof(int) * maxNeighbors);
   }

   numFgNeighbors = new int[fg_atoms->fg_num];
   memset(numFgNeighbors, 0, sizeof(int) * fg_atoms->fg_num);
   fgList = new int*[fg_atoms->fg_num];
   for (int i=0; i<fg_atoms->fg_num; i++)
   {
      fgList[i] = new int[cg_num];
      memset(fgList[i], 0, sizeof(int) * cg_num);
   }

   rcut = 3.0 * cg_sites->rcut;
   double L = fg_atoms->L;
   dimCell = L / rcut;
   numCells = dimCell * dimCell * dimCell;
   cellLength = L / (double)dimCell;

   cellList = new int[fg_atoms->fg_num];
   cellHead = new int[numCells];

   buildCellList();
   buildNeighborList();
   printf("ncells: %d\n", dimCell);
}

void Neighbor::cleanup()
{
   for (int i=0; i<cg_num; i++)
     delete[] list[i];
  
   delete[] numNeighbors;
   delete[] list;
   delete[] cellList;
   delete[] cellHead;
}

/* Building Neighbor List From Cell Partition */
void Neighbor::buildCellList()
{
   for (int i=0; i<numCells; i++)
   {
      cellHead[i] = -1; //empty the head list
   }

   for (int i=0; i<fg_atoms->fg_num; i++)
   {
      int cellIndex = atom2Cell(fg_atoms->r[i]);
      cellList[i] = cellHead[cellIndex];
      cellHead[cellIndex] = i;
   }

}

/* scanning the stencil of cells to generate NN */
void Neighbor::buildNeighborList()
{
   memset(numNeighbors, 0, sizeof(int) * cg_num);
   for (int i=0; i<cg_num; i++)
   {
      memset(list[i], 0, sizeof(int) * maxNeighbors);
   }

   memset(numFgNeighbors, 0, sizeof(int) * fg_atoms->fg_num);
   for (int i=0; i<fg_atoms->fg_num; i++)
   {
      memset(fgList[i], 0, sizeof(int) * cg_num);
   }

   // Build Double Neighbor List from Cell List
   for (int i=0; i<cg_num; i++)
   {
      int cellIndex = atom2Cell(cg_sites->R[i]);
      int x = cellIndex / (dimCell * dimCell);
      int y = (cellIndex - x * dimCell * dimCell) / dimCell;
      int z = cellIndex % dimCell;
 
      int count = 0;
  
      for (int xoff = -1; xoff <= 1; xoff++)
      {
	for (int yoff = -1; yoff <= 1; yoff++)
	{
	   for (int zoff = -1; zoff <= 1; zoff++)
	   {
	      int px = (x + xoff + dimCell) % dimCell;
	      int py = (y + yoff + dimCell) % dimCell;
	      int pz = (z + zoff + dimCell) % dimCell;
	      int proximity = px * dimCell * dimCell + py * dimCell + pz;
	    
    	      int j = cellHead[proximity];
	      while (j != -1)
	      {
		 if (inRange(cg_sites->R[i], fg_atoms->r[j]))
		 {
	           list[i][count] = j;
		   count++;
		   fgList[j][numFgNeighbors[j]] = i;
		   numFgNeighbors[j]++;
		 }
		
		 j = cellList[j];
	      }
	   }
	}	
      }

      numNeighbors[i] = count;
      //printf("size:%d %d\n", i, count);
      /*printf("%d: %d (%d,%d,%d)", i, dimCell, x, y, z);
      for (int k=0; k<count; k++)
        printf("%d ", list[i][k]);
      printf("\n");*/
   }
}

// Helper functions
int Neighbor::atom2Cell(double* r)
{
   int x = r[0] / cellLength;
   int y = r[1] / cellLength;
   int z = r[2] / cellLength;
   return x * dimCell * dimCell + y * dimCell + z;
}

bool Neighbor::inRange(double* R, double* r)
{
   double rsq = 0.0;
   double L = fg_atoms->L;
   double dx = R[0] - r[0];
   if (dx > 0.5*L) dx -= L;
   else if (dx < -0.5*L) dx += L;
   
   double dy = R[1] - r[1];
   if (dy > 0.5*L) dy -= L;
   else if (dy < -0.5*L) dy += L;

   double dz = R[2] - r[2];
   if (dz > 0.5*L) dz -= L;
   else if (dz < -0.5*L) dz += L;

   rsq = dx * dx + dy * dy + dz * dz;
   if (rsq < rcut * rcut)
      return true;
   else
      return false;
}
