#include "neighbor.h"
#include "fg_atoms.h"
#include "cg_sites.h"
#include "mapping.h"
#include "pointers.h"
#include <cstdio>
#include <cstdlib>

Neighbor::Neighbor(Mapping *map) : Pointers(map) {}
Neighbor::~Neighbor()
{
    cleanup();
}

void Neighbor::init()
{
    cg_num = cg_sites->cg_num;
    maxNeighbors = fg_atoms->fg_num;

    numNeighbors = new int[cg_num];
    memset(numNeighbors, 0, sizeof(int) * cg_num);
    list = new int *[cg_num];
    for (int i = 0; i < cg_num; i++)
    {
        list[i] = new int[maxNeighbors];
        memset(list[i], 0, sizeof(int) * maxNeighbors);
    }

    numFgNeighbors = new int[fg_atoms->fg_num];
    memset(numFgNeighbors, 0, sizeof(int) * fg_atoms->fg_num);
    fgList = new int *[fg_atoms->fg_num];
    for (int i = 0; i < fg_atoms->fg_num; i++)
    {
        fgList[i] = new int[cg_num];
        memset(fgList[i], 0, sizeof(int) * cg_num);
    }

    rcut = 2.0 * cg_sites->rcut;
    double L = fg_atoms->L;
    dimCell = L / rcut;
    numCells = dimCell * dimCell * dimCell;
    cellLength = L / (double)dimCell;

    cellList = new int[fg_atoms->fg_num];
    cellHead = new int[numCells];
    memset(cellList, 0, sizeof(int) * fg_atoms->fg_num);
    memset(cellHead, 0, sizeof(int) * numCells);

    buildCellList();
    buildNeighborList();
    printf("ncells: %d\n", dimCell);
}

void Neighbor::cleanup()
{
    for (int i = 0; i < cg_num; i++)
        delete[] list[i];

    delete[] numNeighbors;
    delete[] list;
    delete[] cellList;
    delete[] cellHead;
}

/* Building Neighbor List From Cell Partition */
void Neighbor::buildCellList()
{
    for (int i = 0; i < numCells; i++)
    {
        cellHead[i] = -1; //empty the head list
    }

    for (int i = 0; i < fg_atoms->fg_num; i++)
    {
        int cellIndex = atom2Cell(fg_atoms->r[i]);
        if (cellIndex > numCells)
            printf("%d : Warning!, positions %12.8lf %12.8lf %12.8lf L - %12.8lf\n", i, fg_atoms->r[i][0], fg_atoms->r[i][1], fg_atoms->r[i][2], fg_atoms->L);
        cellList[i] = cellHead[cellIndex];
        cellHead[cellIndex] = i;
    }
    printf("CellList Finished\n");
}

/* scanning the stencil of cells to generate NN */
void Neighbor::buildNeighborList()
{
    //memset(numNeighbors, 0, sizeof(int) * cg_num);
    for (int i = 0; i < cg_num; i++)
    {
        numNeighbors[i] = 0;
        for (int j = 0; j < maxNeighbors; j++)
            list[i][j] = -1;
        //memset(list[i], 0, sizeof(int) * maxNeighbors);
    }

    //memset(numFgNeighbors, 0, sizeof(int) * fg_atoms->fg_num);
    for (int i = 0; i < fg_atoms->fg_num; i++)
    {
        numFgNeighbors[i] = 0;
        for (int j = 0; j < cg_num; j++)
            fgList[i][j] = -1;
        //memset(fgList[i], 0, sizeof(int) * cg_num);
    }

    // Build Double Neighbor List from Cell List
    int sumg = 0;
    for (int i = 0; i < cg_num; i++)
    {
        int cellIndex = atom2Cell(cg_sites->R[i]);
        int x = cellIndex / (dimCell * dimCell);
        int y = (cellIndex - x * dimCell * dimCell) / dimCell;
        int z = cellIndex % dimCell;

        int count = 0;

        if (dimCell >= 3)
        {
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

                        if (proximity >= 0 && proximity < numCells)
                        {
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
            }

            numNeighbors[i] = count;
            sumg += count;
            /*printf("size:%d %d\n", i, count);
            printf("%d: %d (%d,%d,%d)", i, dimCell, x, y, z);
            for (int k=0; k<count; k++)
              printf("%d ", list[i][k]);
            printf("\n");*/
        }
        else
        {
            for (int j = 0; j < fg_atoms->fg_num; j++)
            {
                if (inRange(cg_sites->R[i], fg_atoms->r[j]))
                {
                    list[i][count] = j;
                    count++;
                    fgList[j][numFgNeighbors[j]] = i;
                    numFgNeighbors[j]++;
                }
            }

            numNeighbors[i] = count;
            sumg += count;
        }
    }

    int sumf = 0;
    for (int i = 0; i < fg_atoms->fg_num; i++)
    {
        sumf += numFgNeighbors[i];
        /*printf("Atom %d : ", i+1);
        for (int j=0; j<numFgNeighbors[i]; j++)
        printf("%d ", fgList[i][j]);
            printf("\n");*/
    }
}

// Helper functions
int Neighbor::atom2Cell(double *r)
{
    double L = fg_atoms->L;
    double eps = 1e-5;
    if (r[0] + eps > L) r[0] -= L;
    else if (r[0] + eps < 0) r[0] += L;

    if (r[1] + eps > L) r[1] -= L;
    else if (r[1] + eps < 0) r[1] += L;

    if (r[2] + eps > L) r[2] -= L;
    else if (r[2] + eps < 0) r[2] += L;

    int x = r[0] / cellLength;
    int y = r[1] / cellLength;
    int z = r[2] / cellLength;
    return x * dimCell * dimCell + y * dimCell + z;
}

bool Neighbor::inRange(double *R, double *r)
{
    double rsq = 0.0;
    double L = fg_atoms->L;
    double dx = R[0] - r[0];
    if (dx > 0.5 * L) dx -= L;
    else if (dx < -0.5 * L) dx += L;

    double dy = R[1] - r[1];
    if (dy > 0.5 * L) dy -= L;
    else if (dy < -0.5 * L) dy += L;

    double dz = R[2] - r[2];
    if (dz > 0.5 * L) dz -= L;
    else if (dz < -0.5 * L) dz += L;

    rsq = dx * dx + dy * dy + dz * dz;
    if (rsq < rcut * rcut)
        return true;
    else
        return false;
}
