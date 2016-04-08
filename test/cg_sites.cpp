#include "cg_sites.h"
#include "mapping.h"
#include "fg_atoms.h"
#include "matrix_C.h"
#include "pointers.h"
#include "comm.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>

Cg_sites::Cg_sites(Mapping *map) : Pointers(map) {}
Cg_sites::~Cg_sites() {cleanup();}

void Cg_sites::init()
{
  VAR_BEGIN
    GET_INT(cg_num)
    GET_REAL(tol)
    GET_REAL(timestep)
    GET_INT(niter)
  VAR_END

  R = new double*[cg_num];
  for (int i=0; i<cg_num; i++)
  {
    R[i] = new double[3];
    for (int j=0; j<3; j++) R[i][j] = 0.0;
  }

  V = new double[3*cg_num];
  memset(V, 0, sizeof(double) * 3 * cg_num);

  L = fg_atoms->L;
  int nperdim = (int)cbrt((double)cg_num);
  rcut = 0.5 * L / nperdim;

  printf("fgnum = %d cgnum = %d L = %12.8lf rcut = %12.8lf\n", fg_atoms->fg_num, cg_num, L, rcut);
  //Grid Initialization of CG particle positions
  double size_block = L / nperdim;
  for (int i=0; i<nperdim; i++)
    for (int j=0; j<nperdim; j++)
      for (int k=0; k<nperdim; k++)
      {
	 int ind = nperdim * nperdim * i + nperdim * j + k;
	 R[ind][0] = (i+0.5) * size_block;
	 R[ind][1] = (j+0.5) * size_block;
	 R[ind][2] = (k+0.5) * size_block;
      }

  cgtrj.open("CG_TRJ.lmpstrj", std::ofstream::out);
  cgtrj<<"CG_TRJ : Dynamical Mapping Coarse Graining"<<std::endl;
}

void Cg_sites::cleanup()
{
  cgtrj.close();
  for (int i=0; i<cg_num; i++)
  {
     delete[] R[i];
  }
  delete[] R;
  delete[] V;
}

void Cg_sites::output()
{
  //header
  cgtrj<<"ITEM: TIMESTEP"<<std::endl; 
  cgtrj<<fg_atoms->currentStep<<std::endl;
  cgtrj<<"ITEM: NUMBER OF ATOMS"<<std::endl; 
  cgtrj<<cg_num<<std::endl;
  cgtrj<<"ITEM: BOX BOUNDS pp pp pp"<<std::endl; 
  for (int i=0; i<3; i++) cgtrj<<"0 "<<fg_atoms->L<<std::endl;
  cgtrj<<"ITEM: ATOMS id type x y z vx vy vz"<<std::endl; 

  //body
  for (int i=0; i<cg_num; i++)
  {
    cgtrj<<i+1<<' '<<1<<' '<<R[i][0]<<' '<<R[i][1]<<' '<<R[i][2]<<' '<<V[i]<<' '<<V[i + cg_num]<<' '<<V[i + 2*cg_num]<<std::endl;
  }
 
}

void Cg_sites::firstMapping()
{
  double error = 0.0;
  double** mapMatrix = matrix_C->C;
  double** r = fg_atoms->r;
  int fg_num = fg_atoms->fg_num;

  for (int round = 0; round<niter; round++)
  {
    error = 0.0;
    for (int i=0; i<cg_num; i++)
    {
       double ri0 = R[i][0];
       double ri1 = R[i][1];
       double ri2 = R[i][2];
       R[i][0] = 0.0; R[i][1] = 0.0; R[i][2] = 0.0;
       for (int j=0; j<fg_num; j++)
       {
	  R[i][0] += mapMatrix[i][j] * r[j][0];
	  R[i][1] += mapMatrix[i][j] * r[j][1];
	  R[i][2] += mapMatrix[i][j] * r[j][2];
       }
       //printf("R = %12.8lf, before = %12.8lf\n", R[i][0], ri0);
       error += (R[i][0] - ri0)*(R[i][0] -ri0) + (R[i][1] - ri1)*(R[i][1] - ri1) + (R[i][2] - ri2)*(R[i][2] - ri2);
    }
    
    if (error < tol) return;
    else
      matrix_C->compute();
    if (round%10 == 0) printf("Round %d, error = %e\n", round, error);
  }
}

