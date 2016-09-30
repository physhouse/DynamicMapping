#include "omp.h"
#include "mapping.h"
#include "pointers.h"
#include "matrix_M.h"
#include "matrix_C.h"
#include "cg_sites.h"
#include "fg_atoms.h"
#include "neighbor.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

Matrix_M::Matrix_M(Mapping *map) : Pointers(map) {}
Matrix_M::~Matrix_M() {cleanup();}

void Matrix_M::init()
{
  fg_num  = fg_atoms->fg_num;
  cg_num  = cg_sites->cg_num;
  int dim = 3;

  // Initialize M
  M = new double***[dim];
  for (int idim =0; idim<dim; idim++)
  {
     M[idim] = new double**[dim];
     for (int jdim=0; jdim<dim; jdim++)
     {
	M[idim][jdim] = new double*[cg_num];
	for (int I=0; I<cg_num; I++)
	{
	  M[idim][jdim][I] = new double[cg_num];
	  for (int J=0; J<cg_num; J++) M[idim][jdim][I][J] = 0.0;
	}
     }
  }
}

double Matrix_M::distance(double* R, double* r)
{
   double L = fg_atoms->L;
   double distance = 0.0;
   for (int dim=0; dim<3; dim++)
   {
      double r_dim = R[dim] - r[dim];
      if (r_dim > 0.5 * L) distance += (r_dim - L) * (r_dim - L);
      else if (r_dim < -0.5 * L) distance += (r_dim + L) * (r_dim + L);
      else distance += r_dim * r_dim;
   }
   return distance;
}
// Generating the matrix dc/dR

void Matrix_M::generate_M()
{
   double**   r     = fg_atoms->r;
   double**   R     = cg_sites->R;
   double**   C     = matrix_C->C;
   double***  dw    = matrix_C->dw;
   double**   W     = matrix_C->W;
   double*    w_sum = matrix_C->w_sum;
   int*	      numNeighbors = neighbor->numNeighbors;
   int**      list  = neighbor->list;

   for (int I=0; I<cg_num; I++)
   {
      #pragma omp parallel for
      for (int J=0; J<cg_num; J++)
      {
	for (int jdim = 0; jdim < 3; jdim++)
	{
	  M[jdim][0][I][J] = 0.0;
	  M[jdim][1][I][J] = 0.0;
	  M[jdim][2][I][J] = 0.0;
	}
	
	if (I != J)
	{
	  for (int ii=0; ii<numNeighbors[I]; ii++)
	  {
	    int i = list[I][ii];
	    double dcdx = -dw[J][i][0] / w_sum[i];
	    double dcdy = -dw[J][i][1] / w_sum[i];
	    double dcdz = -dw[J][i][2] / w_sum[i];

	    double sumx = 0.0, sumy = 0.0, sumz = 0.0;

	    for (int index=0; index<numNeighbors[I]; index++)
	    {
	      int j = list[I][index];
	      sumx += C[I][j] * dw[J][j][0] / w_sum[j];
	      sumy += C[I][j] * dw[J][j][1] / w_sum[j];
	      sumz += C[I][j] * dw[J][j][2] / w_sum[j];
	    }

	    dcdx += sumx;
	    dcdy += sumy;
	    dcdz += sumz;

	    for (int jdim = 0; jdim < 3; jdim++)
	    {
	      M[jdim][0][I][J] += C[I][i] * dcdx * (r[i][jdim] - R[I][jdim]);
	      M[jdim][1][I][J] += C[I][i] * dcdy * (r[i][jdim] - R[I][jdim]);
	      M[jdim][2][I][J] += C[I][i] * dcdz * (r[i][jdim] - R[I][jdim]);
	      //M[jdim][0][I][J] += C[I][i] * dcdx * r[i][jdim];
	      //M[jdim][1][I][J] += C[I][i] * dcdy * r[i][jdim];
	      //M[jdim][2][I][J] += C[I][i] * dcdz * r[i][jdim];
	    }	
	  } 
	}
	else
	{
	  for (int ii=0; ii<numNeighbors[I]; ii++)
	  {
	    int i = list[I][ii];
	    double dcdx, dcdy, dcdz;

	    dcdx = dw[I][i][0] / W[I][i] - dw[I][i][0] / w_sum[i];
	    dcdy = dw[I][i][1] / W[I][i] - dw[I][i][1] / w_sum[i];
	    dcdz = dw[I][i][2] / W[I][i] - dw[I][i][2] / w_sum[i];

	    //printf("%e %e %e\n", dcdx,dcdy,dcdz);
	    //printf("%e %e %e\n", dw[I][i][0], W[I][i], w_sum[i]);
	    double sumx = 0.0, sumy = 0.0, sumz = 0.0;

	    for (int index=0; index<numNeighbors[I]; index++)
	    {
	      int j = list[I][index];

	      sumx += C[I][j] * (dw[I][j][0] / W[I][j] - dw[I][j][0] / w_sum[j]);
	      sumy += C[I][j] * (dw[I][j][1] / W[I][j] - dw[I][j][1] / w_sum[j]);
	      sumz += C[I][j] * (dw[I][j][2] / W[I][j] - dw[I][j][2] / w_sum[j]);
	    }

	    dcdx -= sumx;
	    dcdy -= sumy;
	    dcdz -= sumz;

	    //printf("%e %e %e\n", dcdx,dcdy,dcdz);

	    for (int jdim = 0; jdim < 3; jdim++)
	    {
	      M[jdim][0][I][J] += C[I][i] * dcdx * (r[i][jdim] - R[I][jdim]);
	      M[jdim][1][I][J] += C[I][i] * dcdy * (r[i][jdim] - R[I][jdim]);
	      M[jdim][2][I][J] += C[I][i] * dcdz * (r[i][jdim] - R[I][jdim]);
	      //M[jdim][0][I][J] += C[I][i] * dcdx * r[i][jdim];
	      //M[jdim][1][I][J] += C[I][i] * dcdy * r[i][jdim];
	      //M[jdim][2][I][J] += C[I][i] * dcdz * r[i][jdim];
	    }	
	  } 
	    //printf("%f %f %f\n", M[0][0][I][J], M[0][1][I][J], M[0][2][I][J]);
	}
      }
   } 
}

void Matrix_M::cleanup()
{
   for (int idim=0; idim<3; idim++)
   {
     for (int jdim=0; jdim<3; jdim++)
     {
	for (int I=0; I<cg_num; I++) delete[] M[idim][jdim][I];
	delete[] M[idim][jdim];
     }
     delete[] M[idim];
   }
   delete[] M;
   printf("cleaning up M...\n");
}

void Matrix_M::compute()
{
   generate_M();
}
