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
   double     L     = fg_atoms->L;
   double     rcut  = cg_sites->rcut;
   int*	      numNeighbors = neighbor->numNeighbors;
   int**      list  = neighbor->list;

   for (int I=0; I<cg_num; I++)
   {
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

	    for (int j=0; j<numNeighbors[I]; j++)
	    {
	      int index = list[I][j];
	      sumx += C[I][index] * dw[J][index][0] / w_sum[index];
	      sumx += C[I][index] * dw[J][index][1] / w_sum[index];
	      sumx += C[I][index] * dw[J][index][2] / w_sum[index];
	    }

	    dcdx += sumx;
	    dcdy += sumy;
	    dcdz += sumz;

	    for (int jdim = 0; jdim < 3; jdim++)
	    {
	      M[jdim][0][I][J] += C[I][i] * dcdx * r[i][jdim];
	      M[jdim][1][I][J] += C[I][i] * dcdy * r[i][jdim];
	      M[jdim][2][I][J] += C[I][i] * dcdz * r[i][jdim];
	    }	
	  } 
	}
	else
	{
	  for (int ii=0; ii<numNeighbors[I]; ii++)
	  {
	    int i = list[I][ii];
	    double dcdx, dcdy, dcdz;

 	    double rc2 = rcut * rcut;
	    //double deriv = 1.0 + tanh(r_Ii - rcut);
	    double rdim = R[I][0] - r[i][0];
	    if (rdim > 0.5 * L) rdim -= L;
	    else if (rdim < -0.5 * L) rdim += L;
	    dcdx = -2.0 * rdim / rc2 - dw[I][i][0] / w_sum[i];

	    rdim = R[I][1] - r[i][1];
	    if (rdim > 0.5 * L) rdim -= L;
	    else if (rdim < -0.5 * L) rdim += L;
	    dcdy = -2.0 * rdim / rc2 - dw[I][i][1] / w_sum[i];

	    rdim = R[I][2] - r[i][2];
	    if (rdim > 0.5 * L) rdim -= L;
	    else if (rdim < -0.5 * L) rdim += L;
	    dcdz = -2.0 * rdim / rc2 - dw[I][i][2] / w_sum[i];

	    double sumx = 0.0, sumy = 0.0, sumz = 0.0;

	    for (int index=0; index<numNeighbors[I]; index++)
	    {
	      int j = list[I][index];

	      sumx += dw[I][j][0] / W[I][j] - dw[I][j][0] / w_sum[j];
	      sumy += dw[I][j][1] / W[I][j] - dw[I][j][1] / w_sum[j];
	      sumz += dw[I][j][2] / W[I][j] - dw[I][j][2] / w_sum[j];
	    }

	    dcdx -= sumx;
	    dcdz -= sumy;
	    dcdy -= sumz;

	    for (int jdim = 0; jdim < 3; jdim++)
	    {
	      M[jdim][0][I][J] += C[I][i] * dcdx * r[i][jdim];
	      M[jdim][1][I][J] += C[I][i] * dcdy * r[i][jdim];
	      M[jdim][2][I][J] += C[I][i] * dcdz * r[i][jdim];
	    }	
	  } 
	    //printf("%f %f %f\n", M[0][0][I][J], M[0][1][I][J], M[0][2][I][J]);
	
	}
      }
   } 
  
   //test
   /*printf("Testing M...\n");
   for (int i=0; i<cg_num; i++)
   {
      for (int j=0; j<cg_num; j++)
      {
	printf("%12.5lf", M[0][0][i][j]);
      }
      printf("\n");
   }*/

   /*printf("Testing Mapping...\n");
   for (int i=0; i<fg_num; i++)
       printf("%12.8lf\n", w_sum[i]);*/

   /*printf("Testing Mass...\n");
   double* sum_c = new double[fg_num];
   for (int i=0; i<fg_num; i++)
   {
      sum_c[i] = 0.0;
      for (int I=0; I<cg_num; I++)
	sum_c[i] += C[I][i];
   }

   double m_total = 0.0;
   for (int i=0; i<cg_num; i++)
   {
      double mass = 0.0;
      for (int j=0; j<fg_num; j++)
      {
	 mass += C[i][j] / sum_c[j];
      }
      m_total += mass;
      printf("%12.5lf\n", mass);
   } 
  
      printf("TOTAL:  %12.5lf\n", m_total);
   delete[] sum_c;*/
   /*printf("Testing C...\n");
   for (int i=0; i<fg_num; i++)
   {
      for (int j=0; j<cg_num; j++)
      {
	printf("%12.5lf", C[j][i]);
      }
      printf("\n");
   }*/
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
