//2016-04-04
#include "mapping.h"
#include "pointers.h"
#include "matrix_N.h"
#include "matrix_C.h"
#include "cg_sites.h"
#include "fg_atoms.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

// Generating the N matrix, no need to cache the dcdr matrices

Matrix_N::Matrix_N(Mapping *map) : Pointers(map) {}
Matrix_N::~Matrix_N() {cleanup();}

void Matrix_N::init()
{
  fg_num  = fg_atoms->fg_num;
  cg_num  = cg_sites->cg_num;

  // Initialize M
  N = new double***[3];
  for (int idim =0; idim<3; idim++)
  {
     N[idim] = new double**[3];
     for (int jdim=0; jdim<3; jdim++)
     {
	N[idim][jdim] = new double*[cg_num];
	for (int I=0; I<cg_num; I++)
	{
	  N[idim][jdim][I] = new double[fg_num];
	  for (int j=0; j<fg_num; j++) N[idim][jdim][I][j] = 0.0;
	}
     }
  }
}

// Generating the matrix dc/dR
double Matrix_N::distance(double* R, double* r)
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

void Matrix_N::generate_N()
{
   double**   r      = fg_atoms->r;
   double**   R      = cg_sites->R;
   double**   C      = matrix_C->C;
   double***  dw     = matrix_C->dw;
   double**   W      = matrix_C->W;
   double*    w_sum  = matrix_C->w_sum;
   double**   dw_sum = matrix_C->dw_sum;
   double     L	     = fg_atoms->L;

   // Notation dcdr[dim][I][j][i] = /partial C_I_i /partial r_dim_j
   for (int I=0; I<cg_num; I++)
   {
      double sum_denom = 0.0;
      for (int m=0; m<fg_num; m++) sum_denom += W[I][m] / w_sum[m];

      for (int j=0; j<fg_num; j++)
      {
	for (int idim = 0; idim < 3; idim++)
	  for (int jdim = 0; jdim < 3; jdim++)
	    N[idim][jdim][I][j] = 0.0;
	// Note:: dC_I_i/dr_j is not a function of i!
	double dcdx = - (W[I][j] * dw_sum[j][0] - dw[I][j][0] * w_sum[j])/(w_sum[j] * w_sum[j] * sum_denom);
	double dcdy = - (W[I][j] * dw_sum[j][1] - dw[I][j][1] * w_sum[j])/(w_sum[j] * w_sum[j] * sum_denom);
	double dcdz = - (W[I][j] * dw_sum[j][2] - dw[I][j][2] * w_sum[j])/(w_sum[j] * w_sum[j] * sum_denom);

	for (int i=0; i<fg_num; i++)
	{
	  if (i != j)
	  {

	    for (int jdim =0; jdim < 3; jdim++)
	    {
	      /* Notation:: N[d1][d2][I][j], d1->dimension of CG particles, d2->dimension of FG particles 
 	       * N[d1][d2][I][j] = sumover(i){r_i_d1 * dC_I_i/dr_j_d2} */
	      N[jdim][0][I][j] +=  C[I][i] * dcdx * r[i][jdim];
	      N[jdim][1][I][j] +=  C[I][i] * dcdy * r[i][jdim];
	      N[jdim][2][I][j] +=  C[I][i] * dcdz * r[i][jdim];
	    }
	  } 
	  else
	  {
	    double dcdx_diagonal, dcdy_diagonal, dcdz_diagonal;

	    if (W[I][i] < 1E-50)
	    {
	      double r_Ii = distance(R[I], r[i]);	
              double rdim = r[i][0] - R[I][0];
              if (rdim > 0.5 * L) rdim -= L; 
              else if (rdim < -0.5 * L) rdim += L;
	      dcdx_diagonal = -2.0 * rdim / r_Ii + dw_sum[i][0] / w_sum[i] + dcdx;

              rdim = r[i][1] - R[I][1];
              if (rdim > 0.5 * L) rdim -= L; 
              else if (rdim < -0.5 * L) rdim += L;
	      dcdy_diagonal = -2.0 * rdim / r_Ii + dw_sum[i][1] / w_sum[i] + dcdy;

              rdim = r[i][2] - R[I][2];
              if (rdim > 0.5 * L) rdim -= L; 
              else if (rdim < -0.5 * L) rdim += L;
	      dcdz_diagonal = -2.0 * rdim / r_Ii + dw_sum[i][2] / w_sum[i] + dcdz;
	    }
            else
	    {	
	      dcdx_diagonal = -dw[I][i][0] / W[I][i] + dw_sum[i][0] / w_sum[i] + dcdx;
	      dcdy_diagonal = -dw[I][i][1] / W[I][i] + dw_sum[i][1] / w_sum[i] + dcdy;
	      dcdz_diagonal = -dw[I][i][2] / W[I][i] + dw_sum[i][2] / w_sum[i] + dcdz;
	    }

	    for (int jdim =0; jdim < 3; jdim++)
	    {
	      /* Notation:: N[d1][d2][I][j], d1->dimension of CG particles, d2->dimension of FG particles 
 	       * N[d1][d2][I][j] = sumover(i){r_i_d1 * dC_I_i/dr_j_d2} */
	      N[jdim][0][I][j] +=  C[I][i] * dcdx_diagonal * r[i][jdim];
	      N[jdim][1][I][j] +=  C[I][i] * dcdy_diagonal * r[i][jdim];
	      N[jdim][2][I][j] +=  C[I][i] * dcdz_diagonal * r[i][jdim];
	    }
	  }
        }
     } 
   }
}

void Matrix_N::compute()
{
   generate_N();
}

void Matrix_N::cleanup()
{
   for (int idim=0; idim<3; idim++)
   {
     for (int jdim=0; jdim<3; jdim++)
     {
	for (int I=0; I<cg_num; I++) 
	{
	   delete[] N[idim][jdim][I];
	}
	delete[] N[idim][jdim];
     }
     delete[] N[idim];
   }
   delete[] N;
   printf("cleaning up N...\n");
}

