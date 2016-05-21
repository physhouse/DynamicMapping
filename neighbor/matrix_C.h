#ifndef _MATRIX_H
#define _MATRIX_H

#include "pointers.h"

class Matrix_C : protected Pointers
{
public:
   Matrix_C(class Mapping *);
   ~Matrix_C();

   void init();
   double w_Ij(double* R, double* r);
   double weightsVarying(double* R, double* r, double M);
   double weight(double r);
   double massWeights(double M);
   void weight_deriv(double* R, double* r, double* dwvec);
   void weightsVarying_deriv(double* R, double* r, double* dwvec, double M, double alpha);
 
   void matrixGenerator();
   double generateWeights();
   void generateMass();
   void generateMatrixC();
   void iterateSolver();
   void cleanup();
   void compute();

   double**	W;	// cg_num by fg_num weights Matrix W
   double*	w_sum;	// fg vector, W_sum[j] = sum over J W[J][j]
   double***	dw;	// dw/dx matrix, dimension fg by cg by 3
   double**	dw_sum;	// dw/dx sum matrix, dimension fg by 3
   double**	C;	// C matrix, dimension cg*fg 
   double*      sumC;   // sumC[j] = sum over i, C[i][j]
   double*      F;   // F[I]
   int		cg_num;
   int		fg_num;
   double	L;
   double	rcut;

   double	width;  // Width of Gaussian Distribution
   double	M0;	// Constraint Mass
};

#endif
