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
   double weight(double r);
   void weight_deriv(double* R, double* r, double* dwvec);
 
   void matrixGenerator();
   void cleanup();
   void compute();

   double**	W;	// cg_num by fg_num weights Matrix W
   double*	w_sum;	// fg vector, W_sum[j] = sum over J W[J][j]
   double***	dw;	// dw/dx matrix, dimension fg by cg by 3
   double**	dw_sum;	// dw/dx sum matrix, dimension fg by 3
   double**	C;	// C matrix, dimension cg*fg 
   int		cg_num;
   int		fg_num;
   double	L;
   double	rcut;
};

#endif
