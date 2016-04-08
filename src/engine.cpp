#include "engine.h"
#include "pointers.h"
#include "mapping.h"
#include "matrix_C.h"
#include "matrix_M.h"
#include "matrix_N.h"
#include "cg_sites.h"
#include "fg_atoms.h"
#include "lapacke.h"
#include "cblas.h"

extern void cblas_dgemv (const CBLAS_ORDER layout, const CBLAS_TRANSPOSE TransA, 
			 const int M, const int N, const double alpha, const double *A, const int lda, 
			const double *X, const int incX, const double beta, double *Y, const int incY);

Engine::Engine(Mapping* map) : Pointers(map) {init();}
Engine::~Engine() {cleanup();}

//Initilization
void Engine::init()
{
  nsteps = fg_atoms->nframes;
  fg_num = fg_atoms->fg_num;
  cg_num = cg_sites->cg_num;

  int size_cg = 3 * cg_num;
  int size_fg = 3 * fg_num;

  MCG = new double* [size_cg];
  MFG = new double* [size_cg];

  for (int i = 0; i < size_cg; i++)
  {
    MCG[i] = new double[size_cg];
    MFG[i] = new double[size_fg];
  }

  memset(MCG, 0, sizeof(double) * size_cg * size_cg);
  memset(MFG, 0, sizeof(double) * size_cg * size_fg);

  //Filling in the matrices and vectors
  /*double**** 	M = matrix_M->M;
  double**** 	N = matrix_N->N;
  double**   	C = matrix_C->C;

  //Building the MFG and MCG matrices
  for (int idim = 0; idim < 3; idim++)
  {
    int offset_d1 = idim * size_cg;
    for (int i=offset_d1; i<size_cg; j++)
    {
      //MCG = 1 - M
      for (int jdim = 0; jdim < 3; jdim++)
      {
        int offset_d2 = jdim * size_cg;
        for (int j=offset_d2; j<size_cg; j++) 
        {
           double deltaij = (i==j)? 1.0 : 0.0;
	   MCG[i][j] = deltaij - M[idim][jdim][i - offset_d1][j - offset_d2];
        }
      }

      //MFG = C + V
      for (int jdim = 0; jdim < 3; jdim++)
      {
        int offset_d2 = jdim * size_fg;
        for (int j=offset_d2; j<size_fg; j++) 
	{
	  double Cij = (idim == jdim)? C[i-offset_d1][j-offset_d2] : 0.0;
	  MFG[i][j] = Cij + N[idim][jdim][i - offset_d1][j - offset_d2];
        }
      }
    }
  }

  matrixSolver();*/
}

void Engine::cleanup()
{
  int size_cg = 3 * cg_num;

  for (int i=0; i<size_cg; i++)
  {
     delete[] MCG[i];
     delete[] MFG[i];
  }
  delete[] MCG;
  delete[] MFG;
}

//Update Velocity for one frame
void Engine::update()
{
   initFrame();
   matrixSolver();
   integrate();
   endOfFrame();
}

void Engine::exec()
{
   for (int i=0; i<nsteps-1; i++)
   {
      update();
   }
   // The last frame, not loading
   initFrame();
   matrixSolver();
   integrate();
   cg_sites->output();
   fg_atoms->finishReading();
}

void Engine::initFrame()
{
  matrix_C->compute();
  matrix_M->compute();
  matrix_N->compute();

  int size_cg = 3 * cg_num;
  int size_fg = 3 * fg_num;

  memset(MCG, 0, sizeof(double) * size_cg * size_cg);
  memset(MCG, 0, sizeof(double) * size_cg * size_fg);

  //Filling in the matrices and vectors
  double**** 	M = matrix_M->M;
  double**** 	N = matrix_N->N;
  double**   	C = matrix_C->C;

  //Building the MFG and MCG matrices
  for (int idim = 0; idim < 3; idim++)
  {
    int offset_d1 = idim * size_cg;
    for (int i=offset_d1; i<size_cg; i++)
    {
      //MCG = 1 - M
      for (int jdim = 0; jdim < 3; jdim++)
      {
        int offset_d2 = jdim * size_cg;
        for (int j=offset_d2; j<size_cg; j++) 
        {
           double deltaij = (i==j)? 1.0 : 0.0;
	   MCG[i][j] = deltaij - M[idim][jdim][i - offset_d1][j - offset_d2];
        }
      }

      //MFG = C + V
      for (int jdim = 0; jdim < 3; jdim++)
      {
        int offset_d2 = jdim * size_fg;
        for (int j=offset_d2; j<size_fg; j++) 
	{
	  double Cij = (idim == jdim)? C[i-offset_d1][j-offset_d2] : 0.0;
	  MFG[i][j] = Cij + N[idim][jdim][i - offset_d1][j - offset_d2];
        }
      }
    }
  }
    
}

void Engine::matrixSolver()
{
   double* V_CG = cg_sites->V;
   double* v_fg = fg_atoms->v;
   // Calculate MFG*v_fg
   
   enum CBLAS_ORDER order;
   enum CBLAS_TRANSPOSE transa;

   double alpha, beta;
   int m,n,lda,incx,incy;
  
   order = CblasColMajor;
   transa = CblasNoTrans;

   m = 3 * cg_num;
   n = 3 * fg_num;
   lda = m;
   incx = 1;
   incy = 1;
   alpha = 1.0;
   beta = 0.0;

   double* B = new double[m * n];
   for (int i=0; i<m; i++)
   {
      for (int j=0; j<m; j++)
      {
	B[i * lda + j] = MFG[i][j];
      }
   }

   cblas_dgemv(order, transa, m, n, alpha, B, lda, v_fg, incx, beta,
	       V_CG, incy);
   
   delete[] B;
   //Solving Linear Equations MCG*V=y using lapacke
   int size_a = m * m;
   double* A = new double[size_a];
   for (int i=0; i<m; i++)
   {
      for (int j=0; j<m; j++)
      {
	 A[i * lda + j] = MCG[i][j];
      }
   }

   lda = m;
   int ldb = 1;
   int nrhs = 1;
   int ipiv[m];

   int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, m, nrhs, A, lda, ipiv, V_CG, ldb);
   if (info != 0) exit(-1);
   
   delete[] A;
}

void Engine::integrate()
{
   double** R = cg_sites->R;
   double*  V_CG = cg_sites->V;
   double   dtv = cg_sites->timestep;

   for (int i=0; i<cg_num; i++)   
   {
      R[i][0] += V_CG[i] * dtv;
      R[i][1] += V_CG[i + cg_num] * dtv;
      R[i][2] += V_CG[i + 2*cg_num] * dtv;
   }
}

void Engine::endOfFrame()
{
   cg_sites->output();
   fg_atoms->readNextFrame();
}

