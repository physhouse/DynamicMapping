#include "engine.h"
#include "pointers.h"
#include "mapping.h"
#include "matrix_C.h"
#include "matrix_M.h"
#include "matrix_N.h"
#include "cg_sites.h"
#include "fg_atoms.h"
#include "neighbor.h"
#include "lapacke.h"
#include "cblas.h"
#include <ctime>

extern void cblas_dgemv (const CBLAS_ORDER layout, const CBLAS_TRANSPOSE TransA, 
			 const int M, const int N, const double alpha, const double *A, const int lda, 
			const double *X, const int incX, const double beta, double *Y, const int incY);

Engine::Engine(Mapping* map) : Pointers(map) {}
Engine::~Engine() {cleanup();}

//Initilization
void Engine::init(int argc)
{
  if (argc > 3) restartFlag = true;
  else restartFlag = false;

  nsteps = fg_atoms->nframes;
  fg_num = fg_atoms->fg_num;
  cg_num = cg_sites->cg_num;

  int size_cg = 3 * cg_num;
  int size_fg = 3 * fg_num;

  MCG = new double* [size_cg];
  MFG = new double* [size_cg];

  vMap = new double[size_cg * size_fg];
  memset(vMap, 0, sizeof(double) * size_cg * size_fg);

  for (int i = 0; i < size_cg; i++)
  {
    MCG[i] = new double[size_cg];
    memset(MCG[i], 0, sizeof(double) * size_cg);
    MFG[i] = new double[size_fg];
    memset(MFG[i], 0, sizeof(double) * size_fg);
  }

  checkmap.open("CG_check.lmpstrj", std::ofstream::out);
  checkmap<<"CG_CHECK : Checking the Propagation"<<std::endl;

  invMass.open("Mat_miu.dat", std::ofstream::out);
}

void Engine::cleanup()
{
  checkmap.close();
  invMass.close();

  int size_cg = 3 * cg_num;

  for (int i=0; i<size_cg; i++)
  {
     delete[] MCG[i];
     delete[] MFG[i];
  }
  delete[] MCG;
  delete[] MFG;
  delete[] vMap;
}

void Engine::buildNeighbors()
{
   std::clock_t start;
   double cellTime, allTime;
   
   start = std::clock();
   neighbor->buildCellList();
   cellTime = ( std::clock() - start) / (double) CLOCKS_PER_SEC;
   neighbor->buildNeighborList();
   allTime = ( std::clock() - start) / (double) CLOCKS_PER_SEC;
   printf("cell: %12.8lf  neighbor: %12.8lf\n ncells %d\n", cellTime, allTime-cellTime, neighbor->dimCell);
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
   printf("nsteps = %d\n", nsteps);

   if (!restartFlag)
   {
     cg_sites->firstMapping();
     endOfFrame();
   }
   else
   {

   }

   //cg_sites->output();
   for (int i=0; i<nsteps-1; i++)
   {
      update();
      //testing(i);
      if (i % 1 == 0)
      {
        buildNeighbors();
      }
      if (fg_atoms->currentStep % cg_sites->freq == 0)
      {
	printf("Step %d\n", fg_atoms->currentStep);
	//testing(i);
	cg_sites->firstMapping();
      }
   }
   // The last frame, not loading
   initFrame();
   matrixSolver();
   integrate();
   cg_sites->output();
   checker();
   fg_atoms->finishReading();
}

void Engine::testing(int step)
{
   printf("Testing M... Current Step %d\n", step);
   for (int i=0; i<cg_sites->cg_num; i++)
   {
      for (int j=0; j<cg_sites->cg_num; j++)
      {
	printf("%8.8lf\t", matrix_M->M[0][0][i][j]);
      }
      printf("\n");
   }

   printf("Testing N... Current Step %d\n", step);
   for (int i=0; i<cg_sites->cg_num; i++)
   {
      for (int j=0; j<fg_atoms->fg_num; j++)
      {
	printf("%8.8lf\t", matrix_N->N[0][0][i][j]);
      }
      printf("\n");
   }

   printf("Testing C... Current Step %d\n", step);
   for (int i=0; i<cg_sites->cg_num; i++)
   {
      printf("%d : sum = %12.8lf\n", i+1, matrix_C->sumC[i]);
      for (int j=0; j<fg_atoms->fg_num; j++)
      {
	printf("%12.8lf\t", matrix_C->C[i][j]);
      }
      printf("\n");
   }

   printf("Testing dW...\n");
   for (int i=0; i<cg_sites->cg_num; i++)
   {
      printf("%d : \n", i+1);
      for (int j=0; j<fg_atoms->fg_num; j++)
      {
	printf("%e\t", matrix_C->dw[i][j][0]);
      }
      printf("\n");
      for (int j=0; j<fg_atoms->fg_num; j++)
      {
	printf("%e\t", matrix_C->W[i][j]);
      }
      printf("\n");
   }

   printf("Tesing Neighbor...\n");
   int fg2cg = 0; 
   int cg2fg = 0;
   for (int i=0; i<cg_sites->cg_num; i++)
   {
     cg2fg += neighbor->numNeighbors[i];
   }
   for (int i=0; i<fg_atoms->fg_num; i++)
   {
     fg2cg += neighbor->numFgNeighbors[i];
   }
   printf("Neighbors: %d, FGNeighbors: %d\n", cg2fg, fg2cg);

   printf("Testing WSUM...\n");
   for (int i=0; i<fg_atoms->fg_num; i++)
   {
      printf("%d : %e\t", i, matrix_C->w_sum[i]);
   }
   printf("\n");
}

void Engine::initFrame()
{
  matrix_C->compute();
  printf("finishing C\n");
  matrix_M->compute();
  printf("finishing M\n");
  matrix_N->compute();
  printf("finishing N\n");

  //Filling in the matrices and vectors
  double**** 	M = matrix_M->M;
  double**** 	N = matrix_N->N;
  double**   	C = matrix_C->C;
  //Building the MFG and MCG matrices
  for (int idim = 0; idim < 3; idim++)
  {
    int offset_d1 = idim * cg_num;
    for (int i=offset_d1; i < (offset_d1 + cg_num); i++)
    {
      //MCG = 1 - M
      for (int jdim = 0; jdim < 3; jdim++)
      {
        int offset_d2 = jdim * cg_num;
        for (int j=offset_d2; j < (offset_d2 + cg_num); j++) 
        {
           double deltaij = (i==j)? 1.0 : 0.0;
	   //MCG[i][j] = deltaij;
	   MCG[i][j] = deltaij - M[idim][jdim][i - offset_d1][j - offset_d2];
        }
      }
      //MFG = C + V
      for (int jdim = 0; jdim < 3; jdim++)
      {
        int offset_d2 = jdim * fg_num;
        for (int j=offset_d2; j<(offset_d2 + fg_num); j++) 
	{
	  double Cij = (idim == jdim)? C[i-offset_d1][j-offset_d2] : 0.0;
	  //MFG[i][j] = Cij;
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

   int m = 3 * cg_num;
   int lda = m;
   int size_a = m * m;
   int ipiv[m];

   double* A = new double[size_a];
   for (int i=0; i<m; i++)
   {
      for (int j=0; j<m; j++)
      {
	 A[i * lda + j] = MCG[i][j];
      }
   }

   LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, m, A, lda, ipiv);
   LAPACKE_dgetri(LAPACK_ROW_MAJOR, m, A, lda, ipiv);  //Now A stores the inverse of (1-M)

   //Do the matrix Multiplication of (1-M)^(-1) with (C+N)
   enum CBLAS_ORDER order;
   enum CBLAS_TRANSPOSE transa;
   enum CBLAS_TRANSPOSE transb;

   double alpha, beta;
   int n, k, ldb, ldc;

   order = CblasRowMajor;
   transa = CblasNoTrans;
   transb = CblasNoTrans;

   m = 3 * cg_num;
   n = 3 * fg_num;
   k = 3 * cg_num;

   lda = k;
   ldb = n;
   ldc = n;

   alpha = 1.0;
   beta = 0.0;

   double* B = new double[m * n];
   for (int i=0; i<m; i++)
   {
      for (int j=0; j<n; j++)
      {
	B[i * ldb + j] = MFG[i][j];
      }
   }

   cblas_dgemm(order, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, vMap, ldc);
   printf("passing vMap\n");
   
   delete[] B;
   delete[] A;

   // Computing the CG Velocities
   lda = n;
   int incx = 1;
   int incy = 1;
   cblas_dgemv(order, transa, m, n, alpha, vMap, n, v_fg, incx, beta, V_CG, incy);
   printf("passing calculating V_CG\n");
}

// Euler's 1st order integrator for ODE 
void Engine::integrate()
{
   double** R = cg_sites->R;
   double*  V_CG = cg_sites->V;
   double   dtv = cg_sites->timestep;
   double   L = fg_atoms->L;

   for (int i=0; i<cg_num; i++)   
   {
      R[i][0] += V_CG[i] * dtv;
      R[i][1] += V_CG[i + cg_num] * dtv;
      R[i][2] += V_CG[i + 2*cg_num] * dtv;

      if (R[i][0] > L) R[i][0] -= L;
      else if (R[i][0] < 0) R[i][0] += L;
      if (R[i][1] > L) R[i][1] -= L;
      else if (R[i][1] < 0) R[i][1] += L;
      if (R[i][2] > L) R[i][2] -= L;
      else if (R[i][2] < 0) R[i][2] += L;
   }
   
}

void Engine::checker()
{
  /*checkmap<<"ITEM: TIMESTEP"<<std::endl;
  checkmap<<fg_atoms->currentStep<<std::endl;
  checkmap<<"ITEM: NUMBER OF ATOMS"<<std::endl;
  checkmap<<cg_num<<std::endl;
  checkmap<<"ITEM: BOX BOUNDS pp pp pp"<<std::endl;
  for (int i=0; i<3; i++) checkmap<<"0 "<<fg_atoms->L<<std::endl;
  checkmap<<"ITEM: ATOMS id type m x y z vx vy vz"<<std::endl;*/

  double** r = fg_atoms->r;
  double** C = matrix_C->C;
  double   L = fg_atoms->L;
  double   error = 0.0;

  for (int i=0; i<cg_num; i++)
  {
      double Rix = 0.0;
      double Riy = 0.0;
      double Riz = 0.0;

      double CI_square = 0.0;
      for (int j=0; j<fg_num; j++)
      {
	 double r_shift = r[j][0];
	 if ((r_shift - cg_sites->R[i][0]) > 0.5 * L) r_shift -= L;
	 else if ((r_shift - cg_sites->R[i][0]) < -0.5 * L) r_shift += L;
	 Rix += C[i][j] * r_shift;

	 r_shift = r[j][1];
	 if ((r_shift - cg_sites->R[i][1]) > 0.5 * L) r_shift -= L;
	 else if ((r_shift - cg_sites->R[i][1]) < -0.5 * L) r_shift += L;
	 Riy += C[i][j] * r_shift;

	 r_shift = r[j][2];
	 if ((r_shift - cg_sites->R[i][2]) > 0.5 * L) r_shift -= L;
	 else if ((r_shift - cg_sites->R[i][2]) < -0.5 * L) r_shift += L;
	 Riz += C[i][j] * r_shift;

	CI_square += C[i][j] * C[i][j];
      }

      checkmap<<CI_square<<std::endl;

      error += (Rix - cg_sites->R[i][0]) * (Rix - cg_sites->R[i][0]) + (Riy - cg_sites->R[i][1]) * (Riy - cg_sites->R[i][1]) + (Riz - cg_sites->R[i][2]) * (Riz - cg_sites->R[i][2]);
      //checkmap<<i+1<<' '<<1<<' '<<Rix - cg_sites->R[i][0]<<' '<<Riy - cg_sites->R[i][1]<<' '<<Riz - cg_sites->R[i][2]<<' '<<std::endl;
  }

  //checkmap<<error<<std::endl;

  for (int ind=0; ind<3*cg_sites->cg_num; ind++)
  {
  	double sum_bij = 0.0;
  	for (int i=0; i<3*fg_atoms->fg_num; i++)
		sum_bij += vMap[ind * 3 * fg_num + i] * vMap[ind * 3 * fg_num + i];

  	invMass<<sum_bij<<std::endl;
  }

}

void Engine::compute_ke()
{

}


void Engine::endOfFrame()
{
   cg_sites->output();
   //fg_atoms->output();
   checker();
   fg_atoms->readNextFrame();
}

