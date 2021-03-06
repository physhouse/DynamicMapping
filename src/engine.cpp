#include "engine.h"
#include "geom.h"
#include "pointers.h"
#include "mapping.h"
#include "matrix_C.h"
#include "matrix_M.h"
#include "matrix_N.h"
#include "cg_sites.h"
#include "fg_atoms.h"
#include "neighbor.h"
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include "lapacke.h"
#include "cblas.h"
#endif
#include <ctime>

#ifdef __APPLE__
void check_MatrixSolver_lapack(const char * const name, const int info) {
    if (info != 0) {
        if (info < 0) {
            printf("Illegal call to %s in MatrixSolver, argument %d", name, -info);
        } else {
            printf("Singular matrix encountered in %s in MatrixSolver, element (%d, %d)", name, info, info);
        }
        exit(EXIT_FAILURE);
    }
}
#else
extern void cblas_dgemv(const CBLAS_ORDER layout, const CBLAS_TRANSPOSE TransA,
                        const int M, const int N, const double alpha, const double *A, const int lda,
                        const double *X, const int incX, const double beta, double *Y, const int incY);
#endif

Engine::Engine(Mapping *map) : Pointers(map) {}
Engine::~Engine()
{
    cleanup();
}

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

    IMinusM = new double* [size_cg];
    CPlusN = new double* [size_cg];

    vMap = new double[size_cg * size_fg];
    memset(vMap, 0, sizeof(double) * size_cg * size_fg);

    for (int i = 0; i < size_cg; i++)
    {
        IMinusM[i] = new double[size_cg];
        memset(IMinusM[i], 0, sizeof(double) * size_cg);
        CPlusN[i] = new double[size_fg];
        memset(CPlusN[i], 0, sizeof(double) * size_fg);
    }

    checkmap.open("CG_check.lmpstrj", std::ofstream::out);
    checkmap << "CG_CHECK : Checking the Propagation" << std::endl;

    invMass.open("Mat_miu.dat", std::ofstream::out);
}

void Engine::cleanup()
{
    checkmap.close();
    invMass.close();

    int size_cg = 3 * cg_num;

    for (int i = 0; i < size_cg; i++)
    {
        delete[] IMinusM[i];
        delete[] CPlusN[i];
    }
    delete[] IMinusM;
    delete[] CPlusN;
    delete[] vMap;
}

void Engine::buildNeighbors()
{
    std::clock_t start;
    double cellTime, allTime;

    start = std::clock();
    neighbor->buildCellList();
    cellTime = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    neighbor->buildNeighborList();
    allTime = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    printf("cell: %12.8lf  neighbor: %12.8lf\n ncells %d\n", cellTime, allTime - cellTime, neighbor->dimCell);
}

//Update Velocity for one frame
void Engine::update()
{
    computeMatrices();
    matrixSolver();
    integrate();
    cg_sites->output();
    checker();
}

void Engine::exec()
{
    bool do_testing = false;
    printf("nsteps = %d\n", nsteps);

    if (!restartFlag)
    {
        cg_sites->IterateRMapToSelfConsistency();
        cg_sites->output();
        checker();
        fg_atoms->readNextFrame();
    }

    //cg_sites->output();
    for (int i = 0; i < nsteps - 1; i++)
    {
        update();
        fg_atoms->readNextFrame();
        if (do_testing) testing(i);

        // Rebuild the neighbor list every timestep after
        // reading the new FG frame.
        if (i % 1 == 0) buildNeighbors();
        // Enforce self-consistency explicitly at periodic intervals
        // to ensure stable integration on the constraint manifold.
        // This belongs in the integrator.
        if (fg_atoms->currentStep % cg_sites->freq == 0)
        {
            printf("Step %d\n", fg_atoms->currentStep);
            if (do_testing) testing(i);
            cg_sites->IterateRMapToSelfConsistency();
        }
    }
    // The last frame, no need to load next frame.
    update();
    fg_atoms->finishReading();
}

void Engine::testing(int step)
{
    printf("Testing M... Current Step %d\n", step);
    for (int i = 0; i < cg_sites->cg_num; i++)
    {
        for (int j = 0; j < cg_sites->cg_num; j++)
        {
            printf("%8.8lf\t", matrix_M->M[0][0][i][j]);
        }
        printf("\n");
    }

    printf("Testing N... Current Step %d\n", step);
    for (int i = 0; i < cg_sites->cg_num; i++)
    {
        for (int j = 0; j < fg_atoms->fg_num; j++)
        {
            printf("%8.8lf\t", matrix_N->N[0][0][i][j]);
        }
        printf("\n");
    }

    printf("Testing C... Current Step %d\n", step);
    for (int i = 0; i < cg_sites->cg_num; i++)
    {
        printf("%d : sum = %12.8lf\n", i + 1, matrix_C->sumC[i]);
        for (int j = 0; j < fg_atoms->fg_num; j++)
        {
            printf("%12.8lf\t", matrix_C->C[i][j]);
        }
        printf("\n");
    }

    printf("Testing dW...\n");
    for (int i = 0; i < cg_sites->cg_num; i++)
    {
        printf("%d : \n", i + 1);
        for (int j = 0; j < fg_atoms->fg_num; j++)
        {
            printf("%e\t", matrix_C->dw[i][j][0]);
        }
        printf("\n");
        for (int j = 0; j < fg_atoms->fg_num; j++)
        {
            printf("%e\t", matrix_C->W[i][j]);
        }
        printf("\n");
    }

    printf("Tesing Neighbor...\n");
    int fg2cg = 0;
    int cg2fg = 0;
    for (int i = 0; i < cg_sites->cg_num; i++)
    {
        cg2fg += neighbor->numNeighbors[i];
    }
    for (int i = 0; i < fg_atoms->fg_num; i++)
    {
        fg2cg += neighbor->numFgNeighbors[i];
    }
    printf("Neighbors: %d, FGNeighbors: %d\n", cg2fg, fg2cg);

    printf("Testing WSUM...\n");
    for (int i = 0; i < fg_atoms->fg_num; i++)
    {
        printf("%d : %e\t", i, matrix_C->w_sum[i]);
    }
    printf("\n");
}

void Engine::computeMatrices()
{
    matrix_C->compute();
    printf("finishing C\n");
    matrix_M->compute();
    printf("finishing M\n");
    matrix_N->compute();
    printf("finishing N\n");

    //Filling in the matrices and vectors
    double **    **M = matrix_M->M;
    double **    **N = matrix_N->N;
    double      **C = matrix_C->C;
    //Building the CPlusN and IMinusM matrices
    for (int idim = 0; idim < 3; idim++)
    {
        int offset_d1 = idim * cg_num;
        for (int i = offset_d1; i < (offset_d1 + cg_num); i++)
        {
            //IMinusM = 1 - M
            for (int jdim = 0; jdim < 3; jdim++)
            {
                int offset_d2 = jdim * cg_num;
                for (int j = offset_d2; j < (offset_d2 + cg_num); j++)
                {
                    double deltaij = (i == j) ? 1.0 : 0.0;
                    //IMinusM[i][j] = deltaij;
                    IMinusM[i][j] = deltaij - M[idim][jdim][i - offset_d1][j - offset_d2];
                }
            }
            //CPlusN = C + N
            for (int jdim = 0; jdim < 3; jdim++)
            {
                int offset_d2 = jdim * fg_num;
                for (int j = offset_d2; j < (offset_d2 + fg_num); j++)
                {
                    double Cij = (idim == jdim) ? C[i - offset_d1][j - offset_d2] : 0.0;
                    //CPlusN[i][j] = Cij;
                    CPlusN[i][j] = Cij + N[idim][jdim][i - offset_d1][j - offset_d2];
                }
            }
        }
    }

}

void Engine::matrixSolver()
{
    double *V_CG = cg_sites->V;
    double *v_fg = fg_atoms->v;
    // Calculate CPlusN*v_fg

    int m = 3 * cg_num;
    int lda = m;
    int size_a = m * m;
    int ipiv[m];

    double *A = new double[size_a];
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            A[i * lda + j] = IMinusM[i][j];
        }
    }

#ifdef __APPLE__
    int info;
    int lwork = m;
    double work[m];
    dgetrf_(&m, &m, A, &lda, ipiv, &info);
    check_MatrixSolver_lapack("dgetrf", info);
    dgetri_(&m, A, &lda, ipiv, work, &lwork, &info);
    check_MatrixSolver_lapack("dgetri", info);
#else
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, m, A, lda, ipiv);
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, m, A, lda, ipiv);  //Now A stores the inverse of (1-M)
#endif

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

    double *B = new double[m * n];
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            B[i * ldb + j] = CPlusN[i][j];
        }
    }


    // Check CNV vector
    bool check_cnv_vector = true;
    if (check_cnv_vector) {
        cblas_dgemv(order, transa, m, n, alpha, B, n, v_fg, 1, beta, V_CG, 1);
        checkmap << "Printing Nv. Step: " << fg_atoms->currentStep << ", N_CG: "<< cg_num << std::endl;
        for (int i = 0; i < cg_num; i++)
        {
            checkmap << V_CG[i] << '\t' << V_CG[i + cg_num] << '\t' << V_CG[i + 2 * cg_num] << std::endl;
            //checkmap << V_CG[i] - cg_sites->VMAP[i] << '\t' << V_CG[i + cg_num] - cg_sites->VMAP[i + cg_num] << '\t' << V_CG[i + 2 * cg_num] - cg_sites->VMAP[i + 2 * cg_num] << std::endl;
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
    double **R = cg_sites->R;
    double  *V_CG = cg_sites->V;
    double   dtv = cg_sites->timestep;
    double   L = fg_atoms->L;

    for (int i = 0; i < cg_num; i++)
    {
        R[i][0] += V_CG[i] * dtv;
        R[i][1] += V_CG[i + cg_num] * dtv;
        R[i][2] += V_CG[i + 2 * cg_num] * dtv;

        for (int j = 0; j < 3; j++) wrap_coord_into_box(R[i][j], L);
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

    double **C = matrix_C->C;
    double   error = 0.0;

    bool check_mapped_positions = false;
    bool check_sum_c_squared = false;
    for (int i = 0; i < cg_num; i++)
    {
        double recalc_R[3];
        for (int dim = 0; dim < 3; dim++) recalc_R[dim] = 0;

        double CI_square = 0.0;

        // Check the difference between the current and mapped CG site positions.
        // The recalculation of the CG site positions from the map
        // should be in cg_sites.
        if (check_mapped_positions) {
            cg_sites->map_CG_position(i, recalc_R);
            error += (recalc_R[0] - cg_sites->R[i][0]) * (recalc_R[0] - cg_sites->R[i][0]) + (recalc_R[1] - cg_sites->R[i][1]) * (recalc_R[1] - cg_sites->R[i][1]) + (recalc_R[2] - cg_sites->R[i][2]) * (recalc_R[2] - cg_sites->R[i][2]);
            checkmap << i + 1 << ' ' << 1 << ' ' << recalc_R[0] - cg_sites->R[i][0] << ' ' << recalc_R[1] - cg_sites->R[i][1] << ' ' << recalc_R[2] - cg_sites->R[i][2] << ' ' << std::endl;
        }

        // Check the sum of squares of c_Ii coefficients.
        if (check_sum_c_squared) {
            for (int j = 0; j < fg_num; j++) {
                CI_square += C[i][j] * C[i][j];
            }
            checkmap << CI_square << std::endl;
        }
    }

    if (check_mapped_positions) {
        checkmap << error << std::endl;
    }

    // Check the current inverse mass of the CG site.
    for (int i = 0; i < 3 * cg_sites->cg_num; i++)
    {
        double sum_bij_sq = 0.0;
        for (int j = 0; j < 3 * fg_atoms->fg_num; j++)
            sum_bij_sq += vMap[i * 3 * fg_num + j] * vMap[i * 3 * fg_num + j];

        invMass << sum_bij_sq << std::endl;
    }

}
