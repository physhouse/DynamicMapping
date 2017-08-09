#include "cg_sites.h"
#include "mapping.h"
#include "fg_atoms.h"
#include "matrix_C.h"
#include "pointers.h"
#include "comm.h"
#include "assert.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>

Cg_sites::Cg_sites(Mapping *map) : Pointers(map) {}
Cg_sites::~Cg_sites()
{
    cleanup();
}

void Cg_sites::init(int argc, char **argv)
{
    freq = 50;

    VAR_BEGIN
    GET_INT(cg_num)
    GET_REAL(tol)
    GET_REAL(timestep)
    GET_REAL(sigmaOverRcut)
    GET_INT(niter)
    GET_INT(freq)
    VAR_END

    R = new double*[cg_num];
    for (int i = 0; i < cg_num; i++)
    {
        R[i] = new double[3];
        for (int j = 0; j < 3; j++) R[i][j] = 0.0;
    }

    V = new double[3 * cg_num];
    VMAP = new double[3 * cg_num];
    memset(V, 0, sizeof(double) * 3 * cg_num);
    memset(VMAP, 0, sizeof(double) * 3 * cg_num);

    M = new double[cg_num];
    memset(M, 0, sizeof(double) * cg_num);

    double M0 = (double)fg_atoms->fg_num / (double)cg_num; // hard coded! Need to be removed!
    for (int i = 0; i < cg_num; i++)
        M[i] = M0;

    L = fg_atoms->L;
    int nperdim = (int)cbrt((double)cg_num);
    rcut = 0.50 * L / nperdim;

    printf("fgnum = %d cgnum = %d L = %12.8lf rcut = %12.8lf\n", fg_atoms->fg_num, cg_num, L, rcut);
    //Grid Initialization of CG particle positions
    if (argc <= 3)
    {
        double size_block = L / nperdim;
        for (int i = 0; i < nperdim; i++)
            for (int j = 0; j < nperdim; j++)
                for (int k = 0; k < nperdim; k++)
                {
                    int ind = nperdim * nperdim * i + nperdim * j + k;
                    R[ind][0] = (i + 0.5) * size_block;
                    R[ind][1] = (j + 0.5) * size_block;
                    R[ind][2] = (k + 0.5) * size_block;
                }
    }
    else  //read Data from restart files
    {
        assert(argv[4] != NULL);
        int aid, type;
        FILE *fp = fopen(argv[4], "r+");
        for (int i = 0; i < cg_num; i++)
        {
            fscanf(fp, "%d %d %lf %lf %lf %lf", &aid, &type, &M[i], &R[i][0], &R[i][1], &R[i][2]);
        }
        fclose(fp);
        printf("Restart Parsing Finished\n");
    }

    cgtrj.open("CG_TRJ.lmpstrj", std::ofstream::out);
    cgtrj << "CG_TRJ : Dynamical Mapping Coarse Graining" << std::endl;

    // Print out the initial configuration if it needs to be checked.
    if (false) {
        cgtrj<<"ITEM: TIMESTEP"<<std::endl;
        cgtrj<<0<<std::endl;
        cgtrj<<"ITEM: NUMBER OF ATOMS"<<std::endl;
        cgtrj<<cg_num<<std::endl;
        cgtrj<<"ITEM: BOX BOUNDS pp pp pp"<<std::endl;
        for (int i=0; i<3; i++) cgtrj<<"0 "<<fg_atoms->L<<std::endl;
        cgtrj<<"ITEM: ATOMS id type m x y z vx vy vz"<<std::endl;

        for (int i=0; i<cg_num; i++)
        {
          cgtrj<<i+1<<' '<<1<<' '<<M[i]<<' '<<R[i][0]<<' '<<R[i][1]<<' '<<R[i][2]<<' '<<V[i]<<' '<<V[i + cg_num]<<' '<<V[i + 2*cg_num]<<std::endl;
        }
    }
}

void Cg_sites::cleanup()
{
    cgtrj.close();
    for (int i = 0; i < cg_num; i++)
    {
        delete[] R[i];
    }
    delete[] R;
    delete[] V;
    delete[] M;
}

void Cg_sites::output()
{
    double *v = fg_atoms->v;
    double **C = matrix_C->C;

    for (int i = 0; i < cg_num; i++)
    {
        VMAP[i] = 0.0;
        VMAP[i + cg_num] = 0.0;
        VMAP[i + 2 * cg_num] = 0.0;
        int fg = fg_atoms->fg_num;
        for (int j = 0; j < fg; j++)
        {
            VMAP[i]        += C[i][j] * v[j];
            VMAP[i + cg_num]   += C[i][j] * v[j + fg];
            VMAP[i + 2 * cg_num] += C[i][j] * v[j + 2 * fg];
        }
    }
    //header
    cgtrj << "ITEM: TIMESTEP" << std::endl;
    cgtrj << fg_atoms->currentStep << std::endl;
    cgtrj << "ITEM: NUMBER OF ATOMS" << std::endl;
    cgtrj << cg_num << std::endl;
    cgtrj << "ITEM: BOX BOUNDS pp pp pp" << std::endl;
    for (int i = 0; i < 3; i++) cgtrj << "0 " << fg_atoms->L << std::endl;
    cgtrj << "ITEM: ATOMS id type m x y z vx vy vz fx fy fz" << std::endl;

    for (int i = 0; i < cg_num; i++)
    {
        cgtrj << i + 1 << ' ' << 1 << ' ' << M[i] << ' ' << R[i][0] << ' ' << R[i][1] << ' ' << R[i][2] << ' ' << V[i] << ' ' << V[i + cg_num] << ' ' << V[i + 2 * cg_num] << ' ' << VMAP[i] << ' ' << VMAP[i + cg_num] << ' ' << VMAP[i + 2 * cg_num] << std::endl;
    }
}

// Iterate the CG position map to self-consistency.

void Cg_sites::IterateRMapToSelfConsistency()
{
    double sum_of_sq_displacements = 0.0;
    double alpha = 0.70;

    for (int round = 0; round < niter; round++)
    {
        // For each CG site, update the position according to
        // the most recent coefficients, keeping track of how
        // much it was moved and the sum of squares of those
        // displacement vectors.
        sum_of_sq_displacements = 0.0;
        for (int i = 0; i < cg_num; i++)
        {
            double R_initial[3];
            double R_recalc[3];

            for (int dim = 0; dim < 3; dim++) {
                R_initial[dim] = R[i][dim];
                R_recalc[dim] = 0;
            }
            matrix_C->recalc_CG_position(i, R_recalc);

            // The next fixed point iterate is a fraction of the original
            // plus a fraction of the new. 
            sum_of_sq_displacements = 0;
            for (int dim = 0; dim < 3; dim++) {
                R[i][dim] = alpha * R_initial[dim] + (1 - alpha) * R_recalc[dim];
                sum_of_sq_displacements += (R_initial[dim] - R_recalc[dim]) * (R_initial[dim] - R_recalc[dim]);
            }
        }

        if (sum_of_sq_displacements < tol) return;
        else
            matrix_C->compute();
        if (round % 2 == 0)
        {
            printf("Round %d, error = %e\n", round, sum_of_sq_displacements);
        }
    }
    printf("Iterative search for self-consistent CG map failed to converge with %d iterations.", niter);
    exit(EXIT_FAILURE);
}

