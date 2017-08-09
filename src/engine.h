#ifndef ENGINE_H
#define ENGINE_H

#include "pointers.h"
#include <fstream>

class Engine : protected Pointers
{
public:
    Engine(class Mapping *);
    ~Engine();
    void init(int argc);

    void update();
    void exec();

    void computeMatrices();
    void matrixSolver();
    void integrate();
    void buildNeighbors();
    void cleanup();
    void checker();
    void testing(int step);

private:
    double **IMinusM;
    double  *flat_CNv;
    double   ke; //Kinetic Energy

    int     cg_num;
    int     fg_num;
    int     nsteps;
    bool    restartFlag;

    std::ofstream checkmap;
    std::ofstream invMass;
};

#endif
