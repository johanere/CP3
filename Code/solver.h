#ifndef SOLVER_H
#define SOLVER_H

#include "vec3.h"
#include "system.h"
#include "newtoniangravity.h"

class Solver
{
private:
public:
    // For now, public variables...
    System problem;
    NewtonianGravity force;
    // Default, parametrized constructor
    Solver(System newsystem,NewtonianGravity newforce);

    // Destructor
    ~Solver() {}

    void VerletSolve(int T, int N, int printstep);

};

#endif // SOLVER_H
