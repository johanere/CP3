#ifndef SOLVER_H
#define SOLVER_H

#include "vec3.h"
#include "system.h"

class Solver
{
private:

public:
    // For now, public variables...
    System problem;
    // Default, parametrized constructor
    Solver(System newsystem);

    // Destructor
    ~Solver() {}


};

#endif // SOLVER_H
