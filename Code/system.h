#ifndef SYSTEM_H
#define SYSTEM_H

#include <vector>
#include "celestialbody.h"

using std::vector;

//holder Objects,N

class System
{
public:
    // Vector of pointers i.e. memory adresses of solar system objects of type CelestialBody
    std::vector<CelestialBody*> bodies;

    //NewtonianGravity force;
    // Default constructor
    System();

    // Destructor
    ~System() {} // By writing {} after destructor we don't have to write it explicitly in the .cpp file

    void addObject(CelestialBody* newBody) { bodies.push_back(newBody); }
    //void defineForce(Force* newForce) { force=(newForce); }

    void resetForces();
};

#endif // SYSTEM_H
