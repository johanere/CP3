#include <iostream>
#include <cmath>
#include "vec3.h"
#include "celestialbody.h"
#include "system.h"
#include "forwardeuler.h"
#include "newtoniangravity.h"
#include "solver.h"

using std::cout;
using std::endl;

int main()
{

    // Specify discretization paramters through argsys
    //hent fra argsys <- FIKS
    int N = 100;


    // Initializes sun, earth
    double UnitMass=1988500e24; //kg

    double massSun = 1; //
    vec3 posSun(0,0,0);
    vec3 velSun(0,0,0);

    double massEarth = 5.97219e24/UnitMass; // massEarth / UnitMass
    vec3 posEarth(1,0,0);
    vec3 velEarth(0,2*M_PI,0); // v=2pi r/period -> r=1 AU, period=1 yr (standard units)

    //initialize objects
    CelestialBody *sun = new CelestialBody(posSun, velSun, massSun, "sun");
    CelestialBody *earth = new CelestialBody(posEarth, velEarth, massEarth, "earth");

    double G = 4*M_PI*M_PI;

    System SunEarth;  //construct system


    NewtonianGravity force(G, massSun);  //specify force

    //add objects to system class
    SunEarth.addObject(sun);
    SunEarth.addObject(earth);



    Solver Sol(SunEarth,force);


    Sol.EulerSolve(N);




/*
    Solver sol(SunEarth); //construct solver (should be *sol= new Solver(SunEarth,N)? )
    sol.EulerSolve(N); //solves and prints solution to outfile





    cout << "Number of bodies in the solar system is: " << S.bodies.size() << endl;

    // Prints start position
    S.bodies[1]->printObject();

    for (int iStep = 0; iStep < NSteps; iStep++)
    {
        S.resetForces(); // Since we only stores one step at the time...
        force.calculateForces(&S);
        integrator.integrate(&S, h);

        for (int iObj = 0; iObj < S.bodies.size(); iObj++)
        {
            S.bodies[iObj]->printObject();
        }

    }

    // Prints end position, which should be approximately the same as the start
    S.bodies[1]->printObject();
*/

    return 0;
}
