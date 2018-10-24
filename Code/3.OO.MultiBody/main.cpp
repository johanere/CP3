#include <iostream>
#include <cmath>
#include <cstdlib>
#include "vec3.h"
#include "celestialbody.h"
#include "system.h"
#include "newtoniangravity.h"
#include "solver.h"

using std::cout;
using std::endl;

int main(int argc, char * argv[])
{

    // Specify discretization paramters through argsys
    int N = atoi(argv[1]);
    int printstep = atoi(argv[2]);
    int T = atoi(argv[3]);
    
    // Initializes sun, earth
    double UnitMass=1988500e24; //kg

    double veladjust=365.242199;

    double massSun = 1; //
    vec3 posSun(0,0,0);
    vec3 velSun(0,0,0);

    double massEarth = 5.97219e24/UnitMass; // massEarth / UnitMass
    vec3 posEarth(8.679004239479287E-01,4.962296089500106E-01,-1.566030887927330E-04);
    vec3 velEarth(-8.771827063871312E-03, 1.491359022662084E-02, -3.259453001542365E-07); // v=2pi r/period -> r=1 AU, period=1 yr (standard units)


    double massJupiter = 1.9e27/UnitMass; // massEarth / UnitMass
    vec3 posJupiter(-4.560765261803390E+00,-2.957038080810524E+00, 1.142760185974161E-01);
    vec3 velJupiter(4.016730980673702E-03, -5.973206902447022E-03, -6.504563330645101E-05); // v=2pi r/period -> r=1 AU, period=1 yr (standard units)

    velEarth=velEarth*veladjust;
    velJupiter=velJupiter*veladjust;

    //initialize objects
    CelestialBody *sun = new CelestialBody(posSun, velSun, massSun, "sun");
    CelestialBody *earth = new CelestialBody(posEarth, velEarth, massEarth, "earth");
    CelestialBody *jupiter = new CelestialBody(posJupiter, velJupiter, massJupiter, "jupiter");

    double G = 4*M_PI*M_PI;

    System SunEarthJupiter;  //construct system


    NewtonianGravity force(G, massSun);  //specify force

    SunEarthJupiter.addObject(sun);
    SunEarthJupiter.addObject(earth);
    SunEarthJupiter.addObject(jupiter);

    //SunEarth.addObject(mars);


    Solver Sol(SunEarthJupiter, force);

    Sol.VerletSolve(T, N, printstep);



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
