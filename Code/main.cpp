#include <iostream>
#include <cmath>
#include <cstdlib>
#include "vec3.h"
#include "celestialbody.h"
#include "system.h"
#include "forwardeuler.h"
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
    vec3 posEarth(9.079361611992910E-01,4.174196284373897E-01,-9.320364656857048E-05);
    vec3 velEarth(-7.374461784272602E-03,1.561335766240027E-02,1.196212616732829E-07); // v=2pi r/period -> r=1 AU, period=1 yr (standard units)


    double massMars = 6.6e23/UnitMass; // massEarth / UnitMass
    vec3 posMars(1.384662720209662E+00,-8.711968326251433E-02,-3.603051751471911E-02);
    vec3 velMars(1.478464132931038E-03,1.515952151725583E-02,2.813270986502590E-04); // v=2pi r/period -> r=1 AU, period=1 yr (standard units)

    velMars=velMars*veladjust;
    velEarth=velEarth*veladjust;

    //initialize objects
    CelestialBody *sun = new CelestialBody(posSun, velSun, massSun, "sun");
    CelestialBody *earth = new CelestialBody(posEarth, velEarth, massEarth, "earth");
    CelestialBody *mars = new CelestialBody(posMars, velMars, massMars, "mars");
    double G = 4*M_PI*M_PI;

    System SunEarth;  //construct system


    NewtonianGravity force(G, massSun);  //specify force

    //SunEarth.addObject(sun);
    SunEarth.addObject(earth);
    //SunEarth.addObject(mars);


    Solver Sol(SunEarth,force);


    Sol.VerletSolve(T, N, printstep);

    return 0;
}
