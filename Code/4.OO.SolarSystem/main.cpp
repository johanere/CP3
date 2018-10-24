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

    // Initializes planets in solar system
    double UnitMass=1988500e24; //kg

    double veladjust=365.242199;

    double massSun = 1;
    vec3 posSun(2.192279039707197E-03, 5.762703037438353E-03, -1.295390859661947E-04);
    vec3 velSun(-5.271633313328310E-06,  5.466209510902422E-06, 1.241657817440020E-07);

    double massEarth = 5.97219e24/UnitMass;
    vec3 posEarth(8.679004239479287E-01,4.962296089500106E-01,-1.566030887927330E-04);
    vec3 velEarth(-8.771827063871312E-03, 1.491359022662084E-02, -3.259453001542365E-07);

    double massJupiter = 1.9e27/UnitMass;
    vec3 posJupiter(-4.560765261803390E+00,-2.957038080810524E+00, 1.142760185974161E-01);
    vec3 velJupiter(4.016730980673702E-03, -5.973206902447022E-03, -6.504563330645101E-05);

    double massMars = 6.6e23/UnitMass;
    vec3 posMars(-1.586900999129547E+00, 5.001555256353886E-01, 4.922989415831799E-02);
    vec3 velMars(-3.638385221590926E-03, -1.216093048374219E-02, -1.656655640290856E-04);

    double massVenus = 4.9e24/UnitMass;
    vec3 posVenus(-6.859478142184278E-01,2.103240483071906E-01, 4.238650573048407E-02);
    vec3 velVenus(-5.867791574672555E-03, -1.947450171655608E-02, 7.121929284791517E-05);

    double massSaturn = 5.5e26/UnitMass;
    vec3 posSaturn(-3.211170419783875E-01,-1.005045599567478E+01, 1.875288319015644E-01);
    vec3 velSaturn(5.270389048029983E-03, -1.958598991486847E-04, -2.065248307096866E-04);

    double massMercury = 3.3e23/UnitMass;
    vec3 posMercury(-2.327991908789227E-01,-3.904847890125498E-01, -1.095001531240173E-02);
    vec3 velMercury(1.851753087417236E-02, -1.299992586068514E-02, -2.761863418814871E-03);

    double massUranus = 8.8e25/UnitMass;
    vec3 posUranus(1.784901934033485E+01,8.829883531158330E+00, -1.984425267310287E-01);
    vec3 velUranus(-1.772851373238405E-03, 3.341974951808675E-03, 3.527445872207531E-05);

    double massNeptune = 1.03e26/UnitMass;
    vec3 posNeptune(2.861925865229266E+01,-8.803228391659726E+00, -4.782740103805204E-01);
    vec3 velNeptune(9.022088893823943E-04, 3.018763798642282E-03, -8.336640811463031E-05);

    //Adjusting velocities
    velSun = velSun*veladjust;
    velEarth = velEarth*veladjust;
    velJupiter = velJupiter*veladjust;
    velMars = velMars*veladjust;
    velVenus = velVenus*veladjust;
    velSaturn = velSaturn*veladjust;
    velMercury = velMercury*veladjust;
    velUranus = velUranus*veladjust;
    velNeptune = velNeptune*veladjust;

    //initialize objects
    CelestialBody *sun = new CelestialBody(posSun, velSun, massSun, "sun");
    CelestialBody *earth = new CelestialBody(posEarth, velEarth, massEarth, "earth");
    CelestialBody *jupiter = new CelestialBody(posJupiter, velJupiter, massJupiter, "jupiter");
    CelestialBody *mars = new CelestialBody(posMars, velMars, massMars, "mars");
    CelestialBody *venus = new CelestialBody(posVenus, velVenus, massVenus, "venus");
    CelestialBody *saturn = new CelestialBody(posSaturn, velSaturn, massSaturn, "saturn");
    CelestialBody *mercury = new CelestialBody(posMercury, velMercury, massMercury, "mercury");
    CelestialBody *uranus = new CelestialBody(posUranus, velUranus, massUranus, "uranus");
    CelestialBody *neptune = new CelestialBody(posNeptune, velNeptune, massNeptune, "neptune");

    double G = 4*M_PI*M_PI;

    System SolarSystemPlanets;  //construct system



    NewtonianGravity force(G, massSun);  //specify force

    SolarSystemPlanets.addObject(sun);
    SolarSystemPlanets.addObject(earth);
    SolarSystemPlanets.addObject(jupiter);
    SolarSystemPlanets.addObject(mars);
    SolarSystemPlanets.addObject(venus);
    SolarSystemPlanets.addObject(saturn);
    SolarSystemPlanets.addObject(mercury);
    SolarSystemPlanets.addObject(uranus);
    SolarSystemPlanets.addObject(neptune);

    Solver Sol(SolarSystemPlanets, force);

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
