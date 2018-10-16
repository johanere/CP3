#include <iostream>
#include <cmath>
#include "vec3.h"
#include "celestialbody.h"
#include "system.h"
#include "forwardeuler.h"
#include "newtoniangravity.h"

using std::cout;
using std::endl;


class Tclass
{
public:
    // Vector of pointers i.e. memory adresses of solar system objects of type CelestialBody
    double a;
    vec3 kaker;

    // Default constructor
    Tclass(double r, vec3 newkake);

    // Destructor
    ~Tclass() {} // By writing {} after destructor we don't have to write it explicitly in the .cpp file
  };


Tclass::Tclass(double r,vec3 newkake)
{
  a=r;
  kaker=newkake;
}


int main()
{
  vec3 b(0,0,0);
  vec3 c(3,3,3);
  double p=3;
  Tclass ola(p,b);
  ola.kaker=c;

  int N = 10;

  cout<<"huura"<<endl;

  // Initializes sun, earth
  double UnitMass=1988500e24; //kg

  double massSun = 1; //
  vec3 posSun(0,0,0);
  vec3 velSun(0,0,0);

  //initialize objects
  CelestialBody *sun = new CelestialBody(posSun, velSun, massSun, "sun");


}
