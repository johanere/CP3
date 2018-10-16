#include "celestialbody.h"

//initialiseres med N: lager x[n],y[n].... E_k,E_p,mom [n]
//setter x[0]=x_0

CelestialBody::CelestialBody(vec3 newPosition, vec3 newVelocity, double newMass, string newName) :
    position(newPosition), velocity(newVelocity), mass(newMass), name(newName)
{
    force = {0,0,0};
}

void CelestialBody::printObject()
{
    cout << name << " Position: " << position << " Velocity: " << velocity << endl;
}
