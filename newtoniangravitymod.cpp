void GravitationalForce::calculateForcesEinsteinian(System* system)
  {
    /*
    Takes a system object and calculates the einsteinian gravitational force
    between the sun and mercury.
    */

    double c = 63239.7263; //Speed of light [AU/yr]

    double x, y, z, vx, vy, vz, mass;

    /*Unloading positions, velocities and mass*/
    x = system->bodies[1]->position[0];
    y = system->bodies[1]->position[1];
    z = system->bodies[1]->position[2];

    vx = system->bodies[1]->position[0];
    vy = system->bodies[1]->position[1];
    vz = system->bodies[1]->position[2];

    mass = system->bodies[1]->mass;

    double rr, r_cubed; //Distance between sun and mercury and same distance squared
    vec3 r_temp = system->bodies[0]->position - system->bodies[1]->position;
    rr = r_temp.lengthSquared();
    r_cubed = rr*(r_temp.length());

    vec3 rxv(0, 0, 0); //Angular momentum per unit mass

    rxv[0] = y*vz - z*vy;
    rxv[1] = z*vx - x*vz;
    rxv[2] = x*vy - y*vx;

    double ll; // angular momentum per unit mass squared

    ll = rxv.lengthSquared();

    /*Updating force on Merury*/
    system->bodies[1]->force = r_temp * (((m_G * mass)/r_cubed) * (1.0 + (3.0*ll)/(rr*c*c) ));

  }
