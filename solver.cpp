#include "solver.h"
#include "system.h"
//construct TBS TwoBodySolver(system,N)
//har funksjonene EulerSolve og VerletSolve
//printsolution

Solver::Solver(System newsystem)
{
  problem=newsystem;
}
/*
//FIKS H!
//anatar at vi har f og h
    void initprint(string method,int N); (!)
      ofile.open(fileout);
      string fileout=method
        ofile<<setiosflags(ios::showpoint | ios::uppercase);

        outfile line 1: write N gridspount
        outfile line 2: write numer of planets

    void printstep();
    void printend();



    void EulerSolve(int N);
    {
      vec3 NewPos, NewVel;
      double h,NewE_K,NewE_P;
      h=1/(N-1);
      for (int j=0;j<N;j++) //remove j as something to be used?
        {

      //define outfile
      Solve.initializeprint("Euler", string );

      {
        force.calculateForces((!)*system);
        for (CelestialBody *obj : bodies)
          {
            //unload force
            f=obj->force
            m=obj->mass
            //calculate new values

            NewPos=obj->pos+h*(obj->vel)
            NewVel=obj->vel+h*f/m
            NewE_K=obj->1/2*m*len(vel)^2
            NewE_P=obj->
            //update values
            obj->pos=NewPos
            obj->vel=NewVel
            obj->kineticEnergy=NewE_K
            obj->potentialEnergy=NewE_P
            //write to outfile
            write to outfile: obj->pos,vel,kineticEnergy,potentialEnergy //each line after line 2 contains 4 paramters per planet.
          }
          outfile move to next line
      }
    }

    // It is the same as following:
//    for (int iObj = 0; iObj < bodies.size(); iObj++) {
//        bodies[iObj]->force = {0,0,0};
//    }
    void resetForces();*/
