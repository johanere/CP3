#include "solver.h"
#include "system.h"
#include "vec3.h"
//construct TBS TwoBodySolver(system,N)
//har funksjonene EulerSolve og VerletSolve
//printsolution

Solver::Solver(System newsystem,NewtonianGravity newforce) : problem(newsystem), force(newforce)
{
}
//legg til RESET til initial values funksjon!
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
*/

void Solver::EulerSolve(int N)
    {
      //initialize variables
      vec3 NewPos, NewVel, f,x,v;
      vec3 f_old[problem.bodies.size()];
      vec3 x_new[problem.bodies.size()];
      double h,NewE_K,NewE_P,m,h2,two_m;
      h=(double) 1/(N-1);
      cout<<"h"<<h<<endl;
      h2=h*h;
      CelestialBody* obj;

      problem.bodies[1]->printObject();

      //forces at initial positions, f[0]
      force.calculateForces(&problem);
      //calculate x[1] using x[0] v[0] f[0]
      for (int i=0;i<problem.bodies.size();i++)
        {

          obj=problem.bodies[i];
          f=obj->force; //unload force
          two_m=obj->mass; //unload mass
          x=obj->position;//unload position
          v=obj->velocity;//unload velocity
          x_new[i] = x + v*h + f*h2/(two_m);
          f_old[i] = f; //Store current force for body i
        }
      //have x[1],f[0]=f_old
      //loop:
        //find f[i]
        //find v[i] using v[i-1],f_old,f[i]
        //find x[i+1] using x[i],v[i],f[i]
        //set f_oldf[i]
      for (int j=0;j<N-2;j++)
      {
        for (int i=0;i<problem.bodies.size();i++)
          {
            obj->position=x_new[i];
          //obj->force =vec3(0,0,0);
          }
        problem.bodies[1]->printObject();
        problem.resetForces();
        force.calculateForces(&problem);
        for (int i=0;i<problem.bodies.size();i++)
          {
            obj=problem.bodies[i];
            f=obj->force; //unload force
            two_m=obj->mass; //unload mass
            v=obj->velocity;//unload velocity
            obj->velocity= v  + (f_old[i] + f)*h/(two_m);   //calculate new velocity
            v=obj->velocity;//unload velocity
            x_new[i]= x + v*h + f*h2/(two_m);
            f_old[i] = f; //Store current force for body i
          }

          problem.bodies[1]->printObject();
          cout<<endl;

        }




      //find final velocity and position


      force.calculateForces(&problem);
      for (int i=0;i<problem.bodies.size();i++)
        {
          obj=problem.bodies[i];
          f=obj->force; //unload force
          two_m=obj->mass; //unload mass
          x=obj->position;//unload position
          v=obj->velocity;//unload velocity
          obj->velocity=v+(f_old[i]+f)*h/(two_m);
          obj->position=x_new[i];
        }
      }

/*
      //define outfile
      //Solve.initializeprint("Euler", string );

        force.calculateForces(&problem);
        for (int i=0;i<problem.bodies.size();i++)
          {

            obj=problem.bodies[i];
            f=obj->force; //unload force
            m=obj->mass; //unload mass
            x=obj->position;//unload position
            v=obj->velocity;//unload velocity
            f_old[i] = f; //Store current force for body i
            //calculate new values
            NewPos=x+v*h+f*h2/(2*m);

          }

        force.calculateForces(&problem);

        }

      }

          //update values
          //  obj->pos=NewPos
          //  obj->vel=NewVel
            //obj->kineticEnergy=NewE_K
            //obj->potentialEnergy=NewE_P
            //write to outfile
          //  write to outfile: obj->pos,vel,kineticEnergy,potentialEnergy //each line after line 2 contains 4 paramters per planet.

          //outfile move to next line



    // It is the same as following:
//    for (int iObj = 0; iObj < bodies.size(); iObj++) {
//        bodies[iObj]->force = {0,0,0};
//    }
    //void resetForces();*/
