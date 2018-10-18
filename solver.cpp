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

      int Number_of_Bodies = problem.bodies.size();

      vec3 f_i[Number_of_Bodies];
      vec3 f_iplus1[Number_of_Bodies];



      double h, hh;

      //Step-Size
      h=(double) 1/(N-1);
      hh=h*h;

      //Defining celestial body pointer object
      CelestialBody* obj;

      problem.bodies[1]->printObject();

      //forces at initial positions, f[0]
      force.calculateForces(&problem);

      //Populating initial force
      for(int i = 0; i < Number_of_Bodies; i++)
        {
          f_i[i] = problem.bodies[i]->force;
          cout << "Initial forces" << f_i[i] << endl;
        }



      //ALT SKJER FOR TIDEN
      for(int Time = 0; Time < N; Time++)
        {

          for(int i = 0; i < Number_of_Bodies; i++)
            {
              obj = problem.bodies[i];  //Objektet earth, sun fra celestial body

              obj->position = obj->position +  obj->velocity*h +  f_i[i]*(hh/(2*obj->mass));
              cout << "New position of " << problem.bodies[1]->position << " equals " << obj->position << endl;

            }
          problem.resetForces();
          force.calculateForces(&problem);
          for(int i = 0; i < Number_of_Bodies; i++) f_iplus1[i] = problem.bodies[i]->force;  //Ny kraft


        for(int i = 0; i < Number_of_Bodies; i++)
          {

            obj = problem.bodies[i];  //Objektet earth, sun fra celestial body
            obj->velocity = obj->velocity + ( f_iplus1[i] + f_i[i] ) * (h/(2*obj->mass));
            f_i[i] = f_iplus1[i];
          }
        }

    }
/*
        }
      for (int i=0;i<problem.bodies.size();i++)
        {

          obj = problem.bodies[i];  //lagre body nr i inne i obj
          f = obj->force; //unload force
          two_m = obj->mass; //unload mass
          x = obj->position;//unload position
          v = obj->velocity;//unload velocity

          x_new[i] = x + v*h + f*hh/(two_m);
          f_old[i] = f; //Store current force for body i
        }

      }


*/
