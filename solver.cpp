#include "solver.h"
#include "system.h"
#include "vec3.h"

#include <string>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

using namespace std;

ofstream ofile;
//constructor
Solver::Solver(System newsystem,NewtonianGravity newforce) : problem(newsystem), force(newforce)
{
}

//-------------functions-----------
void printstart(int N, int Number_of_Bodies,int printstep,int T) // run once when solving, prints paramters for solving
{
std::string number_of_steps_per_yr = std::to_string(N);
std::string T_tot = std::to_string(T);
std::string number_of_planets = std::to_string(Number_of_Bodies);
string filename;
filename = number_of_planets+"planets_"+number_of_steps_per_yr+"stepspryr_"+T_tot+"yrs.txt";
ofile.open(filename);

ofile<<setw(15)<<setprecision(8)<<N<<" "<< printstep<<" "<<Number_of_Bodies<<" "<< T<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<endl;
}

void printcurrent(double t,vec3 pos,vec3 vel, double m) //run at a given number of steps, prints current values
{
ofile<<setiosflags(ios::showpoint | ios::uppercase);
//finding angular momentum
vec3 L(pos[1]*vel[2]-pos[2]*vel[1],pos[2]*vel[0]-pos[0]*vel[2],pos[0]*vel[1]-pos[1]*vel[0]);
L=L*m;
//writing in file
ofile<<setw(15)<<setprecision(8)<<t<<" "<<pos[0]<<" "<< pos[1]<<" "<<pos[2]<<" "<<vel[0]<<" "<< vel[1]<<" "<<vel[2]<<" "<<L[0]<<" "<< L[1]<<" "<<L[2]<<endl;
}

void printstop() //run once when solving, closing the outfile
{
ofile.close();
}

//solver function
void Solver::EulerSolve(int T, int N, int printstep)
    {

      int Number_of_Bodies = problem.bodies.size();

      //initialize variables
      vec3 f_i[Number_of_Bodies];
      vec3 f_iplus1[Number_of_Bodies];
      double h, hh;
      CelestialBody* obj;

      //Step-Size
      h=(double) 1/(N-1);
      hh=h*h;


      //start print
      if(printstep<=N)
      {
      printstart(N, Number_of_Bodies,printstep,T);
      }
      else
      {
      cout<<"No outfile produced - set printstep<=N for outfile"<<endl;
      }

      //forces at initial positions, f[0]
      force.calculateForces(&problem);

      //set forces at t_0 and print position and velocity at t_0
      for(int i = 0; i < Number_of_Bodies; i++)
        {
          if(printstep<=N)        //check print condition
          {
          printcurrent(0,problem.bodies[i]->position,problem.bodies[i]->velocity,obj->mass);
          }
          f_i[i] = problem.bodies[i]->force;

        }

      //Iterate, through t_1, t_2,...t_n-1, while printing time, position, velocity

      for(int yr=0;yr<T;yr++)
      {

      for(int step = 0; step < (N-1); step++)
        {
          for(int i = 0; i < Number_of_Bodies; i++)
            {
              obj = problem.bodies[i];  //set obj

              obj->position = obj->position +  obj->velocity*h +  f_i[i]*(hh/(2*obj->mass)); //update pos of obj


            }
          problem.resetForces();
          force.calculateForces(&problem);
          for(int i = 0; i < Number_of_Bodies; i++) f_iplus1[i] = problem.bodies[i]->force;  //set new force


        for(int i = 0; i < Number_of_Bodies; i++)
          {

            obj = problem.bodies[i];  //set obj
            obj->velocity = obj->velocity + ( f_iplus1[i] + f_i[i] ) * (h/(2*obj->mass)); //update velocity
            f_i[i] = f_iplus1[i]; //set force _-1 = force_i

            if((step+1)%printstep == 0) //check printperiod
            {
            printcurrent(yr+(step+1)*h,obj->position,obj->velocity,obj->mass);
            }
          }
        if(step == N-2&&(step+1)%printstep != 0) //check that printstep covers T final
        {
          cout<<"Did not print T final! - please adjust printstep"<<endl;
        }

        }
        if(printstep<=N)
        {

        }
      }
      printstop();
    }
