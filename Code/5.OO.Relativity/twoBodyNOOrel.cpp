
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <chrono>
#include <ratio>
#include <cstdlib>

using namespace std;

const double pi = 3.14159265358979323846;
const double c = 63239.7263;

/*Solar System attributes*/
double massSun = 2.0e30;
double massMercury = 3.3e23/massSun;  //Mass of Mercury in solar masses
double GM_O = 4.0*pi*pi;

double AngularMomentumSquared(double rx, double ry, double rz, double vx, double vy, double vz);
double Force_Mercury(double r_dim, double r_norm, double L_squared);

int main(int argc, char * argv[])
  {
    /*Discretization parameters*/
    int N = atoi(argv[1]);        //Number of integration points
    int printstep = atoi(argv[2]);        //Number of integration points
//    double eps = atof(argv[2]);
    double t_max = atoi(argv[3]);    //Solve for 1 year
    double h = t_max/N;   //Step-Size
    double hh = h*h;
    int dimension = 3;  //Solving in two dimensions

    /*Settin initial conditions*/
    double Mercury_x0, Mercury_y0, Mercury_z0, Mercury_v_x0, Mercury_v_y0, Mercury_v_z0;
    double L_squared_old, L_squared_new;

    Mercury_x0   =   0.3075;
    Mercury_y0   =   0.0;
    Mercury_z0   =  0.0;
    Mercury_v_x0 =   0.0;
    Mercury_v_y0 =  12.44;
    Mercury_v_z0 = 0.0;

    double Mercury_perihelion_distance = 0.3075;

    double* r = new double[dimension];
    double* v = new double[dimension];

    r[0] = Mercury_x0;
    r[1] = Mercury_y0;
    r[2] = Mercury_z0;
    v[0] = Mercury_v_x0;
    v[1] = Mercury_v_y0;
    v[2] = Mercury_v_z0;

    /*Velocity-Verlet Method*/
    double r_x_previous, r_y_previous, r_z_previous, v_x_previous, v_y_previous, v_z_previous, r_norm_old, r_norm_new;

    int perihelion_counter = 0;
    double tolerance = 1.e-6;

    ofstream outfile2;
    if(printstep<=N)
      {
        outfile2.open("MercurySunPositionsVerletMethod.dat");
        outfile2 << Mercury_x0 << " " << Mercury_y0 << " " << Mercury_z0 << endl;
      }

    for(int i = 1; i < N; i++)
      {
        r_norm_old = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);

        r_x_previous = r[0];
        r_y_previous = r[1];
        r_z_previous = r[2];
        v_x_previous = v[0];
        v_y_previous = v[1];
        v_z_previous = v[2];

        L_squared_old = AngularMomentumSquared(r_x_previous, r_y_previous, r_z_previous, v_x_previous, v_y_previous, v_z_previous);

        r[0] = r_x_previous + h*v_x_previous - 0.5*hh*Force_Mercury(r_x_previous, r_norm_old, L_squared_old)/massMercury;
        r[1] = r_y_previous + h*v_y_previous - 0.5*hh*Force_Mercury(r_y_previous, r_norm_old, L_squared_old)/massMercury;
        r[2] = r_z_previous + h*v_z_previous - 0.5*hh*Force_Mercury(r_z_previous, r_norm_old, L_squared_old)/massMercury;

        if(fabs(r_norm_old - Mercury_perihelion_distance) < tolerance)
          {
            perihelion_counter++;
          }

        r_norm_new = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        L_squared_new = AngularMomentumSquared(r[0], r[1], r[2], v[0], v[1], v[2]);

        v[0] = v_x_previous - 0.5*h*( Force_Mercury( r[0], r_norm_new, L_squared_new) + Force_Mercury( r_x_previous, r_norm_old, L_squared_old) )/massMercury;
        v[1] = v_y_previous - 0.5*h*( Force_Mercury( r[1], r_norm_new, L_squared_new) + Force_Mercury( r_y_previous, r_norm_old, L_squared_old) )/massMercury;
        v[2] = v_z_previous - 0.5*h*( Force_Mercury( r[2], r_norm_new, L_squared_new) + Force_Mercury( r_z_previous, r_norm_old, L_squared_old) )/massMercury;



        if((i)%printstep == 0) //check printperiod
          {
            outfile2 << r[0] << " " << r[1] << " " << r[2] << endl;
          }
      }

    if(printstep<=N)
      {
        outfile2.close();
      }
    cout << "perihelion_counter = " << perihelion_counter << endl;
    return 0;
  }

double AngularMomentumSquared(double rx, double ry, double rz, double vx, double vy, double vz)
    {
      double L_x, L_y, L_z;
      L_x = ry*vz - rz*vy;
      L_y = rz*vx - rx*vz;
      L_z = rx*vy - ry*vx;

      return L_x*L_x + L_y*L_y + L_z*L_z;
    }

double Force_Mercury(double r_dim, double r_norm, double L_squared)
    /*Calculates gravitational force of the sun on the Mercury*/
    {
      double r_squared = pow(r_norm, 2);
      double r_cubed = pow(r_norm, 3);

      double Newtonian = (GM_O*massMercury)/r_cubed;
//      double relativistic_correction = 1.0 + (3.0*L_squared)/(r_squared*c*c);
      double relativistic_correction = 1;
      return Newtonian*relativistic_correction*r_dim;

    }
