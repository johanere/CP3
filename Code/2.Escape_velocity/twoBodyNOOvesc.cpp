
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <chrono>
#include <ratio>
#include <cstdlib>
#include <ctime>

using namespace std;

const double pi = 3.14159265358979323846;

/*Solar System attributes*/
double massSun = 2.0e30;
double massEarth = 6e24/massSun;  //Mass of Earth in solar masses
double GM_O = 4.0*pi*pi;

double Force_Earth(double r_dim, double r_norm, double);

int main(int argc, char * argv[])
  {
    /*Discretization parameters*/
    int N = atoi(argv[1]);        //Number of integration points
    int printstep = atoi(argv[2]);        //Number of integration points
//    double eps = atof(argv[2]);
    double t_max = atoi(argv[3]);    //Solve for 1 year
    double h = t_max/N;   //Step-Size
    double hh = h*h;
    int dimension = 2;  //Solving in two dimensions

    double beta = atof(argv[4]);  //power of norm of distance between sun and earth

    /*Settin initial conditions*/
    double Earth_x0, Earth_y0, Earth_v_x0, Earth_v_y0;

    double v_esc_analytical = sqrt(2*GM_O);   //Analytical escape velocity Earth 8.88577 AU/yr

    Earth_x0   =   1.0;
    Earth_y0   =   0.0;
    Earth_v_x0 =   0.0;
    Earth_v_y0 =  2*pi+1.0;

    cout << "Testing v_0 = " << Earth_v_y0 << endl;   //Print current initial velocity of Earth

    double* r = new double[dimension];
    double* v = new double[dimension];

    r[0] = Earth_x0;
    r[1] = Earth_y0;
    v[0] = Earth_v_x0;
    v[1] = Earth_v_y0;


    /*Velocity-Verlet Method*/
    auto start2 = chrono::high_resolution_clock::now();  //Timing start
    double r_x_previous, r_y_previous, v_x_previous, v_y_previous, r_norm_old, r_norm_new;


    ofstream outfile2;
    if(printstep<=N)
      {
        outfile2.open("EarthSunPositionsVerletMethod.dat");
        outfile2 << Earth_x0 << " " << Earth_y0 << endl;
      }

    for(int i = 1; i < N; i++)
      {
        r_norm_old = sqrt(r[0]*r[0] + r[1]*r[1]);

        r_x_previous = r[0];
        r_y_previous = r[1];
        v_x_previous = v[0];
        v_y_previous = v[1];

        r[0] = r_x_previous + h*v_x_previous - 0.5*hh*Force_Earth(r_x_previous, r_norm_old, beta)/massEarth;
        r[1] = r_y_previous + h*v_y_previous - 0.5*hh*Force_Earth(r_y_previous, r_norm_old, beta)/massEarth;

        r_norm_new = sqrt(r[0]*r[0] + r[1]*r[1]);

        v[0] = v_x_previous - 0.5*h*( Force_Earth( r[0], r_norm_new, beta ) + Force_Earth( r_x_previous, r_norm_old, beta ) )/massEarth;
        v[1] = v_y_previous - 0.5*h*( Force_Earth( r[1], r_norm_new, beta ) + Force_Earth( r_y_previous, r_norm_old, beta ) )/massEarth;

        if((i)%printstep == 0) //check printperiod
          {
            outfile2 << r[0] << " " << r[1] << endl;
          }
      }

    auto stop2 = chrono::high_resolution_clock::now(); //Timing stop
    auto diff2 = stop2-start2;
    cout  << "N=" << N << " Runtime of Velocity-Verlet algorithm = " << chrono::duration <double,milli> (diff2).count() << "ms" << endl;

    if(printstep<=N)
      {
        outfile2.close();
      }

    return 0;
  }

double Force_Earth(double r_dim, double r_norm, double beta=2)
    /*Calculates gravitational force of the sun on the earth*/
    {
      return (GM_O*massEarth)/pow(r_norm, beta+1)*r_dim;
    }
