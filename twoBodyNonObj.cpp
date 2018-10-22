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
double massSun = 2e30;
double massEarth = 6e24/massSun;  //Mass of Earth in solar masses
double GM_O = 4.0*pi*pi;

double Force_Earth(double dim, double r_cubed_current);
double AngularMomentum(double rx, double ry, double mass, double vx, double vy);
double TotalEnergy(double mass, double rx, double ry, double vx, double vy);
void printEarthAttributes(double L, double E);
bool checkConvervation(double L, double E, double L_init, double E_init, double eps);


int main(int argc, char * argv[])
  {
    /*Discretization parameters*/
    int N = atoi(argv[1]);        //Number of integration points
    int printstep = atoi(argv[2]);        //Number of integration points
//    double eps = atof(argv[2]);
    double t_max = 1.0;    //Solve for 1 year
    double t_min = 0.0;    //Initial time
    double h = (t_max - t_min)/N;   //Step-Size
    double hh = h*h;
    int dimension = 2;  //Solving in two dimensions

    /*Settin initial conditions*/
    double Earth_x0, Earth_y0, Earth_v_x0, Earth_v_y0, r_cubed_old, r_cubed_new;
    double L, E_k, E_p; //Angular Momentum, Kinetic energy and potential energy

    Earth_x0   =  1.0;
    Earth_y0   =  0.0;
    Earth_v_x0 = 0.0;
    Earth_v_y0 =  2*pi;

    double* r = new double[dimension];
    double* v = new double[dimension];

    r[0] = Earth_x0;
    r[1] = Earth_y0;
    v[0] = Earth_v_x0;
    v[1] = Earth_v_y0;

    /*Earth angular momentum and total energy*/
    double Earth_L, Earth_E;
    Earth_L = AngularMomentum(r[0], r[1], massEarth, v[0], v[1]);
    Earth_E = TotalEnergy(massEarth, r[0], r[1], v[0], v[1]);

    /*Euler's method*/
    auto start1 = chrono::high_resolution_clock::now();  //Timing start

    double r_x_previous, r_y_previous, v_x_previous, v_y_previous;

    if(printstep<=N)
    {
    ofstream outfile;
    outfile.open("EarthSunPositionsEulerMethod.dat");
    outfile << Earth_x0 << " " << Earth_y0 << " " << Earth_L << " " << Earth_E << endl;
    }

    /*Current values of angular momentum, potential and kinetic energy */
    double L_curr, E_curr;

    for(int i = 1; i < N; i++)
      {
        r_cubed_old = (r[0]*r[0] + r[1]*r[1])*sqrt(r[0]*r[0] + r[1]*r[1]);

        r_x_previous = r[0];
        r_y_previous = r[1];
        v_x_previous = v[0];
        v_y_previous = v[1];

        r[0] = r_x_previous + h*v_x_previous;
        r[1] = r_y_previous + h*v_y_previous;

        v[0] = v_x_previous - h*Force_Earth(r_x_previous, r_cubed_old)/massEarth;
        v[1] = v_y_previous - h*Force_Earth(r_y_previous, r_cubed_old)/massEarth;

        L_curr = AngularMomentum(r[0], r[1], massEarth, v[0], v[1]);
        E_curr = TotalEnergy(massEarth, r[0], r[1], v[0], v[1]);

        /*Test conservation of angular momentum, potential and kinetic energy */
//        if(checkConvervation(L_curr, E_curr, Earth_L, Earth_E, eps) == 0)
  //        {
            /*Test failed*/
    //        break;
    //      }
        if((i)%printstep == 0) //check printperiod
        {
          outfile << r[0] << " " << r[1] << " " << L_curr << " " << E_curr << endl;
        }
      }

    auto stop1 = chrono::high_resolution_clock::now(); //Timing stop
    auto diff1 = stop1-start1;
    cout  << "N=" << N << " Runtime of Forward-Euler algorithm = " << chrono::duration <double,milli> (diff1).count() << "ms" << endl;

    if(printstep<=N)
    {
      outfile.close();
    }
    /*Reseting intitial conditions*/
    r[0] = Earth_x0;
    r[1] = Earth_y0;
    v[0] = Earth_v_x0;
    v[1] = Earth_v_y0;

    /*Velocity-Verlet Method*/
    auto start2 = chrono::high_resolution_clock::now();  //Timing start

    outfile.open("EarthSunPositionsVerletMethod.dat");
    outfile << Earth_x0 << " " << Earth_y0 << " " << Earth_L << " " << Earth_E << endl;

    for(int i = 0; i < N; i++)
      {
        r_cubed_old = (r[0]*r[0] + r[1]*r[1])*sqrt(r[0]*r[0] + r[1]*r[1]);

        r_x_previous = r[0];
        r_y_previous = r[1];
        v_x_previous = v[0];
        v_y_previous = v[1];

        r[0] = r_x_previous + h*v_x_previous - 0.5*hh*Force_Earth(r_x_previous, r_cubed_old)/massEarth;
        r[1] = r_y_previous + h*v_y_previous - 0.5*hh*Force_Earth(r_y_previous, r_cubed_old)/massEarth;

        r_cubed_new = (r[0]*r[0] + r[1]*r[1])*sqrt(r[0]*r[0] + r[1]*r[1]);

        v[0] = v_x_previous - 0.5*h*( Force_Earth( r[0], r_cubed_new ) + Force_Earth( r_x_previous, r_cubed_old ) )/massEarth;
        v[1] = v_y_previous - 0.5*h*( Force_Earth( r[1], r_cubed_new ) + Force_Earth( r_y_previous, r_cubed_old ) )/massEarth;

        //Computing angular momentum and total energy
        L_curr = AngularMomentum(r[0], r[1], massEarth, v[0], v[1]);
        E_curr = TotalEnergy(massEarth, r[0], r[1], v[0], v[1]);

        outfile << r[0] << " " << r[1] << " " << L_curr << " " << E_curr << endl;

      }

    auto stop2 = chrono::high_resolution_clock::now(); //Timing stop
    auto diff2 = stop2-start2;
    cout  << "N=" << N << " Runtime of Velocity-Verlet algorithm = " << chrono::duration <double,milli> (diff2).count() << "ms" << endl;

    return 0;
  }

double Force_Earth(double dim, double r_cubed_current)
    /*Calculates gravitational force of the sun on the earth*/
    {
      return GM_O*massEarth*(dim/r_cubed_current);
    }
double AngularMomentum(double rx, double ry, double mass, double vx, double vy)
    {
      double L, r, v, m;
      r = sqrt( rx*rx + ry*ry );
      v = sqrt( vx*vx + vy*vy);
      m = mass;
      L = r*m*v;
      return L;
    }
double TotalEnergy(double mass, double rx, double ry, double vx, double vy)
  {
    double U, K, E; //U: potential, K: kinetic, E: total energy
    U = -4*pi*pi*( mass / sqrt( rx*rx + ry*ry) );
    K = 0.5*mass*( vx*vx + vy*vy );
    E = U + K;
    return E;
  }
void printEarthAttributes(double L, double E)
      {
        cout << "--------------------------" << endl;
        cout << "|    Earth Attributes    |" << endl;
        cout << "--------------------------" << endl;
        cout << "Angular Momemntum = " << L << " Js" << endl;
        cout << "Total Energy  = " << E << " J"  << endl;
        cout << "--------------------------" << endl;
      }
bool checkConvervation(double L, double E, double L_init, double E_init, double eps)
        {
          double diff_L, diff_E;
          diff_L = L - L_init;
          diff_E = E - E_init;

          bool L_check, E_check;
          L_check = fabs(diff_L) > eps;
          E_check = fabs(diff_E) > eps;

          if( L_check == 1  || E_check == 1)
            {
              cout << "    Convervation test failed!   "  << endl;
              cout << "---------------------------------------"  << endl;
              cout << "|     |  Initial        |    Current    |" << endl;
              cout << "---------------------------------------"  << endl;
              cout << "|  L  |   "<<L_init<<"   | " << L << " |" << endl;
              cout << "---------------------------------------"  << endl;
              cout << "|  E  |   "<<E_init<<"  | " << E << " |" <<endl;
              cout << "---------------------------------------" << endl;

              return 0>2;
            }
          else
            {
              return 0<2;
            }

        }
