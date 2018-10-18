#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <chrono>
#include <ratio>
#include <ctime>

using namespace std;

const double pi = 3.14159265358979323846;

/*Solar System attributes*/
double massSun = 2e30;
double massEarth = 6e24/massSun;  //Mass of Earth in solar masses
double GM_O = 4.0*pi*pi;

double Force_Earth(double dim, double r_cubed_current);
double AngularMomentum_Earth(double rx, double ry, double mass, double vx, double vy);
double PotentialEnergy_Earth(double mass, double rx, double ry);
double KineticEnergy_Earth(double mass, double vx, double vy);
void printEarthAttributes(double L, double U, double K);
bool checkConvervation(double L, double U, double K, double L_init, double U_init, double K_init);


int main()
  {
    /*Discretization parameters*/
    int N = 100000;        //Number of integration points
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

    /*Earth angular momentum, potential and kinetic energy*/
    double Earth_L, Earth_U, Earth_K;
    Earth_L = AngularMomentum_Earth(r[0], r[1], massEarth, v[0], v[1]);
    Earth_U = PotentialEnergy_Earth(massEarth, r[0], r[1]);
    Earth_K = KineticEnergy_Earth(massEarth, v[0], v[1]);

    /*Euler's method*/
    auto start1 = chrono::high_resolution_clock::now();  //Timing start

    double r_x_previous, r_y_previous, v_x_previous, v_y_previous;

    ofstream outfile;
    outfile.open("EarthSunPositionsEulerMethod.dat");
    outfile << Earth_x0 << " " << Earth_y0 << endl;

    /*Current values of angular momentum, potential and kinetic energy */
    double L_curr, U_curr, K_curr;

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

        outfile << r[0] << " " << r[1] << endl;

        /*Test conservation of angular momentum, potential and kinetic energy */
        L_curr = AngularMomentum_Earth(r[0], r[1], massEarth, v[0], v[1]);
        U_curr = PotentialEnergy_Earth(massEarth, r[0], r[1]);
        K_curr = KineticEnergy_Earth(massEarth, v[0], v[1]);


        if(checkConvervation(L_curr, U_curr, K_curr, Earth_L, Earth_U, Earth_K) == 0)
          {
            /*Test failed*/
            break;
          }

      }

    auto stop1 = chrono::high_resolution_clock::now(); //Timing stop
    auto diff1 = stop1-start1;
    cout  << "N=" << N << " Runtime of Forward-Euler algorithm = " << chrono::duration <double,milli> (diff1).count() << "ms" << endl;

    outfile.close();

    /*Reseting intitial conditions*/
    r[0] = Earth_x0;
    r[1] = Earth_y0;
    v[0] = Earth_v_x0;
    v[1] = Earth_v_y0;

    /*Velocity-Verlet Method*/
    auto start2 = chrono::high_resolution_clock::now();  //Timing start

    outfile.open("EarthSunPositionsVerletMethod.dat");
    outfile << Earth_x0 << " " << Earth_y0 << endl;
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
        outfile << r[0] << " " << r[1] << endl;

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
double AngularMomentum_Earth(double rx, double ry, double mass, double vx, double vy)
    {
      double L, r, v, m;
      r = sqrt( rx*rx + ry*ry );
      v = sqrt( vx*vx + vy*vy);
      m = mass;
      L = r*m*v;
      return L;
    }
double PotentialEnergy_Earth(double mass, double rx, double ry)
    {
      double U;
      U = -4*pi*pi*( mass / sqrt( rx*rx + ry*ry) );
    }
double KineticEnergy_Earth(double mass, double vx, double vy)
    {
      double K;
      K = 0.5*mass*( vx*vx + vy*vy );
    }
void printEarthAttributes(double L, double U, double K)
      {
        cout << "--------------------------" << endl;
        cout << "|    Earth Attributes    |" << endl;
        cout << "--------------------------" << endl;
        cout << "Angular Momemntum = " << L << " Js" << endl;
        cout << "Potential Energy  = " << U << " J"  << endl;
        cout << "Kinetic Energy    = " << K << " J"  << endl;
        cout << "--------------------------" << endl;
      }
bool checkConvervation(double L, double U, double K, double L_init, double U_init, double K_init)
        {
          double eps = 1.e-7;
          double diff_L, diff_U, diff_K;
          diff_L = L - L_init;
          diff_U = U - U_init;
          diff_K = K - K_init;

          bool L_check, U_check, K_check;
          L_check = fabs(diff_L) > eps;
          U_check = fabs(diff_U) > eps;
          K_check = fabs(diff_K) > eps;

          if( L_check == 1  || U_check == 1 || K_check == 1)
            {
              cout << "    Convervation test failed!   "  << endl;
              cout << "---------------------------------------"  << endl;
              cout << "|     |  Initial        |    Current    |" << endl;
              cout << "---------------------------------------"  << endl;
              cout << "|  L  |   "<<L_init<<"   | " << L << " |" << endl;
              cout << "---------------------------------------"  << endl;
              cout << "|  U  |   "<<U_init<<"  | " << U << " |" <<endl;
              cout << "---------------------------------------"  << endl;
              cout << "|  K  |   "<<K_init<<"   | " << K << " |" <<endl;
              cout << "---------------------------------------" << endl;

              return 0>2;
            }
          else
            {
              return 0<2;
            }

        }
