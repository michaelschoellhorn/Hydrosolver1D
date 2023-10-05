#include<iostream>
#include<fstream>
#include<vector>
#include<ctime>
#include<cmath>
#include<string>
#include<sstream>

using namespace std;

/*
Compile and run this script first!
This script solves the euler equations with a classic hydrodynamical solver.
*/

// defining matrix and vector 
typedef vector<double> Vec;
typedef vector<vector<double> > Mat;


class simulation{
    public:
     int N;
     int Ngc;
     Mat Q;
     Mat F; // F_i+1/2
     Vec x; //x in cellmid
     Vec p;
     double delta_x;
     double delta_t;
     double gam;
     string init_flag;

    public:
     simulation(int, int, double, string);
     void iterate(const double, const double, const double, const bool, double);
     void bound();
     void advection();
     double source();
     double source(double); //overload of source to include visc
        
};

/*
Main constructor for simulation class.
*/
simulation::simulation(int Ngridpoints, int Ngc, double gam, string init_flag){

    //set constants
    N = Ngridpoints+2*Ngc;
    this -> Ngc = Ngc;
    this -> gam = gam;
    this -> init_flag = init_flag; //reuse as savedata name
    delta_x = 100.0/(Ngridpoints-1);
    
    //init Matices and vectors
    Q = Mat(3, Vec(N, 0.0));
    F = Mat(3, Vec(N, 0.0));
    x = Vec(N);
    p = Vec(N, 0);
    
    //set x
    for (Vec::size_type i = 0; i!=x.size(); i++){
        x[i] = 0.0+(i-2.0)*delta_x;
        //cout << x[i] << " ";
    }

    //choosing init for Q
    if (init_flag == "gauss")
    {
        for (Vec::size_type i = 0; i!=x.size(); i++){
            Q[0][i] = 1.0 + 0.3*exp(-pow((x[i]-50.0), 2.0)/10.0);
            //cout << Q[0][i];
        }
    }
    else if (init_flag == "jump")
    {
        for (Vec::size_type i = 0; i!=x.size(); i++)
        {
            Q[0][i] = 1.0*(x[i] <= 50.0) + 1.0;
            Q[2][i] = Q[0][i]*1.0;
        }
    }
    else{
        cout << "Unknown initialization flag." <<endl;
    }
}


/*
Method to iterate the system.
*/
void simulation::iterate(const double t_end, const double t_save_step, const double cfl, const bool is_isotherm, double visc_factor){
    double t = 0.0;
    double t_save = t;
    double nu_max = 0.0;

    //calculate nu_max(t=0)
    if (is_isotherm)
    {
        double temp = 0.0;
        double temp2 = 0.0;
        for (Vec::size_type i = 0; i < N; i++)
        {
            temp2 = abs(Q[1][i]/(Q[0][i]));
            if (temp2 > temp){
                temp = temp2;
            }
            
        }
        nu_max = temp + 1.0;
    }
    else
    {
        double temp = 0.0;
        double temp2 = 0.0;
        for (Vec::size_type i = 0; i < N; i++)
        {
            temp2 = sqrt(gam*(gam-1.0)*(Q[2][i]/(Q[0][i]+1E-16) - 0.5*pow(Q[1][i]/(Q[0][i] + 1E-16), 2.0))) + abs(Q[1][i]/(Q[0][i] + 1E-16)); // c_s + |u|
            if (temp2 > temp){
                temp = temp2;
            }
            
        }
        nu_max = temp;
    }
    
    //open filestreams for data output:
    string filenames[3] = {"q1" + init_flag + to_string(visc_factor) + ".txt", "q2" + init_flag + to_string(visc_factor) + ".txt", "q3" + init_flag + to_string(visc_factor)+ ".txt"};
    ofstream outputfile[3];
    for (int i = 0; i < 3; i++)
    {
        outputfile[i].open(filenames[i]);
    }
    
    //loop over t
    while (t < t_end){
        delta_t = cfl*delta_x/nu_max;
        if ((t+delta_t) > t_end){
            delta_t = t_end - t;
        }
        bound();
        advection();
        bound();
        if (is_isotherm)
        {
            nu_max = source();
        }
        else
        {
            nu_max = source(visc_factor);
        }
        //write Q to ofstreams
        if (t > t_save+t_save_step)
        {
            for (Vec::size_type i = Ngc; i < N-Ngc; i++)
            {
                if (i!=Ngc){
                    outputfile[0] << ", ";
                    outputfile[1] << ", ";
                    outputfile[2] << ", ";
                }
                outputfile[0] << Q[0][i];
                outputfile[1] << (Q[1][i]/Q[0][i]);
                outputfile[2] << (Q[2][i]/Q[0][i] - 0.5*pow(Q[1][i], 2.0)/pow(Q[0][i], 2.0));
                if (i == N-1-Ngc){
                    outputfile[0] << "\n";
                    outputfile[1] << "\n";
                    outputfile[2] << "\n";
                }
            }
            t_save = t;
        }
        t+= delta_t;
    }

    //close ofstreams for data output
    for (int i = 0; i < 3; i++)
    {
        outputfile[i].close();
    }  
}


void simulation::advection(){
    //calculate u_i+1/2 and generate flux
    double u_half;
    for (Vec::size_type i = 1; i < N-Ngc+1; i++)
    {
        u_half = 0.5 * (Q[1][i]/(Q[0][i]) + Q[1][i+1]/(Q[0][i+1]));
        if (u_half >= 0)
        {
            F[0][i] = Q[0][i]*u_half;
            F[1][i] = Q[1][i]*u_half;
            F[2][i] = Q[2][i]*u_half;
        }
        else
        {
            F[0][i] = Q[0][i+1]*u_half;
            F[1][i] = Q[1][i+1]*u_half;
            F[2][i] = Q[2][i+1]*u_half;
        }
    }
    //calculate Q_half
    for (size_t k = 0; k < 3; k++)
    {
        for (Vec::size_type i = 2; i < N-Ngc+1; i++)
        {
            Q[k][i] -=delta_t/delta_x*(F[k][i]-F[k][i-1]);
        }
    }
}


/*
This method is used for the non-isothermal case. It includes von Neumann-Richtmeyer viscosity.
----------------------
double visc: 
viscosity factor, usually between 2.0, 3.0. 
Set to 0 for the non-viscous case.
*/
double simulation::source(double visc){
    double nu = 0.0;
    //calculate pressure
    Mat Qtemp(2, Vec(N, 0.0));
    for (Vec::size_type i = Ngc-1; i < N-Ngc+1 ; i++)
    {
        p[i] = (gam - 1.0)*(Q[2][i]-0.5*pow(Q[1][i], 2.0)/(Q[0][i]+1E-16));
        if (Q[1][i+1]/Q[0][i+1] <= Q[1][i-1]/Q[0][i-1]){
            p[i] += 0.25*pow(visc, 2.0) * pow((Q[1][i+1]/Q[0][i+1] - Q[1][i-1]/Q[0][i-1]), 2.0)*Q[0][i];
        }
    }
    //calculate Qnew and nu_max
    for (Vec::size_type i = Ngc; i < N-Ngc; i++)
    {
        Qtemp[0][i] = Q[1][i] - delta_t/(2.0*delta_x)*(p[i+1]-p[i-1]); //rho u
        Qtemp[1][i] = Q[2][i] - delta_t/(2.0*delta_x)*(p[i+1]*Q[1][i+1]/(Q[0][i+1]+1E-16) - p[i-1]*Q[1][i-1]/(Q[0][i-1]+1E-16)); //rho eps
        double temp = sqrt(gam*p[i]/(Q[0][i]+1E-16)) + abs(Qtemp[0][i]/(Q[0][i]+1E-16)); // c_s[i] + |u_[i]|
        if (temp>nu){
            nu = temp; //carry max nu
        }
    }
    //set Q to Qnew
    Q[1] = Qtemp[0];
    Q[2] = Qtemp[1];
    return nu;
}


/*
This method handles the source term in the isothermal case. Sound speed c_s is given by 1.0.
*/
double simulation::source(){
    //calculate pressure
    double nu = 0.0;
    for (Vec::size_type i = Ngc-1; i < N-Ngc+1 ; i++)
    {
        p[i] = pow(1.0, 2.0)*(Q[0][i]);
    }
    //calculate Qnew and updating Q
    double temp = 0.0;
    for (Vec::size_type i = Ngc; i < N-Ngc; i++)
    {
        Q[1][i] = Q[1][i] - delta_t/(2.0*delta_x)*(p[i+1] - p[i-1]); //rho u
        temp = 1.0 + abs(Q[1][i]/(Q[0][i])); // nu atm
        if (temp>nu){
            nu = temp; //set nu max
        }
    }
    return nu;
}


/*
Forces refective boundary conditions
*/
void simulation::bound(){
    /* 
    TODO: Does only accomodate for first order scheme atm
    */
    Q[0][Ngc-1] = Q[0][Ngc];
    Q[0][N-Ngc] = Q[0][N-Ngc-1];

    Q[1][Ngc-1] = -Q[1][Ngc];
    Q[1][N-Ngc] = -Q[1][N-Ngc-1];

    Q[2][Ngc-1] = Q[2][Ngc];
    Q[2][N-Ngc] = Q[2][N-Ngc-1];
}





int main(){
    //1a)
    simulation s1(500, 2, 1.4, "gauss");
    s1.iterate(100.0, 0.2, 0.5, true, 0.0);

    //1b)
    simulation s2(500, 2, 1.4, "jump");
    s2.iterate(40.0, 0.2, 0.5, false, 0.0);

    //1c)
    simulation s3(500, 2, 1.4, "jump");
    s3.iterate(40.0, 0.2, 0.5, false, 3.0);

}

