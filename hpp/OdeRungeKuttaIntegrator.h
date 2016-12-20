#ifndef ODERUNGEKUTTAINTEGRATOR_H
#define	ODERUNGEKUTTAINTEGRATOR_H

#include <iomanip>

#include "TimeIntegrator.h"
using namespace std;

class OdeRungeKuttaIntegrator: public TimeIntegrator
{
public:
    OdeRungeKuttaIntegrator(Ode* ode, bool verb=true, bool dtadap=true,
                         Real atol=1e-2, Real rtol=1e-2,
                         bool scalartol=true, int output_freq=0);
    virtual ~OdeRungeKuttaIntegrator();
                  
    bool integrate(Real& h); //solves problem
    
    int check_correctness(Real& h);             //checks if input parameters are ok
    virtual void print_info();                          //starting info
    virtual void print_statistics();                    //some statistics at the end
    Real get_dt_max();
    int get_n_f_eval();
    int getAcc_steps();
    int getRej_steps();
    void update_rho();      //updates spectral radius
    
    virtual void reinit();
    Real get_elapsed_time();
        
protected:

    virtual void step(const Real t, const Real& h) = 0;     //do the stages
    
    void accepted_step(Real& t, Real& h);
    void rejected_step(Real& t, Real& h);
    virtual void compute_hnew(Real &h);
    virtual void update_n_stages_and_h(Real& h) = 0; //given rho, computes number of stages
    void rho(Real& eigmax, int& iter);//power method approximating rho
    
    Real get_cpu_time();
    string getSolver();
  
public:
    //equation variables
    Vector*& yn;
    Vector* ynpu;
    Vector* eigenvector;
    Vector* integr[5]; //working vectors
    
    //integration settings
    bool dt_adaptivity;  ///<If false the time step is fixed.
    bool internal_rho;   ///<If true the ODE's spectral radius is computed by the internal power like method. Else a method rho is provived by the ODE class.
    bool scalar_tol;     ///<Tells if the tolerances are scalars or if every solution's component has its own.
    Real r_tol;          ///<Scalar relative tolerance.
    Real a_tol;          ///<Scalar absolute tolerance.
    bool verbose;        ///<If true prints info about the timestep.
    
    //integration statistics
    int max_rho;             ///<Maximal spectral radius computed.
    int min_rho;             ///<Minimal spectral radius computed.
    int max_s;               ///<Maximal number of stages used.
    Real dt_max;             ///<Maximal time step used.
    Real dt_min;             ///<Minimal time step used.
    int n_f_eval_rho;        ///<Number of righ hand side evaluations for the spectral radius computation in the power method.
    int n_f_eval;            ///<Number of righ hand side evaluations for time integration.
    int n_steps;             ///<Number of time steps.
    int acc_steps;           ///<Number of acepted steps.
    int rej_steps;           ///<Number of rejected steps.
    
protected:
    //step variables
    Real nu;             ///<Error is rejected if larger than nu
    Real eigmax;         ///<Spectral radius of current time step.
    Real uround;         ///<Minimal time step allowed.
    Real err;            ///<Local error estimation.
    Real errp;           ///<Local error estimated at the previus step.
    Real fac;            ///<Factor used in the choice of the new time step.
    Real facp;           ///<The factor used in the previous step.
    Real facmax;         ///<Maximal allowed factor.
    Real hp;             ///<Previus time step size.
    Real hnew;           ///<Next time step size.
    Real told;           ///<Starting time of previous time step. Used to detect multiple rejections.
    int s;               ///<Number of stages for RKC. Number of stages minus 2 for ROCK2.
    int sold;            ///<Previous number of stages.
    int nrej;            ///<Number of consecutive rejections.
    int nrho;            ///<Number of accpeted steps after last spectral radius estimation.
    bool last;           ///<Is true if the current step is the last one.
    bool reject;         ///<True if last step has been rejected.
    bool Astable;        ///<Tells is the RK method is A stable.
    string solver;       ///<Stores the name of the Runge-Kutta method
    int output_frequency;
    bool firsterrorestimation;
    bool rho_outdated;
    Real elapsed_time;
    string out_file;
        
};

class RungeKutta4: public OdeRungeKuttaIntegrator
{
public:
    RungeKutta4(Ode* ode, bool verb=true, int output_freq=0);
    virtual ~RungeKutta4(){};
    void step(const Real t, const Real& h);
    void update_n_stages_and_h(Real& h);
};

#endif	/* ODERUNGEKUTTAINTEGRATOR_H */