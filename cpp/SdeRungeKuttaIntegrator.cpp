#include <string>
#include <fstream>

#include "SdeRungeKuttaIntegrator.h"

SdeRungeKuttaIntegrator::SdeRungeKuttaIntegrator(Sde* sde, bool verb, 
                         bool dtadap, Real atol, Real rtol, 
                         bool scalartol, int output_freq)
:OdeRungeKuttaIntegrator(sde,verb,dtadap,atol,rtol,scalartol,output_freq)
{
    reinit();
    
    if(sde->isDiagonal())
        Gv.resize(sde->system_size());
    else
        G.resize(sde->system_size(),sde->brownian_size());
    
    needDoubleInt = true;
}

SdeRungeKuttaIntegrator::~SdeRungeKuttaIntegrator()
{
    
}

void SdeRungeKuttaIntegrator::step(const Real t, const Real& h)
{
    Sde* sde = dynamic_cast<Sde*>(ode);
    if(sde->isDiagonal())
        this->step_diagonal_noise(t,h);
    else
        this->step_general_noise(t,h);
}

void SdeRungeKuttaIntegrator::reinit()
{
    OdeRungeKuttaIntegrator::reinit();
    n_g_eval = 0;
}

bool SdeRungeKuttaIntegrator::needDoubleIntegral()
{
    return needDoubleInt;
}

void SdeRungeKuttaIntegrator::print_info()
{
    /**
     * Information about the chosen RungeKuttaIntegrator settings.
     */
    // Printing integration parameters
    
    cout<<"\n-------------------   Sde Solver Info   --------------------"<<endl;
    cout<<"SDE solver: "<<solver<<endl;
    cout<<"Time-step adaptivity "<<(dt_adaptivity ? "enabled.":"disabled.")<<endl;
    cout<<"Spectral radius computed "<<(internal_rho ? "internally.":"externally.")<<endl;
    cout<<"Absolute tolerance = "<<a_tol<<endl;
    cout<<"Relative tolerance = "<<r_tol<<endl;
    cout<<"------------------------------------------------------------\n"<<endl;
}

void SdeRungeKuttaIntegrator::print_statistics()
{
    /**
     * Some statistics about the time integration.
     */
    
    cout<<"\n----------------   Integration Statistics   ----------------"<<endl;
    
    if(Astable)
        cout<<"The spectral radius has not been computed, "<<solver<<" is A-stable."<<endl;
    else
    {
        cout<<"Max estimation of the spectral radius: "<<max_rho<<endl;
        cout<<"Min estimation of the spectral radius: "<<min_rho<<endl;
    }
    cout<<"Number of f eval. for the spectr. radius = "<<n_f_eval_rho<<endl;
    cout<<"Max number of stages used: "<<max_s<<endl;
    cout<<"Number of f evaluations = "<<n_f_eval<<endl;
    cout<<"Number of g evaluations = "<<n_g_eval<<endl;
    cout<<"Maximal time step used: "<<dt_max<<endl;
    cout<<"Steps: "<<n_steps<<endl;
    cout<<"Accepted steps: "<<acc_steps<<endl;
    cout<<"Rejected steps: "<<rej_steps<<endl; 
    cout<<"Elapsed time: "<<elapsed_time<<endl;
    cout<<"------------------------------------------------------------\n"<<endl;
    
    ofstream out(out_file+string("_statistics.txt"), ofstream::out);
    out<<setprecision(16);
    out<<elapsed_time<<", ";
    out<<n_f_eval<<", ";
    out<<n_g_eval<<", ";
    out<<max_s<<", ";
    out<<max_rho<<", ";
    out<<n_steps<<", ";
    out<<acc_steps<<", ";
    out<<rej_steps;
    out.close();
}

int SdeRungeKuttaIntegrator::get_n_g_eval()
{
    return n_g_eval;
}
