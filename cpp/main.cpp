#include "OdeRungeKuttaIntegrator.h"
#include "OdeProblems.h"
#include "Init.h"
#include "MonteCarlo.h"

int main(int argc, char** argv)
{
    // default parameters
    int ntest = 4;
    Real dt = 1.00/16;
    Real rtol = 1e-3;
    Real atol = rtol;
    bool dt_adaptivity = false;
    string solver = "SROCK2";
    //string solver = "ROCK2";
    int output_frequency = 0;
    bool verbose = true;
    string output_file = "solution";
    int iter = 100000000; //monte carlo iterations
    bool continuous_brown = false;
    bool deterministic=false;    
    bool end;
    
    // Initialize parameters
    Init init; // Class for initializing stuff
    GetPot command_line(argc, argv); //command line parser
    //read input parameters and eventually modify defaults
    end=init.read_command_line(command_line, ntest, dt, rtol, atol, dt_adaptivity, 
                      solver, output_frequency, output_file, 
                      iter, continuous_brown,deterministic,verbose);
    if(end)
        return 0;
    
    Ode* ode=0;
    Sde* sde=0;
    TimeIntegrator* integrator=0;
    BrownianMotion* W=0;
    
    if(deterministic)
    {
        init.initOde(ode);
        if(ode==0)
            return 0;
        init.initTimeIntegrator(integrator, ode);
        if(integrator==0)
            return 0;
    }
    else if(iter==1)
    {
        init.initSde(sde);
        if(sde==0)
            return 0;
        init.initTimeIntegrator(integrator, sde);
        if(integrator==0)
            return 0;
        init.initBrownianMotion(W, sde, integrator);
        sde->setW(W);
    }
    
    if(iter ==1)
    {
        integrator->print_info();

        integrator->integrate(dt);

        integrator->print_statistics();
    }
    else
    {
        MonteCarlo mt(iter,init);
        mt.print_info();
        mt.iterate(dt);
        mt.print_statistics();
    }

    if(deterministic)
    {
        delete ode;
        delete integrator;
    }
    else if(iter==1)
    {
        delete sde;
        delete W;
        delete integrator;
    }
    
    
    return 0;
}
