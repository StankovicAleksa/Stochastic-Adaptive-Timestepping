#include "MonteCarlo.h"
#include "OdeRungeKuttaIntegrator.h"
#include <iostream>
#include <fstream>
#include <omp.h>

MonteCarlo::MonteCarlo(int iter, Init init)
:initializer(init)
{
    iterations = iter;
}

void MonteCarlo::iterate(Real h)
{
    Real error=0.;
    Real val=0.;
    Real n_f_eval=0.;
    Real n_g_eval=0.;
    Real acc_steps=0.;
    Real rej_steps=0.;
    string solution_file;
    
    const int ofreq = iterations/10;
    const Real oh=h;
   
    TimeIntegrator* integrator=0;
    Sde* sde=0;
    BrownianMotion* W;
    
    cout<<setprecision(3);
    
    elapsed_time = get_cpu_time();
    
    cout<<"-----------------   Running Monte Carlo   ------------------"<<endl;
//    #pragma omp parallel private(W, integrator, sde, h)
    {
        initializer.initSde(sde);
        initializer.initTimeIntegrator(integrator, sde);
        initializer.initBrownianMotion(W, sde, integrator);
        sde->setW(W);
        
//        #pragma omp master
        {
            solution_file = sde->get_solution_file_name();
        }
        
//        #pragma omp for reduction(+:error,val,n_f_eval,n_g_eval,acc_steps,rej_steps) schedule(dynamic,1000)
        for(int i=1;i<=iterations;i++)
        {
            sde->reinit();
            integrator->reinit();
            h=oh;
            W->clear();

            integrator->integrate(h);

            
            error += (sde->phi()-sde->Exact_phi())/((Real)iterations);
            val += (sde->phi())/((Real)iterations);
            n_f_eval += integrator->get_n_f_eval()/((Real)iterations);
            n_g_eval += integrator->get_n_g_eval()/((Real)iterations);
            acc_steps += integrator->getAcc_steps()/((Real)iterations);
            rej_steps += integrator->getRej_steps()/((Real)iterations);
            
            if(i%ofreq==0)
                cout<<"Monte Carlo iterations at "<<100.*i/(Real)iterations<<"%"<<endl;
        }
        
        delete W;
        delete sde;
        delete integrator;
    }
    // end omp block
    cout<<"------------------------------------------------------------\n"<<endl;
    
    elapsed_time = get_cpu_time()-elapsed_time;
    
    if(error<0.) error=-error;
    
    this->error = error;
    this->val = val;
    this->n_f_eval = n_f_eval;
    this->n_g_eval = n_g_eval;
    this->acc_steps = acc_steps;
    this->rej_steps = rej_steps;
    
    ofstream out(solution_file+"_mt_error.txt", ofstream::out);
    out<<error<<", "<<val<<", "<<n_f_eval<<", "<<n_g_eval<<", "
            <<acc_steps<<", "<<rej_steps<<", "<<elapsed_time<<", "<<iterations;
    out.close();
   
}

Real MonteCarlo::get_cpu_time()
{
    return (Real)clock() / CLOCKS_PER_SEC;
}

void MonteCarlo::print_info()
{
    cout<<"\n-------------------   Monte Carlo Info   -------------------"<<endl;
    cout<<"Number of iterations: "<<iterations<<endl;
    cout<<"SDE solver: "<<initializer.getRk_name()<<endl;
    cout<<"Test number: "<<initializer.getNtest()<<endl;
    cout<<"Time step size: "<<initializer.getDt()<<endl;
    cout<<"Time-step adaptivity "<<(initializer.isDtadap() ? "enabled.":"disabled.")<<endl;
    cout<<"Brownian motion: "<<(initializer.isContinuous()? "continuous.":"discrete.")<<endl;
    cout<<"Absolute tolerance = "<<initializer.getAtol()<<endl;
    cout<<"Relative tolerance = "<<initializer.getRtol()<<endl;
    cout<<"------------------------------------------------------------\n"<<endl;
}

void MonteCarlo::print_statistics()
{
    cout<<setprecision(16);
    cout<<"----------------   Monte Carlo Statistics   ----------------"<<endl;
    cout<<"Running Time: "<<elapsed_time<<" seconds."<<endl;
    cout<<"Mean of f evaluations: "<<n_f_eval<<endl;
    cout<<"Mean of g evaluations: "<<n_g_eval<<endl;
    cout<<"Mean of accepted steps: "<<acc_steps<<endl;
    cout<<"Mean of rejected steps: "<<rej_steps<<endl;
    cout<<"E(phi(X)) = "<<val<<endl;
    cout<<"Error on E(phi(X)): "<<error<<endl;
    cout<<"------------------------------------------------------------\n"<<endl;
}