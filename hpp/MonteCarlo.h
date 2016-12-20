#ifndef MONTECARLO_H
#define	MONTECARLO_H

#include "Init.h"
#include "SdeRungeKuttaIntegrator.h"
#include "Sde.h"

class MonteCarlo
{
public:
    MonteCarlo(int iter, Init init);
    void iterate(Real h);
    
    Real get_cpu_time();
    void print_info();
    void print_statistics();
    
protected:
    Init initializer;
    int iterations;
    Real error;
    Real val;
    Real n_f_eval;
    Real n_g_eval;
    Real acc_steps;
    Real rej_steps;
    Real elapsed_time;
};


#endif	/* MONTECARLO_H */

