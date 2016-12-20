#ifndef TIMEINTEGRATOR_H
#define	TIMEINTEGRATOR_H

#include "Ode.h"

class TimeIntegrator
{
public:
    TimeIntegrator(Ode* ODE) {ode=ODE;};
    virtual ~TimeIntegrator(){};
 
    virtual bool integrate(Real& h) =0; //compute next yn 
    
    virtual int check_correctness(Real& h) =0;             //checks if input parameters are ok
    virtual void print_info() =0;                          //starting info
    virtual void print_statistics() =0;                    //some statistics at the end
    virtual Real get_dt_max() =0;
    virtual bool needDoubleIntegral(){return false;};
    virtual void reinit() = 0;
    virtual int get_n_f_eval(){return 0;};
    virtual int get_n_g_eval(){return 0;};
    virtual int getAcc_steps(){return 0;};
    virtual int getRej_steps(){return 0;};
protected:
    Ode* ode;
};

#endif	/* TIMEINTEGRATOR_H */

