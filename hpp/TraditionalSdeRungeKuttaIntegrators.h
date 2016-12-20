#ifndef TRADITIONALSDERUNGEKUTTAINTEGRATORS_H
#define	TRADITIONALSDERUNGEKUTTAINTEGRATORS_H

#include "SdeRungeKuttaIntegrator.h"

class MilsteinTalay: public SdeRungeKuttaIntegrator
{
public:
    MilsteinTalay(Sde* sde, bool verb=false, bool dtadap=true,
                  Real atol=1e-2, Real rtol=1e-2, 
                  bool scalartol=true, int output_freq=0);
    virtual ~MilsteinTalay();
    
    virtual void step_general_noise(const Real t, const Real& h);
    virtual void step_diagonal_noise(const Real t, const Real& h);
    
protected:
    void update_n_stages_and_h(Real& h);
    
protected:
    Vector* integr5;
        
};

#endif	/* TRADITIONALSDERUNGEKUTTAINTEGRATORS_H */

