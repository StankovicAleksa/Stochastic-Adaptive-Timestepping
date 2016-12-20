#ifndef SDESTABILIZEDINTEGRATORS_H
#define	SDESTABILIZEDINTEGRATORS_H

#include "OdeStabilizedIntegrators.h"
#include "SdeRungeKuttaIntegrator.h"

class SROCK2: public SdeRungeKuttaIntegrator, public DROCK2
{
public:
    SROCK2(Sde* sde, bool verb=true, bool dtadap=true, 
           Real atol=1e-2, Real rtol=1e-2, bool scalartol=true,
           int output_freq=0);
    virtual ~SROCK2();
    
    void step(const Real t, const Real& h);
    void step_diagonal_noise(const Real t, const Real& h);
    void step_general_noise(const Real t, const Real& h);
    
protected:
    Vector* integr5;
};
#endif	/* SDESTABILIZEDINTEGRATORS_H */

