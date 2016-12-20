#ifndef ODESTABILIZEDINTEGRATORS_H
#define	ODESTABILIZEDINTEGRATORS_H

#include "OdeRungeKuttaIntegrator.h"

class RKC: public virtual OdeRungeKuttaIntegrator
{
public:
    RKC(Ode* ode, bool verb=true, bool dtadap=true, 
        Real atol=1e-2, Real rtol=1e-2, bool scalartol=true,
        int output_freq=0);
    virtual ~RKC();
    
protected:    
    void step(const Real t, const Real& h);
    
    void update_n_stages_and_h(Real& h);
    void init_coeffs(Real *w, int s, Real *bj, Real *thj, Real *zj, Real *dzj, 
                     Real *d2zj, Real& kappa);
    void update_coeffs(Real *w, Real *bj, Real *zj, Real *dzj, Real *d2zj, 
                       Real& mu, Real& nu, Real& kappa, Real& ajm1, Real* thj);
    void shift_coeffs(Real *bj, Real *zj, Real *dzj, Real *d2zj, Real* thj);
};


class ROCK2: public virtual OdeRungeKuttaIntegrator
{
public:
    ROCK2(Ode* ode, bool verb=true, bool dtadap=true, 
          Real atol=1e-2, Real rtol=1e-2, bool scalartol=true,
          int output_freq=0);
    virtual ~ROCK2();
    
protected:

    virtual void step(const Real t, const Real& h);

    virtual void update_n_stages_and_h(Real& h);
    void mdegr(int& mdeg, int mp[]);
    
 
protected:
    int mp[2];  ///<It is used in order to find the algorithm precomputed coefficients in the tables.
    
    static int ms[46];      ///<Array of coefficients.
    static Real fp1[46];    ///<Array of coefficients.
    static Real fp2[46];    ///<Array of coefficients.
    static Real recf[4476]; ///<Array of coefficients.
};

class DROCK2: public ROCK2
{
public:
    DROCK2(Ode* ode, bool verb=true, bool dtadap=true, 
          Real atol=1e-2, Real rtol=1e-2, bool scalartol=true,
          int output_freq=0);
    virtual ~DROCK2();
    
protected:

    virtual void step(const Real t, const Real& h);

    virtual void update_n_stages_and_h(Real& h);    
 
protected:
    static Real recf2[184]; ///<Array of coefficients.
    static Real recalph[46];///<Array of coefficients.
};


#endif	/* ODESTABILIZEDINTEGRATORS_H */