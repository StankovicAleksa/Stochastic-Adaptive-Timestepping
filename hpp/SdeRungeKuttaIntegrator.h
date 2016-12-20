#ifndef SDERUNGEKUTTAINTEGRATOR_H
#define	SDERUNGEKUTTAINTEGRATOR_H

#include "OdeRungeKuttaIntegrator.h"
#include "Sde.h"

typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> Matrix;

class SdeRungeKuttaIntegrator: public virtual OdeRungeKuttaIntegrator
{
public:
    SdeRungeKuttaIntegrator(Sde* sde, bool verb=false, 
                            bool dtadap=true, Real atol=1e-2, Real rtol=1e-2, 
                            bool scalartol=true, int output_freq=0);
    virtual ~SdeRungeKuttaIntegrator();
    
    bool needDoubleIntegral();
    
    virtual void reinit();
    void print_info();                     
    void print_statistics();
    int get_n_g_eval();        
    
protected:
    virtual void step(const Real t, const Real& h);
    virtual void step_general_noise(const Real t, const Real& h) =0;
    virtual void step_diagonal_noise(const Real t, const Real& h) =0;
    
protected:
    Matrix G;
    Vector Gv;
    int n_g_eval;
    Real errD;
    Real errS;
    bool needDoubleInt;
};

#endif	/* SDERUNGEKUTTAINTEGRATOR_H */

