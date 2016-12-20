#include "TraditionalSdeRungeKuttaIntegrators.h"

MilsteinTalay::MilsteinTalay(Sde* sde, bool verb, bool dtadap,
                  Real atol, Real rtol, 
                  bool scalartol, int output_freq)
:OdeRungeKuttaIntegrator(sde,verb,dtadap,atol,rtol,scalartol,output_freq),
 SdeRungeKuttaIntegrator(sde,verb,dtadap,atol,rtol,scalartol,output_freq)
{
    integr5 = new Vector(sde->system_size());
    needDoubleInt = true;
    solver = "MilsteinTalay";
}

MilsteinTalay::~MilsteinTalay()
{
    delete integr5;
}

void MilsteinTalay::update_n_stages_and_h(Real& h)
{
    if(h>0.9*2./eigmax)
    {
        h = 0.9*2/eigmax;
        last=false;
    }
}

void MilsteinTalay::step_diagonal_noise(const Real t, const Real& h)
{
    //First convert ode to sde, so that we can call g and the brownian motion
    Sde* sde = dynamic_cast<Sde*>(ode);
    
    Vector*& K1= integr[0];
    Vector*& K2= integr[1];
    Vector*& K3= integr[2];
    Vector*& A = integr[3];
    Vector*& B = integr[4];
    Vector*& C = integr5;
    Vector* swap_ptr=0;
    
    sde->sample(t,h);
    Vector& Ir = sde->getIr(); // Brownian increments
    Vector& Irr = sde->getIrr(); // double stochastic integrals
    Vector& Chir = sde->getChir();
    
    sde->g(t,*yn,Gv); 
    
    //first step
    sde->f(t,*yn,*A);       // A = f(X_0)
    *K1=*yn+h*(*A);
    
    *ynpu = *yn;
    *ynpu += 0.5*h*(*A);

    
    *K2 = *K1 + Gv.cwiseProduct(Ir);
    sde->f(t+h,*K2,*A);
    *ynpu += 0.5*h*(*A);
    
    *K1 *= 0.5;
    *K1 += 0.5*(*yn);
    
    *K2 = *K1+(sqrt(h/2.))*(Gv.cwiseProduct(Chir));
    sde->g(t+h,*K2,*K3);
    *ynpu += 0.5*(K3->cwiseProduct(Ir));
    
    *K2 = *K1-(sqrt(h/2.))*(Gv.cwiseProduct(Chir));
    sde->g(t+h,*K2,*K3);
    *ynpu += 0.5*(K3->cwiseProduct(Ir));
    
    *K2 = Gv.cwiseProduct(Irr);
    *K1 = *yn+(*K2);
    sde->g(t,*K1,*K3);
    *ynpu += 0.5*(*K3);
    
    *K1 = *yn-(*K2);
    sde->g(t,*K1,*K3);
    *ynpu -= 0.5*(*K3);
    
    //update the number of right hand side evaluations
    n_f_eval=n_f_eval + 2;
    n_g_eval= n_g_eval + 5; 
}

void MilsteinTalay::step_general_noise(const Real t, const Real& h)
{
    
}