#include "OdeProblems.h"

// One dimensional test 1
OneDimensionalTest1::OneDimensionalTest1(string output_file)
:Ode(string("Tests/OdeOneDimensionalTest1/")+output_file)
{
    t0=0.0;
    tend=1.0;
    t=t0;
    
    neqn=1;
    cte_rho = false;
    know_rho = true;
    
    init_solution();
}

void OneDimensionalTest1::init_solution()
{
    y = new Vector(neqn);
    
    (*y)(0)=0.;
}

void OneDimensionalTest1::f(Real t, Vector& x, Vector& fx)
{   
    fx(0) = 0.25*x(0)+0.5*sqrt(x(0)*x(0)+1);
}

void OneDimensionalTest1::rho(Real& eigmax)
{
    eigmax = abs(0.25+0.5*(*y)(0)/sqrt((*y)(0)*(*y)(0)+1.))+1;
}

// This is the brusseltaor in 2d

TwoDimensionalBrusselator::TwoDimensionalBrusselator(string output_file)
:Ode(string("Tests/TwoDimensionalBrusselator/")+output_file)
{
    t0=0.0;
    tend=1.0;
    t=t0;
    
    neqn=2;
    cte_rho = false;
    know_rho = true;
    
    init_solution();
}

void TwoDimensionalBrusselator::init_solution()
{
    y = new Vector(neqn);
    
    (*y)(0)=1.5;
    (*y)(1)=3.;
}

void TwoDimensionalBrusselator::f(Real t, Vector& x, Vector& fx)
{   
    fx(0) = 1.+x(0)*x(0)*x(1)-4.*x(0);
    fx(1) = 3.*x(0)-x(0)*x(0)*x(1);
}

void TwoDimensionalBrusselator::rho(Real& eigmax)
{
    Real y1 = (*y)(0);
    Real y2 = (*y)(1);
    Real a = 2.*y1*y2-4.;
    Real b = y1*y1;
    Real c = 3.-2.*y1*y2;
    Real d = -y1*y1;
    b = (a+d)/2.;
    d = a*a+d*d-2.*a*d+4.*b*c;
    d = d/4.;
    if(d>=0)
    {
        d=sqrt(d);
        eigmax = max(abs(b+d),abs(b-d));
    }
    else
    {
        d=sqrt(-d);
        eigmax = sqrt(b*b+d*d);
    }
}


// Neuron cable equation

NeuronCable::NeuronCable(string output_file)
:Ode(string("Tests/NeuronCable/")+output_file)
{
    t0=0.0;
    tend=1.0;
    t=t0;
    nu = 0.01;
    beta= 1.0;
    
    neqn=128;
    cte_rho = true;
    know_rho = true;
    
    init_solution();
}

void NeuronCable::init_solution()
{
    y = new Vector(neqn);
    
    const Real pi=4.*atan(1.);
    for (int j=0;j<neqn;j++)
    {
        Real x=((Real)j)/(neqn-1.);
        (*y)(j)=-70.+20.*cos(15.*pi*x)*(1.-x);
    }
}

void NeuronCable::f(Real t, Vector& x, Vector& fx)
{   

    // Computing diffusion with Neumann bnd conditions
    fx(0)=nu*2.*(x(1)-x(0))*(neqn-1.)*(neqn-1.)-beta*x(0);
    fx(neqn-1)=nu*2.*(x(neqn-2)-x(neqn-1))*(neqn-1.)*(neqn-1.)-beta*x(neqn-1);
    for (int i=1;i<neqn-1;i++)
    {
        fx(i)=nu*(x(i-1)-2.*x(i)+x(i+1))*(neqn-1.)*(neqn-1.)-beta*x(i);
        if(abs(i/(neqn-1.)-0.5)<0.1)
            fx(i) += 5.*exp(1.-1e-2/(1e-2-(i/(neqn-1.)-0-5)*(i/(neqn-1.)-0.5)));
    }
}

void NeuronCable::rho(Real& eigmax)
{
    eigmax = nu*4.*(neqn-1)*(neqn-1)+beta;
}


PDEBrusselator::PDEBrusselator(string output_file)
:Ode(string("Tests/PDEBrusselator/")+output_file)
{
    t0=0.0;
    tend=10.0;
    t=t0;
    
    neqn=2*64;
    cte_rho = false;
    know_rho = false;
    
    A = 1.;
    B = 3.;
    alpha = 1./50.;
    
    init_solution();
}

void PDEBrusselator::init_solution()
{
    y = new Vector(neqn);
    
    int n=neqn/2;
    const Real pi=4.*atan(1.);
    
    for(int i=0;i<n;i++)
    {
        Real x = ((Real)i+1.)/(n+1.);
        (*y)(i) = 1.+sin(2.*pi*x);
        (*y)(i+n) = 3.;
    }
}

void PDEBrusselator::f(Real t, Vector& x, Vector& fx)
{   
    Real ul = 1.;
    Real ur = 1.;
    Real vl = 3.;
    Real vr = 3.;
    int n=neqn/2;
    int npusq = (n+1)*(n+1);
    
    fx(0) = A+ x(0)*x(0)*x(n)-(B+1.)*x(0)+alpha*npusq*(ul-2.*x(0)+x(1));
    fx(n-1) = A+ x(n-1)*x(n-1)*x(n-1+n)-(B+1.)*x(n-1)+alpha*npusq*(x(n-2)-2.*x(n-1)+ur);
    for(int i=1;i<n-1;i++)
        fx(i) = A+ x(i)*x(i)*x(i+n)-(B+1.)*x(i)+alpha*npusq*(x(i-1)-2.*x(i)+x(i+1));
    
    fx(n) = B*x(0) - x(0)*x(0)*x(n)+alpha*npusq*(vl-2.*x(n)+x(n+1));
    fx(2*n-1) = B*x(n-1) - x(n-1)*x(n-1)*x(2*n-1)+alpha*npusq*(x(2*n-2)-2.*x(2*n-1)+vr);
    for(int i=n;i<2*n-1;i++)
        fx(i) = B*x(i-n) - x(i-n)*x(i-n)*x(i)+alpha*npusq*(x(i-1)-2.*x(i)+x(i+1));
}

void PDEBrusselator::rho(Real& eigmax)
{
    eigmax = -1000.;
    cout<<"WARNING: Spectral radius has to be estimated via power method."<<endl;
    cout<<"Choose -intrho 1"<<endl;
}


// Test problem 10 from Krogh
Krogh10::Krogh10(string output_file)
:Ode(string("Tests/Krogh10/")+output_file)
{
    t0= 0.0;
    tend = 6.19216933131963970674;
    t=t0;
    mu = 1./82.45;
    mus = 1.-mu;
    
    neqn=4;
    cte_rho = false;
    know_rho = true;
    
    init_solution();
}

void Krogh10::init_solution()
{
    y = new Vector(neqn);
    
    // These initial conditions give a periodic solution, so they are also
    // exact solution at tend
    (*y)(0)=1.2;
    (*y)(1)=0.0;
    (*y)(2)=0.0;
    (*y)(3)= -1.0493575098031990726;
}

void Krogh10::f(Real t, Vector& x, Vector& fx)
{   
    Real y1 = x(0);
    Real y2 = x(1);
    Real y1p = x(2);
    Real y2p = x(3);
    Real r1 = sqrt(((y1+mu)*(y1+mu)+y2*y2));
    Real r2 = sqrt(((y1-mus)*(y1-mus)+y2*y2));
    r1 = r1*r1*r1;
    r2 = r2*r2*r2;
    
    fx(0) = y1p;
    fx(1) = y2p;
    fx(2) = 2.*y2p+y1-mus*(y1+mu)/r1-mu*(y1-mus)/r2;
    fx(3) = -2.*y1p+y2-mus*y2/r1-mu*y2/r2;
}

void Krogh10::rho(Real& eigmax)
{
    eigmax = abs(0.25+0.5*(*y)(0)/sqrt((*y)(0)*(*y)(0)+1.));
}