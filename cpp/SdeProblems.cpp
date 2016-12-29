#include <cmath>
#include "SdeProblems.h"


// Sde One dimensional test 1 from SROCK2 paper
SdeOneDimensionalTest1::SdeOneDimensionalTest1(string output_file)
:Ode(string("Tests/SdeOneDimensionalTest1/")+output_file),
 Sde(string("Tests/SdeOneDimensionalTest1/")+output_file),
 OneDimensionalTest1(string("Tests/SdeOneDimensionalTest1/")+output_file)       
{
    neqn = 1;
    Wsize = 1;
    
    diagonal = true;
    commutative = true;
    
    cte_rho=false;
    know_rho = true;
}

Real SdeOneDimensionalTest1::phi()
{
    Real x = asinh((*y)(0));
    return x*x;
}

Real SdeOneDimensionalTest1::Exact_phi()
{
    return tend*tend/4.+tend/2.;
}

void SdeOneDimensionalTest1::g(Real t, Vector& x, Vector& G)
{
    G(0) = sqrt((x(0)*x(0)+1.)/2.);
}

void SdeOneDimensionalTest1::g(Real t, Vector& x, Vector& G, int r)
{
    G(0) = sqrt((x(0)*x(0)+1.)/2.);
}

void SdeOneDimensionalTest1::g(Real t, Vector& x, Matrix& G)
{
    G(0,0) = sqrt((x(0)*x(0)+1.)/2.);
}


// Sde One dimensional test 2 from SROCK2 paper
SdeOneDimensionalTest2::SdeOneDimensionalTest2(string output_file)
:Ode(string("Tests/SdeOneDimensionalTest2/")+output_file),
 Sde(string("Tests/SdeOneDimensionalTest2/")+output_file)
{
    neqn = 1;
    Wsize = 10;
    
    diagonal = false;
    commutative = false;
    
    cte_rho = true;
    know_rho = true;
    
    t0=0.0;
    tend=1.0;
    t=t0;
 
    a.resize(Wsize);
    b.resize(Wsize);
    a<< 10., 15., 20., 25., 40., 25., 20., 15., 20., 25.;
    b<< 2., 4., 5., 10., 20., 2., 4., 5., 10., 20.;
    for(int i=0;i<Wsize;i++)
    {
        a(i)=1./a(i);
        b(i)=1./b(i);
    }
    
    init_solution();
}

void SdeOneDimensionalTest2::init_solution()
{
    y = new Vector(neqn);
    
    (*y)(0)=1.;
}

void SdeOneDimensionalTest2::f(Real t, Vector& x, Vector& fx)
{   
    fx(0) = x(0);
}

void SdeOneDimensionalTest2::g(Real t, Vector& x, Vector& G, int r)
{   
    G(0) = a(r)*sqrt(x(0)+b(r));
}

void SdeOneDimensionalTest2::g(Real t, Vector& x, Matrix& G)
{
    for(int r=0;r<Wsize;r++)
        G(0,r) = a(r)*sqrt(x(0)+b(r));
}

Real SdeOneDimensionalTest2::phi()
{
    return (*y)(0)*(*y)(0);
}

Real SdeOneDimensionalTest2::Exact_phi()
{
    return (-68013.-458120.*exp(tend)+14926133.*exp(2.*tend))/14400000.;
}

void SdeOneDimensionalTest2::rho(Real& eigmax)
{
    eigmax = 1.;
}

StochasticTwoBodyProblem::StochasticTwoBodyProblem(string output_file)
:Ode(string("Tests/StochasticTwoBodyProblem/")+output_file),
 Sde(string("Tests/StochasticTwoBodyProblem/")+output_file)
{
    neqn = 4;
    Wsize = 4;
    
    diagonal = true;
    commutative = true;
    
    cte_rho = false;
    know_rho = false;
    
    t0=0.0;
    tend=58.0;
    t=t0;
    
    sigmar = 0.0121;
    sigmaphi = 2.2*1e-4;
    G = 0.0252505160442859;
    m = 0.1*0.1/0.2;
    Gm = G/0.2;
    
    init_solution();
}

void StochasticTwoBodyProblem::init_solution()
{
    y = new Vector(neqn);
    
    (*y)(0)=1.;
    (*y)(1)=1.;
    (*y)(2)=0.01;
    (*y)(3)=1.1;
            
}

void StochasticTwoBodyProblem::f(Real t, Vector& x, Vector& fx)
{   
    // x = (r,phi,v,w)
    fx(0)=x(2);
    fx(1)=x(3);
    fx(2)=x(0)*x(3)*x(3)-Gm/(x(0)*x(0));
    fx(3)=-2.*x(2)*x(3)/x(0);
}

void StochasticTwoBodyProblem::g(Real t, Vector& x, Vector& G)
{
    G(0) = 0.;
    G(1) = 0.;
    G(2) = x(0)*sigmar;
    G(3) = sigmaphi/x(0);
}

Real StochasticTwoBodyProblem::phi()
{
    return m*(*y)(0)*(*y)(0)*(*y)(3);
}

Real StochasticTwoBodyProblem::Exact_phi()
{
    return m*1.1;
}


SdeLinear::SdeLinear(string output_file)
:Ode(string("Tests/SdeLinear/")+output_file),
 Sde(string("Tests/SdeLinear/")+output_file)
{
    neqn = 1;
    Wsize = 1;
    
    diagonal = true;
    commutative = true;
    
    cte_rho = true;
    know_rho = true;
    
    t0=0.0;
    tend=1.0;
    t=t0;
    
    lambda = -1.;
    sigma = 0.1;
    
    init_solution();
}

void SdeLinear::init_solution()
{
    y = new Vector(neqn);
    
    (*y)(0)=1.;
}

void SdeLinear::f(Real t, Vector& x, Vector& fx)
{   
    fx(0) = lambda*x(0);
}

void SdeLinear::g(Real t, Vector& x, Vector& G)
{
    G(0) = sigma*x(0);
}

Real SdeLinear::phi()
{
    return (*y)(0);
}

Real SdeLinear::Exact_phi()
{
    return exp(lambda*tend);
}

void SdeLinear::rho(Real& eigmax)
{
    eigmax = abs(lambda);
}
