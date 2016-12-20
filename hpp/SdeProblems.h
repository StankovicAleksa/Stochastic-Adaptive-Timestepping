#ifndef SDEPROBLEMS_H
#define	SDEPROBLEMS_H

#include "Sde.h"
#include "OdeProblems.h"

// Problem 1 in SROCK2 paper
class SdeOneDimensionalTest1: public Sde, public OneDimensionalTest1
{
public:
    SdeOneDimensionalTest1(string output_file);
    virtual ~SdeOneDimensionalTest1(){};
 
    void g(Real t, Vector& x, Vector& G);
    void g(Real t, Vector& x, Vector& G, int r);
    void g(Real t, Vector& x, Matrix& G);
        
    Real phi();
    Real Exact_phi();
};

// Problem 2 in SROCK2 paper
class SdeOneDimensionalTest2: public Sde
{
public:
    SdeOneDimensionalTest2(string output_file);
    virtual ~SdeOneDimensionalTest2(){};
 
    void init_solution();
    
    void f(Real t, Vector& x, Vector& fx);
    void g(Real t, Vector& x, Vector& G, int r);
    void g(Real t, Vector& x, Matrix& G);
        
    Real phi();
    Real Exact_phi();
    
    void rho(Real& eigmax);
    
protected:
    Vector a;
    Vector b;
};


class StochasticTwoBodyProblem: public Sde
{
public:
    StochasticTwoBodyProblem(string output_file);
    virtual ~StochasticTwoBodyProblem(){};
 
    void init_solution();
    
    void f(Real t, Vector& x, Vector& fx);
    void g(Real t, Vector& x, Vector& G);
        
    Real phi();
    Real Exact_phi();
        
protected:
    Real sigmar;
    Real sigmaphi;
    Real G;
    Real m;
    Real Gm;
};


// Linear 1D  test problem
class SdeLinear: public Sde
{
public:
    SdeLinear(string output_file);
    virtual ~SdeLinear(){};
 
    void init_solution();
    
    void f(Real t, Vector& x, Vector& fx);
    void g(Real t, Vector& x, Vector& G);
        
    Real phi();
    Real Exact_phi();
    
    void rho(Real& eigmax);
    
protected:
    Real lambda;
    Real sigma;

};

#endif	/* SDEPROBLEMS_H */

