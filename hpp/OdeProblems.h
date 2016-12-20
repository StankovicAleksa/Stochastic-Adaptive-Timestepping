#ifndef ODEPROBLEMS_H
#define	ODEPROBLEMS_H

#include "Ode.h"

class OneDimensionalTest1: public virtual Ode
{
public:
    OneDimensionalTest1(string output_file);
    virtual ~OneDimensionalTest1(){};
 
    void init_solution();
    void f(Real t, Vector& x, Vector& fx);
    
    void rho(Real& eigmax);
};

class TwoDimensionalBrusselator: public virtual Ode
{
public:
    TwoDimensionalBrusselator(string output_file);
    virtual ~TwoDimensionalBrusselator(){};
 
    void init_solution();
    void f(Real t, Vector& x, Vector& fx);
    
    void rho(Real& eigmax);
};

class NeuronCable: public virtual Ode
{
public:
    NeuronCable(string output_file);
    virtual ~NeuronCable(){};
 
    void init_solution();
    void f(Real t, Vector& x, Vector& fx);
    
    void rho(Real& eigmax);
    
protected:
    Real nu;
    Real beta;
};

class PDEBrusselator: public virtual Ode
{
public:
    PDEBrusselator(string output_file);
    virtual ~PDEBrusselator(){};
 
    void init_solution();
    void f(Real t, Vector& x, Vector& fx);
    
    void rho(Real& eigmax);
protected:
    Real A;
    Real B;
    Real alpha;
};


class Krogh10: public virtual Ode
{
public:
    Krogh10(string output_file);
    virtual ~Krogh10(){};
 
    void init_solution();
    void f(Real t, Vector& x, Vector& fx);
    
    void rho(Real& eigmax);
protected:
    Real mu;
    Real mus;
};


#endif	/* ODEPROBLEMS_H */