#ifndef BROWNIANMOTION_H
#define	BROWNIANMOTION_H

#include <iostream>
#include <random>
#include <Eigen/Dense>

typedef double Real;
typedef Eigen::VectorXd Vector;
using namespace std;

// This is a very simple class for a Brownian motion. Not possible to 
// sample "past" values yet. Neither double integrals I_rq but only I_rr

class BrownianMotion
{
public:
    BrownianMotion(int size, bool doubleintt, bool commutativee, bool diagonall);
    virtual ~BrownianMotion();
    
    virtual void sample(Real t, Real h) =0;
    Vector& getIr();
    Vector& getChir();
    Vector& getIrr();
    vector<Vector>& getIpq();
    Vector& getOnes();
    
    void resize(int size);
    void clear();
    const bool isContinuous() const;
    int size();
    
protected:
    int m;
    Vector Ir;
    Vector Irr;
    Vector Chir;
    vector<Vector> Ipq;
    Vector Ones;
    
    const bool doubleint;
    const bool commutative;
    const bool diagonal;
    bool continuous;
    
    default_random_engine gen;
};

class ContinuousBrownianMotion: public BrownianMotion
{
public:
    ContinuousBrownianMotion(int size, bool doubleintt, bool commutativee, bool diagonall);
    virtual ~ContinuousBrownianMotion();
    
    void sample(Real t, Real h);
    
protected:
    normal_distribution<Real> normal;
};

class DiscreteBrownianMotion: public BrownianMotion
{
public:
    DiscreteBrownianMotion(int size, bool doubleintt, bool commutativee, bool diagonall);
    virtual ~DiscreteBrownianMotion();
    
    void sample(Real t, Real h);
    
protected:
    uniform_int_distribution<int> xid;
    uniform_int_distribution<int> chid;
    Real sqr3;
};

#endif	/* BROWNIANMOTION_H */

