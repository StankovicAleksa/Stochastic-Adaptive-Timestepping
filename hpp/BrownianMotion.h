#ifndef BROWNIANMOTION_H
#define	BROWNIANMOTION_H

#include <iostream>
#include <random>
#include <Eigen/Dense>
#include <set>

typedef double Real;
typedef Eigen::VectorXd Vector;
using namespace std;


// Required for timestepping 

// Sample of Brownian motion at time t with value val
struct sp_sample{
	Real val;
	Real t;
	sp_sample(Real _t, Real _val){ t=_t; val=_val;}
	sp_sample(){ t=0; val=0;} // default constructor
};
// operator < is needed for ordering by sampling time
bool operator < (sp_sample s1, sp_sample s2);


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
   	virtual void eraseHistory() ; 
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

class ContinuousAdaptedBrownianMotion: public BrownianMotion
{
public:
    ContinuousAdaptedBrownianMotion(int size, bool doubleintt, bool commutativee, bool diagonall);
    virtual ~ContinuousAdaptedBrownianMotion();
    
    void sample(Real t, Real h);
		void eraseHistory();    
protected:
    std::set<sp_sample> *samples;
    normal_distribution<Real> normal;
};

class DiscreteAdaptedBrownianMotion: public BrownianMotion
{
public:
    DiscreteAdaptedBrownianMotion(int size, bool doubleintt, bool commutativee, bool diagonall);
    virtual ~DiscreteAdaptedBrownianMotion();
    
    void sample(Real t, Real h);
		void eraseHistory();    
    
protected:
    std::set<sp_sample> *samples;
    uniform_int_distribution<int> xid;
    uniform_int_distribution<int> chid;
    Real sqr3;
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

