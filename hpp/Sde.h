#ifndef SDE_H
#define	SDE_H

#include "Ode.h"
#include "BrownianMotion.h"

typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> Matrix;

class Sde: public virtual Ode
{
public:    
    Sde(string solution_file);
    virtual ~Sde();
    
    virtual void g(Real t, Vector& x, Vector& G);//only for diagonal noise
    virtual void g(Real t, Vector& x, Matrix& G);//only for general noise
    virtual void g(Real t, Vector& x, Vector& G, int r);//only for general noise
    
    void sample(Real t, Real h);
    Vector& getIr();
    Vector& getChir();
    Vector& getIrr();
    vector<Vector>& getIpq();
    Vector& getOnes();
    
    virtual Real phi() =0;
    virtual Real Exact_phi() =0;
    
    int brownian_size();
    bool is_commutative() const;
    bool isDiagonal() const;
    void reinit();
    void setW(BrownianMotion* W);
        
protected:
    BrownianMotion* W;
    int Wsize;
    bool diagonal;
    bool commutative;
};

#endif	/* SDE_H */

// For diagonal noise
// g(...,Vector) returns a vector where each element 
// is g^i(x^i) for i=1,...,neqn=m

//For general noise
// g(...,Vector,int r) returns the rth diffusion g^r 
// g(...,Matrix) return a matrix with m columns, where each columns is g^r(x))
