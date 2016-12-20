#include <chrono>
#include "BrownianMotion.h"
#include <omp.h>

BrownianMotion::BrownianMotion(int size, bool doubleintt, bool commutativee, 
                               bool diagonall)
:doubleint(doubleintt), commutative(commutativee), diagonal(diagonall)
{
    unsigned seed1 = chrono::system_clock::now().time_since_epoch().count();
//    cout<<"USING FIXED SEED !!"<<endl;
//    gen.seed(1);
    gen.seed(seed1+omp_get_thread_num());
}

BrownianMotion::~BrownianMotion()
{
    
}

Vector& BrownianMotion::getIr()
{
    return Ir;
}

Vector& BrownianMotion::getChir()
{
    return Chir;
}

Vector& BrownianMotion::getIrr()
{
    return Irr;
}

vector<Vector>& BrownianMotion::getIpq()
{
    return Ipq;
}

int BrownianMotion::size()
{
    return m;
}

void BrownianMotion::resize(int size)
{
    m = size;
    Ir.resize(m);
    Ones.resize(m);
    if(!continuous)
        Chir.resize(m);
    for(int i=0;i<m;i++)
        Ones(i) = 1.;
    
    if(doubleint)
        Irr.resize(m);
    if(doubleint && !diagonal)
        Ipq = vector<Vector>(m,Vector(m));
}

void BrownianMotion::clear()
{
    
}

const bool BrownianMotion::isContinuous() const
{
    return continuous;
}

Vector& BrownianMotion::getOnes()
{
    return Ones;
}

ContinuousBrownianMotion::ContinuousBrownianMotion(int size, bool doubleintt, 
                                                   bool commutativee, 
                                                   bool diagonall)
:BrownianMotion(size,doubleintt,commutativee,diagonall), normal(0.0,1.0)
{
    continuous=true;
    resize(size);
}

ContinuousBrownianMotion::~ContinuousBrownianMotion()
{   
}

void ContinuousBrownianMotion::sample(Real t, Real h)
{
    Real sqrh = sqrt(h);
    
    for(int i=0;i<Ir.size();i++)
    {
        Ir(i) = sqrh*normal(gen);
        Irr(i) = 0.5*(Ir(i)*Ir(i)-h);
    }
    
    if(doubleint && !diagonal)
        cout<<"implement double integrals for continuous brownian motion"<<endl;
}

DiscreteBrownianMotion::DiscreteBrownianMotion(int size, bool doubleintt, 
                                               bool commutativee, 
                                               bool diagonall)
:BrownianMotion(size,doubleintt,commutativee,diagonall), xid(1,6), chid(0,1)
{
    sqr3 = sqrt(3.);
    continuous=false;
    resize(size);
}

DiscreteBrownianMotion::~DiscreteBrownianMotion()
{
}

void DiscreteBrownianMotion::sample(Real t, Real h)
{
    Real sqrh = sqrt(h);
    Real xi, chi;
    
    for(int i=0;i<Ir.size();i++)
    {
        xi = xid(gen);
        xi = (xi<=4.) ? 0. : ((xi==5.) ? sqr3 : -sqr3) ;

        Ir(i) = sqrh*xi;
        Irr(i) = 0.5*(Ir(i)*Ir(i)-h);
        
        Chir[i] = chid(gen);
        if(Chir[i]==0.) Chir[i]=-1.;
    }
    
    if(doubleint && !diagonal)
    {
        if(commutative)
        {
            cout<<"Implement commutative sampling"<<endl;
        }
        else
        {
            for(int p=0;p<m;p++) // p=0,...,m-1
            for(int q=0;q<p;q++) // q=0,...,p-1, q<p
            {   
                Ipq[p](q) = 0.5*(Ir[p]*Ir[q]-h*Chir[p]);
                Ipq[q](p) = 0.5*(Ir[p]*Ir[q]+h*Chir[p]);
            }   
            for(int r=0;r<Ir.size();r++)
                Ipq[r](r)=Irr(r);
        }
    }
}