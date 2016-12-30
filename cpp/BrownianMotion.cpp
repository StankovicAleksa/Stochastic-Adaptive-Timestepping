#include <chrono>
#include "BrownianMotion.h"
#include <omp.h>
#include <set>
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
void BrownianMotion::eraseHistory(){
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


// Adapted motions
ContinuousAdaptedBrownianMotion::ContinuousAdaptedBrownianMotion(int size, bool doubleintt, 
                                                   bool commutativee, 
                                                   bool diagonall)
:BrownianMotion(size,doubleintt,commutativee,diagonall), normal(0.0,1.0)
{
    continuous=true;
    resize(size);
		samples=new	std::set<sp_sample>[size];
		//initialize Brownian motion at t=0
		for (int i=0;i<size;i++)
		samples[i].insert(sp_sample(0,0));
}
ContinuousAdaptedBrownianMotion::~ContinuousAdaptedBrownianMotion()
{   
	delete[] samples;
}
void ContinuousAdaptedBrownianMotion::eraseHistory(){
		for (int i=0;i<m;i++){
			samples[i].clear();
		samples[i].insert(sp_sample(0,0));
		}
}

void ContinuousAdaptedBrownianMotion::sample(Real t, Real h)
{

	Real exp; // expectation
	Real std_dev; // standard deviation
	std::set<sp_sample>::iterator sp_is0,sp_is1; 
	Real s_val;
	Real eps=1e-11;
	for ( int i=0;i<Ir.size();i++){
		s_val=upper_bound(samples[i].begin(),samples[i].end(),sp_sample(t-eps,0))->val;	

		// find samples for t0< t < t1
		sp_is1=std::upper_bound(samples[i].begin(),samples[i].end(),sp_sample(t+h-eps,0));
		while ( sp_is1!= samples[i].end() && sp_is1->t < t+ h ){
			sp_is1++;
		}
		sp_is0=prev(sp_is1); 
	
		exp=sp_is0->val-s_val;
		if ( sp_is1 == samples[i].end() ){
			std_dev=sqrt( h-(sp_is0->t - t ));
		} else{
			Real t_t0=t+h-sp_is0->t;
			Real t1_t=sp_is1->t-t-h;
			Real t1_t0=sp_is1->t-sp_is0->t;
			exp+= t_t0*(sp_is1->val-sp_is0->val)/(t1_t0);
			std_dev=sqrt(t_t0*t1_t/t1_t0);
		}
		Ir(i)=exp+std_dev*normal(gen);
		// Biased sampling
		//Ir(i)=sqrt(h)*normal(gen);
    
		Irr(i) = 0.5*(Ir(i)*Ir(i)-h);
		std::cout.flush();
		samples[i].insert(sp_sample(t+h,Ir(i)+s_val));
	}	 
    
    if(doubleint && !diagonal)
        cout<<"implement double integrals for continuous brownian motion"<<endl;
}

DiscreteAdaptedBrownianMotion::DiscreteAdaptedBrownianMotion(int size, bool doubleintt, 
                                               bool commutativee, 
                                               bool diagonall)
:BrownianMotion(size,doubleintt,commutativee,diagonall), xid(1,6), chid(0,1)
{
    sqr3 = sqrt(3.);
    continuous=false;
    resize(size);
		samples=new	std::set<sp_sample>[size];
		//initialize Brownian motion at t=0
		for (int i=0;i<size;i++)
		samples[i].insert(sp_sample(0,0));
}

DiscreteAdaptedBrownianMotion::~DiscreteAdaptedBrownianMotion()
{
}
void DiscreteAdaptedBrownianMotion::eraseHistory(){
		for (int i=0;i<m;i++){
			samples[i].clear();
		samples[i].insert(sp_sample(0,0));
		}
}

void DiscreteAdaptedBrownianMotion::sample(Real t, Real h)
{
    Real sqrh = sqrt(h);
    Real xi, chi;
		Real exp; // expectation
		Real std_dev; // standard deviation
		std::set<sp_sample>::iterator sp_is0,sp_is1; 
		Real s_val;
		Real eps=1e-11;
		for ( int i=0;i<Ir.size();i++){
			s_val=upper_bound(samples[i].begin(),samples[i].end(),sp_sample(t-eps,0))->val;	

			// find samples for t0< t < t1
			sp_is1=std::upper_bound(samples[i].begin(),samples[i].end(),sp_sample(t+h-eps,0));
			while ( sp_is1!= samples[i].end() && sp_is1->t < t+ h ){
				sp_is1++;
			}
			sp_is0=prev(sp_is1); 
		
			exp=sp_is0->val-s_val;
			if ( sp_is1 == samples[i].end() ){
				std_dev=sqrt( h-(sp_is0->t - t ));
			} else{
				Real t_t0=t+h-sp_is0->t;
				Real t1_t=sp_is1->t-t-h;
				Real t1_t0=sp_is1->t-sp_is0->t;
				exp+= t_t0*(sp_is1->val-sp_is0->val)/(t1_t0);
				std_dev=sqrt(t_t0*t1_t/t1_t0);
			}
      xi = xid(gen);
      xi = (xi<=4.) ? 0. : ((xi==5.) ? sqr3 : -sqr3) ;
			Ir(i)=exp+std_dev*xi;
			// Biased sampling
			//Ir(i)=sqrt(h)*xi;
			
			//Irr(i) = 0.5*(Ir(i)*Ir(i)-std_dev*std_dev);
			Irr(i) = 0.5*(Ir(i)*Ir(i)-h);
			std::cout.flush();
			samples[i].insert(sp_sample(t+h,Ir(i)+s_val));
		}	 
    
    for(int i=0;i<Ir.size();i++)
    {
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
bool operator < (sp_sample s1, sp_sample s2){
	return s1.t< s2.t;
}
