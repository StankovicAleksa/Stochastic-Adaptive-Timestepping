#include <limits>
#include <cstdlib>
#include <string>
#include <fstream>

#include "OdeRungeKuttaIntegrator.h"

OdeRungeKuttaIntegrator::OdeRungeKuttaIntegrator(Ode* ode, bool verb, bool dtadap, 
                                    Real atol, Real rtol, bool scalartol, 
                                    int output_freq)
:TimeIntegrator(ode),yn(ode->get_y()),out_file(ode->get_solution_file_name()),
 internal_rho(!ode->do_know_rho())
{
    //auxiliary variables initialization
//    yn = ode->get_y();
    ynpu = new Vector(ode->system_size());
    eigenvector = new Vector(ode->system_size());
    for(int i=0;i<5;i++)
        integr[i]=new Vector(ode->system_size());
    
    uround=numeric_limits<Real>::epsilon();
    
    nu = 1.0;
    
    reinit();
    
    verbose=verb;
    dt_adaptivity = dtadap;
    scalar_tol = scalartol;
    r_tol = rtol;
    a_tol = atol;
    output_frequency=output_freq;
    Astable=false; 
    
    cout<<fixed<<setfill(' ');
}

 OdeRungeKuttaIntegrator::~OdeRungeKuttaIntegrator()
{
    delete ynpu;
    delete eigenvector;
    for(int i=0;i<5;i++)
        delete integr[i];
}

void OdeRungeKuttaIntegrator::reinit()
{
    errp = 0.0;
    err = 0.0;
    firsterrorestimation=true;
    
    nrho = 0;
    facmax=5.0;    
    nrej=0;
    rho_outdated=true;
    
    s=0;
    sold=0;
    
    max_rho=numeric_limits<int>::min();
    min_rho=numeric_limits<int>::max();
    max_s=0;
    n_f_eval_rho=0;
    n_f_eval=0;
    n_steps = 0;
    acc_steps = 0;
    rej_steps = 0;
    dt_max = 0.;
    dt_min = numeric_limits<Real>::max();
}

bool  OdeRungeKuttaIntegrator::integrate(Real& h)
{

    /**
     * Ode is a class describing the Cauchy problem. It has a starting and end 
     * time step, time step, an initial value and a rhs function returning the 
     * right hand side.
     * This function integrates ode with a starting time step h. If one_step is
     * true just one step is executed. Recalling advance with the same arguments
     * does the next step using a new adapted h (id dt_adaptivity is true). 
     * If one_step is false integration is performed until end.
     * The returning value of idid is 1 if we reached the end, 2 if a step
     * has been executed correctly without reaching the end.
     */
    
    Real told;
    Vector *swap_ptr;
    
    Real t = ode->get_t(); 
    Real tend = ode->get_tend();
    
//    Vector*& ;
//    yn = ode->get_y();  //reference to current solution
    
        
// *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** 
//             Initializations
// *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** 
 
    if(last)//we are reusing the integrator
        firsterrorestimation=true;
    told=t-1.;    //different from t
    last=false;   //a priori is not the last step
    reject=false; //if we enter here then the previous step has been accepted    
    
    elapsed_time = get_cpu_time();
    
    for(;;)
    {
        //if we are very close to end or t+h>tend then adjust h
        if(1.1*abs(h)>=abs(tend-t)) //h can be negative (when t>tend)
        {
            h=tend-t;
            last=true; //it will we be last step (if accepted)
        }
        if(h<10.0*uround) //too small h
        {
            cout<<"Tolerances are too small."<<endl;
            return true;  
        }
        
        //compute spectral radius every 25 steps if it isn't constant, or after
        //a step rejection
        if(rho_outdated) 
            update_rho();
            
        //The number of stages chosen depending on h and spectral radius.
        update_n_stages_and_h(h);

        //Computation of an integration step.
        //this function is implemented in derived classes
        step(t,h);
        // at this point yn contains solution at time t, 
        // ynpu is solution at time t+h
        
        n_steps++; //update number of steps
        
        //store last used number of stages
        sold=s;

        if(dt_adaptivity) // if time step adaptivity enabled
            compute_hnew(h);//compute new h
        else //we are not using time step adaptivity
            hnew=h; //keeps same h

        //Accepted step or without time step adaptivity
        if(err<nu || !dt_adaptivity)
        {
            //do some statistics
            dt_max = max(dt_max,h);
            dt_min = min(dt_min,h);
            
            //procedure for accepted step
            accepted_step(t,h);
            
            //ynpu becomes yn
            swap_ptr = yn; 
            yn = ynpu;      //yn takes the new value 
            ynpu = swap_ptr; //and ynpu will be free memory
            
            ode->update_time(t);
            
            if((last && output_frequency>=0) || (output_frequency>0 && acc_steps%output_frequency==0))
                ode->output_solution();
            
            if(last)
                break;
        }
        else
            rejected_step(t,h);
    }
      
    elapsed_time = get_cpu_time()-elapsed_time;
    
    return last;
}     

 void OdeRungeKuttaIntegrator::update_rho()
{
    /**
     * A new spectral radius is computed. Either with the estimation given by
     * ODE::rho or with the RungeKuttaIntegrator::rho internal power method.
     */
    
    if((Astable || ode->constant_rho()) && n_steps>0)
        return;

    if(verbose) cout<<"\n--------------   Spectral Radius Estimation   --------------"<<endl;
    if(!Astable)
    {
        int iter=0;
        //Computed externally by ODE::rho
        if (!internal_rho)
            ode->rho(eigmax);
        //Computed internally by this->rho
        else
            this->rho(eigmax,iter);

        //recover statistics
        if ((int)eigmax+1>max_rho) max_rho=(int)eigmax+1;
        if ((int)eigmax+1<min_rho) min_rho=(int)eigmax+1;

        nrho=0;

        if(verbose) cout<<"Spectral radius estimation: "<<(int)eigmax+1<<endl;
        if(internal_rho && verbose)
            cout<<"Power method converged in "<<iter<<" iterations."<<endl;
    }
    else if(verbose)
        cout<<"Spectral radius not computed since "<<solver<<" is A-stable."<<endl;
    
    if(verbose) cout<<"------------------------------------------------------------\n"<<endl;
    
    rho_outdated=false;
}

 void OdeRungeKuttaIntegrator::compute_hnew(Real& h)
{
    /**
     * Estimates the optimal time step size using the local error estimation
     * of the current and previous step. 
     */
    
    //we use this estimation if the previous step has been accepted and it isn't
    //the first error estimation that we do
    if(!reject && !firsterrorestimation)
    {
        facp=1.0/err;
        fac=facp*(h/hp);
        fac=errp*fac*fac;

        fac=min(facp,fac);
        fac=sqrt(fac);
    }
    else //in first step we go here, or if the previous step has been rejected
        fac=sqrt(1.0/err);

    fac=min(facmax,max(0.1,0.8*fac));
    hnew=h*fac;
}

 void OdeRungeKuttaIntegrator::accepted_step(Real& t, Real& h)
{
    /**
     * Called when the step is accepted. Does output, sets the new time step size
     * and updated some variables.
     */

    //define delta character and red color
    int stages = s+((solver=="ROCK2"||solver=="DampedROCK2") ? 2:0);
    std::string delta = u8"\u0394";
    string bcol="\033[31;1m";
    string ecol="\033[0m";
    if(dt_adaptivity && verbose)
        cout<<"Step t = "<<setw(6)<<setprecision(4)<<t<<", "<<delta<<"t = "<<setw(8)<<setprecision(6)<<h
        <<", s = "<<setw(3)<<stages<<". Accepted with ||e_n||_L2 = "<<setw(6)<<setprecision(4)<<err
        <<" and ||y_n||_Linf = "<<setw(6)<<setprecision(4)<<max(abs(yn->maxCoeff()),abs(yn->minCoeff()))<<endl;
    else if(verbose)
    {
        cout<<"Step t = "<<setprecision(4)<<t<<", "<<delta<<"t = "<<setw(8)<<setprecision(6)<<h
        <<", s = "<<setw(3)<<stages<<","<<(err>1? bcol:"")<<" ||e_n||_L2 = "<<setw(6)<<setprecision(4)<<err<<(err>1? ecol:"")
        <<" and ||y_n||_Linf = "<<setw(6)<<setprecision(4)<<max(abs(yn->maxCoeff()),abs(yn->minCoeff()))<<endl;
    }

    acc_steps++;
    facmax=2.0;  //h can double at most
    t=t+h;
    firsterrorestimation=false;
    errp=err;
    if(reject) //the previous time step has been rejected so a smaller h is chosen.
    {
        hnew = h>0.0 ? min(hnew,h):max(hnew,h); 
        reject=false; //this step has been accepted
        nrej=0;       //set the consecutive rejections to 0
    }
    hp=h;  //previous h
    h=hnew;//next h
    nrho=nrho+1; //consecutive steps without computing the spectral radius
    nrho=(nrho)%25; //set to 0 every 25 steps
    rho_outdated = nrho==0;
}

 void OdeRungeKuttaIntegrator::rejected_step(Real& t, Real& h)
{
    /**
     * Called when the step is rejected. Does output, chooses the new time step
     * and updates some variables.
     */
    
    //define delta character and red color
    int stages = s+((solver=="ROCK2"||solver=="DampedROCK2") ? 2:0);
    std::string delta = u8"\u0394";
    string bcol="\033[31;1m";
    string ecol="\033[0m";
    if(verbose)
        cout<<"Step t = "<<setw(6)<<setprecision(4)<<t<<", "<<delta<<"t = "<<setw(8)<<setprecision(6)<<h
        <<", s = "<<setw(3)<<stages<<"."<<bcol<<" Rejected with ||e_n||_L2 = "<<setw(6)<<setprecision(4)<<err<<ecol
        <<" and ||y_n||_Linf = "<<setw(6)<<setprecision(4)<<max(abs(yn->maxCoeff()),abs(yn->minCoeff()))<<ecol<<endl;
            
    rej_steps++;
    reject=true;
    last=false;     //the step is rejected, it cant be the last one
    h=0.8*hnew;
    facmax=2.0;     //next step step h can at most double

    if(told==t)     //the previous step was also rejected
    {
        nrej=nrej+1;//consecutive rejections
        if (nrej==10)//after 10 consecutive rejections
            h=h*1e-3;
    }
    else
        told=t; //first rejection
    
    //the step has been rejected, so we recompute the spectral radius, unless
    //it has just been computed (nrho==0)
    if(nrho>0)
        rho_outdated=true;
}
 
Real OdeRungeKuttaIntegrator::get_cpu_time()
{
    return (Real)clock() / CLOCKS_PER_SEC;
}

int OdeRungeKuttaIntegrator::getAcc_steps()
{
    return acc_steps;
}

int OdeRungeKuttaIntegrator::getRej_steps()
{
    return rej_steps;
}

string OdeRungeKuttaIntegrator::getSolver()
{
    return solver;
}

Real OdeRungeKuttaIntegrator::get_elapsed_time()
{
    return elapsed_time;
}

int OdeRungeKuttaIntegrator::check_correctness(Real& h)
{
    /**
     * Some parameters checks before starting integration.
     */
    
    //Test the initial step size and tolerances.--------     
    if (h<10.0*uround) 
    {
        cout<<"Initial step-size is too small."<<endl;
        return 0;
    } 
    if (scalar_tol) //scalar tolerances
    {
        if (a_tol<=0.0 || r_tol<=10.0*uround)
        {
            cout<<"Tolerances are too small."<<endl;
            return 0;
        }
    }
    else
    {
        cerr<<"NON SCALAR TOLERANCES NOT IMPLEMENTED YET"<<endl;
        return 0;
        /*
        for(int i=0;i<neqn;i++)
        {
            if (atol[i]<=0.0 || rtol[i]<=10.0*uround)
            {
                cout<<"Tolerances are too small."<<endl;
                idid=-1;
                return;
            }
        }
        */
    }
    
    //everything is ok
    return 1;
}

void OdeRungeKuttaIntegrator::rho(Real& eigmax, int& iter)
{
    /**
    *     rho computes eigmax, a close upper bound of the
    *     spectral radius of the Jacobian matrix using a 
    *     power method (J.N. Franklin (matrix theory)). 
    *     The algorithm used is a small change (initial vector
    *     and stopping criteria) of that of
    *     Sommeijer-Shampine-Verwer, implemented in RKC.
    */

    Real eigmaxo,sqrtu,znor,ynor,quot,dzyn,dfzfn;
    
    const int maxiter=100;
    const Real safe=1.05;
    const Real tol = 1e-3;
    Real t = ode->get_t();
    
    sqrtu= sqrt(uround);

// ------ The initial vectors for the power method are yn --------
//       and yn+c*f(v_n), where vn=f(yn) a perturbation of yn 
//       (if n_steps=0) or a perturbation of the last computed
//       eigenvector (if n_steps!=0). 
    
    Vector*& fn = integr[0];
    Vector* z  = integr[1];
    Vector* swap_ptr;
    
    srand (time(NULL));
    if(n_steps==0)//take a random vector as first vector of power method
    {
        for(int i=0;i<ode->system_size();i++)
            (*z)(i) = ((rand()%2)-0.5)*((double)rand())/((double)RAND_MAX);
    }
    else//take the last vector of last power method call
        *z=*eigenvector;
        
    ode->f(t,*yn,*fn);
    
// ------ Perturbation.--------
    ynor= yn->norm();
    znor= z->norm();
    
    int k;
    
    // Building the vector z so that the difference z-yn is small
    if(ynor!=0.0 && znor!=0.0)
    {
        dzyn=ynor*sqrtu;
        quot=dzyn/znor;
        (*z) *= quot;
        (*z) += *yn;
    }
    else if(ynor!=0.0)
    {
        dzyn=ynor*sqrtu;
        *z=*yn;
        (*z) *= 1.+sqrtu;
    }
    else if(znor!=0.0)
    {
        dzyn=sqrtu;
        quot=dzyn/znor;
        (*z) *= quot;
    }
    else
    {
        dzyn=sqrtu*sqrt(z->size());
        for(int i=0;i<ode->system_size();i++)
            (*z)(i) += sqrtu;
    }
    //here dzyn=||z-yn|| and z=yn+(small perturbation)
    //dzyn=||z-yn|| will be always true, even with new z in the loop
    //Start the power method for non linear operator rhs
        
    eigmax=0.0;
    for(iter=1;iter<=maxiter;iter++)
    {
        ode->f(t,*z,*eigenvector);
        n_f_eval_rho++;

        (*eigenvector) -= *fn; //dz is the new perturbation, not normalized yet
        dfzfn= eigenvector->norm();
                
        eigmaxo=eigmax;
        eigmax=dfzfn/dzyn; //approximation of the Rayleigh quotient (not with dot product but just norms)
        eigmax=safe*eigmax;            
                        
        if (iter>=2 && abs(eigmax-eigmaxo)<= eigmax*tol)
        {
            //The last perturbation is stored. It will very likely be a
            // good starting point for the next rho call.
            *eigenvector=*z;
            (*eigenvector) -= *yn;    
            return;
        }
        if (dfzfn!=0.0)
        {
            quot=dzyn/dfzfn;
            *z=*eigenvector;
            (*z) *= quot;
            (*z) += *yn; //z is built so that dzyn=||z-yn|| still true
        }
        else
            break;
    }
    iter--;
    
    cout<<"ERROR: Convergence failure in the spectral radius computation."<<endl;
}

void OdeRungeKuttaIntegrator::print_info()
{
    /**
     * Information about the chosen RungeKuttaIntegrator settings.
     */
    // Printing integration parameters
    
    cout<<"\n-------------------   Ode Solver Info   --------------------"<<endl;
    cout<<"ODE solver: "<<solver<<endl;
    cout<<"Time-step adaptivity "<<(dt_adaptivity ? "enabled.":"disabled.")<<endl;
    cout<<"Spectral radius computed "<<(internal_rho ? "internally.":"externally.")<<endl;
    cout<<"Absolute tolerance = "<<a_tol<<endl;
    cout<<"Relative tolerance = "<<r_tol<<endl;
    cout<<"------------------------------------------------------------\n"<<endl;
}

void OdeRungeKuttaIntegrator::print_statistics()
{
    /**
     * Some statistics about the time integration.
     */
    
    cout<<"\n----------------   Integration Statistics   ----------------"<<endl;
    
    if(Astable)
        cout<<"The spectral radius has not been computed, "<<solver<<" is A-stable."<<endl;
    else
    {
        cout<<"Max estimation of the spectral radius: "<<max_rho<<endl;
        cout<<"Min estimation of the spectral radius: "<<min_rho<<endl;
    }
    cout<<"Number of f eval. for the spectr. radius = "<<n_f_eval_rho<<endl;
    cout<<"Max number of stages used: "<<max_s<<endl;
    cout<<"Number of f total evaluations = "<<n_f_eval<<endl;
    cout<<"Maximal time step used: "<<dt_max<<endl;
    cout<<"Steps: "<<n_steps<<endl;
    cout<<"Accepted steps: "<<acc_steps<<endl;
    cout<<"Rejected steps: "<<rej_steps<<endl; 
    cout<<"Elapsed time: "<<elapsed_time<<endl;
    cout<<"------------------------------------------------------------\n"<<endl;
    
    ofstream out(out_file+string("_statistics.txt"), ofstream::out);
    out<<setprecision(16);
    out<<elapsed_time<<", ";
    out<<n_f_eval<<", ";
    out<<max_s<<", ";
    out<<max_rho<<", ";
    out<<n_steps<<", ";
    out<<acc_steps<<", ";
    out<<rej_steps;
    out.close();
}

Real OdeRungeKuttaIntegrator::get_dt_max()
{
    return dt_max;
}

int OdeRungeKuttaIntegrator::get_n_f_eval()
{
    return n_f_eval;
}

RungeKutta4::RungeKutta4(Ode* ode, bool verb, int output_freq)
:OdeRungeKuttaIntegrator(ode, verb, false, 1e-2, 1e-2, true, output_freq)
{
    solver = "RK4";
    s = 4;
    max_s = 4;
}

void RungeKutta4::step(const Real t, const Real& h)
{
    
    Vector*& k1= integr[0];
    Vector*& k2= integr[1];
    Vector*& k3= integr[2];
    Vector*& k4= integr[3];
    Vector*& tmp= integr[4];    
        
    ode->f(t,*yn,*k1);
    
    *tmp = *yn;
    *tmp+= h/2.*(*k1);
    ode->f(t+h/2.,*tmp,*k2);
    
    *tmp = *yn;
    *tmp+= h/2.*(*k2);
    ode->f(t+h/2.,*tmp,*k3);
    
    *tmp = *yn;
    *tmp+= h*(*k3);
    ode->f(t+h,*tmp,*k4);
    
    *ynpu = *yn;
    *ynpu += h/6.*(*k1);
    *ynpu += h/3.*(*k2);
    *ynpu += h/3.*(*k3);
    *ynpu += h/6.*(*k4);
    
    n_f_eval += 4;
}

void RungeKutta4::update_n_stages_and_h(Real& h)
{
    if(h>0.9*2.75/(this->eigmax))
    {
        h = 0.9*2.75/(this->eigmax);
        this->last=false;
    }
}
