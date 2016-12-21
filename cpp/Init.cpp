#include "Init.h"


Init::Init()
{
}        

void Init::initOde(Ode*& ode)
{
    if(ntest==1)
        ode = new OneDimensionalTest1(output_file);
    else if(ntest==2)
        ode = new NeuronCable(output_file);
    else if(ntest==3)
        ode = new TwoDimensionalBrusselator(output_file);
    else if(ntest==4)
        ode = new PDEBrusselator(output_file);
    else if(ntest==5)
        ode = new Krogh10(output_file);
    else
        cout<<"Problem not known"<<endl;
}

void Init::initSde(Sde*& sde)
{
    if(ntest==1)
        sde = new SdeOneDimensionalTest1(output_file);
    else if(ntest==2)
        sde = new SdeOneDimensionalTest2(output_file);
    else if(ntest==3)
        sde = new StochasticTwoBodyProblem(output_file);
    else if(ntest==4)
        sde = new SdeLinear(output_file);
    else
        cout<<"Problem not known"<<endl;
}

void Init::initBrownianMotion(BrownianMotion*& W, Sde* sde, TimeIntegrator* integr)
{
    if(continuous){
			if (dtadap)
        W = new ContinuousAdaptedBrownianMotion(sde->brownian_size(), 
                                         integr->needDoubleIntegral(),
                                         sde->is_commutative(),
                                         sde->isDiagonal());
			else
        W = new ContinuousBrownianMotion(sde->brownian_size(), 
                                         integr->needDoubleIntegral(),
                                         sde->is_commutative(),
                                         sde->isDiagonal());
		}
    else{
			if ( dtadap)
        W = new DiscreteAdaptedBrownianMotion(sde->brownian_size(), 
                                       integr->needDoubleIntegral(),
                                       sde->is_commutative(),
                                       sde->isDiagonal());
			else
        W = new DiscreteBrownianMotion(sde->brownian_size(), 
                                       integr->needDoubleIntegral(),
                                       sde->is_commutative(),
                                       sde->isDiagonal());
		}
}

void Init::initTimeIntegrator(TimeIntegrator*& rk, Ode* ode)
{
    if(rk_name=="ROCK2")
        rk = new ROCK2(ode, verbose, dtadap, atol, rtol, 
                       scalartol, output_freq);
    else if(rk_name=="DampedROCK2" || rk_name=="DROCK2")
        rk = new DROCK2(ode, verbose, dtadap, atol, rtol, 
                       scalartol, output_freq);
    else if(rk_name=="RK4")
        rk = new RungeKutta4(ode, verbose, output_freq);
    else if(rk_name=="SROCK2")
    {
        Sde* sde = dynamic_cast<Sde*>(ode);
        rk = new SROCK2(sde, verbose, dtadap, atol, rtol, 
                       scalartol, output_freq);
    }
    else if(rk_name=="MT")
    {
        Sde* sde = dynamic_cast<Sde*>(ode);
        rk = new MilsteinTalay(sde, verbose, dtadap, atol, rtol, 
                       scalartol, output_freq);
    }
    else
        cout<<"Integrator "<<rk_name<<" not known"<<endl;
    // RKC is not working yet
//    else if(rk_name=="RKC")
//        rk = new RKC(ode, onestep, verb, dtadap, atol, rtol, intrho, 
//                       scalartol, output_freq);
}

bool Init::read_command_line(GetPot& command_line, int& ntest, Real& dt, Real& rtol,
                       Real& atol, bool& dt_adaptivity,
                       string& solver, int& output_frequency,
                       string& output_file, int& iter, bool& continuous_brow,
                       bool& deterministic, bool& verbose)
{
    if(command_line.search(2,"-ntest","-test"))
        ntest = command_line.next(ntest);
    if(command_line.search("-dt"))
        dt = command_line.next(dt);
    if(command_line.search("-iter"))
        iter = command_line.next(iter);
    if(command_line.search("-contW"))
        continuous_brow = command_line.next(continuous_brow);
    if(command_line.search("-ode"))
        deterministic = command_line.next(deterministic);
    if(command_line.search("-sde"))
        deterministic = !command_line.next(!deterministic);
    if(command_line.search(2,"-rtol","-rt"))
    {
        rtol = command_line.next(rtol);
        if(!command_line.search(2,"-atol","-at"))
            atol=rtol;
    }
    if(command_line.search(2,"-atol","-at"))
    {
        atol = command_line.next(atol);
        if(!command_line.search(2,"-rtol","-rt"))
            rtol=atol;
    }
    if(command_line.search(3,"-dtadaptivity","-dtadap","-dta"))
        dt_adaptivity = command_line.next(dt_adaptivity);
    if(command_line.search(2,"-verbose","-verb"))
        verbose = command_line.next(verbose);
    if(command_line.search(2,"-outputfile","-ofile"))
        output_file = command_line.next(output_file.c_str());
    if(command_line.search(2,"-outputfreq","-ofreq"))
        output_frequency = command_line.next(output_frequency);
    if(command_line.search(5,"-odesolver","-solver","-integrator","-rk","-sdesolver"))
            solver = command_line.next(solver.c_str());
    
    if(deterministic)//no monte carlo
        iter = 1;
    if(iter>1)//no verbose nor output if doing monte carlo
    {
        verbose=false;
        output_frequency = -1;
    }
    
    if(command_line.search("-help"))
    {
        cout<<"This is a code for running ODE and SDE simulations. A list of problems\n"
            <<"is hardcoded in the executable. Look into OdeList.h and SdeList.h.\n"
            <<"Each problem has a number n depending on the order of apparition in the list."<<endl;
        cout<<"The following options are available:\n"
            <<"    -sde s     : if s=1 run a SDE simulation, if s=0 run a ODE simulation.\n"
            <<"    -ntest n   : chooses the problem. If s=1 it takes the nth problem\n"
            <<"                 from list SdeList.h, else the nth from OdeList.h.\n"
            <<"    -dt h      : time step when running in fixed time step mode or initial time\n"
            <<"                 step when running in adaptive time step mode.\n"
            <<"    -dtadap dta: if dta=1 enables time step adaptivity, else in fixed time step mode.\n"
            <<"    -contW W   : if W=1 uses a continuous brownian motion, else a discrete one.\n"
            <<"    -atol at   : sets the absolute tolerance, by default at=1e-2.\n"
            <<"    -rtol rt   : sets the relative tolerance. If not provided then rtol=atol.\n"
            <<"    -ofile ofi : name of output file. By default ofi=solution.\n"
            <<"    -iter it   : if s=1 then it is the number of Monte Carlo simulations. Else\n"
            <<"                 it is ignored.\n"
            <<"    -ofreq ofr : frequency of write to disk. If ofr=0 writes at end of simulation\n"
            <<"                 only. If ofr=-1 never writes solution. If it>1 we set ofr=-1.\n"
            <<"    -verb v    : enables or disables verbosity. If it>1 then v=0.\n"
            <<"    -solver rk : name of RungeKutta solver to use. The available solvers are:\n"
            <<"                 -ROCK2    Runge-Kutta-Orthogonal-Chebychev order 2,\n"
            <<"                 -DROCK2   ROCK2 with increased damping,\n"
            <<"                 -SROCK2   Stochastic Weak Order 2 ROCK2,\n"
            <<"                 -MT       Weak order 2 Milstein-Talay.\n"<<endl;
        cout<<"Error estimators are available for ROCK2 and SROCK2. Not implemented in DROCK2 and MT."<<endl;
        cout<<"Files are written into folders with same name of the tests."<<endl;
        return true;
    }
    
    this->ntest=ntest;
    this->dt = dt;
    this->rtol = rtol;
    this->atol=atol;
    this->dtadap = dt_adaptivity;
    this->rk_name = solver;
    this->output_freq = output_frequency;
    this->output_file = output_file;
    this->iter = iter;
    this->continuous = continuous_brow;
    this->deterministic = deterministic;
    this->verbose = verbose;
    
    return false;
}

void Init::setOutput_file(string output_file)
{
    this->output_file = output_file;
}

string Init::get_output_filename()
{
    return output_file;
}

string Init::getRk_name()
{
    return rk_name;
}

bool Init::isDtadap() 
{
    return dtadap;
}

Real Init::getAtol() 
{
    return atol;
}

Real Init::getRtol() 
{
    return rtol;
}

Real Init::getDt() const
{
    return dt;
}

int Init::getNtest() const
{
    return ntest;
}

bool Init::isContinuous() const
{
    return continuous;
}
