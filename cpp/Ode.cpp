#include <fstream>

#include "Ode.h"

Ode::Ode(string solutionfile)
:solution_file_name(solutionfile)
{
}

void Ode::rho(Real& eigmax)
{
    cout<<"ERROR: using a non implemented rho function!\nUse -intrho 1."<<endl;
}

bool Ode::do_know_rho() const
{
    return know_rho;
}

Ode::~Ode()
{
    delete y;
}

Vector*& Ode::get_y()
{
    return y;
}

Real Ode::get_t()
{
    return t;
}

Real Ode::get_t0()
{
    return t0;
}

Real Ode::get_tend()
{
    return tend;
}

void Ode::update_time(Real new_time)
{
    t=new_time;
}

int Ode::system_size()
{
    return neqn;
}

bool Ode::constant_rho()
{
    return cte_rho;
}

string Ode::get_solution_file_name()
{
    return solution_file_name;
}

void Ode::output_solution()
{
    
    static int nout=1;
        
    if(nout==1)
    {
        ofstream outfile(solution_file_name+string(".m"), ofstream::out);
				std::cout << solution_file_name;
        outfile<<setprecision(16)<<"t("<<nout<<")="<<t<<";"<<endl;
        outfile<<"y("<<nout<<",:)=[";
        for(int i=0;i<neqn-1;i++)
            outfile<<(*y)(i)<<",";
        outfile<<(*y)(neqn-1)<<"];"<<endl;
        outfile.close();
    }
    else
    {
        ofstream outfile(solution_file_name+string(".m"), ofstream::out | ofstream::app);
        outfile<<setprecision(16)<<"t("<<nout<<")="<<t<<";"<<endl;
        outfile<<"y("<<nout<<",:)=[";
        for(int i=0;i<neqn-1;i++)
            outfile<<(*y)(i)<<",";
        outfile<<(*y)(neqn-1)<<"];"<<endl;
        outfile.close();
    }
    
    nout++;
    
    if(t==tend)
    {
        ofstream outfile(solution_file_name+string(".txt"), ofstream::out);
        outfile<<setprecision(16);
        for(int i=0;i<neqn-1;i++)
            outfile<<(*y)(i)<<", ";
        outfile<<(*y)(neqn-1);
        outfile.close();
    }
    
}
