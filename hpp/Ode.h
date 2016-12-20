#ifndef ODE_H
#define	ODE_H

#include <iostream>
#include <iomanip>
#include <Eigen/Eigen>

using namespace std;
typedef double Real;
typedef Eigen::VectorXd Vector;

class Ode
{
public:
    Ode(string solutionfile);
    virtual ~Ode();
    
    Vector*& get_y();
    Real get_t();
    Real get_t0();
    Real get_tend();
        
    virtual void init_solution() = 0;
    virtual void f(Real t, Vector& x, Vector& fx) =0;
    
    virtual void rho(Real& eigmax);
    bool do_know_rho() const;
    
    void update_time(Real new_time);
    int system_size();
    bool constant_rho();
    
    void output_solution();
    string get_solution_file_name();
    
protected:
    Vector* y;
    Real t;
    Real t0;
    Real tend;
    
    int neqn;
    bool cte_rho;
    bool know_rho;
    
    string solution_file_name;
};


#endif	/* ODE_H */

