#ifndef INIT_H
#define	INIT_H

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <GetPot>

#include "OdeProblems.h"
#include "SdeProblems.h"
#include "OdeStabilizedIntegrators.h"
#include "SdeStabilizedIntegrators.h"
#include "TraditionalSdeRungeKuttaIntegrators.h"

class Init
{
public:
    Init();
    
    void initOde(Ode*& ode);
    void initSde(Sde*& sde);
    void initBrownianMotion(BrownianMotion*& W, Sde* sde, TimeIntegrator* integr);
    void initTimeIntegrator(TimeIntegrator*& rk, Ode* ode);
    bool read_command_line(GetPot& command_line, int& ntest, Real& dt, Real& rtol,
                       Real& atol, bool& dt_adaptivity,
                       string& solver, int& output_frequency,
                       string& output_file, int& iter, bool& continuous_brow,
                       bool& deterministic, bool& verbose);
                       void setOutput_file(string output_file);
   string get_output_filename();
   string getRk_name();
   bool isDtadap();
   Real getAtol();
   Real getRtol();
   Real getDt() const;
   int getNtest() const;
   bool isContinuous() const;
protected:
    int ntest;
    string output_file;
    bool continuous;
    string rk_name;
    bool dtadap;
    Real atol;
    Real rtol;
    bool scalartol;
    int output_freq;
    Real dt;
    int iter;
    bool deterministic;
    bool verbose;
};


#endif	/* INIT_H */

