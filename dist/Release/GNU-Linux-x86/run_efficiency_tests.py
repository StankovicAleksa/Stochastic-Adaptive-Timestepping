#!/usr/bin/env python
import sys
import os
import shutil

cmdargs = sys.argv

tol_pow = [2, 3, 4, 5, 6, 7, 8]
dt_pow = [1, 2, 4, 6, 8, 10, 12]

ntest = 5
ode_solver = "DROCK2"
sde = 0
intrho = 0

# Ignore this if sde=0
contW = 0
mtiter = "1"

program = "stochastic_adaptive_time_stepping"


# -------------DO NOT MODIFY BELOW THIS LINE---------------------
options = "-sde "+str(sde)+" -contW "+str(contW)+" -verb 0 -dtadap 1 -ofreq 0 -onestep 0 -intrho "+str(intrho)+" -rk "+ode_solver\
          + " -ntest " + str(ntest)+" -iter "+mtiter+" -dt 1e-2 -atol "
solutionfile = "EfficiencyTest_"+ode_solver+"_tol_10e_"

for k in tol_pow:
    tol = 1.0/10.0**k
    filename = solutionfile+str(k)
    command = "./"+program+" "+options+str(tol)+" -outputfile "+filename+" > /dev/null &"
    print(command)
    os.system(command)

options = "-sde "+str(sde)+" -contW "+str(contW)+" -verb 0 -dtadap 0 -ofreq 0 -onestep 0 -intrho "+str(intrho)+" -rk "+ode_solver\
          + " -ntest " + str(ntest)+" -iter "+mtiter+" -dt  "
solutionfile = "EfficiencyTest_"+ode_solver+"_dt_2e_"

for k in dt_pow:
    dt = 1.0/2.0**k
    filename = solutionfile+str(k)
    command = "./"+program+" "+options+str(dt)+" -outputfile "+filename+" > /dev/null &"
    print(command)
    os.system(command)


    


