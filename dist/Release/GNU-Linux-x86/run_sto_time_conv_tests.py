#!/usr/bin/env python
import sys
import os
import shutil
import time

cmdargs = sys.argv

p = [0, 1, 2, 3, 4]

ntest = cmdargs[1]     #test number
mtiter = cmdargs[2]    #monte carlo iterations
ode_solver = "SROCK2"
contW = 0              #discrete brownian motion

program = "stochastic_adaptive_time_stepping"


# -------------DO NOT MODIFY BELOW THIS LINE---------------------
options = "-sde 1 -contW "+str(contW) + \
          " -verb 0 -dtadap 0 -ofreq 0 -rk "+ode_solver + \
          " -ntest "+str(ntest)+" -iter "+str(mtiter)+" -dt "
solutionfile = "TimeConvTest_"+ode_solver+"_dt_"


for k in p:
    dt = 1.0/2.0**k
    #print("Running Monte Carlo for dt = 2^(-"+str(k)+")")
    filename = solutionfile+str(k)
    command = "./"+program+" "+options+str(dt)+" -ofile "+filename
    #print(command)
    os.system(command)
    print("\n\n\n")
#print("Simulation Finished.")


    


