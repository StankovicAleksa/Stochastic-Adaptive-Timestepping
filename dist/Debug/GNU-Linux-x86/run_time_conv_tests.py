#!/usr/bin/env python
import sys
import os
import shutil
import time

cmdargs = sys.argv

#p = [2, 4, 6, 8, 10, 12, 14, 16]
p = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]


ntest = cmdargs[1]
ode_solver = "ROCK2"


program = "stochastic_adaptive_time_stepping"


# -------------DO NOT MODIFY BELOW THIS LINE---------------------
options = "-sde 0 -contW 0" + \
          " -verb 0 -dtadap 0 -ofreq 0 -rk "+ode_solver + \
          " -ntest "+str(ntest)+" -dt "
solutionfile = "TimeConvTest_"+ode_solver+"_dt_"

for k in p:
    dt = 1.0/2.0**k
    filename = solutionfile+str(k)
    command = "./"+program+" "+options+str(dt)+" -ofile "+filename+" > /dev/null &"
    print(command)
    os.system(command)

    


