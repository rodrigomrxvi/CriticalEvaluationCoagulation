#
# Compute the time evolution of species given a biochemical scheme given in .xml
# This program includes some comparison with analytical results when available as
# well as the computation of some TGT outcomes if IIa is part of the species
#

import os
from operator import mul
from functools import reduce
from lxml import etree
import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import simps
from scipy.optimize import fsolve
import pylab
from fun_compute import *
from fun_io import *
import time as t


#
# Reading the .in and .xml files - See fun_io.py for description of contains
#
s0, amp0, INPUTFILE, OUTPUTFILE, OUTPUTFILE_GRAD, OUTPUTFILE_ADJOINT,\
final_time, number_of_time, scaling_factor, type_of_plot, \
method, max_dt, grad_type, pnorm, tfixed, valfixed, target, grad_cost_function, grad_check_species, \
grad_check_param, grad_check_eps = read_infile()

number_of_species, spec_name, y0, specs_in, specs_out, \
index_in, index_out, number_of_reactions, reac_name, rate, \
amp0 = read_xmlfile(INPUTFILE,s0,amp0)

print('Biochemistry read from ',INPUTFILE)
print('   ')

#
# Applying the scaling scaling factor to y0 and rate as appropriate
#
y0, rate = scaling(y0,rate,index_in,scaling_factor)

#
# Print the biochemistry details
#
print_biochemistry(spec_name,y0,reac_name,rate,index_in,index_out,scaling_factor)

#
# Solve the biochemical problem: dy/dt = source(y,t)
#
time = np.linspace(0,final_time,number_of_time)
tic = t.time()
out = solve_ivp(source_biochem,[0.0,final_time],y0,                       \
                method=method,                                             \
                t_eval=time,                                              \
                args=(rate, index_in, index_out),                           \
                max_step=max_dt,                                            \
                rtol=1.0e-10,                                              \
                atol=1.0e-10,                                              \
                jac=jac_source_biochem)
y = out.y
toc = t.time()
print("Solving biochemistry in "+"{:e}".format(toc-tic)+" seconds")


#
# Print the computational details
#
print_case(method,max_dt,out,OUTPUTFILE)

#
# Save results in .out file
#
store_results(method,max_dt,scaling_factor,out,spec_name,y0,y,reac_name,rate,\
              index_in,index_out,time,OUTPUTFILE)


#
# Generate analytical solution if available
#
yth = analytical(INPUTFILE,y0,rate,time)

#
# Compute TGT statistics
#
if 'IIa' in spec_name:
    indIIa = spec_name.index("IIa")
    indmIIa = spec_name.index("mIIa")
    totalIIa = simps(y[indIIa],time)
    print(" ")
    print(" TGT statistics:")
    print("Integral of IIa (M.second): ", totalIIa/scaling_factor)
    totalIIa = simps(y[indIIa]+1.2*y[indmIIa],time)
    print("Integral of IIa + 1.2 mIIa (M.second): ", totalIIa/scaling_factor)
    maxtotalIIa = max(y[indIIa]+y[indmIIa])
    #print("Peak value of IIa + mIIa: ", maxtotalIIa/scaling_factor)
    maxtotalIIabis = max(y[indIIa]+1.2*y[indmIIa])
    #print("Peak value of IIa + 1.2 mIIa: ", maxtotalIIabis/scaling_factor)

    tck = interpolate.splrep(time, y[indIIa], k=5,s=0)
    time_of_clot = find_time(time,y[indIIa],tck,2e-9*scaling_factor)
    print("Time of clot: ", time_of_clot," s")
    tmax_rate, max_rate = find_max_rate(time,y[indIIa],tck)
    print("Max of rate: ", max_rate/scaling_factor," at time ",tmax_rate," s")
    tmax, maxIIa = find_max(time,y[indIIa],tck)
    print("Peak value of IIa: ", maxIIa/scaling_factor," at time ",tmax," s")



#
# Generate analytical solution if available
#
yth = analytical(INPUTFILE,y0,rate,time)

#
# Plot species vs time
#
plot_species(type_of_plot,spec_name,time,y,yth,scaling_factor,'false')
pylab.show()
