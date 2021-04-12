#!/usr/bin/python

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.integrate import simps
from matplotlib.backends.backend_pgf import FigureCanvasPgf
matplotlib.backend_bases.register_backend('pdf', FigureCanvasPgf)

pt_IIa = {}

expe = {}

readFile = open("FVIII_FT_Donnees_TGT.csv", 'r')
sepFile_nou = readFile.read().split('\n')
readFile.close()
columns = sepFile_nou[0].split(',')
print(columns)
len_x = len(columns)
ic = 0
for x in range(0,len_x):
    t1 = []
    for plotPair1 in sepFile_nou[1:]:
        xSTRING_1 = plotPair1.split(',')
        if len(xSTRING_1) > 1:
           t1.append(float(xSTRING_1[ic]))
    expe[columns[ic]]= t1
    ic = ic+1

readFile = open("1200sPeak_pt1_IIa.dat", 'r')
sepFile_nou = readFile.read().split('\n')
readFile.close()
columns = sepFile_nou[0].split()
print(columns)
len_x = len(columns)-1
ic = 0
for x in range(0,len_x):
    t1 = []
    for plotPair1 in sepFile_nou[1:]:
        xSTRING_1 = plotPair1.split()
        if ic!=1:
           if len(xSTRING_1) > 1:
              t1.append(float(xSTRING_1[ic]))
    pt_IIa[columns[ic+1]]= t1
    ic = ic+1

IIa = np.asarray(pt_IIa["04:IIa"])*1e9  #+ 1.2*np.asarray(pt_mIIa_0["04:mIIa"]*1e+9)
time_exp = np.asarray(expe["time"])*60.0

font = {'family' : 'serif',
        #'weight' : 'bold',
        'size'   : 25}
#plt.rc('text', usetex=True)
plt.rc('font', **font)

#plt.rcParams.update({'font.size': 25})

plt.plot(pt_IIa["01:time"],IIa,'k-',lw=2.5,label='Numerical')
plt.plot(time_exp,expe["FVIII 100%"],'ko',lw=2.5,label=r'Experimental')
plt.xlim((0,2000))
plt.ylabel(r'$IIa$ nM ')
plt.xlabel(r'$t$ s ')
plt.legend(loc='upper right', prop={'size':25})
plt.show()

