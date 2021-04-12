#!/usr/bin/python

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.integrate import simps
from matplotlib.backends.backend_pgf import FigureCanvasPgf
matplotlib.backend_bases.register_backend('pdf', FigureCanvasPgf)

pt_XII_0 = {}
pt_XII_15 = {}
pt_XII_100 = {}
pt_XI_0 = {}
pt_XI_15 = {}
pt_XI_100 = {}
pt_VIII_0 = {}
pt_VIII_15 = {}
pt_VIII_100 = {}
expe_XII = {}
expe_XI = {}
expe_VIII = {}

# Read File FXII experimental data
readFile = open("FXII_Contact_Donnees_TGT.csv", 'r')
sepFile_nou = readFile.read().split('\n')
readFile.close()
columns = sepFile_nou[0].split(';')
print(columns)
len_x = len(columns)
ic = 0
for x in range(0,len_x):
    t1 = []
    for plotPair1 in sepFile_nou[1:]:
        xSTRING_1 = plotPair1.split(';')
        if len(xSTRING_1) > 1:
           t1.append(float(xSTRING_1[ic]))

    expe_XII[columns[ic]]= t1
    ic = ic+1

# Read File FXI experimental data
readFile = open("FXI_Contact_Donnees_TGT.csv", 'r')
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
    expe_XI[columns[ic]]= t1
    ic = ic+1

# Read File FVIII experimental data
readFile = open("FVIII_Contact_Donnees_TGT.csv", 'r')
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
    expe_VIII[columns[ic]]= t1
    ic = ic+1

# Simulation FVIII
readFile = open("Chatterjee_initial_FVIII/FVIII_100_pt1_IIa.dat", 'r')
sepFile_nou = readFile.read().split('\n')
readFile.close()
columns = sepFile_nou[0].split()
len_x = len(columns)-1
ic = 0
for x in range(0,len_x):
    t1 = []
    for plotPair1 in sepFile_nou[1:]:
        xSTRING_1 = plotPair1.split()
        if ic!=1:
           if len(xSTRING_1) > 1:
              t1.append(float(xSTRING_1[ic]))
    pt_VIII_100[columns[ic+1]]= t1
    ic = ic+1

readFile = open("Chatterjee_initial_FVIII/FVIII_15_pt1_IIa.dat", 'r')
sepFile_nou = readFile.read().split('\n')
readFile.close()
columns = sepFile_nou[0].split()
len_x = len(columns)-1
ic = 0
for x in range(0,len_x):
    t1 = []
    for plotPair1 in sepFile_nou[1:]:
        xSTRING_1 = plotPair1.split()
        if ic!=1:
           if len(xSTRING_1) > 1:
              t1.append(float(xSTRING_1[ic]))
    pt_VIII_15[columns[ic+1]]= t1
    ic = ic+1

readFile = open("Chatterjee_initial_FVIII/FVIII_0_pt1_IIa.dat", 'r')
sepFile_nou = readFile.read().split('\n')
readFile.close()
columns = sepFile_nou[0].split()
len_x = len(columns)-1
ic = 0
for x in range(0,len_x):
    t1 = []
    for plotPair1 in sepFile_nou[1:]:
        xSTRING_1 = plotPair1.split()
        if ic!=1:
           if len(xSTRING_1) > 1:
              t1.append(float(xSTRING_1[ic]))
    pt_VIII_0[columns[ic+1]]= t1
    ic = ic+1

#Simulations Factor FXII
readFile = open("Chatterjee_initial_FXII/FXII_100_pt1_IIa.dat", 'r')
sepFile_nou = readFile.read().split('\n')
readFile.close()
columns = sepFile_nou[0].split()
len_x = len(columns)-1
ic = 0
for x in range(0,len_x):
    t1 = []
    for plotPair1 in sepFile_nou[1:]:
        xSTRING_1 = plotPair1.split()
        if ic!=1:
           if len(xSTRING_1) > 1:
              t1.append(float(xSTRING_1[ic]))
    pt_XII_100[columns[ic+1]]= t1
    ic = ic+1

readFile = open("Chatterjee_initial_FXII/FXII_15_pt1_IIa.dat", 'r')
sepFile_nou = readFile.read().split('\n')
readFile.close()
columns = sepFile_nou[0].split()
len_x = len(columns)-1
ic = 0
for x in range(0,len_x):
    t1 = []
    for plotPair1 in sepFile_nou[1:]:
        xSTRING_1 = plotPair1.split()
        if ic!=1:
           if len(xSTRING_1) > 1:
              t1.append(float(xSTRING_1[ic]))
    pt_XII_15[columns[ic+1]]= t1
    ic = ic+1

readFile = open("Chatterjee_initial_FXII/FXII_0_pt1_IIa.dat", 'r')
sepFile_nou = readFile.read().split('\n')
readFile.close()
columns = sepFile_nou[0].split()
len_x = len(columns)-1
ic = 0
for x in range(0,len_x):
    t1 = []
    for plotPair1 in sepFile_nou[1:]:
        xSTRING_1 = plotPair1.split()
        if ic!=1:
           if len(xSTRING_1) > 1:
              t1.append(float(xSTRING_1[ic]))
    pt_XII_0[columns[ic+1]]= t1
    ic = ic+1

#Simulations Factor FXI
readFile = open("Chatterjee_initial_FXI/FXI_100_pt1_IIa.dat", 'r')
sepFile_nou = readFile.read().split('\n')
readFile.close()
columns = sepFile_nou[0].split()
len_x = len(columns)-1
ic = 0
for x in range(0,len_x):
    t1 = []
    for plotPair1 in sepFile_nou[1:]:
        xSTRING_1 = plotPair1.split()
        if ic!=1:
           if len(xSTRING_1) > 1:
              t1.append(float(xSTRING_1[ic]))
    pt_XI_100[columns[ic+1]]= t1
    ic = ic+1

readFile = open("Chatterjee_initial_FXI/FXI_15_pt1_IIa.dat", 'r')
sepFile_nou = readFile.read().split('\n')
readFile.close()
columns = sepFile_nou[0].split()
len_x = len(columns)-1
ic = 0
for x in range(0,len_x):
    t1 = []
    for plotPair1 in sepFile_nou[1:]:
        xSTRING_1 = plotPair1.split()
        if ic!=1:
           if len(xSTRING_1) > 1:
              t1.append(float(xSTRING_1[ic]))
    pt_XI_15[columns[ic+1]]= t1
    ic = ic+1

readFile = open("Chatterjee_initial_FXI/FXI_0_pt1_IIa.dat", 'r')
sepFile_nou = readFile.read().split('\n')
readFile.close()
columns = sepFile_nou[0].split()
len_x = len(columns)-1
ic = 0
for x in range(0,len_x):
    t1 = []
    for plotPair1 in sepFile_nou[1:]:
        xSTRING_1 = plotPair1.split()
        if ic!=1:
           if len(xSTRING_1) > 1:
              t1.append(float(xSTRING_1[ic]))
    pt_XI_0[columns[ic+1]]= t1
    ic = ic+1

# Build Simulation probes
XII_0 = np.asarray(pt_XII_0["04:IIa"])*1e+9  #+ 1.2*np.asarray(pt_mIIa_0["04:mIIa"])
XII_15 = np.asarray(pt_XII_15["04:IIa"])*1e+9  #+ 1.2*np.asarray(pt_mIIa_0["04:mIIa"])
XII_100 = np.asarray(pt_XII_100["04:IIa"])*1e+9  #+ 1.2*np.asarray(pt_mIIa_0["04:mIIa"])

XI_0 = np.asarray(pt_XI_0["04:IIa"])*1e+9  #+ 1.2*np.asarray(pt_mIIa_0["04:mIIa"])
XI_15 = np.asarray(pt_XI_15["04:IIa"])*1e+9  #+ 1.2*np.asarray(pt_mIIa_0["04:mIIa"])
XI_100 = np.asarray(pt_XI_100["04:IIa"])*1e+9  #+ 1.2*np.asarray(pt_mIIa_0["04:mIIa"])

VIII_0 = np.asarray(pt_VIII_0["04:IIa"])*1e+9  #+ 1.2*np.asarray(pt_mIIa_0["04:mIIa"])
VIII_15 = np.asarray(pt_VIII_15["04:IIa"])*1e+9  #+ 1.2*np.asarray(pt_mIIa_0["04:mIIa"])
VIII_100 = np.asarray(pt_VIII_100["04:IIa"])*1e+9  #+ 1.2*np.asarray(pt_mIIa_0["04:mIIa"])

time_exp_XII = np.asarray(expe_XII["time"])*60.0
time_exp_XI = np.asarray(expe_XI["time"])*60.0
time_exp_VIII = np.asarray(expe_VIII["time"])*60.0

font = {'family' : 'serif',
        'size'   : 25}
plt.rc('font', **font)

a1=plt.subplot(3,1,1)
plt.title('FXII')
plt.plot(pt_XII_100["01:time"],XII_100,'k-',lw=2.5, label= r'Num. $100 \%$')
plt.plot(pt_XII_15["01:time"],XII_15,'k-.',lw=2.5, label= r'Num. $15 \%$')
plt.plot(pt_XII_0["01:time"],XII_0,'k--',lw=2.5, label= r'Num. $0 \%$')
plt.plot(time_exp_XII,expe_XII["FXII 100%"],'kp',lw=2.5,label=r'Exp. $100\%$ ', markevery=1)
plt.plot(time_exp_XII,expe_XII["FXII 15%"],'k>',lw=2.5,label=r'Exp. $15\%$ ', markevery=1)
plt.plot(time_exp_XII,expe_XII["FXII 0%"],'ks',lw=2.5,label=r'Exp. $0\%$ ', markevery=1)
plt.ylabel(r'$II_a$ (nM)')
plt.xlim((0,1000))
plt.legend(loc='upper right', prop={'size':12})
plt.setp(a1.get_xticklabels(), visible=False)

a2=plt.subplot(3,1,2)
plt.title('FXI')
plt.plot(pt_XI_100["01:time"],XI_100,'k-',lw=2.5, label= r'Numerical $FXI_{100 \%}$')
plt.plot(pt_XI_15["01:time"],XI_15,'k-.',lw=2.5, label= r'Numerical $FXI_{15 \%}$')
plt.plot(pt_XI_0["01:time"],XI_0,'k--',lw=2.5, label= r'Numerical $FXI_{0 \%}$')
plt.plot(time_exp_XI,expe_XI["FXI 100%"],'kp',lw=2.5,label=r'Experimental $FXI_{100\%}$ ', markevery=1)
plt.plot(time_exp_XI,expe_XI["FXI 15%"],'k>',lw=2.5,label=r'Experimental $FXI_{15\%}$ ', markevery=1)
plt.plot(time_exp_XI,expe_XI["FXI 0%"],'ks',lw=2.5,label=r'Experimental $FXI_{0\%}$ ', markevery=1)
plt.ylabel(r'$II_a$ (nM)')
plt.xlim((0,1000))
plt.setp(a2.get_xticklabels(), visible=False)

a3=plt.subplot(3,1,3)
plt.title('FVIII')
plt.plot(pt_VIII_100["01:time"],VIII_100,'k-',lw=2.5, label= r'Numerical $FVIII_{100 \%}$')
plt.plot(pt_VIII_15["01:time"],VIII_15,'k-.',lw=2.5, label= r'Numerical $FVIII_{15 \%}$')
plt.plot(pt_VIII_0["01:time"],VIII_0,'k--',lw=2.5, label= r'Numerical $FVIII_{0 \%}$')
plt.plot(time_exp_VIII,expe_VIII["FVIII 100%"],'kp',lw=2.5,label=r'Experimental $FVIII_{100\%}$ ', markevery=1)
plt.plot(time_exp_VIII,expe_VIII["FVIII 15%"],'k>',lw=2.5,label=r'Experimental $FVIII_{15\%}$ ', markevery=1)
plt.plot(time_exp_VIII,expe_VIII["FVIII 0%"],'ks',lw=2.5,label=r'Experimental $FVIII_{0\%}$ ', markevery=1)
plt.ylabel(r'$II_a$ (nM)')
plt.xlabel(r'time (s) ')
plt.xlim((0,1000))

plt.show()
