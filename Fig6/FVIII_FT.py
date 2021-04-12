#!/usr/bin/python

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.integrate import simps
from matplotlib.backends.backend_pgf import FigureCanvasPgf
matplotlib.backend_bases.register_backend('pdf', FigureCanvasPgf)

pt_IIa_0 = {}
pt_IIa_1 = {}
pt_IIa_5 = {}
pt_IIa_15 = {}
pt_IIa_50 = {}
pt_IIa_100 = {}
pt_Initial_IIa_0 = {}
pt_Initial_IIa_1 = {}
pt_Initial_IIa_5 = {}
pt_Initial_IIa_15 = {}
pt_Initial_IIa_50 = {}
pt_Initial_IIa_100 = {}
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

readFile = open("Hockin_modified_FVIII/FVIII_100_pt1_IIa.dat", 'r')
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
    pt_IIa_100[columns[ic+1]]= t1
    ic = ic+1
# Calculate ETP area under the curve
t=pt_IIa_100["01:time"]
dt=t[1]
dt1=(t[2]-t[1])/60
ETP=round(simps(pt_IIa_100["04:IIa"],dx=dt1)*1.0e9,2)
P_ETP=100
# Calculate Lag and max value time
IIa=np.asarray(pt_IIa_100["04:IIa"])*1.0e9
n_max = IIa.argmax()
t_max_100 = round(t[n_max],2)
P_t_max_100=100
IIa_max_100 = round(IIa[n_max],2)
P_IIa_max_100=100
# Calculate LT
IIa=np.asarray(pt_IIa_100["04:IIa"])*1.0e9
n_10nm = np.extract(IIa[0:n_max]<10.0,IIa[0:n_max]).argmax()+1
t_LT_100 = round(t[n_10nm],2)
print('t_lt_10nm=',t_LT_100)
P_t_LT_100=100

readFile = open("Hockin_modified_FVIII/FVIII_0_pt1_IIa.dat", 'r')
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
    pt_IIa_0[columns[ic+1]]= t1
    ic = ic+1
# Calculate ETP area under the curve
t=pt_IIa_0["01:time"]
dt=t[1]
dt1=(t[2]-t[1])/60
ETP_0=round(simps(pt_IIa_0["04:IIa"],dx=dt1)*1.0e9,2)
P_ETP_0=round((ETP_0/ETP)*100,1)
# Calculate Lag and max value time
IIa=np.asarray(pt_IIa_0["04:IIa"])*1.0e9
n_max = IIa.argmax()
t_max_0 = round(t[n_max],2)
P_t_max_0=round((t_max_0/t_max_100)*100,2)
IIa_max_0 = round(IIa[n_max],2)
P_IIa_max_0=round((IIa_max_0/IIa_max_100)*100,2)
# Calculate LT
IIa=np.asarray(pt_IIa_0["04:IIa"])*1.0e9
n_10nm = np.extract(IIa[0:n_max]<10.0,IIa[0:n_max]).argmax()+1
t_LT_0 = round(t[n_10nm],2)
print('t_lt_10nm=',t_LT_0)
P_t_LT_0=round((t_LT_0/t_LT_100)*100,1)

readFile = open("Hockin_modified_FVIII/FVIII_1_pt1_IIa.dat", 'r')
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
    pt_IIa_1[columns[ic+1]]= t1
    ic = ic+1
# Calculate ETP area under the curve
t=pt_IIa_1["01:time"]
dt=t[1]
dt1=(t[2]-t[1])/60
ETP_1=round(simps(pt_IIa_1["04:IIa"],dx=dt1)*1.0e9,2)
P_ETP_1=round((ETP_1/ETP)*100,1)
# Calculate Lag and max value time
IIa=np.asarray(pt_IIa_1["04:IIa"])*1.0e9
n_max = IIa.argmax()
t_max_1 = round(t[n_max],2)
P_t_max_1=round((t_max_1/t_max_100)*100,2)
IIa_max_1 = round(IIa[n_max],2)
P_IIa_max_1=round((IIa_max_1/IIa_max_100)*100,2)
# Calculate LT
IIa=np.asarray(pt_IIa_1["04:IIa"])*1.0e9
n_10nm = np.extract(IIa[0:n_max]<10.0,IIa[0:n_max]).argmax()+1
t_LT_1 = round(t[n_10nm],2)
print('t_lt_10nm=',t_LT_1)
P_t_LT_1=round((t_LT_1/t_LT_100)*100,1)


readFile = open("Hockin_modified_FVIII/FVIII_5_pt1_IIa.dat", 'r')
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
    pt_IIa_5[columns[ic+1]]= t1
    ic = ic+1
# Calculate ETP area under the curve
t=pt_IIa_5["01:time"]
dt=t[1]
dt1=(t[2]-t[1])/60
ETP_5=round(simps(pt_IIa_5["04:IIa"],dx=dt1)*1.0e9,2)
P_ETP_5=round((ETP_5/ETP)*100,1)
# Calculate Lag and max value time
IIa=np.asarray(pt_IIa_5["04:IIa"])*1.0e9
n_max = IIa.argmax()
t_max_5 = round(t[n_max],2)
P_t_max_5=round((t_max_5/t_max_100)*100,2)
IIa_max_5 = round(IIa[n_max],2)
P_IIa_max_5=round((IIa_max_5/IIa_max_100)*100,2)
# Calculate LT
IIa=np.asarray(pt_IIa_5["04:IIa"])*1.0e9
n_10nm = np.extract(IIa[0:n_max]<10.0,IIa[0:n_max]).argmax()+1
t_LT_5 = round(t[n_10nm],2)
print('t_lt_10nm=',t_LT_5)
P_t_LT_5=round((t_LT_5/t_LT_100)*100,1)


readFile = open("Hockin_modified_FVIII/FVIII_15_pt1_IIa.dat", 'r')
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
    pt_IIa_15[columns[ic+1]]= t1
    ic = ic+1
# Calculate ETP area under the curve
t=pt_IIa_15["01:time"]
dt=t[1]
dt1=(t[2]-t[1])/60
ETP_15=round(simps(pt_IIa_15["04:IIa"],dx=dt1)*1.0e9,2)
P_ETP_15=round((ETP_15/ETP)*100,1)
# Calculate Lag and max value time
IIa=np.asarray(pt_IIa_15["04:IIa"])*1.0e9
n_max = IIa.argmax()
t_max_15 = round(t[n_max],2)
P_t_max_15=round((t_max_15/t_max_100)*100,2)
IIa_max_15 = round(IIa[n_max],2)
P_IIa_max_15=round((IIa_max_15/IIa_max_100)*100,2)
# Calculate LT
IIa=np.asarray(pt_IIa_15["04:IIa"])*1.0e9
n_10nm = np.extract(IIa[0:n_max]<10.0,IIa[0:n_max]).argmax()+1
t_LT_15 = round(t[n_10nm],2)
print('t_lt_10nm=',t_LT_15)
P_t_LT_15=round((t_LT_15/t_LT_100)*100,1)


readFile = open("Hockin_modified_FVIII/FVIII_50_pt1_IIa.dat", 'r')
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
    pt_IIa_50[columns[ic+1]]= t1
    ic = ic+1
# Calculate ETP area under the curve
t=pt_IIa_50["01:time"]
dt=t[1]
dt1=(t[2]-t[1])/60
ETP_50=round(simps(pt_IIa_50["04:IIa"],dx=dt1)*1.0e9,2)
P_ETP_50=round((ETP_50/ETP)*100,1)
# Calculate Lag and max value time
IIa=np.asarray(pt_IIa_50["04:IIa"])*1.0e9
n_max = IIa.argmax()
t_max_50 = round(t[n_max],2)
P_t_max_50=round((t_max_50/t_max_100)*100,2)
IIa_max_50 = round(IIa[n_max],2)
P_IIa_max_50=round((IIa_max_50/IIa_max_100)*100,2)
# Calculate LT
IIa=np.asarray(pt_IIa_50["04:IIa"])*1.0e9
n_10nm = np.extract(IIa[0:n_max]<10.0,IIa[0:n_max]).argmax()+1
t_LT_50 = round(t[n_10nm],2)
print('t_lt_10nm=',t_LT_50)
P_t_LT_50=round((t_LT_50/t_LT_100)*100,1)


alpha_0 = np.asarray(pt_IIa_0["04:IIa"])*1e+9  #+ 1.2*np.asarray(pt_mIIa_0["04:mIIa"])
alpha_1 = np.asarray(pt_IIa_1["04:IIa"])*1e+9 #+ 1.2*np.asarray(pt_mIIa_1["04:mIIa"])
alpha_5 = np.asarray(pt_IIa_5["04:IIa"])*1e+9 #+ 1.2*np.asarray(pt_mIIa_5["04:mIIa"])
alpha_15 = np.asarray(pt_IIa_15["04:IIa"])*1e+9 #+ 1.2*np.asarray(pt_mIIa_15["04:mIIa"])
alpha_50 = np.asarray(pt_IIa_50["04:IIa"])*1e+9 #+ 1.2*np.asarray(pt_mIIa_50["04:mIIa"])
alpha_100 = np.asarray(pt_IIa_100["04:IIa"])*1e+9 #+ 1.2*np.asarray(pt_mIIa_100["04:mIIa"])

###########################################################################################
##### Initial Hockin ##################################################################
###########################################################################################

readFile = open("FVIII_Butenas_100_pt1_IIa.dat", 'r')
#readFile = open("FVIII_modified_100_pt1_IIa.dat", 'r')
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
    pt_Initial_IIa_100[columns[ic+1]]= t1
    ic = ic+1
t=pt_Initial_IIa_100["01:time"]
dt=t[1]
dt1=(t[2]-t[1])/60
ETP_100=round(simps(pt_Initial_IIa_100["04:IIa"],dx=dt1)*1.0e9,2)
Initial_P_ETP_100=100.0
# Calculate Lag and max value time
IIa=np.asarray(pt_Initial_IIa_100["04:IIa"])*1.0e9
n_max = IIa.argmax()
t_max_100 = round(t[n_max],2)
Initial_P_t_max_100=100.0
IIa_max_100 = round(IIa[n_max],2)
Initial_P_IIa_max_100=100.0


readFile = open("FVIII_Butenas_0_pt1_IIa.dat", 'r')
#readFile = open("FVIII_modified_0_pt1_IIa.dat", 'r')
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
    pt_Initial_IIa_0[columns[ic+1]]= t1
    ic = ic+1
t=pt_Initial_IIa_0["01:time"]
dt=t[1]
dt1=(t[2]-t[1])/60
ETP_0=round(simps(pt_Initial_IIa_0["04:IIa"],dx=dt1)*1.0e9,2)
Initial_P_ETP_0=round((ETP_0/ETP)*100,1)
# Calculate Lag and max value time
IIa=np.asarray(pt_Initial_IIa_0["04:IIa"])*1.0e9
n_max = IIa.argmax()
t_max_0 = round(t[n_max],2)
Initial_P_t_max_0=round((t_max_0/t_max_100)*100,2)
IIa_max_0 = round(IIa[n_max],2)
Initial_P_IIa_max_0=round((IIa_max_0/IIa_max_100)*100,2)


readFile = open("FVIII_Butenas_1_pt1_IIa.dat", 'r')
#readFile = open("FVIII_modified_1_pt1_IIa.dat", 'r')
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
    pt_Initial_IIa_1[columns[ic+1]]= t1
    ic = ic+1
t=pt_Initial_IIa_1["01:time"]
dt=t[1]
dt1=(t[2]-t[1])/60
ETP_1=round(simps(pt_Initial_IIa_1["04:IIa"],dx=dt1)*1.0e9,2)
Initial_P_ETP_1=round((ETP_1/ETP)*100,1)
# Calculate Lag and max value time
IIa=np.asarray(pt_Initial_IIa_1["04:IIa"])*1.0e9
n_max = IIa.argmax()
t_max_1 = round(t[n_max],2)
Initial_P_t_max_1=round((t_max_1/t_max_100)*100,2)
IIa_max_1 = round(IIa[n_max],2)
Initial_P_IIa_max_1=round((IIa_max_1/IIa_max_100)*100,2)

readFile = open("FVIII_Butenas_5_pt1_IIa.dat", 'r')
#readFile = open("FVIII_modified_5_pt1_IIa.dat", 'r')
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
    pt_Initial_IIa_5[columns[ic+1]]= t1
    ic = ic+1
t=pt_Initial_IIa_5["01:time"]
dt=t[1]
dt1=(t[2]-t[1])/60
ETP_5=round(simps(pt_Initial_IIa_5["04:IIa"],dx=dt1)*1.0e9,2)
Initial_P_ETP_5=round((ETP_5/ETP)*100,1)
# Calculate Lag and max value time
IIa=np.asarray(pt_Initial_IIa_5["04:IIa"])*1.0e9
n_max = IIa.argmax()
t_max_5 = round(t[n_max],2)
Initial_P_t_max_5=round((t_max_5/t_max_100)*100,2)
IIa_max_5 = round(IIa[n_max],2)
Initial_P_IIa_max_5=round((IIa_max_5/IIa_max_100)*100,2)

readFile = open("FVIII_Butenas_15_pt1_IIa.dat", 'r')
#readFile = open("FVIII_modified_15_pt1_IIa.dat", 'r')
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
    pt_Initial_IIa_15[columns[ic+1]]= t1
    ic = ic+1
t=pt_Initial_IIa_15["01:time"]
dt=t[1]
dt1=(t[2]-t[1])/60
ETP_15=round(simps(pt_Initial_IIa_15["04:IIa"],dx=dt1)*1.0e9,2)
Initial_P_ETP_15=round((ETP_15/ETP)*100,1)
# Calculate Lag and max value time
IIa=np.asarray(pt_Initial_IIa_15["04:IIa"])*1.0e9
n_max = IIa.argmax()
t_max_15 = round(t[n_max],2)
Initial_P_t_max_15=round((t_max_15/t_max_100)*100,2)
IIa_max_15 = round(IIa[n_max],2)
Initial_P_IIa_max_15=round((IIa_max_15/IIa_max_100)*100,2)

readFile = open("FVIII_Butenas_50_pt1_IIa.dat", 'r')
#readFile = open("FVIII_modified_50_pt1_IIa.dat", 'r')
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
    pt_Initial_IIa_50[columns[ic+1]]= t1
    ic = ic+1
t=pt_Initial_IIa_50["01:time"]
dt=t[1]
dt1=(t[2]-t[1])/60
ETP_50=round(simps(pt_Initial_IIa_50["04:IIa"],dx=dt1)*1.0e9,2)
Initial_P_ETP_50=round((ETP_50/ETP)*100,1)
# Calculate Lag and max value time
IIa=np.asarray(pt_Initial_IIa_50["04:IIa"])*1.0e9
n_max = IIa.argmax()
t_max_50 = round(t[n_max],2)
Initial_P_t_max_50 =round((t_max_50/t_max_100)*100,2)
IIa_max_50 = round(IIa[n_max],2)
Initial_P_IIa_max_50=round((IIa_max_50/IIa_max_100)*100,2)

Initial_alpha_0 = np.asarray(pt_Initial_IIa_0["04:IIa"])*1e+9  #+ 1.2*np.asarray(pt_mIIa_0["04:mIIa"])
Initial_alpha_1 = np.asarray(pt_Initial_IIa_1["04:IIa"])*1e+9 #+ 1.2*np.asarray(pt_mIIa_1["04:mIIa"])
Initial_alpha_5 = np.asarray(pt_Initial_IIa_5["04:IIa"])*1e+9 #+ 1.2*np.asarray(pt_mIIa_5["04:mIIa"])
Initial_alpha_15 = np.asarray(pt_Initial_IIa_15["04:IIa"])*1e+9 #+ 1.2*np.asarray(pt_mIIa_15["04:mIIa"])
Initial_alpha_50 = np.asarray(pt_Initial_IIa_50["04:IIa"])*1e+9 #+ 1.2*np.asarray(pt_mIIa_50["04:mIIa"])
Initial_alpha_100 = np.asarray(pt_Initial_IIa_100["04:IIa"])*1e+9 #+ 1.2*np.asarray(pt_mIIa_100["04:mIIa"])
time_exp = np.asarray(expe["time"])*60.0


table_vals=[[t_LT_100,P_t_LT_100],[t_LT_50,P_t_LT_50],[t_LT_15,P_t_LT_15],[t_LT_5,P_t_LT_5],[t_LT_1,P_t_LT_1],[t_LT_0,P_t_LT_0]]


col_labels=['LT s','LT %']
row_labels=['$FVIII_{100 \%}$','$FVIII_{50 \%}$','$FVIII_{15 \%}$','$FVIII_{5 \%}$','$FVIII_{1 \%}$','$FVIII_{0 \%}$']

expe_col_labels=['LT s','LT %']
expe_row_labels=['$FVIII_{100 \%}$','$FVIII_{50 \%}$','$FVIII_{15 \%}$','$FVIII_{5 \%}$','$FVIII_{1 \%}$','$FVIII_{0 \%}$']

dt_expe=(time_exp[2]-time_exp[1])/60.0

#100
expe_ETP=round(simps(expe["FVIII 100%"],dx=dt_expe),2)
expe_P_ETP=100
expe_IIa=np.asarray(expe["FVIII 100%"])
expe_n_max_100=expe_IIa.argmax()
expe_t_max_100 = round(time_exp[expe_n_max_100],2)
expe_P_t_max_100=100
expe_IIa_max_100=round(expe_IIa[expe_n_max_100],2)
expe_P_IIa_max_100=100



expe_n_10nm = np.extract(expe_IIa[0:expe_n_max_100]<10.0,expe_IIa[0:expe_n_max_100]).argmax()+1
expe_t_LT_100 = round(time_exp[expe_n_10nm],2)
expe_P_t_LT_100=100

print('expe_IIA(n_10)=',expe_IIa[expe_n_10nm])
print('t(n_10)=',time_exp[expe_n_10nm])

#50
expe_ETP_50=round(simps(expe["FVIII 50%"],dx=dt_expe),2)
expe_P_ETP_50=round((expe_ETP_50/expe_ETP)*100,2)
expe_IIa=np.asarray(expe["FVIII 50%"])
expe_n_max_50=expe_IIa.argmax()
expe_t_max_50 = round(time_exp[expe_n_max_50],2)
expe_P_t_max_50=round((expe_t_max_50/expe_t_max_100)*100,2)
expe_IIa_max_50=round(expe_IIa[expe_n_max_50],2)
expe_P_IIa_max_50=round((expe_IIa_max_50/expe_IIa_max_100)*100,2)

expe_n_10nm = np.extract(expe_IIa[0:expe_n_max_50]<10.0,expe_IIa[0:expe_n_max_50]).argmax()+1
expe_t_LT_50 = round(time_exp[expe_n_10nm],2)
expe_P_t_LT_50=round((expe_t_LT_50/expe_t_LT_100)*100,1)

#15
expe_ETP_15=round(simps(expe["FVIII 15%"],dx=dt_expe),2)
expe_P_ETP_15=round((expe_ETP_15/expe_ETP)*100,2)
expe_IIa=np.asarray(expe["FVIII 15%"])
expe_n_max_15=expe_IIa.argmax()
expe_t_max_15 = round(time_exp[expe_n_max_15],2)
expe_P_t_max_15=round((expe_t_max_15/expe_t_max_100)*100,2)
expe_IIa_max_15=round(expe_IIa[expe_n_max_15],2)
expe_P_IIa_max_15=round((expe_IIa_max_15/expe_IIa_max_100)*100,2)

expe_n_10nm = np.extract(expe_IIa[0:expe_n_max_15]<10.0,expe_IIa[0:expe_n_max_15]).argmax()+1
expe_t_LT_15 = round(time_exp[expe_n_10nm],2)
expe_P_t_LT_15=round((expe_t_LT_15/expe_t_LT_100)*100,1)

#5
expe_ETP_5=round(simps(expe["FVIII 5%"],dx=dt_expe),2)
expe_P_ETP_5=round((expe_ETP_5/expe_ETP)*100,2)
expe_IIa=np.asarray(expe["FVIII 5%"])
expe_n_max_5=expe_IIa.argmax()
expe_t_max_5 = round(time_exp[expe_n_max_5],2)
expe_P_t_max_5=round((expe_t_max_5/expe_t_max_100)*100,2)
expe_IIa_max_5=round(expe_IIa[expe_n_max_5],2)
expe_P_IIa_max_5=round((expe_IIa_max_5/expe_IIa_max_100)*100,2)

expe_n_10nm = np.extract(expe_IIa[0:expe_n_max_5]<10.0,expe_IIa[0:expe_n_max_5]).argmax()+1
expe_t_LT_5 = round(time_exp[expe_n_10nm],2)
expe_P_t_LT_5=round((expe_t_LT_5/expe_t_LT_100)*100,1)

#1
expe_ETP_1=round(simps(expe["FVIII 1%"],dx=dt_expe),2)
expe_P_ETP_1=round((expe_ETP_1/expe_ETP)*100,2)
expe_IIa=np.asarray(expe["FVIII 1%"])
expe_n_max_1=expe_IIa.argmax()
expe_t_max_1 = round(time_exp[expe_n_max_1],2)
expe_P_t_max_1=round((expe_t_max_1/expe_t_max_100)*100,2)
expe_IIa_max_1=round(expe_IIa[expe_n_max_1],2)
expe_P_IIa_max_1=round((expe_IIa_max_1/expe_IIa_max_100)*100,2)

expe_n_10nm = np.extract(expe_IIa[0:expe_n_max_1]<10.0,expe_IIa[0:expe_n_max_1]).argmax()+1
expe_t_LT_1 = round(time_exp[expe_n_10nm],2)
expe_P_t_LT_1=round((expe_t_LT_1/expe_t_LT_100)*100,1)

#0
expe_ETP_0=round(simps(expe["FVIII 0%"],dx=dt_expe),2)
expe_P_ETP_0=round((expe_ETP_0/expe_ETP)*100,2)
expe_IIa=np.asarray(expe["FVIII 0%"])
expe_n_max_0=expe_IIa.argmax()
expe_t_max_0 = round(time_exp[expe_n_max_0],2)
expe_P_t_max_0=round((expe_t_max_0/expe_t_max_100)*100,2)
expe_IIa_max_0=round(expe_IIa[expe_n_max_0],2)
expe_P_IIa_max_0=round((expe_IIa_max_0/expe_IIa_max_100)*100,2)
expe_n_10nm = np.extract(expe_IIa[0:expe_n_max_0]<10.0,expe_IIa[0:expe_n_max_0]).argmax()+1
expe_t_LT_0 = round(time_exp[expe_n_10nm],2)
expe_P_t_LT_0=round((expe_t_LT_0/expe_t_LT_100)*100,1)


expe_table_vals=[[expe_t_LT_100,expe_P_t_LT_100],[expe_t_LT_50,expe_P_t_LT_50],[expe_t_LT_15,expe_P_t_LT_15],[expe_t_LT_5,expe_P_t_LT_5],[expe_t_LT_1,expe_P_t_LT_1],[expe_t_LT_0,expe_P_t_LT_0]]

#ETPFig_7=[["FVIII_(%)","Experiments","Original_Hockin","Modified_Hockin"],
ETPFig_7=[
        [100,expe_P_ETP,Initial_P_ETP_100,P_ETP],
        [50,expe_P_ETP_50,Initial_P_ETP_50,P_ETP_50],
        [15,expe_P_ETP_15,Initial_P_ETP_15,P_ETP_15],
        [5,expe_P_ETP_5,Initial_P_ETP_5,P_ETP_5],
        [1,expe_P_ETP_1,Initial_P_ETP_1,P_ETP_1],
        [0,expe_P_ETP_0,Initial_P_ETP_0,P_ETP_0]]
IIa_max_Fig_7=[
        [100,expe_P_IIa_max_100,Initial_P_IIa_max_100,P_IIa_max_100],
        [50,expe_P_IIa_max_50,Initial_P_IIa_max_50,P_IIa_max_50],
        [15,expe_P_IIa_max_15,Initial_P_IIa_max_15,P_IIa_max_15],
        [5,expe_P_IIa_max_5,Initial_P_IIa_max_5,P_IIa_max_5],
        [1,expe_P_IIa_max_1,Initial_P_IIa_max_1,P_IIa_max_1],
        [0,expe_P_IIa_max_0,Initial_P_IIa_max_0,P_IIa_max_0]]
t_max_Fig_7=[
        [100,expe_P_t_max_100,Initial_P_t_max_100,P_t_max_100],
        [50,expe_P_t_max_50,Initial_P_t_max_50,P_t_max_50],
        [15,expe_P_t_max_15,Initial_P_t_max_15,P_t_max_15],
        [5,expe_P_t_max_5,Initial_P_t_max_5,P_t_max_5],
        [1,expe_P_t_max_1,Initial_P_t_max_1,P_t_max_1],
        [0,expe_P_t_max_0,Initial_P_t_max_0,P_t_max_0]]
ETPFig_7=np.asarray(ETPFig_7)
IIa_max_Fig_7=np.asarray(IIa_max_Fig_7)
print(IIa_max_Fig_7)
t_max_Fig_7=np.asarray(t_max_Fig_7)
np.savetxt('FVIII_ETP.csv',ETPFig_7,delimiter=",")
np.savetxt('FVIII_Peak.csv',IIa_max_Fig_7,delimiter=",")
np.savetxt('FVIII_Ttpeak.csv',t_max_Fig_7,delimiter=",")


font = {'family' : 'serif',
        #'weight' : 'bold',
        'size'   : 25}
#plt.rc('text', usetex=True)
plt.rc('font', **font)

#plt.rcParams.update({'font.size': 25})

a1=plt.subplot(3,1,1)

plt.title('Original Ext')
plt.plot(pt_Initial_IIa_100["01:time"],Initial_alpha_100,'k-',lw=2.5, label= r'$FVIII_{100 \%}$',markevery=500,markersize=15,markerfacecolor='none')
plt.plot(pt_Initial_IIa_50["01:time"],Initial_alpha_50,'k--',lw=2.5, label= r'$FVIII_{50 \%}$',markevery=500,markersize=15,markerfacecolor='none')
plt.plot(pt_Initial_IIa_15["01:time"],Initial_alpha_15,'k-.',lw=2.5, label= r'$FVIII_{15 \%}$',markevery=500,markersize=15,markerfacecolor='none')
plt.plot(pt_Initial_IIa_5["01:time"],Initial_alpha_5,'k:',lw=2.5, label= r'$FVIII_{5 \%}$' ,markevery=50,markersize=15,markerfacecolor='none')
plt.plot(pt_Initial_IIa_1["01:time"],Initial_alpha_1,'kx',lw=2.5, label= r'$FVIII_{1 \%}$' ,markevery=350,markersize=5)
plt.plot(pt_Initial_IIa_0["01:time"],Initial_alpha_0,'ko',lw=2.5, label= r'$FVIII_{0 \%}$' ,markevery=350,markersize=5)
plt.ylabel(r'$IIa$ nM ')
plt.xlim((0,2500))
#plt.legend(loc='upper right', prop={'size':12})
#plt.table(cellText=table_vals,
#                   colWidths = [0.1]*5,
#                   rowLabels=row_labels,
#                   colLabels=col_labels,
#                   loc='center right')

plt.setp(a1.get_xticklabels(), visible=False)

a2=plt.subplot(3,1,2)
plt.title('Modified Ext')
plt.plot(pt_IIa_100["01:time"],alpha_100,'k-',lw=2.5, label= r'$FVIII_{100 \%}$',markevery=500,markersize=15,markerfacecolor='none')
plt.plot(pt_IIa_50["01:time"],alpha_50,'k--',lw=2.5, label= r'$FVIII_{50 \%}$',markevery=500,markersize=15,markerfacecolor='none')
plt.plot(pt_IIa_15["01:time"],alpha_15,'k-.',lw=2.5, label= r'$FVIII_{15 \%}$',markevery=500,markersize=15,markerfacecolor='none')
plt.plot(pt_IIa_5["01:time"],alpha_5,'k:',lw=2.5, label= r'$FVIII_{5 \%}$' ,markevery=500,markersize=15,markerfacecolor='none')
plt.plot(pt_IIa_1["01:time"],alpha_1,'kx',lw=2.5, label= r'$FVIII_{1 \%}$' ,markevery=350,markersize=5)
plt.plot(pt_IIa_0["01:time"],alpha_0,'ko',lw=2.5, label= r'$FVIII_{0 \%}$' ,markevery=350,markersize=5)
plt.ylabel(r'$IIa$ nM ')
#plt.xlabel(r'$t$ s')
plt.xlim((0,2500))
#plt.ylim((0,1.5E-7))
plt.setp(a2.get_xticklabels(), visible=False)

a3=plt.subplot(3,1,3)
plt.title('Experiment')
plt.plot(time_exp,expe["FVIII 100%"],'k-',lw=2.5,label=r'$FVIII_{100\%}$', markevery=500)
plt.plot(time_exp,expe["FVIII 50%"],'k--',lw=2.5,label=r'$FVIII_{50\%} $', markevery=500)
plt.plot(time_exp,expe["FVIII 15%"],'k-.',lw=2.5,label=r'$FVIII_{15\%}$ ', markevery=500)
#plt.plot(time_exp,expe["XV 10%"],'k-.',lw=2.5,label=r'$FVIII_{10\%}$ ', markevery=500)
plt.plot(time_exp,expe["FVIII 5%"],'k:',lw=2.5,label=r'$FVIII_{5\%}$ ', markevery=5)
plt.plot(time_exp,expe["FVIII 1%"],'kx',lw=2.5,label=r'$FVIII_{1\%}$ ', markevery=3)
plt.plot(time_exp,expe["FVIII 0%"],'ko',lw=2.5,label=r'$FVIII_{0\%}$ ', markevery=3)
#plt.table(cellText=expe_table_vals,
#                  colWidths = [0.14,0.14],
#                  rowLabels=row_labels,
#                  colLabels=col_labels,
#                  loc='upper right')
plt.ylabel(r'$IIa$ nM ')
plt.xlabel(r'$t$ s ')
plt.xlim((0,2500))
#plt.ylim((0,150))

#import tikzplotlib
#
#tikzplotlib.save("FVIII.tex")

#plt.savefig('F1a.pgf')

plt.show()
