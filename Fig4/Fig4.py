#!/usr/bin/python
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.integrate import simps
import os 
import csv

DictCases={}
DictVariables={}
path='FVIII'
Files=os.listdir(path)
cases=np.array([0,1,100,15,5,50])
for count, file in enumerate(Files):
    MainDict={}
    Variables={}
    #Read Data
    Filestr=path+'/'+file
    print(Filestr)
    readFile = open(Filestr, 'r')
    sepFile = readFile.read().split('\n')
    readFile.close()
    DataLine= sepFile[3:].index("# Simulation details: -------------------------------------------")+5
    columns = sepFile[DataLine].split() 
    Nspecies = len(columns)-1
    ic=0
    for x in range(0,Nspecies):
        data = []
        for linedata in sepFile[155:]:
            vector = linedata.split()
            if len(vector) > 1:
                data.append(float(vector[ic]))
        MainDict[columns[ic+1]]=data
        ic = ic+1
    print(cases[count])
    DictCases[cases[count]]=MainDict
    # Calculate ETP taumax
    t=MainDict["time"]
    ETP=round(simps(MainDict["IIa"],x=t),2)
    # Calculate IIa_max and tau_max value time
    IIa=np.asarray(MainDict["IIa"])
    n_max = IIa.argmax()
    t_max_ = round(t[n_max],2)
    IIa_max = round(IIa[n_max],2)
    Variables['ETP']=ETP
    Variables['IIaMax']=IIa_max
    Variables['tauMax']=t_max_
    # Save Variables to dictionary
    DictVariables[cases[count]]=Variables

# Write ETP, tau_max and IIa_max to csv
w = csv.writer(open(path+".csv", "w"))
for key, val in DictVariables.items():
    w.writerow([key, val])

###########################################################################################
##### Initial Chatterjee ##################################################################
###########################################################################################

pt_Initial_IIa_0 = {}
pt_Initial_IIa_1 = {}
pt_Initial_IIa_5 = {}
pt_Initial_IIa_15 = {}
pt_Initial_IIa_50 = {}
pt_Initial_IIa_100 = {}

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
    pt_Initial_IIa_100[columns[ic+1]]= t1
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
    pt_Initial_IIa_0[columns[ic+1]]= t1
    ic = ic+1

readFile = open("Chatterjee_initial_FVIII/FVIII_1_pt1_IIa.dat", 'r')
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

readFile = open("Chatterjee_initial_FVIII/FVIII_5_pt1_IIa.dat", 'r')
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
    pt_Initial_IIa_15[columns[ic+1]]= t1
    ic = ic+1

readFile = open("Chatterjee_initial_FVIII/FVIII_50_pt1_IIa.dat", 'r')
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

Initial_alpha_0 = np.asarray(pt_Initial_IIa_0["04:IIa"])*1e+9  #+ 1.2*np.asarray(pt_mIIa_0["04:mIIa"])
Initial_alpha_1 = np.asarray(pt_Initial_IIa_1["04:IIa"])*1e+9 #+ 1.2*np.asarray(pt_mIIa_1["04:mIIa"])
Initial_alpha_5 = np.asarray(pt_Initial_IIa_5["04:IIa"])*1e+9 #+ 1.2*np.asarray(pt_mIIa_5["04:mIIa"])
Initial_alpha_15 = np.asarray(pt_Initial_IIa_15["04:IIa"])*1e+9 #+ 1.2*np.asarray(pt_mIIa_15["04:mIIa"])
Initial_alpha_50 = np.asarray(pt_Initial_IIa_50["04:IIa"])*1e+9 #+ 1.2*np.asarray(pt_mIIa_50["04:mIIa"])
Initial_alpha_100 = np.asarray(pt_Initial_IIa_100["04:IIa"])*1e+9 #+ 1.2*np.asarray(pt_mIIa_100["04:mIIa"])

# Calculate ETP area under the curve 100
t=pt_Initial_IIa_100["01:time"]
dt=t[1]
dt1=(t[2]-t[1])/60
ETP_100=round(simps(pt_Initial_IIa_100["04:IIa"],dx=dt1)*1.0e9,2)
# Calculate IIa max value and tau_max
IIa=np.asarray(pt_Initial_IIa_100["04:IIa"])*1.0e9
n_max = IIa.argmax()
t_max_100 = round(t[n_max],2)
IIa_max_100 = round(IIa[n_max],2)
# Calculate ETP area under the curve 50
t=pt_Initial_IIa_50["01:time"]
dt=t[1]
dt1=(t[2]-t[1])/60
ETP_50=round(simps(pt_Initial_IIa_50["04:IIa"],dx=dt1)*1.0e9,2)
# Calculate IIa max value and tau_max
IIa=np.asarray(pt_Initial_IIa_50["04:IIa"])*1.0e9
n_max = IIa.argmax()
t_max_50 = round(t[n_max],2)
IIa_max_50 = round(IIa[n_max],2)
# Calculate ETP area under the curve 15
t=pt_Initial_IIa_15["01:time"]
dt=t[1]
dt1=(t[2]-t[1])/60
ETP_15=round(simps(pt_Initial_IIa_15["04:IIa"],dx=dt1)*1.0e9,2)
# Calculate IIa max value and tau_max
IIa=np.asarray(pt_Initial_IIa_15["04:IIa"])*1.0e9
n_max = IIa.argmax()
t_max_15 = round(t[n_max],2)
IIa_max_15 = round(IIa[n_max],2)
# Calculate ETP area under the curve 5
t=pt_Initial_IIa_5["01:time"]
dt=t[1]
dt1=(t[2]-t[1])/60
ETP_5=round(simps(pt_Initial_IIa_5["04:IIa"],dx=dt1)*1.0e9,2)
# Calculate IIa max value and tau_max
IIa=np.asarray(pt_Initial_IIa_5["04:IIa"])*1.0e9
n_max = IIa.argmax()
t_max_5 = round(t[n_max],2)
IIa_max_5 = round(IIa[n_max],2)
# Calculate ETP area under the curve 1
t=pt_Initial_IIa_1["01:time"]
dt=t[1]
dt1=(t[2]-t[1])/60
ETP_1=round(simps(pt_Initial_IIa_1["04:IIa"],dx=dt1)*1.0e9,2)
# Calculate IIa max value and tau_max
IIa=np.asarray(pt_Initial_IIa_1["04:IIa"])*1.0e9
n_max = IIa.argmax()
t_max_1 = round(t[n_max],2)
IIa_max_1 = round(IIa[n_max],2)
# Calculate ETP area under the curve 0
t=pt_Initial_IIa_0["01:time"]
dt=t[1]
dt0=(t[2]-t[1])/60
ETP_0=round(simps(pt_Initial_IIa_0["04:IIa"],dx=dt1)*1.0e9,2)
# Calculate IIa max value and tau_max
IIa=np.asarray(pt_Initial_IIa_0["04:IIa"])*1.0e9
n_max = IIa.argmax()
t_max_0 = round(t[n_max],2)
IIa_max_0 = round(IIa[n_max],2)

InitialFVIII=np.array([[ETP_0,IIa_max_0,t_max_0],[ETP_1,IIa_max_1,t_max_1],[ETP_5,IIa_max_5,t_max_5],[ETP_15,IIa_max_15,t_max_15],[ETP_50,IIa_max_50,t_max_50],[ETP_100,IIa_max_100,t_max_100]])
print(InitialFVIII)

#Read Experimental Data
expe = {}
readFile = open(path+"_Contact_Donnees_TGT.csv", 'r')
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
time_exp = np.asarray(expe["time"])*60.0

##Compute ETP, IIa_max and tau_max from experimental data
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

expe_table_valsP=[['100%',expe_P_ETP,expe_P_IIa_max_100,expe_P_t_max_100],
                ['50%',expe_P_ETP_50,expe_P_IIa_max_50,expe_P_t_max_50],
                ['15%',expe_P_ETP_50,expe_P_IIa_max_50,expe_P_t_max_50],
                ['5%',expe_P_ETP_5,expe_P_IIa_max_5,expe_P_t_max_5],
                ['1%',expe_P_ETP_1,expe_P_IIa_max_1,expe_P_t_max_1],
                ['0%',expe_P_ETP_0,expe_P_IIa_max_0,expe_P_t_max_0],
            ]
print('Experimental ETP, IIa_max and tau_max %')
print(expe_table_valsP)

expe_table_vals=[['100%',expe_ETP,expe_IIa_max_100,expe_t_max_100],
                ['50%',expe_ETP_50,expe_IIa_max_50,expe_t_max_50],
                ['15%',expe_ETP_50,expe_IIa_max_50,expe_t_max_50],
                ['5%',expe_ETP_5,expe_IIa_max_5,expe_t_max_5],
                ['1%',expe_ETP_1,expe_IIa_max_1,expe_t_max_1],
                ['0%',expe_ETP_0,expe_IIa_max_0,expe_t_max_0],
            ]
print('Experimental ETP, IIa_max and tau_max')
print(expe_table_vals)

font = {'family' : 'serif',
        'size'   : 25}
plt.rc('font', **font)

plt.title('100 %')
plt.ylabel(r'$II_a$ (nM) ')
plt.xlim((0,1000))
plt.ylim((-30,700))

a1=plt.subplot(3,2,1)
plt.title('100 %')
plt.plot(pt_Initial_IIa_100["01:time"],Initial_alpha_100,'k-',lw=2.5, label= r'$FVIII_{100 \%}$',markevery=500,markersize=15,markerfacecolor='none')
plt.plot(DictCases[100]["time"],DictCases[100]['IIa'],'k--',lw=2.5, label= r'$FVIII_{100 \%}$',markevery=500,markersize=15,markerfacecolor='none')
plt.plot(time_exp,expe["FVIII 100%"],'ko',lw=2.5,label=r'$FVIII_{100\%}$', markevery=1)
plt.ylabel(r'$II_a$ (nM) ')
plt.xlim((0,1000))
plt.ylim((-30,700))
plt.setp(a1.get_xticklabels(), visible=False)

a2=plt.subplot(3,2,3)
plt.title('50 %')
plt.plot(pt_Initial_IIa_50["01:time"],Initial_alpha_50,'k-',lw=2.5, label= r'$FVIII_{50 \%}$',markevery=500,markersize=15,markerfacecolor='none')
plt.plot(DictCases[50]["time"],DictCases[50]['IIa'],'k--',lw=2.5, label= r'$FVIII_{50 \%}$',markevery=500,markersize=15,markerfacecolor='none')
plt.plot(time_exp,expe["FVIII 50%"],'ko',lw=2.5,label=r'$FVIII_{50\%} $', markevery=1)
plt.ylabel(r'$II_a$ (nM)')
plt.xlim((0,1000))
plt.ylim((-30,700))
plt.setp(a2.get_xticklabels(), visible=False)

a3=plt.subplot(3,2,5)
plt.title('15 %')
plt.plot(pt_Initial_IIa_15["01:time"],Initial_alpha_15,'k-',lw=2.5, label= r'$FVIII_{15 \%}$',markevery=500,markersize=15,markerfacecolor='none')
plt.plot(DictCases[15]["time"],DictCases[15]['IIa'],'k--',lw=2.5, label= r'$FVIII_{15 \%}$',markevery=500,markersize=15,markerfacecolor='none')
plt.plot(time_exp,expe["FVIII 15%"],'ko',lw=2.5,label=r'$FVIII_{15\%}$ ', markevery=1)
plt.ylabel(r'$II_a$ (nM)')
plt.xlabel(r'time (s) ')
plt.xlim((0,1000))
plt.ylim((-30,700))

a4=plt.subplot(3,2,2)
plt.title('5 %')
plt.plot(pt_Initial_IIa_5["01:time"],Initial_alpha_5,'k-',lw=2.5, label= r'$FVIII_{5 \%}$' ,markevery=50,markersize=15,markerfacecolor='none')
plt.plot(DictCases[5]["time"],DictCases[5]['IIa'],'k--',lw=2.5, label= r'$FVIII_{5 \%}$',markevery=500,markersize=15,markerfacecolor='none')
plt.plot(time_exp,expe["FVIII 5%"],'ko',lw=2.5,label=r'$FVIII_{5\%}$ ', markevery=1)
plt.xlim((0,1000))
plt.ylim((-30,700))
plt.setp(a4.get_xticklabels(), visible=False)

a5=plt.subplot(3,2,4)
plt.title('1 %')
plt.plot(pt_Initial_IIa_1["01:time"],Initial_alpha_1,'k-',lw=2.5, label= r'$FVIII_{1 \%}$' ,markevery=350,markersize=5)
plt.plot(DictCases[1]["time"],DictCases[1]['IIa'],'k--',lw=2.5, label= r'$FVIII_{1 \%}$',markevery=500,markersize=15,markerfacecolor='none')
plt.plot(time_exp,expe["FVIII 1%"],'ko',lw=2.5,label=r'$FVIII_{1\%}$ ', markevery=1)
plt.xlim((0,1000))
plt.ylim((-30,700))
plt.setp(a5.get_xticklabels(), visible=False)

a6=plt.subplot(3,2,6)
plt.title('0 %')
plt.plot(pt_Initial_IIa_0["01:time"],Initial_alpha_0,'k-',lw=2.5, label= r'$FVIII_{0 \%}$' ,markevery=350,markersize=5)
plt.plot(DictCases[0]["time"],DictCases[0]['IIa'],'k--',lw=2.5, label= r'$FVIII_{0 \%}$',markevery=500,markersize=15,markerfacecolor='none')
plt.plot(time_exp,expe["FVIII 0%"],'ko',lw=2.5,label=r'$FVIII_{0\%}$ ', markevery=1)
plt.xlabel(r'time (s) ')
plt.xlim((0,1000))
plt.ylim((-30,700))
plt.show()
