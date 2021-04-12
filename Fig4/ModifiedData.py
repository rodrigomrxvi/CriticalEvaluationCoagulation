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
cases=np.array([0,1,5,15,50,100])
for count, file in enumerate(cases):
    MainDict={}
    Variables={}
    #Read Data
    #Filestr=path+'/'+file
    Filestr=path+'/Int_pathway_modified_'+path+str(file)+'.out'
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

font = {'family' : 'serif',
        'size'   : 25}
plt.rc('font', **font)

plt.plot(DictCases[1]["time"],DictCases[1]['IIa'],'k-',lw=2.5, label= r'$FVIII_{1 \%}$',markevery=500,markersize=15,markerfacecolor='none')
#plt.plot(DictCases[15["time"],DictCases[15]['IIa'],'k-.',lw=2.5, label= r'$FVIII_{15 \%}$',markevery=500,markersize=15,markerfacecolor='none')
#plt.plot(DictCases[50["time"],DictCases[50]['IIa'],'k--',lw=2.5, label= r'$FVIII_{50 \%}$',markevery=500,markersize=15,markerfacecolor='none')
#plt.plot(DictCases[100["time"],DictCases[100]['IIa'],'k:',lw=2.5, label= r'$FVIII_{100 \%}$',markevery=500,markersize=15,markerfacecolor='none')
plt.ylabel(r'$II_a$ (nM) ')
plt.xlim((0,1000))
plt.ylim((-30,700))

plt.show()

#a1=plt.subplot(3,2,1)
#plt.title('100 %')
#plt.plot(pt_Initial_IIa_100["01:time"],Initial_alpha_100,'k-',lw=2.5, label= r'$FVIII_{100 \%}$',markevery=500,markersize=15,markerfacecolor='none')
#plt.plot(pt_IIa_100["01:time"],alpha_100,'k--',lw=2.5, label= r'$FVIII_{100 \%}$',markevery=500,markersize=15,markerfacecolor='none')
#plt.plot(time_exp,expe["FVIII 100%"],'ko',lw=2.5,label=r'$FVIII_{100\%}$', markevery=1)
#plt.ylabel(r'$II_a$ (nM) ')
#plt.xlim((0,1000))
#plt.ylim((-30,700))
#plt.setp(a1.get_xticklabels(), visible=False)
#
#a2=plt.subplot(3,2,3)
#plt.title('50 %')
#plt.plot(pt_Initial_IIa_50["01:time"],Initial_alpha_50,'k-',lw=2.5, label= r'$FVIII_{50 \%}$',markevery=500,markersize=15,markerfacecolor='none')
#plt.plot(pt_IIa_50["01:time"],alpha_50,'k--',lw=2.5, label= r'$FVIII_{50 \%}$',markevery=500,markersize=15,markerfacecolor='none')
#plt.plot(time_exp,expe["FVIII 50%"],'ko',lw=2.5,label=r'$FVIII_{50\%} $', markevery=1)
#plt.ylabel(r'$II_a$ (nM)')
#plt.xlim((0,1000))
#plt.ylim((-30,700))
#plt.setp(a2.get_xticklabels(), visible=False)
#
#a3=plt.subplot(3,2,5)
#plt.title('15 %')
#plt.plot(pt_Initial_IIa_15["01:time"],Initial_alpha_15,'k-',lw=2.5, label= r'$FVIII_{15 \%}$',markevery=500,markersize=15,markerfacecolor='none')
#plt.plot(pt_IIa_15["01:time"],alpha_15,'k--',lw=2.5, label= r'$FVIII_{15 \%}$',markevery=500,markersize=15,markerfacecolor='none')
#plt.plot(time_exp,expe["FVIII 15%"],'ko',lw=2.5,label=r'$FVIII_{15\%}$ ', markevery=1)
#plt.ylabel(r'$II_a$ (nM)')
#plt.xlabel(r'time (s) ')
#plt.xlim((0,1000))
#plt.ylim((-30,700))
#
#a4=plt.subplot(3,2,2)
#plt.title('5 %')
#plt.plot(pt_Initial_IIa_5["01:time"],Initial_alpha_5,'k-',lw=2.5, label= r'$FVIII_{5 \%}$' ,markevery=50,markersize=15,markerfacecolor='none')
#plt.plot(pt_IIa_5["01:time"],alpha_5,'k--',lw=2.5, label= r'$FVIII_{1 \%}$' ,markevery=350,markersize=5)
#plt.plot(time_exp,expe["FVIII 5%"],'ko',lw=2.5,label=r'$FVIII_{5\%}$ ', markevery=1)
#plt.xlim((0,1000))
#plt.ylim((-30,700))
#plt.setp(a4.get_xticklabels(), visible=False)
#
#a5=plt.subplot(3,2,4)
#plt.title('1 %')
#plt.plot(pt_Initial_IIa_1["01:time"],Initial_alpha_1,'k-',lw=2.5, label= r'$FVIII_{1 \%}$' ,markevery=350,markersize=5)
#plt.plot(pt_IIa_1["01:time"],alpha_1,'k--',lw=2.5, label= r'$FVIII_{1 \%}$' ,markevery=350,markersize=5)
#plt.plot(time_exp,expe["FVIII 1%"],'ko',lw=2.5,label=r'$FVIII_{1\%}$ ', markevery=1)
#plt.xlim((0,1000))
#plt.ylim((-30,700))
#plt.setp(a5.get_xticklabels(), visible=False)
#
#a6=plt.subplot(3,2,6)
#plt.title('0 %')
#plt.plot(pt_Initial_IIa_0["01:time"],Initial_alpha_0,'k-',lw=2.5, label= r'$FVIII_{0 \%}$' ,markevery=350,markersize=5)
#plt.plot(pt_IIa_0["01:time"],alpha_0,'k--',lw=2.5, label= r'$FVIII_{0 \%}$' ,markevery=350,markersize=5)
#plt.plot(time_exp,expe["FVIII 0%"],'ko',lw=2.5,label=r'$FVIII_{0\%}$ ', markevery=1)
#plt.xlabel(r'time (s) ')
#plt.xlim((0,1000))
#plt.ylim((-30,700))
plt.show()
