import numpy as np
import matplotlib.pyplot as plt
import csv

FVIIIETP={}
FVIIIPeak={}
FVIIITtpeak={}

readFile = open("FVIII_ETP.csv", 'r')
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
    FVIIIETP[columns[ic]]= t1
    ic = ic+1

readFile = open("FVIII_Peak.csv", 'r')
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
    FVIIIPeak[columns[ic]]= t1
    ic = ic+1

readFile = open("FVIII_Ttpeak.csv", 'r')
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
    FVIIITtpeak[columns[ic]]= t1
    ic = ic+1

FVIIIrange = [i for i, _ in enumerate(FVIIIETP['FVIII_(%)'])] # Position of plots

#fig, (ax1,ax2,ax3) = plt.subbars(3)

N = 6 # N= number of bars
ind = np.arange(N) # to use in bars position
width = 0.25  # With of bars
plt.style.use('ggplot')

font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 15}
#plt.rc('text', usetex=True)
plt.rc('font', **font)

plt.rcParams.update({'font.size': 15})

ax1=plt.subplot(3,1,1)
ax1.bar(ind+width+width,FVIIIETP['Experiments'],width,color='blue',label='Experiments')
ax1.bar(ind,FVIIIETP['Original_Hockin'],width,color='gray',label='Hockin original')
ax1.bar(ind+width,FVIIIETP['Modified_Hockin'],width,color='red',label='Hockin modified')
ax1.set_title("FVIII")
ax1.set_ylabel("ETP")


ax2=plt.subplot(3,1,2)
ax2.bar(ind,FVIIIPeak['Hockin_initial'],width,color='gray',label='Ext original FVIIa-TF = 1 pM')
ax2.bar(ind+width,FVIIIPeak['Hockin_modifi'],width,color='red',label='Ext modified')
ax2.bar(ind+width+width,FVIIIPeak['Experiments'],width,color='blue',label='Experiments')
ax2.set_ylabel(r"IIa$_{max}$")
plt.legend(loc='best')

ax3=plt.subplot(3,1,3)
ax3.bar(ind+width+width,FVIIITtpeak['Experiments'],width,color='blue',label='Experiments')
ax3.bar(ind,FVIIITtpeak['Hockin_initial'],width,color='gray',label='Hockin original')
ax3.bar(ind+width,FVIIITtpeak['Hockin_modifi'],width,color='red',label='Hockin modified')

plt.xticks(ind + width , ('100%', '50%', '15%', '5%', '1%', '0%'))
ax3.set_ylabel(r"$\tau_{max}$")

plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)

#import tikzplotlib
#tikzplotlib.save("TFVII.tex")
plt.show()
