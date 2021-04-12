import numpy as np
import matplotlib.pyplot as plt
import csv

FXIETP={}
FXIPeak={}
FXITtpeak={}
FXIIETP={}
FXIIPeak={}
FXIITtpeak={}
FVIIIETP={}
FVIIIPeak={}
FVIIITtpeak={}

readFile = open("FXIETP.csv", 'r')
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
    FXIETP[columns[ic]]= t1
    ic = ic+1

readFile = open("FXIPeak.csv", 'r')
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
    FXIPeak[columns[ic]]= t1
    ic = ic+1

readFile = open("FXITtpeak.csv", 'r')
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
    FXITtpeak[columns[ic]]= t1
    ic = ic+1

readFile = open("FXIIETP.csv", 'r')
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
    FXIIETP[columns[ic]]= t1
    ic = ic+1

readFile = open("FXIIPeak.csv", 'r')
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
    FXIIPeak[columns[ic]]= t1
    ic = ic+1

readFile = open("FXIITtpeak.csv", 'r')
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
    FXIITtpeak[columns[ic]]= t1
    ic = ic+1

readFile = open("ContactFVIIIETP.csv", 'r')
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

readFile = open("ContactFVIIIPeak.csv", 'r')
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

readFile = open("ContactFVIIITtpeak.csv", 'r')
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

N = 6 # N= number of plots
ind = np.arange(N) # to use in plots position
width = 0.25  # With of plots
height = 0.6
plt.style.use('ggplot')

font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 15}
#plt.rc('text', usetex=True)
plt.rc('font', **font)

plt.rcParams.update({'font.size': 15})

ax1=plt.subplot(3,3,1)
plt.title('FXI')
ax1.bar(ind + width + width,FXIETP['Experiments'],width,color='blue',label='Experiments')
ax1.bar(ind,FXIETP['Chatterjee_initial'],width,color='gray',label='Chatterjee original')
ax1.bar(ind + width,FXIETP['Chatterjee_modifi'],width,color='red',label='Chatterjee modified')
ax1.set_ylabel("ETP")
plt.setp(ax1.get_xticklabels(),visible=False)

ax2=plt.subplot(3,3,4)
ax2.bar(ind+width+width,FXIPeak['Experiments'],width,color='blue',label='Experiments')
ax2.bar(ind,FXIPeak['Chatterjee_initial'],width,color='gray',label='Chatterjee original')
ax2.bar(ind+width,FXIPeak['Chatterjee_modifi'],width,color='red',label='Chatterjee modified')
ax2.set_ylabel(r"IIa$_{max}$")
plt.setp(ax2.get_xticklabels(),visible=False)

ax3=plt.subplot(3,3,7)
ax3.bar(ind+width+width,FXITtpeak['Experiments'],width,color='blue',label='Experiments')
ax3.bar(ind,FXITtpeak['Chatterjee_initial'],width,color='gray',label='Chatterjee original')
ax3.bar(ind+width,FXITtpeak['Chatterjee_modifi'],width,color='red',label='Chatterjee modified')
plt.xticks(ind + width , ('100%', '50%', '15%', '5%', '1%', '0%'))
ax3.set_ylabel(r"$\tau_{max}$")

ax4=plt.subplot(3,3,2)
plt.title('FXII')
ax4.bar(ind+width+width,FXIIETP['Experiments'],width,color='blue',label='Experiments')
ax4.bar(ind,FXIIETP['Chatterjee_initial'],width,color='gray',label='Chatterjee original')
ax4.bar(ind+width,FXIIETP['Chatterjee_modifi'],width,color='red',label='Chatterjee modified')
plt.setp(ax4.get_xticklabels(),visible=False)

ax5=plt.subplot(3,3,5)
ax5.bar(ind+width+width,FXIIPeak['Experiments'],width,color='blue',label='Experiments')
ax5.bar(ind,FXIIPeak['Chatterjee_initial'],width,color='gray',label='Int original')
ax5.bar(ind+width,FXIIPeak['Chatterjee_modifi'],width,color='red',label='Int modified')
plt.setp(ax5.get_xticklabels(),visible=False)

ax6=plt.subplot(3,3,8)
ax6.bar(ind,FXIITtpeak['Chatterjee_initial'],width,color='gray',label='Int original')
ax6.bar(ind+width,FXIITtpeak['Chatterjee_modifi'],width,color='red',label='Int modified')
ax6.bar(ind+width+width,FXIITtpeak['Experiments'],width,color='blue',label='Experiments')
plt.xticks(ind + width , ('100%', '50%', '15%', '5%', '1%', '0%'))
plt.legend(loc='best')


ax7=plt.subplot(3,3,3)
plt.title('FVIII')
ax7.bar(ind+width+width,FVIIIETP['Experiments'],width,color='blue',label='Experiments')
ax7.bar(ind,FVIIIETP['Chatterjee_initial'],width,color='gray',label='Chatterjee original')
ax7.bar(ind+width,FVIIIETP['Chatterjee_modifi'],width,color='red',label='Chatterjee modified')
plt.setp(ax7.get_xticklabels(),visible=False)

ax8=plt.subplot(3,3,6)
ax8.bar(ind+width+width,FVIIIPeak['Experiments'],width,color='blue',label='Experiments')
ax8.bar(ind,FVIIIPeak['Chatterjee_initial'],width,color='gray',label='Chatterjee original')
ax8.bar(ind+width,FVIIIPeak['Chatterjee_modifi'],width,color='red',label='Chatterjee modified')
plt.setp(ax8.get_xticklabels(),visible=False)

ax9=plt.subplot(3,3,9)
ax9.bar(ind+width+width,FVIIITtpeak['Experiments'],width,color='blue',label='Experiments')
ax9.bar(ind,FVIIITtpeak['Chatterjee_initial'],width,color='gray',label='Chatterjee original')
ax9.bar(ind+width,FVIIITtpeak['Chatterjee_modifi'],width,color='red',label='Chatterjee modified')
plt.xticks(ind + width , ('100%', '50%', '15%', '5%', '1%', '0%'))

#import tikzplotlib
#tikzplotlib.save("FXI.tex")
plt.show()

