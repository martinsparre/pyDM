#!/usr/bin/python
import scipy
import pylab
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import Classes.DM_structure as DM_structure
from scipy import log10
#Fonts:

plt.rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Roman']})
plt.rc('font',**{'family':'serif','serif':['Computer Modern Sans serif']})
plt.rc('text', usetex=True)

plt.rcParams["xtick.major.size"] = 5
plt.rcParams["ytick.major.size"] = 5

#Functions:
def SetLabels(xsize=24,ysize=24):
    ax = pylab.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(xsize)

    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(ysize)


def ReadFile(File):
    f = open(File,'r+')

    betaR = []
    gammaR = []
    kappaR = []

    for line in f:
      if '#' in line:
        continue
      if len(line)<5:
          continue
      
      tmp = line.split()
      betaR.append(float(tmp[4]))
      gammaR.append(float(tmp[5]))
      kappaR.append(float(tmp[6]))
   
   
    return scipy.array(betaR),scipy.array(gammaR),scipy.array(kappaR)


#plot attractor:
PlotAttr = lambda AX,Beta,Gamma,Kappa: AX.fill_between(Beta, (Gamma+Kappa)+0.01, (Gamma+Kappa)-0.01,facecolor = 'gray',edgecolor = 'gray')

beta = []
gamma = []
kappa = []




#for f in ['s1_16.txt','s2_16.txt','s3_16.txt','s4_16.txt']:
#for f in ['00-515_001.txt','0E21_001.txt','hoover29_001.txt','om0-3.5E18_001.txt','omE21_001.txt','s1_16.txt','s2_16.txt','s3_16.txt','s4_16.txt']:
for f in ['00-515_001.txt','0E21_001.txt','hoover29_001.txt','om0-3.5E18_001.txt','omE21_001.txt','s1_16.txt','s2_16.txt','s3_16.txt','s4_16.txt']:
    temp = ReadFile('../HJS2010/'+f)
    beta.append(temp[0])
    gamma.append(temp[1])
    kappa.append(temp[2])


#colors, styles...:
R = scipy.array([1.0,0.0,0.0])
G = scipy.array([0.0,1.0,0.0])
B = scipy.array([0.0,0.0,1.0])
Color = ['black','black','black','black',0.5*(R+G+B),0.5*(R+G+B),0.5*(R+G+B)]
DotStyle = ['-','--','-.',':','-','--','-.']

#Create plot....
plt.subplots_adjust(left=0.06, bottom=0.12, right=0.98, top=0.94,wspace=0.0, hspace=0.5)
ax = plt.subplot(1,3,1)


for i in range(len(beta)):
    ax.fill_between(gamma[i], beta[i]+0.01, beta[i]-0.01,facecolor = 'gray',edgecolor = 'gray')

#plt.title(r'A different projection',fontsize=32)
SetLabels()
plt.xlabel(r'$\gamma$',fontsize=32)
plt.ylabel(r'$\beta$',fontsize=32)






#b,g,k = ReadFile('Infall_200.txt')
#plt.plot(g,b, linestyle='-',marker='o',color='black',label='Infall-simulation',markersize=9)


t=scipy.linspace(-3,-0.8,10)

#print t
plt.plot(t, -0.2*(t+0.8), '--',label=r'$\beta = -0.2 \times (\gamma + 0.8)$',color='black', lw=4)



A = DM_structure.DM_structure('../1HqIso_Impact0_121')
A.FindCenter()
A.FindCenterVel()
A.CenterParticlePositions()
A.CenterParticleVelocities()
AGridSph = A.CreateGridLogBins(NBins=25,Rmin=0.025,Rmax=8.0)
AGridSph.RemoveBadRegions()
plt.plot(AGridSph.Gamma,AGridSph.Beta,'-x',label=r'$\beta=0$-remnant',color='blue',lw=2,ms=9,mew=2)

B = DM_structure.DM_structure('../1HQOM_Impact0_121')
B.FindCenter()
B.FindCenterVel()
B.CenterParticlePositions()
B.CenterParticleVelocities()
BGridSph = B.CreateGridLogBins(NBins=25,Rmin=0.025,Rmax=8.0)
BGridSph.RemoveBadRegions()
plt.plot(BGridSph.Gamma,BGridSph.Beta,'-o',label=r'OM-remnant',color='red',lw=2,ms=9,mew=2)

C = DM_structure.DM_structure('../1HqIso_Impact10_120')
C.FindCenter()
C.FindCenterVel()
C.CenterParticlePositions()
C.CenterParticleVelocities()
CGridSph = C.CreateGridLogBins(NBins=25,Rmin=0.025,Rmax=8.0)
CGridSph.RemoveBadRegions()
plt.plot(CGridSph.Gamma,CGridSph.Beta,'-D',label=r'$\beta=0$, Impact',color='green',lw=2,ms=9,mew=2)


p = Rectangle((0, 0), 1, 1, facecolor="gray",edgecolor="gray")
a,b = ax.get_legend_handles_labels()
a.insert(0,p)
b.insert(0,"Attractor")
ax.legend(a,b,prop=dict(size=20), numpoints=2, ncol=1,frameon=True,loc=1)

plt.xlim((-3.99,-0.01))
plt.ylim((-0.1,0.9))
plt.grid()












ax = plt.subplot(1,3,2)
plt.title(r'$x$-axis',fontsize=32)

for i in range(len(beta)):
    ax.fill_between(gamma[i], beta[i]+0.01, beta[i]-0.01,facecolor = 'gray',edgecolor = 'gray')

#plt.title(r'A different projection',fontsize=32)
SetLabels(ysize=1)
plt.xlabel(r'$\gamma$',fontsize=32)
#plt.ylabel(r'$\beta$',fontsize=32)






#b,g,k = ReadFile('Infall_200.txt')
#plt.plot(g,b, linestyle='-',marker='o',color='black',label='Infall-simulation',markersize=9)


t=scipy.linspace(-3,-0.8,10)

#print t
plt.plot(t, -0.2*(t+0.8), '--',label=r'$\beta = -0.2 \times (\gamma + 0.8)$',color='black', lw=4)



A = DM_structure.DM_structure('../1HqIso_Impact0_121')
A.FindCenter()
A.FindCenterVel()
A.CenterParticlePositions()
A.CenterParticleVelocities()
A.Snapshot.SelectParticlesInCone(1,0,0,3.1415/8.0)
AGridSph = A.CreateGridLogBins(NBins=25,Rmin=0.025,Rmax=8.0)
AGridSph.RemoveBadRegions()
plt.plot(AGridSph.Gamma,AGridSph.Beta,'-x',label=r'$\beta=0$-remnant',color='blue',lw=2,ms=9,mew=2)

B = DM_structure.DM_structure('../1HQOM_Impact0_121')
B.FindCenter()
B.FindCenterVel()
B.CenterParticlePositions()
B.CenterParticleVelocities()
B.Snapshot.SelectParticlesInCone(1,0,0,3.1415/8.0)
BGridSph = B.CreateGridLogBins(NBins=25,Rmin=0.025,Rmax=8.0)
BGridSph.RemoveBadRegions()
plt.plot(BGridSph.Gamma,BGridSph.Beta,'-o',label=r'OM-remnant',color='red',lw=2,ms=9,mew=2)

C = DM_structure.DM_structure('../1HqIso_Impact10_120')
C.FindCenter()
C.FindCenterVel()
C.CenterParticlePositions()
C.CenterParticleVelocities()
C.Snapshot.SelectParticlesInCone(1,0,0,3.1415/8.0)
CGridSph = C.CreateGridLogBins(NBins=25,Rmin=0.025,Rmax=8.0)
CGridSph.RemoveBadRegions()
plt.plot(CGridSph.Gamma,CGridSph.Beta,'-D',label=r'$\beta=0$, Impact',color='green',lw=2,ms=9,mew=2)


#p = Rectangle((0, 0), 1, 1, facecolor="gray",edgecolor="gray")
#a,b = ax.get_legend_handles_labels()
#a.insert(0,p)
#b.insert(0,"Attractor")
#ax.legend(a,b,prop=dict(size=20), numpoints=1, ncol=1,frameon=True,loc='best')

plt.xlim((-3.99,-0.01))
plt.ylim((-0.1,0.9))
plt.grid()










ax = plt.subplot(1,3,3)
plt.title(r'$y$-axis',fontsize=32)

for i in range(len(beta)):
    ax.fill_between(gamma[i], beta[i]+0.01, beta[i]-0.01,facecolor = 'gray',edgecolor = 'gray')

#plt.title(r'A different projection',fontsize=32)
SetLabels(ysize=1)
plt.xlabel(r'$\gamma$',fontsize=32)
#plt.ylabel(r'$\beta$',fontsize=32)






#b,g,k = ReadFile('Infall_200.txt')
#plt.plot(g,b, linestyle='-',marker='o',color='black',label='Infall-simulation',markersize=9)


t=scipy.linspace(-3,-0.8,10)

#print t
plt.plot(t, -0.2*(t+0.8), '--',label=r'$\beta = -0.2 \times (\gamma + 0.8)$',color='black', lw=4)



A = DM_structure.DM_structure('../1HqIso_Impact0_121')
A.FindCenter()
A.FindCenterVel()
A.CenterParticlePositions()
A.CenterParticleVelocities()
A.Snapshot.SelectParticlesInCone(0,1,0,3.1415/8.0)
AGridSph = A.CreateGridLogBins(NBins=25,Rmin=0.025,Rmax=8.0)
AGridSph.RemoveBadRegions()
plt.plot(AGridSph.Gamma,AGridSph.Beta,'-x',label=r'$\beta=0$-remnant',color='blue',lw=2,ms=9,mew=2)

B = DM_structure.DM_structure('../1HQOM_Impact0_121')
B.FindCenter()
B.FindCenterVel()
B.CenterParticlePositions()
B.CenterParticleVelocities()
B.Snapshot.SelectParticlesInCone(0,1,0,3.1415/8.0)
BGridSph = B.CreateGridLogBins(NBins=25,Rmin=0.025,Rmax=8.0)
BGridSph.RemoveBadRegions()
plt.plot(BGridSph.Gamma,BGridSph.Beta,'-o',label=r'OM-remnant',color='red',lw=2,ms=9,mew=2)

C = DM_structure.DM_structure('../1HqIso_Impact10_120')
C.FindCenter()
C.FindCenterVel()
C.CenterParticlePositions()
C.CenterParticleVelocities()
C.Snapshot.SelectParticlesInCone(0,1,0,3.1415/8.0)
CGridSph = C.CreateGridLogBins(NBins=25,Rmin=0.025,Rmax=8.0)
CGridSph.RemoveBadRegions()
plt.plot(CGridSph.Gamma,CGridSph.Beta,'-D',label=r'$\beta=0$, Impact',color='green',lw=2,ms=9,mew=2)


#p = Rectangle((0, 0), 1, 1, facecolor="gray",edgecolor="gray")
#a,b = ax.get_legend_handles_labels()
#a.insert(0,p)
#b.insert(0,"Attractor")
#ax.legend(a,b,prop=dict(size=20), numpoints=1, ncol=1,frameon=True,loc='best')

plt.xlim((-3.99,-0.01))
plt.ylim((-0.1,0.9))
plt.grid()















plt.show()

