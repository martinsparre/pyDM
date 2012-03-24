#/usr/bin/python
import matplotlib.pyplot as plt
import Classes.Gadget2 as Gadget2
import Classes.DM_structure as DM_structure
from scipy import log10
import scipy
import pylab
from math import pi
import copy,sys


ImpactParameter = True 
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
        


FileName = '../1HqIso_Impact10_120'

A = DM_structure.DM_structure(FileName)
A.FindCenter()
A.FindCenterVel()
A.CenterParticlePositions()
A.CenterParticleVelocities()
GridSph = A.CreateGridLogBins(NBins=100,Rmin=0.001)

B = copy.deepcopy(A)
C = copy.deepcopy(A)
D = copy.deepcopy(A)

B.Snapshot.SelectParticlesInCone(1,0,0,3.1415/8.0)
C.Snapshot.SelectParticlesInCone(0,1,0,3.1415/8.0)
D.Snapshot.SelectParticlesInCone(0,0,1,3.1415/8.0)

GridSphA = A.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
GridSphB = B.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
GridSphC = C.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
GridSphD = D.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)



plt.subplots_adjust(left=0.10, bottom=0.13, right=0.95, top=0.92,wspace=0.0, hspace=0.2)
plt.subplot(1,2,1)
plt.title('With impact parameter',fontsize=28)
plt.plot(log10(GridSphA.R),log10(GridSphA.Rho/(GridSphA.Sigma2**1.5) / GridSphB.R**(-1.91) ),'-x',label='Spherical',color='black',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphB.R),log10(GridSphB.Rho/(GridSphB.Sigma2**1.5)/0.03806023375 / GridSphB.R**(-1.91)),'-o',label=r'$x$-axis',color='black',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphC.R),log10(GridSphC.Rho/(GridSphC.Sigma2**1.5)/0.03806023375 / GridSphC.R**(-1.91)),'-D',label=r'$y$-axis',color='black',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphD.R),log10(GridSphD.Rho/(GridSphD.Sigma2**1.5)/0.03806023375 / GridSphD.R**(-1.91)),'-<',label=r'$z$-axis',color='black',lw=2,ms=9,mew=2)
plt.text(0.28,0.5,r'$\rho/\sigma^3\times r^{1.91}$',color='black',fontsize=24)
plt.text(0.28,0.75,r'$2\rho/\sigma_\textrm{rad}^3\times r^{1.91}$ ',color='grey',fontsize=24)

plt.legend(prop=dict(size=18), numpoints=2, ncol=2,frameon=True,loc=2)

plt.plot(log10(GridSphA.R),log10(2*GridSphA.Rho/(GridSphA.Sigma2r**1.5) / GridSphA.R**(-1.91)),'-x',label='Spherical',color='grey',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphB.R),log10(2*GridSphB.Rho/(GridSphB.Sigma2r**1.5)/0.03806023375 / GridSphB.R**(-1.91) ),'-o',label=r'$x$-axis',color='grey',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphC.R),log10(2*GridSphC.Rho/(GridSphC.Sigma2r**1.5)/0.03806023375 / GridSphC.R**(-1.91)),'-D',label=r'$y$-axis',color='grey',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphD.R),log10(2*GridSphD.Rho/(GridSphD.Sigma2r**1.5)/0.03806023375 / GridSphD.R**(-1.91)),'-<',label=r'$z$-axis',color='grey',lw=2,ms=9,mew=2)

plt.ylim((-1.5,1.0))
plt.xlim((-1,0.99))
#plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=3)
plt.grid()
plt.xlabel(r'$\log r$',fontsize=24)
plt.ylabel(r'$\log \rho/\sigma^3 \times r^{1.91}$',fontsize=24)
SetLabels(20,20)


FileName = '../1HqIso_Impact0_121'

A = DM_structure.DM_structure(FileName)
A.FindCenter()
A.FindCenterVel()
A.CenterParticlePositions()
A.CenterParticleVelocities()
GridSph = A.CreateGridLogBins(NBins=100,Rmin=0.001)

B = copy.deepcopy(A)
C = copy.deepcopy(A)
D = copy.deepcopy(A)

B.Snapshot.SelectParticlesInCone(1,0,0,3.1415/8.0)
C.Snapshot.SelectParticlesInCone(0,1,0,3.1415/8.0)
D.Snapshot.SelectParticlesInCone(0,0,1,3.1415/8.0)

GridSphA = A.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
GridSphB = B.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
GridSphC = C.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
GridSphD = D.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)




plt.subplot(1,2,2)
plt.title('Without impact parameter',fontsize=28)
plt.plot(log10(GridSphA.R),log10(GridSphA.Rho/(GridSphA.Sigma2**1.5) / GridSphB.R**(-1.91) ),'-x',label='Spherical',color='black',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphB.R),log10(GridSphB.Rho/(GridSphB.Sigma2**1.5)/0.03806023375 / GridSphB.R**(-1.91)),'-o',label=r'$x$-axis',color='black',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphC.R),log10(GridSphC.Rho/(GridSphC.Sigma2**1.5)/0.03806023375 / GridSphC.R**(-1.91)),'-D',label=r'$y$-axis',color='black',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphD.R),log10(GridSphD.Rho/(GridSphD.Sigma2**1.5)/0.03806023375 / GridSphD.R**(-1.91)),'-<',label=r'$z$-axis',color='black',lw=2,ms=9,mew=2)
plt.text(0.28,0.5,r'$\rho/\sigma^3\times r^{1.91}$',color='black',fontsize=24)
plt.text(0.28,0.75,r'$2\rho/\sigma_\textrm{rad}^3\times r^{1.91}$',color='grey',fontsize=24)

plt.legend(prop=dict(size=18), numpoints=2, ncol=2,frameon=True,loc=2)

plt.plot(log10(GridSphA.R),log10(2*GridSphA.Rho/(GridSphA.Sigma2r**1.5) / GridSphA.R**(-1.91)),'-x',label='Spherical',color='grey',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphB.R),log10(2*GridSphB.Rho/(GridSphB.Sigma2r**1.5)/0.03806023375 / GridSphB.R**(-1.91) ),'-o',label=r'$x$-axis',color='grey',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphC.R),log10(2*GridSphC.Rho/(GridSphC.Sigma2r**1.5)/0.03806023375 / GridSphC.R**(-1.91)),'-D',label=r'$y$-axis',color='grey',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphD.R),log10(2*GridSphD.Rho/(GridSphD.Sigma2r**1.5)/0.03806023375 / GridSphD.R**(-1.91)),'-<',label=r'$z$-axis',color='grey',lw=2,ms=9,mew=2)

plt.ylim((-1.5,1.0))
plt.xlim((-0.99,1))
#plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=3)
plt.grid()
plt.xlabel(r'$\log r$',fontsize=24)
#plt.ylabel(r'$\log \rho/\sigma^3 \times r^{1.91}$',fontsize=24)
SetLabels(20,1)


plt.show()


