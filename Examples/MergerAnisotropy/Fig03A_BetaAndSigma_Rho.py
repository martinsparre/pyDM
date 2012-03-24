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
#E = copy.deepcopy(A)

B.Snapshot.SelectParticlesInCone(1,0,0,3.1415/8.0)
C.Snapshot.SelectParticlesInCone(0,1,0,3.1415/8.0)
D.Snapshot.SelectParticlesInCone(0,0,1,3.1415/8.0)
#E.Snapshot.SelectParticlesInCone(-1,0,0,3.1415/8.0)


PlotParticles = False
if PlotParticles == True:
    plt.title('File'+sys.argv[1])
    plt.plot(B.Snapshot.x,B.Snapshot.y,'.',label='$(1,0,0)$ $\pi/4$')
    plt.plot(C.Snapshot.x,C.Snapshot.y,'.',label='$(0,1,0)$ $\pi/16$')
#    plt.plot(D.Snapshot.x,D.Snapshot.y,'.',label='$(0,0,1)$ $\pi/16$')
#    plt.plot(E.Snapshot.x,E.Snapshot.y,'.',label='$(-1,0,0)$ $\pi/16$')
    plt.legend(loc=2)
    plt.grid()
    plt.xlabel('x')
    plt.xlabel('y')
    SetLabels()
    plt.show()


GridSphA = A.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
GridSphB = B.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
GridSphC = C.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
GridSphD = D.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
#GridSphE = E.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)




plt.subplots_adjust(left=0.1, bottom=0.1, right=0.97, top=0.97,wspace=0.28, hspace=0.25)









plt.subplot(2,2,1)

IDs = scipy.where((GridSphA.R > 0.05 )*(GridSphA.R < 10**(0.7) ))
plt.plot(log10(GridSphA.Rho[IDs]),GridSphA.Beta[IDs],'-x',label=r'Spherical',color='black',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphB.Rho[IDs]/0.03806023375),GridSphB.Beta[IDs],'-o',label=r'$x$-axis',color='blue',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphC.Rho[IDs]/0.03806023375),GridSphC.Beta[IDs],'-D',label=r'$y$-axis',color='red',lw=2,ms=9,mew=2)
if ImpactParameter:
    plt.plot(log10(GridSphD.Rho[IDs]/0.03806023375),GridSphD.Beta[IDs],'-<',label=r'$z$-axis',color='green',lw=2,ms=9,mew=2)
#plt.plot(log10(GridSphE.R),GridSphE.Beta,'->',label='neg. x-axis',color='grey',lw=2,ms=9)
#plt.ylim((-0.3,0.75))
plt.ylim((-0.3,0.5))
plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=3)
plt.grid()
plt.xlabel(r'$\log \rho$',fontsize=24)
plt.ylabel(r'$\beta$',fontsize=24)
SetLabels(20,20)








plt.subplot(2,2,2)

plt.plot(log10(GridSphA.Rho[IDs]),log10(GridSphA.Sigma2[IDs]),'-x',label='Spherical',color='black',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphB.Rho[IDs]/0.03806023375),log10(GridSphB.Sigma2[IDs]),'-o',label=r'$x$-axis',color='black',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphC.Rho[IDs]/0.03806023375),log10(GridSphC.Sigma2[IDs]),'-D',label=r'$y$-axis',color='black',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphD.Rho[IDs]/0.03806023375),log10(GridSphD.Sigma2[IDs]),'-<',label=r'$z$-axis',color='black',lw=2,ms=9,mew=2)
plt.text(-3.8,-0.30,r'$\sigma^2$',color='black',fontsize=24)
plt.text(-3.8,-0.46,r'$\sigma_\textrm{tan}^2$',color='blue',fontsize=24)
plt.text(-3.8,-0.62,r'$\sigma_\textrm{rad}^2$',color='red',fontsize=24)

#plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=4)

plt.plot(log10(GridSphA.Rho[IDs]),log10(GridSphA.Sigma2[IDs]-GridSphA.Sigma2r[IDs]),'-x',label='Spherical',color='blue',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphB.Rho[IDs]/0.03806023375),log10(GridSphB.Sigma2[IDs]-GridSphB.Sigma2r[IDs]),'-o',label=r'$x$-axis',color='blue',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphC.Rho[IDs]/0.03806023375),log10(GridSphC.Sigma2[IDs]-GridSphC.Sigma2r[IDs]),'-D',label=r'$y$-axis',color='blue',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphD.Rho[IDs]/0.03806023375),log10(GridSphD.Sigma2[IDs]-GridSphD.Sigma2r[IDs]),'-<',label=r'$z$-axis',color='blue',lw=2,ms=9,mew=2)

plt.plot(log10(GridSphA.Rho[IDs]),log10(GridSphA.Sigma2r[IDs]),'-x',label='Spherical',color='red',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphB.Rho[IDs]/0.03806023375),log10(GridSphB.Sigma2r[IDs]),'-o',label=r'$x$-axis',color='red',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphC.Rho[IDs]/0.03806023375),log10(GridSphC.Sigma2r[IDs]),'-D',label=r'$y$-axis',color='red',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphD.Rho[IDs]/0.03806023375),log10(GridSphD.Sigma2r[IDs]),'-<',label=r'$z$-axis',color='red',lw=2,ms=9,mew=2)
plt.ylim((-1.3,-0.2))
#plt.xlim((-1,1))
#plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=3)
plt.grid()
plt.xlabel(r'$\log \rho$',fontsize=24)
plt.ylabel(r'$\log \sigma^2$',fontsize=24)
SetLabels(20,20)



plt.subplot(2,2,3)

plt.plot(log10(GridSphA.R[IDs]),log10(GridSphA.Rho[IDs]),'-x',label='Spherical',color='black',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphB.R[IDs]),log10(GridSphB.Rho[IDs]/0.03806023375),'-o',label=r'$x$-axis',color='blue',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphC.R[IDs]),log10(GridSphC.Rho[IDs]/0.03806023375),'-D',label=r'$y$-axis',color='red',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphD.R[IDs]),log10(GridSphD.Rho[IDs]/0.03806023375),'-<',label=r'$z$-axis',color='green',lw=2,ms=9,mew=2)

plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=3)


#plt.xlim((-1,1))
#plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=3)
plt.grid()
plt.xlabel(r'$\log r$',fontsize=24)
plt.ylabel(r'$\log \rho$',fontsize=24)
SetLabels(20,20)

plt.ylim(-4,0.6)

plt.subplot(2,2,4)

plt.plot(log10(GridSphA.R[IDs]),GridSphA.Gamma[IDs],'-x',label='Spherical',color='black',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphB.R[IDs]),GridSphB.Gamma[IDs],'-o',label=r'$x$-axis',color='blue',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphC.R[IDs]),GridSphC.Gamma[IDs],'-D',label=r'$y$-axis',color='red',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphD.R[IDs]),GridSphD.Gamma[IDs],'-<',label=r'$z$-axis',color='green',lw=2,ms=9,mew=2)

plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=3)


#plt.xlim((-1,1))
#plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=3)
plt.grid()
plt.xlabel(r'$\log r$',fontsize=24)
plt.ylabel(r'$\gamma$',fontsize=24)
SetLabels(20,20)

plt.ylim(-4.8,-0.0)


plt.show()


