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


#sys.exit()


plt.subplots_adjust(left=0.09, bottom=0.12, right=0.95, top=0.95,wspace=0.24, hspace=0.2)
plt.subplot(1,2,2)
plt.plot(log10(GridSphA.R),log10(GridSphA.Sigma2),'-x',label='Spherical',color='black',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphB.R),log10(GridSphB.Sigma2),'-o',label=r'$x$-axis',color='black',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphC.R),log10(GridSphC.Sigma2),'-D',label=r'$y$-axis',color='black',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphD.R),log10(GridSphD.Sigma2),'-<',label=r'$z$-axis',color='black',lw=2,ms=9,mew=2)
plt.text(0.5,-0.3,r'$\sigma^2$',color='black',fontsize=24)
plt.text(0.5,-0.42,r'$\sigma_\textrm{tan}^2$',color='blue',fontsize=24)
plt.text(0.5,-0.54,r'$\sigma_\textrm{rad}^2$',color='red',fontsize=24)
#plt.ylim((-1.05,-0.2))
#plt.xlim((-1,1))
plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=3)
#plt.grid()
#plt.xlabel(r'$\log r$',fontsize=24)
#plt.ylabel(r'$\log \sigma^2$',fontsize=24)
#SetLabels(20,20)

#plt.subplot(1,3,2)
plt.plot(log10(GridSphA.R),log10(GridSphA.Sigma2-GridSphA.Sigma2r),'-x',label='Spherical',color='blue',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphB.R),log10(GridSphB.Sigma2-GridSphB.Sigma2r),'-o',label=r'$x$-axis',color='blue',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphC.R),log10(GridSphC.Sigma2-GridSphC.Sigma2r),'-D',label=r'$y$-axis',color='blue',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphD.R),log10(GridSphD.Sigma2-GridSphD.Sigma2r),'-<',label=r'$z$-axis',color='blue',lw=2,ms=9,mew=2)
#plt.ylim((-1.3,-0.4))
#plt.xlim((-1,1))
#plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=3)
#plt.grid()
#plt.xlabel(r'$\log r$',fontsize=24)
#plt.ylabel(r'$\log \sigma^2_\textrm{tan}$',fontsize=24)
#SetLabels(20,20)

#plt.subplot(1,3,3)
plt.plot(log10(GridSphA.R),log10(GridSphA.Sigma2r),'-x',label='Spherical',color='red',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphB.R),log10(GridSphB.Sigma2r),'-o',label=r'$x$-axis',color='red',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphC.R),log10(GridSphC.Sigma2r),'-D',label=r'$y$-axis',color='red',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphD.R),log10(GridSphD.Sigma2r),'-<',label=r'$z$-axis',color='red',lw=2,ms=9,mew=2)
plt.ylim((-1.35,-0.2))
plt.xlim((-1,1))
#plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=3)
plt.grid()
plt.xlabel(r'$\log r$',fontsize=24)
plt.ylabel(r'$\log \sigma^2$',fontsize=24)
SetLabels(20,20)

plt.subplot(1,2,1)
IC = DM_structure.DM_structure('../1HqIso_000')
ICGridSph = IC.CreateGridLogBins(NBins=25,Rmin=0.001,Rmax=9.9)
plt.plot(log10(ICGridSph.R),ICGridSph.Beta,'--',label=r'IC',color='black',lw=2,ms=9,mew=2)

plt.plot(log10(GridSphA.R),GridSphA.Beta,'-x',label=r'Spherical',color='black',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphB.R),GridSphB.Beta,'-o',label=r'$x$-axis',color='blue',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphC.R),GridSphC.Beta,'-D',label=r'$y$-axis',color='red',lw=2,ms=9,mew=2)
if ImpactParameter:
    plt.plot(log10(GridSphD.R),GridSphD.Beta,'-<',label=r'$z$-axis',color='green',lw=2,ms=9,mew=2)
#plt.plot(log10(GridSphE.R),GridSphE.Beta,'->',label='neg. x-axis',color='grey',lw=2,ms=9)




plt.ylim((-0.45,0.75))
plt.xlim((-1,1))
plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=4)
plt.grid()
plt.xlabel(r'$\log r$',fontsize=24)
plt.ylabel(r'$\beta$',fontsize=24)
SetLabels(20,20)
plt.annotate(r'$r_{-2}$', xy=(-0.345, 0.5), xytext=(-0.345, 0.65),fontsize=26,ha='center',arrowprops=dict(facecolor='grey', shrink=0.05))
plt.show()


