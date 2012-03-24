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
        




FileName = '../Hq_0.1_0.5_080'
A = DM_structure.DM_structure(FileName)
#A.Snapshot.SelectIDsMin(1000000)
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


plt.subplots_adjust(left=0.07, bottom=0.15, right=0.97, top=0.91,wspace=0.16, hspace=0.2)


plt.subplot(1,3,1)
plt.title('All particles',fontsize=24)
plt.plot(log10(GridSphA.R),GridSphA.Beta,'-*',label=r'Spherically averaged',color='black',lw=2,ms=9)
plt.plot(log10(GridSphB.R),GridSphB.Beta,'-o',label=r'$x$-axis',color='blue',lw=2,ms=9)
plt.plot(log10(GridSphC.R),GridSphC.Beta,'-D',label=r'$y$-axis',color='red',lw=2,ms=9)

if ImpactParameter:
    plt.plot(log10(GridSphD.R),GridSphD.Beta,'-<',label=r'$z$-axis',color='green',lw=2,ms=9)
#plt.plot(log10(GridSphE.R),GridSphE.Beta,'->',label='neg. x-axis',color='grey',lw=2,ms=9)
plt.ylim((-0.1,0.6))
plt.xlim((-1,1))
#plt.legend(prop=dict(size=18), numpoints=3, ncol=1,frameon=True,loc=4)
plt.grid()
plt.xlabel(r'$\log r$',fontsize=24)
plt.ylabel(r'$\beta$',fontsize=24)
SetLabels(20,20)
plt.annotate(r'$r_{-2}$', xy=(-0.325, 0.2), xytext=(-0.325, 0.35),fontsize=26,ha='center',arrowprops=dict(facecolor='grey', shrink=0.05))

A = DM_structure.DM_structure(FileName)
A.Snapshot.SelectIDsMin(1000000)
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




plt.subplot(1,3,2)
plt.title('Small halo',fontsize=24)
plt.plot(log10(GridSphA.R),GridSphA.Beta,'-*',label=r'Spherically averaged',color='black',lw=2,ms=9)
plt.plot(log10(GridSphB.R),GridSphB.Beta,'-o',label=r'$x$-axis',color='blue',lw=2,ms=9)
plt.plot(log10(GridSphC.R),GridSphC.Beta,'-D',label=r'$y$-axis',color='red',lw=2,ms=9)
if ImpactParameter:
    plt.plot(log10(GridSphD.R),GridSphD.Beta,'-<',label=r'$z$-axis',color='green',lw=2,ms=9)
#plt.plot(log10(GridSphE.R),GridSphE.Beta,'->',label='neg. x-axis',color='grey',lw=2,ms=9)
plt.ylim((-0.0,1.0))
plt.xlim((-1,1))
#plt.legend(prop=dict(size=18), numpoints=3, ncol=1,frameon=True,loc=4)
plt.grid()
plt.xlabel(r'$\log r$',fontsize=24)
#plt.ylabel(r'$\beta$',fontsize=24)
SetLabels(20,20)



A = DM_structure.DM_structure(FileName)
A.Snapshot.SelectIDsMax(1000000)
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




plt.subplot(1,3,3)
plt.title('Main halo',fontsize=24)
plt.plot(log10(GridSphA.R),GridSphA.Beta,'-*',label=r'Spherically averaged',color='black',lw=2,ms=9)
plt.plot(log10(GridSphB.R),GridSphB.Beta,'-o',label=r'$x$-axis',color='blue',lw=2,ms=9)
plt.plot(log10(GridSphC.R),GridSphC.Beta,'-D',label=r'$y$-axis',color='red',lw=2,ms=9)
if ImpactParameter:
    plt.plot(log10(GridSphD.R),GridSphD.Beta,'-<',label=r'$z$-axis',color='green',lw=2,ms=9)
#plt.plot(log10(GridSphE.R),GridSphE.Beta,'->',label='neg. x-axis',color='grey',lw=2,ms=9)
plt.ylim((-0.1,0.4))
plt.xlim((-1,1))
plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=1)
plt.grid()
plt.xlabel(r'$\log r$',fontsize=24)
#plt.ylabel(r'$\beta$',fontsize=24)
SetLabels(20,20)



plt.show()


