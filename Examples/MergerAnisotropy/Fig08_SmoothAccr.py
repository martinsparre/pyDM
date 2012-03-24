#/usr/bin/python
import matplotlib.pyplot as plt
#import Classes.Gadget2 as ReadGadget2
import Classes.DM_structure as DM_structure
from scipy import log10,sqrt,arccos,arcsin
import scipy
import pylab
from math import pi



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
        

        
plt.subplots_adjust(left=0.12, bottom=0.12, right=0.95, top=0.95,wspace=0.24, hspace=0.2)
plt.subplot(1,1,1)
        
        
FileName='../3RemnantsAnd2000000_000'
A = DM_structure.DM_structure(FileName)
A.FindCenter()
A.FindCenterVel()
A.CenterParticlePositions()
A.CenterParticleVelocities()
A.Snapshot.SelectIDsMin(2000000)
GridSph = A.CreateGridLogBins(NBins=20,Rmax=12.59, Rmin=4.01)
plt.plot(log10(GridSph.R[3:-3]),log10(GridSph.Rho[3:-3])  ,'-o',color='black',label='IC, Infall particles',lw=2,ms=9)


FileName='../A3RemnantsAnd2000000_040'
A = DM_structure.DM_structure(FileName)
A.FindCenter()
A.FindCenterVel()
A.CenterParticlePositions()
A.CenterParticleVelocities()
A.Snapshot.SelectIDsMin(2000000)
GridSph = A.CreateGridLogBins(NBins=40)
plt.plot(log10(GridSph.R),log10(GridSph.Rho)  ,'-<',color='black',label='Final, Infall particles',lw=2,ms=9)


FileName='../3RemnantsAnd2000000_000'
A = DM_structure.DM_structure(FileName)
A.FindCenter()
A.FindCenterVel()
A.CenterParticlePositions()
A.CenterParticleVelocities()
A.Snapshot.SelectIDsMax(2000000)
GridSph = A.CreateGridLogBins(NBins=40,Rmin=0.01)
plt.plot(log10(GridSph.R),log10(GridSph.Rho)  ,'-o',color='grey',label='IC, Merger remnant',lw=2,ms=9)


FileName='../A3RemnantsAnd2000000_040'
A = DM_structure.DM_structure(FileName)
A.FindCenter()
A.FindCenterVel()
A.CenterParticlePositions()
A.CenterParticleVelocities()
A.Snapshot.SelectIDsMax(2000000)
GridSph = A.CreateGridLogBins(NBins=40,Rmin=0.017)
plt.plot(log10(GridSph.R),log10(GridSph.Rho)  ,'-<',color='grey',label='Final, Merger remnant',lw=2,ms=9)


plt.legend()

plt.ylim((-6,1))
plt.xlim((-1,1.5))
plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=1)
plt.grid()
plt.xlabel(r'$\log r$',fontsize=24)
plt.ylabel(r'$\log\rho$',fontsize=24)



SetLabels()
plt.show()

