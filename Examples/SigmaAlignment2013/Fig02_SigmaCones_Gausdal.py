#/usr/bin/python
import matplotlib.pyplot as plt
import Classes.Gadget2 as Gadget2
import Classes.DM_structure as DM_structure
from scipy import log10
import scipy
import pylab
from math import pi
import copy,sys
import numpy.random
import healpy as hp
import matplotlib


ImpactParameter = True 
plt.rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Roman']})
plt.rc('font',**{'family':'serif','serif':['Computer Modern Sans serif']})
plt.rc('text', usetex=True)

plt.rcParams["xtick.major.size"] = 5
plt.rcParams["ytick.major.size"] = 5

#Functions:
def SetLabels(xsize=18,ysize=18):
    ax = pylab.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(xsize)

    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(ysize)
        

#Filenames = ['0.1/1HqIso_Impact0_160','0.3/1HqIso_Impact0_160','0.5/1HqIso_Impact0_160','1.0/1HqIso_Impact0_160','1.5/1HqIso_Impact0_160']
Filenames = ['0.5/1HqIso_Impact0_160']
NFiles = len(Filenames)
DIR = '/home/ms/Uni/DarkMatter/AllSimulations/SigmaAlignment2013/HeadonMerger_VaryingVel/'


FileName='/home/ms/Uni/DarkMatter/AllSimulations/SigmaAlignment2013/HeadonMerger_VaryingVel/0.5/1HqIso_Impact0_160'
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
C.Snapshot.SelectParticlesInCone(1.0/scipy.sqrt(2.0),1.0/scipy.sqrt(2.0),0,3.1415/8.0)
D.Snapshot.SelectParticlesInCone(0,1,0,3.1415/8.0)
#E.Snapshot.SelectParticlesInCone(-1,0,0,3.1415/8.0)



GridSphA = A.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
GridSphB = B.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
GridSphC = C.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
GridSphD = D.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)


#IC = DM_structure.DM_structure('../1HqIso_000')
#ICGridSph = IC.CreateGridLogBins(NBins=25,Rmin=0.001,Rmax=9.9)
#plt.plot(log10(ICGridSph.R),ICGridSph.Beta,'--',label=r'IC',color='black',lw=2,ms=9,mew=2)

#plt.title(r'$x+y$-axis cone', fontsize=24)
plt.title('Velocity dispersion',fontsize=24)
plt.plot(log10(GridSphC.R),2*GridSphC.Ekinx,'-o',label=r'$\sigma^2$ (collision axis)',color='blue',lw=3,ms=14,mew=2)
plt.plot(log10(GridSphC.R),2*GridSphC.Ekiny,'-D',label=r'$\sigma^2$ (radial direction)',color='red',lw=3,ms=14,mew=2)
plt.plot(log10(GridSphC.R),2*GridSphC.Ekinz,'-<',label=r'$\sigma^2$ (tangential direction)',color='green',lw=3,ms=14,mew=2)


plt.ylim((-0.0,0.26))
plt.xlim((-1,1.19))
plt.legend(prop=dict(size=24), numpoints=2, ncol=1,frameon=True,loc=1)
plt.grid()
#plt.xlabel(r'$\log r$',fontsize=24)
plt.ylabel(r'$\sigma^2$',fontsize=24)
SetLabels(1,18)


plt.show()


