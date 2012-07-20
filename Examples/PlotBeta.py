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
        
#FileName = '/home/ms/Uni/DarkMatter/SimulationDataAndAnalysis/ROI/OM_ROI00_rAN0.2_HQ_000'
#FileName = '/home/ms/Uni/DarkMatter/AllSimulations/MergerAnisotropy/1HqIso_Impact0_121'
#FileName = '/home/ms/Uni/DarkMatter/AllSimulations/MergerAnisotropy/1HqIso_Impact10_120'
#FileName = '/home/ms/Uni/DarkMatter/AllSimulations/CloneMerge/Mergers7_HQA_081'
#FileName = '/home/ms/Uni/DarkMatter/AllSimulations/MergerAnisotropy/1HQOM_Impact0_121'
#FileName = '/home/ms/Uni/DarkMatter/AllSimulations/MergerAnisotropy/Hq_0.1_0.5_080'
FileName = '/home/ms/Uni/DarkMatter/AllSimulations/MergerAnisotropy/Hq_0.1_0.5Impact_089'
#FileName = '/home/ms/Uni/DarkMatter/SimulationDataAndAnalysis/MergerAnisotropy2012/Simulations/PyDM1/1HqIso_000'

A = DM_structure.DM_structure(FileName)
A.FindCenter()
A.FindCenterVel()
A.CenterParticlePositions()
A.CenterParticleVelocities()
GridSph = A.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50,CalcUncBeta=True)

plt.figure(1)
plt.subplot(1,2,1)


Chi2List = []
Ncones = 0


BetaCones = []

AngleFile = open('../192.txt','r')
#AngleFile = open('../192.txt','r')
for line in AngleFile:
    tmp = line.split()
    if len(tmp) != 5:
        continue
    
    x = float(tmp[0])
    y = float(tmp[1])
    z = float(tmp[2])

    B = copy.deepcopy(A)
    B.Snapshot.SelectParticlesInCone(x,y,z,3.1415/8.0)

    GridSphB = B.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50,CalcUncBeta=True)

    plt.subplot(1,2,1)
    plt.plot(log10(GridSphB.R),GridSphB.Beta,'-',color='black',lw=1)
    BetaCones.append(scipy.mean(GridSphB.Beta[18:23]))
    print scipy.mean(GridSphB.R[18:23])
    Chi2List.append((GridSphB.MeanBetaBootstrap-GridSph.MeanBetaBootstrap)**2 / (GridSphB.UncBetaBootstrap**2+0.025**2))
    Ncones += 1
    print 'Cone number ',Ncones
Chi2List = scipy.array(Chi2List)
AngleFile.close()



plt.plot(log10(GridSph.R),GridSph.Beta,'-o',label=r'Spherical',color='blue',lw=2,ms=9,mew=2)


plt.ylim((-0.45,0.75))
plt.xlim((-1,1))
plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=4)
plt.grid()
plt.xlabel(r'$\log r$',fontsize=24)
plt.ylabel(r'$\beta$',fontsize=24)
SetLabels(20,20)
plt.annotate(r'$r_{-2}$', xy=(-0.345, 0.5), xytext=(-0.345, 0.65),fontsize=26,ha='center',arrowprops=dict(facecolor='grey', shrink=0.05))


plt.subplot(1,2,2)



#Particles[0][numpy.random.randint(0,high=len(Particles[0]),size=len(Particles[0]))]

#calc Chi2:
Chi2Total = Chi2List[0]*0.0
for chi2 in Chi2List:
    Chi2Total += chi2
Chi2Total /= len(Chi2List)
    
#bootstrap to get uncertainty:
BootstrapLists = []
for i in range(100):
    IDs = numpy.random.randint(0,high=len(Chi2List),size=len(Chi2List))
    NewChi2List = Chi2List[[IDs]]
    NewChi2Total = Chi2List[0]*0.0
    for chi2 in NewChi2List:
        NewChi2Total += chi2
    NewChi2Total /= len(NewChi2List)
    BootstrapLists.append(NewChi2Total)

Mean = scipy.mean(BootstrapLists,axis=0)
Std = scipy.std(BootstrapLists,axis=0)

plt.errorbar(log10(GridSphB.R), Mean, yerr=Std, color='black')
plt.plot(log10(GridSphB.R),Chi2Total)
plt.xlim((-1,1))
plt.ylim((0.0,20.0))
plt.grid()
plt.xlabel(r'$\log r$',fontsize=24)
plt.ylabel(r'$\chi^2$',fontsize=24)




import numpy as np
import healpy as hp
plt.figure(2)
NSIDE = 4
m = scipy.array(BetaCones)

hp.mollview(m, title="",min=-0.2,max=0.35)
plt.show()
