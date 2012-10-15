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
def SetLabels(xsize=20,ysize=20):
    ax = pylab.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(xsize)

    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(ysize)
        
BetaNOrder = 5

FileNames = ['00-5G20_001','0G20_001','0GManySmall_39_001','hooverG20_001','om0-3.5G20_001','OMG20_001','s1G20_001','s2G20_001','s3G20_001','s4G20_001']
PATH = '/home/ms/Uni/DarkMatter/AllSimulations/G/'


for f in FileNames:
    A = DM_structure.DM_structure(PATH+f)
    A.FindCenter()
    A.FindCenterVel()
    A.CenterParticlePositions()
    A.CenterParticleVelocities()
    A. SetBetaNOrder(BetaNOrder)
    GridSph = A.CreateGridLogBins(NBins=100,Rmin=0.001)

    GridSphA = A.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)


    plt.plot(GridSphA.Beta[13:-1],GridSphA.BetaN[13:-1]-GridSphA.Beta[13:-1],'-o',color='blue',lw=1,ms=6,mew=1)


plt.text(-0.9,0.8,r'$G$',color='blue',fontsize=18)

        
    
FileNames = ['00-515_001', '0E21_001','hoover29_001','om0-3.5E18_001', 'omE21_001','s1_16','s2_16','s3_16','s4_16']
PATH='/home/ms/Uni/DarkMatter/AllSimulations/HJS2010/'
    
for f in FileNames:
    A = DM_structure.DM_structure(PATH+f)
    A.FindCenter()
    A.FindCenterVel()
    A.CenterParticlePositions()
    A.CenterParticleVelocities()
    A. SetBetaNOrder(BetaNOrder)
    GridSph = A.CreateGridLogBins(NBins=100,Rmin=0.001)

    GridSphA = A.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
    
    plt.plot(GridSphA.Beta,GridSphA.BetaN-GridSphA.Beta,'-x',color='red',lw=1,ms=6,mew=1)

plt.text(-0.9,0.7,'HJS',color='red',fontsize=18)    
    
    
FileNames = ['05ColdCollapse_120','Infall_sub_1e6.bin_201']
PATH='/home/ms/Uni/DarkMatter/AllSimulations/InhomogeneousInfall/'
    
    
for f in FileNames:
    A = DM_structure.DM_structure(PATH+f)
    A.FindCenter()
    A.FindCenterVel()
    A.CenterParticlePositions()
    A.CenterParticleVelocities()
    A. SetBetaNOrder(BetaNOrder)
    GridSph = A.CreateGridLogBins(NBins=100,Rmin=0.001)

    GridSphA = A.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)

    plt.plot(GridSphA.Beta,GridSphA.BetaN-GridSphA.Beta,'-x',color='black',lw=1,ms=6,mew=1)

plt.text(-0.9,0.6,'Cold collapse',color='black',fontsize=18)    
    
    

A = DM_structure.DM_structure('/home/ms/Uni/DarkMatter/AllSimulations/MergerAnisotropy/1HqIso_Impact0_140')
A.FindCenter()
A.FindCenterVel()
A.CenterParticlePositions()
A.CenterParticleVelocities()
A. SetBetaNOrder(BetaNOrder)
GridSph = A.CreateGridLogBins(NBins=100,Rmin=0.001)

GridSphA = A.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)

B = copy.deepcopy(A)
C = copy.deepcopy(A)
D = copy.deepcopy(A)
#E = copy.deepcopy(A)

B.Snapshot.SelectParticlesInCone(1,0,0,3.1415/8.0)
C.Snapshot.SelectParticlesInCone(0,1,0,3.1415/8.0)
D.Snapshot.SelectParticlesInCone(0,0,1,3.1415/8.0)

GridSphB = B.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
GridSphC = C.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
GridSphD = D.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)

plt.plot(GridSphA.Beta[13:-10],GridSphA.BetaN[13:-10]-GridSphA.Beta[13:-10],'-<',color='green',lw=1,ms=8,mew=1)
plt.plot(GridSphB.Beta[13:-10],GridSphB.BetaN[13:-10]-GridSphB.Beta[13:-10],'-<',color='green',lw=1,ms=8,mew=1)
plt.plot(GridSphC.Beta[13:-10],GridSphC.BetaN[13:-10]-GridSphC.Beta[13:-10],'-<',color='green',lw=1,ms=8,mew=1)
plt.plot(GridSphD.Beta[13:-10],GridSphD.BetaN[13:-10]-GridSphD.Beta[13:-10],'-<',color='green',lw=1,ms=8,mew=1)

    
plt.text(-0.9,0.5,r'Merger',color='green',fontsize=18)



FileNames = ['OMG00_001', 'om0-3.5G00_001']
PATH = '/home/ms/Uni/DarkMatter/AllSimulations/G/'



for f in FileNames:
    A = DM_structure.DM_structure(PATH+f)
    A.FindCenter()
    A.FindCenterVel()
    A.CenterParticlePositions()
    A.CenterParticleVelocities()
    A. SetBetaNOrder(BetaNOrder)

    GridSphA = A.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)

    plt.plot(GridSphA.Beta[13:-1],GridSphA.BetaN[13:-1]-GridSphA.Beta[13:-1],'->',color='maroon',lw=1,ms=6,mew=1)


plt.text(-0.9,0.4,r'Osipkov Merritt',color='maroon',fontsize=18) 


x=scipy.arange(-10,10,0.1)
plt.plot(x,0*x,'-',color='grey',lw=2)




plt.xlabel(r'$\beta$',fontsize=24)
plt.ylabel(r'$\beta_5 - \beta$',fontsize=24)
plt.xlim((-1,1))
plt.ylim((-1,1))
plt.grid()
SetLabels()

plt.show()

