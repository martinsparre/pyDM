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
        


ShowSim = 4*[False]+[True]+[False]        
#ShowSim = [False,False,False,False,True,False]             
        
        
SIMULATION_DIR = '/home/ms/Uni/DarkMatter/AllSimulations/G/'
col = ['red', 'blue', 'green', 'maroon','orange','black','pink','yellow','cyan']
mar = ['o','x','+','D','<','>','*','v','p']


ax = plt.subplot(3,2,1)
ax.text(-0.8, 0.9, r'I) $G$-variations', fontsize=20)
plt.subplots_adjust(left=0.11, bottom=0.07, right=0.98, top=0.98,wspace=0.0, hspace=0.0)
Files = [SIMULATION_DIR+'00-5G20_001',SIMULATION_DIR+'s2G20_001',SIMULATION_DIR+'0G20_001',SIMULATION_DIR+'OMG20_001',SIMULATION_DIR+'s3G20_001',SIMULATION_DIR+'s1G20_001',SIMULATION_DIR+'s4G20_001',SIMULATION_DIR+'om0-3.5G20_001']


#i = 1, gives a somewhat unstable beta
for i in range(len(Files)):
    if i!=7:
        continue
    if ShowSim[0] == False:
        break

    FileName = Files[i]

    A = DM_structure.DM_structure(FileName)
    A.FindCenter()
    A.FindCenterVel()
    A.CenterParticlePositions()
    A.CenterParticleVelocities()
    GridSph = A.CreateGridLogBins(NBins=100,Rmin=0.001)


#    B = copy.deepcopy(A)
#    C = copy.deepcopy(A)
#    D = copy.deepcopy(A)


#    B.Snapshot.SelectParticlesInCone(1,0,0,3.1415/8.0)
#    C.Snapshot.SelectParticlesInCone(0,1,0,3.1415/8.0)
#    D.Snapshot.SelectParticlesInCone(0,0,1,3.1415/8.0)

    GridSphA = A.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
#    GridSphB = B.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
#    GridSphC = C.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
#    GridSphD = D.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)



    plt.plot(log10(GridSphA.R),GridSphA.Gamma,'-'+mar[i],label=r'Spherical',color=col[i],lw=1,ms=7,mew=1)

plt.xlim((-0.99,1))
plt.ylim((-3,-1))

#plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=3)

#    plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=4)
#plt.grid()
plt.ylabel(r'$\beta$',fontsize=24)
SetLabels(1,20)
#    plt.annotate(r'$r_{-2}$', xy=(-0.345, 0.5), xytext=(-0.345, 0.65),fontsize=26,ha='center',arrowprops=dict(facecolor='grey', shrink=0.05))


ax = plt.subplot(3,2,2)
ax.text(-0.8, 0.9, r'II) Energy exchange', fontsize=20)


SIMULATION_DIR = '/home/ms/Uni/DarkMatter/AllSimulations/HJS2010/'
Files = [SIMULATION_DIR+'00-515_001',SIMULATION_DIR+'s2_16',SIMULATION_DIR+'0E21_001',SIMULATION_DIR+'omE21_001',SIMULATION_DIR+'s3_16',SIMULATION_DIR+'s1_16',SIMULATION_DIR+'s4_16',SIMULATION_DIR+'om0-3.5E18_001']

for i in range(len(Files)):
    if i!=7:
        continue    
    if ShowSim[1] == False:
        break
    FileName = Files[i]
#    if col[i] not in ['green','red']:
#        continue
    
    

    A = DM_structure.DM_structure(FileName)
    A.FindCenter()
    A.FindCenterVel()
    A.CenterParticlePositions()
    A.CenterParticleVelocities()
    GridSph = A.CreateGridLogBins(NBins=100,Rmin=0.001)


#    B = copy.deepcopy(A)
#    C = copy.deepcopy(A)
#    D = copy.deepcopy(A)


#    B.Snapshot.SelectParticlesInCone(1,0,0,3.1415/8.0)
#    C.Snapshot.SelectParticlesInCone(0,1,0,3.1415/8.0)
#    D.Snapshot.SelectParticlesInCone(0,0,1,3.1415/8.0)

    GridSphA = A.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
#    GridSphB = B.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
#    GridSphC = C.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
#    GridSphD = D.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)



    plt.plot(log10(GridSphA.R),GridSphA.Gamma,'-'+mar[i],label=r'Spherical',color=col[i],lw=1,ms=7,mew=1)
#    plt.plot(log10(GridSphB.R),GridSphB.Gamma,'-'+mar[i],label=r'Spherical',color=col[i],lw=2,ms=9,mew=2)
#    plt.plot(log10(GridSphC.R),GridSphC.Gamma,'-'+mar[i],label=r'Spherical',color=col[i],lw=2,ms=9,mew=2)
#    plt.plot(log10(GridSphD.R),GridSphD.Gamma,'-'+mar[i],label=r'Spherical',color=col[i],lw=2,ms=9,mew=2)    
    
plt.xlim((-0.99,1))
plt.ylim((-3,-1))

#plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=3)

#    plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=4)
#plt.grid()
#plt.xlabel(r'$\log r$',fontsize=24)
#    plt.ylabel(r'$\beta$',fontsize=24)
SetLabels(1,1)
#    plt.annotate(r'$r_{-2}$', xy=(-0.345, 0.5), xytext=(-0.345, 0.65),fontsize=26,ha='center',arrowprops=dict(facecolor='grey', shrink=0.05))



ax = plt.subplot(3,2,3)
ax.text(-0.8, 0.9, r'III) Cold collapse', fontsize=20)

if ShowSim[2] == True:
    FileName = '/home/ms/Uni/DarkMatter/AllSimulations/InhomogeneousInfall/Infall_sub_1e6.bin_201'    

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

    plt.plot(log10(GridSphA.R),GridSphA.Gamma,'-x',label=r'Spherical',color='black',lw=2,ms=9,mew=2)
    plt.plot(log10(GridSphB.R),GridSphB.Gamma,'-o',label=r'$x$-axis',color='blue',lw=2,ms=9,mew=2)
    plt.plot(log10(GridSphC.R),GridSphC.Gamma,'-D',label=r'$y$-axis',color='red',lw=2,ms=9,mew=2)
    plt.plot(log10(GridSphD.R),GridSphD.Gamma,'-<',label=r'$z$-axis',color='green',lw=2,ms=9,mew=2)


plt.xlim((-0.99,1))
plt.ylim((-3,-1))
plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=4)


#plt.grid()
#plt.xlabel(r'$\log r$',fontsize=24)

SetLabels(1,20)
plt.ylabel(r'$\beta$',fontsize=24)

ax = plt.subplot(3,2,4)
ax.text(-0.8, 0.9, r'IV) Major mergers', fontsize=20)

if ShowSim[3] == True: 

    FileName = '/home/ms/Uni/DarkMatter/AllSimulations/CloneMerge/Mergers7_HQA_081'

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

    plt.plot(log10(GridSphA.R),GridSphA.Gamma,'-x',label=r'Spherical',color='black',lw=2,ms=9,mew=2)
    plt.plot(log10(GridSphB.R),GridSphB.Gamma,'-o',label=r'$x$-axis',color='blue',lw=2,ms=9,mew=2)
    plt.plot(log10(GridSphC.R),GridSphC.Gamma,'-D',label=r'$y$-axis',color='red',lw=2,ms=9,mew=2)
    plt.plot(log10(GridSphD.R),GridSphD.Gamma,'-<',label=r'$z$-axis',color='green',lw=2,ms=9,mew=2)


plt.xlim((-0.99,1))
plt.ylim((-3,-1))

plt.legend(prop=dict(size=14), numpoints=2, ncol=2,frameon=True,loc=4)

#plt.grid()

#plt.xlabel(r'$\log r$',fontsize=24)
#plt.ylabel(r'$\beta$',fontsize=24)

SetLabels(1,1)




ax = plt.subplot(3,2,5)
ax.text(-0.5, -0.3, r'V) Unstable model', fontsize=20)

if ShowSim[4] == True: 


    FileName = '/home/ms/Uni/DarkMatter/AllSimulations/ROI/OM_ROI00_rAN0.2_HQ_000'

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

    plt.plot(log10(GridSphA.R),GridSphA.Gamma,'-*',label=r'IC',color='grey',lw=2,ms=9,mew=2)
#    plt.plot(log10(GridSphB.R),GridSphB.Gamma,'-o',color='grey',lw=2,ms=9,mew=2)
#    plt.plot(log10(GridSphC.R),GridSphC.Gamma,'-D',color='grey',lw=2,ms=9,mew=2)
#    plt.plot(log10(GridSphD.R),GridSphD.Gamma,'-<',color='grey',lw=2,ms=9,mew=2)
    
    FileName = '/home/ms/Uni/DarkMatter/AllSimulations/ROI/OM_ROI00_rAN0.2_HQ_081'
    #FileName = '/home/ms/Uni/DarkMatter/AllSimulations/ROI/ROI05_120'
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

    plt.plot(log10(GridSphA.R),GridSphA.Gamma,'-x',label=r'Spherical',color='black',lw=2,ms=9,mew=2)
    plt.plot(log10(GridSphB.R),GridSphB.Gamma,'-o',label=r'$x$-axis',color='blue',lw=2,ms=9,mew=2)
    plt.plot(log10(GridSphC.R),GridSphC.Gamma,'-D',label=r'$y$-axis',color='red',lw=2,ms=9,mew=2)
    plt.plot(log10(GridSphD.R),GridSphD.Gamma,'-<',label=r'$z$-axis',color='green',lw=2,ms=9,mew=2)


plt.xlim((-0.99,1))
plt.ylim((-3,-1))

plt.legend(prop=dict(size=16), numpoints=2, ncol=1,frameon=True,loc=4)

#plt.grid()
plt.ylabel(r'$\beta$',fontsize=24)
plt.xlabel(r'$\log r$',fontsize=24)
SetLabels(20,20)



ax = plt.subplot(3,2,6)
ax.text(-0.8, 0.9, 'VI) Substructure', fontsize=20)

if ShowSim[5] == True: 

    FileName = '/home/ms/Uni/DarkMatter/AllSimulations/ViaLactea/ViaLactea01_1e6_000'

    A = DM_structure.DM_structure(FileName)
    A.FindCenter()
    A.FindCenterVel()
    A.CenterParticlePositions()
    A.CenterParticleVelocities()
    GridSph = A.CreateGridLogBins(NBins=100,Rmin=0.001)

    GridSphA = A.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)

    plt.plot(log10(GridSphA.R),GridSphA.Gamma,'-o',label=r'IC',color='grey',lw=2,ms=9,mew=2)
    
    FileName = '/home/ms/Uni/DarkMatter/AllSimulations/ViaLactea/ViaLactea01_1e6_101'

    A = DM_structure.DM_structure(FileName)
    A.FindCenter()
    A.FindCenterVel()
    A.CenterParticlePositions()
    A.CenterParticleVelocities()
    GridSph = A.CreateGridLogBins(NBins=100,Rmin=0.001)

    GridSphA = A.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)

    plt.plot(log10(GridSphA.R),GridSphA.Gamma,'-x',label=r'Final',color='black',lw=2,ms=9,mew=2)

plt.xlim((-0.99,1))
plt.ylim((-3,-1))

plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=1)

#plt.grid()

plt.xlabel(r'$\log r$',fontsize=24)
SetLabels(20,1)



plt.show()