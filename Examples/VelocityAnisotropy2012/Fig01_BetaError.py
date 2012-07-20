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





    
LogRMinus2 = -1.0    

FileName = '/home/ms/Uni/DarkMatter/AllSimulations/ViaLactea/ViaLactea01_1e6_101'
Cone=False

A = DM_structure.DM_structure(FileName)
A.FindCenter()
A.FindCenterVel()
A.CenterParticlePositions()
A.CenterParticleVelocities()
#GridSph = A.CreateGridLogBins(NBins=100,Rmin=0.001)

B = copy.deepcopy(A)

B.Snapshot.SelectParticlesInCone(1,0,0,3.1415/8.0)

GridSphA = A.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
GridSphB = B.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)

plt.plot(log10(GridSphA.R),GridSphA.Beta,'-x',label=r'Spherical',color='black',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphB.R),GridSphB.Beta,'-o',label=r'$x$-axis',color='blue',lw=2,ms=9,mew=2)



#bootstrap to get uncertainty:

Beta0 = []
BetaPlus = []
BetaMinus = []
for i in range(100):
    if Cone == True:
        GoodIds = numpy.random.randint(0,high=len(B.Snapshot.x),size=len(B.Snapshot.x))
        E=copy.deepcopy(B)
    else:
        GoodIds = numpy.random.randint(0,high=len(A.Snapshot.x),size=len(A.Snapshot.x))
        E=copy.deepcopy(A)
    Snap = E.Snapshot
    Snap.x = Snap.x[GoodIds]
    Snap.y = Snap.y[GoodIds]
    Snap.z = Snap.z[GoodIds]
    Snap.vx = Snap.vx[GoodIds]
    Snap.vy = Snap.vy[GoodIds]
    Snap.vz = Snap.vz[GoodIds]
    Snap.m = Snap.m[GoodIds]
    Snap.ID = Snap.ID[GoodIds]
    Snap.V = Snap.V[GoodIds]
    Snap.NPartTotal = len(Snap.x)
    Snap.NPart = [0,len(Snap.x),0,0,0,0]         
    
    GridSphE = E.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
    
    
    
    ID0 = scipy.fabs(log10(GridSphE.R) - LogRMinus2).argmin()
    IDPlus = scipy.fabs(log10(GridSphE.R) - (LogRMinus2+0.5) ).argmin()    
    IDMinus = scipy.fabs(log10(GridSphE.R) - (LogRMinus2-0.5) ).argmin()
    BetaPlus.append(GridSphE.Beta[IDPlus])
    Beta0.append(GridSphE.Beta[ID0])
    BetaMinus.append(GridSphE.Beta[IDMinus])
#    print log10(GridSphE.R[IDMinus]),log10(GridSphE.R[ID0]),log10(GridSphE.R[IDPlus])
    


print 'BetaMinus',scipy.mean(BetaMinus),scipy.std(BetaMinus)
print 'Beta0',scipy.mean(Beta0),scipy.std(Beta0)
print 'BetaPlus',scipy.mean(BetaPlus),scipy.std(BetaPlus)
print '--------------------'
print log10(GridSphE.R[IDMinus]),scipy.std(BetaMinus)
print log10(GridSphE.R[ID0]),scipy.std(Beta0)
print log10(GridSphE.R[IDPlus]),scipy.std(BetaPlus)
print '--------------------'

plt.errorbar(log10(GridSphE.R[IDMinus]),0.6,yerr=scipy.std(BetaMinus),lw=3,elinewidth=3,color='black')
plt.errorbar(log10(GridSphE.R[ID0]),0.6,yerr=scipy.std(Beta0),lw=3,elinewidth=3,color='black')
plt.errorbar(log10(GridSphE.R[IDPlus]),0.6,yerr=scipy.std(BetaPlus),lw=3,elinewidth=3,color='black')
plt.xlim((-0.99,1))
plt.ylim((-0.6,1.1))

plt.legend(prop=dict(size=14), numpoints=2, ncol=2,frameon=True,loc=4)

plt.show()