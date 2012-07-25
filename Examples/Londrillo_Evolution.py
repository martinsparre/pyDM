#/usr/bin/python
import matplotlib.pyplot as plt
#import Classes.Gadget2 as ReadGadget2
import Classes.DM_structure as DM_structure
from scipy import log10,sqrt,arccos,arcsin
import scipy, scipy.interpolate, scipy.optimize
from math import pi

def FindNearestIDs(x,y,y0):
    ID=(abs(y-y0)).argmin()
    if ID<3 or ID > len(x)-3:
        return [x[ID]],[y[ID]]
    else:
        return x[ID-3:ID+3],y[ID-3:ID+3]

def DoSolve(x,y,y0):
    x1,y1 = FindNearestIDs(x,y,y0)
    if len(x1)>1:
        spline = scipy.interpolate.interp1d(x1,y1-y0,kind='cubic')
        return scipy.optimize.fsolve(spline,x1[len(x1)/2])
    else:
        return x1[0]

FileNumber = scipy.array(range(50))
#FileNumber = FileNumber[::5]
FilePath = '/home/ms/Desktop/LondrilloAnalysis/pyDM/Data/0.15/'


#Color = ['black']*len(FileNames)
plt.subplots_adjust(left=0.06, bottom=0.15,  wspace=0.3,  hspace=0.3)
for i in FileNumber:
    A = DM_structure.DM_structure( FilePath + str(i).zfill(3) )
    A.FindCenter()
    A.FindCenterVel()
    A.CenterParticlePositions()
    A.CenterParticleVelocities()
    A.Snapshot.SelectParticles(A.Snapshot.V+0.5*(A.Snapshot.vx**2+A.Snapshot.vy**2+A.Snapshot.vz**2)<0.0)
    
    CentralPotential = A.Snapshot.V.min()
#    print CentralPotential

    GridSph = A.CreateGridLogBins(NBins=100,Rmin=0.001)
    
#    plt.subplot(1,2,1)
#    plt.plot(log10(GridSph.R),GridSph.CumulativeMass  ,'-',color=Color[i],label=Label[i])        
    
    logrhalf = DoSolve(log10(GridSph.R),GridSph.CumulativeMass, GridSph.CumulativeMass.max()*0.5)
    logr10 = DoSolve(log10(GridSph.R),GridSph.CumulativeMass, GridSph.CumulativeMass.max()*0.1)
    LondrilloPsi = CentralPotential*10.0**logrhalf/GridSph.CumulativeMass.max()
    plt.subplot(1,3,1)
    plt.plot(i,LondrilloPsi,'o',color='blue')

    plt.subplot(1,3,2)
    plt.plot(i,10**(logrhalf-logr10),'o',color='blue')
    
    plt.subplot(1,3,3)
    plt.plot(i,2*scipy.sum(0.5*(A.Snapshot.vx**2+A.Snapshot.vy**2+A.Snapshot.vz**2)) / (0.5*scipy.sum(A.Snapshot.V)),'o',color='blue')    

plt.subplot(1,3,1)
#plt.xlim((-3.1,-0.4))
#plt.ylim(-0.05,1.2)
plt.subplot(1,3,1)
#plt.xlabel(r'$\log_{10} \beta\equiv\left(2T/W\right)_0$',fontsize=18)
plt.ylabel(r'$\psi \equiv r_{0.5}\Phi(0) / M$',fontsize=18)
#plt.legend(loc=2)
plt.grid()

plt.subplot(1,3,2)
#plt.xlim((-3.1,-0.4))
#plt.xlabel(r'$\log_{10} \beta\equiv\left(2T/W\right)_0$',fontsize=18)
plt.ylabel(r'$r_{0.5} / r_{0.1}$',fontsize=18)
plt.grid()

plt.subplot(1,3,3)
#plt.xlim((-3.1,-0.4))
#plt.xlabel(r'$\log_{10} \beta\equiv\left(2T/W\right)_0$',fontsize=18)
plt.ylabel(r'$2T/W$',fontsize=18)
plt.grid()

plt.show()

