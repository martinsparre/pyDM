#/usr/bin/python
import matplotlib.pyplot as plt
#import Classes.Gadget2 as ReadGadget2
import Classes.DM_structure as DM_structure
from scipy import log10,sqrt,arccos,arcsin
import scipy, scipy.interpolate, scipy.optimize
from math import pi
import pyROOT_functions

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


FileNamesShort = ['Hq1e6_rvir5e3_Beta0.001_025',
'Hq1e6_rvir5e3_Beta0.005.bin_025',
'Hq1e6_rvir5e3_Beta0.01.bin_025',
'Hq1e6_rvir5e3_Beta0.025.bin_025',
'Hq1e6_rvir5e3_Beta0.04.bin_025',
'Hq1e6_rvir5e3_Beta0.05.bin_025',
'Hq1e6_rvir5e3_Beta0.075.bin_025',
'Hq1e6_rvir5e3_Beta0.1.bin_025',
'Hq1e6_rvir5e3_Beta0.2.bin_025',
'Hq1e6_rvir5e3_Beta0.15A_025',
'Hq1e6_rvir5e3_Beta0.15A_000',
'Long/Hq1e6_rvir5e3_Beta0.001A_025',
'Long/Hq1e6_rvir5e3_Beta0.005A_025',
'Long/Hq1e6_rvir5e3_Beta0.01A_025',
'Long/Hq1e6_rvir5e3_Beta0.025A_025']

FileNames = FileNamesShort

DataDir = '/home/ms/Uni/DarkMatter/AllSimulations/Londrillo/'

for i in range(len(FileNames)):
    FileNames[i] = DataDir + FileNames[i]


Color = ['red', 'blue', 'green', 'maroon','orange','black','pink','yellow','cyan']*2
Label = [r'$\Delta t = 250$']*len(FileNames)
Beta = [0.001,0.005,0.01,0.025,0.04,0.05,0.075,0.1,0.2,0.15,0.15,0.001,0.005,0.01,0.025]

plt.subplots_adjust(left=0.06, bottom=0.15,  wspace=0.3,  hspace=0.3)

for i in range(len(FileNames)):
    if i>8:
        break
    A = DM_structure.DM_structure(FileNames[i])
    A.FindCenter()
    A.FindCenterVel()
    A.CenterParticlePositions()
    A.CenterParticleVelocities()
    A.Snapshot.SelectParticles(A.Snapshot.V+0.5*(A.Snapshot.vx**2+A.Snapshot.vy**2+A.Snapshot.vz**2)<0.0)
    
    CentralPotential = A.Snapshot.V.min()
    print CentralPotential

    GridSph = A.CreateGridLogBins(NBins=50,Rmin=0.001)
    GridSphInner = A.CreateGridLogBins(NBins=50,Rmin=0.005,Rmax=0.5)    
    
    logrhalf = DoSolve(log10(GridSph.R),GridSph.CumulativeMass, GridSph.CumulativeMass.max()*0.5)
    logr10 = DoSolve(log10(GridSph.R),GridSph.CumulativeMass, GridSph.CumulativeMass.max()*0.1)
    LondrilloPsi = CentralPotential*10.0**logrhalf/GridSph.CumulativeMass.max()
    if i>8:
        Color = 'red'
    else:
        Color = 'blue'
    plt.subplot(2,3,1)
    plt.plot(log10(Beta[i]),LondrilloPsi,'o',color=Color)

    plt.subplot(2,3,2)
    plt.plot(log10(Beta[i]),10**(logrhalf-logr10),'o',color=Color)
    
    plt.subplot(2,3,3)
    plt.plot(log10(Beta[i]),2*scipy.sum(0.5*(A.Snapshot.vx**2+A.Snapshot.vy**2+A.Snapshot.vz**2)) / (0.5*scipy.sum(A.Snapshot.V)),'o',color=Color)    
    
    plt.subplot(2,3,4)
    GammaSmooth = scipy.array(pyROOT_functions.RootSpline(log10(GridSphInner.R[5:-5]),GridSphInner.Gamma[5:-5],BandPass=0.3))[1]
    plt.plot(log10(GridSphInner.R[7:-7]),GridSphInner.Gamma[7:-7]+i,color='blue')
    plt.plot(log10(GridSphInner.R[7:-5]),GammaSmooth[2:]+i,color='red')    
    logrs = DoSolve(log10(GridSphInner.R[7:-5]),GammaSmooth[2:], -2)
    plt.plot(logrs,-2+i,'o',color='green')
    newr = scipy.array([logrs-0.3010])
    newgamma = scipy.array(pyROOT_functions.RootSpline(log10(GridSphInner.R[5:-5]),GridSphInner.Gamma[5:-5],x_out = newr,BandPass=0.3)[1])
    plt.plot(newr,newgamma+i,'o',color='orange')
    
    plt.text(-2.8,-4+1.3*i,'beta='+str(Beta[i]))
    
    plt.subplot(2,3,5)
    plt.plot(log10(Beta[i]),newgamma,'o',color=Color)
    
    plt.subplot(2,3,6)
    plt.plot(LondrilloPsi,newgamma,'o',color=Color)
    
    
plt.subplot(2,3,1)
plt.xlim((-3.1,-0.4))
plt.xlabel(r'$\log_{10} \beta\equiv\left(2T/W\right)_0$',fontsize=18)
plt.ylabel(r'$\psi \equiv r_{0.5}\Phi(0) / M$',fontsize=18)
plt.grid()

plt.subplot(2,3,2)
plt.xlim((-3.1,-0.4))
plt.xlabel(r'$\log_{10} \beta\equiv\left(2T/W\right)_0$',fontsize=18)
plt.ylabel(r'$r_{0.5} / r_{0.1}$',fontsize=18)
plt.grid()

plt.subplot(2,3,3)
plt.xlim((-3.1,-0.4))
plt.xlabel(r'$\log_{10} \beta\equiv\left(2T/W\right)_0$',fontsize=18)
plt.ylabel(r'$2T/W$',fontsize=18)
plt.grid()

plt.subplot(2,3,4)
plt.xlim((-3.1,-0.4))
plt.ylim((-6,12))
plt.xlabel(r'$\log_{10} r$',fontsize=18)
plt.ylabel(r'$\gamma\equiv d\log \rho / d\log r$',fontsize=18)
plt.grid()

plt.subplot(2,3,5)
plt.xlim((-3.1,-0.4))
plt.xlabel(r'$\log_{10} \beta\equiv\left(2T/W\right)_0$',fontsize=18)
plt.ylabel(r'$\gamma( 0.5 r_{-2} )$',fontsize=18)
plt.grid()

plt.subplot(2,3,6)
plt.xlabel(r'$\psi \equiv r_{0.5}\Phi(0) / M$',fontsize=18)
plt.ylabel(r'$\gamma( 0.5 r_{-2} )$',fontsize=18)
plt.grid()
plt.xlim((-4,-2.7))



plt.show()

