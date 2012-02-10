#/usr/bin/python
import matplotlib.pyplot as plt
import Classes.Gadget2 as Gadget2
import Classes.DM_structure as DM_structure
from scipy import log10
import scipy
from math import pi
import copy,sys



if len(sys.argv[1]):
    print 'Give filename, please'


FileName = sys.argv[1]
A = DM_structure.DM_structure(FileName)
A.FindCenter()
A.FindCenterVel()
A.CenterParticlePositions()
A.CenterParticleVelocities()
GridSph = A.CreateGridLogBins(NBins=100,Rmin=0.001)


B = copy.deepcopy(A)
C = copy.deepcopy(A)
D = copy.deepcopy(A)
E = copy.deepcopy(A)

B.Snapshot.SelectParticlesInCone(1,0,0,3.1415/8.0)
C.Snapshot.SelectParticlesInCone(0,1,0,3.1415/8.0)
D.Snapshot.SelectParticlesInCone(0,0,1,3.1415/8.0)
E.Snapshot.SelectParticlesInCone(-1,0,0,3.1415/8.0)


PlotParticles = False
if PlotParticles == True:
    plt.title('File'+sys.argv[1])
    plt.plot(B.Snapshot.x,B.Snapshot.y,'.',label='$(1,0,0)$ $\pi/4$')
    plt.plot(C.Snapshot.x,C.Snapshot.y,'.',label='$(0,1,0)$ $\pi/16$')
    plt.plot(D.Snapshot.x,D.Snapshot.y,'.',label='$(0,0,1)$ $\pi/16$')
    plt.plot(E.Snapshot.x,E.Snapshot.y,'.',label='$(-1,0,0)$ $\pi/16$')
    plt.legend(loc=2)
    plt.grid()
    plt.show()


GridSphA = A.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
GridSphB = B.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
GridSphC = C.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
GridSphD = D.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
GridSphE = E.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)


plt.subplot(2,3,1)
plt.plot(log10(GridSphA.R),GridSphA.Gamma,'-x',label='Spherical',color='black',lw=2)
plt.plot(log10(GridSphB.R),GridSphB.Gamma,'-o',label='x-axis',color='blue',lw=2)
plt.plot(log10(GridSphC.R),GridSphC.Gamma,'-D',label='y-axis',color='red',lw=2)
plt.plot(log10(GridSphD.R),GridSphD.Gamma,'-<',label='z-axis',color='green',lw=2)
plt.plot(log10(GridSphE.R),GridSphE.Gamma,'->',label='neg. x-axis',color='grey',lw=2)
plt.ylim((-5,-0.5))
plt.xlim((-1,1))
plt.legend(prop=dict(size=14), numpoints=2, ncol=1,frameon=True,loc='best')
plt.grid()
plt.xlabel('log r',fontsize=20)
plt.ylabel(r'$\gamma\equiv d\log \rho / d\log r$',fontsize=20)

plt.subplot(2,3,2)
plt.plot(log10(GridSphA.R),GridSphA.Beta,'-x',label='Spherical',color='black',lw=2)
plt.plot(log10(GridSphB.R),GridSphB.Beta,'-o',label='x-axis',color='blue',lw=2)
plt.plot(log10(GridSphC.R),GridSphC.Beta,'-D',label='y-axis',color='red',lw=2)
plt.plot(log10(GridSphD.R),GridSphD.Beta,'-<',label='z-axis',color='green',lw=2)
plt.plot(log10(GridSphE.R),GridSphE.Beta,'->',label='neg. x-axis',color='grey',lw=2)
plt.ylim((-0.2,0.5))
plt.xlim((-1,1))
plt.legend(prop=dict(size=14), numpoints=2, ncol=1,frameon=True,loc='best')
plt.grid()
plt.xlabel('log r',fontsize=20)
plt.ylabel(r'$\beta$',fontsize=20)


#plt.subplot(2,2,3)
#plt.plot(log10(GridSphA.R),log10(GridSphA.R**2*GridSphA.Rho),'-x',label='Spherical',color='black',lw=2)
#plt.plot(log10(GridSphB.R),log10(GridSphB.R**2*GridSphB.Rho),'-o',label='x-axis',color='blue',lw=2)
#plt.plot(log10(GridSphC.R),log10(GridSphC.R**2*GridSphC.Rho),'-D',label='y-axis',color='red',lw=2)
#plt.plot(log10(GridSphD.R),log10(GridSphD.R**2*GridSphD.Rho),'-<',label='z-axis',color='green',lw=2)
#plt.plot(log10(GridSphE.R),log10(GridSphE.R**2*GridSphE.Rho),'->',label='neg. x-axis',color='grey',lw=2)
#plt.ylim((-4.5,-1.0))
#plt.xlim((-1,1))
#plt.legend(prop=dict(size=14), numpoints=2, ncol=1,frameon=True,loc='best')
#plt.grid()
#plt.xlabel('log r',fontsize=20)
#plt.ylabel(r'$\log r^2\rho$',fontsize=20)


plt.subplot(2,3,3)
plt.plot(log10(GridSphA.R),log10(GridSphA.Sigma2-GridSphA.Sigma2r),'-x',label='Spherical',color='black',lw=2)
plt.plot(log10(GridSphB.R),log10(GridSphB.Sigma2-GridSphB.Sigma2r),'-o',label='x-axis',color='blue',lw=2)
plt.plot(log10(GridSphC.R),log10(GridSphC.Sigma2-GridSphC.Sigma2r),'-D',label='y-axis',color='red',lw=2)
plt.plot(log10(GridSphD.R),log10(GridSphD.Sigma2-GridSphD.Sigma2r),'-<',label='z-axis',color='green',lw=2)
plt.plot(log10(GridSphE.R),log10(GridSphE.Sigma2-GridSphE.Sigma2r),'->',label='neg. x-axis',color='grey',lw=2)

plt.ylim((-1.2,-0.2))
plt.xlim((-1,1))
plt.legend(prop=dict(size=14), numpoints=2, ncol=2,frameon=True,loc='best')
plt.grid()
plt.xlabel('log r',fontsize=20)
plt.ylabel(r'$\log \sigma^2_{tan}$',fontsize=20)


plt.subplot(2,3,4)
plt.plot(log10(GridSphA.R),log10(GridSphA.Sigma2r),'--x',label='Spherical',color='black',lw=2)
plt.plot(log10(GridSphB.R),log10(GridSphB.Sigma2r),'--o',label='x-axis',color='blue',lw=2)
plt.plot(log10(GridSphC.R),log10(GridSphC.Sigma2r),'--D',label='y-axis',color='red',lw=2)
plt.plot(log10(GridSphD.R),log10(GridSphD.Sigma2r),'--<',label='z-axis',color='green',lw=2)
plt.plot(log10(GridSphE.R),log10(GridSphE.Sigma2r),'-->',label='neg. x-axis',color='grey',lw=2)
plt.ylim((-1.5,-0.6))
plt.xlim((-1,1))
plt.legend(prop=dict(size=14), numpoints=2, ncol=2,frameon=True,loc='best')
plt.grid()
plt.xlabel('log r',fontsize=20)
plt.ylabel(r'$\log \sigma_r^2$',fontsize=20)


plt.subplot(2,3,5)


#plt.ylim((-1.2,-0.2))
plt.xlim((-1,1))
plt.legend(prop=dict(size=14), numpoints=2, ncol=2,frameon=True,loc='best')
plt.grid()
plt.xlabel('log r',fontsize=20)
plt.ylabel(r'$\log\overline{V}^2$',fontsize=20)


plt.subplot(2,3,6)
plt.plot(log10(GridSphA.R),log10(GridSphA.MeanVr**2),'-x',label='Spherical',color='black',lw=2)
plt.plot(log10(GridSphB.R),log10(GridSphB.MeanVr**2),'-o',label='x-axis',color='blue',lw=2)
plt.plot(log10(GridSphC.R),log10(GridSphC.MeanVr**2),'-D',label='y-axis',color='red',lw=2)
plt.plot(log10(GridSphD.R),log10(GridSphD.MeanVr**2),'-<',label='z-axis',color='green',lw=2)
plt.plot(log10(GridSphE.R),log10(GridSphE.MeanVr**2),'->',label='neg. x-axis',color='grey',lw=2)

#plt.ylim((-1.2,-0.2))
plt.xlim((-1,1))
plt.legend(prop=dict(size=14), numpoints=2, ncol=2,frameon=True,loc='best')
plt.grid()
plt.xlabel('log r',fontsize=20)
plt.ylabel(r'$\log \overline{V_r}^2$',fontsize=20)


plt.show()



