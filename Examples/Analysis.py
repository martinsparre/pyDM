#/usr/bin/python
import matplotlib.pyplot as plt
import ReadGadget2
import DM_structure
from scipy import log10,sqrt,arccos,arcsin
import scipy

from math import pi
def f(E):
    return 1.0/sqrt(2.0) / (2.0*3.14159265359)**3. * sqrt(E) / (1.0-E)**2. * (  (1.-2.*E)* (8.*E**2. - 8.*E    -3.0 ) + 3.*arcsin(sqrt(E))/sqrt(E*(1.-E))      )

def g(E):
    A = 1.0/ sqrt(E*E)
    return (4.*pi)**2*sqrt(2.0*sqrt(E*E)) *(sqrt(A-1.) * (1./8.0 *A**2 - 5.0/12.0 *A - 1.0/3.0   )+ 1.0/8.0 *A * (A**2.-4.*A+8.)*arccos(1./sqrt(A)))

def N(E):
    return f(E)*g(E)

FileName='ISOHQ_001_LOWSOFT_000' #'OM_ROI00_rAN1_HQ_000'
A = DM_structure.DM_structure(FileName)
#Elin = A.CreateEnergyGridLinear(NBins=200)

#Grid = A.CreateGrid(Method='Radial')
GridSph = A.CreateGridLogBins(NBins=100,Rmin=0.001)
#Egrid = A.CreateEnergyGrid()
#V0 = A.GetLowestPotential()
#print  A.Snapshot.vx

#plt.subplot(2,3,1)

#x=scipy.linspace(-1.0,0.0,100)
##plt.plot(x,log10(g(x)),'-',color='black')
##plt.plot(x,log10(f(-x)),'--',color='black')
#plt.plot(-Egrid.E/V0,log10(-Egrid.Mass * V0),'o-',color='blue',alpha=0.5)
#plt.plot(-Elin.E/V0,log10(-Elin.Mass * V0),'<-',color='red',alpha=0.5)
#plt.plot(x,log10(f(-x)*g(x)),'-x',color='green')
##plt.ylim(-1.8,1.0)
#plt.xlim(-1.1,0.1)
#plt.xlabel('E/V0')
#plt.ylabel('N(E)')
#plt.grid()


#def f(x):
    #return 1.0/2.0/3.14159265359 * 1./x/(1.0+x)**3

#x=scipy.arange(0.0,100,0.0001)

#plt.subplot(2,3,2)
#plt.plot(log10(Grid.R),log10(Grid.R**2*Grid.Rho),'o')


#plt.plot(log10(GridSph.R),log10(GridSph.R**2*GridSph.Rho),'<',alpha=0.5)
#plt.plot(log10(x),log10(x**2*f(x)),'-',color='black')
#plt.xlabel('log r')
#plt.ylabel('log r*r*rho')
#plt.xlim(-1.5,2.2)

#plt.subplot(2,3,3)
#plt.xlabel('log r')
#plt.ylabel('log sigma2r')
#plt.plot(log10(GridSph.R),log10(GridSph.Sigma2),'-.',color='red')
#plt.plot(log10(Grid.R),log10(Grid.Sigma2),'--',color='black')
#plt.plot(log10(GridSph.R),log10(GridSph.Sigma2r),'-.',color='red')
#plt.plot(log10(Grid.R),log10(Grid.Sigma2r),'-',color='black')
#plt.xlim(-1.5,2.2)

#plt.subplot(2,3,4)

#plt.xlabel('log r')
#plt.ylabel('beta,gamma,kappa')
#plt.plot(log10(GridSph.R),GridSph.Beta,'-<',color='red',alpha=0.5)
#plt.plot(log10(Grid.R),Grid.Beta,'-o',color='black')
##plt.plot(log10(x), 1./(1.0+1.0/x**2),'--',color='black')
#dx=0.001
#plt.plot(log10(x),x/f(x) * (f(x+dx/2.0)-f(x-dx/2.0))/dx   ,'-',color='grey')
#plt.plot(log10(GridSph.R),GridSph.Gamma  ,'-<',alpha=0.5,color='red')
#plt.plot(log10(Grid.R),Grid.Gamma  ,'-o',color='black')
#plt.plot(log10(GridSph.R),GridSph.Kappa  ,'-v',color='blue',alpha=0.5)
#plt.plot(log10(Grid.R),Grid.Kappa  ,'-x',color='orange')
#plt.xlim(-1.5,2.2)
#plt.ylim(-6,1.5)

#plt.subplot(2,3,5)
#plt.plot(log10(Grid.Rmax),Grid.CumulativeMass  ,'-o',color='black')
#plt.plot(log10(GridSph.Rmax),GridSph.CumulativeMass  ,'->',color='red',alpha=0.5)
#plt.xlim(-1.5,2.2)
#plt.xlabel('log r')
#plt.ylabel('Cumulative Mass')
#plt.show()



plt.plot(log10(GridSph.R),GridSph.Gamma  ,'-',color='orange')

GridSph.RemoveBadRegions()
#GridSph.Smooth(Log10Width=0.2)
GridSph.RemoveBadRegions()
GridSph.SaveToBinary('GrSph.bin')
GridSph.SaveToAscii('GrSph.txt')
plt.plot(log10(GridSph.R),GridSph.Gamma ,'-',color='red')


plt.xlim(-3.1,2.1)
plt.ylim(-5,1)
plt.xlabel('log r')
plt.ylabel('')
plt.plot()
plt.show()
plt.savefig('a0000.ps')
#plt.show()
