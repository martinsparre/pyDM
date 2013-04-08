#/usr/bin/python
import matplotlib.pyplot as plt
import Classes.Gadget2 as Gadget2
import Classes.DM_structure as DM_structure
from scipy import log10,sqrt,arccos,arcsin
import scipy
import scipy.optimize, scipy.interpolate
import pyROOT_functions
from math import pi
import numpy


def f(E):
    E=-E
    return 1.0/sqrt(2.0) / (2.0*3.14159265359)**3. * sqrt(E) / (1.0-E)**2. * (  (1.-2.*E)* (8.*E**2. - 8.*E    -3.0 ) + 3.*arcsin(sqrt(E))/sqrt(E*(1.-E))      )

def g(E):
    A = 1.0/ sqrt(E*E)
    return (4.*pi)**2*sqrt(2.0*sqrt(E*E)) *(sqrt(A-1.) * (1./8.0 *A**2 - 5.0/12.0 *A - 1.0/3.0   )+ 1.0/8.0 *A * (A**2.-4.*A+8.)*arccos(1./sqrt(A)))

def N(E):
    return f(E)*g(E)

FileName='Hernquist100000_000'
FileName='/home/ms/Uni/DarkMatter/AllSimulations/ViaLactea/NoSubTest01_1e6_000'
FileNames=['Hq1e5.bin']


for FileName in FileNames:
    A = DM_structure.DM_structure(FileName)

    plt.subplot(2,2,1)
    v2 = A.Snapshot.vx*A.Snapshot.vx+A.Snapshot.vy*A.Snapshot.vy+A.Snapshot.vz*A.Snapshot.vz
    E = 0.5*v2 + A.Snapshot.V
    HistVal,HistE=numpy.histogram(E,density=1,bins=100)
    BinCenter = (HistE[:-1]+HistE[1:])/2
    plt.plot(BinCenter,log10(HistVal),'-')
    N0 = HistVal[-1]


    plt.subplot(2,2,2)

    IDs = scipy.where((E>-0.9)*(E<-0.8))


    r = scipy.sqrt(A.Snapshot.x*A.Snapshot.x+A.Snapshot.y*A.Snapshot.y+A.Snapshot.z*A.Snapshot.z)
    v2 = A.Snapshot.vx*A.Snapshot.vx+A.Snapshot.vy*A.Snapshot.vy+A.Snapshot.vz*A.Snapshot.vz
    vr = 1.0/(r+1.0e-16)*scipy.array(A.Snapshot.vx*A.Snapshot.x + A.Snapshot.vy*A.Snapshot.y + A.Snapshot.vz*A.Snapshot.z)
    vtheta2 = v2-vr**2 
    L2 = vtheta2*r**2
    L2 = L2[IDs]

    HistVal,HistL2=numpy.histogram(L2/L2.max(),density=1,bins=20)
    BinCenter = (HistL2[:-1]+HistL2[1:])/2
    HistVal /= HistVal[scipy.argmin(scipy.fabs(BinCenter-0.45))]
    plt.plot(BinCenter,HistVal,'-')
    

    

    plt.subplot(2,2,3)
    v2 = A.Snapshot.vx*A.Snapshot.vx+A.Snapshot.vy*A.Snapshot.vy+A.Snapshot.vz*A.Snapshot.vz
    E = 0.5*v2 + A.Snapshot.V
    HistVal,HistE=numpy.histogram(E,density=1,bins=100)
    BinCenter = (HistE[:-1]+HistE[1:])/2
    plt.plot(BinCenter,log10(HistVal),'-')
#    pyROOT_functions.RootFit(BinCenter,HistVal,'[0]*exp( - [1]*(x-1.0) -1   )',Guess={},Limits={}, Fixed={}, FittingRange=[]):
    Guess={}
    Guess[1] = -2.5
    tmp = pyROOT_functions.RootFitErrors(BinCenter,HistVal,BinCenter*0,HistVal*0.001,'[0]* (   exp(-[1]*(x+1.0)) -1   )',Guess=Guess)
    E=scipy.linspace(-0.999,-0.001,100)
    A=tmp[0][0]
    Beta = tmp[0][1]
    print A,Beta
    plt.plot(E,log10(A * ( scipy.exp(-Beta*(E+1.0)) - 1.0  ))  ) 

    print tmp[1]

    


plt.subplot(2,2,1)
E=scipy.linspace(-0.999,-0.001,100)

plt.plot(E,log10(N(E)),'-',color='black')



plt.subplot(2,2,1)

plt.xlabel(r'$E/|\Phi (0)|$',fontsize=16)
plt.ylabel(r'$N(E/|\Phi (0)|)$',fontsize=16)

plt.subplot(2,2,2)

plt.xlabel(r'$L^2/L^2_\max$',fontsize=16)
plt.ylabel(r'$N(L^2)$ (Normalized at $0.45L^2_\max$)',fontsize=16)

plt.subplot(2,2,3)

plt.xlabel(r'$E/|\Phi_0|$',fontsize=16)
plt.ylabel(r'$N(E)$',fontsize=16)






plt.show()


#GridSph = A.CreateGridLogBins()



#x=scipy.logspace(-2.0,2.0,100)
#plt.show()

