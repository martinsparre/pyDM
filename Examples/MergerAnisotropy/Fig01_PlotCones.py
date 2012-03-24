#/usr/bin/python
import matplotlib.pyplot as plt
import Classes.Gadget2 as Gadget2
import Classes.DM_structure as DM_structure
from scipy import log10
import scipy
import pylab
from math import pi
import copy,sys

def Undersample(P, N=25000):
    "P: List with particle positions and others: [x,y,z, vx,vy,vz,V]"
    #Check that the attributes have the same length:

    if len(P)>0:
        Len = len(P[0])
    else:
        print 'P has length 0.'
        return P

    for i in range(len(P)):
        if len(P[i]) != Len:
            print 'P['+str(i)+'] doesnt have the same length as P[0]'
            return P
    
    #Choose points to sample:
    rand = scipy.random.random_integers(0,len(P[0])-1,N)

    newP = []
    for i in range(len(P)):
        newP.append(P[i][rand])

    return newP

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
        

FileName = '../1HqIso_Impact0_121'
A = DM_structure.DM_structure(FileName)
A.FindCenter()
A.FindCenterVel()
A.CenterParticlePositions()
A.CenterParticleVelocities()
#GridSph = A.CreateGridLogBins(NBins=100,Rmin=0.001)


B = copy.deepcopy(A)
C = copy.deepcopy(A)
D = copy.deepcopy(A)
E = copy.deepcopy(A)

B.Snapshot.SelectParticlesInCone(1,0,0,3.1415/8.0)
C.Snapshot.SelectParticlesInCone(0,1,0,3.1415/8.0)

Points = [B.Snapshot.x,B.Snapshot.y,B.Snapshot.z,B.Snapshot.V]
x,y,z,V = Undersample(Points,N=2000)
plt.plot(x,y,'<',label=r'Cone along collision axis',mew=0.0,color='blue',mec='blue',ms=8)
Points = [C.Snapshot.x,C.Snapshot.y,C.Snapshot.z,C.Snapshot.V]
x,y,z,V = Undersample(Points,N=2000)
plt.plot(x,y,'o',label=r'Cone perpendicular to collision axis',mew=0.0,color='red',mec='red',ms=8)

t=scipy.linspace(-10,10,20)
plt.plot(t,0*t,'--',lw=2,label='The collision axis',color='black')
plt.legend(loc=2)
plt.grid()
plt.xlabel('$x$',fontsize=24)
plt.ylabel('$y$',fontsize=24)
plt.xlim((-2,5))
plt.ylim((-2,5))
SetLabels()
plt.legend(prop=dict(size=16), numpoints=1, ncol=1,frameon=True,loc=1)
plt.show()






