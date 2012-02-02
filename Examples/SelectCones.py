#/usr/bin/python
import matplotlib.pyplot as plt
import Classes.Gadget2 as Gadget2
import Classes.DM_structure as DM_structure
from scipy import log10
import scipy
from math import pi
import copy


    

FileName='NoSubTest01_1e6_000'
A = DM_structure.DM_structure(FileName)
GridSph = A.CreateGridLogBins(NBins=100,Rmin=0.001)


B = copy.deepcopy(A)
C = copy.deepcopy(A)
D = copy.deepcopy(A)
E = copy.deepcopy(A)

B.Snapshot.SelectParticlesInCone(1,1,0,3.1415/16.0)
C.Snapshot.SelectParticlesInCone(1.0,-1,0,3.1415/16.0)
D.Snapshot.SelectParticlesInCone(-1,-1,0,3.1415/16.0)
E.Snapshot.SelectParticlesInCone(-1,1,0,3.1415/16.0)

plt.plot(B.Snapshot.x,B.Snapshot.y,'.',label='$(1,1,0)$ $\pi/4$')
plt.plot(C.Snapshot.x,C.Snapshot.y,'.',label='$(1,-1,0)$ $\pi/16$')
plt.plot(D.Snapshot.x,D.Snapshot.y,'.',label='$(-1,-1,0)$ $\pi/16$')
plt.plot(E.Snapshot.x,E.Snapshot.y,'.',label='$(-1,1,0)$ $\pi/16$')
plt.legend(loc=2)
plt.grid()
plt.show()

