#/usr/bin/python
import matplotlib.pyplot as plt
#import Classes.Gadget2 as ReadGadget2
import Classes.DM_structure as DM_structure
from scipy import log10,sqrt,arccos,arcsin
import scipy

from math import pi


FileName='0GManySmall_28_001'
A = DM_structure.DM_structure(FileName)
A.FindCenter()
A.FindCenterVel()
A.CenterParticlePositions()
A.CenterParticleVelocities()
GridSph = A.CreateGridLogBins(NBins=100,Rmin=0.01)


plt.plot(log10(GridSph.R),GridSph.Gamma  ,'-',color='orange')



plt.xlim(-3.1,2.1)
plt.ylim(-5,1)
plt.xlabel('log r')
plt.ylabel('gamma')
plt.plot()
plt.show()

