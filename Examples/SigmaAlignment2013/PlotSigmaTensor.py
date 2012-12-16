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
        

FileName = '/home/ms/Uni/DarkMatter/AllSimulations/SigmaAlignment2013/HeadonMerger_VaryingVel/0.5/1HqIso_Impact0_160'

A = DM_structure.DM_structure(FileName)
A.FindCenter()
A.FindCenterVel()
A.CenterParticlePositions()
A.CenterParticleVelocities()
#GridSph = A.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50,CalcUncBeta=True)



B = copy.deepcopy(A)
B.Snapshot.SelectParticlesInCone(1.0,1,0.0,3.1415/8.0)
R = scipy.sqrt(B.Snapshot.x**2+B.Snapshot.y**2+B.Snapshot.z**2)
B.Snapshot.SelectParticles((R>5) * (R<6))

Tensor = scipy.zeros((3,3))
Tensor[0,0] = scipy.mean(B.Snapshot.vx * B.Snapshot.vx) - scipy.mean(B.Snapshot.vx)  *scipy.mean( B.Snapshot.vx)
Tensor[1,1] = scipy.mean(B.Snapshot.vy * B.Snapshot.vy) - scipy.mean(B.Snapshot.vy)  *scipy.mean( B.Snapshot.vy)               
Tensor[2,2] = scipy.mean(B.Snapshot.vz * B.Snapshot.vz) - scipy.mean(B.Snapshot.vz)  *scipy.mean( B.Snapshot.vz)
Tensor[0,1] = scipy.mean(B.Snapshot.vx * B.Snapshot.vy) - scipy.mean(B.Snapshot.vx)  *scipy.mean( B.Snapshot.vy)
Tensor[1,0] = Tensor[0,1] 
Tensor[0,2] = scipy.mean(B.Snapshot.vx * B.Snapshot.vz) - scipy.mean(B.Snapshot.vx)  *scipy.mean( B.Snapshot.vz)
Tensor[2,0] = Tensor[0,2]  
Tensor[1,2] = scipy.mean(B.Snapshot.vy * B.Snapshot.vz) - scipy.mean(B.Snapshot.vy)  *scipy.mean( B.Snapshot.vz)
Tensor[2,1] = Tensor[1,2]



#Chi2List = []
Ncones = 0


SigmaCones = []

AngleFile = open('Examples/SigmaAlignment2013/192.txt','r')
#AngleFile = open('../192.txt','r')
for line in AngleFile:
    tmp = line.split()
    if len(tmp) != 5:
        continue
    
    x = float(tmp[0])
    y = float(tmp[1])
    z = float(tmp[2])
    
    r = scipy.matrix(scipy.array([x,y,z]))
    T=scipy.matrix(Tensor)
    Sigma = r*T*r.transpose()

    SigmaCones.append(Sigma[0,0])

    Ncones += 1
    print 'Cone number ',Ncones

AngleFile.close()



import numpy as np
import healpy as hp

NSIDE = 4
m = scipy.array(SigmaCones)

hp.mollview(m, title=r"The $\sigma^2$-ellipsoid at (2.5,2.5,0.0)")
plt.text(0,0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(-1.8,0.0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(1.8,0.0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)

#plt.text(1,0,'Halo\ncenter',horizontalalignment='center',verticalalignment='center',fontsize=16)
#plt.text(-1,0,'Halo\ncenter',horizontalalignment='center',verticalalignment='center',fontsize=16)

#plt.text(0,0.9,'Minor axis',horizontalalignment='center',verticalalignment='center',fontsize=16)
#plt.text(0,-0.9,'Minor axis',horizontalalignment='center',verticalalignment='center',fontsize=16)

plt.show()
