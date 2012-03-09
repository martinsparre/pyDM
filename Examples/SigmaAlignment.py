#/usr/bin/python
import matplotlib.pyplot as plt
import Classes.Gadget2 as Gadget2
import Classes.DM_structure as DM_structure
from scipy import log10
import scipy
import pylab
from math import pi
import copy,sys


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
        

if len(sys.argv)<2:
    print 'Give filename, please'


FileName = sys.argv[1]
A = DM_structure.DM_structure(FileName)
A.FindCenter()
#A.FindCenterVel()
A.CenterParticlePositions()
#A.CenterParticleVelocities()
GridSph = A.CreateGridLogBins(Rmin=0.4,Rmax=5.0,NBins=20,CalcVelDispTensor = True)


B = copy.deepcopy(A)

Axis = scipy.array([1.0,0.0,0.0])
Axis /= scipy.sqrt(Axis[0]**2+Axis[1]**2+Axis[2]**2)
print Axis
#LenAxis = scipy.sqrt(Axis[0]**2+Axis[1]**2+Axis[2]**2)

B.Snapshot.SelectParticlesInCone(Axis[0],Axis[1],Axis[2],3.1415/16.0)
GridSphB = B.CreateGridLogBins(NBins=15,CalcVelDispTensor = True,Rmin=0.4,Rmax=5.0)




MinAngle = []
AngleMain = []
for i in range(len(GridSphB.SigmaTensor)):
    if GridSphB.N[i] < 15 or GridSphB.Beta[i]<0.2:
        MinAngle.append(-25)
        AngleMain.append(-25)
        continue
    SigmaTensor = GridSphB.SigmaTensor[i]
    EigenValue,EigenVector = scipy.linalg.eig(SigmaTensor)
    
    print EigenVector[0][0]**2+EigenVector[0][1]**2+EigenVector[0][2]**2
    

    CosAngle0 = abs( scipy.sum(EigenVector[0] * Axis) )
    CosAngle1 = abs( scipy.sum(EigenVector[1] * Axis) )
    CosAngle2 = abs( scipy.sum(EigenVector[2] * Axis) )
#    print CosAngle0, CosAngle1, CosAngle2
    Angle0 = scipy.arccos(CosAngle0)
    Angle1 = scipy.arccos(CosAngle1)
    Angle2 = scipy.arccos(CosAngle2)
    print CosAngle0,CosAngle1,CosAngle2
    
    MinAngle.append(min(abs(Angle0),abs(Angle1),abs(Angle2)))
    a = scipy.where(EigenValue == EigenValue.max())
 
    tmp = scipy.arccos(abs(  scipy.sum(EigenVector[a]*Axis)  ))
    if tmp > 3.1415/2.0:
        tmp=3.1415 - tmp
    AngleMain.append(tmp)


AngleMain = scipy.array(AngleMain)
MinAngle = scipy.array(MinAngle)

ID = scipy.where(MinAngle!=-25)

plt.subplot(1,3,1)
plt.plot(log10(GridSphB.R),GridSphB.Beta,'o')
plt.plot(log10(GridSph.R),GridSph.Beta,'-',color='black')
plt.ylim((-0.4,1.1))
plt.xlim((-0.4,1.0))
plt.xlabel('log r')
plt.ylabel('Beta')

plt.subplot(1,3,2)
plt.plot(log10(GridSphB.R[ID]),MinAngle[ID],'o')
plt.plot(log10(GridSphB.R[ID]),AngleMain[ID],'<')
plt.ylim((-0.1,3.2))
plt.xlim((-0.4,1.0))
plt.xlabel('log r')
plt.ylabel('Angle')

plt.subplot(1,3,3)
plt.plot(log10(GridSphB.R[ID]),GridSphB.N[ID],'o')
plt.ylim((200,800))
plt.xlim((-0.4,1.0))
plt.xlabel('log r')
plt.ylabel('N')

plt.show()
sys.exit()

