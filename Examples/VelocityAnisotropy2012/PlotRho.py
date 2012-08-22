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
import scipy.interpolate,scipy.optimize

def FindNearestIDs(x,y,y0):
    ID=(abs(y-y0)).argmin()
    if ID<3 or ID > len(x)-3:
        return [x[ID]],[y[ID]]
    else:
        return x[ID-3:ID+3],y[ID-3:ID+3]

def DoInterpolate(x,y,x0):
    y1,x1 = FindNearestIDs(y,x,x0)
    if len(x1)>1:
        spline = scipy.interpolate.interp1d(x1,y1,kind='cubic')
        return spline(x0)
    else:
        return x1[0]

def DoSolve(x,y,y0):
    x1,y1 = FindNearestIDs(x,y,y0)
    if len(x1)>1:
        spline = scipy.interpolate.interp1d(x1,y1-y0,kind='cubic')
        return scipy.optimize.fsolve(spline,x1[len(x1)/2],full_output=True)
    else:
        return x1[0]


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
        
#FileName = '/home/ms/Uni/DarkMatter/SimulationDataAndAnalysis/ROI/OM_ROI00_rAN0.2_HQ_000'
#FileName = '/home/ms/Uni/DarkMatter/AllSimulations/MergerAnisotropy/1HqIso_Impact0_121'
#FileName = '/home/ms/Uni/DarkMatter/AllSimulations/MergerAnisotropy/1HqIso_Impact10_120'
#FileName = '/home/ms/Uni/DarkMatter/AllSimulations/CloneMerge/Mergers7_HQA_081'
#FileName = '/home/ms/Uni/DarkMatter/AllSimulations/MergerAnisotropy/1HQOM_Impact0_121'
FileName = '/home/ms/Uni/DarkMatter/AllSimulations/MergerAnisotropy/Hq_0.1_0.5_080'
#FileName = '/home/ms/Uni/DarkMatter/AllSimulations/MergerAnisotropy/Hq_0.1_0.5Impact_089'
#FileName = '/home/ms/Uni/DarkMatter/SimulationDataAndAnalysis/MergerAnisotropy2012/Simulations/PyDM1/1HqIso_000'

A = DM_structure.DM_structure(FileName)
A.FindCenter()
A.FindCenterVel()
A.CenterParticlePositions()
A.CenterParticleVelocities()
GridSph = A.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50,CalcUncBeta=True)

X = []
Y = []
Z = []
Rho26 = []

Cones = []
BetaCones = []

Ncones = 0
AngleFile = open('../12.txt','r')
#AngleFile = open('../192.txt','r')
for line in AngleFile:
    tmp = line.split()
    if len(tmp) != 5:
        continue
    
    x = float(tmp[0])
    y = float(tmp[1])
    z = float(tmp[2])

    B = copy.deepcopy(A)
    B.Snapshot.SelectParticlesInCone(x,y,z,3.1415/8.0)

    GridSphB = B.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50,CalcUncBeta=True)

    plt.plot(log10(GridSphB.R),log10(GridSphB.Rho),'-',color='black',lw=1)
    Cones.append(GridSphB)
    Ncones += 1
    Rho26.append(log10(GridSphB.Rho[26]))
    X.append(x)
    Y.append(y)
    Z.append(z)
    print Ncones,log10(GridSphB.Rho[26])


    
    
    
B = copy.deepcopy(A)
B.Snapshot.SelectParticlesInCone(1,0,0,3.1415/8.0)

GridSphB = B.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50,CalcUncBeta=True)

plt.plot(log10(GridSphB.R),log10(GridSphB.Rho),'-',color='red',lw=2)

B = copy.deepcopy(A)
B.Snapshot.SelectParticlesInCone(-1,0,0,3.1415/8.0)

GridSphB = B.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50,CalcUncBeta=True)

plt.plot(log10(GridSphB.R),log10(GridSphB.Rho),'-',color='blue',lw=2)


ID = scipy.argmin(Rho26)
print 'Min:',X[ID],Y[ID],Z[ID],Rho26[ID]

ID = scipy.argmax(Rho26)
print 'Max:',X[ID],Y[ID],Z[ID],Rho26[ID]

plt.show()


OutputFile = open('xyz.txt','w+')
ids=scipy.where((A.Snapshot.V>-0.8)*(A.Snapshot.V<-0.7))

x = A.Snapshot.x[ids]
y = A.Snapshot.y[ids]
z = A.Snapshot.z[ids]
for i in range(len(x)):
    OutputFile.write(str(x[i])+'\t'+str(y[i])+'\t'+str(z[i])+'\n')
OutputFile.close()


sys.exit()
for rho in log10(GridSphB.Rho)[15:-15:5]:
    r = []
    for cone in Cones:
        
        tmp = DoSolve(log10(cone.R),log10(cone.Rho),rho)
        
        if 'converged' in tmp[3]:
            r.append(tmp[0])
            print tmp[3]
    
    

    r = scipy.array(r)
    AxisRatio = 10.0**(r.max() - r.min())
    plt.plot(r,AxisRatio,'o')

plt.show()