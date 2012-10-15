#/usr/bin/python
import matplotlib.pyplot as plt
import Classes.Gadget2 as Gadget2
import Classes.DM_structure as DM_structure
from scipy import log10,cos,sin
import scipy
import pylab
from math import pi
import copy,sys
import numpy.random
import scipy.interpolate,scipy.optimize


def GetPrincipalAxes(Angle1,Angle2,Angle3):
    "Input: the three euler angles from Tipsy. Output: the three axes..."

    pi = 3.14159265359
    phi =  Angle1 * pi/180.0
    theta = Angle2  * pi/180.0
    psi =  Angle3 * pi/180.0

    a11 = cos(psi) * cos(phi)  - cos(theta) * sin(phi) * sin(psi)   
    a12 = cos(psi) * sin(phi)  + cos(theta) * cos(phi) * sin(psi)   
    a13 = sin(psi) * sin(theta)     
    a21 = -sin(psi) * cos(phi)  - cos(theta) * sin(phi) *  cos(psi)     
    a22 = -sin(psi) * sin(phi)  + cos(theta) * cos(phi) * cos(psi)  
    a23 = cos(psi) * sin(theta)     
    a31 = sin(theta) * sin(phi)     
    a32 = -sin(theta) * cos(phi)    
    a33 = cos(theta)  

    a=scipy.matrix( [[a11,a12,a13],[a21,a22,a23],[ a31,a32,a33]])
    x=scipy.matrix( [[1.0],[0.0],[0.0]])
    y=scipy.matrix( [[0.0],[1.0],[0.0]])
    z=scipy.matrix( [[0.0],[0.0],[1.0]])

#    print a*x
#    print ''
#    print a*y
#    print ''
#    print a*z
    return a*x,a*y,a*z

def ReadTipsyFile(File):
    f = open(File,'r+')

    r = []
    BoverA = []
    CoverA = []
    phi = []
    theta = []
    psi = []
    
    for line in f:
      if '#' in line:
        continue
      if len(line)<5:
          continue
      
      tmp = line.split()
      r.append(float(tmp[0]))
      BoverA.append(float(tmp[1]))
      CoverA.append(float(tmp[2]))
      phi.append(float(tmp[3]))   
      theta.append(float(tmp[4]))   
      psi.append(float(tmp[5]))         
    return scipy.array(r),scipy.array(BoverA),scipy.array(CoverA),scipy.array(phi),scipy.array(theta),scipy.array(psi),



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

Case = 0
if Case == 0:
    FileName = '/home/ms/Uni/DarkMatter/AllSimulations/G/0G20_001'
    ShapeFile = '/home/ms/Uni/DarkMatter/AllSimulations/G/shapes/shape_0G20_001.txt'
    Offset = 0.6
    Color = ['green','orange','red']    
    Label = r'$G$-variation'    
if Case == 1:
    FileName = '/home/ms/Uni/DarkMatter/AllSimulations/HJS2010/0E21_001'
    ShapeFile = '/home/ms/Uni/DarkMatter/AllSimulations/HJS2010/shapes/shape_0E21_001.txt'
    Offset = 1.0
    Color = ['green','orange','red']
    Label = 'Energy exchange'    
if Case == 2:
    FileName = '/home/ms/Uni/DarkMatter/AllSimulations/G/s4G20_001'
    ShapeFile = '/home/ms/Uni/DarkMatter/AllSimulations/G/shapes/shape_s4G20_001.txt'
    Offset = 0.6
    Color = ['green','red','orange']
    Label = r'$G$-variation'
if Case == 3:
    FileName = '/home/ms/Uni/DarkMatter/AllSimulations/InhomogeneousInfall/Infall_sub_1e6.bin_201'
    ShapeFile = '/home/ms/Uni/DarkMatter/AllSimulations/InhomogeneousInfall/shapes/shape_Infall_sub_1e6.bin_201.txt'
    Offset = 0.9
    Color = ['green','orange','red']
    Label = 'Cold collapse'
if Case == 4:
    FileName = '/home/ms/Uni/DarkMatter/AllSimulations/CloneMerge/Mergers7_HQA_081'
    ShapeFile = '/home/ms/Uni/DarkMatter/AllSimulations/CloneMerge/shapes/shape_Mergers7_HQA_081.txt'
    Offset = 0.0
    Color = ['green','orange','red']
    Label = 'Major mergers'
if Case == 5:
    FileName = '/home/ms/Uni/DarkMatter/AllSimulations/ROI/OM_ROI00_rAN0.2_HQ_081'
    ShapeFile = '/home/ms/Uni/DarkMatter/AllSimulations/ROI/shapes/shape_OM_ROI00_rAN0.2_HQ_081.txt'
    Offset = 0.35
    Color = ['green','blue','red']
    Label = 'Unstable model'    
if Case == 6:
    FileName = '/home/ms/Uni/DarkMatter/AllSimulations/ROI/ROI05_120'
    ShapeFile = '/home/ms/Uni/DarkMatter/AllSimulations/ROI/shapes/shape_ROI05_120.txt'
    Offset = 0.2    
    Color = ['red','orange','green']
    Label = 'Unstable model'    
if Case == 7:
    FileName = '/home/ms/Uni/DarkMatter/AllSimulations/G/s1G20_001'
    ShapeFile = '/home/ms/Uni/DarkMatter/AllSimulations/G/shapes/shape_s1G20_001.txt'
    Offset = 0.6
    Color = ['green','orange','red']    
    Label = r'$G$-variation'    

   
   
    #FileName = '/home/ms/Uni/DarkMatter/AllSimulations/CloneMerge/Mergers7_HQA_081'
#ShapeFile = '/home/ms/Uni/DarkMatter/AllSimulations/CloneMerge/shapes/shape_Mergers7_HQA_081.txt'

plt.subplot(1,3,3)
plt.subplots_adjust(left=0.08, bottom=0.14, right=0.97, top=0.97,wspace=0.25, hspace=0.0)
r,BoverA,CoverA,phi,theta,psi = ReadTipsyFile(ShapeFile)
plt.plot(log10(r)+Offset,CoverA,'-o',color='red',label=r'$c/a$',lw=2,ms=7,mew=1)
plt.plot(log10(r)+Offset,BoverA,'-x',color='black',label=r'$b/a$',lw=2,ms=7,mew=1)
plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=4)

plt.ylabel(r'$c/a$ or $b/a$',fontsize=24)
plt.xlabel(r'$\log r/r_{-2}$',fontsize=24)
SetLabels(20,20)
plt.ylim((-0.0,1.05))



Major,Inter,Minor = GetPrincipalAxes(phi[10],theta[10],psi[10])

MajorX = Major[0,0]
MajorY = Major[1,0]
MajorZ = Major[2,0]

InterX = Inter[0,0]
InterY = Inter[1,0]
InterZ = Inter[2,0]

MinorX = Minor[0,0]
MinorY = Minor[1,0]
MinorZ = Minor[2,0]





A = DM_structure.DM_structure(FileName)
A.FindCenter()
A.FindCenterVel()
A.CenterParticlePositions()
A.CenterParticleVelocities()
GridSph = A.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50,CalcUncBeta=True)



Ncones = 0
AngleFile = open('48.txt','r')
#AngleFile = open('../192.txt','r')
for line in AngleFile:
    tmp = line.split()
    if len(tmp) != 5:
        continue
    
    x = float(tmp[0])
    y = float(tmp[1])
    z = float(tmp[2])
    
    MajorAngle = scipy.arccos((x*MajorX + y*MajorY + z*MajorZ) / scipy.sqrt(x**2+y**2+z**2) / scipy.sqrt(MajorX**2+MajorY**2+MajorZ**2))
    InterAngle = scipy.arccos((x*InterX + y*InterY + z*InterZ) / scipy.sqrt(x**2+y**2+z**2) / scipy.sqrt(InterX**2+InterY**2+InterZ**2))    
    MinorAngle = scipy.arccos((x*MinorX + y*MinorY + z*MinorZ) / scipy.sqrt(x**2+y**2+z**2) / scipy.sqrt(MinorX**2+MinorY**2+MinorZ**2))
    
    B = copy.deepcopy(A)
    B.Snapshot.SelectParticlesInCone(x,y,z,3.1415/8.0)

    GridSphB = B.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50,CalcUncBeta=True)
    if MajorAngle < 0.0:
        print '--John john'
    
    plt.subplot(1,3,1)
    if MajorAngle<3.1415/8.0:
        plt.plot(log10(GridSphB.R)+Offset,GridSphB.Beta,'-',color=Color[0],lw=1.5,zorder=6)
        print '---------------------------'
        print '--------Major--------------'
        print '---------------------------'
        print MajorX,MajorY,MajorZ
        print x,y,z
    elif InterAngle<3.1415/8.0:
        plt.plot(log10(GridSphB.R)+Offset,GridSphB.Beta,'-',color=Color[1],lw=1.5,zorder=6)
    elif MinorAngle<3.1415/8.0:
        plt.plot(log10(GridSphB.R)+Offset,GridSphB.Beta,'-',color=Color[2],lw=1.5,zorder=6)
        print '---------------------------'
        print '--------Minor--------------'
        print '---------------------------'
        print MinorX,MinorY,MinorZ
        print x,y,z        
    else:
        plt.plot(log10(GridSphB.R)+Offset,GridSphB.Beta,'-',color='black',lw=0.5,zorder=5)

    plt.subplot(1,3,2)
    if MajorAngle<3.1415/8.0:
        plt.plot(log10(GridSphB.R)+Offset,log10(GridSphB.R**0*GridSphB.Rho),'-',color=Color[0],lw=1.5,zorder=6)
    elif InterAngle<3.1415/8.0:
        plt.plot(log10(GridSphB.R)+Offset,log10(GridSphB.R**0*GridSphB.Rho),'-',color=Color[1],lw=1.5,zorder=6)
    elif MinorAngle<3.1415/8.0:
        plt.plot(log10(GridSphB.R)+Offset,log10(GridSphB.R**0*GridSphB.Rho),'-',color=Color[2],lw=1.5,zorder=6)        
    else:
        plt.plot(log10(GridSphB.R)+Offset,log10(GridSphB.R**0*GridSphB.Rho),'-',color='black',lw=0.5,zorder=5)        
        
    Ncones += 1
    print Ncones

plt.subplot(1,3,1)
plt.xlim((-1,2))
plt.ylim((-1,1))
plt.ylabel(r'$\beta$',fontsize=24)
plt.xlabel(r'$\log r/r_{-2}$',fontsize=24)
SetLabels(20,20)

if Case == 2:
    t=scipy.array([0.14,1.14])
    plt.plot(t,0*t+0.5,'-',color='cyan', lw = 4,zorder=7)
elif Case == 1:
    t=scipy.array([0.31,0.76])
    plt.plot(t,0*t+0.5,'-',color='cyan', lw = 4,zorder=7)

plt.subplot(1,3,2)
plt.xlim((-1,2))

plt.ylabel(r'$\log \rho$',fontsize=24)
plt.xlabel(r'$\log r/r_{-2}$',fontsize=24)
SetLabels(20,20)
if Case in [1,3,4]:
    plt.ylim((-10,2))
    plt.text(-0.25,0.6,Label,color='black',fontsize=22)
else:
    plt.ylim((-10,0))
    plt.text(-0.25,0.6-2,Label,color='black',fontsize=22)

plt.text(-0.75,-6,'Major',color='red',fontsize=22)
plt.text(-0.75,-7,'Intermediate',color='green',fontsize=22)
plt.text(-0.75,-8,'Minor',color='orange',fontsize=22)

if Case == 2:
    t=scipy.array([0.64,0.85])
    plt.plot(t,0*t-4.165,'-',color='cyan', lw = 4,zorder=7)
if Case == 1:
    t=scipy.array([0.33,0.64])
    plt.plot(t,0*t-2.02,'-',color='cyan', lw = 4,zorder=7)

plt.show()