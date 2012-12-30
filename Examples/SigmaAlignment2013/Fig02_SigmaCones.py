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
import healpy as hp



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
        

#Filenames = ['0.1/1HqIso_Impact0_160','0.3/1HqIso_Impact0_160','0.5/1HqIso_Impact0_160','1.0/1HqIso_Impact0_160','1.5/1HqIso_Impact0_160']
Filenames = ['0.5/1HqIso_Impact0_160']
NFiles = len(Filenames)
DIR = '/home/ms/Uni/DarkMatter/AllSimulations/SigmaAlignment2013/HeadonMerger_VaryingVel/'


FileName='/home/ms/Uni/DarkMatter/AllSimulations/SigmaAlignment2013/HeadonMerger_VaryingVel/0.5/1HqIso_Impact0_160'
A = DM_structure.DM_structure(FileName)
A.FindCenter()
A.FindCenterVel()
A.CenterParticlePositions()
A.CenterParticleVelocities()
GridSph = A.CreateGridLogBins(NBins=100,Rmin=0.001)


B = copy.deepcopy(A)
C = copy.deepcopy(A)
D = copy.deepcopy(A)
#E = copy.deepcopy(A)

B.Snapshot.SelectParticlesInCone(1,0,0,3.1415/8.0)
C.Snapshot.SelectParticlesInCone(1.0/scipy.sqrt(2.0),1.0/scipy.sqrt(2.0),0,3.1415/8.0)
D.Snapshot.SelectParticlesInCone(0,1,0,3.1415/8.0)
#E.Snapshot.SelectParticlesInCone(-1,0,0,3.1415/8.0)


PlotParticles = False
if PlotParticles == True:
    plt.title('File'+sys.argv[1])
    plt.plot(B.Snapshot.x,B.Snapshot.y,'.',label='$(1,0,0)$ $\pi/4$')
    plt.plot(C.Snapshot.x,C.Snapshot.y,'.',label='$(0,1,0)$ $\pi/16$')
#    plt.plot(D.Snapshot.x,D.Snapshot.y,'.',label='$(0,0,1)$ $\pi/16$')
#    plt.plot(E.Snapshot.x,E.Snapshot.y,'.',label='$(-1,0,0)$ $\pi/16$')
    plt.legend(loc=2)
    plt.grid()
    plt.xlabel('x')
    plt.xlabel('y')
    SetLabels()
    plt.show()


GridSphA = A.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
GridSphB = B.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
GridSphC = C.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
GridSphD = D.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
#GridSphE = E.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)


#sys.exit()

plt.subplots_adjust(left=0.1, bottom=0.04, right=0.97, top=0.96,wspace=0.0, hspace=0.15)
plt.subplot(3,2,1)
#IC = DM_structure.DM_structure('../1HqIso_000')
#ICGridSph = IC.CreateGridLogBins(NBins=25,Rmin=0.001,Rmax=9.9)
#plt.plot(log10(ICGridSph.R),ICGridSph.Beta,'--',label=r'IC',color='black',lw=2,ms=9,mew=2)
#plt.title(r'$x$-axis cone', fontsize=24)

plt.plot(log10(GridSphB.R),2*GridSphB.Ekinx,'-o',label=r'$\sigma^2_x$',color='blue',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphB.R),2*GridSphB.Ekiny,'-D',label=r'$\sigma^2_y$',color='red',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphB.R),2*GridSphB.Ekinz,'-<',label=r'$\sigma^2_z$',color='green',lw=2,ms=9,mew=2)


#plt.ylim((-0.45,0.75))
plt.ylim((-0.0,0.26))
plt.xlim((-1,1.19))
plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=False,loc=1)
plt.grid()
#plt.xlabel(r'$\log r$',fontsize=24)
plt.ylabel(r'$\sigma^2$',fontsize=24)
SetLabels(1,20)



plt.subplot(3,2,3)
#IC = DM_structure.DM_structure('../1HqIso_000')
#ICGridSph = IC.CreateGridLogBins(NBins=25,Rmin=0.001,Rmax=9.9)
#plt.plot(log10(ICGridSph.R),ICGridSph.Beta,'--',label=r'IC',color='black',lw=2,ms=9,mew=2)

#plt.title(r'$x+y$-axis cone', fontsize=24)
plt.plot(log10(GridSphC.R),2*GridSphC.Ekinx,'-o',label=r'$\sigma^2_x$',color='blue',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphC.R),2*GridSphC.Ekiny,'-D',label=r'$\sigma^2_y$',color='red',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphC.R),2*GridSphC.Ekinz,'-<',label=r'$\sigma^2_z$',color='green',lw=2,ms=9,mew=2)


plt.ylim((-0.0,0.26))
plt.xlim((-1,1.19))
plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=False,loc=1)
plt.grid()
#plt.xlabel(r'$\log r$',fontsize=24)
plt.ylabel(r'$\sigma^2$',fontsize=24)
SetLabels(1,20)



plt.subplot(3,2,5)
#IC = DM_structure.DM_structure('../1HqIso_000')
#ICGridSph = IC.CreateGridLogBins(NBins=25,Rmin=0.001,Rmax=9.9)
#plt.plot(log10(ICGridSph.R),ICGridSph.Beta,'--',label=r'IC',color='black',lw=2,ms=9,mew=2)

#plt.title(r'$y$-axis cone', fontsize=24)
plt.plot(log10(GridSphC.R),2*GridSphC.Ekinx,'-o',label=r'$\sigma^2_x$',color='blue',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphC.R),2*GridSphC.Ekiny,'-D',label=r'$\sigma^2_y$',color='red',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphC.R),2*GridSphC.Ekinz,'-<',label=r'$\sigma^2_z$',color='green',lw=2,ms=9,mew=2)


plt.ylim((-0.0,0.26))
plt.xlim((-1,1.19))
plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=False,loc=1)
plt.grid()
plt.xlabel(r'$\log r$',fontsize=24)
plt.ylabel(r'$\sigma^2$',fontsize=24)
SetLabels(20,20)


#plt.subplot(3,2,2)
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

NSIDE = 4
m = scipy.array(SigmaCones)

hp.mollview(m, title=r"The $\sigma^2$-ellipsoid at $\mathcal{A}$",sub=[3,2,2],min=0.01,max=0.07,flip='geo',margins = (0.02,0.02,0.02,0.02))
plt.text(0,0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(-1.8,0.0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(1.8,0.0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)



#plt.subplot(3,2,4)
R = scipy.sqrt(C.Snapshot.x**2+C.Snapshot.y**2+C.Snapshot.z**2)
C.Snapshot.SelectParticles((R>5) * (R<6))

Tensor = scipy.zeros((3,3))
Tensor[0,0] = scipy.mean(C.Snapshot.vx * C.Snapshot.vx) - scipy.mean(C.Snapshot.vx)  *scipy.mean( C.Snapshot.vx)
Tensor[1,1] = scipy.mean(C.Snapshot.vy * C.Snapshot.vy) - scipy.mean(C.Snapshot.vy)  *scipy.mean( C.Snapshot.vy)               
Tensor[2,2] = scipy.mean(C.Snapshot.vz * C.Snapshot.vz) - scipy.mean(C.Snapshot.vz)  *scipy.mean( C.Snapshot.vz)
Tensor[0,1] = scipy.mean(C.Snapshot.vx * C.Snapshot.vy) - scipy.mean(C.Snapshot.vx)  *scipy.mean( C.Snapshot.vy)
Tensor[1,0] = Tensor[0,1] 
Tensor[0,2] = scipy.mean(C.Snapshot.vx * C.Snapshot.vz) - scipy.mean(C.Snapshot.vx)  *scipy.mean( C.Snapshot.vz)
Tensor[2,0] = Tensor[0,2]  
Tensor[1,2] = scipy.mean(C.Snapshot.vy * C.Snapshot.vz) - scipy.mean(C.Snapshot.vy)  *scipy.mean( C.Snapshot.vz)
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


NSIDE = 4
m = scipy.array(SigmaCones)

hp.mollview(m, title=r"The $\sigma^2$-ellipsoid at  $\mathcal{B}$",sub=[3,2,4],min=0.01,max=0.07,flip='geo', margins = (0.02,0.02,0.02,0.02))
plt.text(0,0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(-1.8,0.0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(1.8,0.0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
















#plt.subplot(3,2,6)
R = scipy.sqrt(D.Snapshot.x**2+D.Snapshot.y**2+D.Snapshot.z**2)
D.Snapshot.SelectParticles((R>5) * (R<6))

Tensor = scipy.zeros((3,3))
Tensor[0,0] = scipy.mean(D.Snapshot.vx * D.Snapshot.vx) - scipy.mean(D.Snapshot.vx)  *scipy.mean( D.Snapshot.vx)
Tensor[1,1] = scipy.mean(D.Snapshot.vy * D.Snapshot.vy) - scipy.mean(D.Snapshot.vy)  *scipy.mean( D.Snapshot.vy)               
Tensor[2,2] = scipy.mean(D.Snapshot.vz * D.Snapshot.vz) - scipy.mean(D.Snapshot.vz)  *scipy.mean( D.Snapshot.vz)
Tensor[0,1] = scipy.mean(D.Snapshot.vx * D.Snapshot.vy) - scipy.mean(D.Snapshot.vx)  *scipy.mean( D.Snapshot.vy)
Tensor[1,0] = Tensor[0,1] 
Tensor[0,2] = scipy.mean(D.Snapshot.vx * D.Snapshot.vz) - scipy.mean(D.Snapshot.vx)  *scipy.mean( D.Snapshot.vz)
Tensor[2,0] = Tensor[0,2]  
Tensor[1,2] = scipy.mean(D.Snapshot.vy * D.Snapshot.vz) - scipy.mean(D.Snapshot.vy)  *scipy.mean( D.Snapshot.vz)
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


NSIDE = 4
m = scipy.array(SigmaCones)

hp.mollview(m, title=r"The $\sigma^2$-ellipsoid at $\mathcal{C}$",sub=[3,2,6],min=0.01,max=0.07,flip='geo',margins = (0.02,0.02,0.02,0.02))
plt.text(0,0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(-1.8,0.0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(1.8,0.0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
















plt.show()


