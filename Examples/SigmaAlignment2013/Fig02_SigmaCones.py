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
import matplotlib

Rmin = 5
Rmax = 6
SigmaMin = 0.01
SigmaMax = 0.07
SigmaMax = 0.06

ImpactParameter = True 
plt.rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Roman']})
plt.rc('font',**{'family':'serif','serif':['Computer Modern Sans serif']})
plt.rc('text', usetex=True)

plt.rcParams["xtick.major.size"] = 5
plt.rcParams["ytick.major.size"] = 5

#Functions:
def SetLabels(xsize=18,ysize=18):
    ax = pylab.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(xsize)

    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(ysize)
        

#FileName='/home/ms/Uni/DarkMatter/AllSimulations/SigmaAlignment2013/HeadonMerger_VaryingVel/0.5/1HqIso_Impact0_160'
FileName='/home/ms/Uni/DarkMatter/AllSimulations/MergerAnisotropy/1HqIso_Impact10_120'


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


B1 = copy.deepcopy(B)
C1 = copy.deepcopy(C)
D1 = copy.deepcopy(D)


GridSphA = A.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
GridSphB = B.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
GridSphC = C.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)
GridSphD = D.CreateGridLogBins(NBins=50,Rmin=0.01,Rmax=50)


plt.subplots_adjust(left=0.08, bottom=0.06, right=0.91, top=0.98,wspace=0.0, hspace=0.16)
plt.subplot(3,3,1)
#IC = DM_structure.DM_structure('../1HqIso_000')
#ICGridSph = IC.CreateGridLogBins(NBins=25,Rmin=0.001,Rmax=9.9)
#plt.plot(log10(ICGridSph.R),ICGridSph.Beta,'--',label=r'IC',color='black',lw=2,ms=9,mew=2)
#plt.title(r'$x$-axis cone', fontsize=24)

plt.plot(log10(GridSphB.R),2*GridSphB.Ekinx,'-o',label=r'$\sigma^2_x$',color='blue',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphB.R),2*GridSphB.Ekiny,'-D',label=r'$\sigma^2_y$',color='red',lw=2,ms=9,mew=2)
plt.plot(log10(GridSphB.R),2*GridSphB.Ekinz,'-<',label=r'$\sigma^2_z$',color='green',lw=2,ms=9,mew=2)


plt.ylim((-0.0,0.26))
plt.xlim((-1,1.19))
plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=False,loc=1)
plt.grid()
#plt.xlabel(r'$\log r$',fontsize=24)
plt.ylabel(r'$\sigma^2$',fontsize=20)
SetLabels(1,18)



plt.subplot(3,3,4)
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
plt.ylabel(r'$\sigma^2$',fontsize=20)
SetLabels(1,18)



plt.subplot(3,3,7)
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
plt.xlabel(r'$\log r$',fontsize=20)
plt.ylabel(r'$\sigma^2$',fontsize=20)
SetLabels(18,18)


#plt.subplot(3,2,2)
R = scipy.sqrt(B.Snapshot.x**2+B.Snapshot.y**2+B.Snapshot.z**2)
B.Snapshot.SelectParticles((R>Rmin) * (R<Rmax))


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

matplotlib.rcParams.update({'font.size': 18})

hp.mollview(m, title="",sub=[3,3,2],min=SigmaMin,max=SigmaMax,flip='geo',margins = (0.02,0.02,0.02,0.02))
plt.text(0,0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(-1.8,0.0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(1.8,0.0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(0,1.2,r'The $\sigma^2$-ellipsoid at $\mathcal{A}$',horizontalalignment='center',verticalalignment='center',fontsize=16)


#plt.subplot(3,2,4)
R = scipy.sqrt(C.Snapshot.x**2+C.Snapshot.y**2+C.Snapshot.z**2)
C.Snapshot.SelectParticles((R>Rmin) * (R<Rmax))

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

hp.mollview(m, title=r"",sub=[3,3,5],min=SigmaMin,max=SigmaMax,flip='geo', margins = (0.02,0.02,0.02,0.02))
plt.text(0,0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(-1.8,0.0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(1.8,0.0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(0,1.2,r'The $\sigma^2$-ellipsoid at  $\mathcal{B}$',horizontalalignment='center',verticalalignment='center',fontsize=16)




#plt.subplot(3,2,6)
R = scipy.sqrt(D.Snapshot.x**2+D.Snapshot.y**2+D.Snapshot.z**2)
D.Snapshot.SelectParticles((R>Rmin) * (R<Rmax))

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

hp.mollview(m, title="",sub=[3,3,8],min=SigmaMin,max=SigmaMax,flip='geo',margins = (0.02,0.02,0.02,0.02))
plt.text(0,0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(-1.8,0.0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(1.8,0.0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(0,1.2,r'The $\sigma^2$-ellipsoid at  $\mathcal{C}$',horizontalalignment='center',verticalalignment='center',fontsize=16)


Rmin = 0.2
Rmax = 0.6
SigmaMin = 0.1
SigmaMax = 0.22
SigmaMax = 0.18

R = scipy.sqrt(B1.Snapshot.x**2+B1.Snapshot.y**2+B1.Snapshot.z**2)
B1.Snapshot.SelectParticles((R>Rmin) * (R<Rmax))


Tensor = scipy.zeros((3,3))
Tensor[0,0] = scipy.mean(B1.Snapshot.vx * B1.Snapshot.vx) - scipy.mean(B1.Snapshot.vx)  *scipy.mean( B1.Snapshot.vx)
Tensor[1,1] = scipy.mean(B1.Snapshot.vy * B1.Snapshot.vy) - scipy.mean(B1.Snapshot.vy)  *scipy.mean( B1.Snapshot.vy)               
Tensor[2,2] = scipy.mean(B1.Snapshot.vz * B1.Snapshot.vz) - scipy.mean(B1.Snapshot.vz)  *scipy.mean( B1.Snapshot.vz)
Tensor[0,1] = scipy.mean(B1.Snapshot.vx * B1.Snapshot.vy) - scipy.mean(B1.Snapshot.vx)  *scipy.mean( B1.Snapshot.vy)
Tensor[1,0] = Tensor[0,1] 
Tensor[0,2] = scipy.mean(B1.Snapshot.vx * B1.Snapshot.vz) - scipy.mean(B1.Snapshot.vx)  *scipy.mean( B1.Snapshot.vz)
Tensor[2,0] = Tensor[0,2]  
Tensor[1,2] = scipy.mean(B1.Snapshot.vy * B1.Snapshot.vz) - scipy.mean(B1.Snapshot.vy)  *scipy.mean( B1.Snapshot.vz)
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

matplotlib.rcParams.update({'font.size': 18})

hp.mollview(m, title="",sub=[3,3,3],min=SigmaMin,max=SigmaMax,flip='geo',margins = (0.02,0.02,0.02,0.02))
plt.text(0,0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(-1.8,0.0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(1.8,0.0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(0,1.2,r'The $\sigma^2$-ellipsoid at $\mathcal{A}$',horizontalalignment='center',verticalalignment='center',fontsize=16)


#plt.subplot(3,2,4)
R = scipy.sqrt(C1.Snapshot.x**2+C1.Snapshot.y**2+C1.Snapshot.z**2)
C1.Snapshot.SelectParticles((R>Rmin) * (R<Rmax))

Tensor = scipy.zeros((3,3))
Tensor[0,0] = scipy.mean(C1.Snapshot.vx * C1.Snapshot.vx) - scipy.mean(C1.Snapshot.vx)  *scipy.mean( C1.Snapshot.vx)
Tensor[1,1] = scipy.mean(C1.Snapshot.vy * C1.Snapshot.vy) - scipy.mean(C1.Snapshot.vy)  *scipy.mean( C1.Snapshot.vy)               
Tensor[2,2] = scipy.mean(C1.Snapshot.vz * C1.Snapshot.vz) - scipy.mean(C1.Snapshot.vz)  *scipy.mean( C1.Snapshot.vz)
Tensor[0,1] = scipy.mean(C1.Snapshot.vx * C1.Snapshot.vy) - scipy.mean(C1.Snapshot.vx)  *scipy.mean( C1.Snapshot.vy)
Tensor[1,0] = Tensor[0,1] 
Tensor[0,2] = scipy.mean(C1.Snapshot.vx * C1.Snapshot.vz) - scipy.mean(C1.Snapshot.vx)  *scipy.mean( C1.Snapshot.vz)
Tensor[2,0] = Tensor[0,2]  
Tensor[1,2] = scipy.mean(C1.Snapshot.vy * C1.Snapshot.vz) - scipy.mean(C1.Snapshot.vy)  *scipy.mean( C1.Snapshot.vz)
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

hp.mollview(m, title=r"",sub=[3,3,6],min=SigmaMin,max=SigmaMax,flip='geo', margins = (0.02,0.02,0.02,0.02))
plt.text(0,0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(-1.8,0.0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(1.8,0.0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(0,1.2,r'The $\sigma^2$-ellipsoid at  $\mathcal{B}$',horizontalalignment='center',verticalalignment='center',fontsize=16)




#plt.subplot(3,2,6)
R = scipy.sqrt(D1.Snapshot.x**2+D1.Snapshot.y**2+D1.Snapshot.z**2)
D1.Snapshot.SelectParticles((R>Rmin) * (R<Rmax))

Tensor = scipy.zeros((3,3))
Tensor[0,0] = scipy.mean(D1.Snapshot.vx * D1.Snapshot.vx) - scipy.mean(D1.Snapshot.vx)  *scipy.mean( D1.Snapshot.vx)
Tensor[1,1] = scipy.mean(D1.Snapshot.vy * D1.Snapshot.vy) - scipy.mean(D1.Snapshot.vy)  *scipy.mean( D1.Snapshot.vy)               
Tensor[2,2] = scipy.mean(D1.Snapshot.vz * D1.Snapshot.vz) - scipy.mean(D1.Snapshot.vz)  *scipy.mean( D1.Snapshot.vz)
Tensor[0,1] = scipy.mean(D1.Snapshot.vx * D1.Snapshot.vy) - scipy.mean(D1.Snapshot.vx)  *scipy.mean( D1.Snapshot.vy)
Tensor[1,0] = Tensor[0,1] 
Tensor[0,2] = scipy.mean(D1.Snapshot.vx * D1.Snapshot.vz) - scipy.mean(D1.Snapshot.vx)  *scipy.mean( D1.Snapshot.vz)
Tensor[2,0] = Tensor[0,2]  
Tensor[1,2] = scipy.mean(D1.Snapshot.vy * D1.Snapshot.vz) - scipy.mean(D1.Snapshot.vy)  *scipy.mean( D1.Snapshot.vz)
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

hp.mollview(m, title="",sub=[3,3,9],min=SigmaMin,max=SigmaMax,flip='geo',margins = (0.02,0.02,0.02,0.02))
plt.text(0,0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(-1.8,0.0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(1.8,0.0,'Collision\naxis',horizontalalignment='center',verticalalignment='center',fontsize=16)
plt.text(0,1.2,r'The $\sigma^2$-ellipsoid at  $\mathcal{C}$',horizontalalignment='center',verticalalignment='center',fontsize=16)



plt.show()


