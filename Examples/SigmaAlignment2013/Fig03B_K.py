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
        

Filenames = ['0.1/1HqIso_Impact0_160','0.3/1HqIso_Impact0_160','0.5/1HqIso_Impact0_160']
NFiles = len(Filenames)
DIR = '/home/ms/Uni/DarkMatter/AllSimulations/SigmaAlignment2013/HeadonMerger_VaryingVel/'

for i in range(NFiles):

    FileName=DIR + Filenames[i]
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
    C.Snapshot.SelectParticlesInCone(0,1,0,3.1415/8.0)
    D.Snapshot.SelectParticlesInCone(0,0,1,3.1415/8.0)
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

    plt.subplots_adjust(left=0.07, bottom=0.12, right=0.97, top=0.93,wspace=0.0, hspace=0.0)
    plt.subplot(NFiles,3,3*i+1)
    #IC = DM_structure.DM_structure('../1HqIso_000')
    #ICGridSph = IC.CreateGridLogBins(NBins=25,Rmin=0.001,Rmax=9.9)
    #plt.plot(log10(ICGridSph.R),ICGridSph.Beta,'--',label=r'IC',color='black',lw=2,ms=9,mew=2)
    plt.title(r'$x$-axis cone', fontsize=24)

    plt.plot(log10(GridSphB.R),GridSphB.Ekinx,'-o',label=r'$K_x$',color='blue',lw=2,ms=9,mew=2)
    plt.plot(log10(GridSphB.R),GridSphB.Ekiny,'-D',label=r'$K_y$',color='red',lw=2,ms=9,mew=2)
    plt.plot(log10(GridSphB.R),GridSphB.Ekinz,'-<',label=r'$K_z$',color='green',lw=2,ms=9,mew=2)


    #plt.ylim((-0.45,0.75))
    plt.ylim((-0.0,0.13))
    plt.xlim((-1,0.99))
    plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=1)
    plt.grid()
    plt.xlabel(r'$\log r$',fontsize=24)
    plt.ylabel(r'Kinetic energy, $K$',fontsize=24)
    SetLabels(20,20)


    plt.subplot(NFiles,3,3*i+2)
    #IC = DM_structure.DM_structure('../1HqIso_000')
    #ICGridSph = IC.CreateGridLogBins(NBins=25,Rmin=0.001,Rmax=9.9)
    #plt.plot(log10(ICGridSph.R),ICGridSph.Beta,'--',label=r'IC',color='black',lw=2,ms=9,mew=2)

    plt.title(r'$y$-axis cone', fontsize=24)
    plt.plot(log10(GridSphC.R),GridSphC.Ekinx,'-o',label=r'$K_x$',color='blue',lw=2,ms=9,mew=2)
    plt.plot(log10(GridSphC.R),GridSphC.Ekiny,'-D',label=r'$K_y$',color='red',lw=2,ms=9,mew=2)
    plt.plot(log10(GridSphC.R),GridSphC.Ekinz,'-<',label=r'$K_z$',color='green',lw=2,ms=9,mew=2)


    plt.ylim((-0.0,0.13))
    plt.xlim((-1,0.99))
    plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=1)
    plt.grid()
    plt.xlabel(r'$\log r$',fontsize=24)
    SetLabels(20,1)

    plt.subplot(NFiles,3,3*i+3)
    #IC = DM_structure.DM_structure('../1HqIso_000')
    #ICGridSph = IC.CreateGridLogBins(NBins=25,Rmin=0.001,Rmax=9.9)
    #plt.plot(log10(ICGridSph.R),ICGridSph.Beta,'--',label=r'IC',color='black',lw=2,ms=9,mew=2)

    plt.title(r'$z$-axis cone', fontsize=24)

    plt.plot(log10(GridSphD.R),GridSphD.Ekinx,'-o',label=r'$K_x$',color='blue',lw=2,ms=9,mew=2)
    plt.plot(log10(GridSphD.R),GridSphD.Ekiny,'-D',label=r'$K_y$',color='red',lw=2,ms=9,mew=2)
    plt.plot(log10(GridSphD.R),GridSphD.Ekinz,'-<',label=r'$K_z$',color='green',lw=2,ms=9,mew=2)



    #plt.ylim((-0.45,0.75))
    plt.ylim((-0.0,0.13))
    plt.xlim((-1,1))
    plt.legend(prop=dict(size=18), numpoints=2, ncol=1,frameon=True,loc=1)
    plt.grid()
    plt.xlabel(r'$\log r$',fontsize=24)
    SetLabels(20,1)


plt.show()


