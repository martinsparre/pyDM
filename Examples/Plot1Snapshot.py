#!/usr/bin/python
import array
import scipy
import os.path
from mayavi import mlab
#from Mayavi2 import DrawSphere
from MartinMisc import Undersample,IndexOfValuesInRange
import Classes.Gadget2 as Gadget2
import Classes.DM_structure as DM_structure




A = DM_structure.DM_structure('/home/ms/Uni/DarkMatter/AllSimulations/HJS2010/s4_16')
A.FindCenter()
A.FindCenterVel()
A.CenterParticlePositions()
A.CenterParticleVelocities()
GridSph = A.CreateGridLogBins(NBins=100,Rmin=0.001)


G0 = A.Snapshot

x = G0.x
y = G0.y
z = G0.z
V = G0.V
P = Undersample([x,y,z,V],N=60000)
#P = [x,y,z,V]

mlab.figure(figure='Snapshot',bgcolor=(0,0,0),size=(600,600))
mlab.clf()
pts = mlab.points3d(P[0],P[1],P[2],P[3],scale_mode='none',scale_factor=1.0,resolution=3,opacity=0.9)


linvalues = scipy.linspace(0.0,10.0,100)
logvalues = scipy.logspace(-1,1.0,100)
MajorX = -0.26871747*logvalues
MajorY = -0.70001169*logvalues
MajorZ = -0.66164534*logvalues
MajorV = linvalues*3.0

logvalues = scipy.logspace(-1,1,100)
InterX = 0.57957833*logvalues
InterY = 0.43113878*logvalues
InterZ = -0.69152608*logvalues
InterV = linvalues
logvalues = scipy.logspace(-1,1.0,100)
MinorX = 0.7693373*logvalues
MinorY = -0.56930044*logvalues
MinorZ = 0.28985708*logvalues
MinorV = linvalues/3.0


mlab.points3d(MajorX,MajorY,MajorZ,MajorV,scale_mode='none',scale_factor=1.0,resolution=3,opacity=0.7,vmin=0.0,vmax=20.0)#extent=[-10,10,-10,10,-10,10]
mlab.points3d(InterX,InterY,InterZ,InterV,scale_mode='none',scale_factor=1.0,resolution=3,opacity=0.7,vmin=0.0,vmax=20.0)#extent=[-10,10,-10,10,-10,10]
mlab.points3d(MinorX,MinorY,MinorZ,MinorV,scale_mode='none',scale_factor=1.0,resolution=3,opacity=0.7,vmin=0.0,vmax=20.0)#extent=[-10,10,-10,10,-10,10]
#DrawSphere(0.6,-2.90751338,  1.84735751, -3.1246891,Color=(0.8,0.5,0.5))

#pts.scene.disable_render = False
#    mlab.savefig('a'+str(N)+'.png')
#    mlab.draw()
mlab.show()
        


