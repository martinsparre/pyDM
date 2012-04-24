#!/usr/bin/python
import scipy
import Gadget2
import copy
from math import sqrt


def CreateClone(OldGadgetFilename, NewFilename, DeltaPos = [0,0,0], DeltaVel = [0,0,0]):
    "This function clones a structure. It can be used to create IC's for merger collisions."
    G = Gadget2.ReadGadget2(OldGadgetFilename)
    #center structure:
    G.x -= G.x.mean()
    G.y -= G.y.mean()
    G.z -= G.z.mean()

    G.vx -= G.vx.mean()
    G.vy -= G.vy.mean()
    G.vz -= G.vz.mean()

    print 'Mean Pos', G.x.mean(),G.y.mean(),G.z.mean()
    print 'Mean Pos',G.vx.mean(),G.vy.mean(),G.vz.mean()

    #create a clone:
    G0 = copy.deepcopy(G)

    deltax=DeltaPos[0]
    deltay=DeltaPos[1]
    deltaz=DeltaPos[2]

    deltavx=DeltaVel[0]
    deltavy=DeltaVel[1]
    deltavz=DeltaVel[2]

    G.x +=deltax
    G0.x -=deltax
    G.y +=deltay
    G0.y -=deltay
    G.z +=deltaz
    G0.z -=deltaz

    G.vx +=deltavx
    G0.vx -=deltavx
    G.vy +=deltavy
    G0.vy -=deltavy
    G.vz +=deltavz
    G0.vz -=deltavz

    x = scipy.array([G.x,G0.x]).flatten()
    y = scipy.array([G.y,G0.y]).flatten()
    z = scipy.array([G.z,G0.z]).flatten()

    vx = scipy.array([G.vx,G0.vx]).flatten()
    vy = scipy.array([G.vy,G0.vy]).flatten()
    vz = scipy.array([G.vz,G0.vz]).flatten()

    m = scipy.array([G.m,G0.m]).flatten()

    #ID = scipy.array([G.ID,G0.ID]).flatten()
    ID = scipy.array(range(len(m)))
    
    NPart = scipy.array(scipy.array([0,len(m),0,0,0,0]))
    NPartTotal = scipy.array([len(m)])
    MassArray = scipy.array(6*[0])
    Time = scipy.array([0.0])

    New = Gadget2.Gadget2Snapshot(x=x, y=y, z=z,vx=vx,vy=vy,vz=vz,ID=ID, m=m,NPartTotal=NPartTotal, NPart=NPart,Time=Time,MassArray=MassArray, V=None)

    Gadget2.WriteGadget2(NewFilename , New)

def MergeHaloes(OldGadgetFilename1,OldGadgetFilename2, NewFilename, DeltaPos = [0,0,0], DeltaVel = [0,0,0],ImportPotential=True):
    "This function clones a structure. It can be used to create IC's for merger collisions."
    G = Gadget2.ReadGadget2(OldGadgetFilename1,ImportPotential=ImportPotential)
    G0 = Gadget2.ReadGadget2(OldGadgetFilename2,ImportPotential=ImportPotential)


    deltax=DeltaPos[0]
    deltay=DeltaPos[1]
    deltaz=DeltaPos[2]

    deltavx=DeltaVel[0]
    deltavy=DeltaVel[1]
    deltavz=DeltaVel[2]


    G0.x += deltax
    G0.y += deltay
    G0.z += deltaz

    G0.vx += deltavx
    G0.vy += deltavy
    G0.vz += deltavz



    x = scipy.append(G.x,G0.x)
    y = scipy.append(G.y,G0.y)
    z = scipy.append(G.z,G0.z)

    vx = scipy.append(G.vx,G0.vx)
    vy = scipy.append(G.vy,G0.vy)
    vz = scipy.append(G.vz,G0.vz)

    m = scipy.append(G.m,G0.m)
    print len(x),len(y),len(z),len(vx),len(vy),len(vz),len(m)

    ID = scipy.array(range(len(m)))
    
    NPart = scipy.array(scipy.array([0,len(m),0,0,0,0]))
    NPartTotal = scipy.array([len(m)])
    MassArray = scipy.array(6*[0])
    Time = scipy.array([0.0])

    New = Gadget2.Gadget2Snapshot(x=x, y=y, z=z,vx=vx,vy=vy,vz=vz,ID=ID, m=m,NPartTotal=NPartTotal, NPart=NPart,Time=Time,MassArray=MassArray, V=None)

    Gadget2.WriteGadget2(NewFilename , New)

    
    
def MergeSnapshots(G,G0,ImportPotential=True):
    "This function clones a structure. It can be used to create IC's for merger collisions."
#    G = Gadget2.ReadGadget2(OldGadgetFilename1,ImportPotential=ImportPotential)
#    G0 = Gadget2.ReadGadget2(OldGadgetFilename2,ImportPotential=ImportPotential)


    x = scipy.append(G.x,G0.x)
    y = scipy.append(G.y,G0.y)
    z = scipy.append(G.z,G0.z)

    vx = scipy.append(G.vx,G0.vx)
    vy = scipy.append(G.vy,G0.vy)
    vz = scipy.append(G.vz,G0.vz)

    m = scipy.append(G.m,G0.m)
    print len(x),len(y),len(z),len(vx),len(vy),len(vz),len(m)

    ID = scipy.array(range(len(m)))
    
    NPart = scipy.array(scipy.array([0,len(m),0,0,0,0]))
    NPartTotal = scipy.array([len(m)])
    MassArray = scipy.array(6*[0])
    Time = scipy.array([0.0])

    New = Gadget2.Gadget2Snapshot(x=x, y=y, z=z,vx=vx,vy=vy,vz=vz,ID=ID, m=m,NPartTotal=NPartTotal, NPart=NPart,Time=Time,MassArray=MassArray, V=None)

    return New
    #Gadget2.WriteGadget2(NewFilename , New)    