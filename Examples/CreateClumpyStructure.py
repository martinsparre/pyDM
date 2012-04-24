#/usr/bin/python
import matplotlib.pyplot as plt
import Classes.DM_structure as DM_structure
import Classes.MergeHaloes
import Classes.Gadget2
import scipy
import copy

from mayavi import mlab

FileName='05small.bin'
Small = DM_structure.DM_structure(FileName, ImportPotential=False)

FileName='05.bin'
Main = DM_structure.DM_structure(FileName, ImportPotential=False)


Small.Snapshot.m = 0.0*Small.Snapshot.m + Main.Snapshot.m[0]

Small.Snapshot.vx = Small.Snapshot.vx*0.0
Small.Snapshot.vy = Small.Snapshot.vy*0.0
Small.Snapshot.vz = Small.Snapshot.vz*0.0

Main.Snapshot.vx = Main.Snapshot.vx*0.0
Main.Snapshot.vy = Main.Snapshot.vy*0.0
Main.Snapshot.vz = Main.Snapshot.vz*0.0

NewPos = copy.deepcopy(Main)
NewPos.Snapshot.SelectIDsMax(5e5+1025)
NewPos.Snapshot.SelectIDsMin(5e5+1000)

Main.Snapshot.SelectIDsMax(5e5)


NewPos = NewPos.Snapshot
New = Main.Snapshot
Small = Small.Snapshot
Main = Main.Snapshot

for i in range(24):
    x,y,z = NewPos.x[i],NewPos.y[i],NewPos.z[i]
    print x,y,z
    
    NewSmall = copy.deepcopy(Small)
    NewSmall.x += x
    NewSmall.y += y
    NewSmall.z += z
    
    New = Classes.MergeHaloes.MergeSnapshots(New,NewSmall,ImportPotential=False)
    print 'Length of new file',len(New.x)

Classes.Gadget2.WriteGadget2('05ColdCollapse.bin', New)