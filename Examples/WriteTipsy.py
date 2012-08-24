#/usr/bin/python
import matplotlib.pyplot as plt
import Classes.Gadget2 as Gadget2
import Classes.DM_structure as DM_structure
import scipy
import pylab
import copy,sys
import scipy.interpolate,scipy.optimize
        
FileName = '/home/ms/Uni/DarkMatter/AllSimulations/HJS2010/s4_16'
FileName = '/home/ms/Uni/DarkMatter/AllSimulations/InhomogeneousInfall/Infall_sub_1e6.bin_201'

A = DM_structure.DM_structure(FileName)
A.FindCenter()
A.FindCenterVel()
A.CenterParticlePositions()
A.CenterParticleVelocities()
A.Snapshot.WriteTipsy('test1.txt')
