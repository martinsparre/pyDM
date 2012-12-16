#/usr/bin/python
import matplotlib.pyplot as plt
import Classes.Gadget2 as Gadget2
import Classes.DM_structure as DM_structure
import scipy
import pylab
import copy,sys
import scipy.interpolate,scipy.optimize



#OutputName=['0.1','0.3','0.5']
OutputName=['1.5']
#Filenames = ['0.1/1HqIso_Impact0_160','0.3/1HqIso_Impact0_160','0.5/1HqIso_Impact0_160']
Filenames = ['1.5/1HqIso_Impact0_160']
NFiles = len(Filenames)
DIR = '/home/ms/Uni/DarkMatter/AllSimulations/SigmaAlignment2013/HeadonMerger_VaryingVel/'

for i in range(NFiles):
    FileName=DIR + Filenames[i]
    A = DM_structure.DM_structure(FileName)
    A.FindCenter()
    A.FindCenterVel()
    A.CenterParticlePositions()
    A.CenterParticleVelocities()
    A.Snapshot.WriteTipsy(OutputName[i] + '.txt')
