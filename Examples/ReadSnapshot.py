#/usr/bin/python
import matplotlib.pyplot as plt
import Classes.DM_structure as DM_structure
import scipy
import copy

FileName = '/home/ms/Uni/DarkMatter/AllSimulations/SigmaAlignment2013/HeadonMerger_VaryingVel/0.1/1HqIso_Impact0_000'
A = DM_structure.DM_structure(FileName)
S = A.Snapshot

print 'Halo 1'
S.SelectIDsMin(1e6)
print '<x>',scipy.median(S.x)
print '<vx>',scipy.median(S.vx)



FileName = '/home/ms/Uni/DarkMatter/AllSimulations/SigmaAlignment2013/HeadonMerger_VaryingVel/0.1/1HqIso_Impact0_000'
A = DM_structure.DM_structure(FileName)
S = A.Snapshot


print 'Halo 2'

S.SelectIDsMax(1e6)
print '<x>',scipy.median(S.x)
print '<vx>',scipy.median(S.vx)
