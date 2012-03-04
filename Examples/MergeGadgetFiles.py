import Classes.MergeHaloes
import scipy
from math import sqrt

def Norm(A):
    return sqrt(A[0]**2+A[1]**2+A[2]**2)
    

filename1 = '1HqIso_Impact0_121Centered.bin'
filename2 = 'Accretion_1600000.bin'
output = 'RemnantsAnd1600000.bin'


DeltaPos = scipy.array([0,0,0])
DeltaVel = scipy.array([0.0,0,0])

Classes.MergeHaloes.MergeHaloes(filename1,filename2,output ,DeltaPos = DeltaPos, DeltaVel=DeltaVel,ImportPotential=False)
