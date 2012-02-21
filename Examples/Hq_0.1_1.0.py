import Classes.MergeHaloes
import scipy
from math import sqrt

def Norm(A):
    return sqrt(A[0]**2+A[1]**2+A[2]**2)
    

filename1 = 'HqIso_000'
filename2 = 'Hq0.1.bin'
output = 'Hq_0.1_1.0.bin'

r=10.0#initial distance is 10
r_vec = scipy.array([1,0,0])
r_vec = r * r_vec/Norm(r_vec)


DeltaPos = r_vec
DeltaVel= scipy.array([-0.4264,0,0])

Classes.MergeHaloes.MergeHaloes(filename1,filename2,output ,DeltaPos = DeltaPos, DeltaVel=DeltaVel,ImportPotential=False)
