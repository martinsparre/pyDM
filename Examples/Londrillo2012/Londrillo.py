#/usr/bin/python
import matplotlib.pyplot as plt
import Classes.DM_structure as DM_structure
import scipy
import Classes.Gadget2 as Gadget2

FileName='Hq1e6_000'
A = DM_structure.DM_structure(FileName)
S = A.Snapshot
R = scipy.sqrt( S.x**2 + S.y**2 + S.z**2 )
Beta = 0.05
S.vx *= scipy.sqrt(Beta)
S.vy *= scipy.sqrt(Beta)
S.vz *= scipy.sqrt(Beta)
T = scipy.sum(0.5*S.m*(S.vx**2+S.vy**2+S.vz**2))
U = 0.5*scipy.sum(S.m*S.V)
print 2*T/scipy.fabs(U)
Gadget2.WriteGadget2(FileName + '_Beta' + str(Beta) , S )
print 'Created file',FileName + '_Beta' + str(Beta)
#plt.semilogx()
#plt.plot(R,S.V)
#plt.show()