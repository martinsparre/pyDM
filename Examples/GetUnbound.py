#/usr/bin/python
import matplotlib.pyplot as plt
import Classes.DM_structure as DM_structure
import scipy
import sys

FileName=sys.argv[1]
A = DM_structure.DM_structure(FileName)

E=0.5 * (A.Snapshot.vx*A.Snapshot.vx+A.Snapshot.vy*A.Snapshot.vy+A.Snapshot.vz*A.Snapshot.vz) + A.Snapshot.V
UnboundIDs = scipy.where((E>0.0))#*(A.Snapshot.ID<1000000))
plt.hist(A.Snapshot.ID[UnboundIDs],bins=11,range=(1,1100000))
plt.show()

plt.hist(scipy.log10(E[UnboundIDs]),bins=100)
plt.semilogy()
plt.show()

plt.plot(A.Snapshot.x[UnboundIDs],A.Snapshot.y[UnboundIDs],'.',color='black')
plt.show()

f=open('ToSteen.txt','w+')

f.write(str(len(UnboundIDs[0]))+'\n')
for i in UnboundIDs[0]:
    f.write(str(A.Snapshot.x[i])+'\t'+str(A.Snapshot.y[i])+'\t'+str(A.Snapshot.z[i])+'\n')    
f.close()

#StupidIDs = scipy.where((E>0.0)*(A.Snapshot.ID<1000000))
#B = DM_structure.DM_structure('../Hq_0.1_1.0/Hq_0.1_1.0.bin_080')
#plt.plot(B.Snapshot.x[StupidIDs],B.Snapshot.y[StupidIDs],'.')
#plt.show()
