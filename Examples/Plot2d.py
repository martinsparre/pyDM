#/usr/bin/python
import matplotlib.pyplot as plt
import Classes.DM_structure as DM_structure
import scipy
import sys
from mayavi import mlab

def Undersample(P, N=25000):
    "P: List with particle positions and others: [x,y,z, vx,vy,vz,V]"
    #Check that the attributes have the same length:

    if len(P)>0:
        Len = len(P[0])
    else:
        print 'P has length 0.'
        return P

    for i in range(len(P)):
        if len(P[i]) != Len:
            print 'P['+str(i)+'] doesnt have the same length as P[0]'
            return P
    
    #Choose points to sample:
    rand = scipy.random.random_integers(0,len(P[0])-1,N)

    newP = []
    for i in range(len(P)):
        newP.append(P[i][rand])

    return newP

if len(sys.argv) < 2:
    print 'No argument'
FileName=sys.argv[1]
A = DM_structure.DM_structure(FileName)
#A.Snapshot.SelectParticlesInCone(1,0,0,3.1415/8.0)

#mlab.figure(figure='Snapshot',bgcolor=(0,0,0),size=(600,600))
#mlab.clf()
Points=[A.Snapshot.x,A.Snapshot.y,A.Snapshot.z,A.Snapshot.V]    
x,y,z,V=Undersample(Points)
#mlab.points3d(x,y,z,V,scale_mode='none',scale_factor=0.119,resolution=3,opacity=0.7)

plt.figure(figsize=(12,12))
plt.subplot(2,2,1)

plt.xlabel('x')
plt.ylabel('y')

plt.axis('equal')
plt.xlim((-10,10))
plt.ylim((-10,10))


plt.plot(Points[0],Points[1],'.',mew=0.0,ms=0.8)


plt.subplot(2,2,2)

plt.plot(Points[0],Points[2],'.',mew=0.0,ms=0.8)
plt.xlabel('x')
plt.ylabel('z')
plt.xlim((-10,10))
plt.ylim((-10,10))


plt.subplot(2,2,3)
plt.plot(Points[1],Points[2],'.',mew=0.0,ms=0.8)
plt.xlabel('y')
plt.ylabel('z')
plt.xlim((-10,10))
plt.ylim((-10,10))

plt.show()
