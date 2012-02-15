#/usr/bin/python
import Classes.Gadget2 as Gadget2
import sys

if len(sys.argv) != 3:
    print 'Proper use:\n python Gadget2ToAscii.py GadgetFileName OutputFile.txt'


f=open(sys.argv[1],'r')

import scipy
N = scipy.fromfile(f, dtype=int, count=1, sep='\n')[0]

x = scipy.fromfile(f, count=N, sep='\n')
if N != scipy.fromfile(f, dtype=int, count=1, sep='\n')[0]:
    print 'Something wrong with format'

y = scipy.fromfile(f, count=N, sep='\n')
if N != scipy.fromfile(f, dtype=int, count=1, sep='\n')[0]:
    print 'Something wrong with format'

z = scipy.fromfile(f, count=N, sep='\n')
if N != scipy.fromfile(f, dtype=int, count=1, sep='\n')[0]:
    print 'Something wrong with format'


vx = scipy.fromfile(f, count=N, sep='\n')
if N != scipy.fromfile(f, dtype=int, count=1, sep='\n')[0]:
    print 'Something wrong with format'

vy = scipy.fromfile(f, count=N, sep='\n')
if N != scipy.fromfile(f, dtype=int, count=1, sep='\n')[0]:
    print 'Something wrong with format'

vz = scipy.fromfile(f, count=N, sep='\n')
if N != scipy.fromfile(f, dtype=int, count=1, sep='\n')[0]:
    print 'Something wrong with format'

m = scipy.fromfile(f, count=N, sep='\n')
if N != scipy.fromfile(f, dtype=int, count=1, sep='\n')[0]:
    print 'Something wrong with format'

ID = scipy.fromfile(f, count=N,dtype=int, sep='\n')
if N != scipy.fromfile(f, dtype=int, count=1, sep='\n')[0]:
    print 'Something wrong with format'

V = scipy.fromfile(f, count=N, sep='\n')
if N != scipy.fromfile(f, dtype=int, count=1, sep='\n')[0]:
    print 'Something wrong with format'

G = Gadget2.Gadget2Snapshot(FileName=None, x=x,y=y,z=z,vx=vx,vy=vy,vz=vz,ID=ID,m=m,V=V,NPart=[0,N,0,0,0,0],NPartTotal=N,MassArray=[0.0]*6,Time=[0.0])


Gadget2.WriteGadget2(sys.argv[2], G)

#tmp = scipy.fromfile(f,)
#G = Gadget2.ReadGadget2(sys.argv[1])#'1HqIso_Impact0_000' is the gadget2-output
#G.WriteToAscii(sys.argv[2])



