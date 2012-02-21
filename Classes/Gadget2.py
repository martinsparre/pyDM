#!/usr/bin/python
import array
import scipy
import os.path

class Gadget2Snapshot:
    """Data container."""
    def __init__(self,FileName=None, x=None,y=None,z=None,vx=None,vy=None,vz=None,ID=None,m=None,V=None,NPart=None,NPartTotal=None,MassArray=None,Time=None):
        self.FileName = FileName
        self.x = scipy.array(x)
        self.y = scipy.array(y)
        self.z = scipy.array(z)

        self.vx = scipy.array(vx)
        self.vy = scipy.array(vy)
        self.vz = scipy.array(vz)

        self.ID = scipy.array(ID)

        self.m = scipy.array(m)

        if V!=None:
            self.V = scipy.array(V)
        else:
            self.V = V

        self.NPart = scipy.array(NPart)
        self.NPartTotal = scipy.array(NPartTotal)
        self.MassArray = scipy.array(MassArray)
        self.Time = scipy.array(Time)
        
    def SelectIDsMax(self,IDMax):
        GoodIds = scipy.where(self.ID<IDMax)
        print GoodIds, len(self.x)
        self.x = self.x[GoodIds]
        self.y = self.y[GoodIds]
        self.z = self.z[GoodIds]
        self.vx = self.vx[GoodIds]
        self.vy = self.vy[GoodIds]
        self.vz = self.vz[GoodIds]
        self.m = self.m[GoodIds]
        self.ID = self.ID[GoodIds]
        self.V = self.V[GoodIds]
        self.NPartTotal = len(self.x)
        self.NPart = [0,len(self.x),0,0,0,0]
        
    
    def SelectIDsMin(self,IDMin):
        GoodIds = scipy.where(self.ID>IDMin)
        self.x = self.x[GoodIds]
        self.y = self.y[GoodIds]
        self.z = self.z[GoodIds]
        self.vx = self.vx[GoodIds]
        self.vy = self.vy[GoodIds]
        self.vz = self.vz[GoodIds]
        self.m = self.m[GoodIds]
        self.ID = self.ID[GoodIds]
        self.V = self.V[GoodIds]
        self.NPartTotal = len(self.x)
        self.NPart = [0,len(self.x),0,0,0,0]
        
    def CenterPos(self,dx,dy,dz):
        self.x -= dx
        self.y -= dy
        self.z -= dz

    def CenterVel(self,dvx,dvy,dvz):
        self.vx -= dvx
        self.vy -= dvy
        self.vz -= dvz
    
    def SelectParticlesInCone(self,x,y,z,OpeningAngle):
        Theta = scipy.arccos((self.x*x+self.y*y+self.z*z) / scipy.sqrt(x**2+y**2+z**2) / scipy.sqrt(self.x**2+self.y**2+self.z**2))
        GoodIds = scipy.where(Theta<OpeningAngle)
        self.x = self.x[GoodIds]
        self.y = self.y[GoodIds]
        self.z = self.z[GoodIds]
        self.vx = self.vx[GoodIds]
        self.vy = self.vy[GoodIds]
        self.vz = self.vz[GoodIds]
        self.m = self.m[GoodIds]
        self.ID = self.ID[GoodIds]
        self.V = self.V[GoodIds]
        self.NPartTotal = len(self.x)
        self.NPart = [0,len(self.x),0,0,0,0]

    def WriteToAscii(self,FileName):
        f=open(FileName,'w+')
        
        self.NPartTotal.tofile(f,sep='\n')
        f.write('\n')
        self.x.tofile(f,sep='\n')
        f.write('\n')
        self.NPartTotal.tofile(f,sep='\n')
        f.write('\n')
        self.y.tofile(f,sep='\n')
        f.write('\n')
        self.NPartTotal.tofile(f,sep='\n')
        f.write('\n')
        self.z.tofile(f,sep='\n')
        f.write('\n')
        self.NPartTotal.tofile(f,sep='\n')
        f.write('\n')
        self.vx.tofile(f,sep='\n')
        f.write('\n')
        self.NPartTotal.tofile(f,sep='\n')        
        f.write('\n')
        self.vy.tofile(f,sep='\n')
        f.write('\n')
        self.NPartTotal.tofile(f,sep='\n')        
        f.write('\n')
        self.vz.tofile(f,sep='\n')
        f.write('\n')
        self.NPartTotal.tofile(f,sep='\n')        
        f.write('\n')
        self.m.tofile(f,sep='\n')
        f.write('\n')
        self.NPartTotal.tofile(f,sep='\n')        
        f.write('\n')
        self.ID.tofile(f,sep='\n')
        f.write('\n')
        self.NPartTotal.tofile(f,sep='\n')        
        f.write('\n')
        self.V.tofile(f,sep='\n')
        f.write('\n')
        self.NPartTotal.tofile(f,sep='\n')
        f.write('\n')
        f.close()
        
        
        

def ReadBinary(GadgetFile, Type, N):
    """GadgetFile: e.g. f1 = open('tmp000.bin','rb')
    Type: e.g. 'f' for float. Look here: http://docs.python.org/library/array.html
    N: Number of numbers/instances to read

    will return the array
    """
    tmp_array = array.array(Type)
    tmp_array.fromfile(GadgetFile, N)
    tmp_array = scipy.array(tmp_array)
    return tmp_array

def WriteBinary(GadgetFile, Type, Array):
    """GadgetFile: e.g. f1 = open('tmp000.bin','wb')
    Type: e.g. 'f' for float. Look here: http://docs.python.org/library/array.html
    Array: scipy-array or list with the values to write...

    will return the array
    """
    tmp_array = array.array(Type, Array)
    tmp_array.tofile(GadgetFile)


def ReadGadget2(FileName,ImportPotential=True):
    if os.path.exists(FileName) == False:
        print 'File: ' + FileName + ' doesnt exist.'

    f = open(FileName,'rb')

    #read (...some parts of...) the header:
    ReadBinary(f,'l',1)

    NPart = ReadBinary(f,'I',6)
    MassArray = ReadBinary(f,'d',6)
    Time = ReadBinary(f,'d',1)

    NPartTotal = sum(NPart)


    #Skip the rest of the header:
    f.seek(0)

    LengthOfBlock0 = ReadBinary(f,'l',1)
    f.read(256)
    LengthOfBlock1 = ReadBinary(f,'l',1)

    if LengthOfBlock0 != LengthOfBlock1:
        print 'Inconsistent block length in file ' + FileName

    #Positions and velocities:
    print 'start pos'
    LengthOfBlock0 = ReadBinary(f,'l',1)

    tmp = ReadBinary(f,'f',3*NPartTotal)
    x=tmp[0::3]
    y=tmp[1::3]
    z=tmp[2::3]

    x=scipy.array(x)
    y=scipy.array(y)
    z=scipy.array(z)

#    x = ReadBinary(f,'f',NPartTotal)
#    y = ReadBinary(f,'f',NPartTotal)
#    z = ReadBinary(f,'f',NPartTotal)


    LengthOfBlock1 = ReadBinary(f,'l',1)

    if LengthOfBlock0 != LengthOfBlock1:
        print 'Inconsistent block length in file ' + FileName


    print 'start vel'
    LengthOfBlock0 =  ReadBinary(f,'l',1)

    tmp = ReadBinary(f,'f',3*NPartTotal)
    vx=tmp[0::3]
    vy=tmp[1::3]
    vz=tmp[2::3]

    vx=scipy.array(vx)
    vy=scipy.array(vy)
    vz=scipy.array(vz)


    LengthOfBlock1 = ReadBinary(f,'l',1)

    if LengthOfBlock0 != LengthOfBlock1:
        print 'Inconsistent block length in file ' + FileName
        
    print 'start ID'
    LengthOfBlock0 =  ReadBinary(f,'l',1)
    ID = ReadBinary(f,'I',NPartTotal)
    LengthOfBlock1 = ReadBinary(f,'l',1)

    if LengthOfBlock0 != LengthOfBlock1:
        print 'Inconsistent block length in file ' + FileName

    print 'start mass'
    LengthOfBlock0 =  ReadBinary(f,'l',1)
    m = ReadBinary(f,'f',NPartTotal)
    LengthOfBlock1 = ReadBinary(f,'l',1)

    #print m[0]

    if LengthOfBlock0 != LengthOfBlock1:
        print 'Inconsistent block length in file ' + FileName

    if  ImportPotential:
        print 'start potential'
        LengthOfBlock0 =  ReadBinary(f,'l',1)
        V = ReadBinary(f,'f',NPartTotal)
        LengthOfBlock1 = ReadBinary(f,'l',1)

        if LengthOfBlock0 != LengthOfBlock1:
            print 'Inconsistent block length in file ' + FileName
    else:
        V=None

    return Gadget2Snapshot(FileName=FileName, x=x,y=y,z=z,vx=vx,vy=vy,vz=vz,ID=ID,m=m,V=V,NPart=NPart,NPartTotal=NPartTotal,MassArray=MassArray,Time=Time)


def WriteGadget2(FileName, G):
    "G:Gadget snapshot class"
    f = open(FileName,'wb+')
    Count = 0
    
    WriteBinary(f, 'l', scipy.array([256]))
    WriteBinary(f, 'I', scipy.array(G.NPart))
    WriteBinary(f, 'd', scipy.array(G.MassArray))
    WriteBinary(f, 'd', scipy.array(G.Time))
    WriteBinary(f, 'l', scipy.array([0]*44) )#44*4 bytes remains in the header...
    WriteBinary(f, 'l', scipy.array([256]))

    if len(G.x) != len(G.y) != len(G.z): 
        print 'something wrong',len(G.x),len(G.y),len(G.z) 
    WriteBinary(f, 'l', scipy.array([3*4*G.NPartTotal]))
    WriteBinary(f, 'f', scipy.array([G.x,G.y,G.z]).transpose().flatten())
    WriteBinary(f, 'l', scipy.array([3*4*G.NPartTotal]))





    WriteBinary(f, 'l', scipy.array([3*4*G.NPartTotal]))
    WriteBinary(f, 'f', scipy.array([G.vx,G.vy,G.vz]).transpose().flatten())
    WriteBinary(f, 'l', scipy.array([3*4*G.NPartTotal]))

    WriteBinary(f, 'l', scipy.array([4*G.NPartTotal]))
    WriteBinary(f, 'I', scipy.array(G.ID))
    WriteBinary(f, 'l', scipy.array([4*G.NPartTotal]))

    WriteBinary(f, 'l', scipy.array([4*G.NPartTotal]))
    WriteBinary(f, 'f', scipy.array(G.m))
    WriteBinary(f, 'l', scipy.array([4*G.NPartTotal]))


    if G.V != None:
#        print 'Writing potentials...'
        WriteBinary(f, 'l', scipy.array([4*G.NPartTotal]))
        WriteBinary(f, 'f', scipy.array(G.V))
        WriteBinary(f, 'l', scipy.array([4*G.NPartTotal]))
    

    f.close()
