#/usr/bin/python
import Gadget2
import os.path
import scipy
import scipy.interpolate,scipy.optimize
import math
import numpy.random

def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]


class EnergyGrid:
    def __init__(self):
        self.E = None
        self.Mass = None
    
    def Smooth(self,Width=0.1):
        self.Mass=scipy.array(pyROOT_functions.RootSpline(scipy.log10(self.E),self.Mass,BandPass=Width))[1]



def FindNearestIDs(x,y,y0):
    ID=(abs(y-y0)).argmin()
    if ID<3 or ID > len(x)-3:
        return [x[ID]],[y[ID]]
    else:
        return x[ID-3:ID+3],y[ID-3:ID+3]

def DoInterpolate(x,y,x0):
    y1,x1 = FindNearestIDs(y,x,x0)
    if len(x1)>1:
        spline = scipy.interpolate.interp1d(x1,y1,kind='cubic')
        return spline(x0)
    else:
        return x1[0]


def DoSolve(x,y,y0):
    x1,y1 = FindNearestIDs(x,y,y0)
    if len(x1)>1:
        spline = scipy.interpolate.interp1d(x1,y1-y0,kind='cubic')
        return scipy.optimize.fsolve(spline,x1[len(x1)/2])
    else:
        return x1[0]



def solve(x,y,y0):
    #"return x0, where y0=y(x0)"
    
    #if y0 < y.min():
        #y0=y.min()+0.05
    
    #if y0 > y.max():
        #y0=y.max()-0.05

    #print x
    #print y
    #print y0
    #spline = scipy.interpolate.interp1d(x,y,kind=3)
    #x = scipy.optimize.fsolve(spline,y0)

    return x[5]

def Interpolate(x,y,x0):
    spline = scipy.interpolate.interp1d(x,y,kind=3, bounds_error=False, fill_value = 0.0)
    return spline(x0)


class Grid:
    def __init__(self):
        self.R = []
        self.V = []
        self.Sigma2 = []
        self.Sigma2r = []
        self.Beta = []
        self.UncBetaBootstrap = []
        self.MeanBetaBootstrap = []        
        self.Gamma = []
        self.Kappa = []
        self.Rho = []
        self.MassInBin = []
        self.N = []
        self.JeansMass = []
        self.Rmax = []
        self.Rmin = []
        self.MeanVr = []
        self.CumulativeMass = []
        self.SigmaTensor = []
        self.Lx = []
        self.Ly = []
        self.Lz = []
        self.dnu2dpsi2_psi = None
        self.dnu2dpsi2 = None
        self.drhodr = None
        self.d2rhodr2 = None
        self.dVdr = None
        self.d2Vdr2 = None
        self.E = None
        self.DF = None
        print "\nGrid initialized\n"
    
    
    def Getd2nudpsi2(self,Psi):

        R = DoSolve(self.R,self.V,-Psi)
        #get r from psi
        M = DoInterpolate(self.R,self.CumulativeMass,R)

        Rho = DoInterpolate(self.R,self.Rho,R)
        drho = DoInterpolate(self.R,self.drhodr,R)
        d2rho = DoInterpolate(self.R,self.d2rhodr2,R)
        #print Psi,R,M,Rho,drho,d2rho
        G=1.0
        dpsi = - G * M / R**2
        d2psi = 2.*G*M / R**3 - 4.0 * 3.14159265359 * G * Rho
        
        temp = d2rho - drho / dpsi * d2psi
        temp /=  dpsi * dpsi
        
        return temp
    
    
    
    def CalcEddington(self,NPoints=50):
        print "\nCalcEddington started"
        PsiMin = (-self.V).min()
        PsiMax = (-self.V).max()
        print PsiMin, PsiMax
        dx = PsiMax-PsiMin
        E = scipy.linspace(PsiMin+dx/200.0,PsiMax-dx/200.0,NPoints)
        DF = []

        for i in range(len(E)):
            print E[i]
            umin = math.atan( 1. / math.sqrt( E[i] ) )
            

            N1 = int(300.0*float(i)/len(E))
            midp = 0.0 #midpoint integration
            
            if i==0:
                du = math.pi/2.0 - umin
                uj = umin+0.5*du
                midp = 2./ math.sin(uj)**2 * self.Getd2nudpsi2( E[i] - 1./math.tan(uj)**2 )
                midp /=math.sqrt(8.)*math.pi*math.pi/du
            else:
                
                du = (math.pi/2.0 - umin) / float(N1)

                for j in range(N1):
                    uj = umin + (j+0.5)*du
                    midp += 2./ math.sin(uj)**2 * self.Getd2nudpsi2( E[i] - 1./math.tan(uj)**2 )
                    
                midp /=math.sqrt(8.)*math.pi*math.pi/du
            
            #midp += - 1./(math.sqrt(8.)*math.pi*math.pi) * 1./sqrt(Cl->E[i]) * drhodr(Gr[NG-1].R,Cl) * Gr[NG-1].R*Gr[NG-1].R / G / Gr[NG-1].Mass;

            DF.append( midp )
        
        self.E = scipy.array(E)
        self.DF = scipy.array(DF)

        print "\nCalcEddington finished"
        return E,DF
        

    def CalculateRhoDerivatives(self,Nderiv=2):
        print "\nCalculateRhoDerivatives started"
        import pyROOT_functions
        drho = []
        for i in range(len(self.Rho)):
            if i<Nderiv or i>len(self.Rho)-1-Nderiv:
                drho.append(1.0)
            else:
                drho.append( (self.Rho[i+Nderiv]-self.Rho[i-Nderiv])/(self.R[i+Nderiv]-self.R[i-Nderiv])    )
        

        self.drhodr = scipy.array(drho)
            
        d2rho = []
        for i in range(len(self.drhodr)):
            if i<Nderiv or i>len(self.drhodr)-1-Nderiv:
                d2rho.append(1.0)
            else:
                d2rho.append( (self.drhodr[i+Nderiv]-self.drhodr[i-Nderiv])/(self.R[i+Nderiv]-self.R[i-Nderiv])    )
        self.d2rhodr2 = scipy.array(d2rho)
        
        print "CalculateRhoDerivatives ended\n"



    def CalculatePotentialDerivatives(self,Nderiv=2,Log10Width=0.2):
        import pyROOT_functions
        dV = []
        for i in range(len(self.V)):
            if i<Nderiv:
                dV.append( (self.V[i+Nderiv]-self.V[i])/(self.R[i+Nderiv]-self.R[i])    )
            elif i>len(self.V)-1-Nderiv:
                dV.append( (self.V[i-Nderiv]-self.V[i])/(self.R[i-Nderiv]-self.R[i])    )
            else:
                dV.append( (self.V[i+Nderiv]-self.V[i-Nderiv])/(self.R[i+Nderiv]-self.R[i-Nderiv])    )
        

        self.dVdr = scipy.array(dV)
            

        d2V = []
        for i in range(len(self.V)):
            if i<Nderiv:
                d2V.append( (self.dVdr[i+Nderiv]-self.dVdr[i])/(self.R[i+Nderiv]-self.R[i])    )
            elif i>len(self.V)-1-Nderiv:
                d2V.append( (self.dVdr[i]-self.dVdr[i-Nderiv])/(self.R[i]-self.R[i-Nderiv])    )
            else:
                d2V.append( (self.dVdr[i+Nderiv]-self.dVdr[i-Nderiv])/(self.R[i+Nderiv]-self.R[i-Nderiv])    )
        self.d2Vdr2 = scipy.array(d2V)


    def SaveToBinary(self,FileName):
        f = open(FileName,'wb+')
        N=len(self.R)
        Gadget2.WriteBinary(f, 'l', scipy.array([N]))
        Gadget2.WriteBinary(f, 'f', scipy.array(self.R))
        Gadget2.WriteBinary(f, 'l', scipy.array([N]))

        Gadget2.WriteBinary(f, 'l', scipy.array([N]))
        Gadget2.WriteBinary(f, 'f', scipy.array(self.V))
        Gadget2.WriteBinary(f, 'l', scipy.array([N]))
        
        Gadget2.WriteBinary(f, 'l', scipy.array([N]))
        Gadget2.WriteBinary(f, 'f', scipy.array(self.Sigma2))
        Gadget2.WriteBinary(f, 'l', scipy.array([N]))

        Gadget2.WriteBinary(f, 'l', scipy.array([N]))
        Gadget2.WriteBinary(f, 'f', scipy.array(self.Sigma2r))
        Gadget2.WriteBinary(f, 'l', scipy.array([N]))
        
        Gadget2.WriteBinary(f, 'l', scipy.array([N]))
        Gadget2.WriteBinary(f, 'f', scipy.array(self.Beta))
        Gadget2.WriteBinary(f, 'l', scipy.array([N]))

        Gadget2.WriteBinary(f, 'l', scipy.array([N]))
        Gadget2.WriteBinary(f, 'f', scipy.array(self.Gamma))
        Gadget2.WriteBinary(f, 'l', scipy.array([N]))
        
        Gadget2.WriteBinary(f, 'l', scipy.array([N]))
        Gadget2.WriteBinary(f, 'f', scipy.array(self.Kappa))
        Gadget2.WriteBinary(f, 'l', scipy.array([N]))

        Gadget2.WriteBinary(f, 'l', scipy.array([N]))
        Gadget2.WriteBinary(f, 'f', scipy.array(self.Rho))
        Gadget2.WriteBinary(f, 'l', scipy.array([N]))
        
        Gadget2.WriteBinary(f, 'l', scipy.array([N]))
        Gadget2.WriteBinary(f, 'f', scipy.array(self.CumulativeMass))
        Gadget2.WriteBinary(f, 'l', scipy.array([N]))
        
        f.close()


    def SaveToAscii(self,FileName):
        f = open(FileName,'w+')
        N=len(self.R)
        f.write(str(N)+'\n')

        for i in self.R:
            f.write(format(i,'.20e')+'\n')
        f.write(str(N)+'\n')
        f.write(str(N)+'\n')
            
        for i in self.Rho:
            f.write(format(i,'.20e')+'\n')
        f.write(str(N)+'\n')
        f.write(str(N)+'\n')
        for i in self.Sigma2:
            f.write(format(i,'.20e')+'\n')
        f.write(str(N)+'\n')
        f.write(str(N)+'\n')
        for i in self.Sigma2r:
            f.write(format(i,'.20e')+'\n')
        f.write(str(N)+'\n')
        f.write(str(N)+'\n')
        for i in self.Beta:
            f.write(format(i,'.20e')+'\n')
        f.write(str(N)+'\n')
        f.write(str(N)+'\n')
        for i in self.Gamma:
            f.write(format(i,'.20e')+'\n') 
        f.write(str(N)+'\n')
        f.write(str(N)+'\n')
        for i in self.Kappa:
            f.write(format(i,'.20e')+'\n')
        f.write(str(N)+'\n')
        f.write(str(N)+'\n')
        for i in self.V:
            f.write(format(i,'.20e')+'\n')
        f.write(str(N)+'\n')
        f.write(str(N)+'\n')
        for i in self.CumulativeMass:
            f.write(format(i,'.20e')+'\n')
        f.write(str(N)+'\n')
        
        if self.d2rhodr2 != None and self.drhodr != None:
            f.write(str(N)+'\n')
            for i in self.drhodr:
                f.write(format(i,'.20e')+'\n')
            f.write(str(N)+'\n')
            f.write(str(N)+'\n')
            for i in self.d2rhodr2:
                f.write(format(i,'.20e')+'\n')
            f.write(str(N)+'\n')

        f.close()


    def Smooth(self,Log10Widths=[0.3,0.1],ChangeN=200):
        print '\nSmoothing started'
        import pyROOT_functions
        
        
        #smoothing in the inner part:
        
        Log10Width=Log10Widths[0]
        V0=scipy.array(pyROOT_functions.RootSpline(scipy.log10(self.R),self.V,BandPass=Log10Width))[1]
        Rho0=scipy.array(pyROOT_functions.RootSpline(scipy.log10(self.R),self.Rho,BandPass=Log10Width))[1]
        Sigma2r0=scipy.array(pyROOT_functions.RootSpline(scipy.log10(self.R),self.Sigma2r,BandPass=Log10Width))[1]
        Sigma20=scipy.array(pyROOT_functions.RootSpline(scipy.log10(self.R),self.Sigma2,BandPass=Log10Width))[1]

        #smoothing in the outer part:
        Log10Width=Log10Widths[1]
        V1=scipy.array(pyROOT_functions.RootSpline(scipy.log10(self.R),self.V,BandPass=Log10Width))[1]
        Rho1=scipy.array(pyROOT_functions.RootSpline(scipy.log10(self.R),self.Rho,BandPass=Log10Width))[1]
        Sigma2r1=scipy.array(pyROOT_functions.RootSpline(scipy.log10(self.R),self.Sigma2r,BandPass=Log10Width))[1]
        Sigma21=scipy.array(pyROOT_functions.RootSpline(scipy.log10(self.R),self.Sigma2,BandPass=Log10Width))[1]
        
        #create final array:
        def Merge(y0,y1,ChangeN):
            y2 = scipy.append(scipy.array(y0[0:ChangeN]),scipy.array(y1[ChangeN:len(y0)]))
            if len(y2) != len(y0) or len(y2) != len(y1):
                print 'something is wrong!123u4njdsf'
            return y2
        
        print 'log-Radius at smooth-change:',scipy.log10(self.R[ChangeN])
        
        self.V = Merge(V0,V1,ChangeN)
        self.Rho = Merge(Rho0,Rho1,ChangeN)
        self.Sigma2r = Merge(Sigma2r0,Sigma2r1,ChangeN)
        self.Sigma2 = Merge(Sigma20,Sigma21,ChangeN)
        
        print self.Rho[ChangeN-4:ChangeN+4]

  #      self.drhodr=scipy.array(pyROOT_functions.RootSpline(scipy.log10(self.R),self.drhodr,BandPass=Log10Width))[1]
  #      self.d2rhodr2=scipy.array(pyROOT_functions.RootSpline(scipy.log10(self.R),self.d2rhodr2,BandPass=Log10Width))[1]

        gamma=[]
        kappa=[]
        gamma.append(0.5*(self.R[0]/self.Rho[0]+self.R[1]/self.Rho[1])*(self.Rho[1]-self.Rho[0])/(self.R[1]-self.R[0]))
        kappa.append(0.5*(self.R[0]/self.Sigma2r[0]+self.R[1]/self.Sigma2r[1])*(self.Sigma2r[1]-self.Sigma2r[0])/(self.R[1]-self.R[0]))
        for i in range(1,len(self.R)-1):
            gamma.append(self.R[i]/self.Rho[i]*(self.Rho[i+1]-self.Rho[i-1])/(self.R[i+1]-self.R[i-1]))
            kappa.append(self.R[i]/self.Sigma2r[i]*(self.Sigma2r[i+1]-self.Sigma2r[i-1])/(self.R[i+1]-self.R[i-1]))

        gamma.append(0.5*(self.R[-2]/self.Rho[-2]+self.R[-1]/self.Rho[-1])*(self.Rho[-1]-self.Rho[-2])/(self.R[-1]-self.R[-2]))
        kappa.append(0.5*(self.R[-2]/self.Sigma2r[-2]+self.R[-1]/self.Sigma2r[-1])*(self.Sigma2r[-1]-self.Sigma2r[-2])/(self.R[-1]-self.R[-2]))
        
        self.Gamma = scipy.array(gamma)
        self.Kappa = scipy.array(kappa)
        self.Beta = 1.0 - (self.Sigma2-self.Sigma2r)/(2.0*self.Sigma2r)
        print '\nSmoothing finished'


    def RemoveBadRegions(self):
        print "\nRemoveBadRegions started"
        
        #Grid-ids to be removed
        ID = [0,1,2,3,len(self.R)-1,len(self.R)-2,len(self.R)-3,len(self.R)-4] #endpoints are annoying....

        #check where gamma is bad:
        for i in scipy.where(scipy.isnan(self.Kappa)==True)[0]:
            ID.append(i)

#        print 'A',len(self.R)

        self.R = scipy.delete(self.R,ID)
        self.V = scipy.delete(self.V,ID)
        self.Sigma2 = scipy.delete(self.Sigma2,ID)
        self.Sigma2r = scipy.delete(self.Sigma2r,ID)
        self.Beta = scipy.delete(self.Beta,ID)
        self.Gamma = scipy.delete(self.Gamma,ID)
        self.Kappa = scipy.delete(self.Kappa,ID)
        self.Rho = scipy.delete(self.Rho,ID)
        self.MassInBin = scipy.delete(self.MassInBin,ID)
        self.N = scipy.delete(self.N,ID)
        self.Rmax = scipy.delete(self.Rmax,ID)
        self.Rmin = scipy.delete(self.Rmin,ID)
        self.MeanVr = scipy.delete(self.MeanVr,ID)
        self.CumulativeMass = scipy.delete(self.CumulativeMass,ID)
        
        if self.drhodr != None:
            self.drhodr =  scipy.delete(self.drhodr,ID)

        if self.d2rhodr2 != None:
            self.d2rhodr2 =  scipy.delete(self.d2rhodr2,ID)
        
        


        print "RemoveBadRegions finished\n"


class DM_structure:
    def __init__(self,Gadget2FileName,Snapshot=None,ImportPotential=True):
        if Snapshot == None:
            print 'initializing DM structure'
            if os.path.exists(Gadget2FileName) == False:
                print '--Error: File '+Gadget2FileName+' does not exist'
            print '-Will now read snapshot '+Gadget2FileName
            self.Snapshot = Gadget2.ReadGadget2(Gadget2FileName,ImportPotential=ImportPotential)
            if ImportPotential==False:
                print 'Warning: potentials are not imported... Some functions will not work'
            print '-Finished reading snapshot '+Gadget2FileName
        else:
            print "Taking snapshot from input..."
            self.Snapshot = Snapshot
        self.EGridLin = None
        self.EGrid = None
        self.GrSph = None
        self.Gr = None
        self.V_index = None
        self.x_C = 0.0
        self.y_C = 0.0
        self.z_C = 0.0
        self.Mean_vx = 0.0
        self.Mean_vy = 0.0
        self.Mean_vz = 0.0
        self.CenterFound = False
        self.CenterVelFound = False


    def GetSnapshot(self):
        if self.Snapshot == None:
            print "-Warning: In GetGadget2Snapshot - snapshot not defined"
        
        return self.Snapshot
        

    def FindCenter(self,Median = True,N_Median=10):
        "define center to be the particle with the lowest potential"
        
        if self.Snapshot.V == None:
            print 'Potentials not defined. Can not find center'
            return None
        self.V_index = scipy.argsort(self.Snapshot.V)
        if Median == False:
            self.x_C = self.Snapshot.x[self.V_index[0]]
            self.y_C = self.Snapshot.y[self.V_index[0]]
            self.z_C = self.Snapshot.z[self.V_index[0]]
        else:
            self.x_C = scipy.median(self.Snapshot.x[self.V_index[0:N_Median]])
            self.y_C =  scipy.median(self.Snapshot.y[self.V_index[0:N_Median]])
            self.z_C = scipy.median( self.Snapshot.z[self.V_index[0:N_Median]])
        
        self.CenterFound = True
       
    def FindCenterVel(self, Median = True,N_Median = 1000):
        if Median == False:
            self.Mean_vx = scipy.mean(self.Snapshot.vx)
            self.Mean_vy = scipy.mean(self.Snapshot.vy)
            self.Mean_vz = scipy.mean(self.Snapshot.vz)
        else:
            if self.Snapshot.V == None:
                print 'Potentials not defined. Can not find center'
                return None

            V_index = scipy.argsort(self.Snapshot.V)
            self.Mean_vx = scipy.median(self.Snapshot.vx[V_index[0:N_Median]])
            self.Mean_vy =  scipy.median(self.Snapshot.vy[V_index[0:N_Median]])
            self.Mean_vz = scipy.median( self.Snapshot.vz[V_index[0:N_Median]])
        
        self.CenterVelFound = True


    def SetCenter(self,x,y,z):
        "Set center manually"
        self.x_C = x
        self.y_C = y
        self.z_C = z
        self.CenterFound = True
    
    def SetCenterVel(self,vx,vy,vz):
        self.Mean_vx = vx
        self.Mean_vy = vy
        self.Mean_vz = vz
        self.CenterVelFound = True       
    
    def CenterParticlePositions(self):
        if self.CenterFound:
            self.Snapshot.x -= self.x_C
            self.Snapshot.y -= self.y_C
            self.Snapshot.z -= self.z_C
        else:
            print 'Warning: Position center not defined'        

    def CenterParticleVelocities(self):
        if self.CenterVelFound:
            self.Snapshot.vx -= self.Mean_vx
            self.Snapshot.vy -= self.Mean_vy
            self.Snapshot.vz -= self.Mean_vz
        else:
            print 'Warning: Velocity center not defined'    
            
    def CenterSnapshotPosAndVel(self):
        self.Snapshot.CenterPos( self.x_C, self.y_C, self.z_C)
        self.Snapshot.CenterVel( self.Mean_vx, self.Mean_vy, self.Mean_vz)
    
    def GetParticlePosVel(self,pID):
        i = scipy.where(self.Snapshot.ID == pID)
        return [self.Snapshot.x[i],self.Snapshot.y[i],self.Snapshot.z[i],self.Snapshot.vx[i],self.Snapshot.vy[i],self.Snapshot.vz[i],self.Snapshot.ID[i] ]
        #for i in range(len(self.Snapshot.y)):
            #if self.Snapshot.ID[i] == pID:
                #return [self.Snapshot.x[i],self.Snapshot.y[i],self.Snapshot.z[i],self.Snapshot.vx[i],self.Snapshot.vy[i],self.Snapshot.vz[i],self.Snapshot.ID[i] ]
                

    def CreateGridLogBins(self,Rmin=False,Rmax=False,NBins=100,CalcVelDispTensor = False,UseSphericalCoordinates=True,CalcUncBeta=False):
        "Spherical bins, distributed equally in logspace"
        print "\nCreateGridLogBins started"
        r = scipy.sqrt(self.Snapshot.x*self.Snapshot.x+self.Snapshot.y*self.Snapshot.y+self.Snapshot.z*self.Snapshot.z)
        
        Index = scipy.argsort(r)
        if Rmin == False:
            Rmin = scipy.mean(r[Index[0:10]])

        if Rmax == False:
            Rmax = scipy.mean(r[Index[-10:-1]])
        
        LeftBinLimit = scipy.logspace(scipy.log10(Rmin),scipy.log10(Rmax),NBins)
        dlogx = scipy.log10(LeftBinLimit[1]/LeftBinLimit[0])
        BinNo = scipy.array(scipy.log10(r/LeftBinLimit[0])/dlogx,int)
        
        BinNo[scipy.where(BinNo<0)]=0
        BinNo[scipy.where(BinNo>NBins-1)]=NBins-1
        
        
        if UseSphericalCoordinates:
            Phi = scipy.arctan2(self.Snapshot.y,self.Snapshot.x)
            Theta = scipy.arccos(self.Snapshot.z/r)
            
            VR = scipy.sin(Theta)*scipy.cos(Phi) * self.Snapshot.vx  + scipy.sin(Theta)*scipy.sin(Phi) * self.Snapshot.vy + scipy.cos(Theta) * self.Snapshot.vz
            VTheta = scipy.cos(Theta)*scipy.cos(Phi) * self.Snapshot.vx + scipy.cos(Theta)*scipy.sin(Phi) * self.Snapshot.vy - scipy.sin(Theta) * self.Snapshot.vz
            VPhi = - scipy.sin(Phi) * self.Snapshot.vx +  scipy.cos(Phi) * self.Snapshot.vy
        else:
            v2 = self.Snapshot.vx*self.Snapshot.vx+self.Snapshot.vy*self.Snapshot.vy+self.Snapshot.vz*self.Snapshot.vz
            vr = 1.0/(r+1.0e-16)*scipy.array(self.Snapshot.vx*self.Snapshot.x + self.Snapshot.vy*self.Snapshot.y + self.Snapshot.vz*self.Snapshot.z)                    
#        vx1 = VR1 * scipy.sin(Theta1) * scipy.cos(Phi1) + VTheta1 * scipy.cos(Theta1) * scipy.cos(Phi1) - VPhi1 * scipy.sin(Phi1)
#        vy1 = VR1 * scipy.sin(Theta1) * scipy.sin(Phi1) + VTheta1 * scipy.cos(Theta1) * scipy.sin(Phi1) + VPhi1 * scipy.cos(Phi1)
#        vz1 = VR1 * scipy.cos(Theta1) - VTheta1 * scipy.sin(Theta1)
#        print self.Snapshot.vx,vx1
#        print 'A',max(abs(self.Snapshot.vx-vx1))
#        print 'B',max(abs(self.Snapshot.vy-vy1))
#        print 'C',max(abs(self.Snapshot.vz-vz1))

        
        self.GrSph = Grid()
        
        for Bin in range(NBins):
            Particles=scipy.where(BinNo==Bin)
            Mass = scipy.sum(self.Snapshot.m[Particles])
            N = len(Particles[0])

            RminBin = LeftBinLimit[Bin]
            
            if Bin != NBins-1:
                RmaxBin = LeftBinLimit[Bin+1]
            else:
                RmaxBin = 2*LeftBinLimit[Bin]-LeftBinLimit[Bin-1]

            Volume = 4.0/3.0 * 3.14159265359 * (RmaxBin**3 - RminBin**3)
            
            if self.Snapshot.V != None:
                V = scipy.mean(self.Snapshot.V[Particles])

                
            
            if UseSphericalCoordinates:
                Sigma2r = scipy.mean(VR[Particles]**2)-scipy.mean(VR[Particles])**2
                Sigma2Theta = scipy.mean(VTheta[Particles]**2)-scipy.mean(VTheta[Particles])**2
                Sigma2Phi = scipy.mean(VPhi[Particles]**2)-scipy.mean(VPhi[Particles])**2
                Sigma2 = Sigma2r + Sigma2Theta + Sigma2Phi
                Beta =  1.0 - (Sigma2Theta+Sigma2Phi)/(2.0*Sigma2r)
                MeanVr = scipy.mean(VR[Particles])


                if CalcUncBeta == True:
                    BetaBootstrap = []
                    for i in range(20):
                        if len(Particles[0])  < 2:
                            break
                            
#                        print 'A',Particles
#                        print 'A1',Particles[0]
#                        print 'B',len(Particles[0])
#                        print 'B1',numpy.random.randint(0,high=len(Particles[0]),size=len(Particles[0]))
#                        print 'B2',Particles[0][scipy.array(numpy.random.randint(0,high=len(Particles[0]),size=len(Particles[0])))]
#                        print 'C',len(Particles)                

                        NewParticles = Particles[0][numpy.random.randint(0,high=len(Particles[0]),size=len(Particles[0]))]

                        NewSigma2r = scipy.mean(VR[NewParticles]**2)-scipy.mean(VR[NewParticles])**2
                        NewSigma2Theta = scipy.mean(VTheta[NewParticles]**2)-scipy.mean(VTheta[NewParticles])**2
                        NewSigma2Phi = scipy.mean(VPhi[NewParticles]**2)-scipy.mean(VPhi[NewParticles])**2
                        NewSigma2 = NewSigma2r + NewSigma2Theta + NewSigma2Phi
                        NewBeta =  1.0 - (NewSigma2Theta+NewSigma2Phi)/(2.0*NewSigma2r)
                        BetaBootstrap.append( NewBeta )
                    
                    if len(Particles[0]) < 2:
                        self.GrSph.UncBetaBootstrap.append(0.0)
                        self.GrSph.MeanBetaBootstrap.append(0.0)
                    else:
                        self.GrSph.UncBetaBootstrap.append(scipy.std(BetaBootstrap))
                        self.GrSph.MeanBetaBootstrap.append(scipy.mean(BetaBootstrap))                
                
            else:
                MeanVr = scipy.mean(vr[Particles])
                Sigma2r = scipy.mean(vr[Particles]**2)
                Sigma2 = scipy.mean(v2[Particles])
                Beta = 1.0 - (Sigma2-Sigma2r)/(2.0*Sigma2r)
            

            if CalcVelDispTensor:
                Tensor = scipy.zeros((3,3))
                Tensor[0,0] = scipy.mean(self.Snapshot.vx[Particles] * self.Snapshot.vx[Particles]) - scipy.mean(self.Snapshot.vx[Particles])  *scipy.mean( self.Snapshot.vx[Particles])
                Tensor[1,1] = scipy.mean(self.Snapshot.vy[Particles] * self.Snapshot.vy[Particles]) - scipy.mean(self.Snapshot.vy[Particles])  *scipy.mean( self.Snapshot.vy[Particles])               
                Tensor[2,2] = scipy.mean(self.Snapshot.vz[Particles] * self.Snapshot.vz[Particles]) - scipy.mean(self.Snapshot.vz[Particles])  *scipy.mean( self.Snapshot.vz[Particles])
                Tensor[0,1] = scipy.mean(self.Snapshot.vx[Particles] * self.Snapshot.vy[Particles]) - scipy.mean(self.Snapshot.vx[Particles])  *scipy.mean( self.Snapshot.vy[Particles])
                Tensor[1,0] = Tensor[0,1] 
                Tensor[0,2] = scipy.mean(self.Snapshot.vx[Particles] * self.Snapshot.vz[Particles]) - scipy.mean(self.Snapshot.vx[Particles])  *scipy.mean( self.Snapshot.vz[Particles])
                Tensor[2,0] = Tensor[0,2]  
                Tensor[1,2] = scipy.mean(self.Snapshot.vy[Particles] * self.Snapshot.vz[Particles]) - scipy.mean(self.Snapshot.vy[Particles])  *scipy.mean( self.Snapshot.vz[Particles])
                Tensor[2,1] = Tensor[1,2]
                self.GrSph.SigmaTensor.append(Tensor)


            Lx = scipy.mean(abs(self.Snapshot.y[Particles]*self.Snapshot.vz[Particles]-self.Snapshot.z[Particles]*self.Snapshot.vy[Particles]))
            Ly = scipy.mean(abs(self.Snapshot.z[Particles]*self.Snapshot.vx[Particles]-self.Snapshot.x[Particles]*self.Snapshot.vz[Particles]))
            Lz = scipy.mean(abs(self.Snapshot.x[Particles]*self.Snapshot.vy[Particles]-self.Snapshot.y[Particles]*self.Snapshot.vx[Particles]))
            
            self.GrSph.Rmin.append(RminBin)
            self.GrSph.Rmax.append(RmaxBin)
            self.GrSph.MassInBin.append(Mass)
            self.GrSph.N.append(N) 
            self.GrSph.Rho.append(Mass / Volume)
            self.GrSph.Sigma2.append(Sigma2 )            
            self.GrSph.Sigma2r.append(Sigma2r )
            self.GrSph.Beta.append(Beta)
            if self.Snapshot.V != None:
                self.GrSph.V.append(V)
            self.GrSph.MeanVr.append(MeanVr)
            self.GrSph.CumulativeMass.append(sum(self.GrSph.MassInBin))
            self.GrSph.Lx.append(Lx)
            self.GrSph.Ly.append(Ly)
            self.GrSph.Lz.append(Lz)

        self.GrSph.R = (scipy.array(self.GrSph.Rmin)+scipy.array(self.GrSph.Rmax))/2.0
        gamma=[]
        kappa=[]
        gamma.append(0.5*(self.GrSph.R[0]/self.GrSph.Rho[0]+self.GrSph.R[1]/self.GrSph.Rho[1])*(self.GrSph.Rho[1]-self.GrSph.Rho[0])/(self.GrSph.R[1]-self.GrSph.R[0]))
        kappa.append(0.5*(self.GrSph.R[0]/self.GrSph.Sigma2r[0]+self.GrSph.R[1]/self.GrSph.Sigma2r[1])*(self.GrSph.Sigma2r[1]-self.GrSph.Sigma2r[0])/(self.GrSph.R[1]-self.GrSph.R[0]))
        for i in range(1,NBins-1):
            gamma.append(self.GrSph.R[i]/self.GrSph.Rho[i]*(self.GrSph.Rho[i+1]-self.GrSph.Rho[i-1])/(self.GrSph.R[i+1]-self.GrSph.R[i-1]))
            kappa.append(self.GrSph.R[i]/self.GrSph.Sigma2r[i]*(self.GrSph.Sigma2r[i+1]-self.GrSph.Sigma2r[i-1])/(self.GrSph.R[i+1]-self.GrSph.R[i-1]))

        gamma.append(0.5*(self.GrSph.R[-2]/self.GrSph.Rho[-2]+self.GrSph.R[-1]/self.GrSph.Rho[-1])*(self.GrSph.Rho[-1]-self.GrSph.Rho[-2])/(self.GrSph.R[-1]-self.GrSph.R[-2]))
        kappa.append(0.5*(self.GrSph.R[-2]/self.GrSph.Sigma2r[-2]+self.GrSph.R[-1]/self.GrSph.Sigma2r[-1])*(self.GrSph.Sigma2r[-1]-self.GrSph.Sigma2r[-2])/(self.GrSph.R[-1]-self.GrSph.R[-2]))
        
        self.GrSph.Gamma=gamma
        self.GrSph.Kappa=kappa
        
        #convert to scipy arrays:
        self.GrSph.Rho = scipy.array(self.GrSph.Rho)
        self.GrSph.MassInBin = scipy.array(self.GrSph.MassInBin)
        self.GrSph.N = scipy.array(self.GrSph.N)
        self.GrSph.Rmin = scipy.array(self.GrSph.Rmin)
        self.GrSph.Rmax = scipy.array(self.GrSph.Rmax)
        self.GrSph.R = scipy.array(self.GrSph.R)
        if self.Snapshot.V != None:
            self.GrSph.V = scipy.array(self.GrSph.V)
        self.GrSph.Sigma2 = scipy.array(self.GrSph.Sigma2)
        self.GrSph.Sigma2r = scipy.array(self.GrSph.Sigma2r)
        self.GrSph.Beta = scipy.array(self.GrSph.Beta)
        self.GrSph.MeanVr = scipy.array(self.GrSph.MeanVr)
        self.GrSph.Gamma = scipy.array(self.GrSph.Gamma)
        self.GrSph.Kappa = scipy.array(self.GrSph.Kappa)
        self.GrSph.CumulativeMass = scipy.array(self.GrSph.CumulativeMass)
        self.GrSph.JeansMassFrac = - self.GrSph.CumulativeMass / self.GrSph.R / self.GrSph.Sigma2r / (self.GrSph.Gamma+self.GrSph.Kappa+2.0*self.GrSph.Beta)
        self.GrSph.Lx = scipy.array(self.GrSph.Lx)
        self.GrSph.Ly = scipy.array(self.GrSph.Ly)
        self.GrSph.Lz = scipy.array(self.GrSph.Lz)
        
        if CalcUncBeta and UseSphericalCoordinates:
            self.GrSph.UncBetaBootstrap = scipy.array(self.GrSph.UncBetaBootstrap)
            self.GrSph.MeanBetaBootstrap = scipy.array(self.GrSph.MeanBetaBootstrap)
        
        print "CreateGridLogBins ended\n"
        
        return self.GrSph
        

    def CreateGrid(self,ParticlesPerBin=30000,Method='Radial',SubtractRadialVelocity = False):
        print "\nCreateGrid started"
        "Creates grid with same number of particles in each bin."
        r = scipy.sqrt(self.Snapshot.x*self.Snapshot.x+self.Snapshot.y*self.Snapshot.y+self.Snapshot.z*self.Snapshot.z)
        v2 = self.Snapshot.vx*self.Snapshot.vx+self.Snapshot.vy*self.Snapshot.vy+self.Snapshot.vz*self.Snapshot.vz
        vr = 1.0/(r+1.0e-16)*scipy.array(self.Snapshot.vx*self.Snapshot.x + self.Snapshot.vy*self.Snapshot.y + self.Snapshot.vz*self.Snapshot.z)

        if Method == 'Radial':
            index = scipy.argsort(r)
        elif Method == 'Potential':
            index = scipy.argsort(self.Snapshot.V)
        
        
        N = len(index)
        BinNo = 0
        self.Gr = Grid()
        while (BinNo+1)*ParticlesPerBin < N:
            Particles = index[ range(BinNo*ParticlesPerBin,(BinNo+1)*ParticlesPerBin) ]
            self.Gr.R.append(scipy.mean(r[Particles]))

            Rmin = r[Particles].min()
            Rmax = r[Particles].max()
            Volume = 4.0/3.0 * 3.14159265359 * (Rmax**3 - Rmin**3)
            Mass = scipy.sum(self.Snapshot.m[Particles])
            self.Gr.Rmin.append(Rmin)
            self.Gr.Rmax.append(Rmax)
            self.Gr.N.append(ParticlesPerBin) 
            self.Gr.MassInBin.append(Mass) 
            self.Gr.Rho.append(Mass / Volume)
            Sigma2 = scipy.mean(v2[Particles])
            MeanVr = scipy.mean(vr[Particles])
            V = scipy.mean(self.Snapshot.V[Particles])
            if SubtractRadialVelocity == True:
                Sigma2r = scipy.mean((vr[Particles]-MeanVr)**2)
            else:
                Sigma2r = scipy.mean(vr[Particles]**2)

            Beta = 1.0 - (Sigma2-Sigma2r)/(2.0*Sigma2r)

            self.Gr.Sigma2.append(Sigma2 )            
            self.Gr.Sigma2r.append(Sigma2r )
            self.Gr.V.append(V)
            self.Gr.Beta.append(Beta)
            self.Gr.MeanVr.append(MeanVr)
            self.Gr.CumulativeMass.append(sum(self.Gr.MassInBin))
            BinNo += 1
        
        
        #derivatives - gamma and kappa:
        gamma=[]
        kappa=[]
        gamma.append(0.5*(self.Gr.R[0]/self.Gr.Rho[0]+self.Gr.R[1]/self.Gr.Rho[1])*(self.Gr.Rho[1]-self.Gr.Rho[0])/(self.Gr.R[1]-self.Gr.R[0]))
        kappa.append(0.5*(self.Gr.R[0]/self.Gr.Sigma2r[0]+self.Gr.R[1]/self.Gr.Sigma2r[1])*(self.Gr.Sigma2r[1]-self.Gr.Sigma2r[0])/(self.Gr.R[1]-self.Gr.R[0]))
        for i in range(1,len(self.Gr.R)-1):
            gamma.append(self.Gr.R[i]/self.Gr.Rho[i]*(self.Gr.Rho[i+1]-self.Gr.Rho[i-1])/(self.Gr.R[i+1]-self.Gr.R[i-1]))
            kappa.append(self.Gr.R[i]/self.Gr.Sigma2r[i]*(self.Gr.Sigma2r[i+1]-self.Gr.Sigma2r[i-1])/(self.Gr.R[i+1]-self.Gr.R[i-1]))

        gamma.append(0.5*(self.Gr.R[-2]/self.Gr.Rho[-2]+self.Gr.R[-1]/self.Gr.Rho[-1])*(self.Gr.Rho[-1]-self.Gr.Rho[-2])/(self.Gr.R[-1]-self.Gr.R[-2]))
        kappa.append(0.5*(self.Gr.R[-2]/self.Gr.Sigma2r[-2]+self.Gr.R[-1]/self.Gr.Sigma2r[-1])*(self.Gr.Sigma2r[-1]-self.Gr.Sigma2r[-2])/(self.Gr.R[-1]-self.Gr.R[-2]))
        
        self.Gr.Gamma=gamma
        self.Gr.Kappa=kappa
        
        #convert to scipy arrays:
        self.Gr.Rho = scipy.array(self.Gr.Rho)
        self.Gr.MassInBin = scipy.array(self.Gr.MassInBin)
        self.Gr.N = scipy.array(self.Gr.N)
        self.Gr.Rmin = scipy.array(self.Gr.Rmin)
        self.Gr.Rmax = scipy.array(self.Gr.Rmax)
        self.Gr.R = scipy.array(self.Gr.R)
        self.Gr.V = scipy.array(self.Gr.V)
        self.Gr.Sigma2 = scipy.array(self.Gr.Sigma2)
        self.Gr.Sigma2r = scipy.array(self.Gr.Sigma2r)
        self.Gr.Beta = scipy.array(self.Gr.Beta)
        self.Gr.MeanVr = scipy.array(self.Gr.MeanVr)
        self.Gr.Gamma = scipy.array(self.Gr.Gamma)
        self.Gr.Kappa = scipy.array(self.Gr.Kappa)
        self.Gr.CumulativeMass = scipy.array(self.Gr.CumulativeMass)
        print "CreateGrid finished\n"
        
        return self.Gr
    
    
    def CreateEnergyGrid(self,ParticlesPerBin=1000):
        v2 = self.Snapshot.vx*self.Snapshot.vx+self.Snapshot.vy*self.Snapshot.vy+self.Snapshot.vz*self.Snapshot.vz
        E = 0.5*v2 + self.Snapshot.V
        
        index = scipy.argsort(E)
        
        tmpE = []
        tmpMass = []
        
        N = len(index)
        BinNo = 0
        self.EGrid = EnergyGrid()
        
        TotalMass = self.Snapshot.m.sum()

        while (BinNo+1)*ParticlesPerBin < N:
            Particles = index[ range(BinNo*ParticlesPerBin,(BinNo+1)*ParticlesPerBin) ]
            Max = E[Particles].max()
            Min = E[Particles].min()
            Mean = E[Particles].mean()
            tmpE.append( Mean )
            tmpMass.append( self.Snapshot.m[Particles].sum() / ( Max - Min )  )
            BinNo += 1
        
        self.EGrid.Mass = scipy.array(tmpMass)
        self.EGrid.E = scipy.array(tmpE)        
        return self.EGrid


    def CreateEnergyGridLinear(self,NBins = 100):
        v2 = self.Snapshot.vx*self.Snapshot.vx+self.Snapshot.vy*self.Snapshot.vy+self.Snapshot.vz*self.Snapshot.vz
        E = 0.5*v2 + self.Snapshot.V
        NPartTotal = len(E)
        
 
        Emin = E.min()
        Emax = E.max()

        LeftBinLimit = scipy.linspace(Emin,Emax,NBins)

        dE = LeftBinLimit[1] - LeftBinLimit[0]
        BinNo = scipy.array((E-LeftBinLimit[0])/dE,int)

        self.EGridLin = EnergyGrid()

        tmpE = []
        tmpMass = []

        for Bin in range(NBins):
            Particles = scipy.where(BinNo==Bin)[0]
            Max = E[Particles].max()
            Min = E[Particles].min()

            tmpE.append( LeftBinLimit[Bin]+dE/2.0 )
            tmpMass.append( len(Particles) / ( Max - Min ) / NPartTotal )

        self.EGridLin.Mass = scipy.array(tmpMass)
        self.EGridLin.E = scipy.array(tmpE)        
        return self.EGridLin

    def GetEGrid(self):
        return self.EGrid
        
    def GetEGridLin(self):
        return self.EGridLin
    
    def GetGrid(self):
        return self.Gr

    def GetGridSph(self):
        return self.GrSph

    def GetLowestPotential(self):
        return self.Snapshot.V.min()
