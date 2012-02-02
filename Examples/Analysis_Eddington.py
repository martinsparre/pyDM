#/usr/bin/python
import matplotlib.pyplot as plt
import ReadGadget2
import DM_structure
from scipy import log10,sqrt,arccos,arcsin
import scipy
import scipy.optimize, scipy.interpolate

from math import pi
def f(E):
    return 1.0/sqrt(2.0) / (2.0*3.14159265359)**3. * sqrt(E) / (1.0-E)**2. * (  (1.-2.*E)* (8.*E**2. - 8.*E    -3.0 ) + 3.*arcsin(sqrt(E))/sqrt(E*(1.-E))      )

def g(E):
    A = 1.0/ sqrt(E*E)
    return (4.*pi)**2*sqrt(2.0*sqrt(E*E)) *(sqrt(A-1.) * (1./8.0 *A**2 - 5.0/12.0 *A - 1.0/3.0   )+ 1.0/8.0 *A * (A**2.-4.*A+8.)*arccos(1./sqrt(A)))

def N(E):
    return f(E)*g(E)

FileName='ISOHQ_001_LOWSOFT_000' #'OM_ROI00_rAN1_HQ_000'
A = DM_structure.DM_structure(FileName)


Grid = A.CreateGrid(Method='Radial')
Grid.CalculateDF()

GridSph = A.CreateGridLogBins()



#x=scipy.logspace(-2.0,2.0,100)
#plt.plot(log10(x),-1.0/(1.0+x),'-',color='black')
#plt.plot(log10(Grid.R),Grid.V,'<-',color='red',alpha=0.5)
#plt.plot(log10(GridSph.R),GridSph.V,'<-',color='blue',alpha=0.5)
#plt.show()

