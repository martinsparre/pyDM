#/usr/bin/python
import matplotlib.pyplot as plt
import Classes.DM_structure as DM_structure
import scipy


from mayavi import mlab

FileName='Hq0.001.bin'
A = DM_structure.DM_structure(FileName, ImportPotential=False)

