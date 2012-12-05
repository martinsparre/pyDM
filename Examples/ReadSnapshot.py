#/usr/bin/python
import matplotlib.pyplot as plt
import Classes.DM_structure as DM_structure
import scipy


FileName = 'Hernquist100000_000'
A = DM_structure.DM_structure(FileName)
S = A.Snapshot
