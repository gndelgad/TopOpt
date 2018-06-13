# -*- coding: utf-8 -*-

from matplotlib import pylab as pl
#import pylab as pl
import numpy as np
import matplotlib.pyplot as plt
import os

Niter=100
it=range(0,Niter-1)

#lire le fichier sortie

Lagrangian=[]
fidJ = open('Lagrangian.txt','r')

fich_rubJ = fidJ.read().split('\n')
for iter in range(0,Niter-1):
	w=1.*float(fich_rubJ[iter])
	Lagrangian = np.append(Lagrangian,w)
fidJ.close()

pl.figure()
pl.plot(it,Lagrangian,linewidth=4)
pl.title('Convergence topology optimization',fontsize=20)
pl.ylabel('Lagrangian',fontsize=20)
pl.xlabel('Iterations',fontsize=20)
pl.xticks(fontsize=10)
pl.yticks(fontsize=10)
#pl.title('Composite optimization for Lambda1')
pl.show()

