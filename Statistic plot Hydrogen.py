from __future__ import division
from sympy import var
from sympy import Symbol, integrate, pi, oo, sin, exp
from sympy.physics.hydrogen import R_nl
import numpy
import matplotlib.pyplot as pyplot
import matplotlib.axis as ax
from scipy.misc import imread
from PIL import Image, ImageDraw


m_b_values = numpy.array([10, 10**2, 10**3, 10**4, 10**5, 10**6])
constraint = numpy.array([10**-11, 10**-9, 10**-7, 10**-5, 10**-3, 10**-1])

pyplot.figure()

img = imread("pythonreadyfade.png")
pyplot.imshow(img,zorder=0, extent= (10**-6, 10**3, 10**-16, 1), aspect='auto')
pyplot.show()

pyplot.plot(m_b_values*(1/(2.6817268*(10**-4)))*(10**-6),constraint, color='limegreen', linewidth=5.0) 
pyplot.xlabel('$m_\phi (MeV)$', fontsize=40)
pyplot.ylabel('$g_{e}g_{n}$', fontsize=40)
pyplot.rc('xtick', labelsize=30) 
pyplot.rc('ytick', labelsize=30) 
pyplot.yscale('log')
pyplot.ylim(10**-16, 1)
pyplot.xlim(10**-6, 10**3)
#matplotlib.pyplot.xticks(ticks=None)
pyplot.xscale('log')
pyplot.xticks([10**-6,10**-3,1,10**3])
pyplot.yticks([10**0,10**-4,10**-8, 10**-12, 10**-16])

pyplot.show()

