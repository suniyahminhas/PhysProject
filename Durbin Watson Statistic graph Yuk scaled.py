from __future__ import division
from sympy import var
from sympy import Symbol, integrate, pi, oo, sin, exp
from sympy.physics.hydrogen import R_nl
import numpy
import matplotlib.pyplot as pyplot
from scipy.misc import imread
from PIL import Image, ImageDraw


m_b_values = numpy.array([10, 1*10**2, 1*10**3, 1*10**4, 1*10**5, 1*10**6])
constraint = numpy.array([1*10**-14,1*10**-12, 1*10**-10, 1.38949549437e-08, 1.93069772888e-06, 0.000117876863479])

#10       1*10**-14
#1*10**2    1*10**-12
#1*10**3    1*10**-10
#1*10**4    1.38949549437e-08
#1*10**5    1.93069772888e-06
#1*10**6    0.000117876863479




m_b_values_H = numpy.array([10, 10**2, 10**3, 10**4, 10**5, 10**6])
constraint_H = numpy.array([10**-11, 10**-9, 10**-7, 10**-5, 10**-3, 0.1]) #Yukawa potential unscaled

pyplot.figure()


img = imread("pythonreadyfade.png") 
pyplot.imshow(img,zorder=0, extent= (10**-6, 10**3, 10**-16, 1), aspect='auto')

pyplot.plot(m_b_values*(1/(2.6817268*(10**-4)))*(10**-6),constraint, color='red', linewidth=4.0, label='1000') 
pyplot.plot(m_b_values_H*(1/(2.6817268*(10**-4)))*(10**-6), constraint_H, color='darkgreen',linewidth=4.0, label='1')
pyplot.legend(loc= 2, fontsize=30)
pyplot.xlabel('$m_\phi (MeV)$', fontsize=40)
pyplot.ylabel('$g_{e}g_{n}$', fontsize=40)
pyplot.rc('xtick', labelsize=30) 
pyplot.rc('ytick', labelsize=30) 
pyplot.yscale('log')
pyplot.ylim(10**-16, 1)
pyplot.xlim(10**-6, 10**3)
pyplot.xscale('log')
pyplot.xticks([10**-6,10**-3,1,10**3])
pyplot.yticks([10**0,10**-4,10**-8, 10**-12, 10**-16])

pyplot.show()
