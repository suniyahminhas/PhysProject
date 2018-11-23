from __future__ import division
import numpy
import matplotlib.pyplot as pyplot
from sympy import var
from sympy import Symbol, integrate, pi, oo, sin, exp
from sympy.physics.hydrogen import R_nl
import numpy

r=Symbol("r", real=True, positive=True)
m_b_values= numpy.array([5,0.5,0.05, 0.005, 0.0005, 0.00005])

pyplot.figure()
def pert_for_function(N,L,Z):
    #this is what would be in the integral (i,e functions to integrate)
    data = integrate((exp(-m_boson*r)/((4*pi*r))*(R_nl(N,L,r,Z))**2)*r**2, (r,0,oo))
    return float(data)

for i in range(len(m_b_values)):
    m_boson=m_b_values[i]/0.510998 # in Mev but this shows in electron mass
    y=numpy.zeros((len(range(12)),1))
    for l in range(12):
        y[l]=pert_for_function(12,l,1)
    pyplot.plot(range(12), y, label=str(m_b_values[i]), marker='.')
    
pyplot.xlabel('l')
pyplot.xlim(0,11)
pyplot.yscale('log')
pyplot.ylim(10**-8, 10**-3)
pyplot.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Boson mass in Mev')
pyplot.ylabel('energy pertubation/(E_hartree*CouplingConstant)')
pyplot.title('Energy perturbation for different m_boson values')
pyplot.show()