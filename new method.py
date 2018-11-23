from __future__ import division
from sympy import var
from sympy import Symbol, integrate, pi, oo, sin, exp
from sympy.physics.hydrogen import R_nl
import numpy
import matplotlib.pyplot as pyplot

r=Symbol("r", real=True, positive=True)
m_b_values= numpy.arange(0,6,0.1)
data_matrix = numpy.zeros((len(m_b_values),1))
pyplot.figure()

def pert_for_function(N,L,Z):
    #this is what would be in the integral (i,e functions to integrate)
    for i in range(len(m_b_values)):
        m_boson = m_b_values[i]
        data = integrate((exp(-m_boson*r)/((4*pi*r))*(R_nl(N,L,r,Z))**2)*r**2, (r,0,oo))
        data_matrix[i] = float(data)
    pyplot.plot(m_b_values,data_matrix, color=C, linestyle=K, label='R_'+str(N)+str(L))
    pyplot.legend()

pyplot.xlabel('m_boson/electron_mass')
pyplot.ylabel('energy pertubation/(E_hartree*CouplingConstant)')
pyplot.ylim(0,0.008)
pyplot.title('Energy perturbation with different boson masses Z=1')
pyplot.show()

C='blue'
K='-'
pert_for_function(2,0,1)
K='--'
pert_for_function(2,1,1)


C='red'
K='-'
pert_for_function(3,0,1)
K='--'
pert_for_function(3,1,1)
K=':'
pert_for_function(3,2,1)

C='green'
K='-'
pert_for_function(4,0,1)
K='--'
pert_for_function(4,1,1)
K=':'
pert_for_function(4,2,1)
K= '-.'
pert_for_function(4,3,1)