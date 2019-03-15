from __future__ import division
from sympy import var
from sympy import Symbol, integrate, pi, oo, sin, exp
from sympy.physics.hydrogen import R_nl
import numpy
import matplotlib.pyplot as pyplot

#this function is an attempt to use degenerate perturbation theory

a_0 = 1 #This can be changed to the genuine bohr radius

#we now define different values for the mass of the boson
m_b_values= numpy.logspace(-3,3)
data_matrix = numpy.zeros((len(m_b_values),1))
r=Symbol("r", real=True, positive=True)

pyplot.figure()

def pert_for_function(N,L,Z):
    #this is what would be in the integral (i,e functions to integrate)
    for i in range(len(m_b_values)):
        m_boson = m_b_values[i]
        data = integrate((exp(-m_boson*r)/((4*pi*r))*(R_nl(N,L,r,Z))**2)*r**2, (r,0,oo))
        data_matrix[i] = float(data)
    pyplot.subplot(1,2,1)
    pyplot.plot(m_b_values,data_matrix, color=C, linestyle=K, label=str(N)+ '    ' + str(L), linewidth=3.0)
    pyplot.xlabel('$m_{\phi} (1/a_0)$', fontsize=30)
    pyplot.ylabel('$dE_{nl} (E_{\mathrm{h}})$', fontsize=30)
    pyplot.yscale('log')
    pyplot.xscale('log')
    pyplot.ylim(10**-6,10**-1)
    pyplot.xticks([10**-3, 10**-1, 10**1, 10**3], fontsize=25)
    pyplot.yticks(fontsize=25)
    pyplot.subplot(1,2,2)
    pyplot.plot(m_b_values*(N**2),data_matrix*(N**2), color=C, linestyle=K, label=str(N)+ '    ' + str(L), linewidth=3.0)
    pyplot.xlabel('$m_{\phi} \cdot n^2 (1/a_0)$', fontsize=30)
    pyplot.ylabel('$dE_{nl} \cdot n^2 (E_{\mathrm{h}})$', fontsize=30)
    pyplot.yscale('log')
    pyplot.xscale('log')
    pyplot.xlim(10**-3, 10**4)
    pyplot.ylim(10**-6,10**-1)
    pyplot.xticks([10**-3, 10**-1, 10**1, 10**3], fontsize=25)
    pyplot.yticks(fontsize=25)
    leg=pyplot.legend(title='        n    l',loc='center left', bbox_to_anchor=(1, 0.5) ,fontsize=25)
    pyplot.setp(leg.get_title(),fontsize='20')
    pyplot.subplots_adjust(hspace=0.3, wspace=0.3)

    
pyplot.xlabel('$m_{\phi} (1/a_0)$', fontsize=30)
pyplot.ylabel('$dE_{nl} (E_{\mathrm{h}})$', fontsize=30)


pyplot.yscale('log')
pyplot.xscale('log')
pyplot.ylim(10**-6,10**-1)
pyplot.xticks([10**-3, 10**-1, 10**1, 10**3], fontsize=25)
pyplot.yticks(fontsize=25)



pyplot.show()

C='purple'
K='-'
pert_for_function(1,0,1)

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