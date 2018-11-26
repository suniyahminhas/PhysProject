from __future__ import division
import numpy
import matplotlib.pyplot as pyplot
from scipy.integrate import quad
import RadialFunctions as f
from sympy import var
from sympy.physics.hydrogen import E_nl
from sympy import Symbol, integrate, pi, oo, sin, exp
from sympy.physics.hydrogen import R_nl

#The aim of this program is to calculate the pertubation differences for each energy level due to the Yukawa potential
#Then to calculate the differences in energy jumps (as seen in atomic spectroscopy)
#Then to compare this value to experimental value to see if there is room for this pertuba
#This will all be evaluated at different m values and for a coupling constant of 1 
r=Symbol("r", real=True, positive=True)
m_b_values= numpy.arange(0,6,0.1)
data_matrix = numpy.zeros((len(m_b_values),1))
pyplot.figure()
#we now define different values for the mass of the boson
def pert_for_function(N,L,Z):
    #this is what would be in the integral (i,e functions to integrate)
    for i in range(len(m_b_values)):
        m_boson = m_b_values[i]
        data = integrate((exp(-m_boson*r)/((4*pi*r))*(R_nl(N,L,r,Z))**2)*r**2, (r,0,oo))
        data_matrix[i] = float(data)
    return data_matrix

#we now define the experimental values and their errors
#all need to be looked at all only approximate for now

h= 6.63*10**-34 #plancks constant, 
k=27.2114 #how many electron volts in one hartree-CHECK ME
e=1.6*10**-19 #electronic charge

def exp_vs_theory(f_exp, grid_1, grid_2):
    Theory = E_nl(grid_1[0], grid_1[2])- E_nl(grid_2[0], grid_2[2]) #The energy level expected by theory
    Enl_exp_h= h*f/(e*k) #the experimental value in terms of hartree for the energy jump En+1-En
    pert_observed = Enl_exp_h - Theory
    pyplot.axhline(y=pert_observed, label='pert observed')
    pert_theory= numpy.subtract(pert_for_function(grid_1[0], grid_1[1], grid_1[2]), pert_for_function(grid_2[0], grid_2[1], grid_2[2]))
    print(pert_theory)
    pyplot.plot(m_b_values, pert_theory, label='pert theory')


f=2466.061*10**12 #this is an example and is simply for the 1s to 2s transition
exp_vs_theory(f,numpy.array([2,0,1]),numpy.array([1,0,1]))

pyplot.legend()
pyplot.xlabel('m_boson/electron_mass')
pyplot.ylabel('energy pertubation/(E_hartree)')
pyplot.title('Energy perturbation with different boson masses, CC=10**-2')
pyplot.show()
