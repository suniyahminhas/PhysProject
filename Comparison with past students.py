from __future__ import division
import numpy
import matplotlib.pyplot as pyplot
from scipy.integrate import quad
import RadialFunctions as k

#this function is an attempt to use degenerate perturbation theory

a_0 = 1 #This can be changed to the genuine bohr radius

print(k.R_10(2))

#we now define different values for the mass of the boson
m_b_values= numpy.arange(10**-4,2,10**-4) #in mass of electron
data_matrix = numpy.zeros((len(m_b_values),1))

pyplot.figure()

def pert_for_function(Radial_function):
    #this is what would be in the integral (i,e functions to integrate)
    def integ_function(r):
        return (numpy.exp(-m_boson*r)/((4*numpy.pi*r))*(Radial_function(r))**2)*r**2
   
    #each mass then works out the integrals 
    for i in range(len(m_b_values)):
        m_boson = m_b_values[i]
        data_matrix[i] =quad(integ_function,0, numpy.inf)[0] #get all answers in hartree units

    #now plot the two energy perturbations against the mass
    pyplot.plot(m_b_values*0.510998,data_matrix*27.2114, label=str(Radial_function.__name__), color=C,linestyle=L)
    #the multiplications are the scalings
    pyplot.legend()

pyplot.xlabel('m_boson/Mev')
pyplot.xscale('log')
pyplot.xlim(10**-4,1)
pyplot.ylabel('energy pertubation/ev')
pyplot.title('Energy perturbation with different boson masses Z=1 logspace')
pyplot.show()

C='green'
L='-'
pert_for_function(k.R_20)
C='blue'
pert_for_function(k.R_10)