from __future__ import division
import numpy
import matplotlib.pyplot as pyplot
from scipy.integrate import quad
import RadialFunctions as f


#we now define different values for the Radial values
R_p_values= numpy.arange(0.1,6,0.1)
data_matrix = numpy.zeros((len(R_p_values),1))

pyplot.figure()

def pert_for_function(Radial_function):
    #this is what would be in the integral (i,e functions to integrate)
    def integ_function(r):
        return (numpy.exp(-r/R_phi)/((4*numpy.pi*r))*(Radial_function(r))**2)*r**2
   
    #each mass then works out the integrals 
    for i in range(len(R_p_values)):
        R_phi = R_p_values[i]
        data_matrix[i] =quad(integ_function,0, numpy.inf)[0]

    #now plot the two energy perturbations against the mass
    pyplot.plot(R_p_values,data_matrix, label=str(Radial_function.__name__), color=C,linestyle=L)
    pyplot.legend()

pyplot.xlabel('R_phi/a_0')
pyplot.ylabel('energy pertubation/(E_hartree*CouplingConstant)')
pyplot.ylim(0,0.008)
pyplot.title('Energy perturbation with different boson masses')
pyplot.show()

C='blue'
L='-'
pert_for_function(f.R_20)
L='--'
pert_for_function(f.R_21)


C='red'
L='-'
pert_for_function(f.R_30)
L='--'
pert_for_function(f.R_31)
L=':'
pert_for_function(f.R_32)

C='green'
L='-'
pert_for_function(f.R_40)
L='--'
pert_for_function(f.R_41)
L=':'
pert_for_function(f.R_42)
L= '-.'
pert_for_function(f.R_43)