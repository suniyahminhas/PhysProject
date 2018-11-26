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

n_values = numpy.arange(2,13)
n_squared= 1/((n_values)**2)

#data taken from Table 5 page 628

# This is data, these are the transiton frequencies given in MHz 
#The transitions are from 1S(1/2) to nP(1/2), with n given by the INDEX
t=numpy.zeros((13,1))
t[2]= 246606800040000
t[3]= 292272900090000 #had to give t[2] and t[3] extra zeros?? Very confused
t[4]= 308257000060000
t[5]= 315656700013000
t[6]= 319676000070000
t[7]= 322098000070000
t[8]= 323670000070000
t[9]= 324749000070000
t[10]= 325518000070000
t[11]= 326092000070000
t[12]= 326525000070000

t=t*10
#NOTE this hasn't been scaled to Hz so works better this way ??
#This is the only scaling that appears correct


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
    return data_matrix

#we now define the experimental values and their errors
#all need to be looked at all only approximate for now

h= 6.626070040*10**-34 #plancks constant, 
k= 27.21138602 #ev in one Hartree
e= 1.6021766208*10**-19 #electronic charge

def data_vs_theory(data):
    E=numpy.zeros((13,3))
    for n in numpy.arange(2,13):
        E[n,0]= h*data[n]/(e*k)
        E[n,1] = E_nl(n,1)-E_nl(1,1)
        E[n,2]= E[n,1]-E[n,0] #Theory mnus experiment
    return E

print data_vs_theory(t)[2:,0]
print data_vs_theory(t)[2:,1]


pyplot.figure()
pyplot.plot(n_squared, data_vs_theory(t)[2:,0], label='exp')
pyplot.plot(n_squared, data_vs_theory(t)[2:,1], label='theory')
pyplot.legend(loc='center left', bbox_to_anchor=(1, 0.5))
pyplot.xlabel('1/n**2')
pyplot.ylabel('Energy/hartree')
pyplot.title('Energy jumps experimental vs theoretical')
pyplot.show()

pyplot.figure()
pyplot.plot(n_squared, data_vs_theory(t)[2:,2], marker='o')
pyplot.legend(loc='center left', bbox_to_anchor=(1, 0.5))
pyplot.xlabel('1/n**2')
pyplot.ylabel('Energy/hartree')
pyplot.title('Energy jumps experimental vs theoretical difference')
pyplot.show()