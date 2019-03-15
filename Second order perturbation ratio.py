#Function to calculate the second order peturbation from Yukawa potential

from __future__ import division
import numpy
import matplotlib.pyplot as pyplot
from scipy.integrate import quad
import RadialFunctions as f
from sympy import var 
from sympy.physics.hydrogen import E_nl
from sympy import var
from sympy import Symbol, integrate, pi, oo, sin, exp
from sympy.physics.hydrogen import R_nl

pyplot.figure()


m_b_values= numpy.logspace(-3,3)
data_matrix = numpy.zeros((len(m_b_values),1))
r=Symbol("r", real=True, positive=True)

def pert_for_function(N,L,Z):
    #this is what would be in the integral (i,e functions to integrate)
    for i in range(len(m_b_values)):
        m_boson = m_b_values[i]
        data = integrate((exp(-m_boson*r)/((4*pi*r))*(R_nl(N,L,r,Z))**2)*r**2, (r,0,oo))
        data_matrix[i] = float(data)
    return data_matrix


#this only works for R_10
def second_order_10(Radial_function):
    m_b_values= numpy.logspace(-3,3)
    other= [f.R_20,f.R_30, f.R_40, f.R_50, f.R_60, f.R_70]
    Energy=[E_nl(2),E_nl(3),E_nl(4), E_nl(5), E_nl(6), E_nl(7)]
    Base_Energy= E_nl(1)
    data_matrix = numpy.zeros((len(m_b_values),len(other)))
    y_matrix = numpy.zeros((len(m_b_values),1))
    #this is what would be in the integral (i,e functions to integrate)
    def integ_function(r):
        return (numpy.exp(-m_boson*r)/((4*numpy.pi*r))*(Radial_function(r))*Radial_function_2(r))*r**2
   
    #each mass then works out the integrals 
    for F in range(len(other)):
        Radial_function_2= other[F]
        H=Energy[F]-Base_Energy
        for i in range(len(m_b_values)):
            m_boson = m_b_values[i]
            K=(quad(integ_function,0, numpy.inf)[0])**2
            y_matrix[i] +=K/H
    ratio=numpy.divide(y_matrix, pert_for_function(1,0,1))
    return ratio 


def second_order_20(Radial_function):
    m_b_values= numpy.logspace(-3,3)
    other= [f.R_30, f.R_40, f.R_50, f.R_60, f.R_70]
    Energy=[E_nl(3),E_nl(4), E_nl(5), E_nl(6), E_nl(7)]
    Base_Energy= E_nl(2)
    data_matrix = numpy.zeros((len(m_b_values),len(other)))
    y_matrix = numpy.zeros((len(m_b_values),1))
    #this is what would be in the integral (i,e functions to integrate)
    def integ_function(r):
        return (numpy.exp(-m_boson*r)/((4*numpy.pi*r))*(Radial_function(r))*Radial_function_2(r))*r**2
   
    #each mass then works out the integrals 
    for F in range(len(other)):
        Radial_function_2= other[F]
        H=Energy[F]-Base_Energy
        for i in range(len(m_b_values)):
            m_boson = m_b_values[i]
            K=(quad(integ_function,0, numpy.inf)[0])**2
            y_matrix[i] +=K/H
    ratio=numpy.divide(y_matrix, pert_for_function(2,0,1))
    return ratio  
    
def second_order_21(Radial_function):
    m_b_values= numpy.logspace(-3,3)
    other= [f.R_31, f.R_41, f.R_51, f.R_61, f.R_71]
    Energy=[E_nl(3),E_nl(4), E_nl(5), E_nl(6), E_nl(7)]
    Base_Energy= E_nl(2)
    data_matrix = numpy.zeros((len(m_b_values),len(other)))
    y_matrix = numpy.zeros((len(m_b_values),1))
    #this is what would be in the integral (i,e functions to integrate)
    def integ_function(r):
        return (numpy.exp(-m_boson*r)/((4*numpy.pi*r))*(Radial_function(r))*Radial_function_2(r))*r**2
   
    #each mass then works out the integrals 
    for F in range(len(other)):
        Radial_function_2= other[F]
        H=Energy[F]-Base_Energy
        for i in range(len(m_b_values)):
            m_boson = m_b_values[i]
            K=(quad(integ_function,0, numpy.inf)[0])**2
            y_matrix[i] +=K/H
    ratio=numpy.divide(y_matrix, pert_for_function(2,1,1))
    return ratio
        
def second_order_30(Radial_function):
    m_b_values= numpy.logspace(-3,3)
    other= [f.R_40, f.R_50, f.R_60, f.R_70]
    Energy=[E_nl(4), E_nl(5), E_nl(6), E_nl(7)]
    Base_Energy= E_nl(3)
    data_matrix = numpy.zeros((len(m_b_values),len(other)))
    y_matrix = numpy.zeros((len(m_b_values),1))
    #this is what would be in the integral (i,e functions to integrate)
    def integ_function(r):
        return (numpy.exp(-m_boson*r)/((4*numpy.pi*r))*(Radial_function(r))*Radial_function_2(r))*r**2
   
    #each mass then works out the integrals 
    for F in range(len(other)):
        Radial_function_2= other[F]
        H=Energy[F]-Base_Energy
        for i in range(len(m_b_values)):
            m_boson = m_b_values[i]
            K=(quad(integ_function,0, numpy.inf)[0])**2
            y_matrix[i] +=K/H
    ratio=numpy.divide(y_matrix, pert_for_function(3,0,1))
    return ratio    
    
def second_order_31(Radial_function):
    m_b_values= numpy.logspace(-3,3)
    other= [f.R_41, f.R_51, f.R_61, f.R_71]
    Energy=[E_nl(4), E_nl(5), E_nl(6), E_nl(7)]
    Base_Energy= E_nl(3)
    data_matrix = numpy.zeros((len(m_b_values),len(other)))
    y_matrix = numpy.zeros((len(m_b_values),1))
    #this is what would be in the integral (i,e functions to integrate)
    def integ_function(r):
        return (numpy.exp(-m_boson*r)/((4*numpy.pi*r))*(Radial_function(r))*Radial_function_2(r))*r**2
   
    #each mass then works out the integrals 
    for F in range(len(other)):
        Radial_function_2= other[F]
        H=Energy[F]-Base_Energy
        for i in range(len(m_b_values)):
            m_boson = m_b_values[i]
            K=(quad(integ_function,0, numpy.inf)[0])**2
            y_matrix[i] +=K/H
    ratio=numpy.divide(y_matrix, pert_for_function(3,1,1))
    return ratio
    
def second_order_32(Radial_function):
    m_b_values= numpy.logspace(-3,3)
    other= [f.R_42, f.R_52, f.R_62, f.R_72]
    Energy=[E_nl(4), E_nl(5), E_nl(6), E_nl(7)]
    Base_Energy= E_nl(3)
    data_matrix = numpy.zeros((len(m_b_values),len(other)))
    y_matrix = numpy.zeros((len(m_b_values),1))
    #this is what would be in the integral (i,e functions to integrate)
    def integ_function(r):
        return (numpy.exp(-m_boson*r)/((4*numpy.pi*r))*(Radial_function(r))*Radial_function_2(r))*r**2
   
    #each mass then works out the integrals 
    for F in range(len(other)):
        Radial_function_2= other[F]
        H=Energy[F]-Base_Energy
        for i in range(len(m_b_values)):
            m_boson = m_b_values[i]
            K=(quad(integ_function,0, numpy.inf)[0])**2
            y_matrix[i] +=K/H
    ratio=numpy.divide(y_matrix, pert_for_function(3,2,1))
    return ratio
    
pyplot.xlabel('$m_{\phi} (1/a_0)$', fontsize=40)
pyplot.ylabel('$dE_n^2/dE_n^1 (E_{\mathrm{h}})$', fontsize=40)
pyplot.rc('xtick', labelsize=30) 
pyplot.rc('ytick', labelsize=30) 

pyplot.yscale('log')
pyplot.xscale('log')
pyplot.show()

    
C='purple'
K='-'
pyplot.plot(m_b_values,second_order_10(f.R_10), color=C, linestyle=K, label='1    0', linewidth= 3.0)


C='blue'
K='-'
pyplot.plot(m_b_values,second_order_20(f.R_20), color=C, linestyle=K, label='2    0', linewidth= 3.0)
K='--'
pyplot.plot(m_b_values,second_order_21(f.R_21), color=C, linestyle=K, label='2    1', linewidth= 3.0)


C='red'
K='-'
pyplot.plot(m_b_values,second_order_30(f.R_30), color=C, linestyle=K, label='3    0', linewidth= 3.0)
K='--'
pyplot.plot(m_b_values,second_order_31(f.R_31), color=C, linestyle=K, label='3    1', linewidth= 3.0)
K=':'
pyplot.plot(m_b_values,second_order_32(f.R_32), color=C, linestyle=K, label='3    2', linewidth= 3.0)



leg=pyplot.legend(title='        n    l', fontsize=30)
pyplot.setp(leg.get_title(),fontsize='30')
pyplot.ylim(10**-4,10**-1)
pyplot.show()


