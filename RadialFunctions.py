from __future__ import division
import numpy
import matplotlib.pyplot as pyplot
from scipy.integrate import quad


#These functions were created using the sympy.physics.hydrogen they were turned into seperate functions so they could be integrated easily

N = numpy.arange(1,8)
data_points=30
range_of_r =20
r_range= numpy.arange(0,range_of_r,range_of_r/data_points)

def R_10(r):
    return 2*numpy.exp(-r)

def R_20(r):
    return numpy.sqrt(2)*(-r + 2)*numpy.exp(-r/2)/4

def R_21(r):
    return numpy.sqrt(6)*r*numpy.exp(-r/2)/12

def R_30(r):
    return 2*numpy.sqrt(3)*(2*r**2/9 - 2*r + 3)*numpy.exp(-r/3)/27

def R_31(r):
    return numpy.sqrt(6)*r*(-2*r/3 + 4)*numpy.exp(-r/3)/81

def R_32(r):
    return 2*numpy.sqrt(30)*r**2*numpy.exp(-r/3)/1215

def R_40(r):
    return (-r**3/768 + r**2/32 - 3*r/16 + 1/4)*numpy.exp(-r/4)

def R_41(r):
    return numpy.sqrt(15)*r*(r**2/8 - 5*r/2 + 10)*numpy.exp(-r/4)/480
    
def R_42(r):
    return numpy.sqrt(5)*r**2*(-r/2 + 6)*numpy.exp(-r/4)/1920

def R_43(r):
    return numpy.sqrt(35)*r**3*numpy.exp(-r/4)/26880

def R_50(r):
    return 2*numpy.sqrt(5)*(2*r**4/1875 - 4*r**3/75 + 4*r**2/5 - 4*r + 5)*numpy.exp(-r/5)/125

def R_51(r):
    return numpy.sqrt(30)*r*(-4*r**3/375 + 12*r**2/25 - 6*r + 20)*numpy.exp(-r/5)/1875

def R_52(r):
    return 2*numpy.sqrt(70)*r**2*(2*r**2/25 - 14*r/5 + 21)*numpy.exp(-r/5)/65625

def R_53(r):
    return numpy.sqrt(70)*r**3*(-2*r/5 + 8)*numpy.exp(-r/5)/328125

def R_54(r):
    return 2*numpy.sqrt(70)*r**4*numpy.exp(-r/5)/4921875

def R_60(r):
    return numpy.sqrt(6)*(-r**5/29160 + r**4/324 - 5*r**3/54 + 10*r**2/9 - 5*r + 6)*numpy.exp(-r/6)/108

def R_61(r):
    return numpy.sqrt(210)*r*(r**4/1944 - 7*r**3/162 + 7*r**2/6 - 35*r/3 + 35)*numpy.exp(-r/6)/11340

def R_62(r):
    return numpy.sqrt(105)*r**2*(-r**3/162 + 4*r**2/9 - 28*r/3 + 56)*numpy.exp(-r/6)/136080

def R_63(r):
    return numpy.sqrt(35)*r**3*(r**2/18 - 3*r + 36)*numpy.exp(-r/6)/1224720

def R_64(r):
    return numpy.sqrt(7)*r**4*(-r/3 + 10)*numpy.exp(-r/6)/7348320

def R_65(r):
    return numpy.sqrt(77)*r**5*numpy.exp(-r/6)/242494560

def R_70(r):
    return 2*numpy.sqrt(7)*(4*r**6/5294205 - 4*r**5/36015 + 2*r**4/343 - 20*r**3/147 + 10*r**2/7 - 6*r + 7)*numpy.exp(-r/7)/343

def R_71(r):
    return numpy.sqrt(21)*r*(-4*r**5/252105 + 16*r**4/7203 - 16*r**3/147 + 16*r**2/7 - 20*r + 56)*numpy.exp(-r/7)/7203

def R_72(r):
    return 2*numpy.sqrt(105)*r**2*(2*r**4/7203 - 12*r**3/343 + 72*r**2/49 - 24*r + 126)*numpy.exp(-r/7)/756315
    
def R_73(r):
    return numpy.sqrt(42)*r**3*(-4*r**3/1029 + 20*r**2/49 - 90*r/7 + 120)*numpy.exp(-r/7)/5294205

def R_74(r):
    return 2*numpy.sqrt(154)*r**4*(2*r**2/49 - 22*r/7 + 55)*numpy.exp(-r/7)/407653785

def R_75(r):
    return 2*numpy.sqrt(231)*r**5*(-2*r/7 + 12)*numpy.exp(-r/7)/8560729485

def R_76(r):
    return 4*numpy.sqrt(3003)*r**6*numpy.exp(-r/7)/779026383135

other= [R_20,R_30, R_40]
K= other[1]
print(K(4))
#pyplot.figure()
#pyplot.xlabel('r/a_0')
#pyplot.ylabel('Wavefunction Rnl')
#pyplot.yscale('log') #This can be removed simply shows what happens with a log scale
#pyplot.legend()
#pyplot.title('Radial Wavefunctions')
#pyplot.show()