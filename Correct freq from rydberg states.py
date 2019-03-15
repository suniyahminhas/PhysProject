from __future__ import division
from sympy import var
from sympy import Symbol, integrate, pi, oo, sin, exp
from sympy.physics.hydrogen import R_nl
import numpy
import matplotlib.pyplot as pyplot

r=Symbol("r", real=True, positive=True)

pyplot.figure()

#m_boson = 10.2388
#C= 0.00031289
#m_boson = 100054
#C= 18212.5   #from 1s to 2p


m_value = numpy.array([1.00000000*10**-4 ,  1.95734178*10**-4  , 3.83118685*10**-4 ,  7.49894209*10**-4,
   1.46779927*10**-3   ,2.87298483*10**-3   ,5.62341325*10**-3   ,1.10069417*10**-2,
   2.15443469*10**-2   ,4.21696503*10**-2   ,8.25404185*10**-2   ,1.61559810*10**-1,
   3.16227766*10**-1  , 6.18965819*10**-1  , 1.21152766*10**0  , 2.37137371*10**0,
   4.64158883*10**0 ,  9.08517576*10**0 ,  1.77827941*10**1 ,  3.48070059*10**1,
   6.81292069*10**1,   1.33352143*10**2,   2.61015722*10**2,   5.10896977*10**2,
   1.00000000*10**3])


 
C_value=numpy.array([     2.54649943*10**-14, 2.54649964*10**-14, 2.54650047*10**-14,2.54650364*10**-14,
 2.54651574*10**-14, 2.54656189*10**-14, 2.54673707*10**-14, 2.54739644*10**-14, 2.54983843*10**-14, 2.55862113*10**-14,
 2.58866745*10**-14,  2.68391539*10**-14, 2.95814615*10**-14, 3.68060397*10**-14,
5.50882499*10**-14, 1.02887682*10**-13,  2.39332274*10**-13, 6.69232893*10**-13,  2.13431465*10**-12,
7.39139819*10**-12, 2.68358504*10**-11, 9.99682138*10**-11, 3.77484266*10**-10, 1.43547857*10**-9, 5.47863211*10**-9
])


h= 6.626070040*10**-34 #plancks constant, 
k= 27.21138602 #ev in one Hartree
e= 1.6021766208*10**-19 #electronic charge

def constraint_compare(n, masses, Cop):
    #this is what would be in the integral (i,e differents to integrate)
    data_1, data_2 =numpy.zeros(len(masses)),numpy.zeros(len(masses))
    for i in range(len(masses)):
        m_boson=masses[i]
        C=Cop[i]
        data_1[i] = C*integrate((exp(-m_boson*r)/((4*pi*r))*(R_nl(1,0,r,1))**2)*r**2, (r,0,oo))
        data_2[i] = C*integrate((exp(-m_boson*r)/((4*pi*r))*(R_nl(n,0,r,1))**2)*r**2, (r,0,oo))
    y = numpy.subtract(data_1, data_2)
    pyplot.plot(masses, ((y*e*k)/h), label=str(n), linestyle='-', linewidth=4.0)
    leg = pyplot.legend(title='n', numpoints=1, loc= 1, fontsize=30)
    pyplot.setp(leg.get_title(),fontsize=30)

constraint_compare(3, m_value, C_value)
constraint_compare(4, m_value, C_value)
constraint_compare(5, m_value, C_value)
constraint_compare(6, m_value, C_value)
constraint_compare(7, m_value, C_value)
constraint_compare(8, m_value, C_value)
constraint_compare(9, m_value, C_value)



pyplot.xlabel('$m_\phi (1/a_\mathrm{0})$', fontsize=40)
pyplot.ylabel('$f (\mathrm{Hz})$', fontsize=40)
pyplot.xlim(10**-4,10**3)
pyplot.xscale('log')
pyplot.xticks(fontsize=30)
pyplot.yticks(fontsize=30)
#pyplot.rc('xtick', labelsize=30)
#pyplot.rc('ytick', labelsize=30)
pyplot.show()