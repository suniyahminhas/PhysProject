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


m_value = numpy.array([ 1.00000000*10**-1  , 1.95734178*10**-1 ,  3.83118685*10**-1  , 7.49894209*10**-1,
1.46779927,  2.87298483 ,  5.62341325 , 1.10069417*10**1,
   2.15443469*10**1 ,  4.21696503*10**1   ,8.25404185*10**1  , 1.61559810*10**2,
   3.16227766*10**2  , 6.18965819*10**2   ,1.21152766*10**3   ,2.37137371*10**3,
   4.64158883*10**3 ,  9.08517576*10**3  , 1.77827941*10**4  , 3.48070059*10**4,
   6.81292069*10**4,   1.33352143*10**5 ,  2.61015722*10**5 ,  5.10896977*10**5,
   1.00000000*10**6])
   

 
C_value=numpy.array([    2.60605270*10**-14,  2.73622125*10**-14, 3.10069202*10**-14, 4.04361988*10**-14, 6.43406953*10**-14, 1.28182801*10**-13, 3.15776788*10**-13,
9.21896981*10**-13, 3.02365885*10**-12, 1.06447527*10**-11,3.89988993*10**-11,1.45977686*10**-10, 5.52601432*10**-10, 2.10412694*10**-9, 8.03592886*10**-9,
3.07375342*10**-8,  1.17664188*10**-7, 4.50603586*10**-7,  1.72597488*10**-6,  6.61180489*10**-6,  2.53296362*10**-5,  9.70397827*10**-5, 
3.71772085*10**-4,  1.42431765*10**-3, 5.45680615*10**-3])




h= 6.626070040*10**-34 #plancks constant, 
k= 27.21138602 #ev in one Hartree
e= 1.6021766208*10**-19 #electronic charge

def constraint_compare(Z, masses, Cop, A):
    #this is what would be in the integral (i,e differents to integrate)
    data_1, data_2 =numpy.zeros(len(masses)),numpy.zeros(len(masses))
    for i in range(len(masses)):
        m_boson=masses[i]
        C=Cop[i]
        data_1[i] = C*integrate((exp(-m_boson*r)/((4*pi*r))*(R_nl(1,0,r,Z))**2)*r**2, (r,0,oo))
        data_2[i] = C*integrate((exp(-m_boson*r)/((4*pi*r))*(R_nl(2,0,r,Z))**2)*r**2, (r,0,oo))
    y = numpy.subtract(data_1, data_2)
    pyplot.plot(masses, A*((y*e*k)/h), label=str(Z), linestyle='-', linewidth=4.0)
    leg = pyplot.legend(title='Z', numpoints=1, loc=2, fontsize=30)
    pyplot.setp(leg.get_title(),fontsize=30)

constraint_compare(2, m_value, C_value,4)
constraint_compare(3, m_value, C_value,7)
constraint_compare(4, m_value, C_value,9)
constraint_compare(5, m_value, C_value,11)
constraint_compare(6, m_value, C_value,12)
constraint_compare(7, m_value, C_value,14)
constraint_compare(8, m_value, C_value,16)



pyplot.xlabel('$m_\phi (1/a_\mathrm{0})$', fontsize=40)
pyplot.ylabel('$f (\mathrm{kHz})$', fontsize=40)
pyplot.xlim(10**-1,10**6)
pyplot.xscale('log')
pyplot.rc('xtick', labelsize=30)
pyplot.rc('ytick', labelsize=30)
pyplot.xticks(fontsize=30)
pyplot.yticks([10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000], [10,20,30,40,50,60,70,80,90], fontsize=30)
pyplot.show()