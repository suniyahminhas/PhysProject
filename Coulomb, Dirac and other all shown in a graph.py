from __future__ import division
import numpy
import matplotlib.pyplot as pyplot
from scipy.integrate import quad
#import RadialFunctions as f
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
t[2]= 2466068000
t[3]= 2922729000 
t[4]= 3082570000
t[5]= 3156567000
t[6]= 3196760000
t[7]= 3220980000
t[8]= 3236700000
t[9]= 3247490000
t[10]= 3255180000
t[11]= 3260920000
t[12]= 3265250000

t=t*10**6
#This as been scaled in Hz

er=numpy.zeros((13,1))

er[2]=4000
er[3]=9000
er[4]=60000
er[5]=13000
er[6]=70000
er[7]=70000
er[8]= 70000
er[9]= 70000
er[10]=70000
er[11]=70000
er[12]=70000

er= er*10**6 #scaling so now in Hz

r=Symbol("r", real=True, positive=True)
m_b_values= numpy.arange(0,6,0.1)
data_matrix = numpy.zeros((len(m_b_values),1))
def pert_for_function(N,L,Z):
    #this is what would be in the integral (i,e functions to integrate)
    for i in range(len(m_b_values)):
        m_boson = m_b_values[i]
        data = integrate((exp(-m_boson*r)/((4*pi*r))*(R_nl(N,L,r,Z))**2)*r**2, (r,0,oo))
        data_matrix[i] = float(data)
    return data_matrix

#CODATA values

h= 6.626070040*10**-34 #plancks constant, 
k= 27.21138602 #ev in one Hartree
e= 1.6021766208*10**-19 #electronic charge
c=1/ (7.2973525664*10**-3)
alp=7.2973525664*10**-3
#all given in hartree units

th =numpy.array([0,0, 2466060355.3412, 2922742963.7933,  3082581430.7380, 3156563616.4629, 3196751390.971,3220983314.718, 3236710746.537, 3247493411.755, 3255206182.761,3260912757.359, 3265253072.974]) 
th= th*10**6

M= 1/ (5.44617021352*10**-4)
mu=M/(1+M)

def Energy_Dirac(n,j):
    return mu*c**2*((1+(alp/(n-(j+1/2)+((j+1/2)**2-alp**2)**(1/2)))**2)**(-1/2)-1)

def data_vs_theory(data):
    E=numpy.zeros((13,3))
    for n in numpy.arange(2,13):
        E[n,0]= h*data[n]/(e*k)
        E[n,1] = (E_nl(n,1)-E_nl(1,1))*(mu)
        E[n,2]= -E[n,1]+E[n,0] #Experiment- theory
    return E
    
def data_vs_theory_dir(data):
    E=numpy.zeros((13,3))
    for n in numpy.arange(2,13):
        E[n,0]= h*data[n]/(e*k)
        E[n,1] = Energy_Dirac(n,1/2)-Energy_Dirac(1,1/2) #experimental differnce
        E[n,2]= -E[n,1]+E[n,0] # experiment minus theory
    return E

def data_vs_quantum(data):
    E=numpy.zeros((13,3))
    for n in numpy.arange(2,13):
        E[n,0]= h*data[n]/(e*k)
        E[n,1] = h*th[n]/(e*k)
        E[n,2]= E[n,0]-E[n,1] #experiment minus theory
    return E
    

errors=h*er/(e*k)
y_er= numpy.zeros(11)
print y_er

pyplot.figure()
pyplot.plot(numpy.arange(2,13), data_vs_theory(t)[2:,2], marker='o',markersize=10, color='red', linestyle="None", label='Coulomb')
pyplot.plot(numpy.arange(2,13), data_vs_theory_dir(t)[2:,2], marker='x',markersize=15, color='blue', linestyle="None", label='Dirac equation')
pyplot.plot(numpy.arange(2,13), data_vs_quantum(t)[2:,2], marker='*', markersize=15, color='green', linestyle="None", label='QED corrections')
pyplot.errorbar(numpy.arange(2,13), y_er, yerr=errors[2:], color='grey', linestyle="None", linewidth=3.0)
leg = pyplot.legend(loc=3, numpoints=1, title= '$\Delta E_{n}^{\mathrm{the}}$ calculated by:', fontsize= 17)
pyplot.setp(leg.get_title(),fontsize=17)

#pyplot.scatter(numpy.arange(2,5), (fs[2:5]*h)/(e*k), label='fine structure')
#pyplot.legend(loc='center left', bbox_to_anchor=(1, 0.5))
pyplot.axhline(y=0, color='black', linestyle=':')
pyplot.xlabel('$n$', fontsize=40)
pyplot.ylabel('$\Delta E_{n}^{\mathrm{exp}}-\Delta E_{n}^{\mathrm{the}}$\n $(10^{-5}E_\mathrm{h})$', fontsize=30)
pyplot.xticks([1,2,3,4,5,6,7,8,9,10,11,12],['$1$','$2$','$3$','$4$','$5$','$6$','$7$','$8$','$9$','$10$','$11$','$12$'])
pyplot.yticks([1.5*10**-5,1*10**-5,5*10**-6, 0, -5*10**-6,-1*10**-5, -1.5*10**-5], 
[r'$1.5$',r'$1.0$', r'$0.5$', '$0.0$', r'$-0.5$',r'$-1.0$', r'$-1.5$'])
pyplot.xlim(1.5,12.5)
pyplot.rc('xtick', labelsize=40) 
pyplot.rc('ytick', labelsize=40) 
pyplot.show()