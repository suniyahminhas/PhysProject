from __future__ import division
import numpy
import matplotlib.pyplot as pyplot
from sympy import var
from sympy import Symbol, integrate, pi, oo, sin, exp
from sympy.physics.hydrogen import R_nl
from sympy.physics.hydrogen import E_nl
import numpy
import scipy.stats
import pandas as pd


#can work with m_b values from 1/100a_0 to (10**4)a_0

h= 6.626070040*10**-34 #plancks constant, 
k= 27.21138602 #ev in one Hartree
e= 1.6021766208*10**-19 #electronic charge
#this needs some kind of order


maas=1.29154967*10**-3
con1=2.54651204*10**-14

n1= numpy.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2])
l1=numpy.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0])
n2=numpy.array([2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,2,2,2,4,4,4,4,8,8,8,10,12,12])
l2=numpy.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,2,0,2,2,2,2,2])
f_exp=numpy.array([2466068000000000,2466068000000000,2922729000000000,2922729000000000,3082570000000000,3082570000000000,3156567000000000,3156567000000000,3196760000000000,3196760000000000,3220980000000000,3220980000000000,3236700000000000,3236700000000000,3247490000000000,3247490000000000,3255180000000000,3255180000000000,3260920000000000,3260920000000000,3265250000000000,3265250000000000,1057847000,10969130000,9911201000,616520017568000,616520150636000,616521388672000,616521843441000,770649350012000,770649504450000,770649561584000,789144886410000,799191710473000,799191727404000])
error=numpy.array([4000000000,4000000000,9000000000,9000000000,60000000000,60000000000,13000000000,13000000000,70000000000,70000000000,70000000000,70000000000,70000000000,70000000000,70000000000,70000000000,70000000000,70000000000,70000000000,70000000000,70000000000,70000000000,9000,100000,12000,15000,10000,10000,24000,9000,8000,7000,40000,10000,7000])
f_th=numpy.array([2466060355341200,2466071324385600,2922742963793300,2922746213883400,3082581430738000,3082582801868100,3156563616462900,3156564318480700,3196751390971000,3196751797231000,3220983314718000,3220983570554000,3236710746537000,3236710917916000,3247493411755000,3247493532128000,3255206182761000,3255206270512000,3260912757359000,3260912823288000,3265253072974000,3265253123756000,1057845900,10969044400,9911198500,616520017550900,616520150628500,616521388681000,616521843427100,770649350016000,770649504449000,770649561578000,789144886426000,799191710482000,799191727409000])

Observed=numpy.subtract(f_exp, f_th)

r=Symbol("r", real=True, positive=True)

def Lag_plot(m_phi, C):
    def pert_for_different(N_1,L_1,Z_1,N_2, L_2, Z_2,):
        #this is what would be in the integral (i,e differents to integrate)
        data_1 = C*integrate((exp(-m_phi*r)/((4*pi*r))*(R_nl(N_1,L_1,r,Z_1))**2)*r**2, (r,0,oo))
        data_2 = C*integrate((exp(-m_phi*r)/((4*pi*r))*(R_nl(N_2,L_2,r,Z_2))**2)*r**2, (r,0,oo))
        return float(data_1 - data_2)

    Expected = numpy.zeros(len(n1))
    for i in range(len(n1)):
        Expected[i] = ((e*k)/h)*(pert_for_different(n1[i], l1[i],1, n2[i], l2[i], 1))

    diff = numpy.subtract(Observed, Expected)
    res_i = numpy.true_divide(diff, error)
    return res_i


def plot_square():
    pyplot.plot((-2,2),(-2,-2),color= 'black', linestyle='--')
    pyplot.plot((-2,2),(2,2), color= 'black', linestyle='--')
    pyplot.plot((-2,-2),(-2,2), color='black', linestyle='--')
    pyplot.plot((2,2),(-2,2), color='black', linestyle='--')


pyplot.subplot(2,3,1)
pyplot.scatter(Lag_plot(maas ,5*(10**-9))[1:], Lag_plot(maas ,5*(10**-9))[:-1])
pyplot.xlabel('$R_{n-1}$', fontsize=20, labelpad=-5)
pyplot.ylabel('$R_{n}$', fontsize=20, labelpad=-5)
plot_square()
pyplot.title(r'$g_{e}g_{N} =5 \times 10^{-9}$', fontsize=20)

pyplot.subplot(2,3,2)
pyplot.scatter((Lag_plot(maas ,10**-9))[1:], Lag_plot(maas ,10**-9)[:-1])
pyplot.xlabel('$R_{n-1}$', fontsize=20, labelpad=-5)
pyplot.ylabel('$R_{n}$', fontsize=20, labelpad=-5)
plot_square()
pyplot.title(r'$g_{e}g_{N} = 10^{-9}$', fontsize=20)

pyplot.subplot(2,3,3)
pyplot.scatter(Lag_plot(maas ,5*10**-10)[1:], Lag_plot(maas ,5*10**-10)[:-1])
pyplot.xlabel('$R_{n-1}$', fontsize=20, labelpad=-5)
pyplot.ylabel('$R_{n}$', fontsize=20, labelpad=-5)
pyplot.title( r'$g_{e}g_{N} =  5 \times 10^{-10}$', fontsize=20)
plot_square()

pyplot.subplot(2,3,4)
pyplot.scatter(Lag_plot(maas ,10**-10)[1:], Lag_plot(maas ,10**-10)[:-1])
pyplot.xlabel('$R_{n-1}$\n' r'$g_{e}g_{N} =  10^{-10}$', fontsize=20, labelpad=-5)
pyplot.ylabel('$R_{n}$', fontsize=20, labelpad=-5)
plot_square()

pyplot.subplot(2,3,5)
pyplot.scatter(Lag_plot(maas ,con1)[1:], Lag_plot(maas ,con1)[:-1])
pyplot.ylabel('$R_{n}$', fontsize=20, labelpad=-5)
pyplot.xlabel('$R_{n-1}$\n' r'$g_{e}g_{N} = 2.55 \times 10^{-14}$', fontsize=20, labelpad=-5)
plot_square()

pyplot.subplot(2,3,6)
res_sm = numpy.true_divide(numpy.subtract(f_exp, f_th), error)
res_sm_im1 = res_sm[1:]
pyplot.scatter(res_sm_im1, res_sm[:-1])
pyplot.xlabel('$R_{n-1}$\n$g_{e}g_{N} = 0$', fontsize=20, labelpad=-5)
pyplot.ylabel('$R_{n}$', fontsize=20, labelpad=-5)
plot_square()




pyplot.show()


