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

h= 6.626070040*10**-34 #plancks constant, 
k= 27.21138602 #ev in one Hartree
e= 1.6021766208*10**-19 #electronic charge

n1= numpy.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2])
l1=numpy.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0])
n2=numpy.array([2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,2,2,2,4,4,4,4,8,8,8,10,12,12])
l2=numpy.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,2,0,2,2,2,2,2])
f_exp=numpy.array([2466068000000000,2466068000000000,2922729000000000,2922729000000000,3082570000000000,3082570000000000,3156567000000000,3156567000000000,3196760000000000,3196760000000000,3220980000000000,3220980000000000,3236700000000000,3236700000000000,3247490000000000,3247490000000000,3255180000000000,3255180000000000,3260920000000000,3260920000000000,3265250000000000,3265250000000000,1057847000,10969130000,9911201000,616520017568000,616520150636000,616521388672000,616521843441000,770649350012000,770649504450000,770649561584000,789144886410000,799191710473000,799191727404000])
error=numpy.array([4000000000,4000000000,9000000000,9000000000,60000000000,60000000000,13000000000,13000000000,70000000000,70000000000,70000000000,70000000000,70000000000,70000000000,70000000000,70000000000,70000000000,70000000000,70000000000,70000000000,70000000000,70000000000,9000,100000,12000,15000,10000,10000,24000,9000,8000,7000,40000,10000,7000])
f_th=numpy.array([2466060355341200,2466071324385600,2922742963793300,2922746213883400,3082581430738000,3082582801868100,3156563616462900,3156564318480700,3196751390971000,3196751797231000,3220983314718000,3220983570554000,3236710746537000,3236710917916000,3247493411755000,3247493532128000,3255206182761000,3255206270512000,3260912757359000,3260912823288000,3265253072974000,3265253123756000,1057845900,10969044400,9911198500,616520017550900,616520150628500,616521388681000,616521843427100,770649350016000,770649504449000,770649561578000,789144886426000,799191710482000,799191727409000])

print len(n1)
Observed=numpy.subtract(f_exp, f_th)

r=Symbol("r", real=True, positive=True)
def Res(m_phi, C):
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

def Stat(m_phi, CON):
    y_values = numpy.zeros(len(CON))
    for i in range(len(CON)):
        C=CON[i]
        print C
        U, D= 0, 0
        R= Res(m_phi,C)
        for j in range(len(R)-1):
            U += (R[j+1]-R[j])**2
        for j in range(len(R)):
            D += (R[j])**2
        y_values[i]= U/D
    return y_values 

maas1=1*10
con1=numpy.arange(10**-12,4*10**-11,5*10**-12)

maas2=1*10**2
con2=numpy.arange(10**-10,4*10**-9,5*10**-10)

maas3=1*10**3
con3=numpy.arange(10**-8,4*10**-7,5*10**-8)

maas4=1*10**4
con4=numpy.arange(10**-6, 4*10**-5, 5*10**-6)

maas5=1*10**5
con5= numpy.arange(10**-4, 4*10**-3, 5*10**-4)

maas6=1*10**6
con6= numpy.arange(10**-2, 4*10**-1, 5*10**-2)

pyplot.figure()

pyplot.subplot(1,3,1)
pyplot.plot(con1, Stat(maas1, con1), linewidth=3.0)
pyplot.title('$m_{\phi}= 10$\n  ', fontsize=35)
pyplot.xlabel('$g_{e}g_{N}$', fontsize=30)
pyplot.ylabel('$\mathrm{DWS}$\n' '$(10^{-6})+1.304513$', fontsize=30)
pyplot.xticks(fontsize=30)
pyplot.yticks([1.304515, 1.304517, 1.304519],[2,4,6], fontsize=30)
pyplot.xscale('log')
pyplot.xlim(10**-12,4*10**-11)

pyplot.subplot(1,3,2)
pyplot.plot(con2, Stat(maas2, con2), linewidth=3.0)
pyplot.title('$m_{\phi}= 100$ \n  ', fontsize=35)
pyplot.xlabel('$g_{e}g_{N}$', fontsize=30)
pyplot.ylabel('$\mathrm{DWS}$\n' '$(10^{-7})+1.3045$', fontsize=30)
pyplot.yticks([1.304506, 1.304512, 1.304518],[0.6,1.2,1.8], fontsize=30)
pyplot.xticks(fontsize=30)
pyplot.xscale('log')
pyplot.xlim(10**-10,4*10**-9)

pyplot.subplot(1,3,3)
pyplot.plot(con3, Stat(maas3, con3), linewidth=3.0)
pyplot.title('$m_{\phi}= 10^{3}$\n   ', fontsize= 35)
pyplot.xlabel('$g_{e}g_{N}$', fontsize=30)
pyplot.ylabel('$\mathrm{DWS}$\n' ' $(10^{-5}) +1.30449$', fontsize=30)
pyplot.yticks([1.30450, 1.30451,1.30452],[1, 2, 3], fontsize=30)
pyplot.xticks(fontsize=30)
pyplot.xscale('log')
pyplot.xlim(10**-8,4*10**-7)

pyplot.subplots_adjust(hspace=0.5, wspace=0.5)

pyplot.show()