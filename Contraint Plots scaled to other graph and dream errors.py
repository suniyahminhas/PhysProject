from __future__ import division
from sympy import var
from sympy import Symbol, integrate, pi, oo, sin, exp
from sympy.physics.hydrogen import R_nl
import numpy
import matplotlib.pyplot as pyplot
from scipy.misc import imread
from PIL import Image, ImageDraw



#attempt to make a constraint plot
h= 6.626070040*10**-34 #plancks constant, 
k= 27.21138602 #ev in one Hartree
e= 1.6021766208*10**-19 #electronic charge
M=0.5109989461 #mass of an electron Mev






r=Symbol("r", real=True, positive=True)
m_b_values= numpy.logspace(-4,6, num=100)

def pert_for_function((N,L,Z)):
    #this is what would be in the integral (i,e functions to integrate)
    data_matrix = numpy.zeros((len(m_b_values),1))
    for i in range(len(m_b_values)):
        m_boson = m_b_values[i]
        data = integrate((exp(-m_boson*r)/((4*pi*r))*(R_nl(N,L,r,Z))**2)*r**2, (r,0,oo))
        data_matrix[i] = float(data)
    return data_matrix
    
def constraint(a,b,error_bar):
    k_new= pert_for_function(a)-pert_for_function(b)
    error_E_bar= error_bar*h/(e*k)
    return error_E_bar/k_new

pyplot.figure()
pyplot.subplot(1,2,1)
img = imread("pythonreadyfade.png")
pyplot.imshow(img,zorder=0, extent= (10**-6, 10**3, 10**-16, 1), aspect='auto')

#pyplot.plot(m_b_values*(1/(2.6817268*(10**-4)))*(10**-6),constraint((1,0,79),(2,1,79),70000), color='darkgreen', label='Au: 1s-2p \n 70,000Hz', linewidth=4.0) 
pyplot.plot(m_b_values*(1/(2.6817268*(10**-4)))*(10**-6),constraint((1,0,1),(2,0,1),10), color='red', label=' 1s-2s', linewidth=4.0) 
#pyplot.plot(m_b_values*(1/(2.6817268*(10**-4)))*(10**-6),constraint((1,0,1),(2,0,1),1), color='blue', label='1Hz', linewidth=4.0) 
#pyplot.plot(m_b_values*(1/(2.6817268*(10**-4)))*(10**-6),constraint((1,0,1),(2,0,1),0.1), color='green', label='0.1Hz', linewidth=4.0) 
pyplot.plot(m_b_values*(1/(2.6817268*(10**-4)))*(10**-6),constraint((1,0,1),(2,1,1),4000*10**6), color='orange', label='1s-2p',linewidth=4.0) #NOTE THIS TIMISING FACTOR MEANS IT IS NOW IN EV
pyplot.plot(m_b_values*(1/(2.6817268*(10**-4)))*(10**-6),constraint((1,0,1),(3,1,1),9000*10**6), color='yellow', label='1s-3p', linewidth=4.0) 
pyplot.plot(m_b_values*(1/(2.6817268*(10**-4)))*(10**-6),constraint((1,0,1),(4,1,1),60000*10**6), color='green', label='1s-4p', linewidth=4.0) 
pyplot.plot(m_b_values*(1/(2.6817268*(10**-4)))*(10**-6),constraint((1,0,1),(5,1,1),13000*10**6), color='blue',label='1s-5p', linewidth=4.0) 
pyplot.plot(m_b_values*(1/(2.6817268*(10**-4)))*(10**-6),constraint((1,0,1),(6,1,1),70000*10**6), color='purple', label='1s-6p', linewidth=4.0) 
pyplot.plot(m_b_values*(1/(2.6817268*(10**-4)))*(10**-6),constraint((1,0,1),(7,1,1),70000*10**6), color='pink', label='1s-7p',linewidth=4.0) 
pyplot.plot(m_b_values*(1/(2.6817268*(10**-4)))*(10**-6),constraint((1,0,1),(8,1,1),70000*10**6), color='deeppink', label='1s-8p', linewidth=4.0) 
pyplot.xlabel('$m_\phi (MeV)$', fontsize=40)
pyplot.ylabel('$g_{e}g_{n}$', fontsize=40)
pyplot.rc('xtick', labelsize=30) 
pyplot.rc('ytick', labelsize=30) 
pyplot.yscale('log')
pyplot.ylim(10**-16, 1)
pyplot.xlim(10**-6, 10**3)
pyplot.xscale('log')
pyplot.xticks([10**-6,10**-3,1,10**3], fontsize=30)
pyplot.yticks([10**0,10**-4,10**-8, 10**-12, 10**-16], fontsize=30)

pyplot.subplot(1,2,2)
img = imread("pythonreadyfade.png")
pyplot.imshow(img,zorder=0, extent= (10**-6, 10**3, 10**-16, 1), aspect='auto')
pyplot.plot(m_b_values*(1/(2.6817268*(10**-4)))*(10**-6),constraint((1,0,1),(2,0,1),10), color='red', label='1s-2s', linewidth=4.0) 
pyplot.plot(m_b_values*(1/(2.6817268*(10**-4)))*(10**-6),constraint((1,0,1),(2,1,1),10), color='orange', label='1s-2p',linewidth=4.0) #NOTE THIS TIMISING FACTOR MEANS IT IS NOW IN MEV
pyplot.plot(m_b_values*(1/(2.6817268*(10**-4)))*(10**-6),constraint((1,0,1),(3,1,1),10), color='yellow', label='1s-3p', linewidth=4.0) 
pyplot.plot(m_b_values*(1/(2.6817268*(10**-4)))*(10**-6),constraint((1,0,1),(4,1,1),10), color='green', label='1s-4p', linewidth=4.0) 
pyplot.plot(m_b_values*(1/(2.6817268*(10**-4)))*(10**-6),constraint((1,0,1),(5,1,1),10), color='blue',label='1s-5p', linewidth=4.0) 
pyplot.plot(m_b_values*(1/(2.6817268*(10**-4)))*(10**-6),constraint((1,0,1),(6,1,1),10), color='purple', label='1s-6p', linewidth=4.0) 
pyplot.plot(m_b_values*(1/(2.6817268*(10**-4)))*(10**-6),constraint((1,0,1),(7,1,1),10), color='pink', label='1s-7p',linewidth=4.0) 
pyplot.plot(m_b_values*(1/(2.6817268*(10**-4)))*(10**-6),constraint((1,0,1),(8,1,1),10), color='deeppink', label='1s-8p', linewidth=4.0) 
pyplot.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 25)
pyplot.xlabel('$m_\phi (MeV)$', fontsize=40)
pyplot.ylabel('$g_{e}g_{n}$', fontsize=40)
pyplot.rc('xtick', labelsize=30) 
pyplot.rc('ytick', labelsize=30) 
pyplot.yscale('log')
pyplot.ylim(10**-16, 1)
pyplot.xlim(10**-6, 10**3)
pyplot.xscale('log')
pyplot.xticks([10**-6,10**-3,1,10**3], fontsize=30)
pyplot.yticks([10**0,10**-4,10**-8, 10**-12, 10**-16], fontsize=30)

#pyplot.tight_layout(pad=0.2, w_pad=0.2, h_pad=0.2)
pyplot.subplots_adjust(hspace=0.3, wspace=0.3)


pyplot.show()
