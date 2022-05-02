import numpy as np
import matplotlib.pyplot as plt
#import delta_ray
import math as m
import mm1

N_A=6.022e23
m_e_ev=0.510998e6 #mass of electron in eV
m_e_kg=9.10938e-31 #mass of electron in kg
r_e=2.8179e-13# in cm
m_p_ev=938.27e6 #mass of proton in eV
m_p = 1.6726219e-27 #mass of proton in kg
epsilon_0=8.854e-12 #F/m
c = 299792458 #speed of light in m/s
#For Argaon
Z=18#atomic number
A=40#atomic mass
I0=15.76#ionisation potential of Argon (in eV)
I=10.4*18
z=1
rho=1.66#kg/m^3

def Beta_P(P,mass):#P=m.c.beta/(1-beta^2)^0.5
	return 1/(1+(mass/P)**2)**0.5
def beta_ke(E,mass):#KE=(1/(1-beta^2)^0.5-1)*m*c^2
	return (1-1/(1+E/mass)**2)**0.5
beta = beta_ke(1e9,m_p_ev)
print(beta)
def dE_dx(E,mass):
	beta = Beta_P(E,mass)
	E_m=2*m_e_ev*beta**2/(1-beta**2)
	W=2*m.pi*N_A*r_e**2*m_e_ev*Z*z**2*1.66e-3/(A*beta**2)
	dE=W*(np.log(2*mass*beta**2*E_m/(I**2*(1-beta**2)))-2*beta**2)
	return dE
E=np.linspace(1.2e9,1.2e12,100000)
I_I_0=[]
for i in range(len(E)):
	I_I_0.append(dE_dx(E[i],m_p_ev))
I_0=min(I_I_0)
print(I_0)
print(sum(I_I_0)/len(I_I_0))
for i in range(len(E)):
	I_I_0[i]=I_I_0[i]/I_0

def readfile(filename,start):
    with open(filename,"r+") as f:
        lines = f.readlines()
        A=[]
        for i in range(start,len(lines)):
            A.append([float(j) for j in lines[i].split(",")])
        del lines
        return A

data=readfile("dedx_data_sauli.csv",0)
data=mm1.transpose(data)




plt.plot(E,I_I_0,label="Calculated")
plt.plot(data[0],data[1],label="Experimental fit data")
#plt.plot(E,dE_dx(E,m_p_ev))
#plt.yscale("log")
plt.xscale("log")
plt.xlabel("P(eV/c)")
plt.ylabel("I/$I_0$")
plt.legend()

plt.xlim(1e9,10e12)
plt.show()