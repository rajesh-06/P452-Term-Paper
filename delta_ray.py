import numpy as np
import matplotlib.pyplot as plt
#import delta_ray
import math as m
import mm1
x = [1178992.1803226308, 1178992.1803226308, 1157278.3958906143, 978978.9645478921, 843689.354676642, 570995.2738808154, 393691.6037296812, 229622.95742485244, 107149.18928470994, 45561.01414228416, 17009.078269334474, 5679.699930004746, 3312.7185920604647, 1236.7233433331205, 470.363539974421, 215.44357086160005, 89.9216644033289, 48.68869058281399, 28.397966882844145, 14.274346247091062, 7.873968625510031, 5.1344736156142945, 3.34809808919345, 3.2864306947469553]
y = [0.000010112397325257372, 0.00001362397504416594, 0.000018700131582663175, 0.000029242839631164663, 0.0000416620591464391, 0.00008777388587317209, 0.00016536605059778108, 0.00034196405431648826, 0.0008056540926458443, 0.0020071917778517845, 0.005488860525268561, 0.01647510863420205, 0.02775656696652856, 0.07312695265401381, 0.19265909257304695, 0.42130381942792033, 0.992575159518511, 1.8016204126818478, 3.270115040593861, 6.762339645729585, 12.047755592835209, 18.49227012339181, 27.86016859643839, 28.91778548026951]
e=1.6e-19 #n Coulomb
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
I=I0*Z
z=1
rho=1.66#kg/m^3
def Beta_P(P,mass):#P=m.c.beta/(1-beta^2)^0.5
	return 1/(1+(mass/P)**2)**0.5

beta=Beta_P(1e9,m_p_ev)

def test(E):
	W=2*m.pi*N_A*r_e**2*m_e_ev*Z*z**2*1.66e-3/(A*beta**2)
	return W/E**2
def beta_ke(E,mass):#KE=(1/(1-beta^2)^0.5-1)*m*c^2
	return (1-1/(1+E/mass)**2)**0.5
def test1(E):
	E_m=2*m_e_ev*beta**2/(1-beta**2)
	W=2*m.pi*N_A*r_e**2*m_e_ev*Z*z**2*1.66e-3/(A*beta**2)
	return W*(1/E-1/E_m)
def dedx(P,mass):
	beta=beta_ke(P,mass)
	E_m=2*m_e_ev*beta**2/(1-beta**2)
	W=2*np.pi*N_A*r_e**2*m_e_ev*Z*z**2*rho*1e-3/(A*beta**2)
	return W*(np.log(2*m_e_ev*beta**2*E_m/(I**2*(1-beta**2)))-2*beta**2)

de_p=[]
de_mu=[]
de_e=[]
P=np.linspace(1,1e4,10000)
#ee=np.linspace(1e3,1e6,10000)
'''
for i in range(len(P)):
	de_p.append(dedx(P[i]*1e6,m_p_ev))
	de_mu.append(dedx(P[i]*1e6,105.658e6))
plt.plot(P,de_p,label="p")
plt.plot(P,de_mu,label="$\mu$")
#plt.yscale('log')
plt.xscale('log')
plt.xlabel("Energy (MeV)")
plt.ylabel("<dE/dX> (eV/cm)")
plt.ylim(0.01e5,0.06e6)
plt.legend()
plt.show()



E_m=2*m_e_ev*beta**2/(1-beta**2) #maximum allowed energy transfer in a single  collision
print(beta, test1(3),E_m)
print(test(1)*(1/3-1/E_m))
E=np.linspace(3,E_m,100000)
#plt.plot(E,mm1.simpson(test,E,E_m,100))
plt.plot(x,y,label='from sauli pdf')
#plt.plot(E,test1(E),label="from my code")
plt.plot(E,test1(E),label="Analytic")

plt.plot(E,mm1.simpson(test,E,E_m,50000),label="50000 points")

plt.plot(E,mm1.simpson(test,E,E_m,25000),label="25000 points")
plt.plot(E,mm1.simpson(test,E,E_m,10000),label="10000 points")
plt.plot(E,mm1.simpson(test,E,E_m,1000),label="1000 points")
plt.plot(E,mm1.simpson(test,E,E_m,100),label="100 points")
#plt.plot(E,mm1.trapezoidal(test,E,E_m,10000),label="trapezoidal")
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$E_0$(eV)')
plt.ylabel('N(E>$E_0$)/cm Ar')
plt.ylim(1e-5,100)
plt.legend()
plt.show()'''
E_m=2*m_e_ev*beta**2/(1-beta**2) #maximum allowed energy transfer in a single  collision
E=np.linspace(3,E_m,100000)
#plt.plot(E,mm1.simpson(test,E,E_m,100000),label="50000 points")
plt.plot(E,test1(E),label="Analytic")
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$E_0$(eV)')
plt.ylabel('N(E>$E_0$)/cm Ar')
plt.ylim(1e-5,100)
#plt.legend()
plt.show()