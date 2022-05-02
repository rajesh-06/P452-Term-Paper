import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 
import mm1
x=np.linspace(0,70,1000)
def landau(x,A,dx,meanE):
	lam=(x-dx)/meanE
	return A*np.exp(-0.5*(lam+np.exp(-lam)))
def gauss(x,A,mu,sig):
	lam=(x-mu)/sig
	return A*np.exp(-0.5*(lam)**2)
def readfile(filename,start):
    with open(filename,"r+") as f:
        lines = f.readlines()
        A=[]
        for i in range(start,len(lines)):
            A.append([float(j) for j in lines[i].split(",")])
        del lines
        return A

data=readfile('landau_data.csv',0)
data=mm1.transpose(data)

la,lb=curve_fit(landau,data[0],data[1])
ga,gb=curve_fit(gauss,data[0],data[1])
N=len(data[0])
df=N-3
print(df)
chi2_l=0
chi2_g=0
for i in range(6,N-13):
    chi2_l+=(data[1][i]-landau(data[0][i],la[0],la[1],la[2]))**2/landau(data[0][i],la[0],la[1],la[2])
    #print(landau(data[0][i],la[0],la[1],la[2]))
print()
for i in range(4,N-20):
    chi2_g+=(data[1][i]-gauss(data[0][i],ga[0],ga[1],ga[2]))**2/gauss(data[0][i],ga[0],ga[1],ga[2])
    #print(gauss(data[0][i],ga[0],ga[1],ga[2]))
print(chi2_l,chi2_g)

plt.plot(data[0],data[1],"o",label="Data")
plt.plot(x,landau(x,*la),label="Landau\n$\chi^2$/ndf="+str(round(chi2_l,2))+"/"+str(df))
plt.plot(x,gauss(x,*ga),label="Gauss\n$\chi^2$/ndf.="+str(round(chi2_g,2))+"/"+str(df))
plt.xlabel("Pulse Height-random unit")
plt.vlines(x = la[1], ymin = 0, ymax = 260,
           colors = 'black', linestyles = 'dashed',
           )
plt.vlines(x = ga[1], ymin = 0, ymax = 250, colors = 'black',linestyles = 'dashed')
plt.ylabel("Counts")
plt.legend()
plt.show()