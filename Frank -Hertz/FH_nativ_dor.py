import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
from sklearn.linear_model import LinearRegression as linreg
from scipy.stats import linregress
from scipy import constants
from uncertainties import unumpy
from uncertainties import ufloat
from uncertainties.umath import *

#%% I - V curve heater curent
def plot_I_V(Heater_Current,Vr,attempt):
    #Heater_Current = 0.28 #A
    fh1 = pd.read_csv(str('FH'+str(attempt)+'.csv'),sep='\t',header=5) # read the data.
    Va1 = np.array(fh1['Va(V)_1']) # accelerating voltage array
    I1 = np.array(fh1['Ia(E-12 A)_1']) # Current array
    T1 = np.array(fh1['T(c)_1']) #temperature array
    plt.plot(Va1,I1,".", label=f'Heater current {Heater_Current}[A], Vr {Vr}[V]',markersize=2)
    plt.ylabel('Current [pA]')
    plt.xlabel('Acceleration voltage [V]')
    plt.grid(True)
    plt.legend()
    return I1,Va1



#%% find picks
def find_local_maximums(array,V):
    n=len(array)
    maximums =[]
    for i in range(1,n-1):
        if array[i] > array[i-1] and array[i] > array[i+1]:
            maximums.append(V[i])
    return maximums

#%% finding values
plt.figure(dpi=300)
I1,Va1=plot_I_V(0.26,1.5,"_1.5_260")
plt.savefig("fig/_1.5_260.jpg")

maximums1=np.array(find_local_maximums(I1,Va1))
maximums1 = maximums1[maximums1>5]
E_exc = np.mean(np.diff(maximums1))
V_contact = ufloat(maximums1[0] - E_exc,0.2)


print("E_exc",E_exc,"\nV_contact",V_contact)

h = constants.physical_constants["Planck constant in eV/Hz"][0]
c = constants.c #meter/sec
Lambda = 253.65e-9 #meter
E_exc_theory = h*c/Lambda
rel_err = abs(E_exc_theory-E_exc)/E_exc_theory

print("rel error:",rel_err*100,"%")
#%% different parmeters - current
plt.figure(dpi=300)
plot_I_V(0.28,1.5,"_1.5_280")
plot_I_V(0.26,1.5,"_1.5_260")
plot_I_V(0.25,1.5,"_1.5_250")
plt.savefig("fig/graph2.jpg")


#%% different parmeters - Vr
plt.figure(dpi=300)
plot_I_V(0.26,1,"_1_260")
plot_I_V(0.26,1.5,"_1.5_260")
plot_I_V(0.26,2,"_2_260")            
plot_I_V(0.26,5.5,"_5.5_260")        
plt.savefig("fig/graph3.jpg")

#%% part 2 
plt.figure(dpi=300)
I,V = plot_I_V(0.28,1,"_1_280_ion_good")
#plt.plot(V[V>12],I[V>12])
plt.legend()
#plt.yscale("log")
TH = 11.7
x,y = V[V>TH],I[V>TH]
#x=x.reshape((-1,1))
reg1 = linregress(x,y)
#model = linreg().fit(x,y)
slope = ufloat(reg1.slope,reg1.stderr)
intercept = ufloat(reg1.intercept,reg1.intercept_stderr)

E_ion = -intercept/slope - V_contact
E_ion_theory = 10.4375 #ev
rel_err = abs(E_ion_theory-E_ion)/E_ion_theory
print("E_ion",E_ion)
print("rel error:",rel_err*100,"%")
print(f"regression: y=({reg1.slope}+-{reg1.stderr})x+({reg1.intercept}+-{reg1.intercept_stderr})")

x = np.linspace(-reg1.intercept/reg1.slope,np.max(V),endpoint=True)
y = reg1.slope * x + reg1.intercept
plt.plot(x,y,"-.",label="regression")
plt.legend()
plt.savefig("fig/graph4.jpg")
