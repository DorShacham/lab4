import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
from sklearn.linear_model import LinearRegression as linreg


#%% I - V curve heater curent
def plot_I_V(Heater_Current,attempt):
    #Heater_Current = 0.28 #A
    fh1 = pd.read_csv(str('FH'+str(attempt)+'.csv'),sep='\t',header=5) # read the data.
    Va1 = np.array(fh1['Va(V)_1']) # accelerating voltage array
    I1 = np.array(fh1['Ia(E-12 A)_1']) # Current array
    T1 = np.array(fh1['T(c)_1']) #temperature array
    plt.figure()
    plt.plot(Va1,I1, label='Heater current {:.2f}[A]'.format(Heater_Current))
    plt.ylabel('Current [pA]')
    plt.xlabel('Acceleration voltage [V]')
    plt.grid()
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

#%% first current
I1,Va1=plot_I_V(0.28,1)
I2,Va2=plot_I_V(0.26,2)
I3,Va3=plot_I_V(0.25,3)
maximums1=find_local_maximums(I1,Va1)
maximums2=find_local_maximums(I2,Va2)
maximums3=find_local_maximums(I3,Va3)
            
        
# finding values
E_exc = np.mean(np.diff(maximums1)[2:])
V_contact = 5.8 - E_exc


#%% part 2 
I,V = plot_I_V(0.28,"_1_280_ion_good")
plt.plot(V[V>12],I[V>12])
#plt.yscale("log")
x,y = V[V>12],I[V>12]
x=x.reshape((-1,1))
model = linreg().fit(x,y)
print(model)
