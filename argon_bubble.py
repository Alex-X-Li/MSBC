# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 14:20:17 2021

@author: Alex Lee
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.interpolate import CubicSpline

h = 111.085e3 ### enthalpy of vap kj/kg  to j/kg
cl = 1.56e3   ### heat capasity kj/...
rhol = 1.035e3  ### g/cm3 to kg/m3
kl = 70.4e-3  ### mw/m/k  to w/m/k
al = kl/rhol/cl  ### thermal diffusivity
rhov = 0.0905e3 ### g/cm3 to kg/m3
Tsat = 90.328  ## K saturation tem
vl = 0.00080382e-4  ## cm2/s to m2/s kinematic viscosity
# vl = 18.3e-6    ##LJ fluid
vl = 3.9*0.00080382e-4    ##test
Tsuper = 130
dT = Tsuper - Tsat
gamma = 3.0919e-3 ### mN/m to N/m

A = np.sqrt(2/3*(h*rhov*dT)/(rhol*Tsat))
Ja = (dT*cl*rhol)/(h*rhov)
B = np.sqrt(12*al/np.pi)*Ja


# A = 17.7
# B = 0.0036
# B = 0.0086 ## LJ fluid

def drdt(w,t, ts, v1, rho1, gamma):
    r = w
    drdt = -(A**2*np.sqrt(t-ts)/B + 2*v1/r) + np.sqrt( A**2-2*gamma/(rho1*r) + (2*v1/r + A**2*np.sqrt( t-ts )/B)**2 )
    return drdt

def d2rdt2(r,drdt, t, ts, v1, rho1, gamma):

    term1 = -A**2/(2*B*np.sqrt(t-ts))
    term2 = 2*v1*drdt/r**2
    term3n = 2*gamma*drdt/(rho1*r**2) + 2*(A**2*np.sqrt(t-ts)/B+2*v1/r)*(-term1-term2)
    term3d = 2*np.sqrt(A**2+(A**2*np.sqrt(t-ts)/B+2*v1/r)**2 - 2*gamma/(rho1*r))
    d2rdt = term1 + term2 + term3n/term3d
    return d2rdt


ts40 = 0.09986e-9
rs40 = 2.821705e-9

ts60 = 0.0797e-9
rs60 = 3.626e-9

ts70 = 0.0563e-9  ### 70eV in ns
rs70 = 3.795e-9   ### 70eV in nm

ts90 = 0.1032e-9
rs90 = 4.951e-9

ts160 = 0.3736846e-9
rs160 = 5.520985e-9


t70 = np.logspace(np.log10(ts70), 1, 1000000)
t40 = np.logspace(np.log10(ts40), 1, 1000000)
t60 = np.logspace(np.log10(ts60), 1, 1000000)
t90 = np.logspace(np.log10(ts90), 1, 1000000)
t160 = np.logspace(np.log10(ts160), 1, 1000000)

sol40 = odeint(drdt, rs40, t40, args=(ts40, vl, rhol, gamma))
sol40 = sol40.flatten()
sol60 = odeint(drdt, rs60, t60, args=(ts60, vl, rhol, gamma))
sol60 = sol60.flatten()
sol70 = odeint(drdt, rs70, t70, args=(ts70, vl, rhol, gamma))
sol70 = sol70.flatten()
sol90 = odeint(drdt, rs90, t90, args=(ts90, vl, rhol, gamma))
sol90 = sol90.flatten()
sol160 = odeint(drdt, rs160, t160, args=(ts160, vl, rhol, gamma))
sol160 = sol160.flatten()

rt40 = CubicSpline(t40, sol40)
rt60 = CubicSpline(t60, sol60)
rt70 = CubicSpline(t70, sol70)
rt90 = CubicSpline(t90, sol90)
rt160 = CubicSpline(t160, sol160)

# plt.figure(dpi = 200)
# plt.plot(t, rt(t))
# plt.xlabel('time [s]')
# plt.ylabel('radius [m]')
# plt.xscale('log')
# plt.yscale('log')
# plt.grid()

data40ev = np.loadtxt('40ev.csv',delimiter = ',')
data60ev = np.loadtxt('60ev.csv',delimiter = ',')
data70ev = np.loadtxt('70ev.csv',delimiter = ',')
data90ev = np.loadtxt('90ev.csv',delimiter = ',')
data160ev = np.loadtxt('160ev.csv',delimiter = ',')   #### 160ev 4rc ~ 80 ev 2rc
# t70ev = data70ev[:, 0]
# r70ev = data70ev[:, 1]

# t70ev = np.array(t70ev)*10**-9
# r70ev = np.array(r70ev)*10**-9

plt.figure(dpi = 200)
plt.plot(t40, rt40(t40))
plt.plot(t60, rt60(t60))
plt.plot(t70, rt70(t70))
plt.plot(t90, rt90(t90))
plt.plot(t160, rt160(t160))
plt.scatter(data40ev[:, 0]*1e-9, data40ev[:, 1]*1e-9, marker='.',label = '40eV')
plt.scatter(data60ev[:, 0]*1e-9, data60ev[:, 1]*1e-9, marker='.',label = '60eV')
plt.scatter(data70ev[:, 0]*1e-9, data70ev[:, 1]*1e-9, marker='.',label = '70eV')
plt.scatter(data160ev[:, 0]*1e-9, data160ev[:, 1]*1e-9, marker='.',label = '80eV')
plt.scatter(data90ev[:, 0]*1e-9, data90ev[:, 1]*1e-9, marker='.',label = '90eV')


plt.xlabel('time [s]')
plt.ylabel('radius [m]')
# plt.xscale('log')
# plt.yscale('log')
plt.grid()
plt.xlim(0, 2e-9)
plt.ylim(0, 20e-9)
plt.legend()

plt.figure(dpi = 200)
plt.plot(t40, rt40(t40)*100)
plt.plot(t60, rt60(t60)*100)
plt.plot(t70, rt70(t70)*100)
plt.plot(t160, rt160(t160)*100)
plt.plot(t90, rt90(t90)*100)

# plt.scatter(data40ev[:, 0]*1e-9, data40ev[:, 1]*1e-9, marker='.',label = '40eV')
# plt.scatter(data60ev[:, 0]*1e-9, data60ev[:, 1]*1e-9, marker='.',label = '60eV')
# plt.scatter(data70ev[:, 0]*1e-9, data70ev[:, 1]*1e-9, marker='.',label = '70eV')
# plt.scatter(data90ev[:, 0]*1e-9, data90ev[:, 1]*1e-9, marker='.',label = '90eV')
# plt.scatter(data160ev[:, 0]*1e-9, data160ev[:, 1]*1e-9, marker='.',label = '160eV')

plt.xlabel('time [s]')
plt.ylabel('radius [cm]')
# plt.xscale('log')
# plt.yscale('log')
plt.grid()
plt.xlim(0, 2)
plt.ylim(0, 0.6)
plt.legend()


def energy_vs_radius(time):
    energy = np.array([40, 60, 70, 80, 90])
    radius = np.array([rt40(time), rt60(time), rt70(time), rt160(time), rt90(time)])
    plt.figure(dpi=200)
    plt.plot(energy, radius)
    plt.scatter(energy, radius) 
    print(radius)
    plt.title('Energy vs radius in '+str(time)+ ' s')
    plt.xlabel('Energy [eV]')
    plt.ylabel('Radius [m]')


vol_data1 = np.loadtxt('comsol_V_0.04.txt', skiprows=5)
vol_data2 = np.loadtxt('comsol_V_0.029.txt', skiprows=5)
vol_data3 = np.loadtxt('comsol_V_0.031.txt', skiprows=5)
vol_data4 = np.loadtxt('comsol_V_0.033.txt', skiprows=5)
vol_data = (vol_data1[:, 2]+vol_data2[:, 2])/2
# r_data = np.loadtxt('edge.txt', skiprows=9)

# real_r = []
# i = 2
# while i < len(r_data[1, :]):
#     idx = np.argmin(r_data[:, i]<0.5)
    
#     real_r = np.append(real_r, r_data[idx, 0])
#     i = i+1


time_v = 0.0003711285159 + vol_data1[:, 1]*1e-3
rad_v29 = (3*2*vol_data2[:, 2]/4/np.pi)**(1/3)
rad_v29 = rad_v29*0.05/rad_v29[0]

rad_v31 = (3*2*vol_data3[:, 2]/4/np.pi)**(1/3)
rad_v31 = rad_v31*0.05/rad_v31[0]

rad_v33 = (3*2*vol_data4[:, 2]/4/np.pi)**(1/3)
rad_v33 = rad_v33*0.05/rad_v33[0]


plt.figure(dpi = 200)
plt.scatter(0.0003711285159, 0.05 , label = 'comsol start')
plt.plot(time_v, rad_v29, label = 'comsol R=0.029')
plt.plot(time_v, rad_v31, label = 'comsol R=0.031')
plt.plot(time_v, rad_v33, label = 'comsol R=0.033')
plt.plot(t160, rt160(t160)*1e3,  label = 'Model')
plt.xlabel('time [s]')
plt.ylabel('radius [mm]')
plt.xlim(0.35e-3, 1e-3)
plt.ylim(0.049, 0.85e-1)
plt.grid()
plt.legend()


