# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 16:00:54 2021

@author: Alex Lee
"""

import numpy as np
import matplotlib.pyplot as plt

def Efield_data(filename, low_z, high_z, zstep, low_r, high_r, rstep, limt ):
    
    coil_data = np.loadtxt(filename, skiprows = 8)
    l_r = low_r
    h_r = high_r
    # rstep = 2
    l_z = low_z
    h_z = high_z
    # zstep = 5
    
    radius = np.arange(l_r,h_r+rstep,rstep)
    zshift = np.arange(l_z,h_z+zstep, zstep)
    comb = np.array(np.meshgrid(radius, zshift))
    comb = comb.T.reshape(-1, 2)
    
    z = coil_data[:, 0]
    
    f_vol = []
    
    i = 1
    while i < len(coil_data[0, :]):
        idx = np.where(  np.abs(coil_data[:, i] - coil_data[int(len(coil_data[:, i])/2), i] ) < limt*coil_data[int(len(coil_data[:, i])/2), i] )
        # print(coil_data[int(len(coil_data[:, i])/2), i])
        z_f = z[idx]
        
        ivol = len(z_f)/len(coil_data[:, i])
        # print(ivol)
        f_vol = np.append(f_vol, ivol)
        i = i + 1
        
    max_idx = np.argmax(np.abs(f_vol))
    print('Maximum volum at [Zshift, radius]=', comb[max_idx])
        
    f_vol_im =  np.ones((len(radius), len(zshift)))
    ir = 0
    ia = 0
    
    while ir< len(radius):
        iz = 0
        while iz < len(zshift):
            f_vol_im [ir, iz] = f_vol[ia]
            # im_std [ir, iz] = std_arr[ia]
            iz = iz +1
            ia = ia +1
    
        ir = ir +1 
    
    return f_vol_im

def angle_data(filename, low_z, high_z, zstep, low_r, high_r, rstep, limt ):
    
    coil_data = np.loadtxt(filename, skiprows = 8)
    l_r = low_r
    h_r = high_r
    # rstep = 2
    l_z = low_z
    h_z = high_z
    # zstep = 5
    
    radius = np.arange(l_r,h_r+rstep,rstep)
    zshift = np.arange(l_z,h_z+zstep, zstep)
    comb = np.array(np.meshgrid(radius, zshift))
    comb = comb.T.reshape(-1, 2)
    
    
    
    z = coil_data[:, 0]
    
    f_vol = []
    
    i = 1
    while i < len(coil_data[0, :]):
        idx = np.where(  np.abs(coil_data[:, i]) < limt )
        z_f = z[idx]
        
        ivol = len(z_f)/len(coil_data[:, i])
        # print(ivol)
        f_vol = np.append(f_vol, ivol)
        i = i + 1
    
    max_idx = np.argmax(np.abs(f_vol))
    print('Maximum volum at [Zshift, radius]=', comb[max_idx])
    
    
    f_vol_im =  np.ones((len(radius), len(zshift)))
    ir = 0
    ia = 0
    
    while ir< len(radius):
        iz = 0
        while iz < len(zshift):
            f_vol_im [ir, iz] = f_vol[ia]
            # im_std [ir, iz] = std_arr[ia]
            iz = iz +1
            ia = ia +1
    
        ir = ir +1 
    
    return f_vol_im


tol_angle = 0.5
tol_enorm = 0.05


angle = angle_data('coils_angle_bounday', 40, 60, 2, -30, 40, 5,  tol_angle)
norm_E = Efield_data('coils_normE_axis', 40, 60, 2, -30, 40, 5, tol_enorm)    

plt.figure()
plt.title(str(tol_angle)+' deg tolerance fiducial volume of coil model')
plt.imshow(np.abs(angle), extent=( 40, 60, 40, -30), aspect='auto')
plt.ylabel('Z Shift [mm]')
plt.xlabel('Radius [mm]')
clb = plt.colorbar()
clb.ax.set_title('Vol ratio')

plt.figure()
plt.title(str(tol_enorm)+'E tolerance Norm E of coil model')
plt.imshow( np.abs(norm_E), extent=( 40, 60, 40, -30), aspect='auto')
plt.ylabel('Z Shift [mm]')
plt.xlabel('Radius [mm]')
clb = plt.colorbar()
clb.ax.set_title('Vol ratio')

plt.figure()
plt.title(str(tol_enorm)+' reweight coil-plane model')
plt.imshow( np.minimum(np.abs(norm_E),np.abs(angle)), extent=( 40, 60, 40, -30), aspect='auto')
plt.ylabel('Z Shift [mm]')
plt.xlabel('Radius [mm]')
clb = plt.colorbar()
clb.ax.set_title('Vol ratio')

# angle_plane = angle_data('plane_angle_bounday', 40, 60, 2, -30, 40, 5,  tol_angle)
# norm_E_plane = Efield_data('plane_normE_axis', 40, 60, 2, -30, 40, 5, tol_enorm)   


# plt.figure()
# plt.title('0.5 deg fiducial volume of plane model')
# plt.imshow(np.abs(angle_plane), extent=( 40, 60, 40, -30), aspect='auto')
# plt.ylabel('Z Shift [mm]')
# plt.xlabel('Radius [mm]')
# clb = plt.colorbar()
# clb.ax.set_title('Vol ratio')

# plt.figure()
# plt.title('99% Norm E of plane model')
# plt.imshow( np.abs(norm_E_plane), extent=( 40, 60, 40, -30), aspect='auto')
# plt.ylabel('Z Shift [mm]')
# plt.xlabel('Radius [mm]')
# clb = plt.colorbar()
# clb.ax.set_title('Vol ratio')

angle_pc = angle_data('plane_coils_angle_bounday', 40, 60, 2, -30, 40, 5,  tol_angle)
norm_E_pc = Efield_data('plane_coils_normE_axis', 40, 60, 2, -30, 40, 5, tol_enorm)   

plt.figure()
plt.title(str(tol_angle)+' deg tolerance fiducial volume of coil-plane model')
plt.imshow(np.abs(angle_pc), extent=( 40, 60, 40, -30), aspect='auto')
plt.ylabel('Z Shift [mm]')
plt.xlabel('Radius [mm]')
clb = plt.colorbar()
clb.ax.set_title('Vol ratio')

plt.figure()
plt.title(str(tol_enorm)+'E tolerance Norm E of coil-plane model')
plt.imshow( np.abs(norm_E_pc), extent=( 40, 60, 40, -30), aspect='auto')
plt.ylabel('Z Shift [mm]')
plt.xlabel('Radius [mm]')
clb = plt.colorbar()
clb.ax.set_title('Vol ratio')

plt.figure(dpi=200)
plt.title(r'0.05|E| and 0.5$^{\circ}$ error tolerance')
plt.imshow( np.minimum(np.abs(norm_E_pc),np.abs(angle_pc)), extent=( 40, 60, 40, -30), aspect='auto')
plt.ylabel('Z Shift [mm]')
plt.xlabel('Radius [mm]')
clb = plt.colorbar()
clb.ax.set_title('Vol ratio')

tol_angle = 0.2
tol_enorm = 0.02

angle_pc = angle_data('plane_coils_angle_bounday', 40, 60, 2, -30, 40, 5,  tol_angle)
norm_E_pc = Efield_data('plane_coils_normE_axis', 40, 60, 2, -30, 40, 5, tol_enorm)   

plt.figure()
plt.title(str(tol_angle)+' deg tolerance fiducial volume of coil-plane model')
plt.imshow(np.abs(angle_pc), extent=( 40, 60, 40, -30), aspect='auto')
plt.ylabel('Z Shift [mm]')
plt.xlabel('Radius [mm]')
clb = plt.colorbar()
clb.ax.set_title('Vol ratio')

plt.figure()
plt.title(str(tol_enorm)+'E tolerance Norm E of coil-plane model')
plt.imshow( np.abs(norm_E_pc), extent=( 40, 60, 40, -30), aspect='auto')
plt.ylabel('Z Shift [mm]')
plt.xlabel('Radius [mm]')
clb = plt.colorbar()
clb.ax.set_title('Vol ratio')

plt.figure(dpi=200)
plt.title(r'0.02|E| and 0.2$^{\circ}$ error tolerance')
plt.imshow( np.minimum(np.abs(norm_E_pc),np.abs(angle_pc)), extent=( 40, 60, 40, -30), aspect='auto')
plt.ylabel('Z Shift [mm]')
plt.xlabel('Radius [mm]')
clb = plt.colorbar()
clb.ax.set_title('Vol ratio')




