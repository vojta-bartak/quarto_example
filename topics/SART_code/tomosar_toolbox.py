# -*- coding: utf-8 -*-
"""
Functions for TomoSAR data processing.
@author: Xiao Liu
E-mail: xiao.liu@mailbox.tu-dresden.de
"""
import os
import pickle
import numpy as np
from tqdm import tqdm
from skimage.transform import resize
from scipy.signal import find_peaks
from math import ceil

#%% tomography
def tomobox(cov,kz_stack,z_vector,outname,tomomethod,overwrite=True):
    if overwrite or not os.path.isfile(outname):
        n_row=np.shape(cov)[0]
        n_col=np.shape(cov)[1]
        n_track=np.shape(cov)[2]    
        n_height=np.shape(z_vector)[0]
    
        tomo=np.zeros((n_row, n_col, n_height), dtype=np.float32)
        
        for n in tqdm(range(n_row)):
            
            for m in range(n_col):
                             
                each_cov = np.squeeze(cov[n,m,:,:])
                each_kz = np.squeeze(kz_stack[n,m,:])               
                         
                # define the steering matrix
                A=np.linalg.multi_dot( [np.matrix(each_kz).T, np.matrix(-1j*z_vector)] )
                A=np.exp(A)
                
                if tomomethod=='capon':#capon
                    inv_temp = np.linalg.inv(each_cov + 1e-2*np.eye(n_track))
                    diag_temp=np.diag(np.linalg.multi_dot( [np.conjugate(np.transpose(A)), inv_temp, A]) )

                    a=1/np.real(diag_temp)
                elif tomomethod=='beamforming':#beamforming
                    diag_temp=np.diag(np.linalg.multi_dot( [np.conjugate(np.transpose(A)), each_cov, A]))

                    a=np.real(diag_temp)/n_track**2
                    
                tomo[n,m,:] = a#power         
                                              
    else:
        with open(outname, 'rb') as f:
            tomo = pickle.load(f)            
        
    return tomo

#%% normalization to 0 and 1
def normalize(layer):
    layer_max = np.amax(layer)
    layer_min = np.amin(layer)
    norm_slice=np.divide((layer - layer_min), (layer_max - layer_min))
    np.seterr(divide='ignore', invalid='ignore')
    return norm_slice

#%% remove residual terrain phase (terrain normalization)
def topo_residual_correction(normalized_stack, kz_stack, z_vector, 
                             hh_tomo_path, terrain_path, overwrite = True):
    
    if os.path.isfile(terrain_path) and overwrite == False:
        terrain = np.load(terrain_path)
    else:
        # use hh tomograms to detect under canopy terrain

        hh_tomo = np.load(hh_tomo_path)        
        terrain = np.zeros((hh_tomo.shape[0], hh_tomo.shape[1]))
        low_bound = 0.3#0.3
        for i in tqdm(range(hh_tomo.shape[0])):
            for j in range(hh_tomo.shape[1]):
                each_tomo = hh_tomo[i,j,:]
                peaks, _ = find_peaks(each_tomo, height=low_bound*np.nanmax(each_tomo)) 
    
                if len(peaks)>0:
                    terrain[i,j] = z_vector[peaks[0]]# height
                    
        del hh_tomo
    
        # save hh_tomo_based undercanopy terrain
        np.save(terrain_path, terrain)
        
    # remove phase casused by difference between dem and dtm (from hh tomosar with capon)
    n_row, n_col, n_track = normalized_stack.shape
    terrain_re = resize(terrain, (n_row, n_col))#.astype('float32')

    for n in tqdm(range(0,n_track)):
        each_kz=np.squeeze(kz_stack[:,:,n])
        each_slc=np.squeeze(normalized_stack[:,:,n])
        normalized_stack[:,:,n]=each_slc*np.exp(1j*each_kz*terrain_re)
        
    return normalized_stack, terrain

#%% down_sampling covariance matrix
# The following codes are based on a lecture from Stefano Tebaldini et.al 
# in PolInSAR 2021 conference

def Filter_Matrix(t_out,t_in,L):

# generates logic filtering matrix

# t_out = output axis in samples
# t_in  = input axis in samples
# L = half a window length in samples
  t_in = t_in.T
  t_out = t_out.T
  T_in, T_out = np.meshgrid(t_in, t_out)
  F = np.abs(T_in-T_out)<(L+1)
  return F

def covmat_downsampling(img_stack, multi_look,  ps_rg, ps_az,
                        coh_matrix_flag = 1):


    n_row, n_col, n_track = img_stack.shape

    #% downsampling because computing capability limit, reference of this code (SARSIM)? better also go to cov function
    # filter and downsample along range
    Lr = ceil(multi_look/2/ps_rg)#ceil:return the nearest and upward integer
    r_out_smpl = np.arange(Lr, n_row-Lr+1, ceil(Lr/2))
    if r_out_smpl.size == 0:
      r_out_smpl = int(np.round(n_row/2))
    Fr = Filter_Matrix(r_out_smpl, np.arange(1,n_row+1,1),Lr)

    # filter and downsample along azimuth
    Lx = ceil(multi_look/2/ps_az)
    x_out_smpl = np.arange(Lx, n_col-Lx+1, ceil(Lx/2))
    if x_out_smpl.size == 0:
      x_out_smpl = int(np.round(n_col/2))
    Fx = Filter_Matrix(x_out_smpl, np.arange(1,n_col+1,1),Lx)

    rows, cols, nTrack = img_stack.shape
    rows_ds = np.shape(Fr)[0]# row number after downsampling
    cols_ds = np.shape(Fx)[0]# column number after downsampling

    cov_matrix = np.ones((rows_ds, cols_ds, nTrack, nTrack), np.complex64)
       
    for i in tqdm(range(0, nTrack)):      
        for j in range(i, nTrack):

            # multiply S1 and complex conjugate of S2
            ms = np.multiply(img_stack[:, :, i], np.conjugate(img_stack[:, :, j]))
            each_cov_matrix = np.linalg.multi_dot([Fr , ms , Fx.conj().T])
            
            cov_matrix[:, :, i, j] = each_cov_matrix
            cov_matrix[:, :, j, i] = np.conjugate(each_cov_matrix)
                 
    if coh_matrix_flag==1:
        for i in tqdm(range(rows_ds)):
            for j in range(cols_ds):
                each=np.squeeze(cov_matrix[ i, j,:, :])
                E_term=np.diag(np.diag(each)**(-1/2))
                cov_matrix[ i, j,:, :] = np.linalg.multi_dot([E_term, each, E_term])

    return cov_matrix, r_out_smpl, x_out_smpl

