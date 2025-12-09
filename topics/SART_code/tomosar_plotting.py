# -*- coding: utf-8 -*-
"""
Functions for plotting TomoSAR results.
@author: Xiao Liu
E-mail: xiao.liu@mailbox.tu-dresden.de
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from tqdm import tqdm
from skimage.transform import resize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.signal import find_peaks

#%% define a function to adjust the colorbar
def color_bar_fit (ax, im, size="5%"):
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size=size, pad=0.05)
    return plt.colorbar(im, cax=cax)

#%% plot InSAR phase with/without flat-earth and topography phase
def insar_quick_look(slc_1_id, slc_2_id, slc_stack, normalized_stack):
    slc_1 = slc_stack[:,:,slc_1_id].squeeze()
    slc_2 = slc_stack[:,:,slc_2_id].squeeze()

    norm_slc_1 = normalized_stack[:,:,slc_1_id].squeeze()
    norm_slc_2 = normalized_stack[:,:,slc_2_id].squeeze()
    
    insar_before = slc_1 * np.conj(slc_2)
    insar_after = norm_slc_1 * np.conj(norm_slc_2)
    
    plt.figure(dpi = 300)
    cmap = 'jet'
    
    ax = plt.subplot(1, 2,1)
    im = ax.imshow( np.angle(insar_before), cmap = cmap )
    color_bar_fit (ax, im)
    ax.set_title('With topographical phase')
    
    ax = plt.subplot(1, 2, 2)
    im = ax.imshow( np.angle(insar_after), cmap = cmap)
    color_bar_fit (ax, im)
    ax.set_title('Without topographical phase')
    
    plt.tight_layout()
    
#%% plot covariance matrix
def cov_mat_plot(cov_mat, img_path):
    n_track = int( cov_mat.shape[2] )
    
    fig = plt.figure(dpi = 300, figsize = (18, 18))
    gs = fig.add_gridspec(n_track, n_track, hspace=0, wspace=0)
    axs = gs.subplots() #sharex='col', sharey='row'
    for r_id in tqdm( range(0,n_track) ):
        axs[r_id, r_id].text(0.5, 0.5, 'Track {}'.format(r_id + 1), horizontalalignment='center',
                             verticalalignment='center', transform=axs[r_id, r_id].transAxes,
                             fontsize=30)
        axs[r_id, r_id].axis("off")
    
        for c_id in range(r_id + 1, n_track):
            img_1 = axs[r_id, c_id].imshow(np.abs(cov_mat[:,:,r_id, c_id].squeeze()),
                                           cmap = 'Greys_r')#cmap = 'viridis')#'Greys_r'
            img_2 = axs[c_id, r_id].imshow(np.angle(cov_mat[:,:,c_id, r_id].squeeze()), cmap = 'jet')
            axs[r_id, c_id].axis("off")
            axs[c_id, r_id].axis("off")
    
    # adding the colorbar
    cax_1 = fig.add_axes([1, 0.49, 0.03, 0.5]) #left, bottom, width, height 
    cax_2 = fig.add_axes([-0.1, 0.01, 0.03, 0.5]) #left, bottom, width, height 
    
    bar_1 = fig.colorbar(img_1, ax=axs.ravel().tolist(), cax = cax_1)
    bar_2 = fig.colorbar(img_2, ax=axs.ravel().tolist(), cax = cax_2)
    
    # change the font size of the colorbar labels
    bar_1.ax.tick_params(labelsize=30)
    bar_2.ax.tick_params(labelsize=30)
    
    bar_1.ax.set_ylabel('Coherence_amplitude',fontsize=30)
    bar_2.ax.set_ylabel('Coherence_phase',fontsize=30)
    
    plt.tight_layout() 
    
    plt.savefig(img_path, dpi=300, bbox_inches='tight')
  
#%% plot tomosar results
def tomo_plot(rg, az, slc_stack, tomo_norm, height, pol, img_path, 
              lidar_rh = np.zeros((1,1))
              ):
    
    # main plot setup
    # set up the subplot layout
    fig = plt.figure(dpi = 300, figsize=(12,10))
    fontsize = 12

    # the horizontal slice plot
    ax1 = fig.add_subplot(221)
    # the vertical profile plot
    ax2 = fig.add_subplot(222)
    # the range slice plot
    ax3 = fig.add_subplot(413)
    # the azimuth slice plot
    ax4 = fig.add_subplot(414)
    plt.subplots_adjust(left=0.1, right=0.2, top=0.3, bottom=0.2)
    
    # set up the plots for range and azimuth slices
    ax3.set_xlim(0, tomo_norm.shape[1]-1)
    ax4.set_xlim(0, tomo_norm.shape[0]-1)
    ax3.set_ylim(0, height * 2)
    ax4.set_ylim(0, height * 2)
    
    # set up the vertical profile plot
    ax2.set_ylabel('Height [m]', fontsize=fontsize)
    ax2.set_xlabel('Reflectivity', fontsize=fontsize)
    ax2.set_title('Vertical point profiles, {}'.format(pol.upper()), fontsize=fontsize)
    ax2.set_xlim(0, 1)
    ax2.set_ylim(-height, height)
    
    # add a cross-hair to the horizontal slice plot
    ax1.axhline(az, linewidth=4, color='w')
    ax1.axvline(rg, linewidth=4, color='w')
    ax3.axvline(rg, linewidth=4, color='w')
    ax4.axvline(az, linewidth=4, color='w')
    
    subset_vertical = tomo_norm[az, rg, :]
    subset_range = np.fliplr(tomo_norm[az, :, :])
    subset_azimuth = np.fliplr(tomo_norm[:, rg, :])

    subset_layer = np.abs( slc_stack[:,:,0].squeeze() )# slc intensity
    subset_layer = resize(subset_layer, (tomo_norm.shape[0], tomo_norm.shape[1]))#.astype('float32')
    
    # plot the vertical profile
    label = 'rg: {}; az: {}'.format(rg, az)
    ax2.plot(subset_vertical, range(-height, height + 1), label=label)
    ax2.legend(loc=0, prop={'size': 12}, markerscale=1)
    
    # plot the range slice
    im3 = ax3.imshow(np.rot90(subset_range, 1), origin='lower', cmap='jet', aspect='auto')
    ax3.set_title('Range slice at azimuth line {}, {}'.format(az, pol.upper()), fontsize=fontsize)
    
    # plot the azimuth slice
    im4 = ax4.imshow(np.rot90(subset_azimuth, 1), origin='lower', cmap='jet', aspect='auto')
    ax4.set_title('Azimuth slice at range line {}, {}'.format(rg, pol.upper()), fontsize=fontsize)
    
    # change yaxis of ax3 and ax4 from bins to meters
    z_vector = np.arange(-height, height + 1, 1)
    ticks = np.linspace(0, 2*height,6)
    ticklabels = ["{:.0f}".format(z_vector[int(i)]) for i in ticks]
    
    ax3.set_yticks(ticks)
    ax3.set_yticklabels(ticklabels)
    ax4.set_yticks(ticks)
    ax4.set_yticklabels(ticklabels)
    ax3.set_ylabel('Height [m]', fontsize=fontsize)
    ax4.set_ylabel('Height [m]', fontsize=fontsize)
    
    im1 = ax1.imshow(subset_layer, origin='upper', cmap='viridis',#cmap='Greys_r',
                    vmin = np.percentile(subset_layer, 10),
                    vmax = np.percentile(subset_layer, 90))# from top to down
    
    ax1.set_xlabel('Range', fontsize=fontsize)
    ax1.set_ylabel('Azimuth', fontsize=fontsize)
    ax1.set_title('Intensity of the reference SLC, {}'.format(pol.upper()), fontsize=fontsize)
    
    color_bar_fit (ax1, im1)
    color_bar_fit (ax3, im3, size='2%')
    color_bar_fit (ax4, im4, size='2%')
    
    if lidar_rh.size > 1:
        lidar_rh = resize(lidar_rh, (tomo_norm.shape[0], tomo_norm.shape[1]))
        lidar_rh_range = lidar_rh[az, :].squeeze() + height
        lidar_rh_azimuth = lidar_rh[:, rg].squeeze() + height
        each_lidar_rh = lidar_rh[az, rg]
        
        ax2.axhline(each_lidar_rh, color='k', label = 'lidar' )#label = 'lidar'

        ax3.plot(lidar_rh_range,linewidth=4, color='k', label = 'lidar' )
        ax4.plot(lidar_rh_azimuth,linewidth=4, color='k', label = 'lidar')

        ax2.legend(loc=0, prop={'size': 12}, markerscale=1)        
        ax3.legend(loc=0, prop={'size': 12}, markerscale=1)
        ax4.legend(loc=0, prop={'size': 12}, markerscale=1)

    plt.tight_layout(pad=1.0, w_pad=0.1, h_pad=0.1)
    
    plt.savefig(img_path, dpi=300, bbox_inches='tight')
    
#%% plot aggregated tomosar profiles in different forest height and agb
def grouped_tomosar_profiles (tomo, lvis_rh, lvis_agb, z_vector, img_path):

    fig = plt.figure(dpi = 300, figsize=(10,6))
    fontsize = 12
    
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
        
    # grouped by lvis height
    h_step = 10
    n_step = 5
    nz=z_vector.shape[0]
     
    for h_id in range(n_step):
        h_start = h_id * h_step
        h_end = (h_id+1) * h_step
        
        mask = np.where((lvis_rh >= h_start) & (lvis_rh < h_end) )
    
        tomo_median = np.zeros_like(z_vector, dtype=np.float32)
        for layer in range(nz):
            each_layer = np.squeeze( tomo[:,:,layer] )[mask]
    
            tomo_median[layer] = np.median(each_layer[each_layer>0])
                
        label = '{}-{} m'.format(h_start, h_end)
        ax1.plot( tomo_median,  z_vector, label = label) 
     
    ax1.set_ylim(0,60)
    ax1.grid(alpha=0.5)
    ax1.legend(bbox_to_anchor=(1.04, 0.25), loc="center left", borderaxespad=0,fontsize=fontsize)
    ax1.set_xlabel('Reflectivity',fontsize=fontsize)
    ax1.set_ylabel('Height (m)',fontsize=fontsize)
    ax1.tick_params(axis='both', labelsize=fontsize)
    ax1.set_title('Grouped by lidar forest height')
    
    # grouped by lvis agb
    agb_step = 100
    n_step = 5
    
    for agb_id in range(n_step):
        agb_start = agb_id * agb_step
        agb_end = (agb_id+1) * agb_step
            
        mask = np.where((lvis_agb >= agb_start) & (lvis_agb < agb_end) )
    
        tomo_median = np.zeros_like(z_vector, dtype=np.float32)
        for layer in range(nz):
            each_layer = np.squeeze( tomo[:,:,layer] )[mask]
            tomo_median[layer] = np.median(each_layer[each_layer>0])
                
        label = '{}-{} Mg/ha'.format(agb_start, agb_end)
        ax2.plot( tomo_median,  z_vector, label = label) 
     
    ax2.set_ylim(0,60)
    ax2.grid(alpha=0.5)
    ax2.legend(bbox_to_anchor=(1.04, 0.25), loc="center left", borderaxespad=0,fontsize=fontsize)
    ax2.set_xlabel('Reflectivity',fontsize=fontsize)
    ax2.set_ylabel('Height (m)',fontsize=fontsize)
    ax2.tick_params(axis='both', labelsize=fontsize)
    ax2.set_title('Grouped by lidar AGB')
  
    plt.tight_layout()
    plt.savefig(img_path, dpi=300, bbox_inches='tight')

#%% compare tomosar phase centre with lidar height
def tomosar_phase_centre(tomo, lvis_rh, z_vector, img_path):
    low_bound = 0.3#0.3
    
    pha_c = np.zeros((tomo.shape[0], tomo.shape[1]))
    
    for i in tqdm(range(tomo.shape[0])):
        for j in range(tomo.shape[1]):
            each_tomo = tomo[i,j,:]
            peaks, _ = find_peaks(each_tomo, height=low_bound*np.nanmax(each_tomo)) 
    
            if len(peaks)>0:
                pha_c[i,j] = z_vector[peaks[-1]]# height
    
    
    mask = np.where((lvis_rh >= 0) & (lvis_rh < 60) )
        
    r = np.corrcoef(lvis_rh[mask], pha_c[mask])[0,1]
    
    res = stats.pearsonr(lvis_rh[mask], pha_c[mask])
    r = res.statistic
    pvalue = res.pvalue

    r = round(r, 2)
    
    fig = plt.figure(dpi = 300, figsize=(12,4))
    
    ax1 = fig.add_subplot(131)#agb
    ax2 = fig.add_subplot(132)#tomo_h
    ax3 = fig.add_subplot(133)#agb vs tomo_h, scatter 
    
    im1 = ax1.imshow(lvis_rh, vmin = 0, vmax = 60)
    ax1.set_title('LVIS RH100 (m)')
    color_bar_fit (ax1, im1)
    
    im2 = ax2.imshow(pha_c, vmin = 0, vmax = 60)
    ax2.set_title('TomoSAR phase centre (m)')
    color_bar_fit (ax2, im2)
    
    im2 = ax3.scatter(lvis_rh[mask], pha_c[mask], s=1)
    if pvalue > 0.01:
        print(ax3.set_title('r={}, p>0.01'.format(r)))
    else:
        print(ax3.set_title('r={}, p<=0.01'.format(r)))    
            
    ax3.axline((0,0), slope=1, color='black')
    ax3.set_xlim(0,60)
    ax3.set_ylim(0,60)
    ax3.set_xlabel('LVIS RH100 (m)')
    ax3.set_ylabel('TomoSAR phase centre (m)')
    
    plt.tight_layout()
    plt.savefig(img_path, dpi=300, bbox_inches='tight')
    
#%% compare agb with tomosar reflectivity at different height layers
def tomosar_layerd_reflectivity (tomo, lvis_agb, min_agb, height, img_path):

    mask = np.where(lvis_agb >= min_agb)

    for h in range(0, 60, 15):
        
        subset_layer = tomo[:,:,height+h].squeeze()
                
        res = stats.pearsonr(lvis_agb[mask], subset_layer[mask])
        r = res.statistic
        pvalue = res.pvalue
        
        r = round(r, 2)
    
        fig = plt.figure(dpi = 300, figsize=(12,4))
    
        ax1 = fig.add_subplot(131)#agb
        ax2 = fig.add_subplot(132)#tomo_h
        ax3 = fig.add_subplot(133)#agb vs tomo_h, scatter 
        
        im1 = ax1.imshow(lvis_agb, vmin = min_agb, vmax = 600)
        ax1.set_title('LVIS AGB (Mg/ha)')
        color_bar_fit (ax1, im1)
    
        im2 = ax2.imshow(subset_layer)
        ax2.set_title('TomoSAR reflectivity at {} m'.format(h))
        color_bar_fit (ax2, im2)
        
        im2 = ax3.scatter(lvis_agb[mask], subset_layer[mask], s=1)
        if pvalue > 0.01:
            print(ax3.set_title('r={}, p>0.01'.format(r)))
        else:
            print(ax3.set_title('r={}, p<=0.01'.format(r)))
            
        ax3.set_xlabel('AGB (Mg/ha)')
        ax3.set_ylabel('Reflectivity')
        
        plt.tight_layout()
        
        img_path = img_path + '_{}m.png'.format(h)
        plt.savefig(img_path, dpi=300, bbox_inches='tight')

#%% look at the data: slc intensity, phase, kz, topographical phase        
def quick_look (slc_id, slc_stack, kz_stack, phase_stack, img_path):
    
    slc_subset = slc_stack[:,:,slc_id].squeeze()
    kz_subset = kz_stack[:,:,slc_id].squeeze()
    phase_subset = phase_stack[:,:,slc_id].squeeze()
      
    slc_amp = np.abs(slc_subset)
    slc_pha = np.angle(slc_subset)
    
    fontsize = 12
    plt.figure(dpi = 300, figsize = (8,8))
    
    ax = plt.subplot(221)
    im = plt.imshow(slc_amp, cmap = 'Greys_r', 
                    vmin = np.percentile(slc_amp, 10),
                    vmax = np.percentile(slc_amp, 90))
    ax.set_title('Track {}, amplitude'.format(slc_id), fontsize=fontsize)
    color_bar_fit (ax, im)
    
    ax = plt.subplot(222)
    im = plt.imshow(slc_pha, cmap = 'jet')
    ax.set_title('Track {}, phase (rad)'.format(slc_id), fontsize=fontsize)
    color_bar_fit (ax, im)
    
    ax = plt.subplot(223)
    im = plt.imshow(kz_subset, cmap = 'jet')
    ax.set_title('Track {}, vertical wavenumber (rad/m)'.format(slc_id), fontsize=fontsize)
    color_bar_fit (ax, im)
    
    ax = plt.subplot(224)
    im = plt.imshow(phase_subset, cmap = 'jet')
    ax.set_title('Track {}, phase correction term (rad)'.format(slc_id), fontsize=fontsize)
    color_bar_fit (ax, im)
    
    plt.tight_layout()
    
    plt.savefig(img_path, dpi=300, bbox_inches='tight')
