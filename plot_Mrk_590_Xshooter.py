#!/usr/bin/env python3
"""
    Author: Freja Amalie Nørby, 27/02/2025

    E-mail: freja.noerby@nbi.ku.dk

"""
import os, sys, time, gc
from datetime import timedelta, datetime
from pathlib import Path
from urllib import request

import numpy as np
import matplotlib.pyplot as plt

from astropy.table import Table
from astropy.units import Unit

from ppxf.ppxf import ppxf, robust_sigma
import ppxf.ppxf_util as util
import ppxf.sps_util as lib
from vorbin.voronoi_2d_binning import voronoi_2d_binning
from plotbin.display_bins import display_bins
from plotbin.display_pixels import display_pixels
from plotbin.plot_velfield import plot_velfield

import functions as f # 

import warnings
warnings.simplefilter("ignore", category=RuntimeWarning)
warnings.simplefilter("ignore", SyntaxWarning)

save = 'output/'
os.makedirs(save, exist_ok=True)

################################################################################################
################################ LOAD DATA
data = Table.read("data/2022_01_27_fluxscaled_spectrum.dat", format="ascii")

lam = data['obs_wl_AA'].data * 10 # for Å
spec = data['flux_cgs'].data 
err = data['err_cgs'].data

z =  0.02647065 # 0.026471 # from ned.ipac.caltech.edu no 19 Ref. https://ned.ipac.caltech.edu/reflookup?refcode=2022ApJS..261....2K
c = 299792.458  # speed of light in km/s

################ CREATE MASK FOR THE THREE SPEC ARMS
#From XSHOOTER Overview https://www.eso.org/sci/facilities/paranal/instruments/xshooter/overview.html
#    UVB, wavelength range 300-559.5 nm
#    VIS, wavelength range 559.5-1024 nm
#    NIR, wavelength range 1024-2480 nm

uv = lam <= 5595
vis = (lam >= 5595) & (lam <= 10240)
nir = lam >= 10240

################################ PLOT IMAGE

plt.figure(figsize=(15,6))
plt.plot(lam[uv], spec[uv], color='b')
plt.plot(lam[vis], spec[vis], color='g')
plt.plot(lam[nir], spec[nir], color='r')

plt.plot(lam, err, color='k')
plt.xlabel('Observed-frame wavelength (nm)', fontsize=16)
plt.ylabel(r"Flux density ($\mathrm{erg\ s^{-1}\ cm^{-2}\ \AA^{-1}}$)", fontsize=16)
plt.xlim(3080, 24810), plt.ylim(-1e-16,6e-15)
plt.savefig(save + 'Mrk590_Xshooter.png')
plt.show()
plt.close()

lam /= (1 + z) # generally simpler to deredshift the spectrum before performing the pPXF fit.

################################################################
em_line_names =            ['[Ne V]', '[Ne V]', '[OII]', '[Ne III]', '[O III]',  'Hbeta',  '[OIII]',       '',  '[OI]',      '', '[NII]', 'Halpha', '[SII]',      '', '[S III]']
em__wavelengths = np.array([3345.821, 3425.881, 3728.82,   3868.760,   4363.21,  4861.33,  4958.911,  5006.84, 6300.30, 6548.05, 6583.46,  6562.82, 6716.47, 6732.67, 9531.100 ])

abs_labels = ['K', 'H', ' G', 'Mg', 'Na', ' ', '', ' CaII']
abs_wavelengths = np.array([3934.777, 3969.588, 4305.61, 5176.7, 5895.6, 8500.36, 8544.44, 8664.52])
################################################################

################################ PLOT IMAGE WITH ABSORPTION LINES AND EMISSION LINES
fig, (ax1, ax2) = plt.subplots(figsize=(15,12), nrows=2, ncols=1)

for ax in [ax1, ax2]:
    
    for i, (wl, label) in enumerate(zip(em__wavelengths, em_line_names)):
        ax.axvline(wl, color='C1', linestyle='dotted')
        
        # Conditional y positioning for even and odd index
        y_pos = 4.5e-15 if i % 2 == 0 else 5e-15
        ax.annotate(label, xy=(wl, 1e-15), xytext=(wl - 40, y_pos),
                    color='orangered', textcoords='data', 
                    horizontalalignment='center', verticalalignment='bottom',
                    fontsize=10, rotation=90)
     
    for i, (wl, label) in enumerate(zip(abs_wavelengths, abs_labels)):
        #ax.axvline(wl, color='r', linestyle='--')
        ax.fill_between([wl - 9, wl + 9], 0, 2, color='gray', alpha=0.5)
        x_pos = wl-40 if i % 2 == 0 else wl + 40
        ax.annotate(label, xy=(wl, 1e-15), xytext=(x_pos, 0.1e-15),
                    color='k', textcoords='data', 
                    horizontalalignment='center', verticalalignment='bottom',
                    fontsize=10, rotation=90)

ax1.plot(lam[uv], spec[uv], color='b')
ax1.set_title("UV spectum", fontsize=18)
ax1.set_xlim(min(lam[uv]-20), max(lam[uv])+20), ax1.set_ylim(-1e-16,6e-15)

ax2.plot(lam[vis], spec[vis], color='g')
ax2.set_xlabel('Redshift corected wavelength (Å)', fontsize=16)
ax2.set_ylabel(r"Flux density ($\mathrm{erg\ s^{-1}\ cm^{-2}\ \AA^{-1}}$)", fontsize=16)
ax2.set_title("VIS sepctrum", fontsize=18)
ax2.set_xlim(min(lam[vis]-20), max(lam[vis])+20), ax2.set_ylim(-1e-16,6e-15)

# Show plot
plt.tight_layout()
plt.savefig(save + 'Abs_Em_lines.png')
plt.show()