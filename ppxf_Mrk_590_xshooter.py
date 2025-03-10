#!/usr/bin/env python3
"""
    Author: Freja Amalie Nørby, 27/02/2025

    E-mail: freja.noerby@nbi.ku.dk
    
    Performing pPXF fitting procedure developed by Cappelari
    "Penalized Pixel-Fitting method by Cappellari & Emsellem (2004)
    as upgraded in Cappellari (2017, 2023)".

    See

    https://ui.adsabs.harvard.edu/abs/2004PASP..116..138C
    https://ui.adsabs.harvard.edu/abs/2017MNRAS.466..798C
    https://ui.adsabs.harvard.edu/abs/2023MNRAS.526.3273C

"""
import os, sys, time
from datetime import timedelta, datetime
from pathlib import Path
from urllib import request

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from astropy.table import Table
from astropy.units import Unit

from ppxf.ppxf import ppxf, robust_sigma
import ppxf.ppxf_util as util
import ppxf.sps_util as lib
from vorbin.voronoi_2d_binning import voronoi_2d_binning
from plotbin.display_bins import display_bins
from plotbin.display_pixels import display_pixels
from plotbin.plot_velfield import plot_velfield

import pickle
import warnings
warnings.simplefilter("ignore", category=RuntimeWarning)
warnings.simplefilter("ignore", SyntaxWarning)

################################################################################################
print('#'*64 + '\n' + f'ppxf for X-shooter Mrk590 (2022)             {
datetime.now().strftime("%Y-%m-%d %H:%M:%S")} '+ '\n' + '#'*64 +'\n')

z = 0.02647065 # from ned ipac no 19 Ref.https://ned.ipac.caltech.edu/reflookup?refcode=2022ApJS..261....2K
c = 299792.458  # speed of light in km/s

################################ LOAD DATA
################
data = Table.read("data/corr_2022_01_27_combined_spectrum.dat", format="ascii")

lam = data['obs_wl_nm'].data * 10 # for Å
spec = data['flux_cgs'].data 
err = data['err_cgs'].data

################
data8 = Table.read("data/corr_2022_01_08_combined_spectrum.dat", format="ascii")

lam8 = data8['obs_wl_nm'].data * 10 # for Å
spec8 = data8['flux_cgs'].data 
err8 = data8['err_cgs'].data

################################

################ CREATE MASK FOR THE THREE SPEC ARMS
#From XSHOOTER Overview https://www.eso.org/sci/facilities/paranal/instruments/xshooter/overview.html
#    UVB, wavelength range 300-559.5 nm
#    VIS, wavelength range 559.5-1024 nm
#    NIR, wavelength range 1024-2480 nm

uv = (lam >= 3000) & (lam <= 5595) 
vis = (lam >= 5595) & (lam <= 10240)
nir = lam >= 10240

#### From XSHOOTER observations headers under Program ID 108.229T.001
#      UVB   VIS   NIR
R = [ 3200, 5000, 5600]

#### DEREDSHIFT THE SPECTRUM
# simpler before the pPXF fit
lam /= (1 + z)
lam8 /= (1 + z)

overlap_uv = 5585/(1+z)
overlap_nir = 10240/(1+z) 

ca_trip = (lam >= 8000) & (lam <= 9000) 
uv_vis = (lam >= 3800) & (lam <= 9000) 
agn = lam <= 3900

nir = lam >= 9000

#### SCALE VIS SPECTRUM
#visual inspection loss of light roughly 10% 
adjust = np.ones_like(spec)
adjust[vis]*=1.1
#spec *= adjust
#spec8 *= adjust
################

################################################################
em_line_names =            ['[Ne V]', '[Ne V]', '[OII]', '[Ne III]', '[S II]', '[O III]',  'He II',  'Hbeta',  '[OIII]',       '', '[N I]', '[Fe VII]', '[N II]',   'He I',  '[Fe VII]', '[OI]',   '[O I]',      '', '[NII]', 'Halpha', '[SII]',      '', '[O II]', '[S III]'] #'[Ar III]', '[Ni II]', '[Ni II]', 
em__wavelengths = np.array([3345.821, 3425.881, 3728.82,   3868.760, 4068.600,   4363.21, 4685.710,  4861.33,  4958.911,  5006.84, 5200.257,  5720.700,  5754.59, 5875.624,    6087.000, 6300.30, 6363.776, 6548.05, 6583.46,  6562.82, 6716.47, 6732.67, 7319.990, 9531.100 ])#  7135.790,  7377.830,  7411.160, 

balmer_labels =           ['H10',     'H9',     'H8',     'Hε',     'Hδ', '    Hγ',     'Hβ', 'Halpha']
balmer_wave = np.array([3797.904, 3835.391, 3889.064, 3970.079, 4101.742, 4340.471, 4861.333, 6562.819])

abs_labels = ['K', 'H', ' G', 'Mg', 'Na', ' ', '', ' CaII']
abs_wavelengths = np.array([3934.777, 3969.588, 4305.61, 5176.7, 5895.6, 8500.36, 8544.44, 8664.52])
#################################################################

################################################################
### Function to iteratively clip the outliers ##################
def clip_outliers(galaxy, bestfit, mask):
    """
    Repeat the fit after clipping bins deviants more than 3*sigma in
    relative error until the bad bins don't change any more. This 
    function uses eq.(34) of Cappellari (2023);
    https://ui.adsabs.harvard.edu/abs/2023MNRAS.526.3273C
    """
    while True:
        scale = galaxy[mask] @ bestfit[mask]/np.sum(bestfit[mask]**2)
        resid = scale*bestfit[mask] - galaxy[mask]
        err = robust_sigma(resid, zero=1)
        ok_old = mask
        mask = np.abs(bestfit - galaxy) < 2.5 * err
        if np.array_equal(mask, ok_old):
            break
            
    return mask

################################################################
### Function to fit the stellar kinematics #####################
def ppxf_fit_and_clean(templates, galaxy, noise, velscale, start, 
                       mask0, lam, lam_temp, apol, mpol, plot=True, quiet=False):
    """
    Simple pPXF wrapper. It performs two pPXF fits: the first one 
    estimates the scatter in the spectrum and identifies the outlier 
    pixels. The second uses the mask obtained from the first fit to 
    exclude the outliers. 
    
    The general approach used in this function is described in 
    Sec.6.5 of Cappellari (2023):
    https://ui.adsabs.harvard.edu/abs/2023MNRAS.526.3273C
    """
    mask = mask0.copy()
    pp = ppxf(templates, galaxy, noise, velscale, start,
              moments=2,  degree=-1, mdegree=-1, lam=lam, lam_temp=lam_temp,
              mask=mask, quiet=quiet)
    if plot:
        #plt.figure(figsize=(8, 8))
        plt.subplot(211)
        pp.plot(color=['black', 'yellowgreen', 'turquoise', 'crimson'])
        plt.title(f"Initial pPXF fit before outliers removal", fontsize=14)
        plt.xlabel('')        
    
    # Add clipped pixels to the original masked emission lines regions and repeat the fit
    mask &= clip_outliers(galaxy, pp.bestfit, mask)
   
    pp = ppxf(templates, galaxy, noise, velscale, start,
              moments=4,  degree=apol, mdegree=mpol, lam=lam, lam_temp=lam_temp,
              mask=mask, quiet=quiet)
    
    pp.optimal_template = templates.reshape(templates.shape[0], -1) @ pp.weights
    
    resid = (pp.galaxy - pp.bestfit)[pp.goodpixels]
    pp.sn = np.nanmedian(pp.galaxy[pp.goodpixels])/robust_sigma(resid)

    if plot:
        plt.subplot(212)
        plt.gca().set_facecolor('white')
        pp.plot()
        plt.title(f'pPXF fit after; $\\sigma$={pp.sol[1]:.0f} km/s; S/N={pp.sn:.1f}', fontsize=14)
    return pp

################################################################
### Function to show the SPS model and templates ###############
def show_sps(title, self, lib, lib_name):
    print('#'*32 + f' SPS templates: {lib_name}')
    plt.figure(figsize=(14,10))
        
    cmap=plt.get_cmap('jet', lib.templates_reshape.shape[1])
    offset = np.array([7,6,5,4,3,2,1,0])*0.5

    idx = np.where(self.weights>0)[0]
    w = lib.age_grid.shape[0]
    n = []
    nn = []
    for i in idx:
        n.append(i % w)  # Store the value of i % 25 in the list n
        if i<=w:nn.append(0)
        else: nn.append(int(i / w) - 1 )
            
    ################
    ax1 = plt.subplot(311)
    plt.title(f'{title} SPS: {lib_name}', fontsize=20)
    if self.goodpixels[0]>=2000: plt.axvspan(3020, 4000, facecolor = 'lightcoral', alpha=0.5, label='UV data not included in fit')
    self.plot(color=['black', 'limegreen', 'darkgreen', 'red'])
    plt.axvspan(overlap_uv-80, overlap_uv+20, facecolor = 'cyan', alpha=0.5) # Overlap between UV and VIS data
    plt.axvline(overlap_uv, color = 'aqua')
    ax1.minorticks_on()
    plt.ylim(0, 3.5e-15)
    
    ################
    ax2 = plt.subplot(312)
    plt.axvspan(self.lam[0], self.lam[-1], facecolor = 'lightgrey', alpha=0.5, label='Data range') # data range
    if self.goodpixels[0]>=2000: plt.axvspan(3020, 4000, facecolor = 'darkgrey', label='UV data not included in fit')
    plt.plot(self.lam_optimal_template, self.optimal_template, label='Combined pPXF Model \n(Sum of Weighted Templates)', color='black', linewidth=1)
    plt.xlim(1600,13000)
    plt.xticks([])
    plt.legend()
    
    if self.apoly is not None:
        # Create first inset plot (Apoly)
        axins1 = inset_axes(ax2, width=0.5, height=0.5, borderpad=4, bbox_to_anchor=(310, 440))
        axins1.plot(self.apoly, color='C1')
        plt.xticks([])  # Hide the x-axis ticks
        axins1.set_title("Apoly")
    
    if self.mpoly is not None:
        # Create second inset plot (Mpoly)
        axins2 = inset_axes(ax2, width=0.5, height=0.5, borderpad=4, bbox_to_anchor=(310, 530))
        axins2.plot(self.mpoly, color='C2')
        plt.xticks([])  # Hide the x-axis ticks
        axins2.set_title("Mpoly")

    ################
    ax3 = plt.subplot(313) 
    plt.axvspan(self.lam[0], self.lam[-1], facecolor = 'lightgrey', alpha=0.5) # data range
    if self.goodpixels[0]>=2000: plt.axvspan(3020, 4000, facecolor = 'darkgrey') # UV data not included in fit
    for i in range(len(idx)):
        plt.plot(lib.lam_temp_full, lib.templates_reshape[:,idx[i]] + offset[i], color=cmap(idx[i]), linewidth=1.5, 
                 label=f'Star: {idx[i]} [Age = {lib.age_grid[n[i], nn[i]]} Gyr, [M/H] = {lib.metal_grid[n[i], nn[i]]}]')
        print(f'Stellar template: {idx[i]}; Weight = {self.weights[idx[i]]};  [Age = {lib.age_grid[n[i],nn[i]]} Gyr, [M/H] = {lib.metal_grid[n[i], nn[i]]}]')
    
    plt.xlabel(r"$\lambda_{\rm rest}$ ($Ångstrøm$)")
    plt.xlim(1600,13000)
    ax3.minorticks_on()
    plt.ylim(1,12)
    plt.legend(fontsize=8)

    ###########  [left, bottom, width, height]
    ax1.set_position([0.1, 0.6, 0.8, 0.33]) 
    ax2.set_position([0.1, 0.3, 0.8, 0.21]) 
    ax3.set_position([0.1, 0.08, 0.8, 0.21]) 

    print('\n')
    return idx, n, nn
    
################################################################
### Function to asve as ASCII file  ############################
import numpy as np
def save_as_ascii(date, self, lib, idx):
    # Open a file to write the output
    with open(f"output/{date}_ppxf_bestfit_ascii.txt", "w") as f:
        f.write("mask   obs_wl_AA    flux_cgs                err_cgs                 bestfit_model         \n")
        f.write("------ -----------  ----------------------  ----------------------  ----------------------\n")
        
        # Write the data in the specified format
        for ma, wl, flux, err, best in zip(mask_combined, self.lam, self.galaxy, self.noise, self.bestfit):
            f.write(f"{ma} {wl:12.2f}  {flux:22.16e}  {err:22.16e}  {best:22.16e}\n")
    
    with open(f"output/{date}_ppxf_templates_ascii.txt", "w") as f:
        # Generate the header dynamically based on the length of idx
        weights_str = 'Weights:\n'
        weights_str += "".join([f"SPS [{idx[i]}]: {self.weights[idx[i]]:22.16e}\n" for i in range(len(idx))]) + "\n"
        f.write(weights_str)
        
        header = f"temp_lam_AA   optimal_model           lib_lam      "
        header += "   ".join([f"SPS [{idx[i]}]             " for i in range(len(idx))]) + " \n"
        f.write(header)
        f.write("------------  ----------------------  -----------  " + "----------------------  "*len(idx) +"\n")
    
        # Write the data in the specified format
        for wl, optemp, lam, *sps_values in zip(self.lam_optimal_template, self.optimal_template, lib.lam_temp_full, *[lib.templates_reshape[:, i] for i in idx]):
            # Join SPS values dynamically based on the length of idx
            sps_values_str = "  ".join([f"{sps:22.16e}" for sps in sps_values])
            f.write(f"{wl:12.2f}  {optemp:22.16e} {lam:12.2f}  {sps_values_str}\n")
    print(f"Saved at output/{date}_ppxf_templates_ascii.txt\n")
    
    
################################################################ 
def runppxf(date, llam, ggalaxy, nnoise, degree, mdegree, sps_name, R, norm_range=[7500, 8400], AGNmask=True, plot=True, quiet=False, save=False):
    """
    Performs two pPXF fits of the galaxy spectrum using stellar population synthesis (SPS) templates.

    Parameters:
    ----------
    date : str
        The date or other identifyer as string.
    llam : array-like
        The wavelength grid of the galaxy spectrum (in Angstroms).
    ggalaxy : array-like
        The galaxy spectrum (flux values corresponding to `llam`).
    nnoise : array-like
        The noise associated with each pixel of the galaxy spectrum (`ggalaxy`).
    sps_name : str
        The name of the stellar population synthesis (SPS) model to be used.
    R : float
        The resolving power of the instrument (e.g., Xshooter).
    norm_range : list, optional, default=[7500, 8400]
        The wavelength range used for normalization, defaulting to [7500, 8400] Å.
        
    Output Parameters:
    -------
    pp : object
        The result of the pPXF fitting procedure, containing the best-fit parameters and other fitting information.
    sps_name : str
        The name of the SPS model used.
    sps : object
        The loaded SPS object corresponding to the specified stellar population synthesis templates.
    velscale : float
        The velocity scale per spectral pixel, calculated from the instrumental dispersion.
        
    Notes:
    ------
    - The function logarithmically rebins the galaxy spectrum and noise to match the instrumental velocity scale.
    - It also normalizes the galaxy spectrum and applies masks for narrow and broad emission lines, noise regions, 
      and overlap between UV and VIS data.
      
    """
    lam_range = [min(llam), max(llam)] 
    if not quiet: print('\n' + '#'*24 + '      ' + sps_name +' SPS templates ['+ f'{min(llam):.0f},{max(llam):.0f}' +'] Å') 
    
    ################ Logarithmically rebin the spectrum and error to a 
    # velocity scale per spectral pixel equal to the instrumental dispersion 
    # of Xshooter. This ensures Nyquist sampling without loss of information.
    
    sigma_inst = c/(R*np.sqrt(4*np.log(4))) 
    if not quiet: print( f"Instrumental dispersion: {sigma_inst:.0f} km/s\n") 
    
    noise, ln_lam, velscale  = util.log_rebin([min(llam),max(llam)], nnoise, velscale=sigma_inst)
    galaxy, ln_lam, velscale  = util.log_rebin([min(llam),max(llam)], ggalaxy, velscale=sigma_inst)
    
    #galaxy /= np.median(galaxy)  # Normalize spectrum to avoid numerical issues
    #noise = np.full_like(galaxy, 0.01) # Use if you want to assume constant noise per pixel; S=1, N=0.01 --> S/N = 100

    ################################ Setup stellar templates
    ppxf_dir = Path('ppxf') # path.dirname(path.realpath(lib.__file__))
    basename = f"spectra_{sps_name}_9.0.npz"
    filename = ppxf_dir / 'sps_models' / basename
    #FWHM_temp = 2.51   # Resolution of E-MILES templates in the fitted range

    libs = lib.sps_lib(filename, velscale, norm_range= norm_range)
    sps = libs
    lam_range_temp = np.exp(libs.ln_lam_temp[[0, -1]])

    ################################ Setup emission line masks width in km/s
    mask_Narrow  = util.determine_mask(ln_lam, lam_range_temp, redshift=0, width=700, custom_lines=em__wavelengths)                 # Emission lines marked from plot_Mrk_590_Xshooter
    mask_Balmer  = util.determine_mask(ln_lam, lam_range_temp, redshift=0, width=1100, custom_lines=balmer_wave)                    # Balmer emission lines
    mask_Broad   = util.determine_mask(ln_lam, lam_range_temp, redshift=0, width=9200, custom_lines=[6548.03, 6583.41, 6562.80])    # [NII]_d & Halpha
    mask_Bbroad  = util.determine_mask(ln_lam, lam_range_temp, redshift=0, width=6000, custom_lines=[4861.33,  4958.911,  5006.84]) # Hbeta & [OIII]
    mask_Bbbroad = util.determine_mask(ln_lam, lam_range_temp, redshift=0, width=3000, custom_lines=[4340.471, 4363.21,  7319.99])  # Hγ & [O III] & [O II]
    mask_Noise   = util.determine_mask(ln_lam, lam_range_temp, redshift=0, width=2000, custom_lines=[7440])                         # Noise near Ca_II trip
    mask_overlap = ~((np.exp(ln_lam) >= (overlap_uv-80)) & (np.exp(ln_lam) <= (overlap_uv+20)))              # Overlap between UV and VIS data see Table 11, XSHOOTER User Manual v. 8 Released On: 2021-06-18
    ################
    if AGNmask: 
        mask_agn = ~(np.exp(ln_lam) <= 4000) 
        mask_combined = mask_Narrow & mask_Broad & mask_Bbroad & mask_Bbbroad & mask_Noise & mask_overlap & mask_Balmer & mask_agn
    else: mask_combined = mask_Narrow & mask_Broad & mask_Bbroad & mask_Bbbroad & mask_Noise & mask_overlap & mask_Balmer
    ################################

    vel0 = 0 # since de-redshifted spectrum else c*np.log(1 + z) aka. (eq. 8 of Cappellari 2017)               
    start = [vel0, 200.]    # (km/s), starting guess for [V,sigma]
    
    pp = ppxf_fit_and_clean(libs.templates, galaxy, noise, velscale, start, mask_combined, np.exp(ln_lam), lam_temp=libs.lam_temp, apol=degree, mpol=mdegree, plot=plot, quiet=quiet)
    if plot: 
        #plt.title(f"pPXF fit with {sps_name} SPS templates") 
        if AGNmask: plt.axvspan(3000, 4000, facecolor = 'lightcoral', alpha=0.5, label='UV data not included in fit') # Overlap between UV and VIS data
        plt.axvspan(overlap_uv-80, overlap_uv+20, facecolor = 'cyan', alpha=0.5) # Overlap between UV and VIS data
        plt.axvline(overlap_uv, color = 'aqua') 
    
        #plt.ylim(-0.4, 1.5)
        plt.minorticks_on()
    if not quiet:
        
        print('\n' + '#'*16 + ' mask delimiters [Å] [start, end, start, end, ...]:\n')
        print(np.insert(pp.lam[np.where(np.diff(mask_combined.astype(int)) !=0)[0]], 0, pp.lam[0]), '\n')
    
    pp.lam_optimal_template = sps.lam_temp
    sps.templates_reshape = sps.templates.reshape(sps.templates.shape[0], -1)
    sps.sps_name = sps_name
    
    if save: 
        packing = {'ppxf_fit': pp, 'sps_lib': sps}

        # Save multiple class instances as dictionary using pickle
        with open(f'output/{date}_ppxffit.pkl', 'wb') as f:
            pickle.dump(packing, f)

    return pp, sps_name, sps, velscale, mask_combined
    
################################################################
print('#'*64 + '\n' + '#'*32 + f'                      27/01/2022\n' +'#'*64)
pp, sps_name, sps, velscale, mask_combined = runppxf('27_01_2022', lam[~nir], spec[~nir], err[~nir], 
                                                     degree=-1, mdegree=-1, sps_name = 'emiles', R=R[1], AGNmask=True, plot=False, save=True)
idx, n, nn = show_sps('27 Jan. X-shooter Mrk590', pp, sps, sps_name)
plt.savefig(f'output/27_01_22_Mrk590_Xshooter_ppxf.pdf', format='pdf')
save_as_ascii('27_01_2022', pp, sps, idx)

print('#'*64 + '\n' + '#'*32 + f'                      08/01/2022\n' +'#'*64)
pp8, sps_name8, sps8, velscale8, mask_combined8 = runppxf('08_01_2022', lam8[~nir], spec8[~nir], err8[~nir],
                                                          degree=-1, mdegree=-1, sps_name = 'emiles', R=R[1], AGNmask=True, plot=False, save=True)
idx8, n8, nn8 = show_sps('08 Jan. X-shooter Mrk590', pp8, sps8, sps_name8)
plt.savefig(f'output/08_01_22_Mrk590_Xshooter_ppxf.pdf', format='pdf')
save_as_ascii('08_01_2022', pp8, sps8, idx8)

