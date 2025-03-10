################################################################

pPXF for Xshooter data of Mrk 590

################################################################

Author: Freja Amalie Nørby

Email: freja.noerby@nbi.ku.dk

Date: 08. March 2025

################################################################

The 0__config.py file ensures that all necessary packages are installed, automatically handling the installation of any missing dependencies. 
Although I provide a modified version of ppxf, it doesn't include all the files. Therefore, you will still need to have the complete ppxf package installed. 

**To run the fitting procedure** execute the ppxf_Mrk_590_xshooter.py script and either uncomment my preset function call or use the following:

pp, sps_name, sps, velscale, mask_combined = runppxf('identifier', lam[range], spec[range], err[range], degree=-1, mdegree=-1, sps_name='emiles',R=R[1], AGNmask=True, plot=False, save=True)

Key Parameters: 

    degree and mdegree are the same as in ppxf, representing the additive Legendre polynomial and multiplicative polynomial to add to the model. 
          A value of -1 means no contribution, and they follow standard polynomial conventions.
          
    sps_name: This refers to the name of the stellar population synthesis (SPS) library. Possible values are: 'fsps', 'galaxev', 'emiles', 'xsl'
          (Note: only MILES and XSL have high resolution beyond 7500Å. For other libraries, restrict the fit to wavelengths beyond 7500Å.)
          
    save: If set to True, a pickle file will be saved. For saving as an ASCII file, use the save_as_ascii('27_01_2022', pp, sps, idx) function.
    
    AGNmask=True: This option will mask the spectrum below 4000Å.

**Note:** The output is not saved to an output log; instead, it is printed to the terminal. Make sure to save the output manually if needed.   

################################################################

Mrk 590 from XSHOOTER observations under Program ID 108.229T.001

################################################################
               Airmass    seeing  
08/01/2022      1.208	   0.51
27/01/2022      1.448	   0.78

UVB slit width = 1.6 , 	R=(λ/Δλ) = 3200, sampling = 9.9
VIS slit width = 1.5 , 	R=(λ/Δλ) = 5000, sampling = 9.5  
NIR slit width = 0.9 , 	R=(λ/Δλ) = 5600, sampling = 3.7 (3.6)

From archive.eso.org/hdr?DpId=XSHOO.2022-01-27T01:00:32.145

################################

chi2/DOF:    Reduced chi-squared 

################################

As a general rule, 
        χ^2 ≫ 1 indicates a poor model fit. 
        
        χ^2 > 1 indicates that the fit has not fully captured the data (or that the error variance has been underestimated)
        
        χ^2 around 1  indicates that the extent of the match between observations and estimates is in accord with the error variance. 
        
        χ^2 < 1 indicates that the model is "overfitting" the data: either the model is improperly fitting noise, or the error variance has been overestimated

See: https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic

############################################

Stellar Population Synthesis (SPS) templates 

############################################

I have Used the Emiles SPS library for my fitting: 
The [E-MILES](http://miles.iac.es/) SPS model templates by [Vazdekis et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016MNRAS.463.3409V) with Padova isochrones and Salpeter IMF.

http://research.iac.es/proyecto/miles/pages/webtools/tune-ssp-models.php

As mentioned, four libraries are distributed with PPXF (or available here: [[micappe/ppxf_data](https://github.com/micappe/ppxf_data)]). However, it's also possible to use other libraries, provided they are packaged correctly as NPZ files
