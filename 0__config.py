import matplotlib.pyplot as plt
# Adjusting default parameters
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.size'] = 8
plt.rcParams['ytick.major.width'] = 1.5
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['xtick.minor.width'] = 1.0
plt.rcParams['ytick.minor.size'] = 4
plt.rcParams['ytick.minor.width'] = 1.0
plt.rcParams['figure.figsize'] = [12, 10] 

# Define default font sizes
font_sizes = {
    'axes.titlesize': 20,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'legend.title_fontsize': 18,
}

plt.rcParams.update(font_sizes)
plt.rc('image', cmap='cubehelix')
plt.rcParams['lines.linewidth'] = 2

print('### Setup initialized ###')

# Only install packages if plotbin is missing
try:
    import plotbin
except ImportError:
    import subprocess, sys
    for package in [ "ppxf", "plotbin", "spectral_cube", "astroML"]:
        subprocess.check_call([sys.executable, "-m", "pip", "install", package], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    print('### Setup complete ###')

'''   

# Only install packages if plotbin is missing
try:
    import plotbin
except ImportError:
    import subprocess, sys
    for package in ["ppxf", "spectral_cube", "vorbin", "plotbin",
import astropy
import numpy as np
import matplotlib
import ppxf

# Check versions
print("Astropy version:", astropy.__version__)
print("Numpy version:", np.__version__)
print("Matplotlib version:", matplotlib.__version__)
print("pPXF version:", ppxf.__version__)
print("Voronoi 2D Binning (vorbin) version:", vorbin.__version__)
print("Plotbin version:", plotbin.__version__)
'''