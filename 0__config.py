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
import matplotlib.pyplot as plt

plt.rcParams.update(font_sizes)
plt.rc('image', cmap='cubehelix')
plt.rcParams['lines.linewidth'] = 2

from matplotlib.colors import LinearSegmentedColormap
Figtree = LinearSegmentedColormap.from_list("Figtree",  ['#904c77', '#e49ab0', '#331F1F', '#c3d350', '#809848', '#2E4D2D']*4)
Minty = LinearSegmentedColormap.from_list("Minty", ["#22577a","#38a3a5","#57cc99","#80ed99","#c7f9cc"]*5)

print('### Setup initialized ###')

# Only install packages if plotbin is missing
try:
    import plotbin
except ImportError:
    import subprocess, sys
    for package in ["ppxf", "spectral_cube", "vorbin", "plotbin", "mpdaf", "tqdm", "astroML", "photutils", "pafit", 'reproject']:
        subprocess.check_call([sys.executable, "-m", "pip", "install", package], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    print('### Setup complete ###')

'''               
"pyarrow", "typhon", "marvin", "photutils", "pafit", "lime-stable", "http.client", "pillow", "imageio"
import astropy
import numpy as np
import matplotlib
import ppxf
import vorbin
import plotbin

# Check versions
print("Astropy version:", astropy.__version__)
print("Numpy version:", np.__version__)
print("Matplotlib version:", matplotlib.__version__)
print("pPXF version:", ppxf.__version__)
print("Voronoi 2D Binning (vorbin) version:", vorbin.__version__)
print("Plotbin version:", plotbin.__version__)
'''