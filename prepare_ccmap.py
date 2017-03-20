"""
Prepare SBC bright star input list for CCMAP

this is saved as a record for reproducibility: it will not do anything, and
must be run interactively.
"""
import aplpy
import matplotlib.pyplot as plt
import numpy as np

from astropy import wcs
from astropy.io import fits
from astropy.table import Table

# reduced photometry table from SBC
f140fn = '13710_M31-CENTER-SBC1_drz_phot.dat'
f140fn = '13710_M31-CENTER-SBC2_drz_phot.dat'
f140fn = '13710_M31-CENTER-SBC3_drz_phot.dat'

# Drizzled F336W image to align SBC to
f336fn = '13710_M31-CENTER_F336W_drz.fits'

# SBC bright limit (plot stars with this mag and brighter)
bright = 21.5
# bright = 20  # SBC3 works better

# Load SBC image and cull bright stars
phot = Table.read(f140fn, format='ascii.commented_header')
isbc, = np.nonzero(phot['vegamag'] < bright)

# Load reference image and grab WCS.
hdu = fits.open(f336fn)
w = wcs.WCS(hdu[1].header)

# Plot reference
gc = aplpy.FITSFigure(f336fn)
# May need to mess with this if changing reference image
gc.show_grayscale(stretch='log', vmid=0.005, vmin=0.01)

# Plot SBC bright stars
gc.show_markers(phot['ra'][isbc], phot['dec'][isbc], layer='140_3',
                marker='o', edgecolor='blue', facecolor='none')


def circ(ra, dec):
    """RA, Dec to X Y"""
    r = phot['ra']
    d = phot['dec']
    idxr = np.argmin(np.sqrt((ra - r) ** 2 + (dec - d) ** 2))
    if np.abs(d[idxr] - dec) > 1e-4:
        print('warning, bad clicking! {} {}: {}'.format(ra, dec,
                                                        np.abs(d[idxr] - dec)))
    return '{0:.6f} {1:.6f}'.format(phot['xcenter'][idxr], phot['ycenter'][idxr])

# Interactively, zoom around and select points keeping in mind which are from
# SBC (red circles) and which are from F336W (image).
# Change n=4 to however many clicks you want to record at a time:
# w.all_pix2world(plt.ginput(n=4), 1)

# Reformat to a file for ccmap (see comment below)
# Last two columns above need to be X,Y in SBC data
# [circ(r,d) for r, d in zip(radecs[:, 0], radecs[:, 1])]


"""
Re-formatted from w.all_pix2world(plt.ginput(n=4), 1) output:
SBC1
10.6922117 ,  41.28278105, 176.150121, 64.330914
10.68996233,  41.2819457 , 205.113797, 335.005810
10.68886333,  41.27958616, 486.484494, 557.811860
10.69000859,  41.27976229, 505.790302, 429.533176
10.69415449,  41.27515858, 1275.027092, 224.291839
10.69240519,  41.27551961, 1164.436891, 386.645743
10.69086282,  41.27587852, 1060.951572, 525.564472
10.68829924,  41.2754164 , 1174.268281, 765.062701
10.68651453,  41.27777828, 648.987562, 882.272635
10.68556523,  41.27758137, 643.042637, 987.127272
10.68500777,  41.27734295, 655.783965, 1057.674264
10.68348769,  41.27914699, 357.548860, 1126.966665

SBC
radecs=np.array([[10.69184468,  41.28281427],
[10.68960635,  41.28199757],
[10.6884577 ,  41.27965463],
[10.68965262,  41.27981417],
[10.69372406,  41.27523536],
[10.69199959,  41.27558809],
[10.69049309,  41.27595735],
[10.68874605,  41.27466367],
[10.68613099,  41.27788406],
[10.68522032,  41.27765189],
[10.68461319,  41.27743006],
[10.68310689,  41.27920714]])


SBC2
radecs = np. array([[ 10.69083534,  41.28663228],
       [ 10.69713143,  41.28426676],
       [ 10.69866368,  41.28359356],
       [ 10.7013921 ,  41.28199667],
       [ 10.6945749 ,  41.28026509],
       [ 10.69541083,  41.28040889],
       [ 10.69842582,  41.28006631],
       [ 10.7001921 ,  41.28006036]])

WFC
10.6912981 ,  41.28659858, 172.091954 1226.432073
10.69757178,  41.28424052, 716.052922 695.352029
10.69910402,  41.2835673 , 858.773051 572.193290
10.70186477,  41.28197975, 1171.583015 365.738666
10.6950376 ,  41.28023137, 1170.700950 1142.755996
10.69587353,  41.28037517, 1179.379238 1052.500566
10.69886613,  41.28004006, 1330.788017 762.161570

SBC3
radecs = np.array([[ 10.70113369,  41.29170334],
       [ 10.70200431,  41.29293723],
       [ 10.70361258,  41.29225521],
       [ 10.70380407,  41.29216492],
       [ 10.7105486 ,  41.29235585],
      [ 10.70866842,  41.2914857 ],
      [ 10.70837086,  41.29059925],
      [ 10.70689013,  41.28987163],
      [ 10.70521082,  41.28804045],
       [ 10.70462685,  41.28870884],
       [ 10.70313559,  41.28947011],
       [ 10.70162273,  41.28927899]])

10.70154761,  41.29172276, 563.700032, 1182.947220
10.70242474,  41.29296762, 424.282107, 1035.884645
10.70404112,  41.29226974, 574.843651, 903.784289
10.70421639,  41.29221117, 591.511248, 890.115320
10.71100635,  41.2923606, 807.867078, 191.191273
10.70911155,  41.29149533, 859.946138, 423.586236
10.70882049,  41.29061985, 969.648366, 494.790084
10.70731054,  41.28990201, 1012.770246, 681.505160
10.70565232,  41.28807692, 1203.490048, 937.491578
10.70506186,  41.28873435, 1092.850357, 967.320018
10.70354139,  41.28950539, 935.108225, 1085.060691
10.70202203,  41.2893033, 909.803903, 1248.572026
"""
