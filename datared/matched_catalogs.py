"""
Prepare SBC bright star input list for CCMAP

this is saved as a record for reproducibility: it will not do anything, and
must be run interactively.
"""
import matplotlib.pyplot as plt
import numpy as np

from astropy import wcs
from astropy.io import fits
from astropy.table import Table
from matplotlib.path import Path
from shapely.geometry import Polygon, Point


def circ(ra, dec, r, d, min=True):
    dist = np.sqrt((ra - r) ** 2 + (dec - d) ** 2)
    if min:
        fn = np.argmin
    else:
        fn = np.argsort
    idx = fn(dist)
    return dist[idx], idx

f140fn = '13710_M31-CENTER-SBC1_drz_phot.dat'
f140fn = '13710_M31-CENTER-SBC2_drz_phot.dat'
f140fn = '13710_M31-CENTER-SBC3_drz_phot.dat'

phot = Table.read(f140fn, format='ascii.commented_header')
radec = np.column_stack([phot['ra'], phot['dec']])
ply = Polygon(radec).convex_hull

# this is very close to exactly the edges of the footprint, FYI
# reg = ply.to_wkt()

reg = ply.envelope.convex_hull.to_wkt()
# verts: [bottomleft, topleft, topright, bottomright, bottomleft]
# (each above is [ra,dec])
verts = parse_poly(reg)

# push the envelope (!) to not be exactly at (4) data points
offset = 0.0001

verts[:, 0][np.array([2, 3])] += offset  # right sides increase
verts[:, 0][np.array([0, 1, 4])] -= offset  # left sides decrease
verts[:, 1][np.array([0, 3, 4])] -= offset  # bottom decrease
verts[:, 1][np.array([1, 2])] += offset  # top increase

# Read UVIS catalog (reference)
f336fn = '13710_M31-CENTER.gst.fits'
hdu = fits.open(f336fn)
w = wcs.WCS(hdu[1].header)
uvis = hdu[1].data
# Cull to SBC region
# radec = np.column_stack([uvis['ra'], uvis['dec']])
# inds, = np.nonzero(Path(verts).contains_points(radec))
# uvis = uvis[inds].copy()

(dists, idxs) = zip(*[circ(uvis['ra'], uvis['dec'], r, d) for r, d in zip(phot['ra'], phot['dec'])])
idxs_4 = np.array(idxs)[np.array(dists) < 0.00003]
pidx, = np.nonzero(np.array(dists) < 0.00003)
ax.plot(phot['vegamag'][pidx] - uvis['F336W_VEGA'][idxs_4], uvis['F336W_VEGA'][idxs_4], 'o')

1. Cull 336 to SBC footprint
1. Find closest matches in ra, dec space
1. Pick tolerance for not matched
1. Check f140lp mag for not matched, make sense? Perhaps go to .phot or .st
1. add f140lp mags to center.gst catalog.
