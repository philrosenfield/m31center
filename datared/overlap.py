import glob

import matplotlib.pyplot as plt
import numpy as np
import seaborn
from astropy.io import fits
from matplotlib.path import Path
from shapely.geometry import Polygon

from fitshelper.footprints import find_footprint, parse_footprint, parse_poly

seaborn.set()

"""
# uvisfile = '13710_M31-CENTER.gst.fits'
# regfile = '13710_M31-CENTER-SBC_drz_footprint_ds9_linear.reg'
# reg = parse_poly(parse_footprint(regfile)[0].split(';')[1])

# uvis = fits.open(uvisfile)
# ura, udec = uvis[1].data['RA'], uvis[1].data['DEC']
# upoints = np.column_stack([ura, udec])

# uinds, = np.nonzero(Path(reg).contains_points(upoints))
"""

all_flc = glob.glob('*flc.fits')
sbcs = glob.glob('jcj*')

pct_tol = 49.
# all_flc = !! ls *flc.fits
# sbcs = !! ls jcj*

wfc3feet = [find_footprint(f) for f in all_flc]
wfc3plys = [Polygon(f) for f in wfc3feet]
sbcfeet = [find_footprint(f) for f in sbcs]
sbcplys = [Polygon(f) for f in sbcfeet]

flcs = []
fig, ax = plt.subplots()
for i, sbcply in enumerate(sbcplys):
    for j, wfc3ply in enumerate(wfc3plys):
        if sbcply.intersects(wfc3ply):
            iply = sbcply.intersection(wfc3ply)
            olap = iply.area / sbcply.area * 100
            if olap > pct_tol:
                # print('{} intersects with {}'.format(sbcs[i], all_flc[j]))
                flcs.append(all_flc[j])
                ax.plot(wfc3feet[j][:, 0], wfc3feet[j][:, 1])
                ax.plot(sbcfeet[i][:, 0], sbcfeet[i][:, 1])

plt.savefig('overlaps.png')
with open('flc.list', 'w') as out:
    [out.write('{:s}\n'.format(flc)) for flc in np.unique(flcs)]
