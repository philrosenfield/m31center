"""Mask away the saturated part of the irac image"""
import numpy as np
from footprints import parse_footprint, write_polygons
from matplotlib.path import Path

incontour = 'irac/contours.reg'
outcontour = 'irac/contours_masked.reg'
irac_mask = 'irac/irac_bulge_mask.dat'

all_regs = parse_footprint(incontour, return_coords=True)
# A little trick from ds9 mess, the longest contours are the ones I want.
inds = np.argsort([len(r) for r in all_regs])[::-1]
regs = np.array(all_regs)[inds][:9]

# I used ginput to select the coords to make this mask.
mask = np.genfromtxt(irac_mask)

# exclude the mask regions.
mregs = [r[~Path(mask).contains_points(r)] for r in regs]

write_polygons(mregs, outcontour)
print('wrote {0:s}'.format(outcontour))
