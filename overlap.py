import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from matplotlib.path import Path
from fitshelper.footprints import parse_poly, parse_footprint

import seaborn
seaborn.set()

uvisfile = '13710_M31-CENTER.gst.fits'
regfile = '13710_M31-CENTER-SBC_drz_footprint_ds9_linear.reg'
reg = parse_poly(parse_footprint(regfile)[0].split(';')[1])

uvis = fits.open(uvisfile)
ura, udec = uvis[1].data['RA'], uvis[1].data['DEC']
upoints = np.column_stack([ura, udec])

uinds, = np.nonzero(Path(reg).contains_points(upoints))
