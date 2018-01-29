import os
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits

hrc = fits.getdata('data/10571_M31NUCLEUS_F330W_F435W.gst.fits')
m1 = hrc['MAG1_ACS']
m2 = hrc['MAG2_STD'] # I know... looks bad. Data has an extra column making MAG2_ACS actually FLAG1
col = m1 - m2
