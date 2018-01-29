"""Make thumbnail of a cutout of P. Barmby's IRAC 3.6 um image"""

import os
from astropy.io import fits
from astropy.wcs import WCS
filename = '../data/m31_irac1.fits'

fn, ext = os.path.splitext(filename)
outname = '{0:s}_thumb{1:s}'.format(fn, ext)

hdu = fits.open(filename)
header = hdu[0].header
wcs = WCS(header)

# center values (int pixels)
x = 2650
y = 8260
# half the thumbnail width (int pixels)
hw = 250

# fig, ax = plt.subplots()
# ax.imshow(hdu[0].data[x-hw:x+hw,y-hw:y+hw])

r, d = wcs.wcs_pix2world(y, x, 1)
keys = {'CRPIX1': hw,
        'CRPIX2': hw,
        'CRVAL1': float(r),
        'CRVAL2': float(d)}

header.update(keys)
thumb = hdu[0].data[x-hw:x+hw, y-hw:y+hw]

fits.writeto(outname, thumb, header=header, overwrite=True)
