"""
Prepare SBC bright star input list for CCMAP

this is saved as a record for reproducibility: it will not do anything, and
must be run interactively.
"""
import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits
from astropy.table import Table
from matplotlib.path import Path


# this could easily be sped up...
def circ(ra, dec, r, d, min=True):
    dist = np.sqrt((ra - r) ** 2 + (dec - d) ** 2)
    if min:
        fn = np.argmin
    else:
        fn = np.argsort
    idx = fn(dist)
    return dist[idx], idx


# Read UVIS catalog (reference). Ok ... this could all be done with
# astropy.table.Table.read however, it wasn't immediately obvious how to
# use astropy.table.Table.write(format='fits') and maintain the fits header.
f336fn = '13710_M31-CENTER.gst.fits'
header, bintab = fits.open(f336fn)

# Aligned and reduced SBC catalogs
f140fns = ['13710_M31-CENTER-SBC1_drz_phot.dat',
           '13710_M31-CENTER-SBC2_drz_phot.dat',
           '13710_M31-CENTER-SBC3_drz_phot.dat']

# Initialize new columns (Ben and Andy like 99 value for NaN)
newcol = np.zeros(len(bintab.data)) + 99.
newcolerr = np.zeros(len(bintab.data)) + 99.

# matching tolerance
tol = 0.00005

for f140fn in f140fns:
    # Read in ra and dec from SBC data
    phot = Table.read(f140fn, format='ascii.commented_header')
    radec = np.column_stack([phot['ra'], phot['dec']])
    # find closest matches to UVIS ra dec
    (dists, idxs) = zip(*[circ(bintab.data['RA'], bintab.data['DEC'], r, d)
                          for r, d in zip(phot['ra'], phot['dec'])])

    # filter by tolerance
    matches = np.array(dists) < tol
    idxs_4 = np.array(idxs)[matches]
    pidx, = np.nonzero(matches)
    print('with tol: {}, matched {} of {}.'.format(tol, len(pidx), len(idxs)))

    # Diag plotting
    # ax.plot(phot['vegamag'][pidx] - bintab.data['F336W_VEGA'][idxs_4],
    #         bintab.data['F336W_VEGA'][idxs_4], 'o')

    # add data to new column.
    # With this loop, farther images from M31 center will overwrite
    # sources from inner, if there are any. (based on the list order)
    newcol[idxs_4] = phot['vegamag'][pidx]
    newcolerr[idxs_4] = phot['vegamagerr'][pidx]

# Make new columns
col = fits.Column(name='F140LP_VEGA', format='E', array=newcol)
colerr = fits.Column(name='F140LP_ERR', format='E', array=newcolerr)
newcols = fits.ColDefs([col, colerr])

# make new HDU
newtab = fits.BinTableHDU.from_columns(bintab.columns + newcols)
newhdu = fits.HDUList([header, newtab])

# filename
filters = '_'.join([f.replace('_VEGA','') for f in bintab.columns.names if
                    f.startswith('F') and f.endswith('VEGA')])
newname = '13710_M31-CENTER_{}.gst.fits'.format(filters)

# save
newhdu.writeto(newname, overwrite=True)
print('wrote {}'.format(newname))
