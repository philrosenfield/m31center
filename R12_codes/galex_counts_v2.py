'''quick script to sum counts in galex images within specified regions'''
import sys
from matplotlib.path import Path
from matplotlib.colors import LogNorm
import matplotlib.pylab as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


'''
import aplpy
from astropy.io import fits
%pylab

gc = aplpy.FITSFigure('brick01_center.fits')
gc.show_colorscale(stretch='log', vmin=1e-2,vmax=0.5, cmap='Greys')


fname = 'b01-6filt-cut.fits'
hdu = fits.open(fname)
mag = hdu[1].data.f275w_vega
recovered, = np.nonzero((mag < 24))
ra = hdu[1].data.ra[recovered]
dec = hdu[1].data.dec[recovered]
mag = mag[recovered]
all_inds, = np.nonzero((ra < 10.73) & (ra > 10.5) & \
                       (dec > 41.23) & (dec < 41.33))
sort = np.argsort(mag[all_inds])[::-1]
s = np.zeros(len(sort)) + 20
bright, = np.nonzero(mag[all_inds][sort] < 22)
s[bright] = 50
gc.show_markers(ra[all_inds][sort][bright], dec[all_inds][sort][bright], s=s)

gc.list_layers()
bright, = np.nonzero(mag[all_inds][sort] < 20)
gc.show_markers(ra[all_inds][sort][bright], dec[all_inds][sort][bright])
gc.set_theme('publication')
cat
daophotcat = 'hst_12058_10_wfc3_uvis_f275w_daophot_trm.cat'
cols = 'X-Center   Y-Center         RA          DEC        ID   MagAp1   MagErr   MagAp2   MagErr     MSky    Stdev       Flux   TotMag   MagErr       CI    Flags'
cols = cols.split()
pcat = np.genfromtxt(daophotcat, names=cols)
pbright2, = np.nonzero(pcat['MagAp2'] < 23.3)
gc.show_markers(pcat['RA'][pbright2], pcat['DEC'][pbright2])
pbright2, = np.nonzero(pcat['MagAp2'] < 22.)
gc.show_markers(pcat['RA'][pbright2], pcat['DEC'][pbright2])

in.fits out.fits ra dec xsize
'X-Center   Y-Center         RA          DEC        ID   MagAp1   MagErr   MagAp2   MagErr     MSky    Stdev       Flux   TotMag   MagErr       CI    Flags'

4348.047, 3503.227, 10.6842213, 41.2689133, 4551.0, 22.668, 0.031, 21.104, 0.021, 0.218, 0.041, 16.909, -999.0, -999.0, 1.564, 1.0)

infits = '/Users/phil/research/PHAT/Data/FullBricks/bens_drz/B01.fits'
outfits = 'brick01_center.fits'
ra = 10.684875
dec = 41.268838
xsize = 10.77-10.71
ysize = 41.335 - 41.3008
montage.mSubimage(infits, outfits, ra, dec, xsize, ysize)
'''

def mag_plot():
    fnames, p1, p2, p3, _, _ = load_data(uvis=True, radec=True)
    hdu = fits.open(fnames[0])
    mag = hdu[1].data.f275w_vega
    recovered, = np.nonzero((mag < 24))

    c1 = np.array([10.69346, 41.26869])
    c2 = np.array([10.70488, 41.26733])
    c3 = np.array([10.71671, 41.26608])

    ra = hdu[1].data.ra[recovered]
    dec = hdu[1].data.dec[recovered]
    mag = mag[recovered]
    grid = np.column_stack((ra, dec))

    cols = ['c', 'm', 'y']
    ax = plt.subplots()[1]

    #img = fits.open('hst_12058_10_wfc3_uvis_f275w_drz.fits')
    #ax.imshow(img[2].data, cmap=plt.cm.Greys)
    #all_inds = np.array([], dtype=np.int)
    for p, c in zip([p1, p2, p3], [c1, c2, c3]):
        ax.plot(p[:, 0], p[:, 1], color='k')
        circ = plt.Circle(c, radius=0.0075, color='k', fill=False)
        ax.add_artist(circ)
        ax.plot(c[0], c[1], '+', color='k', ms=30)
        #path = Path(p)
        #isInside = path.contains_points(grid)
        #inds, = np.nonzero(isInside)
        #all_inds = np.append(all_inds, inds)

    all_inds, = np.nonzero((ra < 10.73) & (ra > 10.5) & \
                           (dec > 41.23) & (dec < 41.33))
    sort = np.argsort(mag[all_inds])[::-1]
    s = np.zeros(len(sort)) + 20
    bright, = np.nonzero(mag[all_inds][sort] < 20)
    s[bright] = 50
    l = ax.scatter(ra[all_inds][sort], dec[all_inds][sort],
                   c=mag[all_inds][sort], cmap=plt.cm.Greys_r, alpha=0.9, s=s,
                   linewidths=0)
    [ax.annotate(r'$%.2f$' % mag[all_inds][sort][b], (ra[all_inds][sort][b],
                dec[all_inds][sort][b])) for b in bright]
    plt.colorbar(l)
    ax.set_xlim(ax.get_xlim()[::-1])
    ticklabel_format(style='scientific', useOffset=False)
    ax.set_xlabel('$RA$')
    ax.set_ylabel('$DEC$')
 
def load_data(fuv=True):
    vmin = 0.001
    if fuv:
        fname = 'PS_M31_MOS00-fd-int.fits'
        vmax = 0.1
    else:
        fname = 'PS_M31_MOS00-nd-int.fits'
        vmax = 7.

    hdu = fits.open(fname)
    w = WCS(hdu[0].header)

    r1 = np.array([[10.69525, 41.27506], [10.68487, 41.27153],
                   [10.68075, 41.28067], [10.69117, 41.28408],
                   [10.69525, 41.27506]])
    r2 = np.array([[10.70446, 41.28147], [10.69404, 41.27792],
                   [10.68992, 41.28706], [10.70029, 41.29053],
                   [10.70446, 41.28147]])
    r3 = np.array([[10.71317, 41.28903], [10.70279, 41.28544],
                   [10.69871, 41.29458], [10.70904, 41.29808],
                   [10.71317, 41.28903]])

    p1 = w.wcs_world2pix(r1, 1)
    p2 = w.wcs_world2pix(r2, 1)
    p3 = w.wcs_world2pix(r3, 1)

    return fname, hdu, p1, p2, p3, vmin, vmax


def main(plot=True):
    '''make plot, print counts in each region'''
    if plot:
        cols = ['c', 'm', 'y']

    for fuv in [True, False]:
        fname, hdu, p1, p2, p3, vmin, vmax = load_data(fuv=fuv)
        # hacky way to get a 2D grid of indices from data
        grid = np.meshgrid(np.arange(hdu[0].shape[0]),
                           np.arange(hdu[0].shape[1]))
        grid = np.concatenate(grid).ravel()
        grid = grid.reshape(2, len(grid)/2).T
        if plot:
            ax = plt.subplots()[1]
            ax.imshow(hdu[0].data, cmap=plt.cm.Greys_r,
                      norm=LogNorm(vmin=vmin, vmax=vmax))
        for k, p in enumerate([p1, p2, p3]):
            path = Path(p)
            isInside = path.contains_points(grid)

            inds, = np.nonzero(isInside)
            counts = np.sum([hdu[0].data[i, j] for i, j in grid[inds]])
            print counts
            if plot:
                ax.plot(p[:, 0], p[:,1], color=cols[k],
                        label=r'$%.6g$' % counts)
                [ax.plot(grid[i][0], grid[i][1], '.', color=cols[k], alpha=0.3)
                 for i in inds]
        if plot:
            ax.set_xlim(1820, 1940)
            ax.set_ylim(1910, 2000)
            ax.set_xlabel(r'$X$', fontsize=20)
            ax.set_ylabel(r'$Y$', fontsize=20)
            ax.legend(loc=0, frameon=False, fontsize=16)
            ax.set_title('%s' % fname, fontsize=16)
            imgname = '%s.png' % fname
            plt.savefig(imgname)
            print 'wrote %s' % imgname
    """
    Is it odd that nuv/fuv counts are almost exactly 10? 
    PS_M31_MOS00-fd-int.fits: 
    # visit counts AB MAG ETC link
    1 4.32889 17.23 http://etc.stsci.edu/etc/results/ACS.im.642422/
    2 3.10702 17.59 http://etc.stsci.edu/etc/results/ACS.im.642423/
    3 2.31266 17.91 http://etc.stsci.edu/etc/results/ACS.im.642424/
    PS_M31_MOS00-nd-int.fits: (all had bright warnings)
    # counts AB MAG ETC link
    1 43.1971 14.73 http://etc.stsci.edu/etc/results/ACS.im.642426/
    2 31.1601 15.09 http://etc.stsci.edu/etc/results/ACS.im.642427/
    3 23.3791 15.40 http://etc.stsci.edu/etc/results/ACS.im.642428/
    
    """
if __name__ == '__main__':
    plot = sys.argv[1]
    main(plot=plot)
