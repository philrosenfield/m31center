'''quick script to sum counts in galex images within specified regions'''
import aplpy
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
 
def load_data(fuv=True)
    vmin = 0.001
    if fuv:
        fname = 'PS_M31_MOS00-fd-int.fits'
        vmax = 0.1
    else:
        fname = 'PS_M31_MOS00-nd-int.fits'
        vmaxs = 7.

    hdu = fits.open(fname)
    
    
    r1 = np.array([10.69525, 41.27506, 10.68487, 41.27153, 10.68075, 41.28067,
                   10.69117, 41.28408, 10.69525, 41.27506])
    r2 = np.array([10.70446, 41.28147, 10.69404, 41.27792, 10.68992, 41.28706,
                   10.70029, 41.29053, 10.70446, 41.28147])
    r3 = np.array([10.71317, 41.28903, 10.70279, 41.28544, 10.69871, 41.29458,
                   10.70904, 41.29808, 10.71317, 41.28903])     
    
    
    
    p1 = r1.reshape(5, 2)
    p2 = r2.reshape(5, 2)
    p3 = r3.reshape(5, 2)


def load_data_old(uvis=True, radec=False):
    '''
    load data, regions hard coded from run with ds9 converted to X,Y
    # Region file format: DS9 version 4.1
    fk5
    polygon(10.722458,41.262667,10.719333,41.271917,10.70875,41.268944,10.711917,41.25975) 
    polygon(10.710708,41.263972,10.707458,41.273333,10.696833,41.270306,10.7,41.261111)
    polygon(10.699208,41.265472,10.695958,41.274528,10.685458,41.271722,10.688625,41.262472)
    '''
    if not uvis:
        fnames = ['PS_M31_MOS00-fd-int.fits',
                  'PS_M31_MOS00-nd-int.fits']
        vmins = [0.001, 0.001]
        vmaxs = [0.1, 7.]

        r1 = np.array([1834.882, 1910.8434,
                       1840.5311, 1933.0404,
                       1859.6186, 1925.8985,
                       1853.8969, 1903.8341,
                       1834.882, 1910.8434])

        r2 = np.array([1856.0811, 1913.9666,
                       1861.9526, 1936.4309,
                       1881.1158, 1929.1583,
                       1875.3973, 1907.0932,
                       1856.0811, 1913.9666])

        r3 = np.array([1876.8284, 1917.5594,
                       1882.6967, 1939.2912,
                       1901.6351, 1932.5545,
                       1895.9195, 1910.3553,
                       1876.8284, 1917.5594])
    elif not radec:
        fnames = ['hst_12058_10_wfc3_uvis_f275w_drz.fits']
        r1 = np.array([1736.8267, 2935.532,
                       1950.3779, 3775.9009,
                       2673.0317, 3505.7225,
                       2456.7051, 2670.4157,
                       1736.8267, 2935.532])

        r2 = np.array([2539.2984, 3053.999,
                       2761.2673, 3904.4815,
                       3486.811, 3629.4961,
                       3270.6034, 2794.0687,
                       2539.2984, 3053.999])

        r3 = np.array([3324.6636, 3190.2898,
                       3546.5179, 4013.091,
                       4263.5444, 3758.2786,
                       4047.4518, 2917.8259,
                       3324.6636, 3190.2898])
    else:
        fnames = ['b01-6filt-cut.fits']
        r1 = np.array([10.722458, 41.262667,
                       10.719333, 41.271917,
                       10.70875, 41.268944,
                       10.711917, 41.25975,
                       10.722458, 41.262667])
        r2 = np.array([10.710708, 41.263972,
                       10.707458, 41.273333,
                       10.696833, 41.270306,
                       10.7, 41.261111,
                       10.710708, 41.263972])
        r3 = np.array([10.699208, 41.265472,
                       10.695958, 41.274528,
                       10.685458, 41.271722,
                       10.688625, 41.262472,
                       10.699208, 41.265472])
        c1 = np.array([10.69346, 41.26869]) 
        c2 = np.array([10.70488, 41.26733])
        c3 = np.array([10.71671, 41.26608])

    p1 = r1.reshape(5, 2)
    p2 = r2.reshape(5, 2)
    p3 = r3.reshape(5, 2)
    vmins = []
    vmaxs = []
    return fnames, p1, p2, p3, vmins, vmaxs

def main(uvis=True, plot=True, radec=False):
    '''make plot, print counts in each region'''
    if plot:
        cols = ['c', 'm', 'y']

    fnames, p1, p2, p3, vmins, vmaxs = load_data(uvis=uvis, ra=radec)
    for l, fname in enumerate(fnames):
        hdu = fits.open(fname)
        if not radec:
            # hacky way to get a 2D grid of indices from data
            grid = np.meshgrid(np.arange(hdu[1].shape[0]),
                               np.arange(hdu[1].shape[1]))
            grid = np.concatenate(grid).ravel()
            grid = grid.reshape(2, len(grid)/2).T
        else:
            grid = np.column_stack((hdu[1].data['ra'], hdu[1].data['dec']))
        if plot:
            ax = plt.subplots()[1]
            if 'uvis' in fname:
                ax.imshow(hdu[2].data, cmap=plt.cm.Greys)
            elif not radec:
                ax.imshow(hdu[0].data, cmap=plt.cm.Greys_r,
                          norm=LogNorm(vmin=vmins[l], vmax=vmaxs[l]))
            else:
                ax.plot(hdu[1].data['ra'], hdu[1].data['dec'], '.', color='k',
                        alpha=0.3) 
        for k, p in enumerate([p1, p2, p3]):
            path = Path(p)
            isInside = path.contains_points(grid)

            inds, = np.nonzero(isInside)
            if not radec:
                counts = np.sum([hdu[1].data[i, j] for i, j in grid[inds]])
            else:
                magin = hdu[1].data['f275w_vega'][inds]
                mag = magin[magin < 99]
                counts = np.sum(10 ** (-0.4 * mag))
            if not uvis:
                print counts
            else:
                print counts * 3.301e-18
            if plot:
                ax.plot(p[:, 0], p[:,1], color=cols[k],
                        label=r'$%.6g$' % counts)
                if not uvis or not radec:
                    [ax.plot(grid[i][0], grid[i][1], '.', color=cols[k], alpha=0.3)
                     for i in inds]
        if plot:
            if not 'uvis' in fname and not ra:
                ax.set_xlim(1820, 1940)
                ax.set_ylim(1960, 1890)
            if radec:
                ax.set_xlabel(r'$\alpha$', fontsize=20)
                ax.set_ylabel(r'$\delta$', fontsize=20)
            else:
                ax.set_xlabel(r'$X$', fontsize=20)
                ax.set_ylabel(r'$Y$', fontsize=20)
            ax.legend(loc=0, frameon=False, fontsize=16)
            imgname = '%s.png' % fname
            plt.savefig(imgname)
            print 'wrote %s' % imgname

if __name__ == '__main__':
    uvis = sys.argv[1]
    plot = sys.argv[2]
    radec = sys.argv[3]
    main(uvis=uvis, plot=plot, radec=radec)
