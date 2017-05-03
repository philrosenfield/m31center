import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np
from astropy import wcs
from astropy.io import fits

# FILTER	PHOTPLAM	PHOTFLAM	STmag	VEGAmag	ABmag
# F140LP	1527.991	2.713e-17	20.316	20.919	23.088
PHOTFLAM = 2.713e-17
VEGAmag0 = 20.919


def apr_red(fitsfile, outfile, makeplots=False):
    """
    aperture photometry on a fitsfild
    """
    from photutils import (Background2D, DAOStarFinder, CircularAperture,
                           aperture_photometry, MeanBackground)
    from photutils.utils import calc_total_error
    from astropy.stats import sigma_clipped_stats
    # Open the fits file
    hdu = fits.open(fitsfile)
    # Take out the data -- this is often not the 0 extention. If this breaks,
    # try this instead:
    # data = hdu[1].data
    data = hdu[0].data
    # the header is almost always the 0th ext.
    w = wcs.WCS(hdu[0].header)
    nanmask = np.isnan(data)
    data[nanmask] = -99
    newdata = data.copy()
    # you might not need to make the image square -- I did it for my mosaic.
    # square-ify: add -99 values to fill a square
    nx, ny = np.shape(data)
    if nx != ny:
        if ny > nx:
            filler = np.zeros((ny - nx, ny)) - 99
            newdata = np.vstack((data, filler))
        else:
            filler = np.zeros((nx, nx - ny)) - 99
            newdata = np.hstack((data, filler))
    else:
        newdata = data.copy()
    # See photutils documentation, box dimension is the nearest factor to 150px
    boxdim = nearest_factor(ny, 150)
    mask = (newdata == -99)
    # sigma_clip = SigmaClip(sigma=3., iters=10)
    bkg = Background2D(newdata, (boxdim, boxdim), filter_size=(9, 9),
                       mask=mask, bkg_estimator=MeanBackground())
    # data - bg
    data_sub = newdata - bkg.background

    mean, median, std = sigma_clipped_stats(data_sub, sigma=8.0, iters=10,
                                            mask=mask)

    daofind = DAOStarFinder(fwhm=4.0, threshold=10. * std + median)
    sources = daofind(data_sub)

    positions = (sources['xcentroid'], sources['ycentroid'])
    # ACS handbook says 4. is inf.
    apertures = CircularAperture(positions, r=4.)

    effective_gain = hdu[0].header['EXPTIME']

    error = calc_total_error(data_sub, bkg.background, effective_gain)
    phot_table = aperture_photometry(data_sub, apertures, error=error,
                                     mask=mask)

    vegamag = -2.5 * np.log10(phot_table['aperture_sum'] * PHOTFLAM) - VEGAmag0
    # vegamagerr = \
    #     -2.5 * np.log10(phot_table['aperture_sum_err'] * PHOTFLAM) - VEGAmag0
    vegamagerr = verrtest(phot_table['aperture_sum'],
                          phot_table['aperture_sum_err'])

    # see astropy.wcs for info:
    ra, dec = w.all_pix2world(phot_table['xcenter'], phot_table['ycenter'], 1)

    phot_table['ra'] = ra
    phot_table['dec'] = dec
    phot_table['vegamag'] = vegamag
    phot_table['vegamagerr'] = vegamagerr

    fmts = {c: '%.6f' for c in ['xcenter', 'ycenter', 'aperture_sum',
                                'aperture_sum_err', 'ra', 'dec', 'vegamag',
                                'vegamagerr']}
    phot_table.write(outfile, format='ascii.commented_header', formats=fmts)
    if makeplots:
        outfig = outfile + '.pdf'
        axs = photplots(newdata, bkg.background, bkg.background_rms, data_sub,
                        fitsfile, outfig)
        apertures.plot(color='red', lw=1.5, alpha=0.5, ax=axs[0, 1])
        plt.savefig(outfig)
    return vegamag, vegamagerr, ra, dec


def verrtest(flux, flerr):
    mag = -2.5 * np.log10(flux * PHOTFLAM) - VEGAmag0
    merr = -2.5 * np.log10((flux - flerr) * PHOTFLAM) - VEGAmag0
    perr = -2.5 * np.log10((flux + flerr) * PHOTFLAM) - VEGAmag0
    return (np.abs(mag - perr) + np.abs(mag - merr)) / 2


def photplots(data, bkg_image, bkg_rms, data_sub, fitsfile, outfig):
    from astropy.visualization import SqrtStretch
    from astropy.visualization.mpl_normalize import ImageNormalize
    norm = ImageNormalize(stretch=SqrtStretch())
    # for plotting
    mask = (data <= 0)
    data[mask] = 0
    data_sub[mask] = 0
    kw = {'norm': norm, 'cmap': 'Greys', 'origin': 'lower'}
    # kw = {'cmap': 'Greys_r', 'origin': 'lower'}
    bkw = {'cmap': 'Greys_r', 'origin': 'lower'}
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
    [ax.grid() for ax in axs.ravel()]
    im = axs[0][0].imshow(data, **kw)
    axs[0][0].set_title(fitsfile)

    im = axs[1][0].imshow(bkg_image, **bkw)
    axs[1][0].set_title('bg')

    im = axs[1][1].imshow(bkg_rms, **bkw)
    axs[1][1].set_title('bg rms')

    # plot background-subtracted image
    ax = axs[0][1]
    im = ax.imshow(data_sub, **kw)
    ax.set_title('bg-sub data')
    plt.tight_layout()
    plt.savefig(outfig)
    return axs



def factors(n):
    """ find the factors of n """
    import itertools
    flatten_iter = itertools.chain.from_iterable
    return np.sort(list(set(flatten_iter((i, n // i)
                                         for i in range(1, int(n**0.5) + 1)
                                         if n % i == 0))))

def nearest_factor(num, nfact):
    """ find the nearest factor to a number"""
    fs = factors(num)
    return fs[np.argmin(np.abs(nfact - fs))]

def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="source extract!")

    parser.add_argument('fitsfiles', nargs='*', type=str,
                        help='fits files')
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    fig, ax = plt.subplots()
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    for fitsfile in args.fitsfiles:
        outfile = fitsfile.replace('.fits', '_phot.dat')
        vegamag, vegamagerr, ra, dec = apr_red(fitsfile, outfile,
                                               makeplots=True)
        inds = np.isfinite(vegamag)
        h, b = np.histogram(vegamag[inds], bins=np.arange(20, 26, 0.1))
        ax.plot(b[1:], h, ls='steps-', label=fitsfile)
        ax.set_yscale('log')
        ax.set_ylim(3, ax.get_ylim()[1])
        ax1.plot(ra, dec, '.', label=fitsfile)
        ax2.plot(vegamag[inds], vegamagerr[inds], '.', label=fitsfile)
    ax.legend(loc='best')
    ax1.legend(loc='best')
    ax2.legend(loc='best')
    fig.savefig('F140LP_LF.pdf')
    fig1.savefig('F140LP_radec.pdf')
    fig2.savefig('F140LP_err.pdf')
    return


if __name__ == "__main__":
    sys.exit(main())
