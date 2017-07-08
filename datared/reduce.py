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
APR_COR = 0.804  # 0.5" Aperture correction

def reduce(fitsfile, outfile, makeplots=False):
    from photutils import (Background2D, DAOStarFinder, CircularAperture,
                           aperture_photometry, MedianBackground,
                           make_source_mask, SigmaClip)
    from photutils.utils import calc_total_error
    from astropy.stats import sigma_clipped_stats

    # background variables
    bkg_boxsize = 100
    bkg_filter_size = (11, 11)
    bkg_edge_method = 'pad'
    bkg_estimator = MedianBackground()

    # sigma clipped stats variables
    sig = 5.
    iters = 5

    # star finder variables
    fwhm = 1.5
    # threshold = 5. * std

    # aperture radius (in pixels)
    # plate scale ~0.034 × 0.030"/pixel
    # infinate aperture is 4"
    # 15 pixels ~ 0.45"
    apr_rad = 5

    # load data and header
    hdu = fits.open(fitsfile)
    data = hdu[0].data
    w = wcs.WCS(hdu[0].header)
    nanmask = np.isnan(data)
    # subtract background
    bkg = Background2D(data, bkg_boxsize, filter_size=bkg_filter_size,
                       edge_method=bkg_edge_method,
                       bkg_estimator=bkg_estimator,
                       mask=nanmask)

    data_sub = data - bkg.background

    # calculate background stats
    mean, median, std = sigma_clipped_stats(data_sub, sigma=sig, iters=iters)
    threshold = 5. * std

    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold, exclude_border=True)
    sources = daofind(data_sub)

    positions = (sources['xcentroid'], sources['ycentroid'])
    #
    apertures = CircularAperture(positions, r=apr_rad)

    effective_gain = hdu[0].header['EXPTIME']

    error = calc_total_error(data_sub, bkg.background, effective_gain)
    phot_table = aperture_photometry(data_sub, apertures, error=error,
                                     mask=nanmask)

    vegamag, vegamagerr = calcmag(phot_table['aperture_sum'],
                                  phot_table['aperture_sum_err'])
    ra, dec = w.all_pix2world(phot_table['xcenter'], phot_table['ycenter'], 1)

    phot_table['ra'] = ra
    phot_table['dec'] = dec
    phot_table['vegamag'] = vegamag
    phot_table['vegamagerr'] = vegamagerr
    fmts = {c: '%.6f' for c in ['xcenter', 'ycenter', 'aperture_sum',
                                'aperture_sum_err', 'ra', 'dec', 'vegamag',
                                'vegamagerr']}
    phot_table.write(outfile, format='ascii.commented_header', formats=fmts,
                     overwrite=True)
    if makeplots:
        outfig = outfile + '.pdf'
        axs = photplots(data, bkg.background, bkg.background_rms, data_sub,
                        fitsfile, outfig)
        apertures.plot(color='r', lw=1.5, alpha=0.5, ax=axs[0, 1])
        plt.savefig(outfig)
    return vegamag, vegamagerr, ra, dec


def calcmag(flux, flerr):
    # flux = flux * APR_COR
    mag = -2.5 * np.log10(flux * PHOTFLAM) - VEGAmag0
    merr = -2.5 * np.log10((flux - flerr) * PHOTFLAM) - VEGAmag0
    perr = -2.5 * np.log10((flux + flerr) * PHOTFLAM) - VEGAmag0
    return mag, (np.abs(mag - perr) + np.abs(mag - merr)) / 2


def photplots(data, bkg_image, bkg_rms, data_sub, fitsfile, outfig):
    from astropy.visualization import LogStretch, SqrtStretch
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
        vegamag, vegamagerr, ra, dec = reduce(fitsfile, outfile, makeplots=True)
        inds = np.isfinite(vegamag)
        h, b = np.histogram(vegamag[inds], bins=np.arange(20, 26, 0.1))
        ax.plot(b[1:], h, ls='steps-', label=fitsfile)
        ax.set_yscale('log')
        if len(args.fitsfiles) > 1:
            ax1.plot(ra, dec, '.', label=fitsfile)
        else:
            cb = ax1.scatter(ra, dec, c=vegamag, cmap=plt.cm.RdBu, label=fitsfile,
                            alpha=0.3, marker='.')
            fig1.colorbar(cb)
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
