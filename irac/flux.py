import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.path import Path
from shapely.geometry import Polygon

regname = 'irac/contours_v3.reg'
fitsname = 'irac/I1_min_cutout_2734.fits'
# fits_name = 'irac/I1_maj_cutout_2752.fits'


def get_areas(verts):
    """
    Calculate the annular area of a list of concentric rings
    Note: assumes inner most ring is first.
    """
    areas = []
    for v in verts:
        pv = Polygon(v)
        areas.append(pv.area)
    # go from out to in
    areas = areas[::-1]
    a = [areas[i] - areas[i + 1] for i in range(len(areas) - 1)]
    # append inner most ring
    a.append(areas[-1])
    # convert from degrees^s to arcsec
    return (a[::-1] * u.degree ** 2).to(u.arcsecond ** 2)


def get_flux(fitsname, regname):
    """
    Calculate the flux of an image within each annular region.
    """
    hdu = fits.open(fits_name)
    header = hdu[0].header
    w = WCS(header)

    # make sure we got the right units.
    q = u.MJy / u.sr
    assert q.to_string().replace(' ', '') == header['BUNIT']

    data = hdu[0].data
    # I don't know how to get an index matching with data.
    # so here is a grid of indices, one for each pixel.
    grid = np.meshgrid(np.arange(data.shape[0]), np.arange(data.shape[1]))
    grid = np.concatenate(grid).ravel()
    grid = grid.reshape(2, len(grid) // 2).T

    # parse footprint and get the area of the annuli.
    verts = parse_footprint(regname, return_coords=True)
    areas = get_areas(verts)
    # Convert the verts to X,Y since its easier than converting the image to
    # ra dec.
    xys = [np.column_stack(w.wcs_world2pix(v[:, 0], v[:, 1], 0))
           for v in verts]

    # grab the pixels in each region
    # note the order of the regions matter.
    # the usual file has the inner first and outter last.
    all_inds = []
    for i, xy in enumerate(xys[::-1]):
        inds, = np.nonzero(Path(xy).contains_points(grid))
        all_inds.append(inds)

    # cut out the inner annuli
    inds_within_annuli = [list(set(all_inds[i]) - set(all_inds[i + 1]))
                          for i in range(len(all_inds) - 1)]
    # don't forget the innermost (that is a filled ellipse)
    inds_within_annuli.append(all_inds[-1])

    counts_reg = []
    flux_reg = []
    # note, we're going back to inner to outer...
    for i, inds in enumerate(inds_within_annuli[::-1]):
        # sum the counts
        counts = np.sum([data[i, j] for i, j in grid[inds]])
        # and convert to arcsecond^2
        counts = (counts * q).to(u.Jy / (u.arcsecond ** 2))
        # integrated flux
        flux = counts * areas[i]

        flux_reg.append(flux)
        counts_reg.append(counts)

        # fig, ax = plt.subplots()
        # ax.plot(grid[:, 0][inds], grid[:, 1][inds], '.', alpha=0.3)
        # print(counts, flux, areas[i])
    return flux_reg, counts_reg, areas
