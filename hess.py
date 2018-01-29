import matplotlib.pyplot as plt
import numpy as np

from .footprints import inside_polygons


def hess(x, y, xbin=None, ybin=None, dx=0.05, dy=0.1, extent=None):
    """Artisan/bespoke np.histogram2d wrapper"""
    if extent is None:
        extent = [x.min(), x.max(), y.min(), y.max()]
    [xmin, xmax, ymin, ymax] = extent

    if xbin is None:
        xbin = np.arange(xmin, xmax, dx)

    if ybin is None:
        ybin = np.arange(ymin, ymax, dy)

    return np.histogram2d(x, y, bins=[xbin, ybin])


def hess_by_reg(reg_name, fitsfile, filter1, filter2, extent=None, dcol=0.1,
                dmag=0.15, index_mag=1, ymag=2, inpoly_kw=None):
    """
    Create hess diagrams for each region in an region file.

    Parameters
    ----------
    reg_name : str
        ds9 formatted region filename

    fitsfile: binary fits table
        input data to make hess diagrams

    filter1, filter2 : str
        column names of the fitsfile

    extent, dcol, dmag:
        extent, dx, dy sent to hess

    ymag : 1 or 2:
        select filter for yaxis of hess.

    index_mag : 1 or 2
        mask the spatial extent of the data by recovered stars in this filter.

    inpoly_kw : dict
        passed to inside_polygons

    Returns
    -------
    hesses, cbin, mbin :
        return of np.histogram 2d

    scaled_hess : Nx2 array
        hesses divided by the area of the region (arcsec^2)
    """
    from astropy.io import fits
    inpoly_kw = inpoly_kw or {}
    if isinstance(fitsfile, fits.fitsrec.FITS_rec):
        uvis = fitsfile
    else:
        uvis = fits.getdata(fitsfile)

    ra = uvis['ra']
    dec = uvis['dec']
    mag1 = uvis[filter1]
    mag2 = uvis[filter2]
    color = mag1 - mag2

    mag = mag1
    if ymag == 2:
        mag = mag2

    imag = filter2
    if index_mag == 1:
        imag = filter1

    detections = get_detections(uvis, filter1, filter2, err=1000.)

    inds, areas = inside_polygons(reg_name, ra[detections], dec[detections],
                                  return_areas=True, **inpoly_kw)
    # degrees^2 to arcsec^2
    areas = np.array(areas) * 3600. * 3600.

    (hesses, cbin, mbin) = zip(*[hess(color[i], mag[i], dx=dmag, dy=dcol,
                                      extent=extent) for i in inds])
    cmds = [np.column_stack([color[detections][i], mag[detections][i]])
            for i in inds]

    # number per arcsec
    # scaled_hess = [hesses[i] / areas[i] for i in range(len(areas))]

    return hesses, areas, cbin, mbin, cmds


def setup_panels(npanels, colorbar=True):
    from mpl_toolkits.axes_grid1 import ImageGrid
    fig = plt.figure(figsize=(npanels * 1.7, 4))
    kw = {'nrows_ncols': (1, npanels),
          'label_mode': 'all',
          'share_all': True}
    if colorbar:
        kw.update({'cbar_location': 'top', 'cbar_mode': 'each'})
    grid = ImageGrid(fig, 111, **kw)
    return fig, grid


def panels(hesses, extent, vmin=None, vmax=None, log=False, xlab=None,
           ylab=None, cmaps=None, cmds=None, plt_kw=None, imshow_kw=None,
           contour_kw=None, contour=False):
    """ImageGrid caller"""
    plt_kw = plt_kw or {}
    imshow_kw = imshow_kw or {}
    contour_kw = contour_kw or {}

    if cmaps is None:
        cmaps = plt.get_cmap(plt.rcParams['image.cmap'])

    if not isinstance(cmaps, list):
        cmaps = [cmaps] * len(hesses)

    if log:
        hesses = np.log10(hesses)

    colorbar = False
    if contour or cmds is None:
        colorbar=True
    fig, grid = setup_panels(len(hesses), colorbar=colorbar)

    for i in range(len(hesses)):
        if cmds is None:
            origin = 'upper'
            img = grid[i].imshow(hesses[i].T, extent=extent,
                                 cmap=cmaps[i], **imshow_kw)
            grid.cbar_axes[i].colorbar(img)
        else:
            origin = 'lower'
            img = grid[i].plot(cmds[i][:, 0], cmds[i][:, 1], '.', **plt_kw)

        if np.sum(hesses[i]) > 0 and contour:
            img = grid[i].contour(hesses[i].T, extent=extent,
                                  origin=origin, zorder=100., **contour_kw)
            grid.cbar_axes[i].colorbar(img)

    if xlab is not None:
        grid.axes_llc.set_xlabel(xlab)
    if ylab is not None:
        grid.axes_llc.set_ylabel(ylab)
    return grid
