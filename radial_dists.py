# from match.scripts.graphics.match_plot import match_plot
from mpl_toolkits.axes_grid1 import ImageGrid
from fitshelper.fitshelper.footprints import parse_footprint, parse_poly
import os
import numpy as np
from matplotlib.path import Path
from astropy.io import fits
from shapely.geometry import Polygon


def stitch_cmap(cmap1, cmap2, stitch_frac=0.5, dfrac=0.001, transparent=False):
    '''
    Code adapted from Dr. Adrienne Stilp
    Stitch two color maps together:
        cmap1 from 0 and stitch_frac
        and
        cmap2 from stitch_frac to 1
        with dfrac spacing inbetween

    ex: stitch black to white to white to red:
    stitched = stitch_cmap(cm.Greys_r, cm.Reds, stitch_frac=0.525, dfrac=0.05)
    '''
    from matplotlib.colors import LinearSegmentedColormap

    def left(seg):
        """left color segment"""
        return [(i * (stitch_frac - dfrac), j, k) for i, j, k in seg]

    def right(seg):
        """right color segment"""
        frac = stitch_frac + dfrac
        return [(i * (1 - frac) + frac, j, k) for i, j, k in seg]

    def new_seg(color):
        """combine left and right segments"""
        seg = left(cmap1._segmentdata[color]) + right(cmap2._segmentdata[color])
        return seg

    rgb = ['blue', 'red', 'green']
    cname = '_'.join((cmap1.name, cmap2.name))
    cdict = dict([(key, new_seg(key)) for key in rgb])
    ncmap = LinearSegmentedColormap(cname, cdict, 1024)

    if transparent:
        # set the middle value to zero transparency.
        # it's probably better if you set alpha on the call using the
        # color map rather than change a single value.
        ncmap._init()
        ind = np.max([np.argmax(ncmap._lut.T[i])
                      for i in range(len(ncmap._lut.T)-1)])
        ncmap._lut[ind][-1] = 0
    return ncmap


def zeroed_cmap(hess, cmap1=plt.cm.Reds_r, cmap2=plt.cm.Blues, dfrac=0.05,
                transparent=False):
    """make a diverging color map with white set to 0.0"""
    fhess = hess[np.isfinite(hess)]
    minfhess = np.abs(np.min(fhess))
    # stitch to make a diverging color map with white set to 0.0
    frac = minfhess / (minfhess + np.abs(np.max(fhess)))
    return stitch_cmap(cmap1, cmap2, stitch_frac=frac, dfrac=dfrac,
                       transparent=transparent)


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


def inside_polygons(reg_name, ra, dec, return_areas=False):
    '''
    list of indices of ra, dec array that are within the bounds of a list of
    polygons in the ds9 formatted file reg_name
    '''
    # load regions
    if os.path.isfile(reg_name):
        feet = parse_footprint(reg_name)
        verts = [parse_poly(f.split(';')[1]) for f in feet]
    else:
        print('{} not found.'.format(reg_name))
        return

    # make masks
    radec = np.column_stack((ra, dec))
    masks = [Path(v).contains_points(radec) for v in verts]

    # make masks into shells, get their indices
    inds = [np.nonzero(masks[i] ^ masks[i+1])[0] for i in range(len(masks)-1)]

    retv = inds
    if return_areas:
        areas = [Polygon(radec[i]).area for i in inds]
        retv = inds, areas

    return retv


def hess_by_reg(reg_name, fitsfile, filter1, filter2, extent, dcol=0.1,
               dmag=0.15, ymag='I'):

    uvis = fits.getdata(fitsfile)

    mag1 = uvis[filter1]
    mag2 = uvis[filter2]
    color = mag1 - mag2

    mag = mag1
    if ymag.upper() == 'I':
        mag = mag2

    inds, areas = inside_polygons(reg_name, uvis['ra'], uvis['dec'],
                                  return_areas=True)

    unit_density =  np.array([len(i) for i in inds]) / areas

    (hesses, cbin, mbin) = zip(*[hess(color[i], mag[i], dx=dmag, dy=dcol,
                                      extent=extent) for i in inds])

    scaled_hess = [hesses[i] * unit_density[i] for i in range(len(hesses))]

    return hesses, scaled_hess, cbin, mbin


def hess_panels(hesses, extent, vmin=None, vmax=None, log=False, xlab=None,
               ylab=None, cmaps=None):
    if cmaps is None:
        cmaps = [plt.get_cmap(plt.rcParams['image.cmap'])] * len(hesses)
    if log:
        hesses = np.log10(hesses)

    fig = plt.figure(figsize=(len(hesses)*1.7, 4))
    grid = ImageGrid(fig, 111, nrows_ncols=(1, len(hesses)), label_mode='L',
                     share_all=True, cbar_location='top', cbar_mode='single')
    for i in range(len(hesses)):
        img = grid[i].imshow(hesses[i].T, extent=extent, vmin=vmin, vmax=vmax,
                             cmap=cmaps[i])
    grid.cbar_axes[0].colorbar(img)

    if xlab is not None:
        grid.axes_llc.set_xlabel(xlab)
    if ylab is not None:
        grid.axes_llc.set_ylabel(ylab)
    return grid


def cmd_labels(filter1, filter2, ymag='I', fextra='_VEGA'):
    yfilter = filter1
    if ymag.upper() == 'I':
        yfilter = filter2

    xlab = (r'$\rm{{{0:s}-{1:s}}}$'.format(filter1, filter2)).replace(fextra, '')
    ylab = (r'$\rm{{{0:s}}}$'.format(yfilter)).replace(fextra, '')
    return xlab, ylab


reg_name = '/Users/rosenfield/Work/m31_center/R12_codes/contours_v3.reg'
fitsfile = '/Users/rosenfield/Desktop/pyraf_workdir/13710_M31-CENTER.gst.fits'
filter2 = 'F336W_VEGA'
extent = [-2, 5, 18, 27]

for filter1 in ['F275W_VEGA', 'F225W_VEGA']:
    xlab, ylab = cmd_labels(filter1, filter2)
    kw = {'xlab': xlab, 'ylab': ylab, 'extent': extent}

    hesses, scaled_hess, cbin, mbin = hess_by_reg(reg_name, fitsfile, filter1,
                                                  filter2, extent,
                                                  dcol=0.2, dmag=0.2)
    hess_panels(hesses, log=True, **kw)

    lscaled_hess = np.log10(scaled_hess)
    for i in range(len(hesses)):
        lscaled_hess[i][np.isneginf(lscaled_hess[i])] = 0

    for i in range(len(hesses)-1):
        h1 = hesses[i]
        h2 = hesses[i+1]
        dif = h1 - h2
        kw['cmaps'] = [plt.cm.Blues, plt.cm.Reds, zeroed_cmap(dif)]
        hess_panels([h1, h2, dif], log=False, **kw)
