# from match.scripts.graphics.match_plot import match_plot
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

from .footprints import parse_footprint
from .hess import hess_by_reg, hess_panels
from .stev import brightest_tracks, load_models
from .utils import zeroed_cmap

# Things to be in the function caller:
# base: '/Users/rosenfield/Library/Mobile Documents/com~apple~CloudDocs/Work/m31center
#reg_name = 'irac/contours_masked.reg'
reg_name = 'irac/contours_v3.reg'
fitsfile = 'data/13710_M31-CENTER_F225W_F275W_F336W_F140LP.gst.fits'
filter2 = 'F336W_VEGA'

uvis = fits.getdata(fitsfile)
verts = parse_footprint(reg_name, return_coords=True)

#inpoly_kw = None
inpoly_kw = {'ax': ax, 'plt_kw': {'alpha': 0.3}}
hbyreg_kw = {'dcol': 0.2, 'dmag': 0.2, 'inpoly_kw': inpoly_kw, 'index_mag': 2}

fig, ax = plt.subplots()
# hepnfn = 'R12_codes/vw_extended/He-pnn.tab.UVbright'
# data = load_model('R12_codes/vw_extended/H-pnn.tab.UVbright', pne=True)
data = load_models(models)
# the brightest track for each metallicity
df = brightest_tracks(data, filter2)

for i, filter1 in enumerate(['F225W_VEGA', 'F140LP_VEGA']):
    xlab, ylab = cmd_labels(filter1, filter2)
    if i < 2:
        extent = [-2., 5, 16, 26.5]
    else:
        extent = [-2.5, 3, 16, 26.5]

    hesses, areas, _, _, cmds = hess_by_reg(verts, uvis, filter1, filter2,
                                            extent, **hbyreg_kw)

    for i, cmd in enumerate(cmds):
        inds = Path(hphb_verts).contains_points(cmd)
        xinds = ~Path(hphb_verts).contains_points(cmd)
        print(len(np.nonzero(inds)[0]) / areas[i])
        print(' ', len(np.nonzero(xinds)) / areas[i])

    kw = {'xlab': xlab, 'ylab': ylab, 'extent': extent, 'cmaps': plt.cm.Blues,
          'plt_kw': {'ms':0.8}, 'contour_kw': {'levels': np.linspace(5,60,10)}}
    grid = panels(hesses, cmds=cmds, **kw)
    #grid = panels(hesses, **kw)

    for ax in grid:
        mf1 = ('m{}'.format(filter1)).replace('_VEGA', '')
        mf2 = ('m{}'.format(filter2)).replace('_VEGA', '')
        ax.plot(d[mf1] - d[mf2], d[mf1])
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:][::-1])
    # sgrid = hess_panels(scaled_hess, **kw)


    for i in range(len(hesses) - 1):
        h1 = hesses[i]
        h2 = hesses[i + 1]
        dif = h1 - h2
        kw['cmaps'] = [plt.cm.Blues, plt.cm.Reds, zeroed_cmap(dif)]
        hess_panels([h1, h2, dif], log=False, **kw)
