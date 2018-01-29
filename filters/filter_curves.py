

def emboss(fg='w', lw=3):
    from matplotlib.patheffects import withStroke
    myeffect = withStroke(foreground=fg, linewidth=lw, alpha=0.5)
    return dict(path_effects=[myeffect])


def discrete_colors(Ncolors, cmap=None):
    '''
    returns list of RGBA tuples length Ncolors

    Parameters
    ----------
    Ncolors: int
        number of discrete colors

    cmap: colormap, optional (defaults to rcParams['image_cmap'])
        colormap to sample

    Returns
    -------
    values: list
        sequence of color values
    '''
    if cmap is None:
        from matplotlib import rcParams
        cmap = rcParams['image.cmap']
    cmap = plt.cm.get_cmap(cmap)
    return [cmap(1. * i / Ncolors) for i in range(Ncolors)]


fd = {'galex_fuv.dat': (1500, 0.4),
      'sbc_F140LP.dat': (1580, 0.03),
      'galex_nuv.dat': (2295, 0.4),
      'uvis1_f225w.dat': (2194, 0.09),
      'uvis1_f275w.dat': (2426, 0.13),
      #'hrc_F330W.dat': (3390, 0.05),
      'uvis1_f336w.dat': (3273, 0.20)}#,
      #'hrc_F475W.dat': (4700, 0.05)}

colors = discrete_colors(8, plt.cm.tab10)

fig, ax0 = plt.subplots()
ax1 = ax0.twinx()
for i, f in enumerate(fd.items()):
    w, t = np.loadtxt(os.path.join('filters', f[0]), unpack=True)
    n = f[0].replace('_', '\n').replace(
        '.dat', '').replace('uvis1', 'uvis').upper()
    ax = ax0
    if 'galex' in f[0]:
        ax = ax1
    ax.plot(w, t, label=n, color=colors[i], zorder=1000 + i)
    ax.text(*f[1], n, color=colors[i], ha='center', zorder=10 - i, **emboss())
ax0.set_xlim(1200, 4000)
ax0.set_ylim(0, 0.26)
ax1.set_ylim(0, 0.65)
ax1.tick_params(labelleft=False, labelright=False, right=False, left=False)
ax0.tick_params(labelleft=False, labelright=False, right=False, left=False)
ax0.set_ylabel('Throughput')
ax0.set_xlabel('Wavelength (Angstroms)')
