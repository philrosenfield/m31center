import numpy as np
import matplotlib.pyplot as plt


def closest_match(num, arr):
    '''index and difference of closet point of arr to num'''
    index = -1
    arr = np.nan_to_num(arr)
    difference = np.abs(num - arr[0])
    for i in range(len(arr)):
        if difference > np.abs(num - arr[i]):
            difference = np.abs(num - arr[i])
            index = i
    return index, difference


def replace_all(text, dic):
    """perfrom text.replace(key, value) for all keys and values in dic"""
    for old, new in dic.items():
        text = text.replace(old, new)
    return text


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
