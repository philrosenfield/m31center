import numpy as np
import pandas as pd
from .utils import closest_match


def add_extinctions(data, logg=4.5, Rv=3.1, filters=None, extfn=None):
    filters = filters or ['F225W', 'F275W', 'F336W', 'F140LP']

    extfn = extfn or 'models/ext_UVbright.dat'
    ext = np.genfromtxt(extfn, names=True)
    afilts = ['A{}'.format(filt) for filt in filters]
    test = [afilt in ext.dtype.names for afilt in afilts]
    if False in test:
        print('Warning, extinction table does not support filters')
        return data
    try:
        approx_te = np.round(10 ** data[logt], decimals=-3)
        unq_te = np.unique(ext['Teff'])
        (matched_inds, diffs) = zip(*[closest_match(t, unq_te) for t in approx_te])
        inds = np.array([np.nonzero((ext['Teff'] == unq_te[i]) &
                                    (ext['logg'] == logg) &
                                    (ext['Rv'] == Rv))[0]
                         for i in matched_inds]).ravel()
        for afilt in afilts:
            data[afilt] = ext[afilt][inds]
    except ValueError:
        for afilt in afilts:
            data[afilt] = ext[afilt]
    return data



def add_apparent_mag(data, dmod=24.47, Av=0.206, filters=None):
    filters = filters or ['F225W', 'F275W', 'F336W', 'F140LP']
    for filt in filters:
        try:
            Al = data['A{}'.format(filt)]
        except KeyError:
            Al = 1.
        data['m{}'.format(filt)] = data[filt] + dmod + (Av * Al)
    return data


def get_detections(uvis, filter1, filter2, rec=40, err=0.15):
    """
    Find the indices of the data that have mag < rec and mag_err < err
    in both filters
    """
    err1 = filter1.replace('VEGA', 'ERR')
    err2 = filter2.replace('VEGA', 'ERR')
    inds, = np.nonzero((uvis[filter1] < rec) & (uvis[filter2] < rec) &
                       (uvis[err1] < err) & (uvis[err2] < err))
    # print('{} and {} have {} stars'.format(filter1, filter2, len(inds)))
    return inds
