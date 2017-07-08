import os

import pandas as pd

from .color_tracks import add_apparent_mag, add_extinctions

mass = 'mass'
z = 'z'
age = 'age'
logt = 'logt'
logl = 'logl'

pne_header = [mass, 'y', z, age, logt, logl, 'mbol', 'F140LP',
              'F175Wfoc', 'F275Wfoc', 'F200LP', 'F218W', 'F225W', 'F275W',
              'F300X', 'F336W', 'F25QTZ', 'FUV', 'NUV', 'U', 'B', 'V']

track_header = [age, logl, logt, 'j', 'Mode', mass, 'XSUP', 'YSUP', 'XCSUP',
                'XOSUP', 'XNSUP', 'RATE', 'XCEN', 'YCEN', 'XC_CEN', 'XN_CEN',
                'XO_CEN', z, 'co']

uvb_header = ['mbol', 'F140LP', 'F175Wfoc', 'F275Wfoc', 'F200LP', 'F218W',
              'F225W', 'F275W', 'F300X', 'F336W', 'F25QTZ', 'FUV', 'NUV',
              'UVW2', 'UVM2', 'UVW1', 'u', 'U', 'B', 'V']

hrc_header = ['mbol', 'F220W', 'F250W', 'F330W', 'F344N', 'F435W', 'F475W',
              'F550M', 'F555W', 'F606W', 'F625W', 'F658N', 'F660N', 'F775W',
              'F814W', 'F850LP', 'F892N']

models = ['parsec_v1.2s_Z0.0001_Y0.249.HB.dat',
          'parsec_v1.2s_Z0.0002_Y0.249.HB.dat',
          'parsec_v1.2s_Z0.0005_Y0.249.HB.dat',
          'parsec_v1.2s_Z0.001_Y0.25.HB.dat',
          'parsec_v1.2s_Z0.002_Y0.252.HB.dat',
          'parsec_v1.2s_Z0.004_Y0.256.HB.dat',
          'parsec_v1.2s_Z0.006_Y0.259.HB.dat',
          'parsec_v1.2s_Z0.008_Y0.263.HB.dat',
          'parsec_v1.2s_Z0.014_Y0.273.HB.dat',
          'parsec_v1.2s_Z0.017_Y0.279.HB.dat',
          'parsec_v1.2s_Z0.01_Y0.267.HB.dat',
          'parsec_v1.2s_Z0.02_Y0.284.HB.dat',
          'parsec_v1.2s_Z0.03_Y0.302.HB.dat',
          'parsec_v1.2s_Z0.04_Y0.321.HB.dat',
          'parsec_v1.2s_Z0.06_Y0.356.HB.dat']

uvbright = [os.path.join('models', 'uvbright', m) for m in models]
acs_hrc = [os.path.join('models', 'acs_hrc', m) for m in models]

hphb_verts = np.array([[-0.44802867,  26.02164819],
                       [-0.43266769,  25.31251606],
                       [-0.41218638,  24.37155226],
                       [-0.38402458,  23.53286713],
                       [-0.32514081,  23.34876552],
                       [-0.17665131,  23.23284969],
                       [0.11264721,  23.28057974],
                       [0.35330261,  23.43740704],
                       [0.62467998,  23.58059718],
                       [0.9421403,  23.89425178],
                       [1.29544291,  24.28291074],
                       [1.78955453,  24.87612705],
                       [1.98412698,  25.48979909]])


def load_models(band='uvbright'):
    return [load_model(os.path.join('models', band, m), band=band) for m in models]


def load_model(filename, dmod=24.47, Av=0.206, logg=4.5, Rv=3.1, pne=False,
               band='uvbright'):
    names = track_header
    filters = ['F225W', 'F275W', 'F336W', 'F140LP']

    if band == 'uvbright':
        names = np.concatenate([track_header, uvb_header])
        extfn = 'models/ext_UVbright.dat'

    if band == 'acs_hrc':
        names = np.concatenate([track_header, hrc_header])
        filters = ['F330W', 'F475W']
        extfn = 'models/ext_acs_hrc.dat'
    if pne:
        names = pne_header

    data = pd.read_table(filename, delim_whitespace=True, names=names,
                         header=0)
    data = add_extinctions(data, logg=logg, Rv=Rv, filters=filters,
                           extfn=extfn)
    data = add_apparent_mag(data, dmod=dmod, Av=Av, filters=filters)
    return data


def select_by_func(data, key, func=np.argmin):
    """
    Select track by function.
    For example:
    brightest: select_by_func(data, filter2, func=np.argmin)
    most massive: select_by_func(data, mass, func=np.argmax)
    """
    return data[data[mass] == data[mass].iloc[func(data[key])]]


def brightest_tracks(data, filt):
    filt = filt.replace('_VEGA', '')
    if not isinstance(data, list):
        datas = [data]
    else:
        datas = data
    return [select_by_func(d, filt) for d in datas]


def cmd_labels(filter1, filter2, ymag='I', fextra='_VEGA'):
    yfilter = filter1
    if ymag.upper() == 'I':
        yfilter = filter2

    xlab = (r'$\rm{{{0:s}-{1:s}}}$'.format(filter1, filter2)
            ).replace(fextra, '')
    ylab = (r'$\rm{{{0:s}}}$'.format(yfilter)).replace(fextra, '')
    return xlab, ylab
