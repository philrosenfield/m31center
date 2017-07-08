import os
from collections import OrderedDict

import numpy as np
from shapely.geometry import Polygon

from .utils import replace_all


def inside_polygons(reg_name, ra, dec, return_areas=False, plot=False, ax=None,
                    plt_kw=None):
    '''
    list of indices of ra, dec array that are within the bounds of a list of
    polygons in the ds9 formatted file reg_name
    '''
    if isinstance(reg_name, list):
        verts = reg_name
    else:
        # load regions
        verts = parse_footprint(reg_name, return_coords=True)

    # make masks
    radec = np.column_stack((ra, dec))
    masks = [Path(v).contains_points(radec) for v in verts]

    # make masks into shells, get their indices
    inds = [np.nonzero(masks[i] ^ masks[i + 1])[0]
            for i in range(len(masks) - 1)]

    retv = inds
    if return_areas:
        areas = []
        for i in inds:
            if len(radec[i]) < 3:
                areas.append(1.)
            else:
                areas.append(Polygon(radec[i]).area)
        retv = inds, areas

    # Annoying to pass ax and then also need to set plot=True...
    if plot or ax is not None or plt_kw is not None:
        plt_kw = plt_kw or {}
        if ax is None:
            fig, ax = plt.subplots()
        [ax.plot(ra[i], dec[i], '.', **plt_kw) for i in inds]

    return retv


def find_footprint(fitsfile, hdrext=1, closed=True):
    """
    wrapper to astropy.wcs.calc_footprint with option to append first coord
    to close the footprint.
    """
    from astropy.io import fits
    from astropy import wcs
    print(fitsfile)
    hdu = fits.open(fitsfile)
    try:
        w = wcs.WCS(hdu[hdrext].header)
    except AssertionError:
        hdrext += 1
        w = wcs.WCS(hdu[hdrext].header)
    foot = w.calc_footprint(hdu[hdrext].header)
    if closed:
        foot = np.vstack((foot, foot[0]))
    return foot


def merge_polygons(polygondic):
    pdict = {k: p for k, p in polygondic.items() if p.is_valid is True}
    first = list(pdict.keys())[0]
    ply = pdict[first]
    for i in pdict.keys():
        if i == first:
            continue
        ply = ply.union(pdict[i])
    ply2 = ply.convex_hull
    return ply2, list(pdict.keys())


def split_polygons(polygonlist, tol=49.):
    """
    Return an OrderedDict of polygons that intersect with the first value and
    an OrderedDict of polygons that do not interesect with the first value.
    """
    ins = OrderedDict()
    outs = OrderedDict()
    if len(polygonlist) > 0:
        first = list(polygonlist.keys())[0]
        ply0 = polygonlist[first].convex_hull
        ins = {first: ply0}
        for i in polygonlist.keys():
            if i == first:
                continue
            ply = polygonlist[i].convex_hull
            if ply0.intersects(ply):
                olap = ply0.intersection(ply).area / ply.area * 100
                if olap > tol:
                    ins[i] = ply
                else:
                    # print(olap)
                    outs[i] = ply
            else:
                outs[i] = ply
    return ins, outs


def group_polygons(polylist):
    """
    group a list of Polygons into ones that intersect with eachother
    returns a list of dictionaries with keys maintaining order of list index
    and values being the associated polygon
    """
    npolys = len(polylist)
    # outs = polylist
    outs = {i: p for i, p in enumerate(polylist)}
    groups = []
    while len(outs) != 0:
        ins, outs = split_polygons(outs)
        if len(ins) > 0:
            groups.append(ins)
    assert npolys == np.sum([len(g) for g in groups]), 'lost polygons'
    return groups


def parse_poly(line, closed=True, return_string=False):
    """
    Convert a polygon into N,2 np array. If closed, repeat first coords at end.
    """
    repd = {'j2000 ': '', 'gsc1 ': '', 'icrs ': '', 'multi': '',
            'polygon': '', ')': '', '(': '', 'other': ''}

    line = replace_all(line.lower(), repd).split('#')[0]
    try:
        # ds9-like format all values separated by ,
        polyline = np.array(line.strip().split(','), dtype=float)
    except:
        try:
            # shapely-like format .to_wkt(): (x0 y0, x1 y1, ...)
            xline = line.strip().replace(' ', ',') \
                                .replace(',,,', ',') \
                                .replace(',,', ',')
            polyline = np.array(xline.strip().split(','), dtype=float)
        except:
            print("Do not know how to parse coords: {}".format(line))
            sys.exit(1)

    if closed:
        if False in polyline[:2] == polyline[-2:]:
            polyline = np.append(polyline, polyline[:2])

    retv = polyline.reshape(len(polyline) // 2, 2)

    if return_string:
        retv = ','.join(['[{:.6f}, {:.6f}]'.format(*c) for c in retv])

    return retv


def write_polygons(poly_list, outfile):
    line = ''
    poly_list = np.asarray(poly_list)

    for poly in poly_list:
        if not isinstance(poly, Polygon):
            poly = Polygon(poly)

        coordstr = ','.join(['{},{}'.format(a, b)
                             for a, b in zip(*poly.exterior.coords.xy)])
        line += 'polygon({})\n'.format(coordstr)

    with open(outfile, 'w') as outp:
        outp.write(line)


def parse_footprint(fname, return_coords=False):
    """
    parse a ds9 linear footprint into a ; separated line:
    filename, polygon, central ra dec
    """
    fmt = '{};{};{}'
    if not os.path.isfile(fname):
        print('{} not found.'.format(fname))
        return

    with open(fname) as inp:
        lines = inp.readlines()
    # filename = fname.replace('_footprint_ds9_linear.reg', '.fits')
    polygons = [p.strip().split('#')[0] for p in lines if 'polygon' in p]
    points = [p.split('#')[0].strip() for p in lines if 'point' in p]
    texts = [p.replace('text', 'point').split('#')[0].strip()
             for p in lines if 'text' in p]
    coords = points
    if len(polygons) != len(coords):
        if len(polygons) != len(texts):
            # rolling our own...
            pcoords = [parse_poly(p) for p in polygons]
            polys = np.array([Polygon(c) for c in pcoords])
            pfmt = 'point {:.6f} {:.6f}'
            coords = [pfmt.format(p.centroid.x, p.centroid.y) for p in polys]
        else:
            coords = texts

    assert len(coords) == len(polygons), 'mismatch'
    retv = pcoords
    if not return_coords:
        retv = [fmt.format(fname, polygons[i], coords[i])
                for i in range(len(coords))]
    return retv
