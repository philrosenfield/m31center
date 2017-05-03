import matplotlib.pyplot as plt
import matplotlib.nxutils as nxutils
import os,sys
import pyfits
import data
from params import *
import graphics
import numpy as np

def radec_to_reg(ra,dec,outfile='radec.reg',shape='circle',size=0.5,header='default',append=False):
    if header == 'default':
        head = 'global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n fk5'
    else:
        head = header
    flag = 'w'
    if append == True: 
        flag = 'a'
        print 'appending',outfile

    f = open(outfile,flag)
    if append != True: f.write('%s'%head)
    if shape == 'circle':
        [f.write('%s( %f, %f,%f")\n'%(shape,ra[i],dec[i],size)) for i in range(len(ra))]
    if shape == 'polygon':
        f.write('%s('%shape)
        line = [''.join('%f,%f'%(r,d)) for r,d in zip(ra,dec)]
        f.write(','.join([str(l) for l in line]))
        f.write(')\n')
    print 'radec_to_reg wrote',outfile
    f.close()

def reg_to_array(reg_file,shape='polygon'):
    lines = open(reg_file,'r').readlines()
    dic = {}
    i = 0
    for line in lines:
        line = line.strip()
        if not line.startswith(shape): continue
        try:
            poly,tag = line.split('# tag=')
        except:
            poly = line
            tag = str(i)
        data = map(float,poly.split('(')[-1].split(')')[0].split(','))
        dic[tag.replace('\'','')] = {'ra': data[::2], 'dec': data[1::2]}
        i+=1
    return dic

def plot_overlap_regions():
    fig = plt.figure()
    ax = plt.axes()
    exposure_reg = os.path.join(PHAT_DATA_BASE,'12058_exposure_coverage_UVIS_F336W_chipswap_bulge.reg')    
    exposure_dict = reg_to_array(exposure_reg)
    for key,dict in exposure_dict.items():
        dict['ra'].append(dict['ra'][0])
        dict['dec'].append(dict['dec'][0])
        ax.plot(dict['ra'],dict['dec'],color='grey')
    
    cols = ['#007929', '#0D56A6', '#4c063e', '#FF2800', 'black']
    v3_contours = PHAT_DIR_EXTRA3+'contours_v3.reg'
    contour_dict = reg_to_array(v3_contours)
    i=0
    for key,dict in contour_dict.items():
        if int(key) ==0: 
            ax.plot(dict['ra'],dict['dec'],'--',color='black')
        else:
            ax.plot(dict['ra'],dict['dec'],color=cols[int(key)-1],lw=2)
    verts = get_overlap_regions()[-1]
    #ax.plot(verts[:,0],verts[:,1],lw=2,color='darkorange')
    #copy and paste from get_overlap_areas....
    #[ax.plot(mv[:,0],mv[:,1],lw=3,color='cyan') for i,mv in enumerate(masked_vertst)]
    #[ax.plot(mv[:,0],mv[:,1],lw=3,color='cyan') for i,mv in enumerate(masked_vertsb)]

    ax.set_xlabel('ra')
    ax.set_ylabel('dec')
    #plt.savefig('overlap_region.png')
    return

def get_overlap_regions(write_reg=False):
    exposure_reg = os.path.join(PHAT_DATA_BASE,'12058_exposure_coverage_UVIS_F336W_chipswap_bulge.reg')
    exposure_dict = reg_to_array(exposure_reg)
    fields = [int(key.split('-')[2].replace('F','')) for key in exposure_dict.keys() if not key.endswith('1')]
    keys =  np.array([key for key in exposure_dict.keys() if not key.endswith('1')])
    fsort = np.sort(fields)
    #order = [3,4,5,6,12,11,10,9]
    ras = []
    decs = []
    order = [5,4,2,3,7,0,1,9]
    for key in keys[order]:
        dict = exposure_dict[key]
        field = int(key.split('-')[2].replace('F',''))
        plt.plot(dict['ra'],dict['dec'],color='grey')
        if field == 3:
            ras.append(dict['ra'][1])
            ras.append(dict['ra'][2])
            decs.append(dict['dec'][1])
            decs.append(dict['dec'][2])
        if 3 < field <= 6:
            ras.append(dict['ra'][3])
            ras.append(dict['ra'][0])
            decs.append(dict['dec'][3])
            decs.append(dict['dec'][0])
        if 9 < field <= 12:
            ras.append(dict['ra'][1])
            ras.append(dict['ra'][2])
            decs.append(dict['dec'][1])
            decs.append(dict['dec'][2])
        if field == 9:
            ras.append(dict['ra'][3])
            ras.append(dict['ra'][0])
            decs.append(dict['dec'][3])
            decs.append(dict['dec'][0])
    
    ras.append(ras[0])
    decs.append(decs[0])
    
    ras2 = []
    decs2 = []
    order2 = [9,1,0,7,11,10,6,8]
    for key in keys[order2]:
        dict = exposure_dict[key]
        field = int(key.split('-')[2].replace('F',''))
        plt.plot(dict['ra'],dict['dec'],color='grey')
        if field == 9:
            ras2.append(dict['ra'][1])
            ras2.append(dict['ra'][2])
            decs2.append(dict['dec'][1])
            decs2.append(dict['dec'][2])
        if 9 < field <= 12:
            ras2.append(dict['ra'][3])
            ras2.append(dict['ra'][0])
            decs2.append(dict['dec'][3])
            decs2.append(dict['dec'][0])
        if 15 < field <= 18:
            ras2.append(dict['ra'][1])
            ras2.append(dict['ra'][2])
            decs2.append(dict['dec'][1])
            decs2.append(dict['dec'][2])
        if field == 15:
            ras2.append(dict['ra'][3])
            ras2.append(dict['ra'][0])
            decs2.append(dict['dec'][3])
            decs2.append(dict['dec'][0])
    
    ras2.append(ras2[0])
    decs2.append(decs2[0])
    
    fields = [int(key.split('-')[2].replace('F','')) for key in exposure_dict.keys() if not key.endswith('1')]
    keys =  np.array(exposure_dict.keys())
    order = [12,3,9,5,1,21]
    ras3 = []
    decs3 =[]
    for key in keys[order]:
        dict = exposure_dict[key]
        field = key.split('-')[2]+'-'+key.split('-')[-1]
        if field == 'F16-2':
            ras3.append(dict['ra'][0])
            decs3.append(dict['dec'][0])        
            ras3.append(dict['ra'][1])
            decs3.append(dict['dec'][1])        
        if field == 'F10-2': 
            ras3.append(dict['ra'][1])
            decs3.append(dict['dec'][1])        
        if field == 'F04-2': 
            ras3.append(dict['ra'][1])
            decs3.append(dict['dec'][1])        
        if field == 'F05-1': 
            ras3.append(dict['ra'][2])
            decs3.append(dict['dec'][2])
            ras3.append(dict['ra'][3])
            decs3.append(dict['dec'][3])
        if field == 'F11-1': 
            ras3.append(dict['ra'][3])
            decs3.append(dict['dec'][3])
        if field == 'F17-1': 
            ras3.append(dict['ra'][3])
            decs3.append(dict['dec'][3])
    ras3.append(ras3[0])
    decs3.append(decs3[0])
    
    verts = np.column_stack((ras,decs))
    verts2 = np.column_stack((ras2,decs2))
    verts3 = np.column_stack((ras3,decs3))
    if write_reg == True:
        for v in [verts,verts2,verts3]:
            radec_to_reg(v[:,0],v[:,1],'overlaps.reg',shape='polygon',append=True)
    
    
    '''
    if field == 'F17-1' or field == 'F16-2' or field == 'F11-1' or field == 'F10-2' or field == 'F05-1' or field == 'F04-2':
        plt.plot(dict['ra'],dict['dec'],color='grey')
        for i in range(len(dict['ra'])):
            plt.plot(dict['ra'][i],dict['dec'][i],'o',color=cols[i])
            if i == 0:
                plt.annotate(field,xy=(dict['ra'][i],dict['dec'][i]))
    for key,dict in contour_dict.items():
        plt.plot(dict['ra'],dict['dec'],color='grey')
    '''
    return verts, verts2, verts3

def get_overlap_areas():
    ast_file = kwargs.get('ast_file')
    ast = ASTs.load_ast_file(ast_file)
    ra = ast['ra']
    dec = ast['dec']
    points = np.column_stack((ra,dec))
    
    v3_contours = PHAT_DIR_EXTRA3+'contours_v3.reg'
    annulus_areas = calcs.area_within_annulus(v3_contours)
    annuli = reg_to_array(v3_contours,shape='polygon')
    aras = [annuli[i]['ra'] for i in np.sort(annuli.keys())]
    adecs = [annuli[i]['dec'] for i in np.sort(annuli.keys())]
    verts = [np.column_stack((aras[i],adecs[i])) for i in range(len(aras))]
    
    center_ra=10.68463
    center_dec=41.26693
    top = np.nonzero(dec>=center_dec)[0]
    bottom =np.nonzero(dec<=center_dec)[0]
    overlaps = reg_to_array('overlaps.reg',shape='polygon')
    mra,mdec = overlaps['2']['ra'],overlaps['2']['dec']
    overlap_verts = np.column_stack((mra,mdec))
    
    overlap_mask = nxutils.points_inside_poly(points, overlap_verts)
    ov_inds = np.nonzero(overlap_mask)[0]
    # much better verts:
    ov_verts = GenUtils.get_verts(ra[ov_inds],dec[ov_inds],nbinsx=300,nbinsy=300)
    
    inds_annuli = calcs.points_inside_ds9_polygon(v3_contours,ra,dec)

    masked_vertst = []
    masked_vertsb = []
    areas = []
    for i,rad_inds in enumerate(inds_annuli):
        tinds = list(set(rad_inds) & set(ov_inds) & set(top))
        binds = list(set(rad_inds) & set(ov_inds) & set(bottom))
        vertst = GenUtils.get_verts(ra[tinds],dec[tinds],nbinsx=50,nbinsy=50)
        vertsb = GenUtils.get_verts(ra[binds],dec[binds],nbinsx=50,nbinsy=50)
        masked_vertst.append(vertst)
        masked_vertsb.append(vertsb)
        areat = GenUtils.area(vertst)
        areab =  GenUtils.area(vertsb)
        areas.append(areat+areab)
    
    # testing
    #[plt.plot(r,d) for r,d in zip(aras,adecs)]
    #[ax.plot(mv[:,0],mv[:,1],lw=3,color='cyan') for i,mv in enumerate(masked_vertst)]
    #[ax.plot(mv[:,0],mv[:,1],lw=3,color='cyan') for i,mv in enumerate(masked_vertsb)]
    # get rid of center exclusion
    areas[0] = areas[0]-annulus_areas[0]
    frac_areas = np.array(areas)/annulus_areas[1:]
    return np.array(areas), frac_areas

def main():
    # /astro/users/ben/ad2/ben/WFC3/12058_M31-B01-F040510111617-UVIS/proc/
    fits_file = os.path.join(PHAT_DATA_BASE,'12058_M31-B01-F040510111617-UVIS_F275W_F336W.st.fits')
    exposure_reg = os.path.join(PHAT_DATA_BASE,'12058_exposure_coverage_UVIS_F336W_chipswap_bulge.reg')
    
    v3_contours = PHAT_DIR_EXTRA3+'contours_v3.reg'
    
    f = pyfits.open(fits_file)
    f_data = f[1].data
    f.close()
    
    Filter1 = 'F275W'
    Filter2 = 'F336W'
    
    mag1 = f_data.field('MAG1_UVIS')
    mag2 = f_data.field('MAG2_UVIS')
    color = mag1-mag2
    ra = f_data.field('RA')
    dec = f_data.field('DEC')
    contour_dict = reg_to_array(v3_contours)
    
    verts, verts2, verts2  = get_overlap_regions()
    
    radec = np.column_stack((ra,dec))
    ind = np.nonzero(nxutils.points_inside_poly(radec, verts))[0]
    ind2 = np.nonzero(nxutils.points_inside_poly(radec, verts2))[0]
    ind3 = np.nonzero(nxutils.points_inside_poly(radec, verts3))[0]
    
    inds = calcs.points_inside_ds9_polygon(v3_contours,ra,dec)
    
    axs = graphics.setup_five_panel_plot(fig=2)
    for i in range(len(inds)):
        bam = list(set(inds[i]) & set(ind3))
        axs[i].plot(color[bam],mag2[bam],'.')
        axs[i].set_ylim(27,19.5)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
