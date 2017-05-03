import os
import sys
import data
import calcs
import graphics
import file_io
#import PHATDataUtils
import text
import diagnostics
import numpy as np
import CrowdTest

def main(default_file):
    home = os.getcwd()
    kwargs = file_io.load_defaults(default_file)
    os.chdir(kwargs['plot_dir'])
    print 'plots are written to',kwargs['plot_dir']
    
    # plot font sizes
    kwargs['ftsize'] = 20
    
    # bin colors
    kwargs['cols'] = ['#007929', '#0D56A6', '#4c063e', '#FF2800', 'black']

    # load bulge photometry file
    kwargs['sup_phot'] = file_io.read_m31_bulge_phot(**kwargs)
    
    # load integrated photometry file
    kwargs['int_flux'] = file_io.load_integrated_phot()
        
    # do completeness tests -- once done, don't need to do them each time.
    # run ASTs.py, update defaults.dat with comp_data to read.
    kwargs['mag1_cut'],kwargs['mag2_cut'] = data.get_mag_cuts(**kwargs)

    # load bulge file
    single = kwargs.get('single')
    bulge_fit = file_io.get_bulge_file(st=bool(abs(single-1)))
    phot_cat = os.path.split(bulge_fit)[1]

    # load contours
    from params import PHAT_DIR_EXTRA3
    contour = PHAT_DIR_EXTRA3 + 'contours_v3.reg'
    
    # areas
    kwargs['areas'] = calcs.area_within_annulus(contour)*(3600.**2)
    
    # overlap regions?
    overlap = kwargs.get('overlap')
    
    # read photometry
    if single==True:
        if overlap == True:
            phot = file_io.read_overlap_data()
            sd1, sd2, both, comp, good = data.photometry_cuts(phot,**kwargs)
            print 'using overlap areas'
            kwargs['areas'] = array([  590.74062232,  2704.71189323,  2911.36133017,  5783.07334901, 9213.35474601])
        else:
            phot = file_io.read_data_diy(bulge_fit)
            sd1, sd2, both, comp, good = data.photometry_cuts(phot,**kwargs)
    else:
        phot = file_io.read_data_st(bulge_fit)
        # sd1 and sd2 are just completeness in mag1, and mag2, both==comp.
        sd1,sd2,both,comp,good = data.completeness_cut(phot,**kwargs)
    
    # divide stars into radial bins [0] = inner [-1] = outter
    rad_inds = calcs.points_inside_ds9_polygon(contour,phot['ra'],phot['dec'])
    
    phot_tests = kwargs.get('phot_tests')
    if phot_tests == True:
        # make photometry plots -- once done, don't need to do them each time
        # can't do them with the full brick (not diy)
        data.make_dolphot_cut_plots(phot,sd1, sd2, both, comp, good, rad_inds,
                                   phot_cat=phot_cat,**kwargs)
        data.make_dolphot_residual_plot(phot, both, rad_inds,
                                   phot_cat='',**kwargs)
    
    # make cmd plot (distances is a kwarg)
    # with completeness lines:
    axs = graphics.five_cmds_complines(phot, rad_inds, both, **kwargs)
    
    # without completeness lines (use good for stars in data)
    # axs = graphics.five_cmds(phot, rad_inds, both, alpha=0.2,**kwargs)
    # axs = graphics.five_cmds(phot, rad_inds, good, **kwargs)
    #sys.exit()

    # make table 1
    kwargs['distances'] = text.make_table1(phot, rad_inds, comp, comp,**kwargs)

    # load padova tracks
    btrack, btrackmore, dtrack, bmstrack = file_io.load_all_tracks(loud=0)
    kwargs['btrack'] = btrack
    kwargs['btrackmore'] = btrackmore
    kwargs['dtrack'] = dtrack
    kwargs['bmstrack'] = bmstrack

    # make cmds of each region
    from params import PHAT_DATA_BASE
    kwargs['disk_fit'] = PHAT_DATA_BASE+'uv-brick09-v4st.fits'
    #calcs.cmds_and_tracks(phot,rad_inds,good,**kwargs)
    
    # get number of each population and the fluxes (as a fn of radius)
    numbers_fluxes = calcs.get_HPHB_PAGB_separation(phot,rad_inds,good,
                                                     **kwargs)
    nhphbs = numbers_fluxes[0]
    npagbs = numbers_fluxes[1]
    flux1_hphbs = numbers_fluxes[2]
    flux2_hphbs = numbers_fluxes[3]
    flux1_pagbs = numbers_fluxes[4]
    flux2_pagbs = numbers_fluxes[5]
    
    calcs.integrated_light_data(*numbers_fluxes,**kwargs)
    # make LF plot, model normalization is done by hard coding in calcs.py
    
    calcs.LF(phot, both, both, rad_inds, **kwargs)    

    # bin by bin comparison
    #CrowdTest.lf_bin_comparison(phot, both, both, rad_inds, **kwargs) 

    # make contributions from hphb and pagb plot
    calcs.hphb_pagb_comp_plot(nhphbs, npagbs, flux1_hphbs, flux2_hphbs,
                        flux1_pagbs, flux2_pagbs,**kwargs)
    
    # make FCT plot
    calcs.time_spent(npagbs,nhphbs,**kwargs)
    
    # make surface density plot
    calcs.numbers_of_pops(npagbs, nhphbs,**kwargs)
    
    # make table 2
    text.make_table2(phot,rad_inds,good,good,nhphbs,npagbs,flux1_hphbs,flux2_hphbs,ff=1e15,
                     **kwargs)
    
    # make hrd plot
    graphics.HRD_Figure(**kwargs)
    
    time_on_cmd = calcs.time_spent_on_cmd(**kwargs)
    #calcs.test_other_tracks(**kwargs)
    #print 'time_on_cmd mean,med:',np.mean(time_on_cmd),np.median(time_on_cmd)
    # do diagnostics -- this is no longer in the paper? (2/10/12)
    #diagnostics.disk_cont(phot,rad_inds,good,**kwargs)
    #f1,f2 = diagnostics.red_star(phot,rad_inds,comp,**kwargs)
    #red_cont1 = np.max(np.array(f1)/(flux1_hphbs+flux1_pagbs))
    #red_cont2 = np.max(np.array(f2)/(flux2_hphbs+flux2_pagbs))    
    #print 'contamination from red stars: max flux/totflux\n'
    #print 'f275w:%f f336w:%f'%(red_cont1,red_cont2)
    os.chdir(home)

def for_talk(defaul_file):
    kwargs = file_io.load_defaults(default_file)
    btrack, btrackmore, dtrack, bmstrack = file_io.load_all_tracks()
    kwargs['btrack'] = btrack
    kwargs['btrackmore'] = btrackmore
    kwargs['yfilter'] = kwargs['Filter1']
    # plot font sizes
    kwargs['ftsize'] = 20
    
    # bin colors
    kwargs['cols'] = ['#007929', '#0D56A6', '#4c063e', '#FF2800', 'black']

    # load bulge photometry file
    kwargs['sup_phot'] = file_io.read_m31_bulge_phot(**kwargs)
    
    kwargs['mag1_cut'],kwargs['mag2_cut'] = data.get_mag_cuts(**kwargs)

    # load bulge file
    single = kwargs.get('single')
    bulge_fit = file_io.get_bulge_file(st=bool(abs(single-1)))
    
    # load contours
    from params import PHAT_DIR_EXTRA3
    contour = PHAT_DIR_EXTRA3 + 'contours_v3.reg'
    
    # areas
    kwargs['areas'] = calcs.area_within_annulus(contour)*(3600.**2)
    
    # read photometry
    if single==True:
        phot = file_io.read_data_diy(bulge_fit)
        sd1, sd2, both, comp, good = data.photometry_cuts(phot,**kwargs)

    kwargs['data_mag1'] = phot['mag1'][good]
    kwargs['data_mag2'] = phot['mag2'][good]
    kwargs['withdata'] = False
    graphics.HRD_Figure(**kwargs)
    rad_inds = calcs.points_inside_ds9_polygon(contour,phot['ra'],phot['dec'])
    kwargs['nomodel'] = True
    calcs.LF(phot, both, both, rad_inds, **kwargs) 
    kwargs['annotate'] =False
    calcs.one_cmd_and_track(phot,rad_inds,good,**kwargs)
    graphics.HRD_Figure_one_by_one(**kwargs)
    numbers_fluxes = calcs.get_HPHB_PAGB_separation(phot,rad_inds,good,
                                                     **kwargs)
    nhphbs = numbers_fluxes[0]
    npagbs = numbers_fluxes[1]
    flux1_hphbs = numbers_fluxes[2]
    flux2_hphbs = numbers_fluxes[3]
    flux1_pagbs = numbers_fluxes[4]
    flux2_pagbs = numbers_fluxes[5]
    
    kwargs['noerror'] = True
    calcs.time_spent(npagbs,nhphbs,**kwargs)
    
if __name__=="__main__":
    default_file = sys.argv[1]
    main(default_file)
    #for_talk(default_file)