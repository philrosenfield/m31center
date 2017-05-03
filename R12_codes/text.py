import numpy as np
import Astronomy
import file_io

def make_table2(phot,rad_inds,cut_inds1,cut_inds2,
                nhphbs,npagbs,flux1_hphbs,flux2_hphbs,ff=1e15,**kwargs):
    
    Filter1 = kwargs.get('Filter1')
    Filter2 = kwargs.get('Filter2')
    areas = kwargs.get('areas')
    sup_phot = kwargs.get('sup_phot')
    rif_36 = sup_phot['rif_36']
    sb_36 = sup_phot['rasb_36']

    # Bulge fraction B/(B+D) 
    frac_bulge = sup_phot['bulge_fraction'] 
    BtoD = frac_bulge/(1.-frac_bulge)

    # integrated flux per bin
    int_flux = kwargs.get('int_flux')
    int_counts1 = int_flux['f275w_cts']
    int_counts2 = int_flux['f336w_cts']
    
    vega1counts = 1.1413e9
    vega1flx = 3.727e-9
    vega2counts = 2.4318e9
    vega2flx = 3.245e-9

    intflx1_cgs = int_counts1*vega1flx/vega1counts
    intflx2_cgs = int_counts2*vega2flx/vega2counts

    # resolved flux per bin
    flux_275w = get_flux_filter(phot,'mag1',rad_inds,cut_inds1,Filter1)
    flux_336w = get_flux_filter(phot,'mag2',rad_inds,cut_inds2,Filter2)

    # number density of hphb and peagb
    number_arcsec2_hphb = nhphbs/areas[1:]
    number_arcsec2_pagb = npagbs/areas[1:]
    
    # flux of hphbs


    # factor
    print 'table 2 flux in units of %g'%(1/ff)
    
    i=0
    print ' 0  & %.1f & %.2f & ... & ... & ... & ... & ... & ...&  ... & ...  \\\ '% (BtoD[i],rif_36[i])
    for i in range(len(areas)-1):
        print ' %i & %.1f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.3f & %.3f \\\ '%(i+1,BtoD[i+1],rif_36[i+1],intflx1_cgs[i]*ff,intflx2_cgs[i]*ff,flux_275w[i]*ff,flux_336w[i]*ff,flux1_hphbs[i]*ff,flux2_hphbs[i]*ff,number_arcsec2_hphb[i],number_arcsec2_pagb[i])
    
    return

def get_flux_filter(phot,mag,rad_inds,cut_inds,filter):
    flux = []
    for i in range(len(rad_inds)):
        inds = list(set(rad_inds[i]) & set(cut_inds))
        fluxs = Astronomy.vegamag2flux(phot[mag][inds],filter,'wfc3uvis')
        flux.append(np.sum(fluxs))
    return flux
            
def get_nstars(phot,mag,rad_inds,cut_inds):
    nstars = []
    for i in range(len(rad_inds)):
        inds = list(set(rad_inds[i]) & set(cut_inds))
        nstars.append(len(phot[mag][inds]))
    return nstars

def make_table1(phot, rad_inds,cut_inds1, cut_inds2,**kwargs):
    # table 1
    print 'table 1'
    areas = kwargs.get('areas')
    sup_phot = kwargs['sup_phot']
    dmod = kwargs.get('dmod')
    
    contour_arcsec = sup_phot['semimajor_axis']
    contour_arcsec_min = np.concatenate(([0],contour_arcsec[:-1]))
    median_arcsec = get_medians(contour_arcsec)
    
    contour_pc = Astronomy.arcsec2kpc(contour_arcsec,dmod=dmod)*1000.
    median_pc = get_medians(contour_pc)
    
    nstars_mag1 = get_nstars(phot,'mag1',rad_inds,cut_inds1)
    nstars_mag2 = get_nstars(phot,'mag2',rad_inds,cut_inds2)
    print 'number of stars:',sum(nstars_mag1)
    comp_data = kwargs.get('comp_data')
    comp50mag1, comp50mag2 = file_io.get_comp50mags(comp_data)
    comp90mag1, comp90mag2 = file_io.get_comp90mags(comp_data)
    
    cont_file = kwargs.get('cont_file')
    contamination = file_io.load_contamination(cont_file)
    i = 0
    #print '0 & %.0f & %.1f & %.1f & %.1f & %.1f & ... & ... & ... & ... & ... \\\ '%(areas[i],contour_arcsec_min[i],contour_arcsec[i],median_arcsec[i],median_pc[i])
    print '0 & %.0f & %.1f & %.1f & %.1f & %.1f & ... & ... & ... & ... & ... \\\ '%(areas[i],contour_arcsec_min[i],contour_arcsec[i],median_arcsec[i],median_pc[i])
    for i in range(len(areas)-1):
        #print '%i & %.0f & %.1f & %.1f & %.1f & %.1f & %i & %.1f & %.1f & %.1f & %.1f \\\ '%(i+1,areas[i+1],contour_arcsec_min[i+1],contour_arcsec[i+1],median_arcsec[i+1],median_pc[i+1],nstars_mag1[i],comp50mag1[i],comp50mag2[i],comp90mag1[i],comp90mag2[i])
        print '%i & %.0f & %.1f & %.1f & %.1f & %.1f & %i & %.1f & %.1f & %.2f & %.2f \\\ '%(i+1,areas[i+1],contour_arcsec_min[i+1],contour_arcsec[i+1],median_arcsec[i+1],median_pc[i+1],nstars_mag1[i],comp90mag1[i],comp90mag2[i],contamination['F275W'][str(i)],contamination['F336W'][str(i)])    
    return median_pc

def get_medians(arr):
    medians = [np.median((arr[i],arr[i+1])) for i in range(len(arr)-1)]
    medians = np.concatenate(([arr[0]/2.],medians))
    return medians
