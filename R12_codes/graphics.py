import matplotlib.nxutils as nxutils
import matplotlib.pyplot as plt
from matplotlib import cm,rc
from matplotlib.patches import FancyArrow
from matplotlib.ticker import NullFormatter, MultipleLocator
nullfmt   = NullFormatter() # no labels
from matplotlib import rcParams
#rcParams['text.usetex']=True
#rcParams['text.latex.unicode']=False
rcParams['axes.linewidth'] = 2
rcParams['ytick.labelsize'] = 'large'
rcParams['xtick.labelsize'] = 'large'
rcParams['axes.edgecolor'] = 'grey'
rcParams['axes.labelsize'] = 'large'
rcParams['savefig.dpi'] = 300
rc('text',usetex=True)
from PadovaTracksUtils import *
from GenUtils import get_afile,tablefile2dict,continuous_hist,continuous_hist_manual,mags2polygon,mag2polygon,manual_hist,closest_match,bin_up,smooth
import GenUtils
import file_io
import os
from params import *

vw_loc = os.path.join(PHAT_TRACKS_BASE,'vw_extended/')
ftsize = 20

def load_ann_kwargs():
    from matplotlib.patheffects import withStroke
    myeffect = withStroke(foreground="w", linewidth=3)
    ann_kwargs = dict(path_effects=[myeffect])
    return ann_kwargs

ann_kwargs = load_ann_kwargs()

def error_wedge(ax,wedge_color,x,p_err,m_err,**pltkwargs):
    return ax.fill_between(x,p_err,m_err,facecolor=wedge_color,**pltkwargs)

def set_up_three_panel_plot():
    fig = plt.figure(figsize=(8,8))
    
    left,width= 0.1,0.4
    bottom, height = 0.1, 0.4
    d=0.01
    lefter = left+width+d
    mid = bottom+height+2*d
    lefts = [left,lefter,left]
    bottoms = [mid,mid,bottom]
    widths = [width,width,2*width+d]
    heights = [height,height,height-0.1]
    
    axs = [plt.axes([l,b,w,h]) for l,b,w,h in zip(lefts,bottoms,widths,heights)]
    return axs

# For HRD_Figure, frac_flux, integrated_LF
def two_panel_plot(sizex,sizey,xlab1,xlab2,ylab,ylab2=None,ftsize=20,mag2_cut=0,mag2_max=1,fig=2):
    fig = plt.figure(fig, figsize=(sizex,sizey))
    left, width = 0.1, 0.4
    bottom, height = 0.12, 0.8  
    
    left2 = left+width + 0.065
    if ylab2 != None: left2 = left+width + 0.08
    axis1 = [left, bottom, width, height]
    axis2 = [left2, bottom, width, height]
    
    ax1 = plt.axes(axis1)
    ax1.set_xlim( (mag2_cut,mag2_max) ) # set all axes limits here
    #ax1.set_ylim( (0.0001, 10.) )
    ax1.set_xlabel(r'%s'%xlab1,fontsize=ftsize)
    ax1.set_ylabel(r'%s'%ylab,fontsize=ftsize)
    
    ax2 = plt.axes(axis2)
    ax2.set_xlim( ax1.get_xlim() )
    #ax2.set_ylim( ax1.get_ylim() )
    ax2.set_xlabel(r'%s'%xlab2,fontsize=ftsize)
    if ylab2 !=None: ax2.set_ylabel(r'%s'%ylab2,fontsize=ftsize)
    return ax1,ax2

def two_panel_plot_v2(sizex,sizey,xlab1,xlab2,ylab,ylab2=None,ftsize=20,mag2_cut=0,mag2_max=1,ylabp2=""):
    fig = plt.figure(2, figsize=(sizex,sizey))
    left, width = 0.07, 0.4
    bottom, height = 0.12, 0.78  
    
    left2 = left+width + 0.09
    if ylab2 != None: left2 = left+width + 0.11
    axis1 = [left, bottom, width, height]
    axis2 = [left2, bottom, width, height]
    
    ax1 = plt.axes(axis1)
    ax1.set_xlim( (mag2_cut,mag2_max) ) # set all axes limits here
    #ax1.set_ylim( (0.0001, 10.) )
    ax1.set_xlabel(r'%s'%xlab1,fontsize=ftsize)
    ax1.set_ylabel(r'%s'%ylab,fontsize=ftsize)
    
    ax2 = plt.axes(axis2)
    ax2.set_xlim( ax1.get_xlim() )
    #ax2.set_ylim( ax1.get_ylim() )
    ax2.set_xlabel(r'%s'%xlab2,fontsize=ftsize)
    if ylab2 !=None: ax2.set_ylabel(r'%s'%ylab2,fontsize=ftsize)
    
    #par2 = ax2.twinx()
    #par2.set_ylabel(ylabp2,fontsize=ftsize)
    return ax1,ax2

def two_panel_plot_vert(oney=True,ftsize=20,mag2_cut=0,mag2_max=1):
    fig = plt.figure(2, figsize=(8,8))
    left, width = 0.13, 0.83
    bottom, height = 0.1, 0.41
    dh = 0.03
    
    axis1 = [left, bottom, width, height]
    axis2 = [left, (bottom+height+dh), width, height]
    
    ax1 = plt.axes(axis1)
    ax1.set_xlim( (mag2_cut,mag2_max) ) # set all axes limits here
    #ax1.set_ylim( (0.0001, 10.) )
    if oney==True: ax1.annotate(r'$\#/ 3.6 \mu \rm{m\ Region\ Integrated\ Flux\ (Jy}^{-1}\rm{)}$',fontsize=ftsize,xy=(0.025,.5),xycoords='figure fraction',va='center',rotation='vertical')
    ax1.set_xlabel(r'$\rm{mag}$',fontsize=ftsize)    
    ax2 = plt.axes(axis2)
    ax2.set_xlim( ax1.get_xlim() )
    #ax2.set_ylim( ax1.get_ylim() )
    ax2.xaxis.set_major_formatter(nullfmt)
    
    return ax1,ax2

def two_panel_plot_vert_bsquat(oney=False,ftsize=20,label=None):
    fig = plt.figure(2, figsize=(8,8))
    left, width = 0.13, 0.83
    bottom, height = 0.1, 0.21
    height2 = .60
    dh = 0.03
    
    axis1 = [left, bottom, width, height]
    axis2 = [left, (bottom+height+dh), width, height2]
    
    ax1 = plt.axes(axis1)
    if oney==True: 
        ax1.annotate(r'$%s$'%label,fontsize=ftsize,xy=(0.025,.5),xycoords='figure fraction',va='center',rotation='vertical')
    ax2 = plt.axes(axis2)
    ax2.xaxis.set_major_formatter(nullfmt)
    
    return ax1,ax2

def two_panel_plot_vert2(oney=True,mag2_cut=30,mag2_max=16,ftsize=20):
    fig = plt.figure(2, figsize=(8,8))
    left, width = 0.15, 0.81
    bottom, height = 0.1, 0.41
    dh = 0.03
    
    axis1 = [left, bottom, width, height]
    axis2 = [left, (bottom+height+dh), width, height]
    
    ax1 = plt.axes(axis1)
    ax1.set_xlim( (mag2_cut,mag2_max) ) # set all axes limits here
    #ax1.set_ylim( (0.0001, 10.) )
    if oney==True: ax1.annotate(r'$\#/ 3.6 \mu \rm{m\ Region\ Integrated\ Flux\ (Jy}^{-1}\rm{)}$',fontsize=ftsize,xy=(0.025,.5),xycoords='figure fraction',va='center',rotation='vertical')
    ax1.set_xlabel(r'$\rm{mag}$',fontsize=ftsize)    
    ax2 = plt.axes(axis2)
    ax2.set_xlim( ax1.get_xlim() )
    #ax2.set_ylim( ax1.get_ylim() )
    ax2.xaxis.set_major_formatter(nullfmt)
    
    return ax1,ax2

# For cmds_and_tracks:
def setup_four_panel(ftsize=20):
    fig = plt.figure(figsize=(8,8))
    left, width = 0.1, 0.4
    bottom, height = 0.1, 0.4
    lefter = left+width+0.01
    higher = bottom+height+0.01
    
    # plot and fig sizes
    fig = plt.figure(1, figsize=(8,8))
    
    ll_axis = [left, bottom, width, height]
    lr_axis = [lefter, bottom, width, height] 
    ul_axis = [left, higher, width, height]
    ur_axis = [lefter, higher, width, height]
    
    ax_ll = plt.axes(ll_axis)
    ax_ll.set_xlim( (-0.75,1.5)) # model and data x limits here
    ax_ll.set_ylim( (25, 18)) # set all y limits here
    
    ax_lr = plt.axes(lr_axis)
    ax_lr.set_xlim( ax_ll.get_xlim() )
    ax_lr.set_ylim( ax_ll.get_ylim() )    
    
    ax_ul = plt.axes(ul_axis)
    ax_ul.set_xlim( ax_ll.get_xlim() )
    ax_ul.set_ylim( ax_ll.get_ylim() )    
    
    ax_ur = plt.axes(ur_axis)
    ax_ur.set_xlim( ax_ll.get_xlim() )
    ax_ur.set_ylim( ax_ll.get_ylim() )    
    
    ax_lr.yaxis.set_major_formatter(nullfmt)
    ax_ur.yaxis.set_major_formatter(nullfmt)
    ax_ur.xaxis.set_major_formatter(nullfmt)
    ax_ul.xaxis.set_major_formatter(nullfmt)
    
    # titles
    #x = fig.text(0.5,0.96,r'$\rm{%s}$' % ('Disk Field'),horizontalalignment='center',verticalalignment='top',size=20)
    ax_ur.set_title(r'$\rm{Disk\ Field\ (\sim 6.5\ kpc)}$',color='black',fontsize=ftsize)
    ax_ul.set_title(r'$\rm{Bulge}\ Field\ (< 120\ \rm{pc})$',color='black',fontsize=ftsize)
    ax_ll.set_ylabel(r'$F336W$',fontsize=ftsize)
    ax_ll.set_xlabel(r'$F275W-F336W$',fontsize=ftsize)
    ax_ul.set_ylabel(ax_ll.get_ylabel(),fontsize=ftsize)    
    ax_lr.set_xlabel(ax_ll.get_xlabel(),fontsize=ftsize)
    
    return ax_ll,ax_lr,ax_ul,ax_ur

def age_ticks(age,ax,color,mag,age_spacing=1e4):
    tick = []
    col_lim = ax.get_xlim()
    mag_lim = ax.get_ylim()
    ind = np.nonzero(
        (color < col_lim[1]) &
        (color > col_lim[0]) &
        (mag < mag_lim[0]) &
        (mag > mag_lim[1])
        )[0]
    if len(ind)<2:
        'no dice.'
    else:
        dt = np.diff(age[ind])
        start = 0
        for i in range(len(dt)):
            so_far = np.sum(dt[start:i])
            if so_far > age_spacing:
                tick.append(i)
                start = i
    return color[ind][tick],mag[ind][tick]
    #age_marks = np.arange(low_age,high_age,age_spacing)
    #tick = []
    #for age in age_marks:
    #    tick.append(closest_match(age,arr))
    #return tick
    
def plot_tracks_on_four_panel_smooth(tracks,ax1,ax2,Av,af275w,af336w,**kwargs):
    annotate = kwargs.get('annotate',True)
    dmod = 24.47
    Filter1 = 'F275W'
    Filter2 = 'F336W'
    from PadovaTracksUtils import get_Hfus_HeB
    for i,track in enumerate(tracks):
        if track == 0: continue
        tm1,tm2,age,HeB,Hfus = get_Hfus_HeB(track,Filter1,Filter2)
        #Mag 2 mag:
        tm1= tm1+dmod+Av*af275w
        tm2= tm2+dmod+Av*af336w
        t_col = tm1-tm2
        #tmag = tm1
        tmag = tm2
        mass = info_from_track_filename(track.name)['M']
        if float(mass) < 1.:
            try:
                ax1.plot(smooth(t_col[HeB]),smooth(tmag[HeB]),
                         color='white',zorder=3,lw=2)
                ax1.plot(smooth(t_col[HeB]),smooth(tmag[HeB]),
                         color='red',zorder=3,lw=1.5,label='HeB')
                tmp_col,tmp_mag = age_ticks(age[HeB],ax1,t_col[HeB],tmag[HeB])
                if float(mass) > 0.5:
                    if annotate == True:
                      ax1.errorbar(tmp_col,tmp_mag,yerr=0.1,
                                     color='red',elinewidth=None,zorder=4,
                                     linewidth=1.5,capsize=0,linestyle='None')
                if (len(HeB) != 0) & (mass != 0.480): 
                    if annotate == True:
                        ax1.annotate('%.2f'%mass,xy=(t_col[np.argmin(tmag[HeB])],
                                                     min(tmag[HeB])),
                                                     color='red',**ann_kwargs)
            except ValueError:
                pass
        
        if (int(mass)/mass == 1.) & (mass < 12):
            try:
                ax2.plot(smooth(t_col[Hfus]),smooth(tmag[Hfus]),
                         color='white',zorder=3,lw=2)
                ax2.plot(smooth(t_col[Hfus]),smooth(tmag[Hfus]),
                         color='blue',zorder=3,lw=1.5)
                tmp_col,tmp_mag = age_ticks(age[Hfus],ax2,t_col[Hfus],
                                            tmag[Hfus])
                if len(Hfus) != 0:
                    if annotate == True:
                        ax2.annotate('%.2f'%mass,xy=(min(t_col[Hfus]),
                                                     tmag[np.argmin(t_col[Hfus])]),
                                                     color='blue',ha='right',**ann_kwargs)
            except ValueError:
                pass

   
    for file in os.listdir(vw_loc):
        if file.endswith('.phat'):
            mass = file.split('z')[0].replace('m','')
            hm = tablefile2dict(vw_loc+file,'#')
            try:
                hm_mag1 =np.array(hm[Filter1])
                hm_mag2 =np.array(hm[Filter2])
            except KeyError:
                hm_mag1 =np.array(hm[Filter1+'1'])
                hm_mag2 =np.array(hm[Filter2+'1'])
            #Mag 2 mag:
            hm_mag1= hm_mag1+dmod+Av*af275w
            hm_mag2= hm_mag2+dmod+Av*af336w
            hm_col = np.array(hm_mag1[1:])-np.array(hm_mag2[1:])
            p = ax1.plot(smooth(hm_col),smooth(hm_mag2[1:]),
                         color = 'white',lw=2)
            tmp_col,tmp_mag = age_ticks(np.array(hm['age'][1:]),
                                        ax1,hm_col,hm_mag2[1:])
            ax1.plot(tmp_col,tmp_mag,'o',color='purple',alpha=0.3,zorder=4)
            p = ax1.plot(smooth(hm_col),smooth(hm_mag2[1:]),
                         color = 'purple',lw=1.5,label='P-AGB')
            if annotate == True:
                ax1.annotate('%.2f'%float(mass),
                             xy=(0.0,min(smooth(hm_mag2))),
                             color='purple',ha='left',**ann_kwargs)
            
    return

def plot_tracks_on_one_panel_smooth(tracks,ax1,Av,**kwargs):
    annotate = kwargs.get('annotate',True)
    dmod = 24.47
    Filter1 = 'F275W'
    Filter2 = 'F336W'
    yfilter = kwargs.get('yfilter')
    af275w = kwargs.get('af275w')
    af336w = kwargs.get('af336w')
    from PadovaTracksUtils import get_Hfus_HeB
    for i,track in enumerate(tracks):
        if track == 0: continue
        tm1,tm2,age,HeB,Hfus = get_Hfus_HeB(track,Filter1,Filter2)
        #Mag 2 mag:
        tm1= tm1+dmod+Av*af275w
        tm2= tm2+dmod+Av*af336w
        t_col = tm1-tm2
        if yfilter == Filter2:
            tmag = tm2
        else:
            tmag = tm1
        mass = info_from_track_filename(track.name)['M']
        if float(mass) < 1.:
            try:
                lw = 1.5
                if float(mass) != 0.58: continue
                lw = 4
                ax1.plot(smooth(t_col[HeB]),smooth(tmag[HeB]),
                         color='white',zorder=3,lw=2)
                ax1.plot(smooth(t_col[HeB]),smooth(tmag[HeB]),
                         color='red',zorder=3,lw=lw,label='HeB')
                tmp_col,tmp_mag = age_ticks(age[HeB],ax1,t_col[HeB],tmag[HeB])
                if float(mass) > 0.5:
                    if annotate == True:
                      ax1.errorbar(tmp_col,tmp_mag,yerr=0.1,
                                     color='red',elinewidth=None,zorder=4,
                                     linewidth=1.5,capsize=0,linestyle='None')
                if (len(HeB) != 0) & (mass != 0.480): 
                    if annotate == True:
                        ax1.annotate('%.2f'%mass,xy=(t_col[np.argmin(tmag[HeB])],
                                                     min(tmag[HeB])),
                                                     color='red',**ann_kwargs)
            except ValueError:
                pass
    
    for file in os.listdir(vw_loc):
        continue    
        if file.endswith('.phat'):
            mass = file.split('z')[0].replace('m','')
            hm = tablefile2dict(vw_loc+file,'#')
            try:
                hm_mag1 =np.array(hm[Filter1])
                hm_mag2 =np.array(hm[Filter2])
            except KeyError:
                hm_mag1 =np.array(hm[Filter1+'1'])
                hm_mag2 =np.array(hm[Filter2+'1'])
            #Mag 2 mag:
            hm_mag1= hm_mag1+dmod+Av*af275w
            hm_mag2= hm_mag2+dmod+Av*af336w
            if yfilter == Filter2:
                hm_mag = hm_mag2
            else:
                hm_mag = hm_mag1
            hm_col = np.array(hm_mag1[1:])-np.array(hm_mag2[1:])
            p = ax1.plot(smooth(hm_col),smooth(hm_mag[1:]),
                         color = 'white',lw=2)
            tmp_col,tmp_mag = age_ticks(np.array(hm['age'][1:]),
                                        ax1,hm_col,hm_mag[1:])
            ax1.plot(tmp_col,tmp_mag,'o',color='purple',alpha=0.3,zorder=4)
            p = ax1.plot(smooth(hm_col),smooth(hm_mag[1:]),
                         color = 'purple',lw=1.5,label='P-AGB')
            if annotate == True:
                ax1.annotate('%.2f'%float(mass),
                             xy=(0.0,min(smooth(hm_mag))),
                             color='purple',ha='left',**ann_kwargs)
            
    return

def plot_tracks_on_four_panel(tracks,ax1,ax2,Av,**kwargs):
    dmod = kwargs.get('dmod')
    af275w = kwargs.get('af275w')
    af336w = kwargs.get('af336w')
    Filter1 = kwargs.get('Filter1')
    Filter2 = kwargs.get('Filter2')
    A1 = Av
    A2 = Av
    from PadovaTracksUtils import get_Hfus_HeB
    for i,track in enumerate(tracks):
        if track == 0: continue
        tm1,tm2,age,HeB,Hfus = get_Hfus_HeB(track,Filter1,Filter2)
        #Mag 2 mag:
        tm1= tm1+dmod+A1*af275w
        tm2= tm2+dmod+A2*af336w
        t_col = tm1-tm2
        #tmag = tm1
        tmag = tm2
        if float(info_from_track_filename(track.name)['M']) < 1.:
            try:
                ax1.plot(t_col[HeB],tmag[HeB],color='darkred',zorder=3,lw=2)#,alpha =1.-info_from_track_filename(phat)['M']/10.)
                ax1.plot(t_col[HeB],tmag[HeB],color='red',zorder=3,lw=1.5,label='HeB')
                #ax1.plot(t_col,tmag,color='darkgreen',zorder=3,lw=1)
                if len(HeB) != 0: 
                    ax1.annotate('%.2f'%info_from_track_filename(track.name)['M'],xy=(t_col[np.argmin(tmag[HeB])],min(tmag[HeB])),color='red',**ann_kwargs)
            except ValueError:
                pass
        
        if (int(info_from_track_filename(track.name)['M'])/info_from_track_filename(track.name)['M'] == 1.) & (info_from_track_filename(track.name)['M'] < 12):
            try:
                #ax2.plot(t_col[Hfus][:-4],tmag[Hfus][:-4],color='pink',zorder=3,lw=2)
                #ax2.plot(t_col[Hfus][4:-1],tmag[Hfus][4:-1],color='blue',zorder=3,lw=1.5)
                #ax2.plot(t_col[Hfus][2],tmag[Hfus][2],'o',color='green',zorder=3,lw=1.5)
                ax2.plot(t_col[Hfus],tmag[Hfus],color='darkblue',zorder=3,lw=2)
                ax2.plot(t_col[Hfus],tmag[Hfus],color='blue',zorder=3,lw=1.5)
                #ax2.plot(t_col,tmag,color='darkgreen',zorder=3,lw=1)
                if len(Hfus) != 0: ax2.annotate('%.2f'%info_from_track_filename(track.name)['M'],xy=(min(t_col[Hfus]),tmag[np.argmin(t_col[Hfus])]),color='blue',ha='right',**ann_kwargs)
            except ValueError:
                pass
                
    for file in os.listdir(vw_loc):
        if file.endswith('.phat'):
            mass = file.split('z')[0].replace('m','')
            hm = tablefile2dict(vw_loc+file,'#')
            try:
                hm_mag1 =np.array(hm[Filter1])
                hm_mag2 =np.array(hm[Filter2])
            except KeyError:
                hm_mag1 =np.array(hm[Filter1+'1'])
                hm_mag2 =np.array(hm[Filter2+'1'])
            #Mag 2 mag:
            hm_mag1= hm_mag1+dmod+A1*af275w
            hm_mag2= hm_mag2+dmod+A2*af336w
            hm_col = np.array(hm_mag1[1:])-np.array(hm_mag2[1:])
            p = ax1.plot(hm_col,hm_mag2[1:],color = 'white',lw=2)
            p = ax1.plot(hm_col,hm_mag2[1:],color = 'purple',lw=1.5,label='P-AGB')
            ax1.annotate('%.2f'%float(mass),xy=(0.0,min(hm_mag2)),color='purple',ha='left',**ann_kwargs)
            
    return
    
def setup_five_panel_plot(ftsize=20,fig=1,noy=True):
    nullfmt   = NullFormatter() # no labels
    left, width = 0.06, 0.175
    bottom, height = 0.14, 0.82
    dl = 0.01
    lefts = [left+(width+dl)*float(i) for i in range(5)]

    # plot and fig sizes
    fig = plt.figure(fig, figsize=(15,4))

    axiss = [[l,bottom,width,height] for l in lefts]

    axs = []
    for i in range(len(axiss)):
        axs.append( plt.axes(axiss[i]))

    for ax in axs:
        ax.set_xlim( (-0.75,2)) # model and data x limits here
        ax.set_ylim( (24.3, 19)) # set all y limits here
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        if axs.index(ax) > 0 and noy==True: ax.yaxis.set_major_formatter(nullfmt)
        ax.set_xlabel(r'$F275W-F336W$',fontsize=ftsize)
        if axs.index(ax) ==0: ax.set_ylabel(r'$F336W$',fontsize=ftsize)
        
    return axs

def setup_ten_panel_plot(ftsize=20,fig=1,noy=True):
    nullfmt   = NullFormatter() # no labels
    # main plot
    left, width = 0.06, 0.175
    bottom, height = 0.43, 0.54
    dl = 0.01
    lefts = [left+(width+dl)*float(i) for i in range(5)]
    
    axiss = [[l,bottom,width,height] for l in lefts]
    
    # bottom plot
    bottom, height = 0.14, 0.25
    baxiss = [[l,bottom,width,height] for l in lefts]
    
    # plot and fig sizes
    fig = plt.figure(fig, figsize=(15,6))

    baxs = []
    axs = []
    for i in range(len(axiss)):
        axs.append( plt.axes(axiss[i]))
        baxs.append( plt.axes(baxiss[i]))
        
    for ax in axs:
        if axs.index(ax) > 0 and noy==True: ax.yaxis.set_major_formatter(nullfmt)
        ax.xaxis.set_major_formatter(nullfmt)
        if axs.index(ax) ==0: ax.set_ylabel(r'$F336W$',fontsize=ftsize)
    for ax in baxs:
        ax.set_xlabel(r'$F275W-F336W$',fontsize=ftsize)
        
    return axs,baxs

def two_panel_squeeze(**kwargs):
    xlab1 = kwargs.get('xlab1','')
    xlab2 = kwargs.get('xlab2','')
    ylab = kwargs.get('ylab','')
    sizex = kwargs.get('sizex',10)
    sizey = kwargs.get('sizey',5)
    ftsize = kwargs.get('ftsize',20)
    fig = kwargs.get('fig',2)
    fig = plt.figure(fig, figsize=(sizex,sizey))
    left, width = 0.1, 0.4
    bottom, height = 0.12, 0.8  
    
    left2 = left+width
    axis1 = [left, bottom, width, height]
    axis2 = [left2, bottom, width, height]
    
    ax1 = plt.axes(axis1)
    ax1.set_xlabel(r'%s'%xlab1,fontsize=ftsize)
    ax1.set_ylabel(r'%s'%ylab,fontsize=ftsize)
    
    ax2 = plt.axes(axis2)
    ax2.set_xlim( ax1.get_xlim() )
    ax2.set_ylim( ax1.get_ylim() )
    ax2.set_xlabel(r'%s'%xlab2,fontsize=ftsize)
    return ax1,ax2
    
def quick_plot(**kwargs):
    sizex = kwargs.get('sizex',8)
    sizey = kwargs.get('sizey',8)
    ftsize = kwargs.get('ftsize',20)
    xlab = kwargs.get('xlab','')
    ylab = kwargs.get('ylab','')
    fig = kwargs.get('fig',2)
    fig = plt.figure(fig, figsize=(sizex,sizey))
    ax = plt.axes()
    ax.set_xlabel(r'%s'%xlab,fontsize=ftsize)
    ax.set_ylabel(r'%s'%ylab,fontsize=ftsize)
    return ax
    
def make_single_plot(arr1,arr2,overpltinds,rad_inds,filter1,filter2,name,xlim,ylim,cols,fig=2):
    fmt = '$%s$ %s'
    kwargs = {'xlab': fmt%(filter2,name), 
              'ylab': fmt%(filter1,name),
              'fig' : fig}
    
    ax = quick_plot(**kwargs)
    ax.plot(arr2,arr1,'.',color='grey')
    for i in range(len(rad_inds)):
        inds = list(set(rad_inds[4-i]) & set(overpltinds))
        ax.plot(arr2[inds],arr1[inds],'.',color=cols[4-i])
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    return ax

def make_double_plot(x1,y1,x2,y2,overpltinds1,overpltinds2,rad_inds,filter1,filter2,ylab,xlim,ylim,cols,fig=2):
    fmt = '$%s$'
    pltkwargs = {'xlab1': fmt%filter1,
                 'xlab2': fmt%filter2,
                 'ylab': ylab,
                 'fig': fig}
    ax1,ax2 = two_panel_squeeze(**pltkwargs)
    ax1.plot(x1,y1,'.',color='grey')
    for i in range(len(rad_inds)):
        inds = list(set(rad_inds[4-i]) & set(overpltinds1))
        ax1.plot(x1[inds],y1[inds],'.',color=cols[4-i])

    ax2.plot(x2,y2,'.',color='grey')
    for i in range(len(rad_inds)):
        inds = list(set(rad_inds[4-i]) & set(overpltinds2))
        ax2.plot(x2[inds],y2[inds],'.',color=cols[4-i])
    ax1.set_xlim(xlim)
    ax1.set_ylim(ylim)
    ax2.set_xlim(xlim)
    ax2.set_ylim(ylim)
    return ax1,ax2

def five_cmds(phot, rad_inds, cut_inds,**kwargs):
    cols = kwargs.get('cols')
    distances = kwargs.get('distances',[34,89,182,353,640])
    alpha = kwargs.get('alpha',1.)
    axs = setup_five_panel_plot()
    for i,ind in enumerate(rad_inds):
        ax = axs[i]
        inds = list(set(cut_inds) & set(ind))
        color = phot['mag1'][inds]-phot['mag2'][inds]
        mag = phot['mag2'][inds]
        #ax.plot(color[ind],mag[ind],'.',ms=6,zorder=2,color='gray')
        ax.plot(color,mag,'.',ms=6,zorder=2,color=cols[i],alpha=alpha)
        ax.set_ylim(25,19)
        ax.annotate(r'%0.f pc'%distances[i],xy=(-.7,19.5),ha='left',fontsize=ftsize-3)
    plt.savefig('five_cmds%s'%EXT)
    plt.close()
    return axs

def five_cmds_complines(phot,rad_inds,both,**kwargs):
    cols = kwargs.get('cols')
    distances = kwargs.get('distances',[34,89,182,353,640])
    alpha = kwargs.get('alpha',1.)
    comp_data = kwargs.get('comp_data')
    
    comp_dat = file_io.load_comp_data(comp_data)
    axs = setup_five_panel_plot()
    for i,ind in enumerate(rad_inds):
        ax = axs[i]
        inds = list(set(ind) & set(both))
        color = phot['mag1'][inds]-phot['mag2'][inds]
        mag = phot['mag2'][inds]
        #ax.plot(phot['mag1'][ind]-phot['mag2'][ind],phot['mag2'][ind],'.',ms=3,zorder=2,color='grey',alpha=0.1)
        ax.plot(color,mag,'.',ms=6,zorder=2,color=cols[i],mec=cols[i],alpha=alpha)
        ax.set_ylim(26,19)
        ax.set_xlim(-.8,2)
        ax.annotate(r'%0.f pc'%distances[i],xy=(-.75,19.5),ha='left',fontsize=ftsize-3)
        col50f2,mag50f2 = put_a_line_on_it(ax,comp_dat['comp50mag2'][i],consty=True)
        col50f1,mag50f1 = put_a_line_on_it(ax,comp_dat['comp50mag1'][i],consty=False,constx=True,filter1=True)
        col90f2,mag90f2 = put_a_line_on_it(ax,comp_dat['comp90mag2'][i],consty=True,color='white',ls='-',lw=2)
        col90f1,mag90f1 = put_a_line_on_it(ax,comp_dat['comp90mag1'][i],consty=False,constx=True,filter1=True,color='white',ls='-',lw=2)
        col50,mag50 = combine_two_lines(col50f2,mag50f2,col50f1,mag50f1)
        col90,mag90 = combine_two_lines(col90f2,mag90f2,col90f1,mag90f1)
        ax.plot(col50,mag50,color='white',lw=2)
        ax.plot(col90,mag90,color='white',lw=2)
        ax.plot(col50,mag50,color='black')
        ax.plot(col90,mag90,color='black',ls='--')
        
    plt.savefig('five_cmds_lines%s'%EXT)
    plt.close()
    return axs

def put_feh_on_top(xdata,ax,twoplot=False,offset_text=0.1,offset_tick=-0.17,shk=0.3):
    #[Fe/H]
    feh = [0.2,0.15,0.06,0.05,0.01]
    if twoplot==True:
        coords = 'axes fraction'
        text = (0.4,1.05)
    else:
        coords = 'figure fraction'
        text=(0.5,0.95)
    for i in range(len(feh)):
        ax.annotate(r'$%.2f$'%feh[i],
                    xytext=(xdata[i],ax.get_ylim()[1]+offset_text),
                    xy=(xdata[i],ax.get_ylim()[1]+offset_tick),
                    arrowprops=dict(facecolor='grey', width=1,shrink=shk,lw=0.0,headwidth=0),
                    fontsize=ftsize-2,
                    ha='center')
    for tick in np.concatenate((ax.xaxis.get_major_ticks(),ax.xaxis.get_minor_ticks())):
        tick.tick2On = False
    ax.annotate(r'$[{\rm Z/H}]$',fontsize=ftsize,xy=text,xycoords=coords)
    return

def HRD_Figure(**kwargs):
    '''
    This is for the HRD figure in the intro.
    '''
    dmod = kwargs.get('dmod')   
    Av = kwargs.get('Av_bulge')
    af275w = kwargs.get('af275w')
    af336w = kwargs.get('af336w')    
    Filter1 = kwargs.get('Filter1')
    Filter2 = kwargs.get('Filter2')
    yfilter = kwargs.get('yfilter',Filter2)
    btrack = kwargs.get('btrack')
    btrackmore = kwargs.get('btrackmore')
    
    xlab1 = r'$\log \ T_{{\rm eff}} \ (K)$'
    ylab = '$\log \ L \ (L_\odot$)'
    xlab2 = '$%s-%s$'%(Filter1,Filter2)
    ylab2 = '$%s$'%yfilter
    
    ax1, ax2 = two_panel_plot(10,5,xlab1,xlab2,ylab,ylab2)
    lmin =1.
    ax1.set_ylim( lmin,4.1 )
    ax1.set_xlim( 5.2,3.3 )
    ax2.set_xlim( (-1.,3.0) )
    ax2.set_ylim( (27., 18.5) )
    thick = 2
    pagb = tablefile2dict(vw_loc+'m0.597z016H.dat.phat','#')
    peagb = get_track_by_mass(btrack,0.58)
    agbm = get_track_by_mass(btrackmore,.48)
    agb = get_track_by_mass(btrack,1.6)
    i=0
    for track in [agbm,peagb,agb]:
        tm1 = Tracks.get_stable_col(track,Filter1+'1')+dmod+Av*af275w
        tm2 = Tracks.get_stable_col(track,Filter2+'1')+dmod+Av*af336w
        t_col = tm1-tm2
        if yfilter == Filter2:
            tmag = tm2
        else:
            tmag = tm1
        logl = Tracks.get_stable_col(track,'LOG_L')
        logte = Tracks.get_stable_col(track,'LOG_TE')
        if track == agbm: 
            ax1.plot(smooth(logte[4:-1]),smooth(logl[4:-1]),color='blue',zorder=3,lw=thick)
            ax2.plot(smooth(t_col[4:-1]),smooth(tmag[4:-1]),color='blue',zorder=3,lw=thick)
            for arrow_pt in [12,25,40]:
                agbm_arr = arrow_on_line(smooth(logte),smooth(logl),arrow_pt,'blue')
                ax1.add_patch(agbm_arr)
            for arrow_pt in [14,20,40]:
                agbm_arr = arrow_on_line(smooth(t_col),smooth(tmag),arrow_pt,'blue',width=.2)
                ax2.add_patch(agbm_arr)
        if track == peagb:
            logl[4:-1][500]
            ax1.plot(smooth(logte[4:-1]),smooth(logl[4:-1]),color='red',zorder=3,lw=thick)
            ax2.plot(smooth(t_col[4:-1]),smooth(tmag[4:-1]),color='red',zorder=3,lw=thick)
            for arrow_pt in [750,850,1200,1400]:
                peagb_arr = arrow_on_line(smooth(logte[4:-1]),smooth(logl[4:-1]),
                                          arrow_pt,'red')
                ax1.add_patch(peagb_arr)
            for arrow_pt in [850,1000,1320]:
                peagb_arr = arrow_on_line(smooth(t_col[4:-1]),smooth(tmag[4:-1]),
                                          arrow_pt,'red',width=.2)
                ax2.add_patch(peagb_arr)
        if track == agb:
            ax1.plot(smooth(logte),smooth(logl),color='black',zorder=3,lw=thick)
            ax2.plot(smooth(t_col),smooth(tmag),color='black',zorder=3,lw=thick)
            for arrow_pt in [700,4000]:
                agb_arr = arrow_on_line(smooth(logte),smooth(logl),arrow_pt,'black')
                ax1.add_patch(agb_arr)
        plt.savefig('hrd%i.png'%i)
        i+=1
    # ZAHB:
    logls,logtes = np.array([]),np.array([])
    t_cols,tmags = np.array([]),np.array([])
    for track in btrackmore:
        tm1 = Tracks.get_stable_col(track,Filter1+'1')+dmod+Av*af275w
        tm2 = Tracks.get_stable_col(track,Filter2+'1')+dmod+Av*af336w
        t_col = tm1-tm2
        if yfilter == Filter2:
            tmag = tm2
        else:
            tmag = tm1
        logl = Tracks.get_stable_col(track,'LOG_L')
        logte = Tracks.get_stable_col(track,'LOG_TE')
        logls = np.append(logls,logl[1:4])
        logtes =np.append(logtes,logte[1:4])
        t_cols = np.append(t_cols,t_col[1:4])
        tmags =np.append(tmags,tmag[1:4])
        
    ax1.plot(smooth(np.sort(logtes)),smooth(logls[np.argsort(logtes)]),
             color='darkgreen',zorder=3,lw=thick)
    ax2.plot(smooth(np.sort(t_cols)),smooth(tmags[np.argsort(t_cols)]),
             color='darkgreen',zorder=3,lw=thick)
    
    # VW Tracks
    hm = pagb
    hm_logl  = np.array(hm['logl'])
    brights = np.nonzero(hm_logl > lmin)
    hm_logte = np.array(hm['logte'])
    try:
        hm_mag1 =np.array(hm[Filter1])
        hm_mag2 =np.array(hm[Filter2])
    except KeyError:
        hm_mag1 =np.array(hm[Filter1+'1'])
        hm_mag2 =np.array(hm[Filter2+'1'])
    hm_mag1= hm_mag1+dmod+Av*af275w
    hm_mag2= hm_mag2+dmod+Av*af336w
    hm_col = np.array(hm_mag1[1:])-np.array(hm_mag2[1:])
    if yfilter == Filter2:
        hm_mag = hm_mag2
    else:
        hm_mag = hm_mag1
    p = ax2.plot(smooth(hm_col[brights]),smooth(hm_mag[1:][brights]),color ='purple',lw=thick)
    p = ax1.plot(smooth(hm_logte[brights]),smooth(hm_logl[brights]),color ='purple',lw=thick)
    for arrow_pt in [10,20]:
        pagb_arr = arrow_on_line(smooth(hm_logte[brights]),smooth(hm_logl[brights]),
                                 arrow_pt,'purple')
        ax1.add_patch(pagb_arr)
    for arrow_pt in [10,20]:
        pagb_arr = arrow_on_line(smooth(hm_col[brights]),smooth(hm_mag[1:][brights]),
                                 arrow_pt,'purple',width=0.2)
        ax2.add_patch(pagb_arr)
    
    fmt = r'$\mathsf{%s}$'
    ftsize = 14
    ax1.annotate(fmt%'AGB',xy=(3.7,3.4),color='black',ha='left',fontsize=ftsize,**ann_kwargs)
    ax1.annotate(fmt%'PE-AGB',xy=(4.5,3.17),color='red',ha='left',fontsize=ftsize,**ann_kwargs)
    ax1.annotate(fmt%'AGB-Manqu\grave e',xy=(4.9,2.2),color='blue',ha='left',fontsize=ftsize,**ann_kwargs)
    ax1.annotate(fmt%'P-AGB',xy=(5,3.8),color='purple',ha='left',fontsize=ftsize,**ann_kwargs)
    ax1.annotate(fmt%'ZAHB',xy=(4,1.3),color='darkgreen',ha='left',fontsize=ftsize,**ann_kwargs)
    
    withdata = kwargs.get('withdata',False)
    if withdata == True:
        mag1 = kwargs['data_mag1']
        mag2 = kwargs['data_mag2']
        if yfilter == Filter2:
            mag = mag2
        else:
            mag = mag1
        ax2.plot(mag1-mag2,mag,'.',ms=1.5,mec=None,color='black')
    
    ax2.annotate(fmt%'AGB',xy=(1.4,29.6),color='black',ha='left',fontsize=ftsize,**ann_kwargs)
    if yfilter == Filter2:
        ax2.annotate(fmt%'PE-AGB',xy=(.88,22.35),color='red',ha='left',fontsize=ftsize,**ann_kwargs)
    else:
        ax2.annotate(fmt%'PE-AGB',xy=(-0.3,21.5),color='red',ha='left',fontsize=ftsize,**ann_kwargs)
    ax2.annotate(fmt%'AGB-Manqu\grave e',xy=(-.24,24.15),color='blue',ha='left',fontsize=ftsize,**ann_kwargs)
    ax2.annotate(fmt%'P-AGB',xy=(-0.3,20.2),color='purple',ha='left',fontsize=ftsize,**ann_kwargs)
    ax2.annotate(fmt%'ZAHB',xy=(-0.07,26.6),color='darkgreen',ha='left',fontsize=ftsize,**ann_kwargs)
    
    xs = [-.8,2,2,-.8,-.8]
    ys = [25,25,19,19,25]
    ax2.plot(xs,ys,'--',lw=3,color='grey')
    plt.savefig('hrd%s'%EXT)
    print 'wrote hrd%s'%EXT
    plt.close()
    return

def HRD_Figure_one_by_one(**kwargs):
    '''
    This is for the HRD figure in the intro but saved at each track draw.
    '''
    fmt = r'$\mathsf{%s}$'
    ftsize = 14
    dmod = kwargs.get('dmod')   
    if dmod == 0:
        offset = 24.47
    else:
        offset = 0.
    Av = kwargs.get('Av_bulge')
    af275w = kwargs.get('af275w')
    af336w = kwargs.get('af336w')    
    Filter1 = kwargs.get('Filter1')
    Filter2 = kwargs.get('Filter2')
    yfilter = kwargs.get('yfilter',Filter2)
    btrack = kwargs.get('btrack')
    btrackmore = kwargs.get('btrackmore')
    if Filter1 == 'F275W':
        file_ext = '.phat'
    else:
        file_ext = '.galex_sloan'
    
    xlab1 = r'$\log \ T_{{\rm eff}} \ (K)$'
    ylab = '$\log \ L \ (L_\odot$)'
    xlab2 = '$%s-%s$'%(Filter1,Filter2)
    ylab2 = '$%s$'%yfilter
    
    ax1, ax2 = two_panel_plot(10,5,xlab1,xlab2,ylab,ylab2)
    lmin =1.
    ax1.set_ylim( lmin,4.1 )
    ax1.set_xlim( 5.2,3.3 )
    ax2.set_xlim( (-1.,3.0) )
    ax2.set_ylim( (27.-offset, 18.5-offset) )
    thick = 2
    pagb = tablefile2dict(vw_loc+'m0.597z016H.dat'+file_ext,'#')
    peagb = get_track_by_mass(btrack,0.58)
    agbm = get_track_by_mass(btrackmore,.48)
    agb = get_track_by_mass(btrack,1.6)
    i=0
    # ZAHB:
    logls,logtes = np.array([]),np.array([])
    t_cols,tmags = np.array([]),np.array([])
    for track in btrackmore:
        if track==0: continue
        try:
            tm1 = Tracks.get_stable_col(track,Filter1+'1')+dmod+Av*af275w
            tm2 = Tracks.get_stable_col(track,Filter2+'1')+dmod+Av*af336w
        except ValueError:
            tm1 = Tracks.get_stable_col(track,Filter1)+dmod+Av*af275w
            tm2 = Tracks.get_stable_col(track,Filter2)+dmod+Av*af336w
        t_col = tm1-tm2
        if yfilter == Filter2:
            tmag = tm2
        else:
            tmag = tm1
        logl = Tracks.get_stable_col(track,'LOG_L')
        logte = Tracks.get_stable_col(track,'LOG_TE')
        logls = np.append(logls,logl[1:4])
        logtes =np.append(logtes,logte[1:4])
        t_cols = np.append(t_cols,t_col[1:4])
        tmags =np.append(tmags,tmag[1:4])
        plt.savefig('hrd%i.png'%i)
    ax1.plot(smooth(np.sort(logtes)),smooth(logls[np.argsort(logtes)]),
             color='darkgreen',zorder=3,lw=thick)
    ax2.plot(smooth(np.sort(t_cols)),smooth(tmags[np.argsort(t_cols)]),
             color='darkgreen',zorder=3,lw=thick)
    ax1.annotate(fmt%'ZAHB',xy=(4,1.3),color='darkgreen',ha='left',fontsize=ftsize,**ann_kwargs)
    plt.savefig('hrd%i.png'%i)
    i+=1
    for track in [agb,peagb,agbm]:
        print track
        try:
            tm1 = Tracks.get_stable_col(track,Filter1+'1')+dmod+Av*af275w
            tm2 = Tracks.get_stable_col(track,Filter2+'1')+dmod+Av*af336w
        except ValueError:
            tm1 = Tracks.get_stable_col(track,Filter1)+dmod+Av*af275w
            tm2 = Tracks.get_stable_col(track,Filter2)+dmod+Av*af336w
        t_col = tm1-tm2
        if yfilter == Filter2:
            tmag = tm2
        else:
            tmag = tm1
        logl = Tracks.get_stable_col(track,'LOG_L')
        logte = Tracks.get_stable_col(track,'LOG_TE')
        if track == agb:
            ax1.plot(smooth(logte),smooth(logl),color='black',zorder=3,lw=thick)
            ax2.plot(smooth(t_col),smooth(tmag),color='black',zorder=3,lw=thick)
            if 'phat' in file_ext:
                for arrow_pt in [700,4000]:
                    agb_arr = arrow_on_line(smooth(logte),smooth(logl),arrow_pt,'black')
                    ax1.add_patch(agb_arr)
            ax1.annotate(fmt%'AGB',xy=(3.7,3.4),color='black',ha='left',fontsize=ftsize,**ann_kwargs)
            # VW Tracks
            hm = pagb
            hm_logl  = np.array(hm['logl'])
            brights = np.nonzero(hm_logl > lmin)
            hm_logte = np.array(hm['logte'])
            try:
                hm_mag1 =np.array(hm[Filter1])
                hm_mag2 =np.array(hm[Filter2])
            except KeyError:
                hm_mag1 =np.array(hm[Filter1+'1'])
                hm_mag2 =np.array(hm[Filter2+'1'])
            hm_mag1= hm_mag1+dmod+Av*af275w
            hm_mag2= hm_mag2+dmod+Av*af336w
            hm_col = np.array(hm_mag1[1:])-np.array(hm_mag2[1:])
            if yfilter == Filter2:
                hm_mag = hm_mag2
            else:
                hm_mag = hm_mag1
            p = ax2.plot(smooth(hm_col[brights]),smooth(hm_mag[1:][brights]),color ='purple',lw=thick)
            p = ax1.plot(smooth(hm_logte[brights]),smooth(hm_logl[brights]),color ='purple',lw=thick)
            ax1.annotate(fmt%'P-AGB',xy=(5,3.8),color='purple',ha='left',fontsize=ftsize,**ann_kwargs)
            if 'phat' in file_ext:
                for arrow_pt in [10,20]:
                    pagb_arr = arrow_on_line(smooth(hm_logte[brights]),smooth(hm_logl[brights]),
                                             arrow_pt,'purple')
                    ax1.add_patch(pagb_arr)
                
                for arrow_pt in [10,20]:
                    pagb_arr = arrow_on_line(smooth(hm_col[brights]),smooth(hm_mag[1:][brights]),
                                             arrow_pt,'purple',width=0.2)
                ax2.add_patch(pagb_arr)
        if track == peagb:
            if yfilter == Filter2:
                ax2.annotate(fmt%'PE-AGB',xy=(.88,22.35-offset),color='red',ha='left',fontsize=ftsize,**ann_kwargs)
            else:
                if 'phat' in file_ext:
                    ax2.annotate(fmt%'PE-AGB',xy=(-0.3,21.5-offset),color='red',ha='left',fontsize=ftsize,**ann_kwargs)
                else:
                    ax2.annotate(fmt%'PE-AGB',xy=(-0.3,-1.1),color='red',ha='left',fontsize=ftsize,**ann_kwargs)
            ax1.annotate(fmt%'PE-AGB',xy=(4.5,3.17),color='red',ha='left',fontsize=ftsize,**ann_kwargs)
            ax1.plot(smooth(logte[4:-1]),smooth(logl[4:-1]),color='red',zorder=3,lw=thick)
            ax2.plot(smooth(t_col[4:-1]),smooth(tmag[4:-1]),color='red',zorder=3,lw=thick)
            if 'phat' in file_ext:
                for arrow_pt in [750,850,1200,1400]:
                    peagb_arr = arrow_on_line(smooth(logte[4:-1]),smooth(logl[4:-1]),
                                              arrow_pt,'red')
                    ax1.add_patch(peagb_arr)
                for arrow_pt in [850,1000,1320]:
                    peagb_arr = arrow_on_line(smooth(t_col[4:-1]),smooth(tmag[4:-1]),
                                              arrow_pt,'red',width=.2)
                    ax2.add_patch(peagb_arr)
        if track == agbm: 
            ax1.annotate(fmt%'AGB-Manqu\grave e',xy=(4.9,2.2),color='blue',ha='left',fontsize=ftsize,**ann_kwargs)
            ax1.plot(smooth(logte[4:-1]),smooth(logl[4:-1]),color='blue',zorder=3,lw=thick)
            ax2.plot(smooth(t_col[4:-1]),smooth(tmag[4:-1]),color='blue',zorder=3,lw=thick)
            if 'phat' in file_ext:
                for arrow_pt in [12,25,40]:
                    agbm_arr = arrow_on_line(smooth(logte),smooth(logl),arrow_pt,'blue')
                    ax1.add_patch(agbm_arr)
                for arrow_pt in [14,20,40]:
                    agbm_arr = arrow_on_line(smooth(t_col),smooth(tmag),arrow_pt,'blue',width=.2)
                    ax2.add_patch(agbm_arr)
        plt.savefig('hrd%i.png'%i)
        i+=1
    
    withdata = kwargs.get('withdata',False)
    if withdata == True:
        mag1 = kwargs['data_mag1']
        mag2 = kwargs['data_mag2']
        if yfilter == Filter2:
            mag = mag2
        else:
            mag = mag1
        ax2.plot(mag1-mag2,mag,'.',ms=1.5,mec=None,color='black')
    
    if 'phat' in file_ext:
        ax2.annotate(fmt%'AGB',xy=(1.4,29.6-offset),color='black',ha='left',fontsize=ftsize,**ann_kwargs)
        ax2.annotate(fmt%'AGB-Manqu\grave e',xy=(-.24,24.15-offset),color='blue',ha='left',fontsize=ftsize,**ann_kwargs)
        ax2.annotate(fmt%'P-AGB',xy=(-0.3,20.2-offset),color='purple',ha='left',fontsize=ftsize,**ann_kwargs)
        ax2.annotate(fmt%'ZAHB',xy=(-0.07,26.6-offset),color='darkgreen',ha='left',fontsize=ftsize,**ann_kwargs)
        xs = [-.8,2,2,-.8,-.8]
        ys = [25,25,19,19,25]
        ax2.plot(xs,ys,'--',lw=3,color='grey')
    else:
        #for galex
        ax2.annotate(fmt%'AGB-Manqu\grave e',xy=(-0.3,3),color='blue',ha='left',fontsize=ftsize,**ann_kwargs)
        ax2.annotate(fmt%'P-AGB',xy=(-0.3,-3),color='purple',ha='left',fontsize=ftsize,**ann_kwargs)
        ax2.annotate(fmt%'ZAHB',xy=(-0.3,4),color='darkgreen',ha='left',fontsize=ftsize,**ann_kwargs)
        ax2.set_ylim( (8, -4) )
        ax2.set_xlim( (-1, 5) )
    plt.savefig('hrd%s'%EXT)
    print 'wrote hrd%s'%EXT
    plt.close()
    return

def arrow_on_line(xarr,yarr,index,color,width=.06):
    x = xarr[index]
    dx = xarr[index+1]-x
    y = yarr[index]
    dy = yarr[index+1]-y
    arr = FancyArrow(x,y,dx,dy,fc=color,
                               ec="none",
                               linewidth=1,
                               length_includes_head=True,
                               head_width=width)
    return arr

def combine_two_lines(x1,y1,x2,y2):
    xy1_0 = np.array([x1[0],y1[0]])
    xy1_1 = np.array([x1[-1],y1[-1]])
    xy2_0 = np.array([x2[0],y2[0]])
    xy2_1 = np.array([x2[-1],y2[-1]])

    xp,yp = seg_intersect(xy1_0,xy1_1,xy2_0,xy2_1)
    
    xind1 = GenUtils.closest_match(xp,x1)
    xind2 = GenUtils.closest_match(xp,x2)
    
    new_xarr = np.hstack((x1[0:xind1],xp,x2[0:xind2][::-1]))
    new_yarr = np.hstack((y1[0:xind1],yp,y2[0:xind2][::-1]))
    return new_xarr,new_yarr

def perp( a ) :
    '''
    line segment intersection using vectors
    see Computer Graphics by F.S. Hill
    http://www.cs.mun.ca/~rod/2500/notes/numpy-arrays/numpy-arrays.html
    '''
    b = np.empty_like(a)
    b[0] = -a[1]
    b[1] = a[0]
    return b

def seg_intersect(a1,a2, b1,b2) :
    '''
    line segment a given by endpoints a1, a2
    line segment b given by endpoints b1, b2
    '''
    da = a2-a1
    db = b2-b1
    dp = a1-b1
    dap = perp(da)
    denom = np.dot( dap, db)
    num = np.dot( dap, dp )
    return (num / denom)*db + b1
    
def put_a_line_on_it(ax,val,consty=True,constx=False,filter1=None,filter2=None,color='black',ls='--',lw=None):
    """
    Trivial if this is just a constant y. But if the constant is actually on a 
    plot of y vs x-Const it's less trivial
    
    if consty == True, then this is just a constant y value.
    
    if yfilter != True, then val is assumed to be on a plot of y vs y-x and val represents
    a value of x.
    """
    xmin,xmax = ax.get_xlim()
    ymin,ymax = ax.get_ylim()

    xarr = np.linspace(xmin,xmax,20)
    # y axis is magnitude...
    yarr = np.linspace(ymax,ymin,20)
    empty_arr = np.zeros(len(xarr))
    
    if consty == True and constx == False:
        # just a contant y value over the plot range of x.
        new_yarr = empty_arr+val
        new_xarr = xarr
    if consty == False and constx == True:
        # a plot of y vs y-x and we want to mark 
        # where a constant value of x is
        new_yarr = yarr
        if filter2 == None:
            new_xarr = val - yarr
        if filter1 == None:
            new_xarr = yarr - val
    #if lw == None:
    #    ax.plot(new_xarr,new_yarr,ls,color=color)
    #else:
    #    ax.plot(new_xarr,new_yarr,ls,color=color,lw=lw)
    #ax.annotate('%.1f'%val,xy=(new_xarr[5],new_yarr[5]),xytext=(new_xarr[5],new_yarr[5]-0.2),
    #            horizontalalignment='left',**graphics.load_ann_kwargs())
    
    return new_xarr,new_yarr