#!/usr/bin/python3
# --------------------------------------------------------------------------------
# amali_raw_to_beta_att
# --------------------------------------------------------------------------------
# read AMALi raw nc data file
#   an calculate beta_att
#
# 
#   timeline
#   
#     29.03.2021 created by Jan Schween  (jschween@uni-koeln.de)
#     ...
#     30.03.2022 removed unnecessary comments
# 
# --------------------------------------------------------------------------------


import os

import xarray as xr

import numpy as np

import datetime

import scipy.constants

import statistics

import time

import matplotlib
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------------

def amali_raw_to_beta_att( 
      filename , # name of netcdf data file
      snr_thres = 4.0, # threshold for snr , value below are set to NAN
      r_gnd_max = 4000., # maximum distance at which ground should appear (Polar 4/5 are not allowed higher than 4000m)
      verbose = 0 
      ) :
    ''' read AMALi raw nc data file and calculate beta_att'''
    
    if verbose > 0 : 
      print( 'try to open "', os.path.basename(filename), '" in "',os.path.dirname(filename),'"... ' )
      t0 = time.time()
    raw_data = xr.open_dataset( filename )
    if verbose > 0 : 
      t1 = time.time()
      print(' ...', t1-t0,'sec')
      print( 'raw_data=' )
      print( raw_data )

    # number of time steps
    N_time = raw_data.dims['time']
    # number of range bins
    N_bins = raw_data.dims['range']
    # number of channels
    N_chnl = raw_data.dims['channel']

    # resolution = delta r ...
    res = raw_data['res'].data

    # range zero:
    #   the laser is triggered a pre-trigger-time t_pretrg after start of the transient recorder
    #   from this time on the laser pulse is on its way, t_pre is in millisec
    t_pretrg = raw_data["pretrigger"].data * 1e-6
    r_pretrg = 1/2 * scipy.constants.speed_of_light * t_pretrg 
    # convert to index
    i_pretrg = int( r_pretrg / res )

    # full overlap between telescope field of view and laser beam after a distance r_ovrlap 
    #   according to Birtes Text Amali_202008.pdf
    #   import from external file ?
    r_ovrlap = 300.
    # convert to index
    i_ovrlap = i_pretrg + int( r_ovrlap / res )
    # ground max = maximum distance at which ground should appear
    i_gnd_max = i_pretrg + int( r_gnd_max / res )

    # number of nodes in atmosphere
    #   from node ia = i_pretrg to ib = i_gnd_max-1 => N = ib-ia+1 =  ...
    N_bin_atm = i_gnd_max - i_pretrg

    if verbose : 
      print( 'Nbins=', N_bins, 'res=', res, 'r_prtrg=', r_pretrg )

    # determine range with r=0 after end of t_pretrg
    r_range     = (0.5 + np.arange(0,N_bins   )) * res - r_pretrg
    r_atm_range = (0.5 + np.arange(0,N_bin_atm)) * res           
    if verbose : 
      print( 'r_range: min, max =', r_range.min(), r_range.max() )
      print( 'r_range =', r_range )  
      print( 'r_atm_range: min, max =', r_atm_range.min(), r_atm_range.max() )


    # aperture = diameter of free telescope opening in meter
    aperture_area = np.pi * (raw_data["aperture"].data / 2)**2 
    if verbose : print( 'aperture_area=', aperture_area )

    # laser nergy per pulse in millijoule from Stachlewska et al 2010
    # pulse_energy = np.array([ 94., X, 94., X, 15, X ])*1e-3 
    pulse_energy = raw_data["pulse_energy"].data * 1e-3
    if verbose > 5 : print( 'pulse_energy=', pulse_energy )

    # wavelength in meter
    wvl = raw_data['channel_wvl'].data * 1e-9 
    if verbose > 5 : print( 'wvl=', wvl )

    # list frequencies etc.
    if verbose : 
      print( 'channels:' )
      for i_ch in range(N_chnl) : 
        print( i_ch, 
               'wvl=', raw_data["channel_wvl"].data[i_ch], 
               'pow=', raw_data["pulse_energy"].data[i_ch], 
               'pol=', raw_data["channel_pol"].data[i_ch], 
               'd/a=', ('a' if raw_data['channel_analog'].data[i_ch] == 1 else 'd'), 
               'Uh min...max=', np.min( raw_data['channel_Uh'].data ), 
                         '...', np.max( raw_data['channel_Uh'].data ), 
               'Ue min...max=', np.min( raw_data['channel_Ue'].data ), 
                         '...', np.max( raw_data['channel_Ue'].data )
             ) 
   

    # select specific channels and times ..
    i_chnl_info = 0
    # i_time_info = 2525
    # i_time_info = 1000
    # i_time_info = 3000
    i_time_info =  0

  
    if verbose : 
      print('read signal from xarray into a np.array. Python needs here some time (> 16sec !) ... ')
      print('  we wait here for ')
      print('    signal = raw_data["signal"].data') 
      print('  to finish !')
      print('  Independent on whether only a sub range or all of signal.data is accessed' ) 
      print('  it takes a long time. Any subsequent access to another part is very fast.' )
      print('  => python reads somehow inefficiently the whole array and then accesses the parts.' )
      print('  => we load here the whole data and enjoy meanwhile this text :-)' )
      print('     (btw. old IDL needs 0.16sec for this operation...)' )
      t0 = time.time()
    signal = raw_data['signal'].data
    # signal = raw_data['signal'].data.to_numpy()  no to_numpy() ! although it is said https://xarray.pydata.org/en/stable/generated/xarray.DataArray.as_numpy.html
    if verbose : 
      t1 = time.time()
      print(' ...', t1-t0,'sec - you see python is slow here !')

    # extract signal in pre trigger range, ground range 
    signal_pretrg = signal[:,:,0:i_pretrg-1]
    signal_gndmax = signal[:,:,i_gnd_max:N_bins-1]
    

    # determine pretrigger background, and its noise = standard deviation 
    # background signal seems to increase linearily during measurement
    # see plot profiles_*.png generated below
    # => determine background noise as mean in pretrigger range in ground range
    # => interpolate linear with height between these two meands
    # => we can calculate without a loop over channel and time by using axis parameter of nanmean function
    # background from mean in pretrigger region (for all channels and all times )
    bkgnd = np.nanmean( signal_pretrg, 2 )
    bg_sd = np.nanstd(  signal_pretrg, 2 )
    # the same for the 'in-ground' signal = bgg = *b*ack*g*round-in-*g*round
    bgg_m = np.nanmean( signal_gndmax, 2 )
    bgg_s = np.nanstd(  signal_gndmax, 2 )
    # mean ranges of both regions
    r_bg_pretrg = np.nanmean( r_range[0:i_pretrg] )
    r_bg_gndmax = np.nanmean( r_range[i_gnd_max:N_bins] )
    # slope of background , for all channels and all times
    bg_slope = (bkgnd-bgg_m)/(r_bg_pretrg-r_bg_gndmax)
   
    # linear increase generates a variance of 1/12 (a dr)^2 
    #   (see https://atmos.meteo.uni-koeln.de/ag_crewell/doku.php?id=internal:administration:jan:linear_regression:variance_of_functions)
    # we subtract this variance to get noise of background signal alone
    bg_noise = np.sqrt( bg_sd**2 - 1./12.* (bg_slope * (r_range[i_pretrg-1]-r_range[0]))**2 )
    if verbose : 
      print( 'bg_stdev=', bg_sd )
      print( 'bg_noise=', bg_noise )

    # median and theoretical Poisson noise = sqrt(signal)
    s_med = np.nanmedian( signal_pretrg, 2 )
    n_psn = np.sqrt(s_med)
   

    # snr
    #   = signal to noise ratio
    #     see eg Heese et al 2010 (https://amt.copernicus.org/articles/3/1763/2010/amt-3-1763-2010.html)
    #     we measure 
    #       signal = atm_bksct + bkgnd 
    #     and calculate 
    #       atm_bsct = signal - bkgnd
    #     signal has noise from bkgnd and from photomultiplier
    #     photon counting of atm_bksct is poisson distributed and thus var = atm_bsckt
    #     but also includes noise of bkgnd. Variances must be added if they are not correlated:
    #       var(signal) = (signal-bkgnd) + bkgnd_noise
    #     we subtract bkgnd i.e. its noise is introduced again thus:
    #       var(signal-bkgnd) = (signal-bkgnd) + 2*bkgnd_noise
    #     signal to noise ratio is accordingly:
    #       snr = (signal-bkgnd) / sqrt(2*bkgnd_noise^2 + (signal-bkgnd))

    # allocate arrays 
    snr_all      = np.full( [N_chnl,N_time,N_bin_atm], float("NaN"), dtype=np.float )
    beta_att_all = np.full( [N_chnl,N_time,N_bin_atm], float("NaN"), dtype=np.float )

    # bkgnd and bkgnd_noise depend only on time and channel but not on height
    # we can do this for all times and all channels at once because we have bakgnd and noise also for all times and channels
    # we only have to cycle about height ...



    # convert shotpower [Joule] in measurement units 
    #   = photon counts for the odd channels = photoncount channels
    #   = in mV for the even channels 
    #       -> we need a conversion mV_to_counts 
    #       ... see plot cnts_vs_mv.png  
    #       -> fit for all data or every x minutes a curve to counts vs mV 
    #       ... slope for small mV gives a relation
    #       x minute could be 2minutes = raw data files = ~120sec = ~120 indices
    #       Fit could be just linear for U < 100 ?
    # ----------------------------------------------------------------
    # for simplicity I estimate here from cnts_vs_mv.png 
    # although this is for only one channel (532nm a p)
    mV_per_photon = 1.0 ; 250./100. 
    # THIS should be adapted by a fit of counts to milliVolts !!!!!!
    # ----------------------------------------------------------------
    for i in range(int(N_chnl/2)) : 
      photon_counts_per_shot = pulse_energy[2*i] / ( scipy.constants.Planck * scipy.constants.speed_of_light / wvl[2*i] )
      pulse_energy[2*i+1] = photon_counts_per_shot
      pulse_energy[2*i  ] = photon_counts_per_shot * mV_per_photon
    if verbose > 5 : print( 'pulse_energy:', pulse_energy )

    # lidar constant according to Stachlevska et al 2010
    C_lidar = pulse_energy * res * aperture_area
    if verbose : print( 'C_lidar:', C_lidar )


    # ----------------------------------------------------------------
    # here comes beta_att ...
    # ----------------------------------------------------------------
    # step 1: subtract background and do range correction per range gate
    for i_range in range(N_bin_atm) :
      # signal_i = signal[:,:,i_range+i_pretrg] - bkgnd[:,:] 
      signal_i = ( signal[:,:,i_range+i_pretrg] - ( bkgnd + bg_slope*(r_atm_range[i_range]-r_bg_pretrg)) )
      # snr_all[     :,:,i_range] = abs( signal_i ) / np.sqrt( 2*bg_sd**2 + abs( signal_i ) )
      snr_all[     :,:,i_range] = abs( signal_i ) / np.sqrt( 2*bg_noise**2 + abs( signal_i ) )
      beta_att_all[:,:,i_range] = signal_i * r_atm_range[i_range]**2

    # step 2: apply lidar constant per channel
    for i_chnl in range( N_chnl ) :
      beta_att_all[i_chnl,:,:] = beta_att_all[i_chnl,:,:] / C_lidar[i_chnl]

    # step 3: remove all beta_att at low snr values
    if snr_thres > 0 : 
      i_low_snr = np.where( snr_all[:,:,:] < snr_thres )
      beta_att_all[i_low_snr] = float('nan')
    # ----------------------------------------------------------------


    # we plot only certain times and channels
    i_chnl = 0
    i_time = int(N_time/2)
    # i_time = 6*60  # 6minutes after start
    # i_time = 60*60 # 60 mins
    # i_time = 65*60 # 65 mins
    # i_time = np.min([ 200*60 , N_time-1 - 100 ] ) # 200mins or 100mins before end

    snr      = snr_all[     i_chnl,i_time,:]
    beta_att = beta_att_all[i_chnl,i_time,:]

    z_km = r_range / 1000. * (-1. if (raw_data['angle'] == 0.) else +1. )


    print(' type(time)         =', type(raw_data['time'].data[i_time]) )
    print(' type(astype(dt.dt))=', type(raw_data['time'].data[i_time].astype(datetime.datetime)) )

    t = raw_data['time'].data[i_time]
    print(' type(t)            =', type(t) )
    print(' type(astype(dt.dt))=', type(t.astype(datetime.datetime)) )


    # time_i_str = datetime.datetime.fromtimestamp( raw_data['time'].data[i_time], tz=datetime.timezone.utc ).strftime('%d.%m.%Y/%H:%M:%S.%f')[:-4]

    # Python and datetime-variables ... we are in HELL !
    #   python wants to enforce clear logic code. 
    #   So here we are:
    #     in our netcdf is time given as 'seconds since 1970-1-1 0:0:0'
    #     xarray is clever and converts this to numpy.datetime64 a fancy datetime type 
    #     Sounds good. But we want to have here a formatted string of the form YYYYMMDD 
    #     and numpy.datetime64 does not know how to freely format a date !
    #     No problem! we convert with numpy.datetime64.astype(datetime.datetime) to type datetime.datetime 
    #     which has function strftime() to create formatted time strings.
    #     this works fine in my test program test_np_datatime.py
    #     But(!) numpy.datetime64 knows different resolutions of time. And if you get as far as nanoseconds
    #     resolution the conversion numpy.datetime64.astype(datetime.datetime) suddenly generates int as type
    #     - no datetime ! - and what does xarray ? it provides my 1sec resolution data in nanosec resoltuion !
    #     No problem: we assume this int is nanoseconds and we can convert it by multiplying with 1e-9 in seconds 
    #     and provide it to function datetime.datetime.fromtimestamp( ) to a type datetime.datetime varaible which knows
    #   
    time_i_str = datetime.datetime.fromtimestamp( raw_data['time'].data[i_time].astype(datetime.datetime)*1e-9,tz=datetime.timezone.utc).strftime('%d.%m.%Y/%H:%M:%S.%f')[:-4]
    time_0_str = datetime.datetime.fromtimestamp( raw_data['time'].data[  0   ].astype(datetime.datetime)*1e-9,tz=datetime.timezone.utc).strftime('%d.%m.%Y/%H:%M:%S.%f')[:-4]
                                              
    # time_start_str = datetime.datetime.fromtimestamp( raw_data['time'].data[0], tz=datetime.timezone.utc )
    # time_start_str = time_start_str.strftime( '%Y%m%d' ) 
    time_start_str = datetime.datetime.fromtimestamp( raw_data['time'].data[0].astype(datetime.datetime)*1e-9, tz=datetime.timezone.utc ).strftime( '%Y%m%d' ) 

    channel_info_str = str(int(raw_data['channel_wvl'].data[i_chnl])) + 'nm_' + raw_data['channel_pol'].data[i_chnl]   + '_'   + ('a' if raw_data['channel_analog'].data[i_chnl] == 1 else 'd')

    

    if verbose : 
      print( 'shape(snr)=', np.shape(snr) )
      print( 'snr[pre_trg]=', snr[0:5], '...', snr[i_pretrg-5:i_pretrg] )
      print( 'snr[pre_trg]: min, max =', snr[0:i_pretrg].min(), snr[0:i_pretrg].max()  )
      print( 'snr[rest   ]=', snr[i_pretrg:i_pretrg+5], '...', snr[N_bins-5:N_bins] )
      print( 'snr[rest   ]: min, max =', snr[i_pretrg:N_bins].min(), snr[i_pretrg:N_bins].max() )
      print( 'snr: min, max =', snr.min(), snr.max() )


    # --------------------------------------------------------------------------------
    # plot
    print( 'plot ...' )
    # we dont want to plot on the screen => do not use x11 server , Agg instead ..
    matplotlib.use('Agg') 


    font = {'family' : 'Latin Modern Roman',
            'weight' : 'normal',
            'size'   : 22}
    matplotlib.rc('font', **font)

  
    # ----------------------------------------
    # plot Nr1: backscatter profile 
    print( 'plot backscatter profile ...' )

    beta_range = 	[ 1e-9, 1e-6 ]
    y_range = np.array([ -r_gnd_max , +200 ]) / 1000.

    fig, ax = plt.subplots( figsize=( 15,9 ) )

    # time_i = datetime.datetime.fromtimestamp( raw_data['time'].data[i_time], tz=datetime.timezone.utc )
    time_i = raw_data['time'].data[i_time]

    ax.set_title( 
      # time_i.strftime( '%d.%m.%Y/%H:%M:%S.%f' )[:-4] + ' ' + 
      time_i_str + ' ' + 
      str(int(raw_data['channel_wvl'].data[i_chnl])) + 'nm '+ 
              raw_data['channel_pol'].data[i_chnl] + ' ' +
      ('anlg.' if raw_data['channel_analog'].data[i_chnl] == 1 else 'digtl.')
      )
    # latex in x-axis title, note the r in front of the string !
    #    if you want to see error lists omit it !!!
    ax.set_xlabel(r'$\beta_{att}$')
    ax.set_xlim( beta_range )
    ax.set_xscale('log')
    ax.set_ylabel('z (km)')
    ax.set_ylim( y_range )
  
    ax.plot( beta_range, [0]*2,        'k'  ) # r=0 line
    ax.plot( beta_range, [z_km[i_ovrlap]]*2, 'k'  ) # r=r_overlap line
    ax.plot( [1        ]*2, y_range,      'k'  ) # beta = 1 line
    ax.plot( [snr_thres]*2, y_range,      'k:' ) # SNR = 5 line

    # plot beta_att
    # use of LaTex for annotations: note the r in front of the string !
    ax.plot( beta_att, z_km[i_pretrg:i_gnd_max] , 'b', label=r'$\beta_{att}$')



    # Noise profile --- derived from SNR
    # SNR = (sig-bkgnd)/noise 
    # beta = (sig-bkgnd)/C*r^2
    #   => beta/SNR = noise/C*r^2
    ax.plot( beta_att/snr, z_km[i_pretrg:i_gnd_max], 'r', label='noise' )

    ax.legend()

    if raw_data['angle'] == 0 : # downward looking lidar
        ax.text( beta_range[0], 0, r' $\uparrow\uparrow$ pretrigger $\uparrow\uparrow$', 
                 verticalalignment ='bottom', 
                 horizontalalignment ='left'
                )
        ax.text( beta_range[0], -r_ovrlap/1000., r' $\downarrow\downarrow$ overlap $\downarrow\downarrow$', 
                 verticalalignment ='top', 
                 horizontalalignment ='left'
               )
    else :
        ax.text( beta_range[0],  0, r' $\downarrow\downarrow$ pretrigger $\downarrow\downarrow$', 
                 verticalalignment ='top', 
                 horizontalalignment ='left'
                )
        ax.text( beta_range[0], +r_ovrlap/1000., r' $\uparrow\uparrow$ overlap $\uparrow\uparrow$', 
                 verticalalignment ='bottom', 
                 horizontalalignment ='left'
               )

    plt.savefig( 'beta_att_profile_'+channel_info_str+'_'+time_start_str+'.png', bbox_inches='tight' )
    # after saving a figure we can close it ...
    plt.close(fig)


    # ----------------------------------------
    # plot Nr 2: counts versus mV 

    print( 'plot counts vs mV ...' )

    # only first N seconds
    N = 60 # N_time-i_time # 60

    i_chnl_0 = int(i_chnl/2)*2
    i_chnl_1 = i_chnl_0+1

    fig, ax = plt.subplots( figsize=( 15,9 ) )
    ax.set_title( 
      # datetime.datetime.fromtimestamp( raw_data['time'].data[i_time  ], tz=datetime.timezone.utc ).strftime( '%d.%m.%Y/%H:%M:%S' ) + 
      # raw_data['time'].data[i_time].strftime( '%d.%m.%Y/%H:%M:%S' ) + 
      time_i_str + 
      ' - ' +
      # datetime.datetime.fromtimestamp( raw_data['time'].data[i_time+N], tz=datetime.timezone.utc ).strftime( '%H:%M:%S' )  +
      # raw_data['time'].data[i_time+N].strftime( '%H:%M:%S' )  +
      datetime.datetime.fromtimestamp( raw_data['time'].data[i_time+N].astype(datetime.datetime)*1e-9,tz=datetime.timezone.utc).strftime( '%H:%M:%S' ) + 
      '  : ' + 
      str(int(raw_data['channel_wvl'].data[i_chnl_1])) + 'nm '+ 
              raw_data['channel_pol'].data[i_chnl_0] + ' ' +
      ('anlg.' if raw_data['channel_analog'].data[i_chnl_1] == 1 else 'digt.') + 
      ' vs ' +
      str(int(raw_data['channel_wvl'].data[i_chnl_0])) + 'nm '+ 
              raw_data['channel_pol'].data[i_chnl_0] + ' ' +
      ('anlg.' if raw_data['channel_analog'].data[i_chnl_0] == 1 else 'digt.')
      )

    # x-axis title and limits    
    ax.set_xlabel('Signal (mV)')
    # ax.set_xlim( [ bkgnd[i_chnl,i_time]-3*bg_sd[i_chnl_0,i_time] , 2000 ] )
    ax.set_xlim( [ -3*bg_sd[i_chnl_0,i_time] , 2000 ] )


    # guiding lines
    # bkgnd and bkgnd+/-noise at time = i_time
    ax.plot( [    0                ]*2, [0,65535], 'k' )  
    ax.plot( [+bg_sd[i_chnl_0,i_time]]*2, [0,65535], 'k:' )
    ax.plot( [-bg_sd[i_chnl_0,i_time]]*2, [0,65535], 'k:' )

    # bkgnd and bkgnd+/-noise at time = i_time
    ax.plot( [0,65535], [        0              ]*2, 'b' )  
    ax.plot( [0,65535], [+bg_sd[i_chnl_1,i_time]]*2, 'b:' )
    ax.plot( [0,65535], [-bg_sd[i_chnl_1,i_time]]*2, 'b:' )


    # y-axis title and limits
    ax.set_ylabel('Signal (counts)')
    ax.set_ylim( [ -10 , 200 ] )
        
    # digital and analog signal in atmosphere
    signal_anlg_i = np.copy( signal[i_chnl_0,i_time:i_time+N, i_pretrg:i_gnd_max  ] )
    signal_cont_i = np.copy( signal[i_chnl_1,i_time:i_time+N, i_pretrg:i_gnd_max  ] )

    print( 'shape(signal_anlg_i)=',np.shape(signal_anlg_i) )
    print( 'shape(signal_cont_i)=',np.shape(signal_anlg_i) )

    for i in np.arange(N_bin_atm) : 
        signal_anlg_i[:,i] -= bkgnd[i_chnl_0,i_time:i_time+N]
        signal_cont_i[:,i] -= bkgnd[i_chnl_1,i_time:i_time+N]

    # points with high SNR
    i_good_snr = np.where( (snr_all[i_chnl_0,i_time:i_time+N,:] >= snr_thres)  )

    print( 'shape(i_good_snr)=',np.shape(i_good_snr) )


    # plot all counts vs analog in black
    ax.plot( signal_anlg_i,             signal_cont_i,             'k.'   )
    # plot only points with high SNR in red
    ax.plot( signal_anlg_i[i_good_snr], signal_cont_i[i_good_snr], 'r.'   )

    # an idea for a 'fit':
    #   two parameter function
    #     P_cnt = P_sat*(1-exp(-ln(2)*P_mV/P_hlf))
    #   reaches at P_mV = P_hlf the value P_cnt = P_sat/2 ...
    # 
    # visual estimates :
    # P_sat = 87.0
    # P_hlf = 65.
    P_sat = 80.
    P_hlf = 70.

    N_fit = 100
    P_anlg_a = bkgnd[i_chnl_0,i_time]
    P_anlg_b = 2000.
    P_fit_in = (P_anlg_b-P_anlg_a)*np.arange(N)/N
    ln2 = np.log(2)

    # try to fit the function:
    #     P_fit = P_sat*( 1 - exp(-ln2*P_sig/P_hlf) )
    #   find P_sat by averaging over high signal points
    #   and then fit 
    #     ln( 1 - P_cnt/P_sat )  to   P_sig 
    #   which should deliver -ln2/P_hlf as slope
    #  
    try_count_mV_fit = 0
    if try_count_mV_fit :
      # find photon counts of saturation 
      i_sat = np.where( (snr_all[i_chnl_0,i_time:i_time+N,:] >= snr_thres) & (signal_anlg_i > 1500.) ) # & (signal_anlg_i < 2000.)  )
      P_sat = np.mean( signal_cont_i[i_sat] )
      y = np.log( 1 - signal_cont_i[i_good_snr]/P_sat )
      x = signal_anlg_i[i_good_snr]
      m,b = np.polyfit( x, y, 1 )

      P_hlf = -ln2/m

      print('==================================================================================')
      print( 'P_sat=', P_sat, 'P_hlf=', P_hlf, 'counts_per_mv=', P_sat*ln2/P_hlf )

      ax.text( 0 , ax.get_ylim()[1]*0.85, str(P_sat)+' '+ ":0.3f".format(P_hlf)+' -> '+':0.3f'.format(P_sat*ln2/P_hlf) ) 

      # end try_count_mV_fit


    P_fit = P_sat * (1.-np.exp(-ln2*P_fit_in/P_hlf))
    
    ax.plot( P_fit_in, P_fit, 'g' )

   

    plt.savefig( 'cnts_vs_mv_'+channel_info_str+'_'+time_start_str+'.png', bbox_inches='tight' )
    plt.close(fig)



    # ----------------------------------------
    # plot Nr. 3: different profiles
    print( 'plot bkgnd, signal and beta profiles ...' )

    fig, ax = plt.subplots( 1,3, figsize=( 16,9 ) )

    z_km_range = np.array( [ z_km.min() , z_km.max() ] )
    bg_range = bkgnd[i_chnl,i_time]+np.array([-1,1])*8*bg_sd[i_chnl,i_time]

    ax[0].set_xlabel('Bkgnd.')
    ax[0].set_xlim( bg_range )
    ax[0].set_ylabel('z (km)')
    ax[0].set_ylim( z_km_range )

    # mark bkgnd mean and range
    ax[0].plot( [bkgnd[i_chnl,i_time]                     ]*2  ,z_km_range , 'k' )
    ax[0].plot( [bkgnd[i_chnl,i_time]-bg_sd[i_chnl,i_time]]*2  ,z_km_range , 'k:' )
    ax[0].plot( [bkgnd[i_chnl,i_time]+bg_sd[i_chnl,i_time]]*2  ,z_km_range , 'k:' )

    ax[0].plot( [bgg_m[i_chnl,i_time]                     ]*2  ,z_km_range , 'r' )
    ax[0].plot( [bgg_m[i_chnl,i_time]-bgg_s[i_chnl,i_time]]*2  ,z_km_range , 'r:' )
    ax[0].plot( [bgg_m[i_chnl,i_time]+bgg_s[i_chnl,i_time]]*2  ,z_km_range , 'r:' )

    # mark r_pretrg=0, r_overlap and r_ground 
    ax[0].plot( bg_range, [0]*2, 'k' )
    ax[0].plot( bg_range, [z_km[i_ovrlap ]]*2, 'k' )
    ax[0].plot( bg_range, [z_km[i_gnd_max]]*2, 'k' )

    # signal in blue
    ax[0].plot( signal[i_chnl,i_time,:], z_km , 'b' )

    # plot straight line connecting both mean backgrounds ...
    zp_m = np.mean(z_km[0:i_pretrg])
    zg_m = np.mean(z_km[i_gnd_max:N_bins])
    m_bg = (bkgnd[i_chnl,i_time] - bgg_m[i_chnl,i_time])/(zp_m-zg_m)
    ax[0].plot( bkgnd[i_chnl,i_time] + m_bg*(z_km_range-zp_m), z_km_range , 'g' )
    ax[0].plot( bkgnd[i_chnl,i_time]-bg_noise[i_chnl,i_time] + m_bg*(z_km_range-zp_m), z_km_range , 'g:' )
    ax[0].plot( bkgnd[i_chnl,i_time]+bg_noise[i_chnl,i_time] + m_bg*(z_km_range-zp_m), z_km_range , 'g:' )



    ax[1].set_title( 
           # time_i.strftime( '%d.%m.%Y/%H:%M:%S.%f' )[:-4] 
           time_i_str + ' ' + channel_info_str )

    sig_range = np.array( [ bg_range[0] , 1e5 ] )

    ax[1].set_xlabel('Signal')
    ax[1].set_xlim( sig_range )
    ax[1].set_xscale('log')
    # ax[1].set_ylabel('z (km)')
    ax[1].set_ylim( z_km_range )


    # lines at bkgnd and bkgnd+/- noise 
    ax[1].plot( [bkgnd[i_chnl,i_time]                       ]*2, z_km_range , 'k' )
    ax[1].plot( [bkgnd[i_chnl,i_time]+bg_sd[i_chnl,i_time]  ]*2, z_km_range , 'k:' )
    ax[1].plot( [bkgnd[i_chnl,i_time]-bg_sd[i_chnl,i_time]  ]*2, z_km_range , 'k:' )

    # plot signal
    ax[1].plot( signal[i_chnl,i_time,:] , z_km, 'b' )

    # linear bkgnd interpolation in green
    ax[1].plot( bkgnd[i_chnl,i_time] + m_bg*(z_km_range-zp_m), z_km_range , 'g' )

    # mark z=0, overlay and maximum ground distance
    ax[1].plot( sig_range, [0]*2 , 'k' )
    ax[1].plot( sig_range, [z_km[i_ovrlap]]*2 , 'k' )
    ax[1].plot( sig_range, [z_km[i_gnd_max]]*2 , 'k' )

  

    beta_range = np.array( [ 1e-10 , 1e-6 ] )
    ax[2].set_xlabel(r'$\beta_{att}$')
    ax[2].set_xlim( beta_range )
    ax[2].set_xscale('log')
    # ax[2].set_ylabel('z (km)')
    ax[2].set_ylim( z_km_range )

    # mark noise line
    ax[2].plot( bg_sd[i_chnl,i_time]/C_lidar[i_chnl] * r_range**2, z_km , 'k:' )

    # mark r_pretrg=0, r_overlap and r_ground 
    ax[2].plot( beta_range, [0]*2, 'k' )
    ax[2].plot( beta_range, [z_km[i_ovrlap ]]*2, 'k' )
    ax[2].plot( beta_range, [z_km[i_gnd_max]]*2, 'k' )

    ax[2].plot( beta_att_all[i_chnl,i_time,:], z_km[i_pretrg:i_gnd_max], 'b' )


    plt.savefig( 'profiles_'+channel_info_str+'_'+time_start_str+'.png', bbox_inches='tight'  )
    plt.close(fig)



    # ----------------------------------------
    # plot nr 4: time height section
    print( 'plot beta t*r section ...' )

    fig, ax = plt.subplots( 3, 1, figsize=( 16,12 ) )

    # plt.subplots_adjust( bottom=0.1, right=0.1, top=0.1)
    # plt.subplots_adjust( left=0.0, right=0.0 )



    # how do i convert the array of timestamps into type datetime ???
    # oh holy python:
    #   because i use datetime as coordinate matplotlib generates a warning
    #   -->"FutureWarning: Using an implicitly registered datetime converter for a matplotlib plotting method. ..."
    #   so i convert to datetime and mathplotlib converts it back ...
    #   but I only want correct nice axis annotations ...
    #   it might be an idea to build a formatter object which can deal with unix epoch (or whatever ...)
    print( 'get time in minutes since start of file ...' )
    # t_plot = (raw_data['time'].data - raw_data['time'].data[0] ) / 60.
    t_plot = raw_data['time'].data
    # t_plot = np.asarray( t_plot, dtype='datetime64[s]')
    

    # t_plot_0 = datetime.datetime.fromtimestamp( raw_data['time'].data[0], tz=datetime.timezone.utc )
    t_plot_0 = raw_data['time'].data[0]




    ax[0].set_title(  
        # t_plot_0.strftime( '%d.%m.%Y/%H:%M:%S.%f' )[:-4] + ' ' + 
        time_0_str + ' ' +
        str(int(raw_data['channel_wvl'].data[i_chnl])) + 'nm '+ 
                raw_data['channel_pol'].data[i_chnl] + ' ' +
        ('anlg.' if raw_data['channel_analog'].data[i_chnl] == 1 else 'digtl.')
        )

    x_tick_format = matplotlib.dates.DateFormatter("%H:%M")
    ax[0].xaxis.set_major_formatter(x_tick_format)

    z_km_atm = z_km[i_pretrg:i_gnd_max]
    z_km_range = np.array( [ np.min(z_km_atm), np.max(z_km_atm) ] )

    ax[0].set_ylabel( 'z (km)' )
    ax[0].set_ylim( z_km_range )

    ax[0].margins(0.0)
    ax[1].margins(0.0)
    ax[2].margins(0.0)

    N_lvl = 20
    # lvl_a = -1.
    # lvl_b = +2.
    lvl_a = -9.0
    lvl_b = -5.0

    log_beta = np.log10( beta_att_all[ i_chnl, :, : ])
    lb_min = np.nanmin(log_beta[np.isfinite(log_beta)])
    lb_max = np.nanmax(log_beta[np.isfinite(log_beta)])
    lb_perc = np.nanpercentile( log_beta, [ 0, 1, 5, 10, 25, 50, 75, 90, 95, 99, 100  ] )
    print( 'min(lb)=', lb_min )
    print( 'max(lb)=', lb_max )
    print( 'perc(lb)=', lb_perc )
    # lvl_a = lb_min
    lvl_a = lb_perc[1] # 1% percentile
    # lvl_a = lb_perc[2] # 5% percentile
    # lvl_a = lb_perc[3] # 10% percentile
    # lvl_b = lb_perc[7] # 90% percentile
    # lvl_b = lb_perc[8] # 95% percentile
    lvl_b = lb_perc[9] # 99% percentile
    # lvl_b = lb_max

    levels = lvl_a + (lvl_b-lvl_a)*np.arange( N_lvl )/(N_lvl-1)
    colors = np.arange( N_lvl )/(N_lvl-1)

    print('contour or colormesh or ...' )
    t0 = time.time()

    # with the dimensions here (large N_time, N_height, N_lvl) contourf becomes VERY SLOW (33sec !)
    # cs = ax[0].contourf( t_plot, r_range/1000, np.transpose( log_beta ) , 
    #                      levels=levels 
    #                      # colors=colors
    #                     )

    # an alternativ could be imshow( ) ... which does boxes for every value
    # it is really fast: 0.4sec
    # but it takes different paramters and we have to adapt ...
    #    imshow( Z, extent=[xmin,xmax,ymin,ymax], ... )
    # cs = ax[0].imshow( np.transpose( log_beta ) ) , 
    #                    # extent = [ min(t_plot), min(r_range/1000), max(t_plot), max(r_range/1000) ] ,
    #                    extent = [ min(t_plot), max(t_plot), min(r_range/1000), max(r_range/1000) ] ,
    #                    aspect = 'auto', 
    #                    vmin = lvl_a, 
    #                    vmax = lvl_b 
    #                  )

    # another alternative pcolormesh (or eventually pcolor)
    # fast: 0.6sec
    # parameters similar (but not identical) to contourf 
    cs = ax[0].pcolormesh( t_plot, z_km_atm, np.transpose( log_beta ) , 
                           vmin = lvl_a, 
                           vmax = lvl_b 
                          )

    t1 = time.time()
    print(' ...', t1-t0,'sec')


    cb = fig.colorbar( cs, ax=ax[:], shrink=0.5 )
    # cb.set_label( r'$\log_{10}(\beta_{att}/\beta_{ovl})$' )
    cb.set_label( r'$\log_{10}(\beta_{att})$' )


    ax[1].xaxis.set_major_formatter(x_tick_format)
    # ax[1].set_xlabel( 'time (Unix epoch :-) )' )

    bg_ns_perc = np.nanpercentile( bg_sd[i_chnl,:], [ 0, 1, 5, 10, 25, 50, 75, 90, 95, 99, 100  ] )
    bg_mean = np.nanmean( bkgnd[i_chnl,:] )

    print( 'back ground noise perc.=', bg_ns_perc )

    bg_min = bg_mean-3*bg_ns_perc[9] # 99% percentile 
    bg_max = bg_mean+3*bg_ns_perc[9] # 99% percentile


    ax[1].set_ylabel( 'bkgnd '+('(mV)' if raw_data['channel_analog'].data[i_chnl] == 1 else '(counts)') )

    ax[1].set_ylim( [ bg_min, bg_max ] )


    ax[1].plot( t_plot, bkgnd[i_chnl,:],                   'k'   )
    ax[1].plot( t_plot, bkgnd[i_chnl,:] + bg_sd[i_chnl,:], 'g:'  )
    ax[1].plot( t_plot, bkgnd[i_chnl,:] - bg_sd[i_chnl,:], 'g:'  )
    

    ax[2].set_ylabel( r'$log_{10}(\beta_{ovl})$' )
    ax[2].plot( t_plot, log_beta[:,i_ovrlap-i_pretrg],        'k'   )

    ax[2].set_ylim( [ lvl_a, lvl_b ] )

    ax[2].xaxis.set_major_formatter(x_tick_format)

    # ax[2].set_xlabel('minutes since start' )
    ax[2].set_xlabel('time (UTC)' )


    print('  save fig ...')
    t0 = time.time()

    plt.savefig( 'beta_att_'+channel_info_str+'_'+time_start_str+'.png', bbox_inches='tight'  )
    plt.close(fig)

    t1 = time.time()
    print('...', t1-t0, 'sec' )


    # ----------------------------------------
    # plot nr 5: contours of background
    print( 'plot background contours (colourmesh) ...' )
    t0 = time.time()

    fig, ax = plt.subplots( figsize=( 16,8 ) )

    ax.set_title(  
        # t_plot_0.strftime( '%d.%m.%Y/%H:%M:%S.%f' )[:-4] + ' ' + 
        time_0_str + ' ' +        
        str(int(raw_data['channel_wvl'].data[i_chnl])) + 'nm '+ 
                raw_data['channel_pol'].data[i_chnl] + ' ' +
        ('anlg.' if raw_data['channel_analog'].data[i_chnl] == 1 else 'digtl.')
        )

    # ax.set_xlabel('minutes since start')
    ax.set_xlabel('time (UTC)')
    ax.xaxis.set_major_formatter(x_tick_format)

    ax.set_ylabel('z (km)')

    cs = ax.pcolormesh( t_plot, z_km, np.transpose( signal[i_chnl,:,:] ) , 
                           vmin = bkgnd[i_chnl,i_time]-8*bg_sd[i_chnl,i_time], 
                           vmax = bkgnd[i_chnl,i_time]+8*bg_sd[i_chnl,i_time]
                          )

    cb = fig.colorbar( cs )
    cb.set_label( 'signal around bkgnd '+('(mV)' if raw_data['channel_analog'].data[i_chnl] == 1 else '(counts)') )

    plt.savefig( 'bkgnd_'+channel_info_str+'_'+time_start_str+'.png', bbox_inches='tight'  )
    plt.close(fig)

    t1 = time.time()
    print('...', t1-t0, 'sec' )


    # ----------------------------------------
    # plot nr 6: contours of SNR
    print( 'plot SNR contours (colormesh) ...' )

    t0 = time.time()

    fig, ax = plt.subplots( figsize=( 16,8 ) )

    ax.set_title(  
        # t_plot_0.strftime( '%d.%m.%Y/%H:%M:%S.%f' )[:-4]  + ' ' + 
        time_0_str + ' ' +
        str(int(raw_data['channel_wvl'].data[i_chnl])) + 'nm '+ 
                raw_data['channel_pol'].data[i_chnl] + ' ' +
        ('anlg.' if raw_data['channel_analog'].data[i_chnl] == 1 else 'digtl.')
        )

    # ax.set_xlabel('minutes since start')
    ax.set_xlabel('time (UTC)')
    ax.xaxis.set_major_formatter(x_tick_format)

    ax.set_ylabel('z (km)')

    cs = ax.pcolormesh( t_plot, z_km_atm, np.transpose( snr_all[i_chnl,:,:] ) , 
                           vmin = 0, 
                           vmax = 10
                          )

    cb = fig.colorbar( cs )
    cb.set_label( 'SNR' )

    plt.savefig( 'snr_'+channel_info_str+'_'+time_start_str+'.png', bbox_inches='tight'  )
    plt.close(fig)

    t1 = time.time()
    print('...', t1-t0, 'sec' )


    # ----------------------------------------
    # plot nr 7: depolarization, and color ratio 
    print( 'plot depolarization and color ratio (colormesh) ...' )
    t0 = time.time()

    fig, ax = plt.subplots( 2, 1, figsize=( 16,8 ) )

    ax[0].set_title(  
        # t_plot_0.strftime( '%d.%m.%Y/%H:%M:%S.%f' )[:-4]
        time_0_str
        )

    ax[0].xaxis.set_major_formatter(x_tick_format)

    ax[0].set_ylabel('z (km)')


    depol = np.log10(beta_att_all[ 2,:,: ]) - np.log10(beta_att_all[ 0,:,: ])
    depol_perc = np.nanpercentile( depol, [ 0, 1, 5, 10, 25, 50, 75, 90, 95, 99, 100  ] )
    print( 'depol percentiles:', depol_perc )
    depol_min = depol_perc[1]
    depol_max = depol_perc[9]

    cs = ax[0].pcolormesh( t_plot, z_km_atm, np.transpose( depol ) , 
                           vmin = depol_min, 
                           vmax = depol_max
                          )

    cb = fig.colorbar( cs, ax=ax[0] )
    # cb.set_label( 'depol.ratio '+str(int(raw_data["channel_wvl"].data[2]))+'nm '+raw_data["channel_pol"].data[2]+
    #                       '  / '+str(int(raw_data["channel_wvl"].data[0]))+'nm '+raw_data["channel_pol"].data[0] )
    cb.set_label( r'$\log_{10}(\beta_{'+str(int(raw_data['channel_wvl'].data[2]))+r','+raw_data["channel_pol"].data[2]+r'}/\beta_{'+raw_data["channel_pol"].data[0]+r'})$' )

    # cb.set_label( 'depol.ratio' )


    col_ratio = np.log10(beta_att_all[ 0,:,: ]) - np.log10(beta_att_all[ 5,:,: ])
    col_ratio_perc = np.nanpercentile( col_ratio, [ 0, 1, 5, 10, 25, 50, 75, 90, 95, 99, 100  ] )
    print( 'col.ratio percentiles:', col_ratio_perc )
    col_ratio_min = col_ratio_perc[1]
    col_ratio_max = col_ratio_perc[9]

    cs = ax[1].pcolormesh( t_plot, z_km_atm, np.transpose( col_ratio  ) , 
                           vmin = col_ratio_min, 
                           vmax = col_ratio_max
                          )

    cb = fig.colorbar( cs, ax=ax[1]  )
    # cb.set_label( 'color.ratio '+str(int(raw_data["channel_wvl"].data[0]))+'nm /'+
    #                              str(int(raw_data["channel_wvl"].data[5]))+'nm' )
    cb.set_label( r'$\log_{10}(\beta_{'+str(int(raw_data["channel_wvl"].data[0]))+r'}/\beta_{'+str(int(raw_data["channel_wvl"].data[5]))+r'})$' )
    # cb.set_label( 'color.ratio' )
    
    # ax[1].set_xlabel('minutes since start')
    ax[1].set_xlabel('time (UTC)')
    ax[1].xaxis.set_major_formatter(x_tick_format)

                        
    ax[1].set_ylabel('z (km)')



    plt.savefig( 'depol_ratio_'+time_start_str+'.png', bbox_inches='tight'  )
    plt.close(fig)

    t1 = time.time()
    print('...', t1-t0, 'sec' )


    print('amali_raw_to_beta_att: done')
    print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ')




# --------------------------------------------------------------------------------
# test


# amali_raw_to_beta_att( 'amali_L00_20170523_1000.nc' ) # last ~15minutes noise (laser off)
# amali_raw_to_beta_att( 'amali_L00_20170523_1100.nc' ) # this is only white noise (laser off)
# amali_raw_to_beta_att( 'amali_L00_20170523_1200.nc' )
# amali_raw_to_beta_att( 'amali_L00_20170523_1300.nc' )

# amali_raw_to_beta_att( 'amali_L00_20170525_0901.nc' )
# amali_raw_to_beta_att( 'amali_L00_20170525_1001.nc' )
# amali_raw_to_beta_att( 'amali_L00_20170525_1101.nc' )
# amali_raw_to_beta_att( 'amali_L00_20170525_1201.nc' )

# AC3 campaign March 2022

# first data from 09.03.2022  ... test run of instrument
# amali_raw_to_beta_att( '/data/obs/campaigns/halo-ac3/polar5_tmp/amali/l00/2022/03/09/amali_l00_20220309_1144.nc', snr_thres=0.0, r_gnd_max = 3500., verbose=1 )


# first flight from 20.03.2022 ... long section of closed cloud deck 2-2.5km below plane
# instruments data need to be included in raw data netcdf !
# amali_raw_to_beta_att( '/data/obs/campaigns/halo-ac3/polar5_tmp/amali/l00/2022/03/20/amali_l00_20220320_1142.nc', snr_thres=0.0, r_gnd_max = 3500., verbose=1 )

# second flight from 22.03.2022 ... big variety of clouds, ci, ac, cs?, aerosol ...
# amali_raw_to_beta_att( '/data/obs/campaigns/halo-ac3/polar5_tmp/amali/l00/2022/03/22/amali_l00_20220322_1304.nc', snr_thres=0.0, verbose=1 )

# third flight from 25.03.2022 ... big gap with no data, ~30min of strongly varying clout tops
# amali_raw_to_beta_att( '/data/obs/campaigns/halo-ac3/polar5_tmp/amali/l00/2022/03/25/amali_l00_20220325_1227.nc', snr_thres=0.0, r_gnd_max = 3500., verbose=1 )

# fourth flight from 28.03.2022 ... stronllgy varying clout tops at 2+/-0.5km
# amali_raw_to_beta_att( '/data/obs/campaigns/halo-ac3/polar5_tmp/amali/l00/2022/03/28/amali_l00_20220328_1239.nc', snr_thres=0.0, verbose=1 )

# fifth flight from 29.03.2022 ... low clouds 500m above surface
# amali_raw_to_beta_att( '/data/obs/campaigns/halo-ac3/polar5_tmp/amali/l00/2022/03/29/amali_l00_20220329_1230.nc', snr_thres=0.0, verbose=1 )

# sixth flight from 30.03.2022 ... broken clouds between 2.3+/-0.3km below plane
# amali_raw_to_beta_att( '/data/obs/campaigns/halo-ac3/polar5_tmp/amali/l00/2022/03/30/amali_l00_20220330_1148.nc', snr_thres=0.0, verbose=1 )

# seventh flight from 1.4.2022 ... stronlgy varying cloud tops, repeatedly crossed
# amali_raw_to_beta_att( '/data/obs/campaigns/halo-ac3/polar5_tmp/amali/l00/2022/04/01/amali_l00_20220401_0925.nc', snr_thres=0.0, verbose=1 )
# amali_raw_to_beta_att( '/data/obs/campaigns/halo-ac3/polar5_tmp/amali/l00/2022/04/01/amali_l00_20220401_0925.nc', snr_thres=1.5, r_gnd_max=3500., verbose=1 )

# eigth flight from 4.4.2022 ... clout tops increasing/decreasing from 0.3...1km, several times crossed
# amali_raw_to_beta_att( '/data/obs/campaigns/halo-ac3/polar5_tmp/amali/l00/2022/04/04/amali_l00_20220404_1019.nc', snr_thres=0.0, verbose=1 )
# amali_raw_to_beta_att( '/data/obs/campaigns/halo-ac3/polar5_tmp/amali/l00/2022/04/04/amali_l00_20220404_1019.nc', snr_thres=1.5, r_gnd_max=3500., verbose=1 )

# ninth flight from 5.4.2022 ... hazy(?) clouds ~1km above ocean
# amali_raw_to_beta_att( '/data/obs/campaigns/halo-ac3/polar5_tmp/amali/l00/2022/04/05/amali_l00_20220405_1025.nc', snr_thres=0.0, verbose=1 )
# amali_raw_to_beta_att( '/data/obs/campaigns/halo-ac3/polar5_tmp/amali/l00/2022/04/05/amali_l00_20220405_1025.nc', snr_thres=0, r_gnd_max=3500., verbose=1 )


# tenth flight from 7.4.2022 ... 
# amali_raw_to_beta_att( '/data/obs/campaigns/halo-ac3/polar5_tmp/amali/l00/2022/04/07/amali_l00_20220407_0950.nc', snr_thres=0, r_gnd_max=3500., verbose=1 )

# eleventh flight from 10.4.2022 ... 
amali_raw_to_beta_att( '/data/obs/campaigns/halo-ac3/p5/polar5_tmp/amali/l00/2022/04/01/amali_l00_20220401_0925.nc', snr_thres=0, r_gnd_max=3500., verbose=1 )

