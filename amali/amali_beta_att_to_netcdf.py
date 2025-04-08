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
import pandas as pd
import datetime

import scipy.constants

import statistics

import time
import datetime as dt

import matplotlib
import matplotlib.pyplot as plt
import subprocess

# --------------------------------------------------------------------------------

def amali_beta_att_to_netcdf( 
      filename , # name of netcdf data file
      dest_path,
      mission,
      gpsp,
      snr_thres = 4.0, # threshold for snr , value below are set to NAN
      r_gnd_max = 4000., # maximum distance at which ground should appear (Polar 4/5 are not allowed higher than 4000m)
      verbose = 0
      ) :
    #/data/obs/campaigns/halo-ac3/p5/polar5_tmp/amali/l00/2022/03/29/amali_beta_att_l00_20220329_1230.nc
    ''' read AMALi raw nc data file and calculate beta_att'''
    
    if verbose > 0 : 
      print( 'try to open "', os.path.basename(filename), '" in "',os.path.dirname(filename),'"... ' )
      t0 = time.time()
    raw_data = xr.open_dataset(filename)
    print(raw_data.variables)

    if verbose > 0 : 
      t1 = time.time()
      print(' ...', t1-t0,'sec')
      print( 'raw_data=' )
      print( raw_data['channel_wvl'].data )

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
    if 'aperture' not in list(raw_data.variables): 
        raw_data["aperture"]= np.array(0.099)
    aperture_area = np.pi * (raw_data["aperture"].data / 2)**2 
    if verbose : print( 'aperture_area=', aperture_area )

    # laser energy per pulse in millijoule from Stachlewska et al 2010
    # pulse_energy = np.array([ 94., X, 94., X, 15, X ])*1e-3 
    if 'pulse_energy' not in list(raw_data.variables): 
        pulse_energy = np.array([ 94., 94., 94., 94., 15., 15.])*1e-3 
    else: 
        pulse_energy = raw_data["pulse_energy"].data * 1e-3
    if verbose > 5 : print( 'pulse_energy=', pulse_energy )

    # wavelength in meter
    if 'channel_wvl' not in list(raw_data.variables): 
        wvl = np.array([532., 532., 532., 532., 355., 355.])* 1e-9 
    else:
        wvl = raw_data['channel_wvl'].data * 1e-9 
    if verbose > 5 : print( 'wvl=', wvl )

    # list frequencies etc.
    #if verbose : 
      #print( 'channels:' )
      #for i_ch in range(N_chnl) : 
        #print( i_ch, 
               #'wvl=', raw_data["channel_wvl"].data[i_ch], 
               #'pow=', raw_data["pulse_energy"].data[i_ch], 
               #'pol=', raw_data["channel_pol"].data[i_ch], 
               #'d/a=', ('a' if raw_data['channel_analog'].data[i_ch] == 1 else 'd'), 
               #'Uh min...max=', np.min( raw_data['channel_Uh'].data ), 
               #          '...', np.max( raw_data['channel_Uh'].data ), 
               #'Ue min...max=', np.min( raw_data['channel_Ue'].data ), 
               #          '...', np.max( raw_data['channel_Ue'].data )
             #) 
   

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
    snr_all      = np.full( [N_chnl,N_time,N_bin_atm], float("NaN"), dtype=float )
    beta_att_all = np.full( [N_chnl,N_time,N_bin_atm], float("NaN"), dtype=float )

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
    #i_chnl = 0
    #i_time = int(N_time/2)
    # i_time = 6*60  # 6minutes after start
    # i_time = 60*60 # 60 mins
    # i_time = 65*60 # 65 mins
    # i_time = np.min([ 200*60 , N_time-1 - 100 ] ) # 200mins or 100mins before end

    #snr      = snr_all[     i_chnl,i_time,:]
    #beta_att = beta_att_all[i_chnl,i_time,:]

    z_km = r_range / 1000. * (-1. if (raw_data['angle'] == 0.) else +1. )

    #print(' type(time)         =', type(raw_data['time'].data[i_time]) )
    #print(' type(astype(dt.dt))=', type(raw_data['time'].data[i_time].astype(datetime.datetime)) )

    #t = raw_data['time'].data[i_time]
    #print(' type(t)            =', type(t) )
    #print(' type(astype(dt.dt))=', type(t.astype(datetime.datetime)) )


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
    #time_i_str = datetime.datetime.fromtimestamp( raw_data['time'].data[i_time].astype(datetime.datetime)*1e-9,tz=datetime.timezone.utc).strftime('%d.%m.%Y/%H:%M:%S.%f')[:-4]
    #time_0_str = datetime.datetime.fromtimestamp( raw_data['time'].data[  0   ].astype(datetime.datetime)*1e-9,tz=datetime.timezone.utc).strftime('%d.%m.%Y/%H:%M:%S.%f')[:-4]
                                              
    # time_start_str = datetime.datetime.fromtimestamp( raw_data['time'].data[0], tz=datetime.timezone.utc )
    # time_start_str = time_start_str.strftime( '%Y%m%d' ) 
    #time_start_str = datetime.datetime.fromtimestamp( raw_data['time'].data[0].astype(datetime.datetime)*1e-9, tz=datetime.timezone.utc ).strftime( '%Y%m%d' ) 

    #channel_info_str = str(int(raw_data['channel_wvl'].data[i_chnl])) + 'nm_' + raw_data['channel_pol'].data[i_chnl]   + '_'   + ('a' if raw_data['channel_analog'].data[i_chnl] == 1 else 'd')

    

    if verbose>=0 : 
      #print( 'shape(snr)=', np.shape(snr) )
      #print( 'snr[pre_trg]=', snr[0:5], '...', snr[i_pretrg-5:i_pretrg] )
      #print( 'snr[pre_trg]: min, max =', snr[0:i_pretrg].min(), snr[0:i_pretrg].max()  )
      #print( 'snr[rest   ]=', snr[i_pretrg:i_pretrg+5], '...', snr[N_bins-5:N_bins] )
      #print( 'snr[rest   ]: min, max =', snr[i_pretrg:N_bins].min(), snr[i_pretrg:N_bins].max() )
      #print( 'snr: min, max =', snr.min(), snr.max() )
      time_first_str=filename[-16:]
      nc_filename = 'amali_beta_att_l00_'+time_first_str
      # extract date from time_str
      YYYY = time_first_str[0:4]
      MM   = time_first_str[4:6]
      DD   = time_first_str[6:8]

        # if dest_path is given parse for <YYYY>, <MM>, etc. and evtl. create it
      if len(dest_path) > 0 :

            if verbose >= 10 : 
                print( 'time_first_str = ', time_first_str )
                print( 'YYYY MM DD = ', YYYY, MM, DD )
                print( '   dest_path=', dest_path )
       
            # replace <YYYY>, <MM>, <DD> by respective parts of the date
            dest_path = dest_path.replace('<YYYY>', YYYY ).replace('<MM>', MM ).replace('<DD>', DD )
            
            if verbose >= 10 : print( '=> dest_path=', dest_path )

            # create destination directory
            subprocess.run( [ 'mkdir', '-p', dest_path ] )

            # append '/' to path for use below
            if dest_path[-1] != '/' : dest_path += '/'
            
      if time_first_str=='20220309_1144.nc':flight_id='test flight'
      if time_first_str!='20220309_1144.nc':# this was just the test flight for which GPS does not exist      
          #load GPS_INS file for lon, lat , alt and get flight_id  
          GPS_data=[]
          GPS_filename= glob(gpsp+'/'+mission+'_P5_GPS_INS_'+YYYY+MM+DD+'_RF*')
          print(gpsp+'/'+mission+'_P5_GPS_INS_'+YYYY+MM+DD+'_RF*','FILENAME',GPS_filename)
          for i in range(0,len(GPS_filename)):
              GPS = xr.open_dataset(GPS_filename[i])
              #GPS_data = xr.open_mfdataset(GPS_filename)
              GPS_data.append(GPS)
              print('GPS',GPS)
              print('GPS_data',GPS_data)
          GPS_data=xr.concat(GPS_data, dim='time')    
          flight_id=GPS_filename[-1][-7:-3]
         
    #import pdb; pdb.set_trace()
    #t_plot = raw_data['time'].data
    #t_plot_0 = raw_data['time'].data[0]
    z_km_atm = z_km[i_pretrg:i_gnd_max]
    #z_km_range = np.array( [ np.min(z_km_atm), np.max(z_km_atm) ] )
    log_beta = np.log10( beta_att_all[:, :, : ])
    raw_data.time.attrs['units']='Seconds since 01.01.1970 00:00:00'
    raw_data=xr.decode_cf(raw_data)
    log_beta=xr.DataArray(data=log_beta,name='log_beta',dims=["i_channel", "time","height"], coords={'time':raw_data['time'],'height':z_km_atm*1000.,'i_channel':raw_data.channel.values})#dict(time=(["t"],raw_data['time'].data),height=(["z"], z_km_atm*1000.),i_channel=(["i"])
    snr=xr.DataArray(data=snr_all,name='snr',dims=["i_channel", "time","height"], coords={'time':raw_data['time'],'height':z_km_atm*1000.,'i_channel':raw_data.channel.values})#dict(time=(["t"],raw_data['time'].data),height=(["z"], z_km_atm*1000.),i_channel=(["i"])

    # --------------------------------------------------------------------------------
    # write to netcdf
    if verbose >= 1 : print( 'write data to "'+dest_path+nc_filename+'"' )
    da=xr.merge([log_beta,snr,raw_data.rename({'channel':'i_channel'})['channel_wvl'],raw_data.rename({'channel':'i_channel'})['channel_pol'], raw_data.rename({'channel':'i_channel'})['channel_analog']])
    if time_first_str!='20220309_1144.nc':    
        GPS_data=GPS_data.sel(time=pd.DatetimeIndex(da['time'].values), method='nearest')
        #coordinate von da
        GPS_data=GPS_data.assign_coords(time=da.time)
        #print(GPS_data.alt.data)
        da=xr.merge([da,GPS_data.lon,GPS_data.lat,GPS_data.alt ])
        #da['lon']=GPS_data.lon
        #da['lat']=GPS_data.lat
        #da['alt']=GPS_data.alt
        da.lon.attrs['standard_name']='longitude'
        da.lon.attrs['long_name']='WGS84 datum/longitude'
        da.lon.attrs['units']='degrees_east'
        
        da.lat.attrs['standard_name']='latitude'
        da.lat.attrs['long_name']='WGS84 datum/latitude'
        da.lat.attrs['units']='degrees_north'

        da.alt.attrs['standard_name']='altitude'
        da.alt.attrs['long_name']='aircraft flight altitude above mean sea level'
        da.alt.attrs['units']='m'    
    
    #set attributes
    da.attrs['institution']='Institute of Geophysics and Meteorology (IGM), University of Cologne; Alfred Wegener Institute - Research Unit Potsdam'
    da.attrs['source']='airborne observation'
    da.attrs['author']= 'Jan Schween (jschween@uni-koeln.de); Imke Schirmacher (imke.schirmacher@uni-koeln.de)'
    da.attrs['convention']='CF-1.8'
    da.attrs['featureType']='trajectory'
    da.attrs['mission']=mission
    da.attrs['platform']='Polar 5'
    da.attrs['flight_id']=flight_id
    da.attrs['title']='logarithmic volume attenuated backscatter coefficient based on AMALi observations onboard Polar 5 during '+mission
    da.attrs['instrument']='AMALi (Airborne Mobile Aerosol Lidar for Arctic research)'
    da.attrs['history']='corrected for background, range and incomplete overlap and reformatted (amali_raw_to_netcdf.py, amali_beta_att_to_netcdf.py)'
    da.attrs['contact']='imke.schirmacher@uni-koeln.de'
    da.attrs['created']=dt.datetime.utcnow().strftime('%d.%m.%Y/%H:%M/%Z')
    
    da.log_beta.attrs['standard_name']='logarithmic_volume_attenuated_backwards_scattering_function_in_air'
    da.log_beta.attrs['long_name']='logarithmic volume attenuated backscatter coefficent'
    da.log_beta.attrs['units']='m-1 sr-1'
    da.log_beta.attrs['description']='corrected for background, range and incomplete overlap and reformatted'

    da.snr.attrs['standard_name'] = 'signal_to_noise_ratio'
    da.snr.attrs['long_name'] = 'signal to noise ratio'
    da.snr.attrs['units'] = '-'
    da.snr.attrs['description'] = 'Signal to noise ratio'

    del(da.time.attrs["units"])
    da.time.attrs['standard_name']='time'
    da.time.attrs['long_name']='time in seconds since epoch'
    # da.time.attrs['units']='Seconds since 01.01.1970 00:00:00'
    
    da.height.attrs['standard_name']='altitude'
    da.height.attrs['long_name']='distance to aircraft'
    da.height.attrs['units']='m'
    
    da.i_channel.attrs['standard_name']='channel'
    da.i_channel.attrs['long_name']='AMALi channel number'
    da.i_channel.attrs['comment']='AMALi transmitts radiation with 532 nm (parallel and perpendicular polarized, both analog and digital) and 355 nm (non-polarized, analog and digital) wavelength'
    
    da.channel_wvl.attrs['standard_name']='radiation_wavelength'
    da.channel_wvl.attrs['long_name']='wavelength of transmitted radiation'
    da.channel_wvl.attrs['units']='nm'
    da.channel_wvl.attrs['comment']='AMALi transmitts radiation with 532 nm (parallel and perpendicular polarized, both analog and digital) and 355 nm (non-polarized, analog and digital) wavelength'
    
    da.channel_pol.attrs['standard_name']='radiation_polarization'
    da.channel_pol.attrs['long_name']='polarization of transmitted radiation'
    da.channel_pol.attrs['comment']='AMALi transmitts radiation with 532 nm (parallel and perpendicular polarized, both analog and digital) and 355 nm (non-polarized, analog and digital) wavelength'
    da.channel_pol.attrs['description']='parallel polarized: p; perpendicular polarized: s; non-polarized: n'
    
    da.channel_analog.attrs['standard_name']='signal_kind'
    da.channel_analog.attrs['long_name']='analog or digital signal'
    da.channel_analog.attrs['comment']='AMALi transmitts radiation with 532 nm (parallel and perpendicular polarized, both analog and digital) and 355 nm (non-polarized, analog and digital) wavelength'
    da.channel_analog.attrs['description']='analog: 1; digital: 0'
    da.to_netcdf( dest_path+nc_filename )
    # --------------------------------------------------------------------------------
    if verbose >= 1 : print( 'amali_raw_to_netcdf: done.' )      


from glob import glob 

# --------------------------------------------------------------------------------

# data base is in:
# acloud campaign = 19.5. - 26.6.2017 
# db = '/data/obs/campaigns/acloud/p5/amali/l00'
# gpsp='/data/obs/campaigns/acloud/p5/gps_ins'
# mission='ACLOUD'

# aflux campaign 21.3. - 11.3.2019
#db = '/data/obs/campaigns/aflux/p5/amali/l00'
#gpsp='/data/obs/campaigns/aflux/p5/gps_ins'
#mission='AFLUX'

# mosaic-aca campagin 31.8. - 13.9.2019
# db = '/data/obs/campaigns/mosaic-aca/p5/amali/l00'
# gpsp='/data/obs/campaigns/mosaic-aca/p5/gps_ins'
# mission='MOSAiC-ACA'

# halo-ac3 campaign 20.3.2022 - ... preliminary directories ?
# db = '/data/obs/campaigns/halo-ac3/p5/amali/l00'
# gpsp='/data/obs/campaigns/halo-ac3/p5/gps_ins'
# mission='HALO-AC3'

# compex-ec campaign 2.4.2025 - ... preliminary directories ?
db = '/data/obs/campaigns/compex-ec/p5/amali/l00'
gpsp='/data/obs/campaigns/compex-ec/p5/gps_ins'
mission='COMPEX-EC'

search_path = db 

# destination is sub directory l00 (level zero zero)
dest_base = db

# search days = directories
print( 'search days in', search_path )

# all days
# fp = glob( search_path+'/????/??/??/amali_l00_*' )
# fp = glob( search_path+'/2017/05/27/amali_l00_*' )
# fp = glob( search_path+'/2020/09/10/amali_l00_*' )
# fp = glob( search_path+'/2022/04/01/amali_l00_*' )
# fp = glob( search_path+'/2019/03/21/amali_l00_*' )
fp = glob( search_path+'/2025/04/07/amali_l00_*' )

print( 'found ', len(fp),'day directories' )


i = 0
if '/data/obs/campaigns/aflux/p5/amali/l00/2019/03/20/amali_l00_20190320_1701.nc' in fp:
    fp.remove('/data/obs/campaigns/aflux/p5/amali/l00/2019/03/20/amali_l00_20190320_1701.nc')

  
for fp_i in fp:

  i += 1

  print( '--------------------------------------------------' )
  print( 'processing dir#',i,'=', fp_i )

  #amali_raw_to_netcdf( fp_i, dest_path=dest_base+'/<YYYY>/<MM>/<DD>', verbose=5 )
  amali_beta_att_to_netcdf( fp_i,dest_base+'/<YYYY>/<MM>/<DD>', mission,gpsp,snr_thres=0, r_gnd_max=4500., verbose=0 )

# def amali_beta_att_to_netcdf( 
#       filename , # name of netcdf data file
#       dest_path,
#       mission,
#       gpsp,
#       snr_thres = 4.0, # threshold for snr , value below are set to NAN
#       r_gnd_max = 4000., # maximum distance at which ground should appear (Polar 4/5 are not allowed higher than 4000m)
#       verbose = 0
#       ) :
