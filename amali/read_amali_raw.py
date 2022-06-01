# --------------------------------------------------------------------------------
# read_amali_raw
# --------------------------------------------------------------------------------
# read AMALi raw data files
# 
#   timeline
#   
#     18.03.2021 created by Jan Schween  (jschween@uni-koeln.de my first python program !)
#     ...
#     02.04.2021 inverted order of dimensions: channel always first 
#     30.03.2022 changed to two date and time digits in time.unit 
#                  -> seconds since 1970-01-01 00:00:00
#                                        ^  ^  ^  ^  ^
#                                          
# 
# --------------------------------------------------------------------------------

import numpy as np

import os

import datetime

import xarray as xr


import matplotlib
import matplotlib.pyplot as plt

import scipy.constants # physical constants: speed_of_light, ...



# --------------------------------------------------------------------------------
def read_amali_raw( filename, plot_every=0, verbose=0 ):
    """read AMALi raw data files, return a structure with the data"""
#     filename :-)
#     verbose  message from program
#       0 -> no messages
#       1 -> start message : routine name and filename
#       5 -> details 
#     plot_every = N
#       plot every N'th set of profiles -> YYYYmmdd_HHMMSS.png
#       a file has 2min of data every sec or less => 120 profile sets 
#       => be aware that this may trash your disk ...
#       plotting is relativley slow and will reduce speed significantly !
#       use this only if you have problems with the data ... ;-)
# --------------------------------------------------------------------------------

    if verbose >= 1 : print( 'read_amali_raw: try to read "'+filename+'"' )

    # we need to convert mac_epoch (sec. since t0_mac=1.1.1904) to unix epoch (sec.since t0lx=1.1.1970)
    t0_mac = datetime.datetime( 1904, 1, 1, tzinfo=datetime.timezone.utc ).timestamp()

    # file size in bytes
    file_size = os.path.getsize( filename )


    # in case we plot : we dont want to plot on the screen => do not use x11 server , Agg instead ..
    if plot_every > 0 : matplotlib.use('Agg') 

    
    # open file for *r*eading as *b*inary
    with open( filename, "rb",  ) as f :
        
        # read first header as text line
        #   although the file is opened binary we can use .readline()
        #   but we get a bytearray and we have to convert it with .decode()
        #   then strip CR.LF (.rstrip()) and split along spaces (.split())
        if verbose >= 5 : print( 'try to read line1' )
        header1 = f.readline(  ).decode().rstrip().split()
        if verbose >= 5 : print( 'header1=', header1 )
        # write into dictionary
        header = { "op": header1[0] ,  # operator
                   "hw": header1[1] ,  # hardware
                   "sw": header1[2] }  # software
        if verbose >= 5 : print( 'header1=', header )
    
        if verbose >= 5 : print( 'try to read line2' )
        header2 = f.readline().decode().rstrip().split()
        if verbose >= 5 : print( header2 )
        header.update( 
                   { "res"       : float(header2[0]),  # resolution = range gatelength in meter
                     "qsdelay"   : float(header2[1]),  # q-switch delay = time after which laser is q-switched in microsec, should be 135microsec
                     "rate"      : float(header2[2]),  # pulse repetition rate in 1/sec
                     "pretrigger": float(header2[3]),  # time in microseconds before laser shot => use for noise level determination 
                     "angle"     : float(header2[4]),  # zenith angle in degrees 0 or 180deg for nadir or zenith measurement
                     "ntracks"   : int(  header2[5])   # number of channels
                   } )

        if verbose >= 5 : print( 'header=', header )
        
        # number of channels as handy int variable
        N_chnl = header["ntracks"]
        
        # read channel information   
        # define arrays for channel info
        channel_wvl = np.full(   [N_chnl], float("NaN"), dtype=np.float )

        # ncview does not like string type in file ... => should not use the following ...
        channel_pol = np.full(   [N_chnl], 'x', dtype=np.str   )
        # try type conversion of string to byte -> python believes at some point it must convert to int ???
        # channel_pol = np.full(   [N_chnl], bytes('x','utf-8'), dtype=np.byte   )
        # try encode  ... same error !?!?!?
        # channel_pol = np.full(   [N_chnl], 'x'.encode(), dtype=np.byte   )

        channel_Uh  = np.full(   [N_chnl], '-999',       dtype=np.int   )
        channel_Ue  = np.full(   [N_chnl], '-999',       dtype=np.int   )
        channel_analog = np.full([N_chnl], 255,          dtype=np.byte  )

        channel_info_str = []

        for i in range( N_chnl ) :
            # read channel information convert to string and split into substrings
            s = f.readline(  ).decode().rstrip()
            channel_info_str.append(s)
            # split into substrings
            s = s.split()
            # write into arrays 
            channel_wvl[i] = float(s[0]) 
            # ncview does not like string variables - we should not use this ... but cant avoid ...
            channel_pol[i] =       s[1]
            # pythont thinks at some point it has to use here int('p') - we cant use this ...
            # channel_pol[i] = bytes(s[1], 'utf-8' )
            # try encode ()
            # channel_pol[i] = s[1].encode()

            channel_Uh [i] = int(  s[2])
            channel_Ue [i] = int(  s[3])
            channel_analog[i] = int(  s[4])

            # end of read channel info loop

    

        # number of bins to estimate amount of data in file
        #   this number is given for every data record 
        #   it could in principle vary but probably not 
        #   Stachlewska et al 2010 mention that it should be 1000 and 1700 for Nadir and Zenith measurements, respectively
        N_bins = 1700 
  
        # data block size
        # every file should contain 2 min of data = 120sec ... or less 
        #   we try to estimate number of times from known sizes
        # ascii header has 27 characters plus two times CR+LF
        db_hdr_str_len = 2 + 27 + 2 
        db_sngl_chnl_len = db_hdr_str_len + 2*N_bins
        db_len = N_chnl * db_sngl_chnl_len
        #   number of records in file is thus remaining file size divided by size of one data record
        N_time_est = int( (file_size - f.tell()) / db_len )
        # we increase by 20% to be on the save side (wrong N_bins or header size ... )
        N_time_max = int( 1.2 * (N_time_est+1) )

        if verbose >= 5 : print( 'N_time_est=', N_time_est, ' -> we use:', N_time_max )

        # the core data 
        time    = np.full( [N_chnl,N_time_max], float("NaN"), dtype=np.float64 )
        N_shot  = np.full( [N_chnl,N_time_max], -999,         dtype=np.int16   )
        N_range = np.full( [N_chnl,N_time_max], -999,         dtype=np.int16   )
        res     = np.full( [N_chnl,N_time_max], float("NaN"), dtype=np.float )

        # received signal is stored binary as array of unsigned 2byte integer in little endian
        #   we ignore here the endianess as it is only important while reading
        #   we also use an order of dimensions which allows profile by profile saving
        signal  = np.full( [N_chnl,N_time_max,N_bins], 0, dtype=np.uint16 )

    
        # init indices to data 
        i_time = 0
        i_chnl = 0

        # buffer to read text lines
        s = 'not empty'

        # read data
        #   ... as long as there could be a complete data block
        while file_size - f.tell() >= db_sngl_chnl_len :

            if verbose >= 10 : print( (file_size - f.tell()) / db_sngl_chnl_len, ' blocks to go ...',  )
    
            # one empty line ... readline will return CR+LF thus it will have length 2
            s = f.readline(  )

            if len(s) > 0 :            
                # read data block header, convert to string, remove trailing CR LF and split along ' '
                db_hdr_str = f.readline(  ).decode().rstrip().split()

                if verbose >= 5 : print( f.tell(), '/', file_size, i_time, i_chnl, db_hdr_str, 
                                         datetime.datetime.fromtimestamp(float(db_hdr_str[0])+t0_mac, tz=datetime.timezone.utc ).strftime('%d.%m.%Y/%H:%M:%S.%f')[:-4]     )
        
                time[     i_chnl,i_time] = float(db_hdr_str[0])
                N_shot[   i_chnl,i_time] = int(  db_hdr_str[1])
                N_range[  i_chnl,i_time] = int(  db_hdr_str[2])
                res[      i_chnl,i_time] = float(db_hdr_str[3])
   
                N_words = N_range[ i_chnl, i_time ]

                # read BINARY data
                #   as array[N_words] of little endian unsigned 2 byte integer -> "<u2"
                signal_i = np.fromfile( f, dtype="<u2", count=N_words )
       
                # place signal data in array
                signal[ i_chnl, i_time, 0:N_words ] = signal_i

 
                # plot
                if plot_every > 0 :
                  if (i_time % plot_every == 0) and (i_chnl == N_chnl-1) :
                    if verbose >= 10 : print( 'plot' )
                    z = range( N_range[i_chnl,i_time] ) * res[i_chnl,i_time]
                    z_pt= 1./2. * scipy.constants.speed_of_light * header["pretrigger"]*1e-6 
                    color = [ '#00ff00', '#009900', '#0000ff', '#000099',  '#ff0000', '#990000' ]
                    xrange = [ 0.8 , 65536 ]
                    t_unx_i = datetime.datetime.fromtimestamp( time[i_chnl,i_time] + t0_mac, tz=datetime.timezone.utc )
                    plt.ioff()
                    fig, ax = plt.subplots( figsize=( 12,9 ) )
                    ax.set_title( t_unx_i.strftime( '%d.%m.%Y/%H:%M:%S.%f' )[:-4] )
                    ax.set_xlabel('P (mV , counts)')
                    ax.set_xscale('log')
                    ax.set_xlim( xrange )
                    ax.set_ylabel('z (m)')
                    for i in range(N_chnl) :
                        ax.plot( signal[i,i_time,:] , z, color[i], label=channel_info_str[i]  )
                    ax.plot( xrange  , 2 * [z_pt], 'k:' )
                    ax.legend()
                    plt.savefig( t_unx_i.strftime( '%Y%m%d_%H%M%S' )+'.png' )

                    # end plot script

        
                # increase channel index
                i_chnl += 1
                # if all channels are done increase time counter and reset channel index
                if i_chnl == N_chnl :
                    i_chnl  = 0
                    i_time += 1
                    # here could go a check whether size of arrays is going to be exceeded

                # end of len(s) > 0 block

            else: 
                print( 'empty data block ???\n', '- - - - - - - - - - - - - - - - - - - - - - - - - ' )
    
            # end of data read loop

        # end of with f do ...


    # post processing
        
    # number of records read
    N_time = i_time
    if verbose >= 5 : print('read', N_time, 'records ...' )

    # only if we found data
    if N_time > 0 : 

        # adapt size of arrays if less data was read
        if N_time < N_time_max : 
            if verbose >= 5 : 
                print( 'N_time < N_time_max  => resize ...' )
                print( 'signal.shape = ', signal.shape )

            time      = time[    :, 0:N_time ]
            N_shot    = N_shot[  :, 0:N_time ]
            N_range   = N_range[ :, 0:N_time ]
            res       = res[     :, 0:N_time ]
            signal    = signal[  :, 0:N_time, : ]

            if verbose >= 5 : print( 'signal.shape = ', signal.shape )
            # end if N_time < N_time_max

        # convert time to unix epoch, use only time for first channel
        t_unx = time[0,:] + t0_mac

        # range of lidar in meter 
        #   derived quantity  -> level 1 data
        #   it misses a correct zero point (determine pretrigger ...)
        # r_range = (0.5 + np.arange(0,N_bins)) * header["res"]

        # report maxima and minima
        if verbose >= 5 : 
            i = 0
            print( 't_beg=', time[0,i], t_unx[i], datetime.datetime.fromtimestamp( t_unx[i], tz=datetime.timezone.utc ).strftime( '%d.%m.%Y %H:%M:%S.%f' )[:-4] )

            i = N_time-1
            print( 't_end=', time[0,i], t_unx[i], datetime.datetime.fromtimestamp( t_unx[i], tz=datetime.timezone.utc ).strftime( '%d.%m.%Y %H:%M:%S.%f' )[:-4] )

            print( 't_unx : min, max, mean =', t_unx.min(), t_unx.max(), t_unx.mean()  )
            print( 't_unx : min, max, mean =', datetime.datetime.fromtimestamp( t_unx.min( ), tz=datetime.timezone.utc ).strftime( '%d.%m.%Y %H:%M:%S.%f' )[:-4], 
                                               datetime.datetime.fromtimestamp( t_unx.max( ), tz=datetime.timezone.utc ).strftime( '%d.%m.%Y %H:%M:%S.%f' )[:-4], 
                                               datetime.datetime.fromtimestamp( t_unx.mean(),tz=datetime.timezone.utc ).strftime( '%d.%m.%Y %H:%M:%S.%f' )[:-4]
                                               )

            for i_chnl in range(N_chnl) : 
                print( channel_info_str[i_chnl],': min, max, mean, median = ', 
                       signal[i_chnl,:,:].min(), 
                       signal[i_chnl,:,:].max(), 
                       signal[i_chnl,:,:].mean(), 
                       np.median( signal[i_chnl,:,:]) 
                       )


            print( 'variable shapes:' )
            print( 'wvl:', np.shape(channel_wvl) )
            print( 'pol:', np.shape(channel_pol) ) 
            print( 'Uh :', np.shape(channel_Uh ) ) 
            print( 'Ue :', np.shape(channel_Ue ) ) 
            print( 'analog:', np.shape(channel_analog ) ) 
            print( 't_unx:', np.shape(t_unx) ) 
            # print( 'r_range:', np.shape(r_range) ) 
            print( 'signal:', np.shape(signal) ) 

            # end report min and max ...


        if verbose > 5 : 
            print( 'time :', time.shape,  time )
            print( 't_unx:', t_unx.shape, t_unx )
            print( 'signal:', signal.shape )

        # create a datastructure of xarrays representing content of this file
        data = xr.Dataset( 
            data_vars = dict( 

                # instruments parameters taken from Stachlevska et al 2010
                aperture       = ([],  0.099,  dict( longname = 'telescope clear aperture diameter', units = 'm', comment = 'From Stachlevska et. al 2010' ) ),
                beam_diameter  = ([],  0.006,  dict( longname = 'laser beam diameter',               units = 'm', comment = 'From Stachlevska et. al 2010' ) ),
                beam_divergence= ([],  2.59,   dict( longname = 'laser beam divergence',          units = 'mrad', comment = 'From Stachlevska et. al 2010' ) ),
                pulse_duration = ([], 11.38,   dict( longname = 'laser pulse duration',             units = 'ns', comment = 'From Stachlevska et. al 2010' ) ),
                pulse_energy   = ( ["channel"], np.array([ 94., 94., 94., 94., 15., 15. ]), 
                          dict( longname = 'pulse energy per laser shot', unit = 'mJ', comment = 'From Stachlevska et. al 2010, if you have newer numbers replace these.' )
                          ), 

                # from file header
                res     = ( [],    header["res"],        dict( longname = 'resolution',            units = 'm' ) ), 
                qsdelay = ( [],    header["qsdelay"],    dict( longname = 'q-switch delay',        units = '1e-6s' ) ), 
                rate    = ( [],    header["rate"],       dict( longname = 'pulse repetition rate', units = '1/s' ) ), 
                pretrigger = ( [], header["pretrigger"], dict( longname = 'pretrigger time',       units = '1e-6s' ) ), 
                angle   = ( [],    header["angle"],      dict( longname = 'nadir distance angle',  units = 'deg', comment = '0 for nadir looking, 180 for zenith looking.\n The angle just gives the rough orientaton, for an exact angle you need the airplane data.' ) ), 

                # information about channels
                channel_wvl  = ( ["channel"], channel_wvl,    dict( longname='wavelength',   units='nm'       ) ),
                channel_pol  = ( ["channel"], channel_pol,    dict( longname='polarization',                      comment='p=parallel, s=senkrecht, n=none' ) ),
                channel_analog=( ["channel"], channel_analog, dict( longname='Flag for analog (mV) or digital (counts) measurement', comment='1=analog, 0=digital', _FillValue=255 ) ),
                channel_Uh   = ( ["channel","file_time"], np.transpose([channel_Uh]) ,     dict( longname='U_high',       units='V',            comment='Voltage for photomultiplier', _FillValue=-999 ) ),
                channel_Ue   = ( ["channel","file_time"], np.transpose([channel_Ue]) ,     dict( longname='U_e',          units='mV',           comment='Voltage for amplifier', _FillValue=-999 ) ),
                               
                # time coordinate of full dataset
                time       = ( ["time"      ],   t_unx,       dict( longname='time', units='seconds since 1970-01-01 00:00:00', calendar='standard', _FillValue=float("NaN") ) ),

                # first time in this file
                file_time  = ( ["file_time" ], [ t_unx[0] ],  dict( longname='start time of 2minute raw data files',   units='seconds since 1970-01-01 00:00:00', calendar='standard', _FillValue=float("NaN") ) ),

                # original filename
                file_name  = ( ["file_time" ], [ os.path.basename(filename).encode() ],   dict( longname='name of raw data file' ) ),


                # we do not provide range as it is given by res and pretrigger ...      
                # range = ( ["range"], r_range, dict( longname='distance', units='m', _FillValue=float("NaN") ) ),

                # information about beams
                #   a beam means the average of N_shots 
                #   these informations are given in the raw files for every data block 
                #   accordingly they have shape=["time", "channel"]
                #   but most of them depend not on channel and even not on time

                #   N_shot varies only with time but not with channel  (->shape=["time"])
                # N_shot_beam = ( ["time", "channel"], N_shot,   dict( longname='number of laser shots averaged', units='1', _FillValue=-999 ) ),
                N_shot_beam = ( ["time"], N_shot[0,:],   dict( longname='number of laser shots averaged', units='1', _FillValue=-999 ) ),


                #   res_beam varies not at all and it is also given above as scalar taken from the ascii header 
                # we drop it here
                # res_beam    = ( ["channel","time"], res,      dict( longname='range gate length', units='m', _FillValue=float("NaN") ) ),

                #   N_range may vary with [file_]time ... i would not recommend to change during flight 
                #   ... and this reading routine can not deal with that
                # N_range_beam= ( ["time"], [N_range[:,0]], dict( longname='range gate length', units='m', _FillValue=-999 ) ),


                # the core data
                #   we take max value of its type (=65535) as Fillvalue 
                signal = ( ["channel", "time", "range"], signal, dict( longname='received signal', units='mV or counts', comment = 'For identification of channels see variables channel_*.', _FillValue=np.iinfo(np.uint16).max ) )
                
                ),

            # global attributes
            attrs = dict( 
                description = 'AMALi Lidar raw data',
                source      = filename ,  
                operator    = header["op"], 
                hardware    = header["hw"], 
                software    = header["sw"],
                time_range  = ( datetime.datetime.fromtimestamp( t_unx[   0    ], tz=datetime.timezone.utc ).strftime( '%d.%m.%Y %H:%M:%S.%f' )[:-4]
                                +' - '+
                                datetime.datetime.fromtimestamp( t_unx[N_time-1], tz=datetime.timezone.utc ).strftime( '%d.%m.%Y %H:%M:%S.%f' )[:-4]
                              ),
                creator     = 'read_amali_raw.py  by  Jan Schween (jschween@uni-koeln.de)', 
                creation_date = datetime.datetime.today().strftime('%d.%m.%Y/%H:%M')
                )
            )




        # data.to_netcdf( 'test.nc' )

        if plot_every > 0 :
            fig, ax = plt.subplots( )
            ax.set_title( t_unx_i.strftime( '%d.%m.%Y/%H:%M:%S.%f' )[:-4] )
            ax.set_xlabel('time')
            ax.set_ylabel('signal')
            ax.set_ylim([ -20, 1200 ])
            for i in range(N_chnl) :
              ax.plot( t_unx, signal[i,:,600] , label=channel_info_str[i]  )
            ax.legend()
            plt.savefig( 'amali_raw_timeseries_'+datetime.datetime.fromtimestamp( t_unx[0], tz=datetime.timezone.utc ).strftime( '%Y%m%d_%H%M' )+'.png' )

        # end if N_time > 0
        
    else : # no data
        if verbose >= 5 : print( 'Found no data in file.' )
        data = []

    return data


    if verbose >= 1 : print( 'read_amali_raw: done' )

    # end read_amali_raw









# --------------------------------------------------------------------------------
# test

# filename = '/data/obs/campaigns/acloud/amali/raw/2017/05/23/a17523082037.491' # laser not running - dark currents
# filename = '/data/obs/campaigns/acloud/amali/raw/2017/05/23/a17523082437.439' # dark currents only
# filename = '/data/obs/campaigns/acloud/amali/raw/2017/05/23/a17523091457.715' # dark currents , changing ... ? 
# filename = '/data/obs/campaigns/acloud/amali/raw/2017/05/23/a17523091657.737' # dark currents only

# filename = '/data/obs/campaigns/acloud/amali/raw/2017/05/23/a17523135856.851'  # laser on, until 12:59:56+ , clear sky , dark currents stable 

# filename = '/data/obs/campaigns/acloud/amali/raw/2017/05/23/a17523102057.543' # laser on, varying dark currents, ...

# read_amali_raw( filename, plot_every=10, verbose=5 )


# middle of first ac3 flight on 20.3.2022
# a raw data file is 2minutes or 120 profiles => plot every 60 :-)
# filename = '/data/obs/campaigns/halo-ac3/polar5_tmp/amali/2022/03/20/a22320132832.266' 
# read_amali_raw( filename, plot_every=60, verbose=5 )


# --------------------------------------------------------------------------------
