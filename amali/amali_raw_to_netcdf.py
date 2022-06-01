# --------------------------------------------------------------------------------
# amali_raw_to_netcdf
# --------------------------------------------------------------------------------
# read AMALi raw data files based on search path, 
#   combine data into one xarray datset
# 
#   timeline
#   
#     29.03.2021 created by Jan Schween  (jschween@uni-koeln.de)
#      2.04.2021 use combine_by_coords( ... , data_vars='minimal' )
#     11.05.2021 added dest_path
#     18.05.2021 use shell command find with regular expression
# 
# 
# --------------------------------------------------------------------------------


import numpy as np

import xarray as xr

import datetime as dt

# from glob import glob

import subprocess

from  read_amali_raw  import read_amali_raw  


# --------------------------------------------------------------------------------

def amali_raw_to_netcdf( 
    search_path,    # search path of the raw data files
    dest_path = '', # path to which netcdf is written, you may use <YYYY>, <MM>, <DD> to include the first date
    verbose = 0     # 
    ) :
    '''read AMALi raw data files, combine data into one x-array, and write to netcdf '''

    if verbose >= 1 : print( 'amali_raw_to_netcdf: try to find "'+search_path+'"' )
  
    # glob is not specific enough
    # filenames = glob( search_path )
    # Amali raw files have very special names:
    #   a19325113600.748
    #   IYYMDDhhmmss.nnn
    #     I = identifier: a,b,c,...
    #     YYMDD = date as two digit year, one(!) digit month, two digit day
    #     hhmmss = time as two digit hour minute seconds
    #     nnn = three digit milliseonds (i guess) can even be four digits but then only 1000
    # => use shell command find with a regular expression - see test_find_regex.py
    find_files = subprocess.run( 
        [ 'find', search_path, '-regextype', 'sed', '-regex', '.*/[a-z][0-9]\{11\}\.[0-9]\{3,4\}' ],
        capture_output=True, 
        text=True
        )
    filenames = find_files.stdout.splitlines()

    # sort filenames
    filenames.sort(  )

  
    N_file = len( filenames )

    if verbose >= 5 : print( 'found', N_file, 'files, try to read ...' )

    if N_file > 0 :

        i_file = 0

        for fn_i in filenames :

            # if verbose >= 10 : print( 'call read_amali_raw(', fn_i, ')' )
            
            # read new data_set
            data_i = read_amali_raw( fn_i, verbose=verbose/2 )
     
            if len(data_i) > 0 :
                # build data_set 
                if i_file == 0 :
                    if verbose >= 10 : print( '  found first file ...' )
                    data = data_i
                
                    global_attrs = data.attrs
                else :
                    # combine new data_i with total data 
                    if verbose >= 10 : 
                        print( '  append data and data_i :'	 )
                        print( '    dims[time]=', data.dims['time'], data_i.dims['time'] )
                        print( '    dims[file_time] =', data.dims['file_time'], data_i.dims['file_time'] )
                        print( '    file_time = ...',            
                                                                 data[  'file_time'].data[-1],                  
                                                                 data_i['file_time'].data[ 0]
                                      )
                        print( '    file_time = ...', 
                                      dt.datetime.fromtimestamp( data[  'file_time'].data[-1], tz=dt.timezone.utc ).strftime( '%d.%m.%Y/%H:%M:%S.%f' )[:-4], 
                                      dt.datetime.fromtimestamp( data_i['file_time'].data[ 0], tz=dt.timezone.utc ).strftime( '%d.%m.%Y/%H:%M:%S.%f' )[:-4]
                                      )
                        print( '    file_name = ...', data['file_name'].data[-1], data_i['file_name'].data[0] )
                        print( '    file_name.len =', len(data['file_name'].data[-1]), len(data_i['file_name'].data[0]) )
                        print( '    signal.shape =', np.shape(data['signal']), np.shape(data_i['signal'].data) )

                        print( '    data.coords   =', data.coords   )
                        print( '    data_i.coords =', data_i.coords )
                        fnn = [ data['file_name'] , data_i['file_name'] ]


                    # join data and data_i
                    data = xr.combine_by_coords( [ data, data_i ] , data_vars='minimal' )


                if verbose >= 10 : print( "data.dims['time']=",data.dims['time'] )

                i_file += 1

                # end if len(data_i) > 0

            # end for fn_i in filenames

        if verbose >= 10 : 
            N_time = data.dims['time']
            N_chnl = data.dims['channel']
            print( 'read ',N_file,' files with in total ', N_time, 'records' )
            print( 'time: [0] ... [N-1] =', 
                   dt.datetime.fromtimestamp( data['time'].data[  0     ], tz=dt.timezone.utc ).strftime( '%d.%m.%Y/%H:%M:%S.%f' )[:-4], 
                   '...', 
                   dt.datetime.fromtimestamp( data['time'].data[N_time-1], tz=dt.timezone.utc ).strftime( '%d.%m.%Y/%H:%M:%S.%f' )[:-4]
                 )
            print( 'time: min ... max =', 
                   dt.datetime.fromtimestamp( np.min(data['time'].data[  0     ]), tz=dt.timezone.utc ).strftime( '%d.%m.%Y/%H:%M:%S.%f' )[:-4], 
                   '...', 
                   dt.datetime.fromtimestamp( np.max(data['time'].data[N_time-1]), tz=dt.timezone.utc ).strftime( '%d.%m.%Y/%H:%M:%S.%f' )[:-4]
                 )
            print( 'N_channels=', N_chnl  )
            for i_chnl in np.arange(N_chnl) :
              print( data['channel_wvl'].data[i_chnl], data['channel_pol'].data[i_chnl], data['channel_analog'].data[i_chnl], 
                     'min, max, mean, median=', 
                     np.min(    data['signal'].data[i_chnl,:,:]), 
                     np.max(    data['signal'].data[i_chnl,:,:]), 
                     np.mean(   data['signal'].data[i_chnl,:,:]), 
                     np.median( data['signal'].data[i_chnl,:,:])
                     )
            print( 'N_shots: min, max, mean, median=', 
                   np.min(    data['N_shot_beam'].data ), 
                   np.max(    data['N_shot_beam'].data ), 
                   np.mean(   data['N_shot_beam'].data ), 
                   np.median( data['N_shot_beam'].data ) 
                 )
                     

        # we need to take the attributes of the very first file
        data.attrs = global_attrs

        data.attrs['source'] = search_path

        data.attrs['time_range'] = (
            dt.datetime.fromtimestamp( min(data['time'].data), tz=dt.timezone.utc ).strftime( '%d.%m.%Y %H:%M:%S.%f' )[:-4] 
            + ' - ' + 
            dt.datetime.fromtimestamp( max(data['time'].data), tz=dt.timezone.utc ).strftime( '%d.%m.%Y %H:%M:%S.%f' )[:-4]
            )

        data.attrs['creator'] = 'amali_raw_to_netcdf.py  by  Jan Schween (jschween@uni-koeln.de)', 

        data.attrs['creation_date'] = dt.datetime.utcnow().strftime('%d.%m.%Y/%H:%M/%Z')


        time_first_str = dt.datetime.fromtimestamp( data['time'].data[0], tz=dt.timezone.utc ).strftime( '%Y%m%d_%H%M' )
        nc_filename = 'amali_l00_'+time_first_str+'.nc'

        # if dest_path is given parse for <YYYY>, <MM>, etc. and evtl. create it
        if len(dest_path) > 0 :

            # extract date from time_str
            YYYY = time_first_str[0:4]
            MM   = time_first_str[4:6]
            DD   = time_first_str[6:8]

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


        # --------------------------------------------------------------------------------
        # write to netcdf
        if verbose >= 1 : print( 'write data to "'+dest_path+nc_filename+'"' )
        data.to_netcdf( dest_path+nc_filename )
        # --------------------------------------------------------------------------------

        # end if N_file > 0

    if verbose >= 1 : print( 'amali_raw_to_netcdf: done.' )


    # end def amali_raw_to_netcdf

# --------------------------------------------------------------------------------


# test:
# YYYY = '2017'
# YY = YYYY[2:4]
# mm   = '05'
# m1   = mm[1]
# # dd   = '23'
# dd   = '25'
# HH   = '10'


# amali_raw_to_netcdf( '/data/obs/campaigns/acloud/amali/raw/'+YYYY+'/'+mm+'/'+dd+'/a'+YY+m1+dd+HH+'*.*' , verbose=5 )



# for h in np.arange( 8, 16 ) :
#      HH = "{:02d}".format(h)
#      #                     /data/obs/campaigns/acloud/amali/raw/  2017  /  05  /  23  /a  17  5 23 hh ...
#      amali_raw_to_netcdf( '/data/obs/campaigns/acloud/amali/raw/'+YYYY+'/'+mm+'/'+dd+'/a'+YY+m1+dd+HH+'*.*' , verbose=5 )




