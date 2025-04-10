#!/usr/bin/python3
# --------------------------------------------------------------------------------
# make_amali_daily_raw_netcdf.py
# --------------------------------------------------------------------------------
# 
# search amali files and convert to one flight-netcdf files
# 
# --------------------------------------------------------------------------------


from glob import glob 

from amali_raw_to_netcdf import amali_raw_to_netcdf

import xarray as xr
# --------------------------------------------------------------------------------

# data base is in:
# acloud campaign = 19.5. - 26.6.2017 
# db = '/data/obs/campaigns/acloud/p5/amali'

# aflux campaign 21.3. - 11.3.2019
#db = '/data/obs/campaigns/aflux/p5/amali'

# mosaic-aca campagin 31.8. - 13.9.2019
# db = '/data/obs/campaigns/mosaic-aca/p5/amali'

# ac3 campaign 20.3.2022 - ... preliminary directories ?
# db = '/data/obs/campaigns/halo-ac3/p5/amali'
# db = '/tmp/amali/'

# compex-ec campagin 2.4. - 16.4.2025
db = '/data/obs/campaigns/compex-ec/p5/amali'

search_path = db 
# for acloud data is in subdir raw for aflux there is no subdir ...
search_path = db+'/raw'

# destination is sub directory l00 (level zero zero)
dest_base = db+'/l00'
# dest_base = '/tmp'+'/l00'

# search days = directories
print( 'search days in', search_path )

# all days
fp = glob( search_path+'/????/??/??' )

# redo June 2017
# fp = glob( search_path+'/2017/06/??' )

# redo March 2019
# fp = glob( search_path+'/2019/03/??' )

# redo August 31 2020
# fp = glob( search_path+'/2020/08/31' )

# redo September 9 2020
# fp = glob( search_path+'/2020/09/10' )


# test single days
# fp = glob( search_path+'/2019/03/21' )
# intermittent a and b files on 2.6.2017 and 14.6.2017
# fp = glob( search_path+'/2017/06/02' )
# fp = glob( search_path+'/2017/06/14' )

# files from 31.3.2019 in 1.4. directory ... sorted away ... redo
# fp = glob( search_path+'/2019/04/01' )



# AC3 campaign March 2022
# March 9 2022 ... first data ?
# fp = glob( search_path+'/2022/03/09' )
# March 20 2022 ... first flight 
# fp = glob( search_path+'/2022/03/20' )
# March 22 2022 ... second flight 
# fp = glob( search_path+'/2022/03/22' )
# March 25 2022 ... third flight 
# fp = glob( search_path+'/2022/03/25' )

# March 28 2022 ... fourth flight 
# fp = glob( search_path+'/2022/03/28' )
# March 29 2022 ... fifth flight 
# fp = glob( search_path+'/2022/03/29' )
# March 30 2022 ... sixth flight 
# fp = glob( search_path+'/2022/03/30' )
# Aril 01 2022 ... seventh flight 
# fp = glob( search_path+'/2022/04/01' )
# Aril 04 2022 ... eigth flight 
# fp = glob( search_path+'/2022/04/04' )
# Aril 05 2022 ... ninth flight 
# fp = glob( search_path+'/2022/04/05' )
# Aril 07 2022 ... tenth flight 
# fp = glob( search_path+'/2022/04/07' )
# Aril 10 2022 ... eleventh flight 
# fp = glob( search_path+'/2022/04/10' )




print( 'found ', len(fp),'day directories' )


i = 0
for fp_i in fp :

  i += 1

  print( '--------------------------------------------------' )
  print( 'processing dir#',i,'=', fp_i )
#   data_files = 'a' + fp_i.split('/')[-3][2::] +fp_i.split('/')[-2].lstrip('0') +fp_i.split('/')[-1] + '*.*'
#   print(data_files)
  destination = dest_base + '/' + '/'.join(fp_i.split('/')[-3::])
#   print(destination)
#   import subprocess
#   # subprocess.run(['mkdir','-p',destination])
  amali_raw_to_netcdf( fp_i, dest_path=destination, verbose=0 )
  # amali_raw_to_netcdf( fp_i, dest_path='/tmp/', verbose=1 )
#   amali_raw_to_netcdf( fp_i, data_files, dest_path=destination, verbose=5 )
#   # amali_raw_to_netcdf( '/data/obs/campaigns/acloud/p5/amali/raw/2017/05/27/', '20170527a.mat',dest_path=dest_base+'/<YYYY>/<MM>/<DD>', verbose=5 )

#   ncfiles = glob( ).sort( )
#   for ncfile in ncfiles

