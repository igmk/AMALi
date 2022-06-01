

Date.from.matlab<- function(val){as.POSIXct((val - 719529)*86400, origin = "1970-01-01", tz = "UTC")}

#read nc file from matlab processing
readAmalidataStep1<- function(Pfad){
  library(ncdf4)
  lidar.nc <- nc_open(Pfad)
  
  names(lidar.nc$var)
  Time <- ncvar_get(lidar.nc, "Time")
  Distance <- ncvar_get(lidar.nc, "Distance")
  LidarRatio <- ncvar_get(lidar.nc, "LidarRatio")
  
  P532p_raw <- ncvar_get(lidar.nc, "P532p_raw")
  P532s_raw<- ncvar_get(lidar.nc, "P532s_raw")
  P355_raw<- ncvar_get(lidar.nc, "P355_raw")
  
  P532p_bg<- ncvar_get(lidar.nc, "P532p_bg")
  P532s_bg<- ncvar_get(lidar.nc, "P532s_bg")
  P355_bg<- ncvar_get(lidar.nc, "P355_bg")
  
  P532p_noise<- ncvar_get(lidar.nc, "P532p_noise")
  P532s_noise<- ncvar_get(lidar.nc, "P532s_noise")
  P355_noise<- ncvar_get(lidar.nc, "P355_noise")
  
  B532p<- ncvar_get(lidar.nc, "B532p")
  BSR532p<- ncvar_get(lidar.nc, "BSR532p")
  A532p<- ncvar_get(lidar.nc, "A532p")
  
  Error_Code<- ncvar_get(lidar.nc, "Error_Code")
  Position_UBC<- ncvar_get(lidar.nc, "Position_UBC")
  LC532p<- ncvar_get(lidar.nc, "LC532p")
  
  d <- list(
  
  TIMEDIST = data.table(
    matlabtime = rep(Time, each=length(Distance)),
    dist = rep(Distance, length((Time))),
    P532p_raw = as.vector(P532p_raw),
    P532s_raw = as.vector(P532s_raw),
    P355_raw = as.vector(P355_raw),
    P532p_noise = as.vector(P532p_noise),
    P532s_noise = as.vector(P532s_noise),
    P355_noise = as.vector(P355_noise),
    B532p_LR1 = as.vector(B532p[,,1]),
    B532p_LR2 = as.vector(B532p[,,2]),
    B532p_LR3 = as.vector(B532p[,,3]),
    BSR532p_LR1 = as.vector(BSR532p[,,1]),
    A532p_LR1 = as.vector(A532p[,,1]),
    A532p_LR2 = as.vector(A532p[,,2]),
    A532p_LR3 = as.vector(A532p[,,2]),
    Position_UBC_LR1 = as.vector(Position_UBC[,,1]),
    LC532p = as.vector(LC532p[,,1])
  ),
  
  TIME = data.table(
    matlabtime = Time,
    P532p_bg = P532p_bg,
    P532s_bg = P532s_bg,
    P355_bg  = P355_bg,
    Error_Code_LR1 = Error_Code[,1],
    Error_Code_LR2 = Error_Code[,2],
    Error_Code_LR3 = Error_Code[,3]
    )
  )
  
  d$TIMEDIST$datetime<- Date.from.matlab(d$TIMEDIST$matlabtime)
  d$TIME$datetime <- Date.from.matlab(d$TIME$matlabtime)
  
  return(d)
  
}


readGPS <- function(GPSpfad){
  header <- read.table(GPSpfad, nrows = 1, header = FALSE, sep ='\t', stringsAsFactors = FALSE)
  header2 <-read.table(GPSpfad, nrows = 1, header=FALSE, sep='\t', stringsAsFactors = FALSE, skip=1)
  GPS <- read.table(GPSpfad, header = F, skip = 4, sep = '\t')
  colnames(GPS) <- c(unlist(header)[1:6],unlist(header2)[7:13])
  rm(header, header2)
  
  GPS$datetime <- as.POSIXct( paste(GPS$YYYY,'-',GPS$MM,'-', GPS$DD,' ', GPS$HH,':',GPS$mm,':',GPS$ss, sep=''),tz ='UTC' )
  GPS<- GPS[7:ncol(GPS)]
  
  GPS$Lon_dezdeg <- floor(GPS$Lon/100)+ ((GPS$Lon/100-floor(GPS$Lon/100))/60*100)
  GPS$Lon_dezdeg[which(GPS$`Lon dir` == "W")] <- GPS$Lon_dezdeg[which(GPS$`Lon dir` == "W")] * -1
  GPS$Lat_dezdeg <- floor(GPS$Lat/100)+ ((GPS$Lat/100-floor(GPS$Lat/100))/60*100)
  
  GPS$datetime_fac  <- as.integer( GPS$datetime)
  GPS$datetime <- NULL 
  return(GPS)
}


readINS <- function(INSpfad){
  header <- read.table(INSpfad, nrows = 1, header = FALSE, sep ='\t', stringsAsFactors = FALSE)
  header2 <-read.table(INSpfad, nrows = 1, header=FALSE, sep='\t', stringsAsFactors = FALSE, skip=1)
  INS <- read.table(INSpfad, header = F, skip = 4)
  colnames(INS) <- c(unlist(header)[1:6],unlist(header2)[7:13])
  rm(header, header2)
  
  INS$datetime <- as.POSIXct( paste(INS$YYYY,'-',INS$MM,'-', INS$DD,' ', INS$HH,':',INS$mm,':',INS$ss, sep=''),tz ='UTC' )
  INS<- INS[7:ncol(INS)]
  
  INS$datetime_fac  <- as.integer( INS$datetime)
  INS$Latitude  <- NULL
  INS$Longitude <- NULL
  INS$datetime  <- NULL
  INS$`Inertial Altitude` <- NULL
  INS$`Ground Speed` <- NULL
  return(INS)
}




readmirac <- function(Pfad){
  radar.nc <- nc_open(Pfad)
  
  planealt <- ncvar_get(radar.nc, "alt_platform")
  radarheight <- ncvar_get(radar.nc, "alt")
  usefullheight <- radarheight<max(planealt) & radarheight>-50
  usefulltime <- planealt > 2000
  
  planealt <- ncvar_get(radar.nc, "alt_platform")[usefulltime]
  radarheight <- ncvar_get(radar.nc, "alt")[usefullheight]
  
  ze <- ncvar_get(radar.nc, "Ze")[usefullheight,usefulltime]
  radarflag <- ncvar_get(radar.nc, "Ze_flag")[usefullheight,usefulltime]
  
  vm_raw<-  ncvar_get(radar.nc, "vm_raw")[usefullheight,usefulltime]
  vm_sensor<-  ncvar_get(radar.nc, "v_sensor_r")[usefulltime]
  v_sigma<-  ncvar_get(radar.nc, "v_sigma")[usefullheight,usefulltime]
  datetime<-  ncvar_get(radar.nc, "time_target")[usefullheight,usefulltime]
  time_diff_target <- ncvar_get(radar.nc, "time_diff_target")[usefullheight,usefulltime]
  measurementtime <- ncvar_get(radar.nc, "time")[usefulltime] #die Zeit ist in ganzen Sekunden abgelegt, aber der Chirp vielleicht nicht unbedingt ganze Sekunden lang
  
  #miliseconds <- ncvar_get(radar.nc, "sampleTms") !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  #wenn er den offset von der zeit nicht mitnimmt: 
  if(measurementtime[1]< 1e5){
    measurementtime <- measurementtime+datetime[1,1]-measurementtime[1]
  }
  
  nc_close(radar.nc)
  
  #v_par_r <- vm_raw - as.matrix(rep(vm_sensor, length(planealt)))
  #v_par_r <- sweep(vm_raw, 2, vm_sensor, FUN = "-") #sensorr geschwindigkeit von der mean doppler abziehen 
  #das ist quatsch, weil es noch nicht entfaltet ist
  gc()
  print('reading done')
  #miliseconds.nc <- nc_open('/atm_meas/polar_5_6/mirac_radar/data/2017/milliseconds/170602milisecond.nc')
  #miliseconds.nc <- nc_open('/atm_meas/polar_5_6/mirac_radar/data/2017/milliseconds/170527_01milisecond.nc')
  #miliseconds <- ncvar_get(miliseconds.nc, "sampleTms")[usefulltime]
  #nc_close(miliseconds.nc)
  
  measurementtime <- measurementtime #+ miliseconds/1000
  datetime <- sweep(time_diff_target*(-1),2,measurementtime, FUN = "+")
  # wie mach ich das in der Schleife ?
  #Datatable mit Daten erstellen
  radar <- data.table(radarheight=rep(radarheight, length(planealt)),
                      #altitude=rep(planealt, each=length(radarheight)),
                      datetime  = as.vector(datetime),
                      measurementtime = rep(measurementtime,each=length(radarheight)),
                      #miliseconds = rep(miliseconds,each=length(radarheight)),
                      ze   = as.vector(ze),
                      v_sigma   = as.vector(v_sigma),
                      vm_raw   = as.vector(vm_raw),
                      vm_sensor = rep(vm_sensor, each=length(radarheight)),
                      #v_par_r   = as.vector(v_par_r),
                      ze_flag = as.vector(radarflag)
                       ) #asvector geht spaltenweise,
  radar[is.nan(ze), ze:= -9999]
  
  radar$datetime <- as.POSIXct(radar$datetime, origin = '1970-01-01')
  
  print('saved as dt')
  #radar$datetime_fac <- as.factor(radar$datetime)
  #rm(ze, radarflag, radar.nc,miliseconds)
  
  return(radar)
}

