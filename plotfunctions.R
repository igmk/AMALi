library(ggplot2)
library(ggmap)
library(data.table)


colourchoice <- function(farbe, p) {
  #passenden colourscale auswaehlen
  if (farbe == "log(bsr532)"|
      farbe == "log(BSR532p_LR2)") {
    p <- p + fill$bsr_VIS
  }
  else if (farbe == "beta532s/beta532") {
    p <- p + fill$depol
  }
  else if (farbe == "beta355/beta532"|
           farbe == "(P355_raw)/(P532p_raw)*5") {
    #plotdatatime$beta355 <- NA
    p <- p + fill$CR#+
    #geom_line(data=plotdatatime, aes(y=cloudtop1), fill="transparent")
  }
  else if (farbe == "lidarratio532") {
    p <- p + fill$LR
  }
  else if (farbe == "log(beta532)" |
           farbe == "log(B532p_LR1)" |
           farbe == "log(B532p_LR2)" |
           farbe == "log(B532p_LR3)" |
           farbe == "log(attbeta532)"|
           farbe == "log(average_attbeta532)"|
           farbe == "log(attbeta532_CALIOP)"|
           farbe == "10*log10(attbeta532)"|
           farbe == "10*log10(attbeta532_CALIOP)") {
    p <- p + fill$beta_VIS
  }
  else if (farbe == "p532s_bc/p532_bc"|
           farbe == "voldepol") {
    p <- p + fill$volumedepol2
  }
  else if (farbe == "log(p532_bc)" |
           farbe == "log(p532_bc*dist^2)") {
    p <- p + fill$signal_VIS
  }
  else if (farbe == "10*log10(alpha532)") {
    p <- p + fill$alpha_VIS2
  }
  else if (farbe == "log(bsr355)") {
    p <- p + fill$bsr_UV
  }
  else if (farbe == "10*log10(ze)"|
           farbe == "10*log10(average_ze)"|
           farbe == "average_10log10ze"|
           farbe == "ze_CLOUDSAT"|
           farbe == "ze_CS") {
    p <- p + fill$radar +
      geom_tile(data=plotradar[10*log10(ze)>20], fill='violetred4', aes(y=radarheight), height=5.1)
      # geom_tile(
      #   data = plotdata[!duplicated(plotdata, by = c("datetime_fac", "height_rounded")) &
      #                     log(bsr532) > 3 &
      #                     p532snr > 3 & eval(parse(text = cond))],
      #   aes(alpha = log(bsr532)),
      #   width = 0.006,
      #   fill = 'black'
      # ) +
      # scale_alpha_continuous(
      #   expression(paste('log(', BSR [532], ')')),
      #   range = c(0, 0.9),
      #   expand = c(0, 0),
      #   limits = c(3, 10)
      # )
    #print("addBSR")
  }
  else if (farbe == "attbeta532_CALIOP - attbeta532") {
    p <- p + fill$dif_attbeta
  }  
  else if (farbe == "ze_flag") {
    p <- p + fill$LR
  }
  else if (farbe == "v_sigma") {
    p <- p + fill$spectralwidth
  }
  else if (farbe == "vm_raw") {
    p <- p + scale_fill_gradientn(expression("vm_raw"),
                                  colours =c(rainbow(12)[3:12],'darkred'),
                                  #limits = c(20,50),
                                  na.value="transparent")
  }
  else if (farbe == "vm_dev") {
    p <- p + fill$dopplerabw
  }
  else if (farbe == "vm_dif" | farbe == "vm_smoothdif") {
    p <- p + fill$dopplerdz
  }
  
  else if (farbe == "ice") {
    p <- p + scale_fill_manual(values = c('navy', 'cyan'), na.value = "grey")
  }
  else {
    print("kenne farbe nicht")
    break 
  }
}

nameplot <- function(farbe) {
  if (farbe == "beta532s/beta532") {
    farbe2 <- "beta532s_beta532"
  }
  else if (farbe == "beta355/beta532"|
           farbe =="(P355_raw)/(P532p_raw)*5") {
    farbe2 <- "beta355_beta532"
  }
  else if (farbe == "p532s_bc/p532_bc") {
    farbe2 <- "p532s_bc_p532_bc"
  }
  else if (farbe == "log(p532_bc*dist^2)") {
    farbe2 <- "p532_bc_dist2"
  }
  else{
    farbe2 <- farbe
  }
  return(farbe2)
}

lonplot <- function(farbe_input=farbe, cond_input=cond, data_input=plotdata, cloudtop = F, hmax=3200) {
  p <-
    ggplot(data_input[eval(parse(text = cond_input))], aes(Lon_dezdeg, height)) + tl +
    geom_line(data=plotdatatime, colour= 'grey70', aes(y=Altitude))+
    scale_y_continuous('Altitude [m]',
                       expand = c(0, 0),
                       limits = c(-50, hmax)) +
    scale_x_continuous("Longitude", 
                       expand = c(0, 0)) +
    #coord_fixed(ratio = 2e-4)+
    theme(axis.text.x = element_text(angle = 90))+
    theme(legend.position = 'bottom')

  p <- colourchoice(farbe_input, p)
  farbe2 <- nameplot(farbe_input)
  
  if ("p532_bc" %in% colnames(data_input) | 'P532p_raw' %in% colnames(data_input)) {
    p <- p + geom_tile(aes(width = quantile(diff(plotdatatime$Lon_dezdeg),0.98),
                           fill = eval(parse(text = farbe_input)),
                           height = 7.6 #abs(diff(unique(height)))*1.1
                           ))  #so ist der balken sso breiter als 90% der Werte, was macht die 1.5 da?
  }else{
    if("width_lon"%in% colnames(data_input)){
      p <- p + geom_rect(aes( x = NULL, y = NULL,
                              xmin=Lon_dezdeg,
                             xmax=Lon_dezdeg+width_lon,
                             ymin=radarheight,
                             ymax=radarheight+5.1,
                             #height=5.1,
                             fill = eval(parse(text = farbe_input))
                             
      )) 
    }else{
    p <- p + geom_tile(aes(width = quantile(diff(plotdatatime$Lon_dezdeg),0.98)*2,
                           fill = eval(parse(text = farbe_input)),
                           height = 5.1,
                           y= radarheight
                           )) 
  }}
  
  if ("section" %in% colnames(data_input)) {
    p <- p + facet_grid(rows = vars(section))
  }
  
  # if (exists('timestamps')) {
  #   p <- p + geom_text(data = timestamps,
  #                      aes(y=2000, label=substr( datetime_fac, 12, 19)), angle =90, colour="grey20")
  # }
  
  if (cloudtop==T) {
    p <- p + geom_point(data=plotdata[liquidtoptop ==T], aes(y=height), colour='grey5', size=0.005)+
      geom_point(data=plotdata[liquidtoptop ==F], aes(y=height), colour='grey80', size=0.005)
  }
 
  ggsave(
    paste(strftime(data_input[eval(parse(text = cond_input))]$datetime[1], format = '%Y%m%d_%H%M'),
          '_',farbe2,'_lon.png',sep = ""),
    width = 1.7 * klein, #poster 1.7
    height = 1 * klein,
    unit = "cm" ,
    limitsize = FALSE,
    bg = "transparent"
  )
  print(paste(strftime(data_input[eval(parse(text = cond_input))]$datetime[1], format = '%Y%m%d_%H%M'),
              '_',farbe2,'_lon.png',sep = ""))
  return(p)
}

latplot <- function(farbe, cond) {
  p <-
    ggplot(plotdata[eval(parse(text = cond))], aes(Lat_dezdeg, height_rounded, 
                                                   fill = eval(parse(text = farbe)))) + 
    tl +
    geom_tile(width = 0.001) +
    scale_y_continuous('Altitude [m]',
                       expand = c(0, 0),
                       limits = c(-50, 3500)) +
    scale_x_continuous("Lattitude", expand = c(0, 0)) +
    theme(legend.position = 'bottom')
  if ("section" %in% colnames(plotdata)) {
    p <- p + facet_grid(rows = vars(section))
  }
  
  p <- colourchoice(farbe, p)
  farbe2 <- nameplot(farbe)
  
  ggsave(
    paste(
      strftime(plotdata[eval(parse(text = cond))]$datetime[1], format = '%Y%m%d_%H%M'),
      '_',
      farbe2,
      '_lat.png',
      sep = ""
    ) ,
    width = 2 * klein,
    height = 2 * klein,
    unit = "cm" ,
    limitsize = FALSE
  )
}

timeplot <- function(farbe_input=farbe, cond_input=cond, data_input=plotdata, cloudtop = F) {
  p <-
    ggplot(data_input[eval(parse(text = cond_input))] ) + 
    tl +
    scale_y_continuous('Altitude [m]',
                       expand = c(0, 0),
                       limits = c(-50, 3500)) +
    scale_x_datetime("Time [UTC]", 
                     expand = c(0, 0)) +
    #facet_grid(rows=vars(section))+
    theme(legend.position = 'bottom')
  
  p <- colourchoice(farbe_input, p)
  
  if (farbe_input == "10*log10(ze)" |
      farbe_input == "ze_flag"|
      farbe_input == "vm_dev") {
    p <-
      p + geom_tile(width = 2,
                    na.rm = T,
                    aes(datetime, radarheight, 
                            fill = eval(parse(text = farbe_input))))#+ + miliseconds / 1000 - 8
    # geom_tile(data= data_input[!duplicated(data_input, by=c("datetime_fac", "height_rounded"))&log(bsr532)>3 & p532snr>3 & eval(parse(text=cond_input))],
    #           aes(alpha=log(bsr532)),width=1, fill='black')
  }else {
    p <- p + geom_tile(width = 1,aes(datetime, height_rounded, 
                                     fill = eval(parse(text = farbe_input))))
  }
  
  farbe2 <- nameplot(farbe_input)
  ggsave(
    paste(
      strftime(plotdata[eval(parse(text = cond_input))]$datetime[1], format = '%Y%m%d_%H%M'),
      '_',
      farbe2,
      '_time.png',
      sep = ""
    ) ,
    width = 3 * klein,
    height = 1 * klein,
    unit = "cm" ,
    limitsize = FALSE,
    bg = "transparent"
  )
  return(p)
}

altitudedensity <- function(farbe=farbe, cond=cond) {
  p <-
    ggplot(plotdata[eval(parse(text = cond))], aes(eval(parse(text = farbe)), height)) + tl +
    geom_bin2d(bins = 80) +
    fill$density1e4 +
    scale_x_continuous(farbe) +
    scale_y_continuous("Altitude [m]", limits = c(0, 3000)) +
    geom_hline(yintercept = 75,
               lintype = "dashed",
               alpha = 0.8)
  if (farbe == "p532s_bc/p532_bc") {
    p <- p + scale_x_continuous(farbe, limits = c(0, 5))
  }
  
  farbe2 <- nameplot(farbe)
  ggsave(
    paste(
      strftime(plotdata[eval(parse(text = cond))]$datetime[1], format = '%Y%m%d_%H%M'),
      '_',
      farbe2,
      '_heigth_density.png',
      sep = ""
    ) ,
    width = klein,
    height = klein,
    unit = "cm" ,
    limitsize = FALSE,
    bg = "transparent"
  )
  
}

#map+Overview
overviewmap<- function(cond_var =cond) {
  library(mapproj)
  library(mapdata)
  #install.packages('rworldxtra')
  library(rworldxtra)
  #world <- data(countries)
  world <- map_data("world")
  spitzbergen <- subset(world, subregion == "Svalbard")
  rm(world)
  #map of cuurent plotdata
  karte <- ggplot(data = plotdatatime[datetime_fac %in% plotdata[eval(parse(text = cond_var))]$datetime_fac], 
                  aes(x = Lon_dezdeg, y = Lat_dezdeg)) + tl + #geom_osm()+
    #geom_path(data=spitzbergen, aes(x=long, y=lat, group=group))+
    geom_polygon(data = spitzbergen,
                 aes(x = long, y = lat, group = group),
                 fill = 'cornsilk') +
    geom_point(size = 0.5, aes(colour = datetime)) +
    #geom_point(data=plotdatatime2,aes(colour=datetime_fac[1]), size=0.5)+
    coord_map("ortho") +
    scale_color_gradientn(trans = "time",colours = rainbow(10))+
    labs(x = "", y = "",  colour = '') +
    theme(
      panel.background = element_rect(
        fill = "azure2",
        colour = "azure2",
        size = 2,
        linetype = "solid"
      ),
      legend.position = 'right'
    )
  
  if (exists("DSplot")){
    karte <- karte+
      geom_point(data = DSplot, aes(x=GPS.Longitude..deg., y=GPS.Latitude..deg.), colour="orange")
  }
  
  if ("section" %in% colnames(plotdatatime)){
    karte <- karte+
      facet_wrap(~section)
  }
  ggsave(
    paste(
      strftime(plotdata[eval(parse(text = cond_var))]$datetime[1], format = '%Y%m%d_%H%M'),
      '_map.png',
      sep = ""
    ) ,
    width = 2 * klein,
    height = 2 * klein,
    unit = "cm" ,
    limitsize = FALSE,
    bg = "transparent"
  )
}


# scale layout etc================================================================

## 
#### plotting libraries

library(ggplot2)
library(RColorBrewer)
library(viridis)
#library(egg)
#library(gridExtra) #to merge severalplots into one
#library(metR) #to get multiple colourscales


## plot layout -----

#setwd("/atm_meas/awipev/lidartmp/lidar4/bkulla/Amali/R/Quicklooks")


#general Layout
tl <- theme_light(base_size = 16)+
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        # panel.grid.major = element_blank(), # get rid of major grid for Elena
        # panel.grid.minor = element_blank(), # get rid of minor grid
        plot.caption=element_text(margin=margin(t=15),
                                  face="italic", size=8, colour='grey'),#bemerkung/caption nach unten rechts schieben
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", linetype = 0) # get rid of legend panel bg
        
  )

#cloudcols <- c("no visible cloud" = "cyan", "high clouds" = "gray60", "low clouds" = "black")
klein=15
normal=30

####colourscales ----

jet.colors <-
  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))


fill <- list()
fill$LR<- scale_fill_gradientn(expression(paste(LR [532])),
                               colours =c(rainbow(12)[1:9]),
                               na.value="transparent",
                               #values = c(0,       0.33,     0.34,    0.36,     0.39,   0.42,          0.46,    0.5,    0.6,      0.75,     0.9),
                               limits = c(3,28))
fill$beta_VIS<-scale_fill_gradientn(expression(paste('log(',beta [532],')')),
                                    colours =c('white','grey95','grey80','grey65','yellow','springgreen3',brewer.pal(6,'YlGnBu')[4:6],'darkorchid4', 'black'),
                                    na.value="transparent",
                                    values = c(0,       0.33,     0.34,    0.36,     0.39,   0.42,          0.46,    0.5,    0.6,      0.75,     0.9),
                                    limits = c(-20,-2))

fill$beta_VIS_calipso<-scale_fill_gradientn(expression(paste('log(',beta [532],')')),
                                            colours =c('white','grey95','grey80','grey65','yellow','springgreen3',brewer.pal(6,'YlGnBu')[4:6],'darkorchid4', 'black'),
                                            na.value="transparent",
                                            values = c(0,       0.33,     0.34,    0.36,     0.39,   0.42,          0.46,    0.5,    0.6,      0.75,     0.9),
                                            limits = c(-16,5))

fill$bsr_VIS<-scale_fill_gradientn(expression(paste('log(',BSR [532],')')),
                                   colours =c('white','grey95','grey80','grey65','yellow','springgreen3',brewer.pal(6,'YlGnBu')[4:6],'darkorchid4', 'black'),
                                   na.value="transparent",
                                   values = c(0,       0.01,    0.04,    0.07,     0.1,    0.2,          0.3,    0.5,    0.6,      0.75,     0.9),
                                   limits = c(0,12))

fill$alpha_VIS<-scale_fill_gradientn(expression(paste('10log(',alpha [532],')')),
                                     colours =c('white','grey95','grey80','grey65','yellow','springgreen3',brewer.pal(6,'YlGnBu')[4:6],'darkorchid4', 'black'),
                                     na.value="transparent",
                                     values = c(0,       0.1,     0.2,    0.3,     0.39,   0.42,          0.46,    0.5,    0.6,      0.75,     0.9),
                                     #limits = c(-15,1)) log(alpha)
                                     limits = c(-65,1))#10log10(alpha)

fill$alpha_VIS2<-scale_fill_gradientn(expression(paste('10log(',alpha [532],')')),
                                      colours =c('white','grey98','grey84','grey70','honeydew3','springgreen3',brewer.pal(6,'YlGnBu')[4:6],'darkorchid4', 'black'),
                                      na.value="transparent",
                                      values = c(0,      0.25 ,   0.35,     0.42,    0.5,     0.55,   0.6,          0.65,    0.7,    0.8,      0.9),
                                      #limits = c(-15,1)) log(alpha)
                                      limits = c(-65,1))#10log10(alpha)

fill$signal_VIS<-scale_fill_gradientn(expression(paste('log(',sigma [532],')')),
                                      colours =c('white','grey95','grey80','grey65','yellow','springgreen3',brewer.pal(6,'YlGnBu')[4:6],'darkorchid4', 'black'),
                                      na.value="transparent",
                                      values = c(0,       0.2,     0.25,    0.28,     0.35,   0.42,          0.46,    0.5,    0.6,      0.75,     0.9),
                                      limits = c(17,28))

fill$beta_UV <- scale_fill_gradientn(expression(paste(beta [355])),
                                     colours = c('grey95', 'grey90','grey80','pink','plum','purple', 'navy', 'black'),
                                     values =  c(0,         0.15,   0.2,     0.25,    0.3,    0.4,       0.5,        0.75),
                                     na.value="transparent",
                                     limits=c(-20,0))

fill$bsr_UV <-         scale_fill_gradientn(expression(paste(BSR [355])),
                                            colours = c('grey95', 'grey90','grey80','pink','plum','purple', 'navy', 'black'),
                                            values = c(0,       0.01,    0.04,    0.07,     0.1,    0.2,          0.3,    0.5,    0.6,      0.75,     0.9),
                                            na.value="transparent",
                                            limits=c(0,10))


fill$depol <-          scale_fill_gradientn('depolarisation\nratio',
                                            colours = c('grey10','navy','steelblue', 'cyan', 'grey90'), 
                                            values = c(0,0.001,0.1,0.15,0.2,1),
                                            na.value="transparent",
                                            limits=c(0,0.2))
fill$volumedepol <-    scale_fill_gradientn('volume\ndepolarisation\nratio',
                                            colours = c('grey10','navy','steelblue', 'cyan', 'grey90'), 
                                            values = c(0,0.01,0.2,0.3,0.4,1),
                                            na.value="transparent",
                                            limits=c(0,3))
fill$volumedepol2 <-    scale_fill_gradientn('volume\ndepolarisation\nratio',
                                             colours = c('grey10','navy','steelblue', 'cyan', 'grey90'), 
                                             values = c(0,0.01,0.05,0.1,0.2,1),
                                             na.value="transparent",
                                             limits=c(0,1))

fill$CR <-             scale_fill_gradientn('colour\nratio',
                                            colours = c('grey90','yellow','orange','tomato','darkred', 'black'),
                                            na.value="transparent", 
                                            limits=c(0,7))

fill$radar <-          scale_fill_gradientn('radar\nreflectivity\ndBZ',
                                            colours = c(jet.colors(10),'darkred', 'violetred4'),
                                            na.value="transparent",
                                            limits= c(-40,20))
fill$radarcloudsat <-  scale_fill_gradientn('radar\nreflectivity\ndBZ',
                                            colours = c(jet.colors(10),'darkred', 'violetred4'),
                                            na.value="transparent",
                                            limits= c(10,38))

fill$radar2 <-         scale_fill_gradientn('radar\nreflectivity\ndBZ',
                                            colours = rev(c(inferno(12))),
                                            na.value="transparent",
                                            limits= c(-40,20))
fill$radarcolour <-  scale_colour_gradientn('radar\nreflectivity\ndBZ',
                                            colours = c(rev(rainbow(12)[1:10]),'darkred'),
                                            na.value="transparent",
                                            limits= c(-40,20))

fill$spectralwidth <-   scale_fill_gradientn(expression("spectral\nwidth"),
                                             colours =c('white','grey85', rev(inferno(12))[2:12], 'black', 'black'),
                                             values = c(0,0.22,seq(0.32,0.68,length.out = 12),1),
                                             limits = c(0,1.1),
                                             na.value="transparent")

fill$dopplerabw <- scale_fill_gradientn(expression(paste(v ['d'] ,' - ',tilde(v)[column])),
                                        colours =c('royalblue', 'grey98','darkred'),
                                        limits = c(-3,3),
                                        na.value="transparent")

fill$dopplerdzabs<-   scale_fill_gradientn(expression(paste('| v '['n'] ,' - ',v['n+1'],'|')),
                                           colours=c("grey95",rainbow(16)[9:16]),
                                           #colours =c('white', 'grey50','black'),
                                           values = c(0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
                                           limits = c(0,0.075),
                                           na.value="transparent")
fill$dopplerdz <- scale_fill_gradientn(expression(paste(v ['n'] ,' - ',v['n+1'])),
                                       colours =c('steelblue', 'cornsilk1','firebrick4'),
                                       limits = c(-0.075,0.075),
                                       na.value="transparent")


fill$dif_attbeta <- scale_fill_gradientn(expression(paste(beta[532],"'II","(CALIOP) - ",beta[532],"'II","(AMALi)")),
                                         colours =c('navy','royalblue', 'grey98','darkred', 'maroon'),
                                         limits = c(-1.5e-3,1.5e-3),
                                         na.value="transparent")

fill$dif_attbeta2 <- scale_fill_gradientn(expression(paste(beta[532],"'II","(CALIOP) - ",beta[532],"'II","(AMALi)")),
                                         colours =c('navy','royalblue', 'grey98','darkred', 'maroon'),
                                         limits = c(-30,30),
                                         na.value="transparent")

fill$dif_ze <- scale_fill_gradientn("Ze[dB](CLOUDSAT CPR) - Ze[dB](MiRAC)",
                                         colours =c('navy','royalblue', 'grey98','darkred','maroon'),
                                         limits = c(-80,80),
                                         na.value="grey")


fill$density1e4 <- scale_fill_gradientn("",colours=c("grey95",rainbow(16)[9:16]), na.value="red", 
                                        limits=c(0,10000))

fill$density1e3 <- scale_fill_gradientn("",colours=c("grey95",rainbow(16)[9:16]), na.value="red", 
                                        limits=c(0,1000))

fill$density2e2 <- scale_fill_gradientn("",colours=c("grey95",rainbow(16)[9:16]), na.value="red", 
                                        limits=c(0,200))

fill$density1e2 <- scale_fill_gradientn("",colours=c("grey95",rainbow(16)[9:16]), na.value="red", 
                                        limits=c(0,100))

fill$density1e_3 <- scale_fill_gradientn("",colours=c("grey95",rainbow(16)[9:16]), na.value="red", 
                                        limits=c(0,1e-3))

set_panel_heights <- function(g, heights){
  g$heights <- grid:::unit.list(g$heights) # hack until R 3.3 comes out
  id_panels <- unique(g$layout[g$layout$name=="panel", "t"])
  g$heights[id_panels] <- heights
  g
}


#mogelfunktion weil meine R-version alt ist:
isFALSE <- function (x) {  is.logical(x) && length(x) == 1L && !is.na(x) && !x}

