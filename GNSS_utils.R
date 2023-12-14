
# GNSS_utils written for PRJ3280


library(tidyverse)
library(sf)
library(readxl)
library(mapview)
library(plotly)
library(ggforce)

# to include this set of utilities in other code:
#source("C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/Code/R/PRJ3820-DOC-GNSS/GNSS_utils.R")


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### crs_convert
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

crs_convert = function(csv_ffn = NULL, crs_src, crs_dst, plot_dst=T){
  
  
  # 
  
  # EPSG:9988     ITRF2020      Cartesian/Geocentric     Ellipsoid: GRS 1980
  # EPSG:9990     ITRF2020      Geodetic  (degrees)      Ellipsoid: GRS 1980
  
  # EPSG:7789     ITRF2014      Cartesian/Geocentric     Ellipsoid: GRS 1980
  # EPSG:7912     ITRF2014      Geodetic                 Ellipsoid: GRS 1980
  
  # EPSG:32760    UTM Zone 60S  Geodetic                 Ellipsoid: GRS 1980
  
  
  csv_ffn = NULL
  crs_src = 7912
  crs_dst = 2193
  
  if (is.null(csv_ffn)){csv_ffn = "C:/temp/coords2convert.csv"}
  
  src_data = read.csv("C:/temp/coords2convert.csv")
  
  if (is.na(src_data$LAT_SRC[1] )){
    
    #this means that the input data is in DMS format
    src_data = src_data %>% 
      mutate(
        LAT_SRC = sign(LAT_DEG) *(abs(LAT_DEG) + LAT_MIN/60 + LAT_SEC/3600),
        LON_SRC = sign(LON_DEG) *(abs(LON_DEG) + LON_MIN/60 + LON_SEC/3600)
      )
  }
  
  
  crs_conv2 = function(DMS_LATLON = NULL, DD_LATLON = NULL, src_crs, dst_crs){
    
    
    DMS_LATLON = c(-40,20,10.51917, 175,48,45.27946)
    DMS_LATLON = c(-40,20,10.5481, 175,48,45.4390)
    
    src_crs = 9990
    dst_crs = 32760
    
    
    if (!is.null(DMS_LATLONG)){
      
      DD_LAT = sign(DMS_LATLON[1]) * (abs(DMS_LATLON[1]) + DMS_LATLON[2]/60 + DMS_LATLON[3]/3600)
      DD_LON = sign(DMS_LATLON[4]) * (abs(DMS_LATLON[4]) + DMS_LATLON[5]/60 + DMS_LATLON[6]/3600)
      
      DD_LATLON = c(DD_LAT, DD_LON)
      
    }
    
    DF = data.frame(X = DD_LATLON[2], Y = DD_LATLON[1])
    SF_src = st_as_sf(DF, coords = c("X","Y"), crs = "EPSG:7912")
    SF_dst = st_transform(SF_src, crs = dst_crs)
    SF_2193 = st_transform(SF_dst, crs = "EPSG:2193")
    
    sprintf("%.12f",st_coordinates(SF_dst))
    sprintf("%.12f",st_coordinates(SF_2193))
  }
  
  
  src_sf = st_as_sf(src_data, coords = c("LON_SRC", "LAT_SRC"), crs = crs_src) 
  
  
  
  # src_sf = st_as_sf(src, coords = c(-4855980.940, 355080.488), crs = 7789) 
  mapview(src_sf)
  
  int_sf = st_transform(src_sf, "epsg:4326") 
  dst_sf = st_transform(int_sf, "epsg:2193") 
  
  # 
  
  dst_coords = st_coordinates(dst_sf) %>% as.data.frame()
  
  dst = src_data %>% 
    mutate(
      SRC_CRS = crs_src,
      LAT_DST = dst_coords$Y,
      LON_DST = dst_coords$X,
      DST_CRS = crs_dst
    )
  
  if (plot_dst){mapview(dst_sf)}
  
  sprintf("%.3f",dst$LAT_DST) 
  sprintf("%.3f",dst$LON_DST) 
  
  results_out_ffn = paste0(tools::file_path_sans_ext(csv_ffn),"_converted.csv")
  write_csv(dst, results_out_ffn)
  
  
}


draw_circle=function(g, RADIUS, CX, CY, NPOINTS = 100, COL = "red", LTYPE = "dashed"){
  
  CIRC = circleFun(center = c(CX, CY), diameter = RADIUS*2, npoints = NPOINTS)
  
  g = g + geom_path(data = CIRC, aes(x,y), color = COL, linetype = LTYPE)
  g
  
}

circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  circ = data.frame(x = xx, y = yy) 
  # circ_pos = circ %>% filter(y>=center[2]) %>% dplyr::arrange(x)
  # circ_neg = circ %>% filter(y<=center[2]) %>% dplyr::arrange(x, "desc")
  # circ_srt = bind_rows(circ_pos, circ_neg)
  return(circ)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### rotate coordinate system to find angle of rotation that
# minimizes correlation between EASTING AND NORTHING Coordinates
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+


find_semi_major = function(D, known_loc){
  
  # D = read_csv("c:/temp/CRSD.csv") %>% 
  #   mutate(
  #     DX = X-known_loc$X,
  #     DY = Y-known_loc$Y,
  #   )
  # known_loc = data.frame(PEG = 6, X = 1838913.032, Y = 5531118.301)
  
  D = D %>% 
    mutate(
      DX = X - known_loc$X,
      DY = Y - known_loc$Y
    )
  
  interv_deg = 1 # degrees
  interv = interv_deg*pi/180
  angles = seq(-pi/2, pi/2, interv)
  
  for (i in 1:length(angles)){
    
    itheta = angles[4]
    itheta = angles[i]  
    Drot = D %>% 
      mutate(
        
        DXROT = DX * cos(itheta) + DY * sin(itheta),
        DYROT = -1 * DX * sin(itheta) + DY * cos(itheta)
        
      )
    
    COV_orig = cov(Drot$DX, Drot$DY)
    COV = cov(Drot$DXROT, Drot$DYROT)
    
    DROW = data.frame(i, itheta, COV)
    
    if (i==1){
      DOUT = DROW
    } else {
      DOUT = DOUT %>% bind_rows(DROW)
    }
  }
  
  print(DOUT)
  
  
  min_cov_theta = DOUT$itheta[which.min(abs(DOUT$COV))]
  min_cov_theta_deg = min_cov_theta*180/pi
  
  g = ggplot(DOUT) + geom_point(aes(itheta, COV))
  g = g + geom_vline(aes(xintercept = min_cov_theta))
  g
  
  Drot_final = D %>% 
    mutate(
      DXROT = DX * cos(min_cov_theta) + DY * sin(min_cov_theta),
      DYROT = -1 * DX * sin(min_cov_theta) + DY * cos(min_cov_theta)
      
    )
  
  
  
  g = ggplot(Drot_final) + geom_point(aes(DX, DY))
  g = g + geom_point(aes(DXROT, DYROT), color = "red")
  g = g + lims(x = c(-.75,.75), y = c(-.75,.75))
  g
  
  
  
  DX_SIGMA_95 = 2 * sd(Drot_final$DXROT)
  DY_SIGMA_95 = 2 * sd(Drot_final$DYROT)
  
  DX_SIGMA_95_ROT = DX_SIGMA_95 * cos(min_cov_theta) + 0 * sin(min_cov_theta)
  DY_SIGMA_95_ROT = 0 * cos(min_cov_theta) + DY_SIGMA_95 * sin(min_cov_theta)
  
  if (DX_SIGMA_95 > DY_SIGMA_95){
    Semimajor_azimuth = min_cov_theta
    Semiminor_azimuth = min_cov_theta - pi/2
    Semimajor_2SD = DX_SIGMA_95
    Semiminor_2SD = DY_SIGMA_95
  }else{
    Semimajor_azimuth = min_cov_theta - pi/2
    Semiminor_azimuth = min_cov_theta 
    Semimajor_2SD = DY_SIGMA_95
    Semiminor_2SD = DX_SIGMA_95
  }
  
  DOUT = data.frame(
    Semimajor_azimuth = Semimajor_azimuth,
    Semiminor_azimuth = Semiminor_azimuth,
    Semimajor_2SD = Semimajor_2SD,
    Semiminor_2SD = Semiminor_2SD)
  
  
  # NOW ROTATE THE sigma values
  
  # ggplot(DOUT) + geom_ellipse(aes(x0 = 0,y0 = 0, a = Semimajor_2SD, b = Semiminor_2SD, angle = 90 ))
  
  
  return(DOUT)
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### plot ellipse
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
# g = ggplot object that the ellipse will be plotted to
# XC = x coordinate of centre of ellipse
# XY = y coordinate of centre of ellipse
# A = semi-minor axis of ellipse
# B = semi-major axis of ellipse
# THETA = angle (measured anticlockwise from positive horisontal axis) of ellipse

plot_ellipse = function(g, XC, YC, A, B, THETA){
  
  # 
  # XC = 0
  # YC = 0
  # A = 4
  # B = 2
  # THETA = pi/4
  
  t = seq(0, 2*pi, 0.01) 
  x = XC + A*cos(t)*cos(THETA) - B*sin(t)*sin(THETA)
  y = YC + A*cos(t)*cos(THETA) + B*sin(t)*cos(THETA)
  
  DF = data.frame(t,x,y)
  
  g2 = ggplot(DF) + geom_path(aes(x=x,y=y))
  g2 = g2 + geom_point()
  
  g2 = g2 + geom_point(data = DF, aes(x,y))
  
  g3 = ggplot(DF) + geom_point(aes(x=x,y=y))
  g3 = g3 + geom_path(aes(x,y), color="red")
  return(g)
  
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### calc_acc_and_prec
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+

calc_acc_and_prec = function(GNSS_data, known_loc){
  
  n = nrow(GNSS_data)
  melap = max(GNSS_data$melap)
  
  GNSS_data = GNSS_data %>% 
    mutate(
      X_DEV = (X - known_loc$X),
      Y_DEV = (Y - known_loc$Y),
      
      X_DEV2 = X_DEV^2,
      Y_DEV2 = Y_DEV^2
    )
  
  GNSS_data_smy = GNSS_data %>% 
    summarise(
      mean_X_DEV = mean(X_DEV),
      mean_Y_DEV = mean(Y_DEV),
      
      M_X_tmp = sum(X_DEV2),
      M_Y_tmp = sum(Y_DEV2),
      SIGMA_X = sd(X_DEV),
      SIGMA_Y = sd(Y_DEV)
    ) %>% 
    mutate(
      M_X = sqrt(M_X_tmp/n),
      M_Y = sqrt(M_Y_tmp/n),
      # M_ACC = sqrt(M_X^2 + M_Y^2),
      M_ACC = sqrt(mean_X_DEV^2 + mean_Y_DEV^2),
      M_PREC = sqrt(SIGMA_X^2 + SIGMA_X^2),  #this is DRMS
      DRMS_65 = sqrt(SIGMA_X^2 + SIGMA_X^2),
      DRMS_95 = 2*DRMS_65,
      SIG_Y_X_RATIO = SIGMA_Y/SIGMA_X,
      CEP_50 = ifelse(SIG_Y_X_RATIO > 0.3, 0.62*SIGMA_Y + 0.56*SIGMA_X, DRMS_95/2.4),#accurate when SIGMA
      CEP_95 = CEP_50*2.08,
      N_MEAS = n,
      minutes_elapsed = melap)
  
} 



### ===> Function analyze_POS\\


# calls 

# 
# anl_pos_fig_1
# calc_acc_and_prec

#-----------------------------------------------------------------------#
# read_POS: read an RTKLIB output file
#-----------------------------------------------------------------------#


read_POS = function(POS_ffn,FMT_OPT=NULL){
  
  lines = readLines(POS_ffn,100)
  NSKIP = min(which(!( str_detect(lines, "^\\%"))))-1
  
  if (is.null(FMT_OPT)){
    
    D = read_csv(POS_ffn, skip = NSKIP, col_names =F, show_col_types = FALSE)
    
    
    
    HDRS = c("GPS_Week", "GPS_Time",  "Lat_DD",  "Lon_DD",  "Ht_ell", "Q",
             "nvsv","sdn", "sde", "sdu", "sdne" , "sdeu","sdun", "age" ,"ratio")
    names(D) <- HDRS
    D %>% pull(GPS_Week) %>% unique()
    D = D %>% 
      mutate(
        across(Lat_DD:ratio, as.numeric),
        melap = (GPS_Time - GPS_Time[1])/60,
        UTC_TIME_PX = GPS_tm_2_UTC(GPS_Week , GPS_Time )) 
    
    
  }else{
    
    leapseconds_curr = 18
    
    if (FMT_OPT==2){
      
      D = read_csv(POS_ffn, skip = NSKIP, col_names =F)
      HDRS = c("GPS_DateTime",   "Lat_DD",  "Lon_DD",  "Ht_ell", "Q",
               "nvsv","sdn", "sde", "sdu", "sdne" , "sdeu","sdun", "age" , "ratio")
      names(D) <- HDRS
      D = D %>%
        dplyr::filter(!is.na(GPS_DateTime) | GPS_DateTime != "" | GPS_DateTime != "NA") %>% 
        mutate(
          GPS_Time = as.POSIXct(GPS_DateTime, format = "%Y/%m/%d %H:%M:%S"),
          melap = (as.numeric(GPS_Time) - as.numeric(GPS_Time[1]))/60,
          UTC_TIME_PX = GPS_Time + seconds(leapseconds_curr)
        )
      
      
    } else {
      
      D = read_csv(POS_ffn, skip = NSKIP, col_names =F)
      HDRS = c("GPS_Week", "GPS_Time",  "Lat_DD",  "Lon_DD",  "Ht_ell", "Q",
               "nvsv","sdn", "sde", "sdu", "sdne" , "sdeu","sdun", "age" ,"ratio")
      names(D) <- HDRS
      D %>% pull(GPS_Week) %>% unique()
      D = D %>% 
        mutate(
          across(Lat_DD:ratio, as.numeric),
          melap = (GPS_Time - GPS_Time[1])/60,
          UTC_TIME_PX = GPS_tm_2_UTC(GPS_Week , GPS_Time )) 
      
      
    }
  }
  return(D)
}



#======================================================================================#
# figure 1. the xy plot of the solution
#======================================================================================#

anl_pos_fig_1 = function(D_sf, RCVR, OCC, PEG, PERIOD, METHOD, known_loc, plotdir, plotdim, limax, CEP=NULL){
  
  CRDS = D_sf %>% 
    st_coordinates() %>%  as.data.frame() %>% as_tibble() %>% 
    mutate(melap = D_sf$melap)
  
  CRDS_last5 = CRDS %>% tail(10) %>% summarise(across(everything(), mean))
  
  # DEL = sqrt( (known_loc$X - CRDS$X)^2 + (known_loc$Y - CRDS$Y)^2 )
  
  AccPrec = calc_acc_and_prec(CRDS, known_loc)
  
  if (is.null(CEP)){CEP = c(AccPrec$CEP_50, AccPrec$CEP_95)}
  
  minx = min(CRDS$X)
  maxx = max(CRDS$X)
  miny = min(CRDS$Y)
  maxy = max(CRDS$Y)
  rngx = maxx - minx
  rngy = maxy - miny
  rng_largest = max(c(rngx,rngy))
  # midx = mean(c(minx, maxx))
  # midy = mean(c(miny, maxy))
  midx = known_loc$X
  midy = known_loc$Y
  
  # plotdim = rng_largest * 1.3
  # plotdim = 3
  xlim = c(midx - 0.5*plotdim, midx + 0.5*plotdim)
  ylim = c(midy - 0.5*plotdim, midy + 0.5*plotdim)
  
  
  
  if (limax){
    g  = ggplot(D_sf) + geom_sf(aes(fill = melap, color = melap)) + coord_sf(datum = 2193, xlim = xlim, ylim = ylim) + geom_point(data = known_loc, aes(X,Y), color = "red", size =5) + scale_fill_continuous(type = "viridis") + scale_color_continuous(type = "viridis") 
  }else{
    g  = ggplot(D_sf) + geom_sf(aes(fill = melap, color = melap)) + coord_sf(datum = 2193) + geom_point(data = known_loc, aes(X,Y), color = "red", size =5) + scale_fill_continuous(type = "viridis") + scale_color_continuous(type = "viridis") 
    
  }
  
  if (!is.null(CEP)){
    # 
    # CEP_50 = CEP[1]
    # XNDIVS_50 = 100
    # DIVSZ_50 = 2*CEP_50/XNDIVS_50
    # XSERIES_50 = seq(-1*CEP_50, CEP_50, by=DIVSZ_50)
    # YSERIES_50_POS =     sqrt(CEP_50^2 - ((XSERIES_50 - known_loc$X)^2)) + known_loc$Y
    # YSERIES_50_NEG = -1*(sqrt(CEP_50^2 - ((XSERIES_50 - known_loc$X)^2)) + known_loc$Y)
    # CEP_50_df = data.frame(X = c(XSERIES_50, XSERIES_50), Y  = c(YSERIES_50_POS, YSERIES_50_NEG))
    # 
    # CEP_95 = CEP[2]
    # XNDIVS_95 = 100
    # DIVSZ_95 = 2*CEP_95/XNDIVS_95
    # XSERIES_95 = seq(-1*CEP_95, CEP_95, by=DIVSZ_95)
    # YSERIES_95_POS =     sqrt(CEP_95^2 - ((XSERIES_95 - known_loc$X)^2)) + known_loc$Y
    # YSERIES_95_NEG = -1*(sqrt(CEP_95^2 - ((XSERIES_95 - known_loc$X)^2)) + known_loc$Y)
    # CEP_95_df = data.frame(X = c(XSERIES_95, XSERIES_95), Y  = c(YSERIES_95_POS, YSERIES_95_NEG))
    # 
    CEP_50 = CEP[1]
    g = draw_circle(g, CEP_50, CX = known_loc$X, CY = known_loc$Y, NPOINTS = 100)
    CEP_95 = CEP[2]
    g = draw_circle(g, CEP_95, CX = known_loc$X, CY = known_loc$Y, NPOINTS = 100)
    
    
  }
  
  
  
  g = g + scale_x_continuous(breaks = seq(xlim[1], xlim[2]))
  g = g + scale_y_continuous(breaks = seq(ylim[1], ylim[2]))
  
  g = g + labs(x= "EASTING (m)", y = "NORTHING (m)", color = "Occupation\nTime\n(minutes)")
  g = g + scale_fill_continuous(guide = 'none')
  
  
  g = g + annotate("text", x = xlim[1] + diff(xlim)*.05, y = ylim[2] - diff(ylim)*.05, label = RCVR, hjust=0)
  g
  
  
  
  DESCR = paste("Peg", PEG," - Occupation", OCC, " - ", PERIOD)
  g = g + annotate("text", x = xlim[1] + diff(xlim)*.05, y = ylim[2] - diff(ylim)*.15, label = DESCR, hjust=0)
  
  
  ACC_PREC_METRICS = paste0("M_Acc = ", formatC(AccPrec$M_ACC, digits = 11) , ", M_Prec = ", formatC(AccPrec$M_PREC, digits = 11) )
  g = g + annotate("text", x = xlim[1] + diff(xlim)*.05, y = ylim[2] - diff(ylim)*.25, label = ACC_PREC_METRICS, hjust=0)
  g
  
  g = g + annotate("text", x = xlim[1] + diff(xlim)*.05, y = ylim[2] - diff(ylim)*.55, label = paste0("DRMS(65%) = ", formatC(AccPrec$DRMS_65, digits = 3), "m"), hjust=0)
  g = g + annotate("text", x = xlim[1] + diff(xlim)*.05, y = ylim[2] - diff(ylim)*.60, label = paste0("DRMS(95%) = ", formatC(AccPrec$DRMS_95, digits = 3), "m"), hjust=0)
  g = g + annotate("text", x = xlim[1] + diff(xlim)*.05, y = ylim[2] - diff(ylim)*.65, label = paste0("CEP(50%) = ", formatC(AccPrec$CEP_50, digits = 3), "m"), hjust=0)
  g = g + annotate("text", x = xlim[1] + diff(xlim)*.05, y = ylim[2] - diff(ylim)*.70, label = paste0("CEP(95%) = ", formatC(AccPrec$CEP_95, digits = 3), "m"), hjust=0)
  g
  
  g = g + annotate("text", x = xlim[1] + diff(xlim)*.05, y = ylim[2] - diff(ylim)*.85, label = paste0("Known Peg Location = \nX: ",known_loc$X , ", Y: ", known_loc$Y), hjust=0)
  g
  
  g = g + annotate("text", x = xlim[1] + diff(xlim)*.05, y = ylim[2] - diff(ylim)*.95, label = paste0("Measured Peg Location = \nX: ",formatC(CRDS_last5$X, digits=10) , ", Y: ", formatC(CRDS_last5$Y, digits=10)), hjust=0)
  g
  
  
  g1 = g
  
  # plot_subdn = paste0(plotdir, RCVR, "/")
  # if (!dir.exists(plot_subdn)){dir.create(plot_subdn)}
  
  plot_fname_g1 = paste0(RCVR,"-OCC#", OCC, "-PEG#", PEG,"-Method-",METHOD,"-", PERIOD,".png" )
  plot_ffname_g1 = paste0(plotdir, plot_fname_g1)
  ggsave(plot_ffname_g1, g1, width = 18, height =18, units = 'cm')
  
  return(g1)
}

#======================================================================================#
# figure 2. DISTRIBUTION OF ERROR
#======================================================================================#

anl_pos_quantify_error = function(known_loc, D_sf){
  
  
  
  
  CRDS = D_sf %>% 
    st_coordinates() %>%  as.data.frame() %>% as_tibble() %>% mutate(melap = D$melap)
  
  AccPrec = calc_acc_and_prec(CRDS, known_loc)
  
  
  
  CRDS_last5 = CRDS %>% tail(10) %>% summarise(across(everything(), mean))
  
  # DEL = sqrt( (known_loc$X - CRDS_last5$X)^2 + (known_loc$Y - CRDS_last5$Y)^2 )
  
  
  
  minx = min(CRDS$X)
  maxx = max(CRDS$X)
  miny = min(CRDS$Y)
  maxy = max(CRDS$Y)
  rngx = maxx - minx
  rngy = maxy - miny
  rng_largest = max(c(rngx,rngy))
  # midx = mean(c(minx, maxx))
  # midy = mean(c(miny, maxy))
  midx = known_loc$X
  midy = known_loc$Y
  
  # plotdim = rng_largest * 1.3
  # plotdim = 3
  xlim = c(midx - 0.5*plotdim, midx + 0.5*plotdim)
  ylim = c(midy - 0.5*plotdim, midy + 0.5*plotdim)
  
  errors_df = data.frame(x_meas = CRDS$X, y_meas = CRDS$Y) %>% 
    as_tibble() %>% 
    mutate(x_dev = known_loc$X - x_meas, y_dev = known_loc$Y)
  
  
  
  
}


#======================================================================================#
# figure 3. DISTRIBUTION OF ERROR
#======================================================================================#

# anl_pos_quantify_error = function(known_loc, D_sf){

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### wineq_path
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Check to make sure these three files exist
wineq_path = function(PATH){
  PATHOUT = str_replace(PATH, "\\/mnt\\/c", "C\\:\\")
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### read_rinex
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Check to make sure these three files exist
read_rinex = function(PATH){
  library(asteRisk)
  testGPSnav <- readGPSNavigationRINEX("C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/data/Field-Test-Ballance-20231108/R10-Base/RINEX211/01393111.23n")
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### GPS_tm_2_UTC
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Check to make sure these three files exist
GPS_tm_2_UTC = function(gps_week, gps_sec_of_current_week){
  # gps_week = c(2287,2288)
  # gps_sec_of_current_week = c(250846,250846)
  
  gps_week_in_sec = 7*86400*(gps_week)
  sec_since_rollover = gps_week_in_sec+gps_sec_of_current_week
  rollover_date = as.POSIXct("1980-1-6",  origin = "1970-01-01 00:00:00", tz = "UTC")
  rollover_date_nmc = as.numeric(rollover_date)
  leapseconds = 18
  current_sec = rollover_date_nmc + sec_since_rollover + leapseconds
  current_datetime_UTC = as.POSIXct(current_sec,  origin = "1970-01-01 00:00:00", tz = "UTC")
}




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### proc_rtk
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# InputFilesOption:
#   0 = Autonomous
#   1 = PPK (Base Station Obs and Nav)
#   2 = PPK (Base station Obs + Nav) + Precise EPH + Precise Clock


proc_rtk = function(Rover, Occ, Method, Period, InputFilesOption = 2, plotdim = 3, replot_only = F){
  
  
  
  testmode = F
  
  if (testmode){
    
    Rover = "JVD"
    Occ = 2
    Method = "RTK-R10"
    Period = "60MIN"
    # plotdim=3
    # Method = "RTK-CORS"
  }
  
  
  
  datadir_w10 = "C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/data/"
  datadir_wsl = "/mnt/c/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/data/"
  
  fielddatadir_w10 = paste0(datadir_w10, "Field-Test-Ballance-20231108/")
  fielddatadir_wsl = paste0(datadir_wsl, "Field-Test-Ballance-20231108/")
  
  solution_dn_w10 = paste0(fielddatadir_w10, "SOLUTIONS/")
  solution_dn_wsl = paste0(fielddatadir_wsl, "SOLUTIONS/")
  
  #read the metadata table
  
  
  field_metadata_binary_ffn = paste0(datadir_w10, "GNSS_DATA_COLLATED.RDS")
  reread_xls=F
  
  if (reread_xls){
    field_metadata_xls_ffn = paste0(datadir_w10, "GNSS_DATA_COLLATED.xlsx")
    
    field_metadata = read_excel(field_metadata_xls_ffn, sheet = "BAL-FIELD-TEST-231108", range = "A6:AR24") %>% 
      filter(!is.na(RECEIVER))
    
    saveRDS(field_metadata, field_metadata_binary_ffn)
    
    
  } else {
    
    field_metadata = readRDS(field_metadata_binary_ffn)
    
    
  }
  
  # get the relevant row of metadata
  current_metadata_row = field_metadata %>% filter(RECEIVER  == Rover & OCC_NUM == Occ)
  pegid = current_metadata_row$PEG_ID
  pegnum = as.numeric(substr(pegid,6,6))
  
  rover_fn = current_metadata_row$RAW_RINEX_FNAME
  rover_sdn = paste0(Rover, "/BAL-", Rover, "-OCC", Occ,"-P", pegnum, "-RINEX/" )
  rover_ffn = paste0(fielddatadir_wsl, rover_sdn, rover_fn)
  
  
  print(paste("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"))
  print(paste("++++   Start processing data from Rover:", Rover," - Occ#", Occ,"at peg # ", pegnum, "using method", Method))
  print(paste("++++   Time Period: ", Period))
  print(paste("--------------------------------------------------------------------------------------"))
  
  proc_metadata = list(Rover = Rover, Occ = Occ, Peg = pegnum, Method = Method, Period=Period)
  
  
  
  base_obs_ffn = case_when(
    Method == "RTK-R10" & Occ == 1 ~ paste0(fielddatadir_wsl, "R10-Base/RINEX211/01393111.23o"),
    Method == "RTK-R10" & Occ == 2 ~ paste0(fielddatadir_wsl, "R10-Base/RINEX211/01393120.23o"),
    Method == "RTK-CORS-DNVK" & Occ == 1 ~ paste0(fielddatadir_wsl, "POSITIONZ-DNVK/OCC1/dnvk3110_20.23o"),
    Method == "RTK-CORS-DNVK" & Occ == 2 ~ paste0(fielddatadir_wsl, "POSITIONZ-DNVK/OCC2/dnvk3120.23o"),
    Method == "RTK-CORS-WANG" & Occ == 1 ~ paste0(fielddatadir_wsl, "POSITIONZ-WANG/OCC1/wang3110_20.23o"),
    Method == "RTK-CORS-WANG" & Occ == 2 ~ paste0(fielddatadir_wsl, "POSITIONZ-WANG/OCC2/wang3120.23o"),
    Method == "RTK-CORS-WRPA" & Occ == 1 ~ paste0(fielddatadir_wsl, "POSITIONZ-WRPA/OCC1/wrpa3110_20.23o"),
    Method == "RTK-CORS-WRPA" & Occ == 2 ~ paste0(fielddatadir_wsl, "POSITIONZ-WRPA/OCC2/wrpa3120.23o"),
    Method == "RTK-CORS-RDLV" & Occ == 1 ~ paste0(fielddatadir_wsl, "POSITIONZ-RDLV/OCC1/rdlv3110_20.23o"),
    Method == "RTK-CORS-RDLC" & Occ == 2 ~ paste0(fielddatadir_wsl, "POSITIONZ-RDLV/OCC2/rdlv3120.23o"),
    
  )
  
  base_nav_ffn = case_when(
    Method == "RTK-R10" & Occ == 1 ~ paste0(fielddatadir_wsl, "R10-Base/RINEX211/01393111.23n"),
    Method == "RTK-R10" & Occ == 2 ~ paste0(fielddatadir_wsl, "R10-Base/RINEX211/01393120.23n"),
    Method == "RTK-CORS-DNVK" & Occ == 1 ~ paste0(fielddatadir_wsl, "POSITIONZ-DNVK/OCC1/auto3110_20.23n"),
    Method == "RTK-CORS-DNVK" & Occ == 2 ~ paste0(fielddatadir_wsl, "POSITIONZ-DNVK/OCC2/auto3120.23n"),
    Method == "RTK-CORS-WANG" & Occ == 1 ~ paste0(fielddatadir_wsl, "POSITIONZ-WANG/OCC1/auto3110_20.23n"),
    Method == "RTK-CORS-WANG" & Occ == 2 ~ paste0(fielddatadir_wsl, "POSITIONZ-WANG/OCC2/auto3120.23n"),
    Method == "RTK-CORS-WRPA" & Occ == 1 ~ paste0(fielddatadir_wsl, "POSITIONZ-WRPA/OCC1/auto3110_20.23n"),
    Method == "RTK-CORS-WRPA" & Occ == 2 ~ paste0(fielddatadir_wsl, "POSITIONZ-WRPA/OCC2/auto3120.23n"),
    Method == "RTK-CORS-RDLV" & Occ == 1 ~ paste0(fielddatadir_wsl, "POSITIONZ-RDLV/OCC1/auto3110_20.23n"),
    Method == "RTK-CORS-RDLV" & Occ == 2 ~ paste0(fielddatadir_wsl, "POSITIONZ-RDLV/OCC2/auto3120.23n"),
    
  )
  
  precise_eph_ffn = case_when(
    Occ == 1 ~ paste0(fielddatadir_wsl, "PREC-EPH-CLK/CODE Precise/COD0OPSFIN_20233110000_01D_05M_ORB.SP3"),
    Occ == 2 ~ paste0(fielddatadir_wsl, "PREC-EPH-CLK/CODE Precise/COD0OPSFIN_20233120000_01D_05M_ORB.SP3"),
    
    
  )
  
  precise_clk_ffn = case_when(
    Occ == 1 ~ paste0(fielddatadir_wsl, "PREC-EPH-CLK/CODE Precise/COD0OPSFIN_20233110000_01D_05S_CLK.CLK"),
    Occ == 2 ~ paste0(fielddatadir_wsl, "PREC-EPH-CLK/CODE Precise/COD0OPSFIN_20233120000_01D_05S_CLK.CLK"),
    
    
  )
  
  
  FileTypeVec =     c("Rover_obs", "Base_obs", "Base_nav", "Precise_eph", "Precise_clk")
  FilesReqd_IFO_0 = c(1          ,  0        ,      0    ,    0         ,     0        )
  FilesReqd_IFO_1 = c(1          ,  1        ,      1    ,    0         ,     0        )
  FilesReqd_IFO_2 = c(1          ,  1        ,      1    ,    1         ,     1        )
  
  FilesReqd = case_when(
    InputFilesOption == 0 ~ FilesReqd_IFO_0,
    InputFilesOption == 1 ~ FilesReqd_IFO_1,
    InputFilesOption == 2 ~ FilesReqd_IFO_2
  )
  
  FilesList = c(rover_ffn,base_obs_ffn, base_nav_ffn, precise_eph_ffn, precise_clk_ffn )
  
  InputFileTable = data.frame(
    FileType = FileTypeVec, 
    Reqd = FilesReqd, 
    Found = 0, 
    FileName = FilesList) %>% 
    filter(Reqd==1) 
  
  NR = nrow(InputFileTable)
  
  for (i in 1:NR){
    
    cfilename = InputFileTable$FileName[i]
    if (file.exists(wineq_path(cfilename))){InputFileTable$Found[i]=1}
  }
  
  
  ALL_FOUND = sum(InputFileTable$Found) == NR
  
  if (ALL_FOUND){
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print("!!!......found all files........!!!!!!!!!")
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
  } else {
    
    InputFileTable_badrows = InputFileTable %>% 
      filter(Found==0)
    
    for (ii in 1:nrow(InputFileTable_badrows)){
      
      print(paste("Input File", InputFileTable_badrows$FileName[ii],"not found ..."))
      
    }
    
    
  }
  
  
  # 
  # if (file.exists(wineq_path(rover_ffn))){
  #   print(paste("Rover file:", rover_ffn, "exists"))
  #   files_fnd = files_fnd + 1
  # }else{
  #   print(paste("Rover file:", rover_ffn, "not found"))
  # } 
  # 
  # if (file.exists(wineq_path(base_obs_ffn))){
  #   print(paste("Base obs file:", base_obs_ffn, "exists"))
  #   files_fnd = files_fnd + 1
  # }else{
  #   print(paste("Base obs file:", base_obs_ffn, "not found"))
  # } 
  # 
  # if (file.exists(wineq_path(base_nav_ffn))){
  #   print(paste("Base nav file:", base_nav_ffn, "exists"))
  #   files_fnd = files_fnd + 1
  # }else{
  #   print(paste("Base nav file:", base_nav_ffn, "not found"))
  # } 
  # 
  # if (file.exists(wineq_path(precise_eph_ffn))){
  #   print(paste("Base nav file:", precise_eph_ffn, "exists"))
  #   files_fnd = files_fnd + 1
  # }else{
  #   print(paste("Preceise ephemeris data file:", precise_eph_ffn, "not found"))
  # } 
  # 
  # if (file.exists(wineq_path(precise_clk_ffn))){
  #   print(paste("Base nav file:", precise_clk_ffn, "exists"))
  #   files_fnd = files_fnd + 1
  # }else{
  #   print(paste("Base nav file:", precise_clk_ffn, "not found"))
  # } 
  # 
  # print("=================================================================")
  # files_fnd
  
  if (ALL_FOUND){
    
    # Build the name of the solution file
    Solution_fn = paste0("BAL_OCC",Occ,"_", Rover,"_", Method, "_", Period, ".pos" )
    Solution_ffn = paste0(solution_dn_wsl, Solution_fn)
    Solution_ffn_w10 = paste0(solution_dn_w10, Solution_fn)
    #now build the cmd
    
    qm = "\""
    cmd_stub = paste0("wsl rnx2rtkp -o ", qm, Solution_ffn, qm," -p 3 -f 3 -s \",\" -c")
    
    # buld the time mark portion of the command if it is used
    if (Period != "FULL" ){
      
      # 
      # Period = "30MIN_TO_60MIN"
      # Period = "30MIN"
      
      
      start_AND_stop = str_detect(Period, "_TO_")|str_detect(Period, "_to_")
      
      if (start_AND_stop){
        start_MIN = as.numeric(str_extract(Period, "^\\d{1,4}"))
        end_MIN = as.numeric(str_extract(Period, "(?<=\\d{1,4}MIN\\_TO\\_)\\d{1,4}"))
      } else {
        start_MIN = 0
        end_MIN = as.numeric(str_extract(Period, "^\\d{1,4}"))
      }
      
      
      #get the start time of the measurement
      
      tm_start_obs = current_metadata_row$OCC_START
      tm_end_obs = current_metadata_row$OCC_END
      
      tm_start_proc = tm_start_obs + minutes(start_MIN)
      tm_end_proc = tm_start_obs + minutes(end_MIN)
      
      fmt = "%Y/%m/%d %H:%M:%S"
      
      tm_start_proc_txt = strftime(tm_start_proc, fmt, tz = "UTC")
      tm_end_proc_txt = strftime(tm_end_proc, fmt, tz = "UTC")
      
      timeSliceString = paste0("-ts ", tm_start_proc_txt, " -te ", tm_end_proc_txt) 
      
      cmd_stub_mdf = paste0(cmd_stub, " ", timeSliceString)
      
    } else {
      
      cmd_stub_mdf = cmd_stub
      
    }
    
    #now create full CMD string
    
    CMD = paste0(cmd_stub_mdf,  " ",
                 qm, rover_ffn, qm, " ",
                 qm, base_obs_ffn, qm, " ",
                 qm, base_nav_ffn, qm, " ",
                 qm, precise_eph_ffn, qm, " ",
                 qm, precise_clk_ffn, qm)
    
    
    
    print(CMD)
    
    
    
    
    if (!replot_only){ 
      system(CMD)
      
      print(paste("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"))
      print(paste("++++   Finished processing data from Rover:", Rover," - Occ#", Occ,"at peg # ", pegnum, "using method", Method))
      print(paste("++++   Time Period: ", Period))
      print(paste("++++   Solution saved to: ", basename(Solution_ffn_w10)))
      print(paste("++++   in directory: ", dirname(Solution_ffn_w10)))
      print(paste("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"))
      # print(paste("Now plotting results"))
      # print(paste("Done"))
      # print(paste("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"))
      
    } 
    
    # PLOTDN = paste0("C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/plots/Solution-Plots/", Rover, "/", Method, "/")
    # if (!dir.exists(PLOTDN)){dir.create(PLOTDN, recursive = T)}
    # 
    # 
    # POS_ANL_RESULTS = analyze_POS(Rover, OCC=Occ, PEG=pegnum, PERIOD=Period, Solution_ffn_w10, plotdir = PLOTDN, plotdim = plotdim)
    
    
    
    
    
    
    
    
    
    POS_DATA = read_POS(Solution_ffn_w10)
    
    OUT = list(proc_metadata, POS_DATA, POS_FFN = Solution_ffn_w10, CMD)
    
  } else {
    
    print("A file is missing ....")
    
    
  }
  
  
  
}

#-----------------------------------------------------------------------#
# analyze_POS: Analyses an RTKLIB Solution "POS" file (calls function read_POS)
#-----------------------------------------------------------------------#

analyze_POS = function(RCVR, OCC, PEG, METHOD = NULL, PERIOD, POS_ffn, plotdir = "C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/plots/Solution-Plots/" , plotdim = 3, limax = T, tm_win = NULL, FMT_OPT=1, qualfilt = NULL){
  # 
  # RCVR = "JVD"
  # OCC = 2
  # PEG = 4
  # PERIOD = "60 MIN"
  # METHOD = "RTK-R10"
  # FMT_OPT = 1
  # plotdim = 3
  # limax = T
  # tm_win = NULL
  # 
  # #
  # # POS_ffn = "C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/data/Field-Test-Ballance-20231108/JVD/SOLUTIONS/OCC1-PPK-R10BASE/SOLN_101.pos"
  # POS_ffn = "C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/data/Field-Test-Ballance-20231108/JVD/SOLUTIONS/OCC1-PPK-R10BASE/SOLN_102.pos"
  
  # POS_ffn = "C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/
  # data/Field-Test-Ballance-20231108/SOLUTIONS/BAL_OCC1_JVD_RTK-R10_60MIN.pos"
  # 
  
  
  # plotdir = "C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/plots/Solution-Plots/"
  plotdir = "C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/data/Field-Test-Ballance-20231108/SOLUTIONS/plots/"
  
  KnownPegPositions_binary_ffn = "C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/data/KnownPegPositions.RDS"
  
  KnownPegPositions = readRDS(KnownPegPositions_binary_ffn)
  
  known_loc = data.frame(Peg = PEG, X = KnownPegPositions$EASTING[PEG], Y = KnownPegPositions$NORTHING[PEG])
  
  
  
  
  D = read_POS(POS_ffn, FMT_OPT = FMT_OPT)
  
  if (!is.null(qualfilt)){
    
    D = D %>% filter(Q %in% qualfilt)
    
    
    
  }
  
  
  if (!is.null(tm_win)){
    D = D %>% filter(melap >= tm_win[1] & melap <= tm_win[2])
  }
  
  # D_sf = D %>% 
  #   st_as_sf(coords = c("Lon_DD", "Lat_DD"), crs = 7912) %>% 
  #   st_transform(4326) %>% 
  #   st_transform(2193)  
  # 
  D_sf = D %>% 
    st_as_sf(coords = c("Lon_DD", "Lat_DD"), crs = 4326) %>% 
    st_transform(2193)  %>% 
    mutate(
      melap = D$melap
    )
  
  # D_sf = D %>% 
  #  st_as_sf(coords = c("Lon_DD", "Lat_DD"), crs = 2193)
  # 
  
  
  
  
  
  CRDS = D_sf %>% 
    st_coordinates() %>%  
    as.data.frame() %>% 
    as_tibble() %>% 
    mutate(melap = D_sf$melap)
  
  
  
  
  OUTPUT_fig1_g = anl_pos_fig_1(D_sf, RCVR, OCC, PEG, METHOD, PERIOD,known_loc,  plotdim = 3, plotdir, limax)
  
  #get the semimajor angle
  # SM = find_semi_major(CRDS, known_loc)
  
  #now draw an ellipse
  # g = plot_ellipse(OUTPUT_fig1_g, known_loc$X,known_loc$Y, SM$Semiminor_2SD, SM$Semimajor_2SD, SM$Semimajor_azimuth)
  # g  
  # 
  # g = ggplot(CRDS, aes(X,Y)) + geom_point()
  # g = ggMarginal(g,CRDS, X,Y )
  # g = plot_ellipse(g, known_loc$X,known_loc$Y, SM$Semiminor_2SD, SM$Semimajor_2SD, SM$Semimajor_azimuth)
  # g 
  
  OUTPUT_AP_table = calc_acc_and_prec(CRDS, known_loc)
  
  
  
  #-----------------------------------------------------------------------#
  # collate_results
  #-----------------------------------------------------------------------#
  
  CollateDF_ffn =  "C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/data/Field-Test-Ballance-20231108/SOLUTIONS/GNSS_PERFORMANCE_STATS_COLLATED.csv"
  
  
  PX_UTC = GPS_tm_2_UTC(D$GPS_Week, D$GPS_Time)
  TM_STT = strftime(min(PX_UTC), format = "%Y/%m/%d %H:%M:%S", tz="UTC")
  TM_ENN = strftime(max(PX_UTC), format = "%Y/%m/%d %H:%M:%S", tz="UTC")
  NSV_AVG = mean(D$nvsv)
  PDOP_AVG = NA
  X_AVG= mean(CRDS$X)
  Y_AVG= mean(CRDS$Y)
  
  DATAROW = data.frame(
    RCVR = RCVR,
    PEG = PEG,
    KNOWN_X = known_loc$X,
    KNOWN_Y = known_loc$Y,
    METHOD = METHOD,
    TM_STT = TM_STT,
    TM_ENN = TM_ENN,
    NMEAS = OUTPUT_AP_table$N_MEAS,
    MIN_ELAPSED = OUTPUT_AP_table$minutes_elapsed,
    X_AVG = X_AVG,
    Y_AVG = Y_AVG,
    SIGMA_X = OUTPUT_AP_table$SIGMA_X,
    SIGMA_Y = OUTPUT_AP_table$SIGMA_Y,
    M_ACC = OUTPUT_AP_table$M_ACC,
    M_PREC = OUTPUT_AP_table$M_PREC,
    DRMS_65 = OUTPUT_AP_table$DRMS_65,
    DRMS_95 = OUTPUT_AP_table$DRMS_95,
    CEP_50 = OUTPUT_AP_table$CEP_50,
    CEP_95 = OUTPUT_AP_table$CEP_95,
  )
  
  
  
  
  COLLATED_DATA = read_csv(CollateDF_ffn)
  COLLATED_DATA_NEW = COLLATED_DATA %>% bind_rows(DATAROW)
  
  write_csv(COLLATED_DATA_NEW, CollateDF_ffn)
  
  
  OUTACC = paste("Difference between Known position and Measured Position = ", formatC(DATAROW$M_ACC, digits = 7), " +/- ", formatC(DATAROW$M_PREC, digits = 7), "Std Dev")
  print(OUTACC)
  
  OUT = list(OUTACC, DATAROW, CRDS, OUTPUT_fig1_g)
  return(OUT)
  
}




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### rtk_proc_anal
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# This is a wrapper function to process raw rinex data using RTKLIB (via the R function proc_rtk) and then analyse it
# (using the R function analyze_POS)


#Check to make sure these three files exist
rtk_proc_anal = function(Rover, Occ, Method, Period, plotdim = 3, replot_only = F){
  # gps_week = c(2287,2288)
  # gps_sec_of_current_week = c(250846,250846)
  
  
  print(paste("==================================================================="))
  print(paste("==================================================================="))
  print(paste("New Analysis: Rover =",Rover,"Occ = ", Occ, "Method = ",Method, "Period = ", Period))
  print(paste("==================================================================="))
  print(paste("==================================================================="))
  #set the plot directory
  PLOTDN = "C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/plots/Solution-Plots/"
  
  #process the raw GNSS data
  RTKLIB_RESULTS = proc_rtk(Rover, Occ, Method, Period, plotdim = 3, replot_only = replot_only)
  
  
  #grab the metadata
  pegnum = RTKLIB_RESULTS[[1]]$Peg
  method = RTKLIB_RESULTS[[1]]$Method
  period = RTKLIB_RESULTS[[1]]$Period
  #grab the file name of the solution file (POS file)
  POS_ffn = RTKLIB_RESULTS[[3]]
  
  #analyse the solution file
  POS_ANAL_RESULTS = analyze_POS(RCVR=Rover, OCC = Occ, PEG=pegnum, METHOD = method, PERIOD=period, POS_ffn, plotdir = PLOTDN , plotdim = 3, limax = T, tm_win = NULL, FMT_OPT=1)
  
  
  
  
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### get_qual_stats
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

get_qual_stats = function(POS_data_input = NULL, POS_ffn = NULL){
  
  # POS_ffn = "C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/data/Field-Test-Ballance-20231108/SOLUTIONS/BAL_OCC1_R10-Rover_RTK-CORS-WRPA_180MIN.pos"
  # POS_data_input = NULL
  
  if (is.null(POS_data_input) & !is.null(POS_ffn)){
    
    POS_data = read_POS(POS_ffn)
    
  }
  
  
  qual_stats = POS_data %>% group_by(Q) %>%
    summarise(qual_n = length(Q)) %>% 
    mutate(qual_perc = 100 * qual_n/sum(qual_n))
  
  return(qual_stats)
  
  
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### qual_based_soln
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

qual_based_soln = function(POS_data_input = NULL, POS_ffn = NULL, known_coords = NULL){
  
  testmode = F
  
  if (testmode){
    POS_ffn = "C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/data/Field-Test-Ballance-20231108/SOLUTIONS/BAL_OCC1_R10-Rover_RTK-CORS-WRPA_180MIN.pos"
    POS_data_input = NULL
    known_coords = c(1838917.665,	5531119.762,	102.372)
    known_coords = c(1839374.859,	5531494.916,	79.451)
  }
  
 
  
  
  #work out if the input is a file name or a data framne
  if (is.null(POS_data_input) & !is.null(POS_ffn)){
    POS_data_mdf = read_POS(POS_ffn)
  } else if (!is.null(POS_data_input) & is.null(POS_ffn)){
    POS_data_mdf = POS_data_input
  }
  
  #convert coordinates
  POS_data_with_COORDS = convert_crs_df(POS_data_mdf, src_crs = 4326, dst_crs = 2193, 
                            input_coord_labs = c("Lon_DD", "Lat_DD"),
                            output_coord_labs = c("LON", "LAT"))
  
  
  # if known coordinates are supplied then create a boolean called pos_known
  if (!is.null(known_coords)){pos_known = T}else{pos_known = T}
  
  # if known coordinates are supplied then add difference columns
  if (pos_known){POS_data_with_COORDS = add_diff_from_known(POS_data_with_COORDS, known_coords =   known_coords) }
  
  
  
 
  POS_data_stats = calc_POS_data_stats(POS_data_with_COORDS, known_coords = known_coords) 
  
  print(t(POS_data_stats))
  
  POS_data_fixed = POS_data_all %>% filter(Q==1)
  
  
  return(t(POS_data_stats))
  
  
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### calc_POS_data_stats
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


calc_POS_data_stats = function(POS_data, known_coords=NULL){
  
  if (!is.null(known_coords)){
    POS_data_all_stats = POS_data %>% 
      summarise(
        LAT_avg = mean(LAT),
        LON_avg = mean(LON),
        ELE_avg = mean(Ht_ell),
        LAT_sd = sd(LAT),
        LON_sd = sd(LON),
        ELE_sd = sd(Ht_ell),
        DEL_LON_avg = mean(DEL_LON),
        DEL_LAT_avg = mean(DEL_LAT),
        DEL_LON_sd = sd(DEL_LON),
        DEL_LAT_sd = sd(DEL_LAT),
        KNOWN_POS_X = known_coords[1],
        KNOWN_POS_Y = known_coords[2],
        M_ACC = sqrt( (LON_avg-KNOWN_POS_X)^2 + (LAT_avg-KNOWN_POS_Y)^2  ),
        M_PREC = sqrt(DEL_LON_sd^2 + DEL_LAT_sd^2),
        
        LON_avg_fixed = mean(LON[Q==1]),
        LAT_avg_fixed = mean(LAT[Q==1]),
        
        LON_avg_float = mean(LON[Q==2]),
        LAT_avg_float = mean(LAT[Q==2]),
        
        LON_avg_fixed_or_float = mean(LON[Q==1|Q==2]),
        LAT_avg_fixed_or_float = mean(LAT[Q==1|Q==2]),
        
        M_ACC_all = sqrt( (LON_avg-KNOWN_POS_X)^2 + (LAT_avg-KNOWN_POS_Y)^2  ),
        M_ACC_fixed = sqrt( (LON_avg_fixed - KNOWN_POS_X)^2 + (LAT_avg_fixed - KNOWN_POS_Y)^2),
        M_ACC_float = sqrt( (LON_avg_float - KNOWN_POS_X)^2 + (LAT_avg_float - KNOWN_POS_Y)^2),
        M_ACC_fixed_or_float = sqrt( (LON_avg_fixed_or_float - KNOWN_POS_X)^2 + (LAT_avg_fixed_or_float - KNOWN_POS_Y)^2),
        PERCENT_FIXED = 100 * length(which(Q==1))/length(Q),
        PERCENT_FLOAT = 100 * length(which(Q==2))/length(Q)
        
      )
    
  }else{
    POS_data_all_stats = POS_data %>% 
      summarise(
        LAT_avg = mean(LAT),
        LON_avg = mean(LON),
        ELE_avg = mean(Ht_ell),
        LAT_sd = sd(LAT),
        LON_sd = sd(LON),
        ELE_sd = sd(Ht_ell),
        M_PREC = sqrt(LON_sd^2 + LAT_sd^2)
      )
  }
  return(POS_data_all_stats)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### add_diff_from_known
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# a funtion to add columns based on known coordinates
add_diff_from_known = function(POS_DF, known_coords){
  
  elev_known = length(known_coords)==3
  
  POS_DF_mdf01 = POS_DF %>% 
    mutate(
      KNOWN_POS_X = known_coords[1],
      KNOWN_POS_Y = known_coords[2],
      KNOWN_POS_Z = ifelse(elev_known, known_coords[3], NA),
      DEL_LON = LON - KNOWN_POS_X,
      DEL_LAT = LAT - KNOWN_POS_Y,
      DEL_ELE = ifelse(elev_known,Ht_ell - KNOWN_POS_Z,NA),
      DEL_LON_SQ =  DEL_LON^2,
      DEL_LAT_SQ =  DEL_LAT^2,
      DEL_LATLON = sqrt(DEL_LON_SQ + DEL_LAT_SQ)
    )
  
  return(POS_DF_mdf01)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### convert_crs_df: input a data frame with X and Y in one projection
###        output data frame with transformed coordinates
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+

convert_crs_df = function(df, src_crs = 4326, dst_crs = 2193, 
                          input_coord_labs = c("Lon_DD", "Lat_DD"),
                          output_coord_labs = c("LON", "LAT")){
  
  # df = SOLN_DATA
  
  transformed_coords = df %>% 
    st_as_sf(coords = coord_labs, crs = src_crs) %>% 
    st_transform(crs = dst_crs) %>% 
    st_coordinates() %>% 
    as_tibble() 
  
  names(transformed_coords) = output_coord_labs
  
  df_out = df %>% bind_cols(transformed_coords)
  
  return(df_out)
  
  
}



proc_rtk_generic = function(InputFile_ffn = NULL, Output_ffn,Rover_obs_ffn, Base_obs_ffn, Base_nav_ffn, Precise_eph, Precise_clk){
  
  testmode = F
  if (testmode){
    
    
    wsl_rootDrive = "/mnt/c/"
    w10_rootDrive = "C:/"
    
    
    
    root_wsl = paste0(wsl_rootDrive, "Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/data/Field-Test-Ballance-20231108/")
    root_w10 = paste0(w10_rootDrive, "Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/data/Field-Test-Ballance-20231108/")
    
    rover_obs_ffn     = paste0(root_wsl, "R10-Base/Bens-Rinex-Files/01393111.23o")
    base_obs_ffn  = paste0(root_wsl,"POSITIONZ-DNVK/OCC1/dnvk3110_20.23o")
    base_nav_ffn  = paste0(root_wsl,"POSITIONZ-DNVK/OCC1/auto3110_20.23n")
    precise_eph_ffn = paste0(root_wsl)
    precise_clk_ffn = paste0(root_wsl, "PREC-EPH-CLK/CODE Precise/COD0OPSFIN_20233110000_01D_05S_CLK.CLK")
    Solution_ffn  = paste0(root_wsl,"SOLUTIONS/R10-BASE-RTKLIB_SOLUTION.pos")
    
    
    w10_Solution_ffn = paste0(root_w10,"SOLUTIONS/R10-BASE-RTKLIB_SOLUTION.pos")
    
    InputFile_ffn = "C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/data/proc_rtk_generic_filelists/RTKLIB-inputfilelist-R10-Base-01393111.txt"
    
    
  }
  
  if(!is.null(InputFile_ffn)){
    
    lines = readLines(InputFile_ffn)
    remove("config_ffn")
    for (iline in lines){
      
      
      
      varname = str_trim(str_split(iline,"\\=")[[1]][1])
      ffn_value = str_trim(str_split(iline,"\\=")[[1]][2])
      assign(varname, ffn_value)
      print(varname)
      print(ffn_value)
    }
  }
  
  if (exists("config_ffn")){config=T}else(config=F)
  
  files2find = c(rover_obs_ffn, base_obs_ffn, precise_eph_ffn, precise_clk_ffn)
  for (ifile in files2find){
    
    ifile_w10 = wineq_path(ifile)
    if (!file.exists(ifile_w10)){
      print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print(paste("Could not find ", wineq_path(ifile_w10)))
      stop()
    }
    
  }
  
  
  
  qm = "\""
  
  if (config){
    
    cmd_stub = paste0("wsl rnx2rtkp -k ", qm, config_ffn, qm, " -o ", qm,Solution_ffn,qm)
  } else {
     cmd_stub = paste0("wsl rnx2rtkp -o ", qm, Solution_ffn, qm," -p 3 -f 1 -s \",\" -c")
  }
  # 
  # 
  # start_AND_stop = str_detect(Period, "_TO_")|str_detect(Period, "_to_")
  # 
  # if (start_AND_stop){
  #   start_MIN = as.numeric(str_extract(Period, "^\\d{1,4}"))
  #   end_MIN = as.numeric(str_extract(Period, "(?<=\\d{1,4}MIN\\_TO\\_)\\d{1,4}"))
  # } else {
  #   start_MIN = 0
  #   end_MIN = as.numeric(str_extract(Period, "^\\d{1,4}"))
  # }
  # 
  # 
  #get the start time of the measurement
  # 
  # tm_start_obs = current_metadata_row$OCC_START
  # tm_end_obs = current_metadata_row$OCC_END
  # 
  # tm_start_proc = tm_start_obs + minutes(start_MIN)
  # tm_end_proc = tm_start_obs + minutes(end_MIN)
  # 
  # fmt = "%Y/%m/%d %H:%M:%S"
  # 
  # tm_start_proc_txt = strftime(tm_start_proc, fmt, tz = "UTC")
  # tm_end_proc_txt = strftime(tm_end_proc, fmt, tz = "UTC")
  
  # timeSliceString = paste0("-ts ", tm_start_proc_txt, " -te ", tm_end_proc_txt) 
  
  # cmd_stub_mdf = paste0(cmd_stub, " ", timeSliceString)
  #   
  # } else {
  #   
  #   cmd_stub_mdf = cmd_stub
  #   
  # }
  # 
  #now create full CMD string
  
  #delete the old solution file so that there is no confusion
  if (file.exists(wineq_path(Solution_ffn))){
    file.rename(wineq_path(Solution_ffn), paste0(tools::file_path_sans_ext(wineq_path(Solution_ffn)), "_old_", strftime(Sys.time(), format = "%Y-%m-%d_%H%M%S"), ".txt"))
}
  
  cmd_stub_mdf = cmd_stub 
  
  CMD = paste0(cmd_stub_mdf,  " ",
               qm, rover_obs_ffn, qm, " ",
               qm, base_obs_ffn, qm, " ",
               qm, base_nav_ffn, qm, " ",
               qm, precise_eph_ffn, qm, " ",
               qm, precise_clk_ffn, qm)
  
  
  
  print(CMD)
  
  system(CMD)
  
  SOLN_DATA = read_POS(wineq_path(Solution_ffn))
  
  QUAL_STATS = get_qual_stats(SOLN_DATA)
  
  
  
  print(paste("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"))
  print(paste("++++   Finished processing data from Rover file:", rover_obs_ffn))
  # print(paste("++++   Time Period: ", Period))
  print(QUAL_STATS)
  print(paste("++++   Solution saved to: ", basename(Solution_ffn)))
  print(paste("++++   in directory: ", dirname(Solution_ffn)))
  print(paste("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"))
  # print(paste("Now plotting results"))
  # print(paste("Done"))
  # print(paste("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"))
  
  
  w10_Solution_ffn = wineq_path(Solution_ffn)
  
  return(list(POS_ffn = w10_Solution_ffn, SOLN_DATA = SOLN_DATA,QUAL_STATS=QUAL_STATS))
  
} 



anal_POS_generic = function(POS_ffn, known_coords = NULL){
  
  # POS_ffn = "c:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/data/Field-Test-Ballance-20231108/SOLUTIONS/R10-BASE-RTKLIB_SOLUTION-01393120.pos"
  
  known_coords = c(1839374.859,	5531494.916)

  # SOLN_DATA = read_POS(POS_ffn,FMT_OPT=NULL) 
  # 
  
  # SOLN_DATA_mdf = convert_crs_df(SOLN_DATA, src_crs = 4326, dst_crs = 2193, 
  #                             input_coord_labs = c("Lon_DD", "Lat_DD"),
  #                             output_coord_labs = c("LON", "LAT"))
  # 
  
  
  SOLN_QBASED = qual_based_soln(POS_data_input = NULL, POS_ffn = POS_ffn, known_coords = known_coords)
  
  SOLN_QBASED_df = as.data.frame(t(SOLN_QBASED))
  print(paste("LAT = ",formatC(SOLN_QBASED_df$LAT_avg,digits = 11)))
  print(paste("LON = ",formatC(SOLN_QBASED_df$LON_avg,digits = 11)))
  print(paste("ELEV = ",formatC(SOLN_QBASED_df$ELE_avg,digits = 6)))
  
  print(paste("ACCURACY = ",formatC(SOLN_QBASED_df$M_ACC,digits = 4)))
  print(paste("PREC = ",formatC(SOLN_QBASED_df$M_PREC,digits = 4)))
  print(paste("PERCENT FIXED = ",formatC(SOLN_QBASED_df$PERCENT_FIXED,digits = 6)))
  
  return(SOLN_QBASED)
  
}







