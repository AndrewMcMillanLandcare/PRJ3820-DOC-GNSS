
# GNSS_utils written for PRJ3280


library(tidyverse)
library(sf)
library(mapview)
library(plotly)


# to include this set of utilities in other code:
#source("C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/Code/R/PRJ3820-DOC-GNSS/GNSS_utils.R")


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### crs_convert
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

crs_convert = function(csv_ffn = NULL, crs_src, crs_dst, plot_dst=T){
  
  
  # 
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




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### calc_acc_and_prec
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+

calc_acc_and_prec = function(GNSS_data, known_loc){
  
  n = nrow(GNSS_data)
  melap = max(GNSS_data$melap)
  
  GNSS_data = GNSS_data %>% 
    mutate(
      X_DEV2 = (X - known_loc$X)^2,
      Y_DEV2 = (Y - known_loc$Y)^2
    )
  
  GNSS_data_smy = GNSS_data %>% 
    summarise(
      M_X_tmp = sum(X_DEV2),
      M_Y_tmp = sum(Y_DEV2),
      SIGMA_X = sd(X_DEV2),
      SIGMA_Y = sd(Y_DEV2),
    ) %>% 
    mutate(
      M_X = sqrt(M_X_tmp/n),
      M_Y = sqrt(M_Y_tmp/n),
      M_ACC = sqrt(M_X^2 + M_Y^2),
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



analyze_POS = function(RCVR, OCC, PEG, PERIOD, POS_ffn, plotdir = "C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/plots/Solution-Plots/" , plotdim = 3, limax = T, tm_win = NULL){

  RCVR = "JVD"
  OCC = 1
  PEG = 6
  PERIOD = "60 MIN"
  # 
  POS_ffn = "C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/data/Field-Test-Ballance-20231108/JVD/SOLUTIONS/OCC1-PPK-R10BASE/SOLN_101.pos"
  POS_ffn = "C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/data/Field-Test-Ballance-20231108/JVD/SOLUTIONS/OCC1-PPK-R10BASE/SOLN_102.pos"

  FMT_OPT = 2
  
    
  # POS_ffn = "C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/data/Field-Test-Ballance-20231108/BAL_OCC1_JVD_RTK-R10_60MIN.pos"
  # plotdim = 3
  # limax = T
  # tm_win = NULL
  # plotdir = "C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/plots/Solution-Plots/"
  
  KnownPegPositions_binary_ffn = "C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/data/KnownPegPositions.RDS"
  
  KnownPegPositions = readRDS(KnownPegPositions_binary_ffn)
  
  known_loc = data.frame(Peg = PEG, X = KnownPegPositions$EASTING[PEG], Y = KnownPegPositions$NORTHING[PEG])
  
  
  if (FMT_OPT==2){
    NSKIP = 25
    D = read_csv(POS_ffn, skip = NSKIP, col_names =F)
    HDRS = c("GPS_DateTime",   "Lat_DD",  "Lon_DD",  "Ht_ell", "Q",
             "nvsv","sdn", "sde", "sdu", "sdne" , "sdeu","sdun", "age" , "ratio")
    names(D) <- HDRS
    D = D %>%
      dplyr::filter(!is.na(GPS_DateTime) | GPS_DateTime != "" | GPS_DateTime != "NA") %>% 
      mutate(
        GPS_Time = as.POSIXct(GPS_DateTime, format = "%Y/%m/%d %H:%M:%S"),
        melap = (as.numeric(GPS_Time) - as.numeric(GPS_Time[1]))/60)
    
  } else {
    
    NSKIP=11
    D = read_csv(POS_ffn, skip = NSKIP, col_names =F)
    HDRS = c("GPS_Week", "GPS_Time",  "Lat_DD",  "Lon_DD",  "Ht_ell", "Q",
             "nvsv","sdn", "sde", "sdu", "sdne" , "sdeu","sdun", "age" ,"ratio")
    names(D) <- HDRS
    D %>% pull(GPS_Week) %>% unique()
    D = D %>% 
      mutate(
        across(Lat_DD:ratio, as.numeric),
        melap = (GPS_Time - GPS_Time[1])/60) 
    
    
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
    st_transform(2193)  
  
  # D_sf = D %>% 
  #  st_as_sf(coords = c("Lon_DD", "Lat_DD"), crs = 2193)
  # 
  
  
  Fig1 = anl_pos_fig_1(D_sf, RCVR, OCC, PEG, PERIOD,known_loc,  plotdir)
  
  
  
  
  CRDS = D_sf %>% 
    st_coordinates() %>%  as.data.frame() %>% as_tibble() %>% mutate(melap = D$melap)
  
  
  AccPrec = calc_acc_and_prec(CRDS, known_loc)
  
  print(AccPrec)

  
  
  
  
  
  errors_df = data.frame(x_meas = CRDS$X, y_meas = CRDS$Y) %>% 
    as_tibble() %>% 
    mutate(x_dev = known_loc$X - x_meas, y_dev = known_loc$Y)
  
  
  # average of last 5 positions
  
  # figure 2. the change in the errors as a function of minutes elapsed
  
  g  = ggplot(D_sf) + geom_point(aes(x = melap, y = sdne))
  g
  ggplotly(g)
  # average of last 5 positions
  
  
  # ggplotly(g)
  
  # 
  # mapview(D_sf, zcol = "melap")
  # 
  # ggplot(D) + geom_point(aes(Lon_DD, Lat_DD))
  # 
  # errors = D %>% 
  #   gather(key = errortype, value = value, sdn:sdeu)
  # 
  # ggplot(errors) + geom_point(aes(melap, value, colour =errortype))
  
  
  
  print(paste("Difference between Known position and Measured Position = ", DEL))
  
  return(g1)
  
}



#======================================================================================#
# figure 1. the xy plot of the solution
#======================================================================================#

anl_pos_fig_1 = function(D_sf, RCVR, OCC, PEG, PERIOD, known_loc, plotdir){
  
  CRDS = D_sf %>% 
    st_coordinates() %>%  as.data.frame() %>% as_tibble() %>% mutate(melap = D$melap)
  
  CRDS_last5 = CRDS %>% tail(10) %>% summarise(across(everything(), mean))
  
  DEL = sqrt( (known_loc$X - CRDS$X)^2 + (known_loc$Y - CRDS$Y)^2 )
  
  AccPrec = calc_acc_and_prec(CRDS, known_loc)
  
  
  
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
  g = g + scale_x_continuous(breaks = seq(xlim[1], xlim[2]))
  g = g + scale_y_continuous(breaks = seq(ylim[1], ylim[2]))
  
  g = g + labs(x= "EASTING (m)", y = "NORTHING (m)", color = "Occupation\nTime\n(minutes)")
  g = g + scale_fill_continuous(guide = 'none')
  
  
  g = g + annotate("text", x = xlim[1] + diff(xlim)*.05, y = ylim[2] - diff(ylim)*.05, label = RCVR, hjust=0)
  g
  
  
  
  DESCR = paste("Peg", PEG," - Occupation", OCC, " - ", PERIOD)
  g = g + annotate("text", x = xlim[1] + diff(xlim)*.05, y = ylim[2] - diff(ylim)*.15, label = DESCR, hjust=0)
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
  
  plot_fname_g1 = paste0(RCVR,"-OCC#", OCC, "-PEG#", PEG,"-", PERIOD,".png" )
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
  
  DEL = sqrt( (known_loc$X - CRDS_last5$X)^2 + (known_loc$Y - CRDS_last5$Y)^2 )
  
  
  
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


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### proc_rtk
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

proc_rtk = function(Rover, Occ, Method, Period, plotdim = 3, replot_only = F){
  
  # Rover = "JVD"
  # Occ = 2
  # Method = "RTK-R10"
  # Period = "60MIN"
  # plotdim=3
  # Method = "RTK-CORS"
  
  print(paste("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"))
  print(paste("++++   Start processing data from Rover:", Rover," - Occ#", Occ,"at peg # ", pegnum, "using method", Method))
  print(paste("++++   Time Period: ", Period))
  
  print(paste("--------------------------------------------------------------------------------------"))
  
  
  
  datadir_w10 = "C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/data/"
  datadir_wsl = "/mnt/c/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/data/"
  
  fielddatadir_w10 = paste0(datadir_w10, "Field-Test-Ballance-20231108/")
  fielddatadir_wsl = paste0(datadir_wsl, "Field-Test-Ballance-20231108/")
  
  
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
  
  
  base_obs_ffn = case_when(
    Method == "RTK-R10" & Occ == 1 ~ paste0(fielddatadir_wsl, "R10-Base/RINEX211/01393111.23o"),
    Method == "RTK-R10" & Occ == 2 ~ paste0(fielddatadir_wsl, "R10-Base/RINEX211/01393120.23o"),
    Method == "RTK-CORS" & Occ == 1 ~ paste0(fielddatadir_wsl, "POSITIONZ-DNVK/OCC1/dnvk3110_20.23o"),
    Method == "RTK-CORS" & Occ == 2 ~ paste0(fielddatadir_wsl, "POSITIONZ-DNVK/OCC1/dnvk3120.23o"),
    
  )
  
  base_nav_ffn = case_when(
    Method == "RTK-R10" & Occ == 1 ~ paste0(fielddatadir_wsl, "R10-Base/RINEX211/01393111.23n"),
    Method == "RTK-R10" & Occ == 2 ~ paste0(fielddatadir_wsl, "R10-Base/RINEX211/01393120.23n"),
    Method == "RTK-CORS" & Occ == 1 ~ paste0(fielddatadir_wsl, "POSITIONZ-DNVK/OCC1/auto3110_20.23n"),
    Method == "RTK-CORS" & Occ == 2 ~ paste0(fielddatadir_wsl, "POSITIONZ-DNVK/OCC1/auto3120.23n"),
    
  )
  
  
  
  
  files_fnd = 0
  
  if (file.exists(wineq_path(rover_ffn))){
    print(paste("Rover file:", rover_ffn, "exists"))
    files_fnd = files_fnd + 1
  }else{
    print(paste("Rover file:", rover_ffn, "not found"))
  } 
  
  if (file.exists(wineq_path(base_obs_ffn))){
    print(paste("Base obs file:", base_obs_ffn, "exists"))
    files_fnd = files_fnd + 1
  }else{
    print(paste("Base obs file:", base_obs_ffn, "not found"))
  } 
  
  if (file.exists(wineq_path(base_nav_ffn))){
    print(paste("Base nav file:", base_nav_ffn, "exists"))
    files_fnd = files_fnd + 1
  }else{
    print(paste("Base nav file:", base_nav_ffn, "not found"))
  } 
  
  
  
  if (files_fnd==3){
    
    # Build the name of the solution file
    Solution_fn = paste0("BAL_OCC",Occ,"_", Rover,"_", Method, "_", Period, ".pos" )
    Solution_ffn = paste0(fielddatadir_wsl, Solution_fn)
    Solution_ffn_w10 = paste0(fielddatadir_w10, Solution_fn)
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
                 qm, base_nav_ffn, qm)
    
    
    
    print(CMD)
    
    
    
    
    if (!replot_only){ 
      system(CMD)
      
      print(paste("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"))
      print(paste("++++   Finished processing data from Rover:", Rover," - Occ#", Occ,"at peg # ", pegnum, "using method", Method))
      print(paste("++++   Time Period: ", Period))
      print(paste("++++   Solution saved to: ", basename(Solution_ffn_w10)))
      print(paste("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"))
      print(paste("Now plotting results"))
      print(paste("Done"))
      print(paste("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"))
      
    }
    
    
    PLOTDN = paste0("C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/plots/Solution-Plots/", Rover, "/", Method, "/")
    if (!dir.exists(PLOTDN)){dir.create(PLOTDN, recursive = T)}
    
    
    analyze_POS(Rover, OCC=Occ, PEG=pegnum, PERIOD=Period, Solution_ffn_w10, plotdir = PLOTDN, plotdim = plotdim)
    
    
  }
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### wineq_path
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Check to make sure these three files exist
wineq_path = function(PATH){
  PATHOUT = str_replace(base_obs_ffn, "\\/mnt\\/c", "C\\:\\")
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### read_rinex
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Check to make sure these three files exist
read_rinex = function(PATH){
  library(asteRisk)
  testGPSnav <- readGPSNavigationRINEX("C:/Users/McMillanAn/OneDrive - MWLR/Projects/PRJ3820-DOC-GNSS/data/Field-Test-Ballance-20231108/R10-Base/RINEX211/01393111.23n")
}
