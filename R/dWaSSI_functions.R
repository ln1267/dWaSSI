# dWaSSI-C Model functions

# Distributed wrap of dWaSSI model----
#' @title Distribution wrap of dWaSSI-C
#' @description FUNCTION_DESCRIPTION
#' @param sim.dates  list of all dates
#' @return outputs potential evapotranspiration (mm day-1)
#' @details For details see Haith and Shoemaker (1987)
#' @examples
#' \dontrun{
#' distHydroSim(str.date,
#'              hru.par = hru_par, hru.info = hru_info, hru.elevband = NULL,
#'              clim.prcp = stcroix$prcp.grid, clim.tavg = stcroix$tavg.grid,
#'              snow.flag = 0, progress.bar = TRUE)
#' }
#' @rdname distHydroSim
#' @export
#'
distHydroSim <- function(sim.dates,
                         par.sacsma = NULL,par.petHamon=NULL, par.snow17=NULL,par.routLah=NULL,
                         hru.info = NULL, hru.elevband = NULL,
                         clim.prcp = NULL, clim.tavg = NULL,
                         hru.lai=NULL,hru.lc.lai=NULL,huc.lc.ratio=NULL,
                         snow.flag = 0, progress.bar = TRUE,month=F) {

  # date vectors
  sim_date <- seq.Date(sim.dates[["Start"]]-years(warmup), sim.dates[["End"]], by = "month")
  #jdate <- as.numeric(format(sim_date, "%j"))
  jdate <-  c(15,46,76,107,137,168,198,229,259,290,321,351)
  ydate <- as.POSIXlt(sim_date)$year + 1900
  ddate <- as.POSIXlt(as.Date(paste0(ydate,"/12/31")))$yday + 1
  days_mon<-sapply (sim_date,numberOfDays)

  # Simulation parameters
  sim_num  <- length(sim_date)  # number of simulation steps (months)
  hru_num  <- nrow(hru.info)  # number of HRUs
  # Extract the sim period from the climate data for each HRU
  climate_date <- seq.Date(sim.dates[["Start_climate"]], sim.dates[["End_climate"]], by = "month")
  grid_ind <-which(climate_date %in% sim_date)
  hru_prcp <- clim.prcp[grid_ind,]
  hru_tavg <- clim.tavg[grid_ind,]

  # HRU parameters
  hru_num   <- nrow(hru.info)  # total number of hrus in the watershed
  hru_lat     <- hru.info$HRU_Lat # Latitude of each HRU (Deg)
  hru_lon     <- hru.info$HRU_Lon # Longitude of each HRU (Deg)
  hru_area    <- hru.info$HRU_Area # Area of each HRU (as %)
  hru_elev    <- hru.info$`HRU_Elev(m)` # Elevation of each HRU (m)
  hru_flowlen <- hru.info$`HRU_FlowLen(m)` # Flow length for each HRU (m)
  tot_area    <- sum(hru_area) # Total area of the watershed system

  #  Elevation band parameters
  if (!is.null(hru.elevband)) {
    elevband_num         <- (dim(hru.elevband)[2]-2)/2
    hru_elevband_frac     <- hru.elevband[,3:(elevband_num +2)]
    hru_elevband_elev_avg <- hru.elevband[,(elevband_num +3):ncol(hru.elevband)]
  }

  #Final surfaceflow and baseflow at the outlet
  FLOW_SURF <- vector(mode = "numeric", length = sim_num)
  FLOW_BASE <- vector(mode = "numeric", length = sim_num)

  # Loop through each HRU and calculate adjust temp and pet values
  if(progress.bar == TRUE) pb <- txtProgressBar(min=0, max=hru_num, style=3)
  for (h in 1:hru_num) {

    # Vectors for surface and baseflow generated at each hru
    flow_direct_tot <- vector(mode = "numeric", length = sim_num)
    flow_base_tot   <- vector(mode = "numeric", length = sim_num)

    # If elevation bands are not specified
    if (is.null(hru.elevband)) {

      # If snow module is active, calculate snow
      if (snow.flag == 1) {
        hru_mrain <- snow17(par_snow17[h,], hru_prcp[,h], hru_tavg[,h], hru_elev[h], jdate)
      } else {
        hru_mrain <- snow_melt(ts.prcp = hru_prcp[,h],ts.temp = hru_tavg[,h],snowrange = c(-5,1))$prcp
      }
      # PET using Hamon equation
      hru_pet <- hamon(par = par_petHamon[h], tavg = hru_tavg[,h], lat = hru_lat[h], jdate = jdate)
      hru_pet<-hru_pet*days_mon
      # Calculate hru flow for each elevation band
      ## calculate flow for each land cover type
      if(is.null(hru.lc.lai) | is.null(huc.lc.ratio)){

        if(!is.null(hru.lai)) hru_pet<- hru_pet*hru.lai[,h]*0.0222+0.174*hru_mrain+0.502*hru_pet+5.31*hru.lai[,h]

        if( month){
          out <- sacSma_mon(par = par_sacsma[h,], prcp = hru_mrain, pet = hru_pet)
        }else{
          out <- sacSma(par = par_sacsma[h,], prcp = hru_mrain, pet = hru_pet)
        }

      }else{
        # if lai data is not enough
        if(is.null(sim.dates[["Start_lai"]])) str.date.lai<-sim.dates[["Start_climate"]]
        if(is.null(sim.dates[["End_lai"]])) end.date.lai<-sim.dates[["End_climate"]]

        ## set the first year lai to before
        if(sim.dates[["Start_lai"]]>sim.dates[["Start"]]-years(warmup)){
          lack_lai_years<-floor(as.numeric(difftime(sim.dates[["Start_lai"]],sim.dates[["Start"]]-years(warmup),units="days")/365))
          hru.lc.lai<-lapply(hru.lc.lai,function(x) rbind(apply(x[1:12,],2,rep,lack_lai_years),x))
          sim.dates[["Start_lai"]]<-sim.dates[["Start"]]-years(warmup)
        }
        ## set the last year lai to after
        if(sim.dates[["End_lai"]]<sim.dates[["End"]]){
          lack_lai_years<-floor(as.numeric(difftime(sim.dates[["End"]],sim.dates[["End_lai"]],units="days")/365))
          start1<-length(hru.lc.lai[[1]][,1])-11
          hru.lc.lai<-lapply(hru.lc.lai,function(x) rbind(x,apply(x[start1:(start1+11),],2,rep,lack_lai_years)))
          sim.dates[["End_lai"]]<-sim.dates[["End"]]
        }

        # filter lai data by the simulation date
        lai_date <- seq.Date(sim.dates[["Start_lai"]], sim.dates[["End_lai"]], by = "month")
        grid_ind_lai <-which(lai_date %in% sim_date)
        hru.lc.lai<-lapply(hru.lc.lai,function(x) x[grid_ind_lai,] )

        # Calculate PET based on LAI
        hru_lc_pet<-hru_pet*hru.lc.lai[[h]]*0.0222+0.174*hru_mrain+0.502*hru_pet+5.31*hru.lc.lai[[h]]
        hru_lc_out <-apply(hru_lc_pet,2,sacSma_mon,par = par_sacsma[h,], prcp = hru_mrain )

        # Weight each land cover type based on it's ratio
        hru_out<-lapply(c(1:length(huc.lc.ratio[h,])),function(x) lapply(hru_lc_out[[x]],function(y) huc.lc.ratio[h,x]*y))
        out<-lapply(names(hru_lc_out[[1]]), function(var) apply(sapply(hru_out, function(x) x[[var]]),1,sum))
        names(out)<-names(hru_lc_out[[1]])
      }

      flow_direct_tot <- out$surf
      flow_base_tot   <- out$base

    } else {

      #Loop through each elevation band
      for (nn in 1:elevband_num) {

        if(hru_elevband_frac[h, nn] > 0) {

          # Adjusted temperature
          hru_tavg_adj <- hru_tavg[[h]] - 6.1 * (hru_elev[[h]] - hru_elevband_elev_avg[h, nn]) / 1000

          # Adjusted PET
          hru_pet_adj <- hamon(par = par_petHamon[h], tavg = hru_tavg_adj, lat = hru_lat[h], jdate = jdate)

          # If snow module is active, calculate snow
          if (snow.flag == 1) {
            hru_mrain <- snow17(par_snow17, hru_prcp[[h]], hru_tavg_adj, hru_elev[h], jdate)
          } else {
            hru_mrain <- hru_prcp[[h]]
          }

          # Calculate hru flow for each elevation band
          out <- sacSma(par = par_sacsma[h,], prcp = hru_mrain, pet = hru_pet_adj)
          flow_direct_tot <- flow_direct_tot + out$surf * hru_elevband_frac[h,nn]
          flow_base_tot   <- flow_base_tot   + out$base * hru_elevband_frac[h,nn]

        }

      } # close loop elev. bands

    }

    #---------- Channel flow routing from Lohmann routing model ---------------#
    out2 <- lohamann(par = par_routLah[h,], inflow.direct = flow_direct_tot,
                     inflow.base = flow_base_tot, flowlen = hru_flowlen[h])

    FLOW_SURF <- FLOW_SURF + out2$surf * hru_area[h] / tot_area
    FLOW_BASE <- FLOW_BASE + out2$base * hru_area[h] / tot_area

    if(progress.bar == TRUE) setTxtProgressBar(pb, h)
  }

  if(progress.bar == TRUE) close(pb)
  # get the real simulate period result
  FLOW_SURF<-FLOW_SURF[(warmup*12+1):length(FLOW_SURF)]
  FLOW_BASE<-FLOW_BASE[(warmup*12+1):length(FLOW_BASE)]
  out<- lapply(out, function(x) x[(warmup*12+1):length(x)])
  return(list(SURF = FLOW_SURF, BASE = FLOW_BASE,OUT=out))
}


# Parallel wrap of WaSSI-C model----
#' @title Parallel wrap of WaSSI-C model
#' @description FUNCTION_DESCRIPTION
#' @param sim.dates  list of all dates
#' @param warmup years for warming up WaSSI-C model
#' @param mcores how many cores using for simulation
#' @return outputs potential evapotranspiration (mm day-1)
#' @details For details see Haith and Shoemaker (1987)
#' @examples
#' \dontrun{
#' distHydroSim(sim.dates = sim_dates,
#'              hru.par = hru_par, hru.info = hru_info, hru.elevband = NULL,
#'              clim.prcp = stcroix$prcp.grid, clim.tavg = stcroix$tavg.grid,
#'              snow.flag = 0, progress.bar = TRUE)
#' }
#' @rdname dWaSSIC
#' @export
#'
dWaSSIC<- function(sim.dates, warmup=3,mcores=1,
                   par.sacsma = NULL,par.petHamon=NULL,par.routing=NULL,
                   hru.info = NULL,
                   clim.prcp = NULL, clim.tavg = NULL,
                   hru.lai=NULL,hru.lc.lai=NULL,huc.lc.ratio=NULL)
{

  # date vectors
  sim_date <- seq.Date(sim.dates[["Start"]]-years(warmup), sim.dates[["End"]], by = "month")
  #jdate <- as.numeric(format(sim_date, "%j"))
  jdate <-  c(15,46,76,107,137,168,198,229,259,290,321,351)
  ydate <- as.POSIXlt(sim_date)$year + 1900
  ddate <- as.POSIXlt(as.Date(paste0(ydate,"/12/31")))$yday + 1
  days_mon<-sapply (sim_date,numberOfDays)

  # Simulation parameters
  sim_num  <- length(sim_date)  # number of simulation steps (months)
  hru_num  <- nrow(hru.info)  # number of HRUs
  # Extract the sim period from the climate data for each HRU
  climate_date <- seq.Date(sim.dates[["Start_climate"]], sim.dates[["End_climate"]], by = "month")
  grid_ind <-which(climate_date %in% sim_date)
  hru_prcp <- clim.prcp[grid_ind,]
  hru_tavg <- clim.tavg[grid_ind,]

  # HRU parameters
  hru_num   <- nrow(hru.info)  # total number of hrus in the watershed
  hru_lat     <- hru.info$HRU_Lat # Latitude of each HRU (Deg)
  hru_lon     <- hru.info$HRU_Lon # Longitude of each HRU (Deg)
  hru_area<-NULL;hru_elev<-NULL;hru_flowlen<-NULL
  if("HRU_Area" %in% names(hru.info)) hru_area <- hru.info$HRU_Area # Area of each HRU (as %)
  if("HRU_Elev(m)" %in% names(hru.info)) hru_elev    <- hru.info$`HRU_Elev(m)` # Elevation of each HRU (m)
  if("HRU_FlowLen(m)" %in% names(hru.info)) hru_flowlen    <- hru.info$`HRU_FlowLen(m)` # Elevation of each HRU (m)
  # WaSSI for each HRU and calculate adjust temp and pet values

  WaSSI<-function(h){

    # Using snowmelt function to calculate effective rainfall

    hru_mrain <- snow_melt(ts.prcp = hru_prcp[,h],ts.temp = hru_tavg[,h],snowrange = c(-5,1))$prcp

    # PET using Hamon equation
    hru_pet <- hamon(par = par_petHamon[h], tavg = hru_tavg[,h], lat = hru_lat[h], jdate = jdate)
    hru_pet<-hru_pet*days_mon
    # Calculate hru flow for each elevation band
    ## calculate flow for each land cover type

    # if lai data is not enough
    if(is.null(sim.dates[["Start_lai"]])) str.date.lai<-sim.dates[["Start_climate"]]
    if(is.null(sim.dates[["End_lai"]])) end.date.lai<-sim.dates[["End_climate"]]

    ## set the first year lai to before
    if(sim.dates[["Start_lai"]]>sim.dates[["Start"]]-years(warmup)){
      lack_lai_years<-floor(as.numeric(difftime(sim.dates[["Start_lai"]],sim.dates[["Start"]]-years(warmup),units="days")/365))
      hru.lc.lai<-lapply(hru.lc.lai,function(x) rbind(apply(x[1:12,],2,rep,lack_lai_years),x))
      sim.dates[["Start_lai"]]<-sim.dates[["Start"]]-years(warmup)
    }
    ## set the last year lai to after
    if(sim.dates[["End_lai"]]<sim.dates[["End"]]){
      lack_lai_years<-floor(as.numeric(difftime(sim.dates[["End"]],sim.dates[["End_lai"]],units="days")/365))
      start1<-length(hru.lc.lai[[1]][,1])-11
      hru.lc.lai<-lapply(hru.lc.lai,function(x) rbind(x,apply(x[start1:(start1+11),],2,rep,lack_lai_years)))
      sim.dates[["End_lai"]]<-sim.dates[["End"]]
    }

    # filter lai data by the simulation date
    lai_date <- seq.Date(sim.dates[["Start_lai"]], sim.dates[["End_lai"]], by = "month")
    grid_ind_lai <-which(lai_date %in% sim_date)
    hru.lc.lai<-lapply(hru.lc.lai,function(x) x[grid_ind_lai,] )

    # Calculate PET based on LAI
    hru_lc_pet<-hru_pet*hru.lc.lai[[h]]*0.0222+0.174*hru_mrain+0.502*hru_pet+5.31*hru.lc.lai[[h]]
    hru_lc_out <-apply(hru_lc_pet,2,sacSma_mon,par = par_sacsma[h,], prcp = hru_mrain )

    # Weight each land cover type based on it's ratio
    hru_out<-lapply(c(1:length(huc.lc.ratio[h,])),function(x) lapply(hru_lc_out[[x]],function(y) huc.lc.ratio[h,x]*y))
    out<-lapply(names(hru_lc_out[[1]]), function(var) apply(sapply(hru_out, function(x) x[[var]]),1,sum))
    names(out)<-names(hru_lc_out[[1]])

    return(out)
  }

  huc_routing<-function(h){
    out2<-lohamann(par = par.routing[h,], inflow.direct = out[[h]][["surf"]],
                   inflow.base = out[[h]][["base"]], flowlen = hru_flowlen[h])
    FLOW_SURF <- out2$surf * hru_area[h] / sum(hru_area)
    FLOW_BASE <- out2$base * hru_area[h] / sum(hru_area)
    return(list(surf=FLOW_SURF,base=FLOW_BASE))
  }

  # get the real simulate period result
  if(mcores>1){
    result<-mclapply(c(1:hru_num), WaSSI,mc.cores = mcores)
  }else{
    result<-lapply(c(1:hru_num), WaSSI)
  }

  out<- lapply(result, function(x) as.data.frame(x)[(warmup*12+1):length(x[[1]]),])

  # routing based on catchment area
  if(is.null(par.routing & is.null(hru_flowlen) & is.null(hru_area))){

    return(out)

  }else{
    # Channel flow routing from Lohmann routing model
    out2 <- mclapply(c(1:hru_num), huc_routing,mc.cores = mcores)
    out3<-lapply(names(out2[[1]]), function(var) apply(sapply(out2, function(x) x[[var]]),1,sum))

  }

  return(list(FLOW_SURF = out3[[1]], FLOW_BASE=out3[[2]],HUC=out))

}

# Hamon Potential Evapotranspiration Equation----
#' @title Hamon Potential Evapotranspiration Equation
#' @description The Hamon method is also considered as one of the simplest estimates
#'     of potential Evapotranspiration.
#' @param par proportionality coefficient (unitless)
#' @param tavg  vector of mean daily temperature (deg C)
#' @param lat latitude ()
#' @param jdate a day number of the year (julian day of the year)
#' @return outputs potential evapotranspiration (mm day-1)
#' @details For details see Haith and Shoemaker (1987)
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1deg
#'  }
#' }
#' @rdname hamon
#' @export
hamon <- function(par, tavg, lat, jdate) {

  var_theta <- 0.2163108 + 2 * atan(0.9671396 * tan(0.0086 * (jdate - 186)))
  var_pi <- asin(0.39795 * cos(var_theta))
  daylighthr <- 24 - 24/pi * acos((sin(0.8333 * pi/180) + sin(lat * pi/180) * sin(var_pi))/(cos(lat *
                                                                                                  pi/180) * cos(var_pi)))

  esat <- 0.611 * exp(17.27 * tavg/(237.3 + tavg))

  return(par * 29.8 * daylighthr * (esat/(tavg + 273.2)))

}


# SNOW Model based on Tavg----
#' @title SNOW Model
#' @description Peter's Snow model
#' @param ts.prcp numeric vector of precipitation time-series (mm)
#' @param ts.temp numeric vector of average temperatuer time-series (Deg C)
#' @param snowrange (-3,1) the temperature range between snow and rain
#' @param meltmax  maximium ratio of snow melt (default: 0.5)
#' @param snowpack  initial snowpack (default: 0)
#' @return a dataframe of c("prcp","snowpack","snowmelt"), effective rainfall, snowpack and snowmelt
#' @rdname snowmelt
#' @details This is based on Peter's Fortran code
#' @examples
#' \dontrun{
#' mellt<-snow_melt(ts.prcp,ts.temp,meltmax=0.5,snowrange=c(-3,1),snowpack=0)
#' }
#' @export
snow_melt<-function(ts.prcp,ts.temp,meltmax=0.5,snowrange=c(-3,1),snowpack=0) {

  # Define outputs
  ts.snowpack <- vector(mode = "numeric", length = length(ts.prcp))
  ts.snowmelt <- vector(mode = "numeric", length = length(ts.prcp))
  ts.prcp.eff<- vector(mode = "numeric", length = length(ts.prcp))
  # loop each time step
  for (i in c(1:length(ts.prcp))){
    tavg<-ts.temp[i];prcp<-ts.prcp[i]
    # ---- FOR TEMPERATURE LESS THAN SNOW TEMP, CALCULATE SNOWPPT
    if (tavg <= snowrange[1]) {
      snow<-prcp
      rain<-0
      # ---- FOR TEMPERATURE GREATER THAN RAIN TEMP, CALCULATE RAINPPT
    }else if (tavg >= snowrange[2]){
      rain<-prcp
      snow<-0
      # ---- FOR TEMPERATURE BETWEEN RAINTEMP AND SNOWTEMP, CALCULATE SNOWPPT AND RAINPPT
    }else{
      snow<- prcp*((snowrange[2] - tavg)/(snowrange[2] - snowrange[1]))
      rain<-prcp-snow
    }
    # Calculate the snow pack
    snowpack<-snowpack+snow
    # ---- CALCULATE SNOW MELT FRACTION BASED ON MAXIMUM MELT RATE (MELTMAX) AND MONTHLY TEMPERATURE

    snowmfrac = ((tavg-snowrange[1])/(snowrange[2] - snowrange[1]))*meltmax

    if (snowmfrac >= meltmax) snowmfrac <- meltmax
    if (snowmfrac <0) snowmfrac <- 0

    # ---- CALCULATE AMOUNT OF SNOW MELTED (MM) TO CONTRIBUTE TO INFILTRATION & RUNOFF

    # ---- IF SNOWPACK IS LESS THAN 10.0 MM, ASSUME IT WILL ALL MELT (MCCABE AND WOLOCK, 1999)
    # ---- (GENERAL-CIRCULATION-MODEL SIMULATIONS OF FUTURE SNOWPACK IN THE WESTERN UNITED STATES)

    if (snowpack <= 10.0) {
      snowm<-snowpack
      snowpack<- 0

    }else{
      snowm = snowpack * snowmfrac
      snowpack = snowpack - snowm
    }

    # -- COMPUTE THE INFILTRATION FOR A GIVEN MONTH FOR EACH LAND USE
    ts.snowpack[i]<-snowpack
    ts.snowmelt[i]<-snowm
    ts.prcp.eff[i]<-rain+snowm
  }
  return(data.frame(prcp=ts.prcp.eff,snowpack=ts.snowpack,snowmelt=ts.snowmelt))
}


# Routing calculate average flow for hucs----
#' hucrouting
#'
#' Calculate accumulated flow for each huc, based on the flow direction
#' @param water this is the original flow data for each huc
#' @param routpar flow direction between hucs
#' @export
#' @examples
#' hrurouting(water,outpar)

hrurouting<-function(water,routpar,mc_cores=1){
  library(parallel)
  max_level<-max(routpar$LEVEL)

  hru_accm<-function(hru,water,routpar){
    hru<-as.numeric(hru)
    water$flow[water$HUC12==hru] +sum(water$flow[water$HUC12 %in% routpar$FROM[which(routpar$TO==hru)]])
  }

  for (level in c(max_level:1)){
    hrus<-unique(routpar$TO[routpar$LEVEL==level])
    #print(paste0("There are ",length(hrus)," hrus in level ",level))

    if(length(hrus)>100) {
      flowaccu<-mclapply(hrus,hru_accm,water=water,routpar=routpar,mc.cores = mc_cores)
      #print("using paralell")
    }else{
      flowaccu<-lapply(hrus,hru_accm,water=water,routpar=routpar)
    }
    for (i in c(1:length(hrus))) water$flow[water$HUC12==hrus[i]]<- flowaccu[[i]]
  }
  return(water)
}


# A SUN_ET function old----
#' A SUN_ET function
#'
#' Calculate PET of each pixel/watershed with monthly temperature data by using Hamon's PET formular and then calculate SUN_ET with monthly ET~f(P,LAI,PET)
#' @param data_monthly this is the monthly dataframe data, which includes all main variabels (P,T,LAI, etc.)
#' @param pars A names vector inlcudes "ALT","LAT","LONG",and "VEG" infomation
#' @export
#' @examples
#' PET<-f_ET_SUN(data_monthly_frame,pars)

f_ET_SUN<-function(data_monthly,pars){
  ## Hamon's PET ----
  # calculate PET in monthly with T
  f_Hamon<-function (par=1.2, tavg, lat, jday)
  {
    var_theta <- 0.2163108 + 2 * atan(0.9671396 * tan(0.0086 *
                                                        (jday - 186)))
    var_pi <- asin(0.39795 * cos(var_theta))
    daylighthr <- 24 - 24/pi * acos((sin(0.8333 * pi/180) + sin(lat *
                                                                  pi/180) * sin(var_pi))/(cos(lat * pi/180) * cos(var_pi)))
    esat <- 0.611 * exp(17.27 * tavg/(237.3 + tavg))
    return(par * 29.8 * daylighthr * (esat/(tavg + 273.2)))
  }
  # monthly time seqeue
  data_daily<-month2daily(data_monthly,ts = T)
  data_daily$Jday<-as.numeric(format(index(data_daily), "%j"))
  #calculate PET with Hamon's PET formular
  data_daily$PET<-f_Hamon(tavg=data_daily$T,lat=pars["LAT"],jday = data_daily$Jday)
  data_mon<-daily2monthly(data_daily,FUN = sum)
  data_monthly$PET <- data_mon$PET

  # calculate ET based on SUN GE's ET and LAI, PET relationship
  data_monthly$E <- -4.79 + 0.75*data_monthly$PET + 3.29*data_monthly$LAI+0.04*data_monthly$P

  data_monthly
}



#' A f_WaSSI function
#'
#' Calculate Q with SMA-SMC and snow model
#' @param data_in this is the monthly dataframe data, which includes all main variabels (P, T, E, Q etc.)
#' @param pars A names vector inlcudes "ALT","LAT","LONG",and "VEG" infomation
#' @param soil_in A names vector inlcudes all soil parameters
#' @param calibrate Default =NA, defines wheter calibrate model in watershed scale
#' @param daily A logic object defined the scale of input data "default is daily=F monthly"
#' @param scale Simulation scale c("month","ann")
#' @param y_s,y_e The start and end year for simulation
#' @export
#' @examples
#' WaS<-f_WaSSI(data_in,pars,soil_in,calibrate=F,daily=F,y_s=NA,y_e=NA,scale="month")
#'

f_WaSSI<-function(data_in,
                  pars,
                  soil_in,
                  calibrate=F,
                  daily=F,
                  y_s=NA,
                  y_e=NA,
                  scale="month"){

  require(hydromad)
  require(hydroTSM)
  require(hydroGOF)
  require(SPEI)
  require(zoo)
  ##############################################################################
  # Prepare input data
  ##############################################################################
  # if the input is a zoo
  if(is.zoo(data_in)){
    y_start<- as.Date(index(data_in)[1])
    y_end<-as.Date(index(data_in)[length(data_in[,1])])
    HydroTestData<-data_in
    data_in_all<-data_in
  }else{

    ## creat the target data zoo
    # defined the time range base on the input data frame
    y_start<-paste(min(data_in$YEAR),"-01-01",sep="")
    y_end<-paste(max(data_in$YEAR),"-12-31",sep="")

    # daily time seqeue
    timeseq <- seq(as.POSIXct(y_start, tz = "GMT"),as.POSIXct(y_end, tz = "GMT"),by = "day")

    # if input is daily data
    if(daily){
      #print("Daily simulation")
      data_daily<-data_in

    }else{

      # monthly time seqeue
      timeseq_month <- seq(as.POSIXct(y_start, tz = "GMT"),as.POSIXct(y_end, tz = "GMT"),by = "month")

      ## get the input monthly data
      Date_monthly<-data.frame(YEAR=substr(as.character(timeseq_month),1,4),MONTH=substr(as.character(timeseq_month),6,7))

      # get the default P, T data
      data_monthly<-data.frame(Date_monthly,P=data_in$P,T=data_in$T)
      data_in_all<-as.zooreg(zoo(data_in[,c(-1,-2)], order.by = timeseq_month))
      # get the Q validation data if exit
      if (length(which(names(data_in) %in% "Q"))>0){data_monthly$Q<-data_in$Q}

      # get the PET data if exit otherwise calculate it by "thornthwaite" with T and latitude
      if (length(which(names(data_in) %in% "E"))>0){data_monthly$E<-data_in$E} else{data_monthly$E<-as.vector(thornthwaite(data_in$T,lat = pars["LAT"]))}

      # get Date info for all daily data
      DATE_daily<-data.frame(YEAR=substr(as.character(timeseq),1,4),MONTH=substr(as.character(timeseq),6,7),DAY=substr(as.character(timeseq),9,10))

      # get the number of days for each month
      numberOfDays <- function(date) {
        m <- format(date, format="%m")

        while (format(date, format="%m") == m) {
          date <- date + 1
        }

        return(as.integer(format(date - 1, format="%d")))
      }

      days<-sapply(timeseq_month,function(x) numberOfDays(as.Date(x)))

      # averaged the monthly summed data and join it to daily Date
      ## get mean daily value
      # index needed to averaged
      index_ave<-which(!names(data_monthly) %in% c("ID","YEAR","MONTH","Tmin","T","Tmax","Tavg_C"))

      data_monthly[,index_ave]<-data_monthly[,index_ave]/days

      data_daily<-join(DATE_daily,data_monthly,by=c("YEAR","MONTH"))
    }

    HydroTestData <- as.zooreg(zoo(data_daily[,c(-1,-2,-3)], order.by = timeseq))
  }
  head(HydroTestData)

  ##############################################################################
  # Simulate snow----
  ##############################################################################
  ##  Simulate snow

  # print default parameters for snowmelt calculation
  #print(str(hydromad.options("snow")))

  # Use T mean instead of E to calculate snow
  names(HydroTestData)[c(which(names(HydroTestData)=="E"),which(names(HydroTestData)=="T"))]<-c("E0","E")

  #data(HydroTestData)
  mod0 <- hydromad(HydroTestData, sma = "snow", routing = "expuh")
  ## simulate with some arbitrary parameter values
  mod1 <- update(mod0, Tmax = 1, Tmin = -1, cr = 1, cs = 1,
                 kd = 3, kf = 1, rcap = 0.5,
                 d = 200, f = 0.5, e = 0.1, tau_s = 10)

  ## plot results with state variables
  snow_out <- predict(mod1, return_state = T)
  #xyplot(cbind(HydroTestData[,1:2], snow = snow_out[,c(1,2,4)]))

  # give effective rainfall to origal data
  HydroTestData[,"P"]<-snow_out[,"U"]
  names(HydroTestData)[c(which(names(HydroTestData)=="E0"),which(names(HydroTestData)=="E"))]<-c("E","T")

  ##############################################################################
  # Simulate Q and AET based on SMA-SMC
  ##############################################################################

  # calculate Q and actual ET based on SMA-SAC

  # define simulation time period
  if(is.na(y_s)){y_start_sim<-y_start}else{y_start_sim<-paste(y_s,"-01-01",sep="")}
  if(is.na(y_e)){y_end_sim<-y_end}else{y_end_sim<-paste(y_e,"-12-31",sep="")}

  # select particular time period for calibration and validation
  index_sim<-which(index(HydroTestData)<= as.POSIXct(y_end_sim) & index(HydroTestData)>= as.POSIXct(y_start_sim) )

  ## if Q exist ### Fit SAM-SAC model soil parameters based on “Hydromad” fiting function in Watershed scale with Q validation data
  if (calibrate & length(which(names(HydroTestData)=="Q"))>0){

    print("calibration")
    ## an unfitted model, with ranges of possible parameter values
    modx <- hydromad(HydroTestData[index_sim,], sma = "sacramento")

    ## now try to fit it
    fitx <- fitByOptim(modx)

    # print the summary info for the fitted model
    print(summary(fitx))
    print("NSE=")
    print(NSE(fitx$data[,"Q"],fitx$U[,"U"]))
    fitx$fit.result$NSE=NSE(fitx$data[,"Q"],fitx$U[,"U"])
    # plot fitted result with P
    # pdf(paste("images/Calibration_",name,".pdf",sep=""))
    # xyplot(fitx, with.P = TRUE, type = c("l", "g"))
    # dev.off()
    # transfer soil input
    soil_pars<-list("uztwm" = 1, "uzfwm" = 150, "uzk" = 0.1, "pctim" = 1e-06, "adimp" = 0,
                    "zperc" = 28.8997, "rexp" = 5, "lztwm" = 205.652, "lzfsm" = 758.774,
                    "lzfpm" = 1000, "lzsk" = 0.149213, "lzpk" = 0.00607691, "pfree" = 0.582714)
    soil_pars[c("uztwm" ,"uzfwm","uzk","zperc", "rexp", "lztwm", "lzfsm","lzfpm", "lzsk", "lzpk", "pfree")]<-fitx$parlist[c("uztwm" ,"uzfwm","uzk","zperc", "rexp", "lztwm", "lzfsm","lzfpm", "lzsk", "lzpk", "pfree")]

    # simulate Q with SMA-SAC model get all soil moisture data
    out<-sacramento.sim(HydroTestData[index_sim,],uztwm = soil_pars["uztwm"], uzfwm = soil_pars["uzfwm"], uzk = soil_pars["uzk"], pctim = soil_pars["pctim"], adimp = soil_pars["adimp"],
                        zperc = soil_pars["zperc"], rexp = soil_pars["rexp"], lztwm = soil_pars["lztwm"], lzfsm = soil_pars["lzfsm"],
                        lzfpm = soil_pars["lzfpm"], lzsk = soil_pars["lzsk"], lzpk = soil_pars["lzpk"], pfree = soil_pars["pfree"],return_state = T)

    # print the summary info for the fitted model
    #print(summary(out))
    # print(head(out) )


  }else{
    ### Calculate ET based on SAM-SAC with soil parameters inputs
    print("simulation without calibration")
    # transfer soil input
    soil_pars<-list("uztwm" = 1, "uzfwm" = 150, "uzk" = 0.1, "pctim" = 1e-06, "adimp" = 0,
                    "zperc" = 28.8997, "rexp" = 5, "lztwm" = 205.652, "lzfsm" = 758.774,
                    "lzfpm" = 1000, "lzsk" = 0.149213, "lzpk" = 0.00607691, "pfree" = 0.582714)
    soil_pars[c("uztwm" ,"uzfwm","uzk","zperc", "rexp", "lztwm", "lzfsm","lzfpm", "lzsk", "lzpk", "pfree")]<-soil_in[c("uztwm" ,"uzfwm","uzk","zperc", "rexp", "lztwm", "lzfsm","lzfpm", "lzsk", "lzpk", "pfree")]

    # simulate Q with SMA-SAC model get all soil moisture data
    out<-sacramento.sim(HydroTestData[index_sim,],uztwm = soil_pars["uztwm"], uzfwm = soil_pars["uzfwm"], uzk = soil_pars["uzk"], pctim = soil_pars["pctim"], adimp = soil_pars["adimp"],
                        zperc = soil_pars["zperc"], rexp = soil_pars["rexp"], lztwm = soil_pars["lztwm"], lzfsm = soil_pars["lzfsm"],
                        lzfpm = soil_pars["lzfpm"], lzsk = soil_pars["lzsk"], lzpk = soil_pars["lzpk"], pfree = soil_pars["pfree"],return_state = T)

    # print the summary info for the fitted model
    #print(summary(out))
    # print(head(out) )
  }

  ##############################################################################
  # Merge simulated result
  ##############################################################################

  # merge all data
  snow_out_al<-snow_out[index_sim,]
  out_result<-merge(snow_out_al,out)
  #print(head(out_result) )
  #print(summary(out_result))

  ##############################################################################
  # Tansfer daily data to monthly or annual
  ##############################################################################

  # input data
  if(daily){

    index_in<-which(index(data_in_all)<= as.POSIXct(y_end_sim) & index(data_in_all)>= as.POSIXct(y_start_sim) )
    data_in_all<-data_in_all[index_in,]
    # index_T<-which(names(data_in_all) %in% c("T","Tmin","Tmax"))
    # index_NT<-which(! names(data_in_all) %in% c("T","Tmin","Tmax"))
    #
    # Water_mon<-daily2monthly(data_in_all[,index_NT],FUN=sum,na.rm = T)
    # Temp_mon<-daily2monthly(data_in_all[,index_T],FUN=mean,na.rm = T)
    # HydroTestData_mon<-cbind(Water_mon,Temp_mon)
    result_daily<-merge(data_in_all,out_result)
    result_mon<-daily2monthly(result_daily,FUN = sum,na.rm = T,out.fmt = "%Y-%m-%d")

  }else{

    result_mon<-daily2monthly(out_result,FUN = sum,na.rm = T,out.fmt = "%Y-%m-%d")
    #timeseq_month
    index_in<-which(index(data_in_all)<= as.POSIXct(y_end_sim) & index(data_in_all)>= as.POSIXct(y_start_sim) )
    data_in_all<-data_in_all[index_in,]
    # index_T<-which(names(data_in_all) %in% c("T","Tmin","Tmax"))
    # index_NT<-which(! names(data_in_all) %in% c("T","Tmin","Tmax"))
    #
    # Water_mon<-daily2monthly(data_in_all[,index_NT],FUN=sum,na.rm = T)
    # Temp_mon<-daily2monthly(data_in_all[,index_T],FUN=mean,na.rm = T)
    # HydroTestData_mon<-cbind(Water_mon,Temp_mon)
    index(result_mon)<-index(data_in_all)
    result_mon<-merge(data_in_all,result_mon)
  }

  # calculate GEP based on AET and WUE
  GEP<-as.numeric(result_mon[,"AET"])* pars["WUE"]
  # calculate RE based on relationship between RE and ET for vegetation types
  ER<-as.numeric(result_mon[,"AET"])*pars["ER_m"]+pars["ER_n"]
  GEP<-as.zooreg(zoo(GEP, order.by = index(result_mon)))
  ER<-as.zooreg(zoo(ER, order.by = index(result_mon)))
  # merge carbon variables
  result_mon<-merge(result_mon,GEP,ER)

  summary(result_mon)


  ##############################################################################
  # Output result
  ##############################################################################
  if (calibrate & length(which(names(HydroTestData)=="Q"))>0){
    if( scale=="MONTH" | scale=="month" ){

      list("Result"=result_mon,"Fit"=fitx)

    }else if(scale=="ANN" | scale=="ann" | scale=="annual" | scale=="ANNUAL"| scale=="Annual"){

      result_ann<-monthly2annual(result_mon,FUN = sum,na.rm = T)
      list("Result"=result_ann,"Fit"=fitx)
    }else{
      result_daily
      list("Result"=result_daily,"Fit"=fitx)
    }

  }else{

    if( scale=="MONTH" | scale=="month" ){

      result_mon

    }else if(scale=="ANN" | scale=="ann" | scale=="annual" | scale=="ANNUAL"| scale=="Annual"){

      result_ann<-monthly2annual(result_mon,FUN = sum,na.rm = T)
      result_ann
    }else{
      result_daily
    }
  }
}



#' Calculate WaSSI for each pixel function
#'
#' Calculate WaSSI for each pixel function
#' @param lin this is vector data, which includes all main variabels (P, T, E, Q etc.)
#' @param pars A names vector inlcudes "ALT","LAT","LONG",and "VEG" infomation
#' @param S_y,E_y The start and end year for input climate data
#' @param S_y_LAI,E_y_LAI The start and end year for LAI data
#' @param watershed whether calculate in watershed scale
#' @param calibrate whether calibrate first
#' @param y_s,y_e the simulate start and end year
#' @param Q A dataframe of Q data with c["YEAR","Month","Q"]
#' @export
#' @examples
#' WaS<-f_cal_WaSSI(lin,S_y,E_y,S_y_LAI,E_y_LAI,watershed=F,Q=NA,calibrate=NA,y_s=NA,y_e=NA)
#'
f_cal_WaSSI<-function(lin,S_y,E_y,S_y_LAI,E_y_LAI,watershed=F,Q=NA,calibrate=NA,y_s=NA,y_e=NA){

  Year_C<-rep(c(S_y:E_y), each=12)
  Year_LAI<-c(S_y_LAI:E_y_LAI)
  Month_C<-rep(c(1:12),E_y-S_y+1)
  Date.frame<-data.frame(YEAR=Year_C,Month=Month_C)

  data_monthly_frame<-data.frame(Date.frame,P=lin[1:length(Month_C)],T=lin[(length(Month_C)+1):(2*length(Month_C))],LAI=NA)

  # in watershed scale merge Q data by date
  if( watershed){ data_monthly_frame<-merge(data_monthly_frame,Q,by=c("YEAR","Month"),all.x=T) }

  data_monthly_frame$LAI[(Date.frame$YEAR %in% Year_LAI)]<-lin[(2*length(Month_C)+1):(2*length(Month_C)+12*length(Year_LAI))]

    # give the oldest LAI to the history
  if (S_y_LAI > S_y){
    data_monthly_frame$LAI[(Date.frame$YEAR < S_y_LAI)]<-data_monthly_frame$LAI[(Date.frame$YEAR == S_y_LAI)]
  }

  # give the newest LAI to the future
  if (E_y_LAI > E_y){
    data_monthly_frame$LAI[(Date.frame$YEAR > E_y_LAI)]<-data_monthly_frame$LAI[(Date.frame$YEAR == E_y_LAI)]
  }

  # read soil pars
  soil_in<-list("uztwm" = 1, "uzfwm" = 150, "uzk" = 0.1, "zperc" = 28.8997, "rexp" = 5, "lztwm" = 205.652, "lzfsm" = 758.774,
                "lzfpm" = 1000, "lzsk" = 0.149213, "lzpk" = 0.00607691, "pfree" = 0.582714)

  soil_in[1:11]<-lin[(length(lin)-16):(length(lin)-6)]

  pars<-c("ALT"=lin[(length(lin)-5)],"LAT"=lin[(length(lin)-4)],"LONG"=lin[(length(lin)-3)],"WUE"=lin[(length(lin)-2)],"ER_m"=lin[(length(lin)-1)],"ER_n"=lin[(length(lin)-0)])

  # calculate SUN-ET
  data_monthly_frame<-f_ET_SUN(data_monthly_frame,pars=pars)

  result<-f_WaSSI(data_monthly_frame,pars,soil_in,calibrate=calibrate,y_s=y_s,y_e=y_e,scale="month")
  result
  #print(head(result))

}

