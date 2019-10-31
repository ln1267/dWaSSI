# This is the global Rscript for this app

Ning<-"Ning Liu"
testfile<-"tmp/test.csv"
# Load all libraries----
library("dplyr")
library("plyr")
library("zoo")
# library("hydromad")
# library("hydroGOF")
library("xts")

uploadfilenames<-"test.csv"

if(file.exists("www/outp.RData")){
  load("www/outp.RData")
  result_vars<-names(outp$Out)
  input_vars<-names(outp$Input)
}
# Read csv and add timestamp to the dataframe based on "Year" and "Month" columns ----
f_read_csv<-function(filename,filesep){

  df <- read.csv(
    paste0("tmp/",filename),
    header = T,
    sep = filesep
  )
  times<-sum(c("Year","Month","Day") %in% names(df))
  if(times==1){
    df$Timestamp<-as.Date(paste0(df$Year,"0101"),format="%Y%m%d")
  }else if(times==2){
    df$Timestamp<-as.Date(paste0(df$Year,substr(df$Month+100,2,3),"01"),format="%Y%m%d")
  }else if(times==3){
    df$Timestamp<-as.Date(paste0(df$Year,substr(df$Month+100,2,3),substr(df$Day+100,2,3)),format="%Y%m%d")
  }
  return(df)
}

# Subset the input data for simulation ----
f_subset<-function(station,daterange,filename){
  # Select station and dateperiod
  a<-lapply(c("climate","lai","lucc","soil","flow"),function(x) regexpr(x,tolower(filename)))

  if(!file.exists(paste0("tmp/",filename))) return(print(paste0("There is no uploaded file named as ",filename)))
  test<-f_read_csv(filename,filesep = ",")
  if(!station %in% unique(test$Station)) return(print(paste0("There is no Station ID named as ",station)))
  if(a[[1]][1]>0){
    Climate<-f_read_csv(filename,filesep = ",")
    if(range(Climate$Timestamp)[1]>daterange[1]) return(print(paste0("The start date should be later than ",range(Climate$Timestamp)[1])))
    if(range(Climate$Timestamp)[2]<daterange[2]) return(print(paste0("The end date should be earler than ",range(Climate$Timestamp)[2])))
    Climate<-subset(Climate,Station == station & Timestamp >=daterange[1] & Timestamp <=daterange[2])
    write.csv(Climate[-1],"./www/Climate.csv",row.names = F)
  }else if(a[[2]][1]>0){
    LAI<-f_read_csv(filename,filesep = ",")
    if(range(LAI$Timestamp)[1]>daterange[1]) return(print(paste0("The start date should be later than ",range(LAI$Timestamp)[1])))
    if(range(LAI$Timestamp)[2]<daterange[2]) return(print(paste0("The end date should be earler than ",range(LAI$Timestamp)[2])))
    LAI<-subset(LAI,Station == station & Timestamp >=daterange[1] & Timestamp <=daterange[2])
    write.csv(LAI[-1],"./www/LAI.csv",row.names = F)
  }else if(a[[3]][1]>0){
    LUCC<-f_read_csv(filename,filesep = ",")
    LUCC<-subset(LUCC,Station == station)
    write.csv(LUCC[-1],"./www/LUCC.csv",row.names = F)
  }else if(a[[4]][1]>0){
    SOIL<-f_read_csv(filename,filesep = ",")
    SOIL<-subset(SOIL,Station == station)
    write.csv(SOIL[-1],"./www/SOIL.csv",row.names = F)
  }else if(a[[5]][1]>0){
    Flow<-f_read_csv(filename,filesep = ",")
    if(range(Flow$Timestamp)[1]>daterange[1]) return(print(paste0("The start date should be later than ",range(Flow$Timestamp)[1])))
    if(range(Flow$Timestamp)[2]<daterange[2]) return(print(paste0("The end date should be earler than ",range(Flow$Timestamp)[2])))
    Flow<-subset(Flow,Station == station & Timestamp >=daterange[1] & Timestamp <=daterange[2])
    Flow<-subset(Flow,Station == station)
    write.csv(Flow[-1],"./www/Flow.csv",row.names = F)
  }
  print("Finished subsetting data and saved to www folder")
}


# Merge all input data ----

f_merge_input<-function(){

  for(inputfile in c("Climate.csv","LAI.csv","SOIL.csv","LUCC.csv")){
  if(!file.exists(paste0("www/",inputfile))) return(print(paste0("There is no ", input$inputfiles," input")))
  }

  if(!file.exists("www/Flow.csv")) print("There is no flow input")

  climate<-read.csv("www/Climate.csv",sep=",",stringsAsFactors = F)
  lai<-read.csv("www/LAI.csv",sep=",",stringsAsFactors = F)

  data_TS<-merge(climate,lai,by=c("Timestamp"),all.x=T)
  data_TS<-data_TS[c("Timestamp","P","T","LAI")]

  if(file.exists("www/Flow.csv")){
    Streamflow<-read.csv("www/Flow.csv",sep=",",stringsAsFactors = F)
    data_TS<-merge(data_TS,Streamflow,by=c("Timestamp"),all.x=T)
    data_TS<-data_TS[c("Timestamp","P","T","LAI","Q")]
  }

  data_TS$Timestamp<-as.Date(data_TS$Timestamp,format="%Y-%m-%d")
  data_TS[2:length(data_TS)]<-round(data_TS[2:length(data_TS)],2)
  save(data_TS,file="www/data_TS.Rdata")
  print("Summary informaiton")
  print(summary(data_TS))

}

# function for WaSSI-C----
## Calculate wue ----
f_WaSSI_WUE<-function(veg_ratio){
  WUEs<-data.frame(
    IGBP=c("ENF","EBF","DNF","DBF","MF","CSH","OSH","WSA","SAV","GRA","WET","CRO"),
    WUE_IGBP_SUN=c(2.46,2.59,2.46,3.2,2.74,1.37,1.33,1.26,1.26,2.12,0,3.13),  # from SUN 2010
    WUE_IGBP_ZH=c(2.710,2.675,2.71,3.280,3.349,1.734,1.679,1.689,2.032,2.333,2.595,2.533), # from \cite{Zhang2016d}
    ER_IGBP_SUN_n=c(9.9,19.6,9.9,30.8,24.4,11.4,9.7,14.7,25.2,18.9,7.8,40.6),
    ER_IGBP_SUN_m=c(0.68,0.61,0.68,0.45,0.62,0.69,0.56,0.63,0.53,0.64,0.56,0.43),
    WUE_IGBP_OzFlux=c(1.07,3.18,2.71,1.79,3.349,1.734,1.66,1.77,1.52,1.79,0.72,1.14)
  )

  veg_ratio<-merge(veg_ratio,WUEs,by="IGBP",all.x=T)
  .ER_n<-sum(veg_ratio$ER_IGBP_SUN_n*veg_ratio$Ratio)
  .ER_m<-sum(veg_ratio$ER_IGBP_SUN_m*veg_ratio$Ratio)
  .WUE<-sum(veg_ratio$WUE_IGBP_ZH*veg_ratio$Ratio)
  return(c("WUE"=.WUE,"ER_m"=.ER_m,"ER_n"=.ER_n))
}

# veg_ind$IGBP<-WUEs$IGBP[veg_ind$Levels]
# f_WaSSI_WUE(veg_ratio =veg_ind )

## This function allows you to get the number of days for a specific month.----
numberOfDays<-function(date) {
  m <- format(date, format="%m")

  while (format(date, format="%m") == m) {
    date <- date + 1
  }

  return(as.integer(format(date - 1, format="%d")))
}


daily2monthly<-function (x, FUN, na.rm = TRUE, ...)
{
  library(xts)
  library(zoo)
  if (missing(FUN))
    stop("Missing argument value: 'FUN' must contain a valid function for aggregating the values")
  if (sfreq(x) %in% c("monthly", "quarterly", "annual"))
    stop("Invalid argument: 'x' is not a (sub)daily/weekly ts. 'x' is a ",
         sfreq(x), " ts")
  dates <- time(x)
  months <- as.Date(as.yearmon(time(x)))
  tmp <- aggregate(x, by = months, FUN, na.rm = na.rm)
  nan.index <- which(is.nan(tmp))
  if (length(nan.index) > 0)
    tmp[nan.index] <- NA
  inf.index <- which(is.infinite(tmp))
  if (length(inf.index) > 0)
    tmp[inf.index] <- NA
  if (NCOL(tmp) == 1)
    tmp <- zoo(as.numeric(tmp), time(tmp))
  return(tmp)
}

sfreq <- function(x, min.year=1800) {

  # Checking that 'class(x)'
  valid.class <- c("xts", "zoo")
  if (length(which(!is.na(match(class(x), valid.class )))) <= 0)
    stop("Invalid argument: 'x' must be in c('xts', 'zoo')" )

  out <- periodicity(x)$scale # xts::periodicity

  if (out == "yearly") out <- "annual"

  return(out)

} # 'sfreq' END

## Hamon's PET ----
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

## Calculate ET based on SUN ----
f_SUN_ET<-function(data_monthly,pars){
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

## month2daily for converting monthly data to daily----
month2daily<-function(data_in,ts=F,fields=NA){

  data_in$Year<-as.integer(format(data_in$Timestamp,"%Y"))
  data_in$Month<-as.integer(format(data_in$Timestamp,"%m"))
  data_in<-data_in[-1]
  ## creat the target data zoo
  # defined the time range base on the input data frame
  y_start <- paste(min(data_in$Year), "-01-01", sep = "")
  y_end <- paste(max(data_in$Year), "-12-31", sep = "")

  # daily time seqeue
  timeseq <- seq(as.POSIXct(y_start), as.POSIXct(y_end), by = "day")

  # monthly time seqeue
  timeseq_month <-seq(as.POSIXct(y_start), as.POSIXct(y_end), by = "month")

  ## get the input monthly data
  Date_monthly <-  data.frame(Year = as.integer(format(timeseq_month,"%Y")),
                              Month = as.integer(format(timeseq_month,"%m")))

  # get the default P, T data
  data_monthly <- merge(Date_monthly,data_in,by=c("Year","Month"),all.x=T)
  data_monthly<-arrange(data_monthly,Year,Month)

  # get Date info for all daily data
  DATE_daily <- data.frame(Year = as.integer(format(timeseq,"%Y")),
                           Month = as.integer(format(timeseq,"%m")),
                           Day=as.integer(format(timeseq,"%d")))

  # Averge the monthly data to daily
  if(!is.na(fields)){
    days <-sapply(timeseq_month, function(x)
      numberOfDays(as.Date(x)))
    # monthly data to daily data
    col_index<-which(names(data_monthly) %in% fields)
    data_monthly[col_index]<- data_monthly[col_index]/days
  }

  data_daily <-
    merge(DATE_daily,
          data_monthly,
          by = c("Year", "Month"),
          all.x = T)
  data_daily<-arrange(data_daily,Year,Month,Day)
  data_daily[-c(1:3)]<-round(data_daily[-c(1:3)],2)
  if(ts)  data_daily<-as.zooreg(zoo(data_daily[,c(-1,-2,-3)], order.by = timeseq))
  print(head(data_daily))
  data_daily
}
#data_TS_daily<-month2daily(ET,fields = c("P","Q","PET","E"),ts = T)


## simulated snow----
f_snow_calculation <- function(HydroData) {
  # Simulate snow to get the effective rainfall using "p" and "Tmean"
  require(hydromad)
  # print default parameters for snowmelt calculation
  print(str(hydromad.options("snow")))

  # Use T mean instead of E to calculate snow
  if(which(names(HydroData) =="E")>0) names(HydroData)[which(names(HydroData) =="E")] <- c("E0")
  names(HydroData)[which(names(HydroData) =="T")] <- c("E")

  mod0 <- hydromad(HydroData, sma = "snow", routing = "expuh")
  ## simulate with some arbitrary parameter values
  mod1 <- update(mod0,Tmax = 1,Tmin = -1,cr = 1,cs = 1,kd = 3,kf = 1,rcap = 0.5,d = 200,f = 0.5,e = 0.1,tau_s = 10) #

  ## plot results with state variables
  snow_out <- predict(mod1, return_state = T)
  #xyplot(cbind(HydroTestData[,1:2], snow = snow_out[,c(1,2,4)]))

  # give effective rainfall to origal data
  HydroData$P_original<-HydroData$P
  HydroData$P <- snow_out$U
  names(HydroData)[which(names(HydroData) =="E")] <- c("T")
  if(which(names(HydroData) =="E0")>0) names(HydroData)[which(names(HydroData) =="E0")] <- c("E")
  print(head(HydroData))
  HydroData
}

#snow_out<-f_snow_calculation(data_TS_daily)

## Simulate streamflow ----
f_WaSSI_Q<-function(HydroData,soil_in=NA,y_s=NA,y_e=NA,calibrate=F){
  # calculate Q and actual ET based on SMA-SAC

  # define simulation time period
  if(is.na(y_s)){y_start_sim<-y_start}else{y_start_sim<-paste(y_s,"-01-01",sep="")}
  if(is.na(y_e)){y_end_sim<-y_end}else{y_end_sim<-paste(y_e,"-12-31",sep="")}

  # select particular time period for calibration and validation
  index_sim<-which(index(HydroData)<= as.POSIXct(y_end_sim) & index(HydroData)>= as.POSIXct(y_start_sim) )

  ## if Q exist ### Fit SAM-SAC model soil parameters based on “Hydromad” fiting function in Watershed scale with Q validation data
  if (calibrate & length(which(names(HydroData)=="Q"))>0){

    print("calibration")
    ## an unfitted model, with ranges of possible parameter values
    modx <- hydromad(HydroData[index_sim,], sma = "sacramento")

    ## now try to fit it
    fitx <- fitByOptim(modx)

    # print the summary info for the fitted model
    print(summary(fitx))
    # print("NSE=")
    #print(NSE(fitx$data[,"Q"],fitx$U[,"U"]))
    #fitx$fit.result$NSE=NSE(fitx$data[,"Q"],fitx$U[,"U"])
    # plot fitted result with P
    cc<- hydromad::xyplot(fitx, with.P = T, type = c("l", "g"))
    print(cc)
    # transfer soil input
    soil_pars<-list("uztwm" = 1, "uzfwm" = 150, "uzk" = 0.1, "pctim" = 1e-06, "adimp" = 0,
                    "zperc" = 28.8997, "rexp" = 5, "lztwm" = 205.652, "lzfsm" = 758.774,
                    "lzfpm" = 1000, "lzsk" = 0.149213, "lzpk" = 0.00607691, "pfree" = 0.582714)
    soil_pars[c("uztwm" ,"uzfwm","uzk","zperc", "rexp", "lztwm", "lzfsm","lzfpm", "lzsk", "lzpk", "pfree")]<-fitx$parlist[c("uztwm" ,"uzfwm","uzk","zperc", "rexp", "lztwm", "lzfsm","lzfpm", "lzsk", "lzpk", "pfree")]

    # simulate Q with SMA-SAC model get all soil moisture data
    out<-sacramento.sim(HydroData[index_sim,],uztwm = soil_pars["uztwm"], uzfwm = soil_pars["uzfwm"], uzk = soil_pars["uzk"], pctim = soil_pars["pctim"], adimp = soil_pars["adimp"],
                        zperc = soil_pars["zperc"], rexp = soil_pars["rexp"], lztwm = soil_pars["lztwm"], lzfsm = soil_pars["lzfsm"],
                        lzfpm = soil_pars["lzfpm"], lzsk = soil_pars["lzsk"], lzpk = soil_pars["lzpk"], pfree = soil_pars["pfree"],return_state = T)

    # print the summary info for the fitted model
    print(summary(out))
    return(list("Input"=HydroData[index_sim,],"Out"=out,"Fitx"=fitx,"Index_sim"=index_sim,"Start_year"=y_start_sim,"End_year"=y_end_sim))
  }else{
    ### Calculate ET based on SAM-SAC with soil parameters inputs
    print("simulation without calibration")
    # transfer soil input
    soil_pars<-list("uztwm" = 1, "uzfwm" = 150, "uzk" = 0.1, "pctim" = 1e-06, "adimp" = 0,
                    "zperc" = 28.8997, "rexp" = 5, "lztwm" = 205.652, "lzfsm" = 758.774,
                    "lzfpm" = 1000, "lzsk" = 0.149213, "lzpk" = 0.00607691, "pfree" = 0.582714)
    soil_pars[c("uztwm" ,"uzfwm","uzk","zperc", "rexp", "lztwm", "lzfsm","lzfpm","lzsk", "lzpk", "pfree")]<-soil_in[c("uztwm" ,"uzfwm","uzk","zperc", "rexp", "lztwm", "lzfsm","lzfpm", "lzsk", "lzpk", "pfree")]

    # simulate Q with SMA-SAC model get all soil moisture data
    out<-sacramento.sim(HydroData[index_sim,],uztwm = soil_pars["uztwm"], uzfwm = soil_pars["uzfwm"], uzk = soil_pars["uzk"], pctim = soil_pars["pctim"], adimp = soil_pars["adimp"],
                        zperc = soil_pars["zperc"], rexp = soil_pars["rexp"], lztwm = soil_pars["lztwm"], lzfsm = soil_pars["lzfsm"],
                        lzfpm = soil_pars["lzfpm"], lzsk = soil_pars["lzsk"], lzpk = soil_pars["lzpk"], pfree = soil_pars["pfree"],return_state = T)

    # print the summary info for the fitted model
    print(summary(out))
    #print(head(out) )
    return(list("Input"=HydroData[index_sim,],"Out"=out,"Index_sim"=index_sim,"Start_year"=y_start_sim,"End_year"=y_end_sim))
  }
}
#
soil<-read.csv("www/SOIL.csv")
#outp<-f_WaSSI_Q(snow_out,soil,calibrate =T,y_s = 2000,y_e = 2006)
# NSE(outp$Input$Q[13:84],outp$Out$U[13:84])
# coords_frame<-read.csv("www/Coords_basin.csv")
# .alt<-100;.lat<-coords_frame$Lat[coords_frame$Station==station];.long<-coords_frame$Long[coords_frame$Station==station]
# pars<-c("ALT"=.alt,"LAT"=.lat,"LONG"=.long,"WUE"=NA,"ER_m"=NA,"ER_n"=NA)
# pars[c(4:6)]<-f_WaSSI_WUE(LUCC_catchment)
# out_mon<-f_WaSSI_output(outp,pars)
# NSE(out_mon$Annual$Q,out_mon$Annual$U)
# plot(out_mon$Annual$Q,out_mon$Annual$U,lty=2)
# plot(out_mon$Monthly$Q,out_mon$Monthly$U,lty=2)
# summary(lm(out_mon$Monthly$Q~out_mon$Monthly$U))

## output_waSSI and caculate carbon----
f_WaSSI_output<-function(out,pars,snow_out=NA,daily=T,scale="mon"){
  library(hydroTSM)
  # Merge simulated result----
  if(is.na(snow_out)){
    out_result<-out$Out
  }else{
    # merge all data
    snow_out_al<-snow_out[out$Index_sim,]
    out_result<-merge(snow_out_al,out)
  }

  # Tansfer daily data to monthly or annual----

  # merge input data to daily
  result_daily<-merge(outp$Input,out_result)

  # aggregate to monthly
  mean_index<-c("T","Tmin","Tmax","LAI","uztwc","uzfwc","lztwc","lzfsc","lzfpc","adimc")
  index_mean<-which(names(result_daily) %in% mean_index)
  index_sum<-which(!names(result_daily) %in% mean_index)

  Sum_variable_mon<-daily2monthly(result_daily[,index_sum],FUN=sum,na.rm = T)
  Mean_variable_mon<-daily2monthly(result_daily[,index_mean],FUN=mean,na.rm = T)
  result_mon<-merge(Sum_variable_mon,Mean_variable_mon)

  # calculate GEP based on AET and WUE
  GEP<-as.numeric(result_mon[,"AET"])* pars["WUE"]
  # calculate RE based on relationship between RE and ET for vegetation types
  ER<-as.numeric(result_mon[,"AET"])*pars["ER_m"]+pars["ER_n"]
  GEP<-as.zooreg(zoo(GEP, order.by = index(result_mon)))
  ER<-as.zooreg(zoo(ER, order.by = index(result_mon)))
  # merge carbon variables
  result_mon<-merge(result_mon,GEP,ER)

  summary(result_mon)

  ## aggregate to annual
  Sum_variable_ann<-daily2annual(result_daily[,index_sum],FUN=sum,na.rm = T)
  Mean_variable_ann<-daily2annual(result_daily[,index_mean],FUN=mean,na.rm = T)
  result_ann<-merge(Mean_variable_ann,Sum_variable_ann)

  # Output result----

  return(list("Annual"=result_ann,"Monthly"=result_mon,"Daily"=result_daily))
}
