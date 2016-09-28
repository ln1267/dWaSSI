# WaSSI-Model functions (f_SUN_ET,f_WaSSI,f_cal_WaSSI)

#' A SUN_ET function
#'
#' Calculate PET of each pixel/watershed with monthly temperature data by using Hamon's PET formular and then calculate SUN_ET with monthly ET~f(P,LAI,PET)
#' @param data_monthly this is the monthly dataframe data, which includes all main variabels (P,T,LAI, etc.)
#' @param pars A names vector inlcudes "ALT","LAT","LONG",and "VEG" infomation


# calculate PET in monthly with T
f_ET_SUN<-function(data_monthly_frame,pars){
  y_start<-paste(min(data_monthly_frame$YEAR),"-01-01",sep="")
  y_end<-paste(max(data_monthly_frame$YEAR),"-12-31",sep="")
  # monthly time seqeue
  timeseq_mon <- seq(as.POSIXct(y_start, tz = "GMT"),as.POSIXct(y_end, tz = "GMT"),by = "month")

  days<-sapply(timeseq_mon,function(x) numberOfDays(as.Date(x)))

  #calculate PET with Hamon's PET formular
  MJD<-rep(c(15,46,76,107,137,168,198,229,259,290,321,351),length(days)/12)

  data_monthly_frame$PET <- days*0.1651*1.2*216.7*6.108*exp(17.27*data_monthly_frame$T/(data_monthly_frame$T+237.3))/(data_monthly_frame$T+237.3)*2*acos(-1*tan(pars["LAT"]*.0174529)*tan(0.4093*sin(2*3.1415*MJD/365.-1.405)))/3.1415

  # calculate ET based on SUN GE's ET and LAI, PET relationship
  data_monthly_frame$E <- -4.79 + 0.75*data_monthly_frame$PET + 3.29*data_monthly_frame$LAI+0.04*data_monthly_frame$P

  data_monthly_frame
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

f_WaSSI<-function(data_in,pars,soil_in,calibrate=NA,daily=F,y_s=NA,y_e=NA,scale="month"){

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
  head(HydroTestData)

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

  # calculate Q and actual ET based on SMA-SAC

  # define simulation time period
  if(is.na(y_s)){y_start_sim<-y_start}else{y_start_sim<-paste(y_s,"-01-01",sep="")}
  if(is.na(y_e)){y_end_sim<-y_end}else{y_end_sim<-paste(y_e,"-12-31",sep="")}

  # select particular time period for calibration and validation
  index_sim<-which(index(HydroTestData)<= as.POSIXct(y_end_sim) & index(HydroTestData)>= as.POSIXct(y_start_sim) )

  ## if Q exist ### Fit SAM-SAC model soil parameters based on “Hydromad” fiting function in Watershed scale with Q validation data
  if (! is.na(calibrate) & length(which(names(HydroTestData)=="Q"))>0){

    print("calibration")
    ## an unfitted model, with ranges of possible parameter values
    modx <- hydromad(HydroTestData[index_sim,], sma = "sacramento")

    ## now try to fit it
    fitx <- fitByOptim(modx)

    # print the summary info for the fitted model
    print(summary(fitx))
    print("NSE=")
    print(NSE(fitx$data[,"Q"],fitx$U[,"U"]))

    # plot fitted result with P
    #xyplot(fitx$parlist, with.P = TRUE, type = c("l", "g"))

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

  # merge all data
  snow_out_al<-snow_out[index_sim,]
  out_result<-merge(snow_out_al,out)
  #print(head(out_result) )
  #print(summary(out_result))

  # output result
  if( scale=="MONTH" | scale=="month" ){
    result_mon<-daily2monthly(out_result,FUN = sum,na.rm = T,out.fmt = "%Y-%m-%d")

    # input data
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
    summary(result_mon)
    head(result_mon)
    result_mon
  }else if(scale=="ANN" | scale=="ann"){
    result_ann<-daily2annual(out_result,FUN = sum,na.rm = T)
    summary(result_ann)
    head(result_ann)
    result_ann
  }else{
    out_result
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

f_cal_WaSSI<-function(lin,S_y,E_y,S_y_LAI,E_y_LAI,watershed=F,Q=NA,calibrate=NA,y_s=NA,y_e=NA){

  Year_C<-rep(c(S_y:E_y), each=12)
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

  soil_in[1:11]<-lin[(length(lin)-14):(length(lin)-4)]

  pars<-c("ALT"=lin[(length(lin)-3)],"LAT"=lin[(length(lin)-2)],"LONG"=lin[(length(lin)-1)],"VEG"=lin[(length(lin)-0)])

  # calculate SUN-ET
  data_monthly_frame<-f_ET_SUN(data_monthly_frame,pars=pars)

  result<-f_WaSSI(data_monthly_frame,pars,soil_in,calibrate=calibrate,y_s=y_s,y_e=y_e,scale="month")
  result
  #print(head(result))

}

