##################################################

# This file includs all functions used for studying the relationship between VEG and PRE in AU

###---------list of functions
# f_lib_check(libs) ##libs is a vector of required libraries
#	f_parallel(data=null, fun=null, type="parallel")		for setting up parallel methods
#f_cor(x,y) ## Calculate the correlationship between two vectors ("x" and "y") and return ("r" and "p")
#	f_m2y(data, fun="mean")			for transfering monthly frame data to annual data by using fun="sum" or "mean"
#	f_summary()					for outputing summary infos of all "data frame objects" in memory
#	f_dp(data,seasonal=TRUE,year_start,year_end)	for seasonal or annual changepoint and MK trend analysis
#	f_plot<-function(data,info,annual=FALSE,monthly=FALSE) for monthly or annual grid data plot
#	f_grid2basin(data,type="annual",fun="mean")	#Transfer grid frame data to basin data by fun="mean"

##################################################


#' A requied libraries Load and check Function
#'
#' This function allows you to check whether the required library has been installed, otherwise it will be installed and load.
#' @param libs A character vector of names of required libraries.
#' @keywords libraries
#' @export
#' @examples
#' libs<-c("ggplot2","caTools")
#' f_lib_check(libs)

# Load and check libraries
f_lib_check<-function(libs){
  for (lib in libs ){
    if(lib %in% rownames(installed.packages())){

    }else{
      install.packages(lib,repos='http://cran.us.r-project.org')
    }
  }

  a<-lapply(libs, require, character.only = TRUE)
}


# TRIM THE ANOMALIES FOR A VARIABLE
### treat +0.5% and -0.5% value as anomaly
##########################################################

cutAnomalies <- function(x){
  # Cut the anomolies
  toPlot <- c(x)
  toPlot1<- toPlot[!is.na(toPlot)]
  toPlot <-toPlot1
  sortedLow<-sort(toPlot)
  lowCut<-sortedLow[length(sortedLow)*0.005]
  sortedHigh<-sort(toPlot, decreasing=T)
  highCut<-sortedHigh[length(sortedHigh)*0.005]
  x[which(x>highCut)]<-highCut
  x[which(x<lowCut)]<-lowCut
  x
}

f_cut<-function(x){

  low<-quantile(x,0.05,na.rm=T)
  high<-quantile(x,0.95,na.rm=T)
  x[x>high]<-high
  x[x<low]<-low
  x
}



# setup parallel in Magnus or Zeus
f_Parallel_set<-function(name="zeus"){
  if (name=="magnus"){
    print("using Magnus for processing")
    cl<<-makeCluster(40, type="FORK", outfile = "parallel_mag.txt")  # set up parallel
    print(mem_used())#detectCores()-1
    print(detectCores())

  }else{
    print("using Zeus for processing")
    cl<<-makeCluster(detectCores()-1, type="FORK", outfile = "parallel_zeus.txt")  # set up parallel
    print(mem_used())#
    print(detectCores())
  }
}

# this plot theme is for ggplot
theme_grid <- function(base_size = 12, base_family = "Times"){
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      #line = element_line(colour="black"),
      #text = element_text(colour="black"),
      axis.title = element_text(size = 14,face="bold"),
      axis.text = element_text(colour="black", size=12),
      #strip.text = element_text(size=12),
      legend.key=element_rect(colour=NA, fill =NA),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black", size=2),
      panel.background = element_rect(fill = "grey70", colour = "black"),
      strip.background = element_rect(fill = NA),
      legend.position="right",
      legend.background = element_blank()

    )
}

# function for MK trend analysis and change points detection ("trend" and "changepoint" packages)
## ts_in is the input time serie; name is the output pdf name; seasonal is wether for seasonal data; plot is whether plot result; main is the title for plot; Y_name is the title for y_axiel; sig is the sig threhold
f_MK_CP<-function(ts_in,name="",seasonal=F,plot=F,main="",Y_name="Streamflow (mm)",sig=0.05){
  require(trend)
  require(changepoint)
  ##  MK and CP analysis
  if (! seasonal){

    #    print("Annual trend analysis")
    #    print(name)
    # changepoint detect
    cp_mean<-cpt.mean(ts_in)
    means<-cp_mean@param.est$mean
    cp_year<-cp_mean@cpts[1]
    ## get changepoint info if it eixts
    if(length(means)>1){

      cp_year<-time(ts_in)[cp_year]
      change<-round((means[2]-means[1])/means[1],3)*100

    }else{

      means<-NA
      cp_year<-NA
      change<-NA
    }

    # MK test
    mk.t<-mk.test(ts_in)
    sen.res <- sens.slope(ts_in)

    # plot for selected stations
    if(plot & sen.res$b.sen !=0 & mk.t$pvalg< sig){
      #print("plot data")
      pdf(paste("result/pdf/",name,".pdf",sep=""),family="Times",width=10)
      t <- (1:(length(ts_in)))
      s.pred <- sen.res$intercept + sen.res$b.sen * t
      s.pred.ts <- ts(s.pred)
      tsp(s.pred.ts) <- tsp(ts_in)
      plot(cp_mean,xlab="Year",ylab=Y_name,main=main)
      lines(s.pred.ts, lty=2)
      dev.off()
    }

    # return
    return(list(CP_M=means,CP_Pec=change,CP_Y=cp_year,MK_P=round(mk.t$pvalg,3),MK_Slope=round(sen.res$b.sen,2)))

  }else{

    print("Seasonal trend analysis")

    # MK test
    mk.t<-smk.test(ts_in)
    #cmk.t<-csmk.test(ts_in)
    sen.res <- sea.sens.slope(ts_in)
    print(names[i])
    print(sen.res)
    print(mk.t)

  }
}


## Calculate annual mean and anomaly (SAI or pecentage change) of a dataset (in array) in parallel (return a list with ("MEAN", "ANOM"))
f_SAI<-function(data=data,method="SAI",mask=NA,plot=F,anom_ab=5){

  #  mask data base on mask
  if(length(mask)>1){
    for (i in 1:dim(data)[3]){
      data[,,i][mask]<-NA
    }
  }

  # Get the mean and sd for the dataset for each pixel
  data_mean<- parApply(cl,data,c(1,2),mean,na.rm=T)
  data_sd<- parApply(cl,data,c(1,2),sd,na.rm=T)

  print("Start anomaly")
  #  annomly
  anom<-array(0,c(nrow(data),ncol(data),dim(data)[3]))
  if(method=="SAI"){
    for (i in 1:dim(data)[3]){anom[,,i]<-(data[,,i]-data_mean)/data_sd}
  }else{
    for (i in 1:dim(data)[3]){anom[,,i]<-(data[,,i]-data_mean)/data_mean*100}
  }

  # abnormal control
  anom[abs(anom)>anom_ab]<-NA
  anom[is.infinite(anom)]<-NA
  print("anom range")
  print(range(anom,na.rm=T))

  # plot annomly of dataset
  if(plot==T){
    require(ggplot2)
    nrows<-nrow(data_mean)
    ncols<-ncol(data_mean)
    a<-raster(data_mean, xmn=112.9, xmx=154, ymn=-43.74, ymx=-8.98)

    # with raster plot or by dataframe
    if(1){
      for (i in 1:dim(anom)[3]){

        nrows<-nrow(anom)
        ncols<-ncol(anom)
        a<-raster(anom[,,i], xmn=112.9, xmx=154, ymn=-43.74, ymx=-8.98)
        plot(a)
        pdf(paste("result/images/anom_",i+1999,".pdf",sep=""))
        plot(a)
        dev.off()
      }
    }else{
      LAT<-rep(seq(-10, by=-0.05, length.out = nrows),ncols)
      LONG<-rep(seq(112.5, by=0.05, length.out = ncols),each=nrows)
      anom_frame<-data.frame(ID=rep(c(1:(nrows*ncols)),15),YEAR=rep(c(2000:2014),each=(nrows*ncols)),LONG=rep(LONG,15),LAT=rep(LAT,15),ANOM=as.vector(anom))

      gplot<-ggplot(aes(x = LONG, y = LAT), data = anom_frame) +
        geom_raster(aes(LONG, LAT, fill=anom_frame[[5]]))+
        facet_wrap( ~ YEAR, ncol=5)+
        scale_fill_gradient(low = 'red', high = 'green',name=names(anom_frame)[5],na.value = "white") +
        coord_equal()+ #xlim=c(102.4, 104.1),ylim = c(30.72,33.18)
        labs(x="Latitude",y="Longitude")+
        theme_grid()

      ggsave(gplot,file =paste("Ann_",names(anom_frame),".pdf",sep="")[5],dpi = 300)

    }
  }
  # return result
  return(list(MEAN=data_mean,ANOM=anom))
}


## Calculate the correlationship between two vectors ("x" and "y") and return ("r" and "p")
f_cor <- function(x,y) {

  res <- try(cor.test(x,y ), silent=TRUE)
  if (class(res)=="try-error") {
    res <- setNames(c(NA, NA), c("estimate","p.value"))

  }else{
    .res<-unclass(res)[c("estimate","p.value")]
    res<-unlist(.res)
  }
  return(res)
}


## setup parallel for "parallel" or "doParallel" or "foreach" or "snow"

f_parallel<-function(data=null, fun=null, type="parallel"){

  if (type=="parallel"){

    library(parallel)
    print("using parallel package to simulate")
    print(paste("Num of Cores=", detectCores()))

    ## set up parallel type to "FORK", all cluster share the global variables
    cl<-makeCluster(detectCores()-1, type="FORK")

  }else if (type=="doParallel"){

    library(doParallel)
    print("using doParallel package to simulate")
    print(paste("Num of Cores=", detectCores()))

    cl<-makeCluster(detectCores()-1, type="FORK")  # set up parallel
    clusterEvalQ(cl, library(rms)) # load required packages "rms"

    print(mem_used())
    #	cl<-makeCluster(detectCores()-1)  # set up parallel
    print(detectCores())
    #	clusterExport(cl,c("x"))    # share default data for all threads

  }else if (type=="foreach"){

    library(foreach)
    library(doParallel)
    print("using foreach package to simulate")

    cl<-makeCluster(detectCores()-1,outfile = "foreach_debug.txt")  # set up parallel
    registerDoParallel(cl)

    foreach(exponent = 2:4,
            .combine = c,
            .export = "base",
            .packages = c("rms", "mice")
    )  %dopar%  {
      tryCatch({
        c(1/x, x, 2^x)
      }, error = function(e) return(paste0("The variable '", x, "'", " caused the error: '", e, "'")))
    }

    stopCluster(cl)

  }else if (type=="snow"){



  }else if (type=="snow"){

    print("using snow package to simulate")

    lnxOptions <-list(host = "itasca", rscript = "/group/director1234/software/zeus/apps/gcc/4.8.3/r/3.2.3/lib64/R/bin/Rscript", snowlib = "/home/nliu/R/x86_64-pc-linux-gnu-library/3.2")
    cl <- makeCluster(c( rep(list(lnxOptions), 2)), type = "SOCK")
    x<-NDVI_mon_82_13$NDVI

    nc<-length(cls)

    clusterExport(cl,c("x"))    # share default data for all threads

    system.time(
      STA<-parLapply(cl,seq(1,564400),f_dp,data=x,year_start=1982,year_end=2013) # using parLapply to simulate data in parallel way
    )
    ## combin general returned data
    STA<-do.call(rbind,STA)
    save(STA,file="STA.RData")
    stopCluster(cl)

  }else{
    print("using snowfall and ff package to simulate")
    x<-as.ffdf(NDVI_mon_82_13)
    cores<-detectCores()-1
    sfInit(parallel=TRUE, cpus=cores, type="SOCK")
    sfLibrary(ff)
    sfLibrary(bfast)
    sfLibrary(trend)
    sfExport("x")
    sfClusterSetupRNG()
    #system.time(ls<-sfLapply(1:564400, f_change,data=x,year_start=1982,year_end=2013,variable="NDVI"))
    system.time(ls<-sfLapply(1:564400, f_mk,data=x,year_start=1982,year_end=2013,variable="NDVI"))
    la<-do.call("rbind",ls)
    save(la,file="mk.RData")
    sfStop()
  }



}

##Transfer monthly frame data to annual data by fun="sum" ot "mean"
f_m2y<-function(data, fun="mean"){

  .linshi<-melt(data,id=c(1,2,3))
  .out<-dcast(.linshi, ID+YEAR~variable, get(fun), na.rm=TRUE)
  return(.out)

}

##Transfer grid frame data to basin data by fun="mean"
f_grid2basin<-function(data,type="annual",fun="mean"){
  if(type=="HUC"){

    .linshi<-melt(data,id=c(1,2))
    .out<-dcast(.linshi, BASIN~variable, get(fun), na.rm=TRUE)
    return(.out)
  }else if(type=="annual"){

    .linshi<-melt(data,id=c(1,2,3))
    .out<-dcast(.linshi, BASIN+YEAR~variable, get(fun), na.rm=TRUE)
    return(.out)

  }else if(type=="month"){

    .linshi<-melt(data,id=c(1,2,3,4))
    .out<-dcast(.linshi, BASIN+YEAR+MONTH~variable, get(fun), na.rm=TRUE)
    return(.out)
  }
}
## summary funtion which can output summary information for all data frame objects in memory
f_summary<-function(){
  print("print info for all dataframe objects")
  a<-ls(envir=.GlobalEnv)
  print(a)
  for (i in c(1:length(a))){
    if ( is.data.frame(get(a[i]))){
      print(a[i])
      str(get(a[i]))
      print(summary.data.frame(get(a[i])))
    }
  }

}

## summary funtion which can output summary information for all list objects in memory
f_list_summary<-function(){
  print("print info for all list objects")
  a<-ls(envir=.GlobalEnv)
  #print(a)
  for (i in c(1:length(a))){
    if (is.list(get(a[i]))){
      print(a[i])
      str(get(a[i]))
      len<-length(get(a[i]))
      for (j in 1:len){
        if (is.data.frame(get(a[i])[[j]])){
          print(names(get(a[i])[j]))
          print(summary.data.frame(get(a[i])[[j]]))
        }
      }
    }
  }
}

####################################################################
## changepoint detection using "bfast" package and MK test using "trend" package
## in seasonal (default) and annual scale
## changepoint detection using "bfast" package
## http://www.sciencedirect.com/science/article/pii/S003442570900265X
######################################################################

f_dp<-function(data,seasonal=TRUE,year_start,year_end){
  require(bfast)
  require(trend)
  require(pryr)
  #.start<-(n-1)*((year_end-year_start+1)*12)+1
  #.end<-n*((year_end-year_start+1)*12)
  #print(n)
  print(mem_used())

  if(seasonal){
    # seasonal changepoint and trend detection
    print("seasonal scale process")
    # IF all data are NA or equal, there would be no trend
    if (!(any(is.na(data)) | isTRUE(all.equal(data, rep(data[1], length(data)))))){
      .linshi<-ts(data,frequency = 12,start = c(year_start,1))

      rdist<-12/length(.linshi)
      # ratio of distance between breaks (time steps) and length of the time series

      fit <- bfast(.linshi,h=rdist, season="harmonic", max.iter=1,breaks=2)
      .out<-fit$output[[1]]
      if ( is.list(.out$bp.Vt)){.trend_change<-.out$bp.Vt$breakpoints}else{.trend_change<-NA}
      if ( is.list(.out$ci.Wt)){.season_change<-.out$ci.Wt[[1]][2]}else{.season_change<-NA}

      ## MK trend detection using "trend" package

      .outmk<-smk.test(.linshi)
      .outslope<-sea.sens.slope(.linshi)

    }else{

      ##---change point detect result

      .trend_change<-NA
      .season_change<-NA

      ## MK test result

      .outmk<-data.frame(tautot=NA,pvalue=NA)
      .outslope<-data.frame(b.sen=NA)

    }
    return(list(CP_trend=.trend_change,CP_season=.season_change,TAU=.outmk$tautot,PMK=.outmk$pvalue,SLOPE=.outslope$b.sen))
  }else{
    # annual changepoint and trend detection
    print("annual scale process")
    # IF all data are NA or equal, there would be no trend
    if (!(any(is.na(data)) | isTRUE(all.equal(data, rep(data[1], length(data)))))){
      .linshi<-ts(data,frequency = 1, start = year_start)

      rdist<-3/length(.linshi)
      #print(.linshi)
      fit <- bfast(.linshi,h=rdist,season = "none", max.iter=1,breaks=2)

      .out<-fit$output[[1]]

      if ( is.list(.out$bp.Vt)){.trend_change<-.out$bp.Vt$breakpoints}else{.trend_change<-NA}

      ## MK trend detection using "trend" package

      .outmk<-mk.test(.linshi)
      .outslope<-sens.slope(.linshi)

    }else{

      ##---change point detect result

      .trend_change<-NA

      ## MK test result

      .outmk<-data.frame(tautot=NA,pvalue=NA)
      .outslope<-data.frame(b.sen=NA)

    }
    return(list(CP_trend=.trend_change,TAU=.outmk$tautot,PMK=.outmk$pvalue,SLOPE=.outslope$b.sen))
  }
}

####################################################################
## this function is used for plot Brick result
## 1, the dataframe data will be transfer to brick
## 2, using ggplot to plot brick
#####################################################################
f_grid_plot<-function(data,info,annual=FALSE,monthly=FALSE,plot=FALSE){
  #info<-c("latmin"=1, "latmax"=1, "longmin"=1, "longmax"=1, "ncols"=1, "nrows"=1, "nbands"=1,"year_start"=1, "year_end"=1,"annual"=0,"monthly"=0)
  ## data is the original data frame and info consists the required infor mation for transfering data from frame to brick
  require(raster)
  require(plyr)
  require(rasterVis)
  require(ggplot2)

  .brick.plot<-function(bricks,name){

    .gplot<-gplot(bricks) + geom_tile(aes(fill = value)) +
      facet_wrap(~ variable) + #,ncol=3
      #scale_fill_gradient(low = 'white', high = 'blue') +
      scale_fill_gradientn(colours=c("blue","green","yellow","red"),name=name,na.value = "grey70")+ #limits=c(500,1000),
      coord_equal()+
      theme_set(theme_bw())

    #ggsave(filename, plot = last_plot(), device = NULL, path = NULL, scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"), dpi = 300, limitsize = TRUE, ...)
    #ggsave(gplot,file =paste("trend/","P_Q.pdf",sep=""),width = 10,  units = c("cm"),dpi = 300)
    filename<-paste("images/",name,".pdf",sep="")
    if (plot){print(.gplot)}
    ggsave(.gplot,file=filename,width = 20,  units = c("cm"),dpi=300)

    #print(class(a))
  }

  if(annual){

    if(monthly){
      ## this is for monthly result
      data<-arrange(data,YEAR,MONTH,ID)
      for (var in 4:length(data)){
        nbands=(info["year_end"]-info["year_start"]+1)*12
        .array<-array(data[[var]],c(info["ncols"],info["nrows"],nbands))
        .brick<-brick(.array,xmn=info["longmin"],xmx=info["longmax"],ymn=info["latmin"],ymx=info["latmax"],crs="+proj=longlat+datum=WGS84")
        names(.brick)<-paste("Monthly",names(data)[var],paste(rep(c(info["year_start"]:info["year_end"]),each=12),c(1:12)))
        print(paste("Monthly",names(data)[var]))
        .name<-paste("Monthly",names(data)[var],sep="_")
        .brick.plot(.brick,.name)

        #print(.brick)
      }
    }else{
      ## this is for annual result
      data<-arrange(data,YEAR,ID)
      for (var in 3:length(data)){
        nbands=info["year_end"]-info["year_start"]+1
        .array<-array(data[[var]],c(info["ncols"],info["nrows"],nbands))
        .brick<-brick(.array,xmn=info["longmin"],xmx=info["longmax"],ymn=info["latmin"],ymx=info["latmax"],crs="+proj=longlat+datum=WGS84")
        names(.brick)<-paste("Annual",names(data)[var],c(info["year_start"]:info["year_end"]))
        #print(.brick)
        print(paste("Annual",names(data)[var]))
        .name<-paste("Annual",names(data)[var],sep="_")
        .brick.plot(.brick,.name)
      }
    }
  }else{
    ## this is for HUC result
    for (var in 2:length(data)){
      nbands=1
      .array<-array(data[[var]],c(info["ncols"],info["nrows"],nbands))
      .brick<-brick(.array,xmn=info["longmin"],xmx=info["longmax"],ymn=info["latmin"],ymx=info["latmax"],crs="+proj=longlat+datum=WGS84")
      names(.brick)<-paste("HUC",names(data)[var])
      #print(.brick)
      .name<-paste("HUC",names(data)[var],sep="_")
      print(paste("HUC",names(data)[var]))
      .brick.plot(.brick,.name)

    }
  }

}

####################################################################
## this function is used for plot scatter valid result
## 1, the dataframe data will be transfer to brick
## 2, using ggplot to plot brick
#####################################################################
f_scatter_plot<-function(data,info,annual=FALSE,monthly=FALSE){
  cof<-coef(lm(Q[3:9] ~ Observed[3:9], data = ann_mean_MJ))

  ggplot(ann_mean_MJ, aes(x=Observed[3:9], y=Q[3:9])) +
    geom_point(size=4) +    # Use hollow circles
    geom_abline(intercept = cof[1], slope = cof[2]) +   # Don't add shaded confidence region
    #geom_abline(intercept = 0, slope = 1,linetype="dashed") +
    scale_x_continuous(name="Observed annual runoff (mm)") +
    scale_y_continuous(name="Simulated annual runoff (mm)")+#limits=c(300, 700)
    theme(axis.title.x = element_text(family="Times",face="bold", colour="black", size=12),
          axis.title.y  = element_text(family="Times",face="bold", colour="black", size=12),
          axis.text.x  = element_text(family="Times",face="bold",size=10),
          axis.text.y  = element_text(family="Times",face="bold",size=10))+
    annotate("text",family="Times", x = 500, y = 525, label = "Y = 0.94 * X - 93.8", fontface="italic",size=8)+
    annotate("text",family="Times", x = 500, y = 500, label="R^2 = 0.75\n RMSE = 135 mm", size=6)

  #ylab(expression("today's temperature is "*-5~degree*C))
  #qplot(1,1) + ylab(expression(Temp^2))

  ####-------plot veg water balance box
}

f_box_plot<-function(name1){
  g_plot<-ggplot(data = ann_mean_main_veg, aes(x = VEG, y = ann_mean_main_veg[[a]])) +
    stat_boxplot(geom = "errorbar", stat_params = list(width = 0.5), geom_params = list()) +
    geom_boxplot() + xlab("Vegetation types") +
    ylab(name1) + theme_bw(base_size = 16, base_family = "Times")
  ggsave(g_plot,file =paste("box/",names(ann_mean_MJ),".pdf",sep="")[a],dpi = 300)
  #print(g_plot)
}

###   plot annual mean line
f_line_plot<-function(name1){
  r<-coef(lm(mean_ann_MJ_Y[[a]] ~ YEAR, data = mean_ann_MJ_Y))
  print(r[2])
  l_plot<- ggplot(data = mean_ann_MJ_Y, aes(x = YEAR, y = mean_ann_MJ_Y[[a]])) + geom_point(size=4,shape=21, fill="white") +
    geom_line(size = 1) + scale_x_continuous(breaks=2002:2014)+
    xlab("YEAR") + ylab(name1) + theme_bw(base_size = 14, base_family = "Times") +
    geom_abline(intercept = r[1], slope = r[2])
  ggsave(l_plot,file =paste("line/",names(mean_ann_MJ_Y),".pdf",sep="")[a],dpi = 300)
  print(names(mean_ann_MJ_Y)[a])
}


#################################################
#multiplot for ploting multi row and cols ggplot function

################################################
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
# convert daily data to annual (zoo)
daily2annual<-function (x, FUN, na.rm = TRUE, out.fmt = "%Y-%m-%d", ...)
{
  if (missing(FUN))
    stop("Missing argument value: 'FUN' must contain a valid function for aggregating the values")
  if (sfreq(x) %in% c("annual"))
    stop("Invalid argument: 'x' is already an annual ts !!")
  if (is.na(match(out.fmt, c("%Y", "%Y-%m-%d"))))
    stop("Invalid argument: 'out.fmt' must be in c('%Y', '%Y-%m-%d')")
  dates <- time(x)
  y <- as.numeric(format(dates, "%Y"))
  years <- factor(y, levels = unique(y))
  tmp <- aggregate(x, by = years, FUN, na.rm = na.rm)
  nan.index <- which(is.nan(tmp))
  if (length(nan.index) > 0)
    tmp[nan.index] <- NA
  inf.index <- which(is.infinite(tmp))
  if (length(inf.index) > 0)
    tmp[inf.index] <- NA
  if (out.fmt == "%Y") {
    time(tmp) <- format(time(tmp), "%Y")
  }
  else time(tmp) <- as.Date(paste(time(tmp), "-01-01", sep = ""))
  if (NCOL(tmp) == 1)
    tmp <- zoo(as.numeric(tmp), time(tmp))
  return(tmp)
}

# convert daily data to monthly (zoo)
daily2monthly<-function (x, FUN, na.rm = TRUE, ...)
{
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

}

#' Transfering matrix/array to raster/brick
#' @param data A matrix or array object.
#' @param infonc A filename of a raster object for getting the raster extent info.
#' @keywords cats
#' @export
#' @examples
#' rc<-f_2raster(darray,infonc="/Dataset/backup/CABLE/ET_ann_82_14.nc")
#'

# This is for transfering array to raster
f_2raster<-function(data,infonc=NA){
  #infonc is a target raster
  require(raster)
  if(is.na(infonc)){
    info<-raster("/Dataset/backup/CABLE/ET_ann_82_14.nc")
  }else{
    info<-raster(infonc)
  }

  if(is.matrix(data)){
    .grid<-raster(data,xmn=info@extent@xmin,xmx=info@extent@xmax,ymn=info@extent@ymin,ymx=info@extent@ymax,crs=crs(info))

  }else{
    .grid<-brick(data,xmn=info@extent@xmin,xmx=info@extent@xmax,ymn=info@extent@ymin,ymx=info@extent@ymax,crs=crs(info))
  }
  return(.grid)
}

#' Zonal raster/brick based on a shapefile
#' @param ncfilename A filename for a "raster*" type object.
#' @param basin A ploygon object.
#' @param fun function for doing zonal (fun="mean"/"sum")
#' @param varname define the variable name
#' @param zonal_field select the field from shapefile file for naming the result
#' @param start the start year for the time series of the input raster
#' @param scale the time step the the input raster
#' @param weight Whether weight polygon for mean
#' @keywords zonal
#' @export
#' @examples
#' sta_shp<-f_sta_shp_nc(ncfilename="/Dataset/backup/CABLE/ET_ann_82_14.nc",
#' basin,fun="mean",varname="ET",zonal_field="Station",start=1982,scale="annual")
#'

# zonal brick

f_sta_shp_nc<-function(ncfilename,basin,fun="mean",varname,zonal_field,start,scale="month",df=T,weight=T){
  require(reshape2)
  require(raster)
  da<-brick(ncfilename)
  NAvalue(da)<- 0
  if(fun=="mean" | fun=="Mean" | fun=="MEAN"){
    ex <- raster::extract(da, basin, fun=mean, na.rm=TRUE, df=df,weights=weight,normalizWeights=TRUE, small=TRUE)
  }else{
    ex <- raster::extract(da, basin, fun=sum, na.rm=TRUE, df=df,normalizWeights=TRUE, small=TRUE)
  }

  if (df){
    sta_catchment<-as.data.frame(t(ex[-1]))

    if(scale=="month" | scale=="Month" | scale=="MONTH"){
      dates<-seq(as.Date(paste(start,"01-01",sep="-")),by="1 month",length.out = dim(da)[3])
      sta_catchment$Year<-as.integer(format(dates,"%Y"))
      sta_catchment$Month<-as.integer(format(dates,"%m"))
      names(sta_catchment)<-c(as.character(basin[[zonal_field]]),"Year","Month")
      sta_catchment<-melt(sta_catchment,id=c("Year","Month"))
      names(sta_catchment)<-c("Year","Month",zonal_field,varname)

    }else if(scale=="annual" | scale=="Annual" | scale=="ANNUAL"){
      dates<-seq(as.Date(paste(start,"01-01",sep="-")),by="1 year",length.out = dim(da)[3])
      sta_catchment$Year<-as.integer(format(dates,"%Y"))
      names(sta_catchment)<-c(as.character(basin[[zonal_field]]),"Year")
      sta_catchment<-melt(sta_catchment,id=c("Year"))
      names(sta_catchment)<-c("Year",zonal_field,varname)

    }else{
      dates<-seq(as.Date(paste(start,"01-01",sep="-")),by="1 day",length.out = dim(da)[3])
      sta_catchment$Year<-as.integer(format(dates,"%Y"))
      sta_catchment$Month<-as.integer(format(dates,"%m"))
      sta_catchment$Day<-as.integer(format(dates,"%d"))
      names(sta_catchment)<-c(as.character(basin[[zonal_field]]),"Year","Month","Day")
      sta_catchment<-melt(sta_catchment,id=c("Year","Month","Day"))
      names(sta_catchment)<-c("Year","Month","Day",zonal_field,varname)
    }
  }else{
    sta_catchment<-ex

    if(scale=="month" | scale=="Month" | scale=="MONTH"){
      # names(sta_catchment)<-as.character(seq(as.Date(paste(start,"01-01",sep="-")),by="1 month",length.out = dim(da)[3]))
      rownames(sta_catchment)<-as.character(basin[[zonal_field]])

    }else if(scale=="annual" | scale=="Annual" | scale=="ANNUAL"){
      # names(sta_catchment)<-as.character(seq(as.Date(paste(start,"01-01",sep="-")),by="1 year",length.out = dim(da)[3]))
      rownames(sta_catchment)<-as.character(basin[[zonal_field]])

    }else{
      #  names(sta_catchment)<-as.character(seq(as.Date(paste(start,"01-01",sep="-")),by="1 day",length.out = dim(da)[3]))
      rownames(sta_catchment)<-as.character(basin[[zonal_field]])
    }

  }
  sta_catchment
}
