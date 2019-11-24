# This is the global Rscript for this app

# function for checking libs
f_lib_check<-function(libs){
  for (lib in libs ){
    if(lib %in% rownames(installed.packages())){

    }else{
      install.packages(lib,repos='http://cran.us.r-project.org')
    }
  }

  a<-lapply(libs, require, character.only = TRUE)
}

theme_ning<-function(size.axis=5,size.title=6){
  theme_bw(base_family = "serif") %+replace%
    theme(axis.title = element_text(face="bold", colour="black", size=size.title),
          axis.text= element_text(angle=0, vjust=0.3, size=size.axis),
          legend.title = element_text(colour="black", size=size.axis, face="bold"),
          legend.text = element_text(colour="black", size = size.axis),
          strip.text.x = element_text(size = size.axis,margin=margin(4, 2, 6, 2), face="bold"),
          strip.text.y = element_text(size = size.axis,margin=margin(4, 2, 4, 6), face="bold",angle=-90),
          legend.key.size=unit(1.2, "lines"),
          legend.box.spacing=unit(1, "mm"),
          strip.background = element_blank(),
          plot.title = element_text(vjust = 2.5,hjust = 0.5,face="bold")
    )
}

## Paste one to one for two vectors, matrixes or arrays----
#' Paste by value one to one for two vectors, matrixes or arrays
#' @param x The first object, which can be vector, matrix or array.
#' @param y The second object, which can be vector, matrix or array.
#' @param sep The separate letter
#' @keywords paste
#' @export
#' @examples
#' x<-c(1,2,3)
#' y<-c("A","B","C")
#' f_paste(x,y,sep="-")
f_paste<-function(x,y,sep=""){
  dimx<-dim(x)
  if(is.null(dimx)){
    sapply(c(1:length(x)),function(a) paste(x[a],y[a],sep=sep))
  }else{
    pas<-sapply(c(1:length(x)),function(a) paste(as.vector(x)[a],as.vector(y)[a],sep=sep))
    if(length(dimx)==2){
      matrix(pas,dimx)
    }else{
      array(pas,dimx)
    }
  }
}

f_sta_shp_nc_new<-function(ncfilename,basin,fun="mean",varname,zonal_field,yr.start,scale="1 month",weight=T,plot=T){
  require(dplyr)
  require(raster)
  require(tidyr)
  require(sp)
  da<-brick(ncfilename)
  if (!compareCRS(basin,da))  basin<-spTransform(basin,crs(da))
  da<-crop(da,basin)

  if(plot) {
    plot(da[[1]],basin)
    plot(basin,add=T)
  }

  beginCluster()
  basin_r<-rasterize(basin,da[[1]],field=zonal_field)
  endCluster()

  dates<-seq(as.Date(paste0(yr.start,"-01-01")),by=scale,length.out = dim(da)[3])

  da_matrix<-cbind(values(basin_r),values(da))
  colnames(da_matrix)<-c("BasinID",as.character(dates))

  da_sta<-da_matrix%>%
      as.data.frame()%>%
      filter(!is.na(BasinID))%>%
      group_by(BasinID)%>%
      summarise_all(.funs = fun,na.rm=T)%>%
      melt(id="BasinID")%>%
      rename(Date=variable)%>%
      mutate(Date=as.Date(Date))
  names(da_sta)[3]<-varname
  names(da_sta)[1]<-zonal_field
  return(da_sta)
}

f_sta_shp_nc<-function(ncfilename,basin,fun="mean",varname,zonal_field,yr.start,scale="month",weight=T,plot=T){
  require(dplyr)
  require(raster)
  require(tidyr)
  require(sp)
  da<-brick(ncfilename)
  if (!compareCRS(basin,da))  basin<-spTransform(basin,crs(da))
  da<-crop(da,basin)
  #NAvalue(da)<- 0
  if(plot) {
    plot(da[[1]],basin)
    plot(basin,add=T)
  }

  beginCluster()
  if(fun=="mean" | fun=="Mean" | fun=="MEAN"){
    ex <- raster::extract(da, basin, fun=mean, na.rm=TRUE, weights=weight)
  }else{
    ex <- raster::extract(da, basin, fun=sum, na.rm=TRUE)
  }
  endCluster()

  if(scale=="month" | scale=="Month" | scale=="MONTH"){
    dates<-seq(as.Date(paste0(yr.start,"-01-01")),by="1 month",length.out = dim(da)[3])
    sta_catchment<-t(ex)%>%
      round(digits = 3)%>%
      as.data.frame()%>%
      mutate(Year=as.integer(format(dates,"%Y")),
             Month=as.integer(format(dates,"%m")))%>%
      gather(BasinID,values,1:length(basin))%>%
      mutate(BasinID=rep(basin[[zonal_field]],each=length(dates)))%>%
      dplyr::select(BasinID,Year,Month,values)
    names(sta_catchment)<-c(zonal_field,"Year","Month",varname)

  }else if(scale=="annual" | scale=="Annual" | scale=="ANNUAL"){
    dates<-seq(as.Date(paste0(yr.start,"-01-01")),by="1 year",length.out = dim(da)[3])
    sta_catchment<-t(ex)%>%
      round(digits = 3)%>%
      as.data.frame()%>%
      mutate(Year=as.integer(format(dates,"%Y")))%>%
      gather(BasinID,values,1:length(basin))%>%
      mutate(BasinID=rep(basin[[zonal_field]],each=length(dates)))%>%
      dplyr::select(BasinID,Year,values)

    names(sta_catchment)<-c(zonal_field,"Year",varname)

  }else{
    dates<-seq(as.Date(paste0(yr.start,"-01-01")),by="1 day",length.out = dim(da)[3])

    sta_catchment<-t(ex)%>%
      round(digits = 3)%>%
      as.data.frame()%>%
      mutate(Year=as.integer(format(dates,"%Y")),
             Month=as.integer(format(dates,"%m")),
             Day=as.integer(format(dates,"%d")))%>%
      gather(BasinID,values,1:length(basin))%>%
      mutate(BasinID=rep(basin[[zonal_field]],each=length(dates)))%>%
      dplyr::select(BasinID,Year,Month,Day,values)
    names(sta_catchment)<-c(zonal_field,"Year","Month","Day",varname)

  }

  sta_catchment
}
# this function is used for zonal LAI of each HRU in by a shp file
hru_lc_zonal<-function(classname,daname,shp,fun='mean',field=NULL,plot=T){
  require(raster)
  require(sp)
  # read the class and data by their names
  class<-raster(classname)
  da<-brick(daname)
  if(sum(res(da)==res(class))==2) crs(da)<-crs(class)
  # crop data based on the input shp
  if (!compareCRS(shp,da)) shp<-spTransform(shp,crs(da))
  da<-crop(da,shp)

  if (!compareCRS(shp,class)) shp<-spTransform(shp,crs(class))
  class<-crop(class,shp)

  # resample class data based on the input data (Their geo reference could be different)
  #class<-projectRaster(class,da[[1]],method='ngb')
  #da<-projectRaster(da,class,method='ngb')
  # class<-raster::mask(class,shp)
  # da<-raster::mask(da,shp)

  shp@data[,field]<-as.character(shp@data[,field])

  # get the number class and their percentage and plot some base map
  print(raster::unique(class))
  if(plot){
    nclass<-raster::unique(class)
    print("percentage for each lc")
    print(round(table(matrix(class))/sum(table(matrix(class)))*100,2))
    plot(class)
    plot(da[[1]],add=T,alpha=0.5)
    plot(shp,add=T)
  }

  # funtion for zonal each polygon
  f_zonal<-function(i){
    #print(i)
    polygon1<-shp[i,]
    class1<-crop(class,polygon1)
    if (!compareCRS(polygon1,da)) polygon1<-spTransform(polygon1,crs(da))
    da1<-crop(da,polygon1)
    if (!compareCRS(polygon1,class1)) polygon1<-spTransform(polygon1,crs(class1))
    polygon1<-spTransform(polygon1,crs(class1))
    if(!sum(res(da1)==res(class1))==2){
      beginCluster()
      da1<-projectRaster(da1,class1,method='ngb')
      endCluster()
    }
    if(!extent(class1)==extent(da1)) extent(class1)=extent(da1)
    class1<-raster::mask(class1,polygon1)
    da1<-raster::mask(da1,polygon1)
    if(sum(unique(class1))<1){
      da_zonal1<-NA
      return(da_zonal1)
    }else if(ncell(class1)==1){
      da_zonal1<-data.frame("lai"=values(da1)[1,])
      colnames(da_zonal1)<-paste0("Lc_",unique(class1))
      return(da_zonal1)
    }
    da_zonal1<-zonal(da1, class1, fun)

    if(nrow(da_zonal1)<2){
      namesls<-paste0("Lc_",da_zonal1[,1])
      da_zonal1<-data.frame(namesls=da_zonal1[1,-1])
      colnames(da_zonal1)<-namesls
    }else{
      namesls<-paste0("Lc_",da_zonal1[,1])
      da_zonal1<-t(da_zonal1[,-1])
      colnames(da_zonal1)<-namesls
    }

    return(da_zonal1)
  }

  # Run sta
  if(length(shp)>1){
    da_zonal<- lapply(c(1:length(shp)), f_zonal)
    names(da_zonal)<-shp@data[,field]
  }else{
    da_zonal<-zonal(da, class, fun)
    namesls<-paste0("Lc_",da_zonal)
    da_zonal<-t(da_zonal[,-1])
    colnames(da_zonal)<-namesls
  }

  return(da_zonal)
}

hru_lc_ratio_new<-function(classname,shp,field=NULL){
  library(raster)
  require(sp)

  class<-raster(classname)
  if (!compareCRS(shp,class)) shp<-spTransform(shp,crs(class))
  beginCluster()
  shp_r<-rasterize(shp,class,field="BasinID")
  endCluster()

  class_ratio<-data.frame("BasinID"=values(shp_r),"Class"=values(class))%>%
    filter(!is.na(BasinID))%>%
    group_by(BasinID,Class)%>%
    summarise(n=n())%>%
    dcast(BasinID~Class)
  names(class_ratio)<-c("BasinID",paste0("Lc_",names(class_ratio)[-1]))
  class_ratio[is.na(class_ratio)]<-0
  Sum_n<-apply(class_ratio[,-1], 1,sum)
  for (i in c(2:ncol(class_ratio))) class_ratio[,i]<-round(class_ratio[,i]/Sum_n,3)

  return(class_ratio)
}

hru_lc_ratio<-function(classname,shp,field=NULL,mcores=1){
  library(raster)
  require(sp)

  class<-raster(classname)
  if (!compareCRS(shp,class)) shp<-spTransform(shp,crs(class))
  class<-crop(class,shp)
  class<-mask(class,shp)
  print(table(matrix(class)))

  # zonal_for each polygon
  f_zonal<-function(i,shp=shp,class=class){
    polygon1<-shp[i,]
    class1<-crop(class,polygon1)
    class_ratio<-as.data.frame(table(matrix(class1))/sum(table(matrix(class1))))
    names(class_ratio)<-c("Class","Ratio")
    class_ratio$Ratio<-round(class_ratio$Ratio,3)
    class_ratio$Count<-table(matrix(class1))
    class_ratio[field]<-polygon1@data[field]
    return(class_ratio)
  }

  # Run sta
  if(length(shp)>1){
    if(.Platform$OS.type=="windows"){
      # cl<-makeCluster(mcores, type="SOCK")
      # clusterExport(cl, c("shp","class","f_zonal"))
      # result<-clusterApply(cl,c(1:length(shp)), f_zonal,shp=shp,class=class)
      # stopCluster(cl)
      lc_ratio<- lapply(c(1:length(shp)), f_zonal,shp=shp,class=class)
    }else{
    lc_ratio<- mclapply(c(1:length(shp)), f_zonal,shp=shp,class=class,mc.cores=mcores)
    }
    lc_ratio<-do.call(rbind,lc_ratio)
  }else{
    class_ratio<-as.data.frame(table(matrix(class))/sum(table(matrix(class))))
    names(class_ratio)<-c("Class","Ratio")
    class_ratio$Ratio<-round(class_ratio$Ratio,3)
    class_ratio$Count<-table(matrix(class))
    class_ratio[field]<-polygon1@data[field]
    lc_ratio<-class_ratio
  }

  return(lc_ratio)
}

f_landlai_new<-function(lcfname,laifname,Basins,byfield,yr.start,scale="1 month"){

  require(raster)
  require(sp)
  # read the class and data by their names
  lcclass<-raster(lcfname)
  da<-brick(laifname)
  if(sum(res(da)==res(lcclass))==2) crs(da)<-crs(lcclass)
  # crop data based on the input shp
  if (!compareCRS(Basins,da)) Basins<-spTransform(Basins,crs(da))
  da<-crop(da,Basins)

  if (!compareCRS(Basins,lcclass)) shp<-spTransform(Basins,crs(lcclass))
  lcclass<-crop(lcclass,Basins)

  beginCluster()
  Basins_r<-rasterize(Basins,lcclass,field=byfield)
  endCluster()

  da_matrix<-cbind(values(Basins_r),values(lcclass),values(da))
  dates<-seq(as.Date(paste0(yr.start,"-01-01")),by=scale,length.out = dim(da)[3])
  colnames(da_matrix)<-c("BasinID","Class",as.character(dates))

  hru_lais<-da_matrix%>%
    as.data.frame()%>%
    filter(!is.na(BasinID))%>%
    mutate(Class=paste0("Lc_",as.integer(as.character(Class))))%>%
    group_by(BasinID,Class)%>%
    summarise_all(.funs = "mean",na.rm=T)%>%
    melt(id=c("BasinID","Class"))%>%
    rename(Date=variable,LAI=value)%>%
    mutate(Date=as.Date(Date))%>%
    mutate(Year=as.integer(format(Date,"%Y")),
         Month=as.integer(format(Date,"%m")))%>%
    dplyr::select(BasinID,Year,Month,Class,LAI)%>%
    dcast(BasinID+Year+Month~Class)%>%
    arrange(BasinID,Year,Month)

  hru_lais[is.na(hru_lais)]<-0
  names(hru_lais)[1]<-byfield
  return(hru_lais)

}


# this function is used for zonal LAI of each lc in the HRU
f_landlai<-function(lcfname,laifname,Basins,byfield,yr.start){
  hru_lai<-hru_lc_zonal(classname = lcfname,
                        daname =laifname,
                        shp = Basins,
                        field = "BasinID")
  lcs<-paste0("Lc_",unique(raster(lcfname)))

  f_fillmatix<-function(a,lcs){
    a<-round(a,3)
    prel<-length(a[1,])+1
    if(prel< length(lcs)+1){
      lacks<-lcs[which(! lcs %in%  colnames(a))]
      for (i in c(1:length(lacks))) a<-cbind(a,lcadd=0)
      colnames(a)[prel:length(a[1,])]<-lacks
      a<-a[,lcs]
    }
    a
  }
  ha<-lapply(hru_lai, f_fillmatix,lcs)
  hru_lais<-do.call(rbind,ha)
  hru_lais<-cbind("BasinID"=rep(as.integer(names(hru_lai)),each=length(hru_lai[[1]][,1])),
                  "Year"=rep(c(as.integer(yr.start):(length(hru_lai[[1]][,1])/12+as.integer(yr.start)-1)),each=12),
                  "Month"=c(1:12),
                  hru_lais)
  hru_lais<-as.data.frame(hru_lais)
  hru_lais<-arrange(hru_lais,BasinID,Year,Month)
  hru_lais[is.na(hru_lais)]<-0
  return(hru_lais)
}

f_cellinfo<-function(classfname,Basins,byfield="BasinID",demfname=NULL){
  require(tidyr)
  require(dplyr)
  require(raster)
  # hru_lcs<-hru_lc_ratio(classname =classfname,
  #                       shp = Basins,
  #                       field = byfield)%>%
  #   mutate(Class=paste0("Lc_",Class))%>%
  #   dplyr::select(-Count)%>%
  #   spread(Class, Ratio,fill=0)
  hru_lcs<-hru_lc_ratio_new(classname =classfname,
                        shp = Basins,
                        field = byfield)
  print("getting elevation for each HUC")
  if(!"Elev_m" %in% names(Basins) & !is.null(demfname)){
    dem<-raster(demfname)
    # beginCluster()
    # .a<-raster::extract(dem,Basins,df=T,fun=mean,na.rm=T,weight=T)
    # endCluster()
    #
    if (!compareCRS(Basins,dem)) Basins<-spTransform(Basins,crs(dem))
    beginCluster()
    Basins_r<-rasterize(Basins,dem,field="BasinID")
    endCluster()

    .a<-data.frame("BasinID"=values(Basins_r),"dem"=values(dem))%>%
      filter(!is.na(BasinID))%>%
      group_by(BasinID)%>%
      summarise(dem=mean(dem,na.rm=T))
    Basins$Elev_m<-round(.a[,2],2)
  }

  cellinfo<-Basins@data%>%
    dplyr::select(one_of(c(byfield,"Area_m2","Latitude","Longitude","Elev_m","Flwlen_m")))%>%
    arrange(get(byfield))%>%
    mutate(Area_m2=round(Area_m2,1))%>%
 #   mutate(Elev_m=round(Elev_m,1))%>%
    mutate(Latitude=round(Latitude,3))%>%
    mutate(Longitude=round(Longitude,3))%>%
    merge(hru_lcs,by=byfield)

  if("Elev_m" %in% names(cellinfo)) cellinfo[["Elev_m"]]<-round(cellinfo[["Elev_m"]],1)

  return(cellinfo)
}
f_soilinfo<-function(soilfname,Basins){
  SOIL<-brick(soilfname)
  Basins<-spTransform(Basins,crs(SOIL))
  SOIL_catchment<-raster::extract(SOIL,Basins,fun=mean,na.rm=T,weights=T)
  # fill NA values
  SOIL_catchment[is.infinite(SOIL_catchment)]<-NA
  SOIL_catchment[is.na(SOIL_catchment)]<-0
  SOIL_catchment<-round(SOIL_catchment,4)

  colnames(SOIL_catchment)<-c("uztwm", "uzfwm" , "uzk", "zperc" , "rexp" , "lztwm" , "lzfsm",
                              "lzfpm", "lzsk" , "lzpk" , "pfree")

  SOIL_catchment<-as.data.frame(cbind(BasinID=Basins$BasinID,SOIL_catchment))
  return(SOIL_catchment)
}

f_crop_roi<-function(da,roi_shp,plot=T){
  if (!compareCRS(roi_shp,da))  roi_shp<-spTransform(roi_shp,crs(da))
  da<-crop(da,roi_shp)
  da<-mask(da,roi_shp)
  if(plot) plot(da[[1]]);plot(roi_shp,add=T)
  return(da)
}

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

f_stream_level_pete<-function(streamfile=NA,mc_cores=1){

  stream<-read.csv(streamfile)
  stream_level<-stream[c("FROM","TO")]
  stream_level$LEVEL<-NA

  lev<-1
  for (i in c(1:500)){

    if(lev==1){
      index_lev_down<-which(!stream_level$TO %in% stream_level$FROM)
      stream_level$LEVEL[index_lev_down]<-lev
      lev<-lev+1
    }

    index_lev_up<-which(stream_level$TO %in% stream_level$FROM[index_lev_down])
    stream_level$LEVEL[index_lev_up]<-lev
    index_lev_down<-index_lev_up
    lev<-lev+1
  }

  return(stream_level)
}

f_upstreamHUCs<-function(BasinID,routpar){
  level_to<-routpar$LEVEL[routpar$TO==BasinID]
  upids<-NA
  To<-BasinID
  if(length(level_to)>0){

    while (length(level_to)>0){
      FROM_HUCs<-routpar$FROM[routpar$TO %in% To]

      upids<-c(upids,FROM_HUCs)

      To<-routpar$FROM[routpar$TO %in% To]

      level_to<-routpar$LEVEL[routpar$TO %in% To]

    }

    return(upids[-1])
  }else{
    return(NULL)
  }

}

hrurouting<-function(Flwdata,routpar,mc_cores=1){
  library(parallel)
  max_level<-max(routpar$LEVEL)

  hru_accm<-function(hru,water,routpar){
    hru<-as.numeric(hru)
    water$flow[water$BasinID==hru] +sum(water$flow[water$BasinID %in% routpar$FROM[which(routpar$TO==hru)]])
  }

  for (level in c(max_level:1)){
    hrus<-unique(routpar$TO[routpar$LEVEL==level])
    #print(paste0("There are ",length(hrus)," hrus in level ",level))

    if(length(hrus)>100) {
      flowaccu<-mclapply(hrus,hru_accm,water=Flwdata,routpar=routpar,mc.cores = mc_cores)
      #print("using paralell")
    }else{
      flowaccu<-lapply(hrus,hru_accm,water=Flwdata,routpar=routpar)
    }
    for (i in c(1:length(hrus))) Flwdata$flow[Flwdata$BasinID==hrus[i]]<- flowaccu[[i]]
  }
  return(water)
}

#' This function allows you to get the number of days for a specific month.
#' @param date A date object.
#' @keywords cats
#' @export
#' @examples
#' date<-as.Date("2001-01-01")
#' numberOfDays(date)
numberOfDays <- function(date) {
  m <- format(date, format="%m")

  while (format(date, format="%m") == m) {
    date <- date + 1
  }

  return(as.integer(format(date - 1, format="%d")))
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


# Parallel wrap of WaSSI-C model----
#' @title Parallel wrap of WaSSI-C model
#' @description FUNCTION_DESCRIPTION
#' @param sim.dates  list of all dates
#' @param warmup years for warming up WaSSI-C model
#' @param mcores how many cores using for simulation
#' @param WUE.coefs the input parameters of WUE
#' @param ET.coefs the input parameter of ET
#' @return outputs potential evapotranspiration (mm day-1)
#' @details For details see Haith and Shoemaker (1987)
#' @examples
#' \dontrun{
#' dWaSSIC(sim.dates = sim_dates,
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
                   hru.lai=NULL,hru.lc.lai=NULL,huc.lc.ratio=NULL,WUE.coefs=NULL,ET.coefs=NULL)
{
  require(lubridate)
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

  if (!is.null(huc.lc.ratio)) NLCs<-length(huc.lc.ratio[1,])
  # WaSSI for each HRU and calculate adjust temp and pet values
  WaSSI<-function(h){

    # Using snowmelt function to calculate effective rainfall

    snow.result<-snow_melt(ts.prcp = hru_prcp[,h],ts.temp = hru_tavg[,h],snowrange = c(-5,1))
    hru_rain <- snow.result$prcp

    # PET using Hamon equation
    hru_pet <- hamon(par = par.petHamon[h], tavg = hru_tavg[,h], lat = hru_lat[h], jdate = jdate)
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
    if(is.null(ET.coefs)){
      hru_lc_pet<-hru_pet*hru.lc.lai[[h]]*0.0222+0.174*hru_mrain+0.502*hru_pet+5.31*hru.lc.lai[[h]]
    }else{
      # Calculate the actual ET based on SUN's ET model
      hru_lc_pet<-matrix(NA,nrow =length(hru_pet), ncol=NLCs)
      for (i in c(1:NLCs)){
        hru_lc_pet[,i]<- ET.coefs[i,"Intercept"] +
          ET.coefs[i,"P_coef"]*hru_mrain+
          ET.coefs[i,"PET_coef"]*hru_pet+
          ET.coefs[i,"LAI_coef"]*hru.lc.lai[[h]][i]+
          ET.coefs[i,"P_PET_coef"]*hru_mrain*hru_pet +
          ET.coefs[i,"P_LAI_coef"]*hru_mrain*hru.lc.lai[[h]][i] +
          ET.coefs[i,"PET_LAI_coef"] *hru_pet*hru.lc.lai[[h]][i]
      }
    }

    # combin all the input
    hru_in<-list(prcp=hru_prcp[,h],
                    temp = hru_tavg[,h],
                    rain=hru_rain,
                    snowpack=snow.result$snowpack,
                    snowmelt=snow.result$snowmelt,
                    PET_hamon=hru_pet)

    # calculate flow based on PET and SAC-SMA model
    hru_lc_out <-apply(hru_lc_pet,2,sacSma_mon,par = par.sacsma[h,], prcp = hru_mrain )

    if(!is.null(WUE.coefs)){
      # Calculate Carbon based on WUE for each vegetation type
      for (i in c(1:NLCs)){
        hru_lc_out[[i]][["GEP"]]<-hru_lc_out[[i]][["totaet"]] *WUE.coefs$WUE[i]
        hru_lc_out[[i]][["RECO"]]<-WUE.coefs$RECO_Intercept[i] + hru_lc_out[[i]][["GEP"]] *WUE.coefs$RECO_Slope[i]
        hru_lc_out[[i]][["NEE"]] <-hru_lc_out[[i]][["RECO"]]-hru_lc_out[[i]][["GEP"]]
      }
    }
    for (i in c(1:NLCs)){
      hru_lc_out[[i]][["LAI"]]<-hru.lc.lai[[h]][i]
      hru_lc_out[[i]][["PET"]]<-hru_lc_pet[,i]
    }

    # Weight each land cover type based on it's ratio
    hru_out<-lapply(c(1:length(huc.lc.ratio[h,])),function(x) lapply(hru_lc_out[[x]],function(y) huc.lc.ratio[h,x]*y))
    out<-lapply(names(hru_lc_out[[1]]), function(var) apply(sapply(hru_out, function(x) x[[var]]),1,sum))
    names(out)<-names(hru_lc_out[[1]])

    return(list(input=hru_in,out=out,lc_out=hru_lc_out))
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
  out<-list()
  out[["input"]]<- lapply(result[["input"]], function(x) as.data.frame(x)[(warmup*12+1):length(x[[1]]),])

  # routing based on catchment area
  if(is.null(par.routing) | is.null(hru_flowlen) | is.null(hru_area)){

    return(list(HUC=out))

  }else{
    # Channel flow routing from Lohmann routing model
    out2 <- mclapply(c(1:hru_num), huc_routing,mc.cores = mcores)
    out3<-lapply(names(out2[[1]]), function(var) apply(sapply(out2, function(x) x[[var]]),1,sum))

  }

  return(list(FLOW_SURF = out3[[1]], FLOW_BASE=out3[[2]],HUC=out))

}


librs<-c("dplyr","zip","lubridate","raster","ggplot2","leaflet","rgdal","rgeos","leaflet.extras","parallel","shinyFiles","tidyr","reshape2")
f_lib_check(librs)

if(!exists("data_input")) data_input<-list()
if(!exists("BasinShp")) BasinShp<-NULL
Ning<-"Ning Liu"
if(!exists("data_simulation")) data_simulation<-list()
if(!exists("Sim_dates")) Sim_dates<-list()
if(!exists("resultOutput")) resultOutput<-list()
mcores=detectCores()-1
#if(.Platform$OS.type=="windows") mcores=1
warmup<-2
