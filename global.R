# This is the global Rscript for this app

Ning<-"Ning Liu"

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

f_sta_shp_nc<-function(ncfilename,basin,fun="mean",varname,zonal_field,start,scale="month",weight=T,plot=T){
  require(dplyr)
  require(raster)
  require(tidyr)
  da<-brick(ncfilename)
  da<-crop(da,basin)
  #NAvalue(da)<- 0
  if(plot) {
    plot(da[[1]],basin)
    plot(basin,add=T)
  }
  if(fun=="mean" | fun=="Mean" | fun=="MEAN"){
    ex <- raster::extract(da, basin, fun=mean, na.rm=TRUE, weights=weight)
  }else{
    ex <- raster::extract(da, basin, fun=sum, na.rm=TRUE)
  }

  if(scale=="month" | scale=="Month" | scale=="MONTH"){
    dates<-seq(as.Date(paste0(start,"-01-01")),by="1 month",length.out = dim(da)[3])
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
    dates<-seq(as.Date(paste0(start,"-01-01")),by="1 year",length.out = dim(da)[3])
    sta_catchment<-t(ex)%>%
      round(digits = 3)%>%
      as.data.frame()%>%
      mutate(Year=as.integer(format(dates,"%Y")))%>%
      gather(BasinID,values,1:length(basin))%>%
      mutate(BasinID=rep(basin[[zonal_field]],each=length(dates)))%>%
      dplyr::select(BasinID,Year,values)

    names(sta_catchment)<-c(zonal_field,"Year",varname)

  }else{
    dates<-seq(as.Date(paste0(start,"-01-01")),by="1 day",length.out = dim(da)[3])

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
  # read the class and data by their names
  class<-raster(classname)
  da<-brick(daname)

  # crop data based on the input shp
  class<-crop(class,shp)
  da<-crop(da,shp)

  # resample class data based on the input data (Their geo reference could be different)
  class<-resample(class,da[[1]],method='ngb')
  class<-raster::mask(class,shp)
  da<-raster::mask(da,shp)

  shp@data[,field]<-as.character(shp@data[,field])

  # get the number class and their percentage and plot some base map
  print(raster::unique(class))
  if(plot){
    nclass<-raster::unique(class)
    print(table(matrix(class)))
    plot(class)
    plot(da[[1]],add=T,alpha=0.5)
    plot(shp,add=T)
  }

  # funtion for zonal each polygon
  f_zonal<-function(i){
    polygon1<-shp[i,]
    class1<-crop(class,polygon1)
    da1<-crop(da,polygon1)
    class1<-raster::mask(class1,polygon1)
    da1<-raster::mask(da1,polygon1)
    da_zonal1<-zonal(da1, class1, fun)
    namesls<-paste0("Lc_",da_zonal1[,1])
    da_zonal1<-t(da_zonal1[,-1])
    colnames(da_zonal1)<-namesls
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

hru_lc_ratio<-function(classname,shp,field=NULL,mcores=1){
  library(raster)
  class<-raster(classname)
  class<-crop(class,shp)
  class<-mask(class,shp)
  print(table(matrix(class)))

  # zonal_for each polygon
  f_zonal<-function(i){
    polygon1<-shp[i,]
    class1<-crop(class,polygon1)
    class_ratio<-as.data.frame(table(matrix(class1))/sum(table(matrix(class1))))
    names(class_ratio)<-c("Class","Ratio")
    class_ratio$Ratio<-round(class_ratio$Ratio,2)
    class_ratio$Count<-table(matrix(class1))
    class_ratio[field]<-polygon1@data[field]
    return(class_ratio)
  }

  # Run sta
  if(length(shp)>1){
    lc_ratio<- mclapply(c(1:length(shp)), f_zonal,mc.cores=mcores)
    lc_ratio<-do.call(rbind,lc_ratio)
  }else{
    class_ratio<-as.data.frame(table(matrix(class))/sum(table(matrix(class))))
    names(class_ratio)<-c("Class","Ratio")
    class_ratio$Ratio<-round(class_ratio$Ratio,2)
    class_ratio$Count<-table(matrix(class))
    class_ratio[field]<-polygon1@data[field]
    lc_ratio<-class_ratio
  }

  return(lc_ratio)
}

# this function is used for zonal LAI of each lc in the HRU
f_landlai<-function(lcfname,laifname,Basins,byfield,yr.start,yr.end){
  hru_lai<-hru_lc_zonal(classname = lcfname,
                        daname =laifname,
                        shp = Basins,
                        field = "BasinID")
  lcs<-paste0("Lc_",unique(raster(lcfname)))

  f_fillmatix<-function(a,lcs){
    a<-round(a,2)
    prel<-length(a[1,])+1
    if(prel< length(lcs)+1){
      lacks<-lcs[which(! lcs %in%  colnames(a))]
      for (i in c(1:length(lacks))) a<-cbind(a,lcadd=0)
      colnames(a)[prel:length(a[1,])]<-lacks
      a<-a[,lcs]
    }
  }

  ha<-lapply(hru_lai, f_fillmatix,lcs)
  hru_lais<-do.call(rbind,ha)
  hru_lais<-cbind("BasinID"=rep(as.integer(names(hru_lai)),each=length(hru_lai[[1]][,1])),
                  "Year"=rep(c(yr.start:yr.end),each=12),
                  "Month"=c(1:12),
                  hru_lais)
  hru_lais<-as.data.frame(hru_lais)
  hru_lais<-arrange(hru_lais,BasinID,Year,Month)
  return(hru_lais)
}

f_cellinfo<-function(classfname,Basins,byfield="BasinID",demfname=NULL){
  require(tidyr)
  hru_lcs<-hru_lc_ratio(classname =classfname,
                        shp = Basins,
                        field = byfield)%>%
    mutate(Class=paste0("Lc_",Class))%>%
    select(-Count)%>%
    spread(Class, Ratio,fill=0)

  if(!"Elev_m" %in% names(Basins) & !is.null(demfname)){
    dem<-raster(demfname)
    .a<-raster::extract(dem,Basins,df=T,fun=mean,na.rm=T,weight=T)
    Basins$Elev_m<-round(.a[,2],2)
  }

  cellinfo<-Basins@data%>%
    select(one_of(c(byfield,"Area_m2","Latitude","Longitude","Elev_m","Flwlen_m")))%>%
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

  SOIL_catchment<-raster::extract(SOIL,Basins,fun=mean,na.rm=T,weights=T)
  # fill NA values
  SOIL_catchment[is.infinite(SOIL_catchment)]<-NA
  SOIL_catchment[is.na(SOIL_catchment)]<-0
  SOIL_catchment<-round(SOIL_catchment,2)

  colnames(SOIL_catchment)<-c("uztwm", "uzfwm" , "uzk", "zperc" , "rexp" , "lztwm" , "lzfsm",
                              "lzfpm", "lzsk" , "lzpk" , "pfree")

  SOIL_catchment<-cbind(BasinID=Basins$BasinID,SOIL_catchment)
  return(SOIL_catchment)
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

librs<-c("dplyr","raster","ggplot2","leaflet","rgdal","rgeos","leaflet.extras","parallel")
f_lib_check(librs)

