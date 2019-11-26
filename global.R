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

f_read_basin<-function(fname){
  require(rgdal)
  Basins<- readOGR(fname)
  if (!"BasinID" %in% names(Basins)) Basins$BasinID<-c(1:length(Basins[,1]))
  if (is.factor(Basins$BasinID)) Basins$BasinID<-as.integer(as.character(Basins$BasinID))
  # Add latitude and longitude infor to the Basin
  if(!"Latitude"  %in% names(Basins)) {

    # get the contral coordinates of each polygon
    if(is.na(crs(Basins))) proj4string(Basins)<-CRS("+init=epsg:4326")
    Basins_wgs<-spTransform(Basins,CRS("+init=epsg:4326"))
    Basin_coords<-gCentroid(Basins_wgs, byid=TRUE)
    rownames(Basin_coords@coords)<-Basins$BasinID
    Basin_coords<-as.data.frame(Basin_coords)

    Basins[["Latitude"]]=Basin_coords$y
    Basins[["Longitude"]]=Basin_coords$x

    # Get the area of each polygon
    if (!is.projected(Basins)){
      Basins_pro<-spTransform(Basins,CRS("+init=epsg:32648"))
      Basins[["Area_m2"]]=round(gArea(Basins_pro,byid = T),2)
    }else{
      Basins[["Area_m2"]]=round(gArea(Basins,byid = T),2)
    }

  }
  return(Basins)
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
  require(reshape2)
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
      summarise_all(.funs = fun,na.rm=T)

  # Check whether it has missing HUCs
  missing_indx<-which(!basin[[zonal_field]] %in% unique(da_sta$BasinID) )
  if(length(missing_indx)>0){
    shp<-basin[missing_indx,]
    beginCluster()
    b<-raster::extract(da,shp,df=T,fun=mean,na.rm=T,weight=T)
    endCluster()
    b[,1]<-shp[[zonal_field]]
    colnames(b)<-c("BasinID",as.character(dates))
    da_sta<-rbind(da_sta,b)
  }
  da_sta<-da_sta%>%
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
    Basins_r<-rasterize(Basins,dem,field=byfield)
    endCluster()

    .a<-data.frame("BasinID"=values(Basins_r),"Elev_m"=values(dem))%>%
      filter(!is.na(BasinID))%>%
      group_by(BasinID)%>%
      summarise(Elev_m=mean(Elev_m,na.rm=T))

    # Check whether it has missing HUCs
    missing_indx<-which(!Basins$BasinID %in% .a$BasinID)
    if(length(missing_indx)>0){
      shp<-Basins[missing_indx,]
      beginCluster()
      .b<-raster::extract(dem,shp,df=T,fun=mean,na.rm=T,weight=T)
      endCluster()
      .b[,1]<-shp[[byfield]]
      names(.b)<-c("BasinID","Elev_m")
      .a<-rbind(.a,.b)
    }
    names(.a)[1]<-byfield
    Basins<-merge(Basins,.a,by=byfield)
    Basins$Elev_m<-round(Basins$Elev_m,2)
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

hru_lc_imp<-function(impname,classname,shp,byfield=NULL){
  library(raster)
  require(sp)

  class<-raster(classname)
  imp<-raster(impname)
  if (!compareCRS(shp,class)) shp<-spTransform(shp,crs(class))
  class<-crop(class,shp)

  beginCluster()
  imp<- projectRaster(imp,class)
  shp_r<-rasterize(shp,class,field=byfield)
  endCluster()

  lc_imp<-data.frame("BasinID"=values(shp_r),"Class"=values(class),imp=values(imp))%>%
    filter(!is.na(BasinID))%>%
    group_by(BasinID,Class)%>%
    summarise(imp=mean(imp,na.rm=T))%>%
    dcast(BasinID~Class)


  names(lc_imp)<-c(byfield,paste0("Lc_",names(lc_imp)[-1]))
  lc_imp[is.na(lc_imp)]<-0

  return(lc_imp)
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


#' @title Sacremento Soil Moisture Accounting Model SAC-SMA
#' @description revised based on sacsmaR package
#' @param par model parameters (11 soil parameters)
#' @param ini.states initial parameters
#' @param prcp daily precipitation data
#' @param pet potential evapotranspiration, in mm
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' sacSma_mon(pet, prcp,par)
#' }
#' @rdname sacSim_mon
#' @export
sacSma_mon <- function(pet, prcp,par,ini.states = c(0,0,500,500,500,0)) {
  if(sum(names(par) %in% c("UZTWM","UZFWM","UZK", "ZPERC",  "REXP", "LZTWM", "LZFSM", "LZFPM",  "LZSK",  "LZPK", "PFREE"))==11){
    uztwm  <-  par["UZTWM"]    # Upper zone tension water capacity [mm]
    uzfwm  <-  par["UZFWM"]    # Upper zone free water capacity [mm]
    lztwm  <-  par["LZTWM"]    # Lower zone tension water capacity [mm]
    lzfpm  <-  par["LZFPM"]    # Lower zone primary free water capacity [mm]
    lzfsm  <-  par["LZFSM"]    # Lower zone supplementary free water capacity [mm]
    uzk    <-  par["UZK"]    # Upper zone free water lateral depletion rate [1/day]
    lzpk   <-  par["LZPK"]    # Lower zone primary free water depletion rate [1/day]
    lzsk   <-  par["LZSK"]    # Lower zone supplementary free water depletion rate [1/day]
    zperc  <-  par["ZPERC"]    # Percolation demand scale parameter [-]
    rexp   <-  par["REXP"]   # Percolation demand shape parameter [-]
    pfree  <-  par["PFREE"]   # Percolating water split parameter (decimal fraction)
    pctim  <- 0 #   par[12]   # Impervious fraction of the watershed area (decimal fraction)
    adimp  <- 0 #  par[13]   # Additional impervious areas (decimal fraction)
    riva   <- 0 #  par[14]   # Riparian vegetation area (decimal fraction)
    side   <- 0 # par[15]   # The ratio of deep recharge to channel base flow [-]
    rserv  <- 0 #par[16]   # Fraction of lower zone free water not transferrable (decimal fraction)
  }else{
    print("Input soil parameter is missing")
  }

  # Initial Storage States (SAC-SMA)
  uztwc <- uztwm # Upper zone tension water storage
  uzfwc <- 1.0 # Upper zone free water storage
  lztwc <- lztwm # Lower zone tension water storage
  lzfsc <- lzfsm*0.75 # Lower zone supplementary free water storage
  lzfpc <- lzfpm*0.75 # Upper zone primary free water storage
  adimc <- 0 # Additional impervious area storage

  # CONVERT FLOW RATES IN UNITS OF 1/D TO MONTHLY RATES

  # uzk=uzk*30
  # lzpk=lzpk*30
  # lzsk=lzsk*30

  # RESERVOIR STATE ARRAY INITIALIZATION
  simaet  <- vector(mode = "numeric", length = length(prcp))
  simflow   <- vector(mode = "numeric", length = length(prcp))
  base_tot  <- vector(mode = "numeric", length = length(prcp))
  surf_tot  <- vector(mode = "numeric", length = length(prcp))
  uztwc_ts   <- vector(mode = "numeric", length = length(prcp))
  uzfwc_ts  <- vector(mode = "numeric", length = length(prcp))
  lztwc_ts  <- vector(mode = "numeric", length = length(prcp))
  lzfpc_ts  <- vector(mode = "numeric", length = length(prcp))
  lzfsc_ts  <- vector(mode = "numeric", length = length(prcp))

  thres_zero  <- 0.00001 # Threshold to be considered as zero
  parea       <- 1 - adimp - pctim

  for (i in 1:length(prcp)) {

    ### Set input precipitation and potential evapotranspiration
    pr = prcp[i] # This could be effective rainfall, a sum of rainfall and snowmelt
    edmnd = pet[i]

    ## Compute for different compnents...
    # ET(1), ET from Upper zone tension water storage
    et1 <- edmnd * uztwc/uztwm
    red <- edmnd - et1  # residual ET demand
    uztwc <- uztwc - et1

    # ET(2), ET from upper zone free water storage
    et2 <- 0
    #print(paste0("I=",i," uztwm= ",uztwm," uztwc= ",uztwc," et1= ", et1, " pr= ",pr," pet= ",edmnd))
    # in case et1 > uztws, no water in the upper tension water storage
    if (uztwc <= 0) {
      et1 <- et1 + uztwc #et1 = uztwc
      uztwc <- 0
      red <- edmnd - et1

      # when upper zone free water content is less than residual ET
      if (uzfwc < red) {

        # all content at upper zone free water zone will be gone as ET
        et2 <- uzfwc
        uzfwc <- 0
        red <- red - et2
        if (uztwc < thres_zero) uztwc <- 0
        if (uzfwc < thres_zero) uzfwc <- 0

        # when upper zone free water content is more than residual ET
      } else {
        et2 <- red  # all residual ET will be gone as ET
        uzfwc <- uzfwc - et2
        red <- 0
      }

      # in case et1 <= uztws, all maximum et (et1) are consumed at uztwc,
      # so no et from uzfwc (et2=0)
    } else {

      # There's possibility that upper zone free water ratio exceeds
      #upper zone tension water ratio. If so, free water is transferred to
      #tension water storage

      if((uztwc / uztwm) < (uzfwc / uzfwm)) {
        uzrat = (uztwc + uzfwc) / (uztwm + uzfwm)
        uztwc = uztwm * uzrat
        uzfwc = uzfwm * uzrat
      }

      if(uztwc < thres_zero) uztwc = 0
      if(uzfwc < thres_zero) uzfwc = 0

    }

    # ET(3), ET from Lower zone tension water storage when residual ET > 0
    et3 <- red * lztwc / (uztwm + lztwm) #residual ET is always bigger than ET(3)
    lztwc <- lztwc - et3

    # if lztwc is less than zero, et3 cannot exceed lztws
    if(lztwc < 0) {
      et3   <- et3 + lztwc  # et3 = lztwc
      lztwc <- 0
    }

    # Water resupply from Lower free water storages to Lower tension water storage
    saved  <- rserv * (lzfpm + lzfsm)
    ratlzt <- lztwc / lztwm
    ratlz  <- (lztwc + lzfpc + lzfsc - saved) / (lztwm + lzfpm + lzfsm - saved)

    # water is first taken from supplementary water storage for resupply
    if (ratlzt < ratlz) {

      del <- (ratlz - ratlzt) * lztwm
      lztwc <- lztwc + del  # Transfer water from lzfss to lztws
      lzfsc <- lzfsc - del

      # if tranfer exceeds lzfsc then remainder comes from lzfps
      if(lzfsc < 0) {
        lzfpc <- lzfpc + lzfsc
        lzfsc <- 0
      }
    }

    if(lztwc < thres_zero) {lztwc <- 0}
    # Comment for additional imprevious ET
    # # ET(5), ET from additional impervious (ADIMP) area
    # # ????? no idea where this come from, I think there's a possibility that et5 can be negative values
    et5   <- et1 + (red + et2) * (adimc - et1 - uztwc) / (uztwm + lztwm)
    adimc <- adimc - et5
    if(adimc < 0) {
      #et5 cannot exceed adimc
      et5 <- et5 + adimc # et5 = adimc
      adimc <- 0
    }
    et5 <- et5 * adimp

    # Time interval available moisture in excess of uztw requirements
    twx <- pr + uztwc - uztwm

    # all moisture held in uztw- no excess
    if(twx < 0) {
      uztwc <- uztwc + pr
      twx <- 0
      # moisture available in excess of uztw storage
    } else {
      uztwc = uztwm
    }
    #
    # for now twx is excess rainfall after filling the uztwc
    #
    adimc <- adimc + pr - twx

    # Compute Impervious Area Runoff
    roimp <- pr * pctim

    # Initialize time interval sums
    sbf   <- 0  # Sum of total baseflow(from primary and supplemental storages)
    ssur  <- 0  # Sum of surface runoff
    sif   <- 0  # Sum of interflow
    sperc <- 0  # Time interval summation of percolation
    sdro  <- 0  # Sum of direct runoff from the additional impervious area

    # Determine computational time increments for the basic time interval
    ninc <- floor(1.0 + 0.2*(uzfwc+twx))  # Number of time increments that interval is divided into for further soil-moisture accountng

    dinc <- 1.0 / ninc                    # Length of each increment in days
    pinc <- twx / ninc                    # Amount of available moisture for each increment

    # Compute free water depletion fractions for the time increment
    #(basic depletions are for one day)
    duz   <- 1 - (1 - uzk)^dinc
    dlzp  <- 1 - (1 - lzpk)^dinc
    dlzs  <- 1 - (1 - lzsk)^dinc
    # This is the version of Peter's
    # duz   <- uzk*dinc
    # dlzp  <- lzpk*dinc
    # dlzs  <- lzsk*dinc

    #print(paste0("ninc=", str(ninc)))

    # Start incremental for-loop for the time interval
    for (n in 1:ninc){

      adsur <- 0 # Amount of surface runoff. This will be updated.
      excess<- 0  # the excess of LZ soil water capacity

      # Compute direct runoff from adimp area
      ratio <- (adimc - uztwc) / lztwm
      if(ratio < 0) ratio <- 0

      # Amount of direct runoff from the additional impervious area
      addro <- pinc*(ratio^2)

      # Compute baseflow and keep track of time interval sum
      # Baseflow from free water primary storage
      bf_p <- lzfpc * dlzp
      lzfpc <- lzfpc - bf_p
      if(lzfpc <= 0.0001) {
        bf_p  <- bf_p + lzfpc
        lzfpc <- 0
      }

      sbf <- sbf + bf_p
      spbf<- sbf + bf_p
      # Baseflow from free water supplemental storage
      bf_s  <- lzfsc * dlzs
      lzfsc <- lzfsc - bf_s
      if (lzfsc <= 0.0001) {
        bf_s <- bf_s + lzfsc
        lzfsc <- 0
      }

      # Total Baseflow from primary and supplemental storages
      sbf <- sbf + bf_s

      # Compute PERCOLATION- if no water available then skip.
      if((pinc + uzfwc) <= 0.01) {
        uzfwc <- uzfwc + pinc
      } else {

        # Limiting drainage rate from the combined saturated lower zone storages
        percm <- lzfpm * dlzp + lzfsm * dlzs
        perc <- percm * uzfwc / uzfwm

        # DEFR is the lower zone moisture deficiency ratio
        defr <- 1.0 - (lztwc + lzfpc + lzfsc)/(lztwm + lzfpm + lzfsm)

        if(defr < 0) {defr <- 0}

        perc <- perc * (1.0 + zperc * (defr^rexp))

        # Note. . . percolation occurs from uzfws before pav is added

        # Percolation rate exceeds uzfws
        if(perc >= uzfwc) {perc <- uzfwc}

        uzfwc <- uzfwc - perc    # Percolation rate is less than uzfws.

        # Check to see if percolation exceeds lower zone deficiency.
        check <- lztwc + lzfpc + lzfsc + perc - lztwm - lzfpm - lzfsm
        if(check > 0) {
          perc <- perc - check
          uzfwc <- uzfwc + check
        }

        # SPERC is the time interval summation of PERC
        sperc <- sperc + perc

        # Compute interflow and keep track of time interval sum. Note that PINC has not yet been added.
        del <- uzfwc * duz # The amount of interflow

        ## Check whether interflow is larger than uzfwc
        if (del > uzfwc) {
          del<-uzfwc
          uzfwc<-0.0
        }else{
          uzfwc <- uzfwc - del
        }

        sif <- sif + del

        # Distribute percolated water into the lower zones. Tension water
        # must be filled first except for the PFREE area. PERCT is
        # percolation to tension water and PERCF is percolation going to
        # free water.

        perct <- perc * (1.0 - pfree)  # Percolation going to the tension water storage
        if((perct + lztwc) <= lztwm) {

          lztwc <- lztwc + perct
          percf <- 0 # Pecolation going to th lower zone free water storages

        } else {

          percf <- lztwc + perct - lztwm
          lztwc <- lztwm

        }

        # Distribute percolation in excess of tension requirements among the free water storages.
        percf <- percf + (perc * pfree)

        if(percf != 0) {

          # Relative size of the primary storage as compared with total lower zone free water storages.
          hpl <- lzfpm / (lzfpm + lzfsm)

          # Relative fullness of each storage.
          ratlp <- lzfpc / lzfpm
          ratls <- lzfsc / lzfsm

          # The fraction going to primary
          fracp <- hpl * 2 * (1 - ratlp) / (2 - ratlp - ratls)

          if(fracp > 1.0) {fracp <- 1.0}

          percp <- percf * fracp # Amount of the excess percolation going to primary
          percs <- percf - percp # Amount of the excess percolation going to supplemental
          lzfsc <- lzfsc + percs


          if(lzfsc > lzfsm) {
            percs <- percs - lzfsc + lzfsm
            lzfsc <- lzfsm
          }

          lzfpc <- lzfpc + percf - percs

          # This is different to Peter's
          #
          # Check to make sure lzfps does not exceed lzfpm
          if(lzfpc >= lzfpm) {
            excess <- lzfpc - lzfpm
            lztwc <- lztwc + excess
            lzfpc <- lzfpm
            if(lztwc >= lztwm) {
              excess <- lztwc - lztwm
              lztwc <- lztwm
            }
          }

        }

        #

        # Distribute PINC between uzfws and surface runoff
        if((pinc+excess) != 0) {

          # check if pinc exceeds uzfwm
          if((pinc + uzfwc+excess) <= uzfwm) {

            uzfwc <- uzfwc + pinc+excess  # no surface runoff
          } else {
            sur <- pinc + uzfwc + excess - uzfwm # Surface runoff
            uzfwc <- uzfwm

            ssur = ssur + (sur * parea)

            # ADSUR is the amount of surface runoff which comes from
            # that portion of adimp which is not currently generating
            # direct runoff. ADDRO/PINC is the fraction of adimp
            # currently generating direct runoff.
            adsur = sur * (1.0 - addro / pinc)
            ssur = ssur + adsur * adimp

          }
        }
      }

      adimc <- adimc + pinc - addro - adsur
      if(adimc > (uztwm + lztwm)) {
        addro = addro + adimc - (uztwm + lztwm)
        adimc = uztwm + lztwm
      }

      # Direct runoff from the additional impervious area
      sdro  = sdro + (addro * adimp)

      if(adimc < thres_zero) {adimc <- 0}

    } # END of incremental for loop

    # Compute sums and adjust runoff amounts by the area over which they are generated.

    # EUSED is the ET from PAREA which is 1.0 - adimp - pctim
    eused <- et1 + et2 + et3
    sif <- sif * parea

    # Separate channel component of baseflow from the non-channel component
    tbf <- sbf * parea   # TBF is the total baseflow
    bfcc <- tbf / (1 + side)    # BFCC is baseflow, channel component

    bfp = (spbf * parea) / (1.0 + side)
    bfs = bfcc - bfp
    if (bfs < 0.) bfs = 0
    bfncc = tbf - bfcc # BFNCC IS BASEFLOW, NON-CHANNEL COMPONENT

    # Ground flow and Surface flow
    base <- bfcc                       # Baseflow and Interflow are considered as Ground inflow to the channel
    surf <- roimp + sdro + ssur + sif  # Surface flow consists of Direct runoff and Surface inflow to the channel

    # ET(4)- ET from riparian vegetation.
    et4 <- (edmnd - eused) * riva  # no effect if riva is set to zero

    # Compute total evapotransporation - TET
    eused <- eused * parea
    tet <- eused + et4 + et5

    # Check that adimc >= uztws
    # This is not sure?
    #if(adimc > uztwc) adimc <- uztwc

    # Total inflow to channel for a timestep
    tot_outflow <- surf + base - et4;

    ### ------- Adjustments to prevent negative flows -------------------------#

    # If total outflow <0 surface and baseflow needs to be updated
    if (tot_outflow < 0) {

      tot_outflow = 0; surf = 0; base = 0;

    } else {

      surf_remainder = surf - et4
      surf <- max(0,surf_remainder)

      if (surf_remainder < 0) { # In this case, base is reduced

        base = base + surf_remainder
        if (base < 0) base = 0
      }
    }

    # Total inflow to channel for a timestep
    simaet[i]  <- tet
    simflow[i]  <- tot_outflow
    surf_tot[i] <- surf
    base_tot[i] <- base
    uztwc_ts[i] <- uztwc
    uzfwc_ts[i]  <- uzfwc
    lztwc_ts[i]  <- lztwc
    lzfpc_ts[i]  <- lzfpc
    lzfsc_ts[i]  <- lzfsc
  } #close time-loop
  return(data.frame("aetTot" = simaet,"flwTot" = simflow, "flwSurface" = surf_tot, "flwBase" = base_tot,
                    "uztwc"=uztwc_ts,"uzfwc"=uzfwc_ts,
                    "lztwc"=lztwc_ts,"lzfpc"=lzfpc_ts,"lzfsc"=lzfsc_ts))
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
  # WaSSI<-function(h){
  #
  #   # Using snowmelt function to calculate effective rainfall
  #
  #   snow.result<-snow_melt(ts.prcp = hru_prcp[,h],ts.temp = hru_tavg[,h],snowrange = c(-5,1))
  #   hru_rain <- snow.result$prcp
  #
  #   # PET using Hamon equation
  #   hru_pet <- hamon(par = par.petHamon[h], tavg = hru_tavg[,h], lat = hru_lat[h], jdate = jdate)
  #   hru_pet<-hru_pet*days_mon
  #   # Calculate hru flow for each elevation band
  #   ## calculate flow for each land cover type
  #
  #   # if lai data is not enough
  #   if(is.null(sim.dates[["Start_lai"]])) str.date.lai<-sim.dates[["Start_climate"]]
  #   if(is.null(sim.dates[["End_lai"]])) end.date.lai<-sim.dates[["End_climate"]]
  #
  #   ## set the first year lai to before
  #   if(sim.dates[["Start_lai"]]>sim.dates[["Start"]]-years(warmup)){
  #     lack_lai_years<-floor(as.numeric(difftime(sim.dates[["Start_lai"]],sim.dates[["Start"]]-years(warmup),units="days")/365))
  #     hru.lc.lai<-lapply(hru.lc.lai,function(x) rbind(apply(x[1:12,],2,rep,lack_lai_years),x))
  #     sim.dates[["Start_lai"]]<-sim.dates[["Start"]]-years(warmup)
  #   }
  #   ## set the last year lai to after
  #   if(sim.dates[["End_lai"]]<sim.dates[["End"]]){
  #     lack_lai_years<-floor(as.numeric(difftime(sim.dates[["End"]],sim.dates[["End_lai"]],units="days")/365))
  #     start1<-length(hru.lc.lai[[1]][,1])-11
  #     hru.lc.lai<-lapply(hru.lc.lai,function(x) rbind(x,apply(x[start1:(start1+11),],2,rep,lack_lai_years)))
  #     sim.dates[["End_lai"]]<-sim.dates[["End"]]
  #   }
  #
  #   # filter lai data by the simulation date
  #   lai_date <- seq.Date(sim.dates[["Start_lai"]], sim.dates[["End_lai"]], by = "month")
  #   grid_ind_lai <-which(lai_date %in% sim_date)
  #   hru.lc.lai<-lapply(hru.lc.lai,function(x) x[grid_ind_lai,] )
  #
  #   # Calculate PET based on LAI
  #   if(is.null(ET.coefs)){
  #     hru_lc_pet<-hru_pet*hru.lc.lai[[h]]*0.0222+0.174*hru_mrain+0.502*hru_pet+5.31*hru.lc.lai[[h]]
  #   }else{
  #     # Calculate the actual ET based on SUN's ET model
  #     hru_lc_pet<-matrix(NA,nrow =length(hru_pet), ncol=NLCs)
  #     for (i in c(1:NLCs)){
  #       hru_lc_pet[,i]<- ET.coefs[i,"Intercept"] +
  #         ET.coefs[i,"P_coef"]*hru_mrain+
  #         ET.coefs[i,"PET_coef"]*hru_pet+
  #         ET.coefs[i,"LAI_coef"]*hru.lc.lai[[h]][i]+
  #         ET.coefs[i,"P_PET_coef"]*hru_mrain*hru_pet +
  #         ET.coefs[i,"P_LAI_coef"]*hru_mrain*hru.lc.lai[[h]][i] +
  #         ET.coefs[i,"PET_LAI_coef"] *hru_pet*hru.lc.lai[[h]][i]
  #     }
  #   }
  #
  #   # combin all the input
  #   hru_in<-list(prcp=hru_prcp[,h],
  #                   temp = hru_tavg[,h],
  #                   rain=hru_rain,
  #                   snowpack=snow.result$snowpack,
  #                   snowmelt=snow.result$snowmelt,
  #                   PET_hamon=hru_pet)
  #
  #   # calculate flow based on PET and SAC-SMA model
  #   hru_lc_out <-apply(hru_lc_pet,2,sacSma_mon,par = par.sacsma[h,], prcp = hru_mrain )
  #
  #   if(!is.null(WUE.coefs)){
  #     # Calculate Carbon based on WUE for each vegetation type
  #     for (i in c(1:NLCs)){
  #       hru_lc_out[[i]][["GEP"]]<-hru_lc_out[[i]][["totaet"]] *WUE.coefs$WUE[i]
  #       hru_lc_out[[i]][["RECO"]]<-WUE.coefs$RECO_Intercept[i] + hru_lc_out[[i]][["GEP"]] *WUE.coefs$RECO_Slope[i]
  #       hru_lc_out[[i]][["NEE"]] <-hru_lc_out[[i]][["RECO"]]-hru_lc_out[[i]][["GEP"]]
  #     }
  #   }
  #   for (i in c(1:NLCs)){
  #     hru_lc_out[[i]][["LAI"]]<-hru.lc.lai[[h]][i]
  #     hru_lc_out[[i]][["PET"]]<-hru_lc_pet[,i]
  #   }
  #
  #   # Weight each land cover type based on it's ratio
  #   hru_out<-lapply(c(1:length(huc.lc.ratio[h,])),function(x) lapply(hru_lc_out[[x]],function(y) huc.lc.ratio[h,x]*y))
  #   out<-lapply(names(hru_lc_out[[1]]), function(var) apply(sapply(hru_out, function(x) x[[var]]),1,sum))
  #   names(out)<-names(hru_lc_out[[1]])
  #
  #   return(list(input=hru_in,out=out,lc_out=hru_lc_out))
  # }

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


# WaSSI for each HRU and calculate adjust temp and pet values
WaSSI<-function(hru,datain,sim.dates){
  require(dplyr)
  library(tidyr)

  # Set the input variables
  hru_prcp<-subset(datain$Climate,BasinID==hru)$Ppt_mm[sim.dates$climate_index]
  hru_tavg<-subset(datain$Climate,BasinID==hru)$Tavg_C[sim.dates$climate_index]
  NoLcs<-length(names(datain[["LAI"]]))-3
  hru_lat<-subset(datain$Cellinfo,BasinID==hru)$Latitude
  jdate<-sim.dates$jdate
  days_mon<-sim.dates$dateDays
  hru.lc.lai<-subset(datain$LAI,BasinID==hru)[sim.dates$lai_index,-c(1:3)]
  par.sacsma<-unlist(subset(datain$Soilinfo,BasinID==hru)[-c(1)])
  names(par.sacsma)<-toupper(names(par.sacsma))
  huc.lc.ratio<-unlist(subset(datain$Cellinfo,BasinID==hru)[grep("Lc",names(datain$Cellinfo))])

  ET.coefs<-datain$ET_coefs
  WUE.coefs<-datain$WUE_coefs

  # Using snowmelt function to calculate effective rainfall and snowmelt
  snow.result<-snow_melt(ts.prcp = hru_prcp,ts.temp = hru_tavg,snowrange = c(-5,1))
  hru_rain <- snow.result$prcp

  # PET using Hamon equation
  hru_pet <- hamon(par = 1, tavg = hru_tavg, lat = hru_lat, jdate = jdate)
  hru_pet<-hru_pet*days_mon

  # Calculate hru flow for each elevation band
  # Calculate PET_Sun based on LAI
  if(is.null(ET.coefs)){
    hru_lc_pet<-hru_pet*hru.lc.lai*0.0222+0.174*hru_rain+0.502*hru_pet+5.31*hru.lc.lai
  }else{
    # Calculate the PET_Sun based on user defined model
    hru_lc_pet<-matrix(NA,nrow =length(hru_pet), ncol=NoLcs)
    for (i in c(1:NoLcs)){
      hru_lc_pet[,i]<- ET.coefs[i,"Intercept"] +
        ET.coefs[i,"P_coef"]*hru_rain+
        ET.coefs[i,"PET_coef"]*hru_pet+
        ET.coefs[i,"LAI_coef"]*hru.lc.lai[[i]]+
        ET.coefs[i,"P_PET_coef"]*hru_rain*hru_pet +
        ET.coefs[i,"P_LAI_coef"]*hru_rain*hru.lc.lai[[i]] +
        ET.coefs[i,"PET_LAI_coef"] *hru_pet*hru.lc.lai[[i]]
    }
  }

  # Set the min PET of Hamon or Sun PET as the PET
  #hru_lc_pet<- apply(hru_lc_pet,2,function(x) apply(cbind(x,hru_pet),1,min))

  # calculate flow based on PET and SAC-SMA model
  hru_lc_out <-apply(hru_lc_pet,2,sacSma_mon,par = par.sacsma, prcp = hru_rain)

  # Calculate Carbon based on WUE for each vegetation type
  if(!is.null(WUE.coefs)){
    for (i in c(1:NoLcs)){
      hru_lc_out[["GEP"]][[i]]<-hru_lc_out[["totaet"]][[i]] *WUE.coefs$WUE[i]
      hru_lc_out[["RECO"]][[i]]<-WUE.coefs$RECO_Intercept[i] + hru_lc_out[[i]][["GEP"]] *WUE.coefs$RECO_Slope[i]
      hru_lc_out[["NEE"]][[i]] <-hru_lc_out[["RECO"]][[i]]-hru_lc_out[["GEP"]][[i]]
    }
  }
  # Add LAI and PET to each land cover and subset the target period
  for (i in c(1:NoLcs)){
    hru_lc_out[[i]][["LAI"]]<-hru.lc.lai[[i]]
    hru_lc_out[[i]][["PET"]]<-hru_lc_pet[,i]
    hru_lc_out[[i]][["Date"]]<-sim.dates$Seq_date
    hru_lc_out[[i]]<-hru_lc_out[[i]]%>%filter(Date>=sim.dates$Start) %>% dplyr::select(-Date)
  }

  # combin all the input
  hru_in<-data.frame(Date=sim.dates$Seq_date,
                     prcp=hru_prcp,
                     temp = hru_tavg,
                     rain=hru_rain,
                     snowpack=snow.result$snowpack,
                     snowmelt=snow.result$snowmelt,
                     PET_hamon=hru_pet)%>%
    filter(Date>=sim.dates$Start)%>%
    dplyr::select(-Date)

  # process data for each lc
  vars<-names(hru_lc_out[[1]])
  #vars<-vars[1:(length(vars)-1)]
  hru_lc_out<-lapply(vars, function(x) sapply(hru_lc_out, "[[", x))
  names(hru_lc_out)<-vars
  hru_out<-sapply(hru_lc_out, function(x) apply(x, 1, weighted.mean,huc.lc.ratio) )
  #colnames(hru_out)<-vars
  #hru_lc_out$Date<-hru_in$Date
  hru_in<-cbind(hru_in,hru_out)
  return(list(output=hru_in,lc_output=hru_lc_out))
}



librs<-c("dplyr","sf","zip","lubridate","raster","ggplot2","leaflet","rgdal","rgeos","cartography","leaflet.extras","parallel","shinyFiles","tidyr","reshape2")
f_lib_check(librs)

if(!exists("data_input")) data_input<-list()
if(!exists("BasinShp")) BasinShp<<-NULL
Ning<-"Ning Liu"
if(!exists("data_simulation")) data_simulation<-list()
if(!exists("Sim_dates")) Sim_dates<-list()
if(!exists("resultOutput")) resultOutput<-list()
mcores=detectCores()-1
#if(.Platform$OS.type=="windows") mcores=1
warmup<-2
