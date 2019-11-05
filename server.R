#
# This is the server logic of a dWaSSI model's web application. You can run the
# application by clicking 'Run App' above.
#

library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output,session) {
  options(shiny.maxRequestSize=30*1024^2)

  # Define the info printing
  inforprint<-reactiveValues(reading=NULL,
                            processing="Data process log:",
                            plotting="Plotting log: ",
                            simulating="Simulationg log: ")
  f_addinfo<-function(tab,content=NULL){
    inforprint[[tab]]<<-paste(inforprint[[tab]],content,sep="\n")
  }

  # Tab: Upload data----
  ## Action: Print reading infor ----
  output$printreadinginfo<-renderPrint({
    if(!is.null(inforprint$reading)) print(inforprint$reading)
  })
  if (!exists("data_input")) data_input<-list()
  ## Action: Print the uploaded file names ----
  observeEvent(input$readdata,{
    print(input$Input_climate)
    if(is.null(c(input$Input_climate,input$Input_LAI,input$Input_cellinfo,input$Input_Soilinfo))) {
      inforprint$reading<-"!!! Warining: You need to select at least one file!"
      return()
    }
    inforprint$reading<-"Reading input data ..."
    if (!is.null(input$Input_climate)) data_input[["Climate"]]<<-read.csv(input$Input_climate$datapath)
    if (!is.null(input$Input_LAI)) data_input[["LAI"]]<<-read.csv(input$Input_LAI$datapath)
    if (!is.null(input$Input_cellinfo)) data_input[["Cellinfo"]]<<-read.csv(input$Input_cellinfo$datapath)
    if (!is.null(input$Input_Soilinfo)) data_input[["Soilinfo"]]<<-read.csv(input$Input_Soilinfo$datapath)
    inforprint$reading<-"Finished reading all the input data ..."
    output$printsummary<-renderPrint({
      for (i in 1:length(data_input)){
        var<-names(data_input)[i]
        df<-data_input[[i]]
        print(paste0(var," data----"))
        print(head(df))
        print(summary(df))
      }
    })
  })

  # Tab: Process data----

# # Select file by name
#   filepaths<-reactiveValues(laifname=NULL)
#   volumes = getVolumes()
#   observe({
#     shinyFileChoose(input, "Input_lai_fpath", roots = volumes, session = session)
#
#     if(!is.null(input$Input_lai_fpath)){
#       # browser()
#       file_selected<-parseFilePaths(volumes, input$Input_lai_fpath)
#       filepaths$laifname <- renderText(as.character(file_selected$datapath))
#     }
#   })

#print(filepaths$laifname)
  ## Print: Print the processig infor ----
  output$printprocessinginfo<-renderPrint({
    if(!is.null(inforprint$processing)) writeLines(inforprint$processing)
  })

  ## Action: Read Basin shapefile ----
  if (!exists("BasinShp")) BasinShp<-NULL
  observeEvent(input$Input_basin,{
    inforprint$processing<-"Processing log: "
    f_addinfo("processing","Reading shapfile ...")
    Basins<- readOGR(input$Input_basin$datapath)
    Basins$BasinID<-c(1:length(Basins[,1]))

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
       Basins_pro<-spTransform(Basins,CRS("+init=epsg:3857"))
       Basins[["Area_m2"]]=round(gArea(Basins_pro,byid = T),2)
     }else{
       Basins[["Area_m2"]]=round(gArea(Basins,byid = T),2)
     }

    }
    BasinShp<<-Basins
    f_addinfo("processing","Finished processing the uploaded shapefile.")
  })


  ## Action: Read raster data  ----
  if (!exists("input_rasters")) input_rasters<-list()
  if (!exists("data_input")) data_input<-list()

  observeEvent(input$processrasters,{

    inforprint$processing<-"Processing log: "
    if (is.null(BasinShp)) {
      f_addinfo("processing","You need to upload a shapefile including the basins.")
      return()
    }

    # process Climate input
    if(is.null(input$Input_temp_raster)| is.null(input$Input_precp_raster)){
      f_addinfo("processing","Warning: You need to provide the climate raster time series!")

    }else if("Climate" %in% names(data_input) & ! input$updateClimate){
      f_addinfo("processing","Warning: Climate data has been processed!")

    }else{

      f_addinfo("processing","Processing climate data ...")
      Tmean_catchment<-f_sta_shp_nc(input$Input_temp_raster$datapath,BasinShp,varname = "Tavg_C",yr.start = input$yearStartClimate,zonal_field = "BasinID")
      Pre_catchment<-f_sta_shp_nc(input$Input_precp_raster$datapath,BasinShp,varname = "Ppt_mm",yr.start = input$yearStartClimate,zonal_field = "BasinID")
      climate<-Pre_catchment%>%
        mutate(Tavg_C=Tmean_catchment$Tavg_C)%>%
        arrange(BasinID,Year,Month)%>%
        mutate(Ppt_mm=round(Ppt_mm,3))%>%
        mutate(Tavg_C=round(Tavg_C,3))%>%
        dplyr::select(BasinID,Year,Month,Ppt_mm,Tavg_C)
      data_input[["Climate"]]<<-climate
      f_addinfo("processing","Finished processing climate input!")
      print("Finished processsing Climate data ...")
      #print(summary(climate))
    }

    # Process cellinfo
    if(is.null(input$Input_lc_raster)){
      f_addinfo("processing","Warning: There is no land cover input!")

    }else if("Cellinfo" %in% names(data_input) & ! input$updateLc){
      f_addinfo("processing","Warning: Cellinfo data has been processed!")

    }else{
      f_addinfo("processing","processing Cellinfo data ...")
      print(paste0("processing Cellinfo ..."))
      cellinfo<-f_cellinfo(classfname=input$Input_lc_raster$datapath,
                           Basins=BasinShp,
                           byfield="BasinID",
                           demfname=input$Input_dem_raster$datapath)
      data_input[["Cellinfo"]]<<-cellinfo
      print("Finished processsing Cellinfo data ...")
      f_addinfo("processing","Finished processing cellinfo!")
    }
      # process LAI input
    if(is.null(input$Input_lc_raster) | (is.null(input$Input_lai_raster) & input$Input_lai_fpath=="~" )){
      f_addinfo("processing","Warning: There is no land cover or LAI data!")

    }else if("LAI" %in% names(data_input) & ! input$updateLai){
      f_addinfo("processing","Warning: LAI data has been processed!")

    }else{
      print("Processing LAIinfo")
      laifname<-input$Input_lai_fpath
      if(laifname=="~")  laifname<-input$Input_lai_raster$datapath
      laiinfo<-f_landlai(lcfname=input$Input_lc_raster$datapath,
                          laifname=laifname,
                          Basins=BasinShp,
                          byfield="BasinID",
                          yr.start=input$yearStartLai)
      data_input[["LAI"]]<<-laiinfo
      print("Finished processsing LAI data ...")
      f_addinfo("processing","Finished processing Lai data!")
    }
      # process SOIL input
    if(is.null(input$Input_soil_raster)){
      f_addinfo("processing","Warning: There is no soil data!")

      }else if("Soilinfo" %in% names(data_input) & ! input$updateSoil){
      f_addinfo("processing","Warning: Soilinfo data has been processed!")

      }else{
      print("Processing Soilinfo")
      soilinfo<-f_soilinfo(soilfname=input$Input_soil_raster$datapath,
                          Basins=BasinShp)
      data_input[["Soilinfo"]]<<-soilinfo
      print("Finished processsing LAI data ...")
      f_addinfo("processing","Finished processing Soil data!")
    }
    # Print the summary of the input
      output$prntraster<-renderPrint({
        if(length(data_input)<1) return("No data!")
        for (i in 1:length(data_input)){
          var<-names(data_input)[i]
          df<-data_input[[i]]
          print(paste0(var," data----"))
          print(head(df))
          print(summary(df))
        }
      })
  })


  ## PlotAction: Plot Basin and raster data  ----
      observeEvent(input$plotrasterdata,{

        Tmean_avg<-NULL
        Pre_avg<-NULL
        Lai_avg<-NULL
        dem<-NULL

        if(!is.null(input$Input_dem_raster$datapath)){
          print("reading dem")
          dem<-raster(input$Input_dem_raster$datapath)
          f_addinfo("processing","Finished reading dem data!")
        }
        if(!is.null(input$Input_lc_raster$datapath)){
          print("reading lc")
          lc<-raster(input$Input_lc_raster$datapath)
          f_addinfo("processing","Finished reading ladn cover data!")
        }
        if(!is.null(input$Input_temp_raster$datapath)){
          print("reading temperature data")
          f_addinfo("processing","Calculate average for temperature time series!")
          Tmean_br<-brick(input$Input_temp_raster$datapath)
          beginCluster()
          Tmean_avg <- clusterR(Tmean_br, calc, args=list(fun=mean,na.rm=T))
          endCluster()
          }
        if(!is.null(input$Input_precp_raster$datapath)){
          print("reading precipitation data")
          Pre_br<-brick(input$Input_precp_raster$datapath)
          beginCluster()
          Pre_br <- clusterR(Pre_br, calc, args=list(fun=mean,na.rm=T))
          Pre_avg<-Pre_br*12
          endCluster()
        }

        if(!is.null(input$Input_lai_raster$datapath)){
          print("reading LAI data")
          Lai_br<-brick(input$Input_lai_raster$datapath)
          beginCluster()
          Lai_avg <- clusterR(Lai_br, calc, args=list(fun=mean,na.rm=T))
          endCluster()
        }


      output$basinrastermap <- renderLeaflet({
        input$plotrasterdata
        # Plot the BasinShp
        print("plotting basic map")
        input_leaflet<-leaflet() %>%
            addTiles(group = "OSM (default)") %>%
            addProviderTiles('Esri.WorldImagery',group = "Esri.Imagery")%>%
            addFullscreenControl()
        grps<-c("OSM (default)", "Esri.Imagery")
        ovlgrps<-NULL

        if (!is.null(BasinShp)) {
          print("Add basin to basic map")
          popups<-paste0("BasinID = ",BasinShp$BasinID)
          if("Elev_m" %in% names(BasinShp)) popups<-f_paste(popups,paste0("; Elevatation = ",round(BasinShp$Elev_m,0),"m"))

          ovlgrps<-c(ovlgrps,"Watershed")
          input_leaflet<-input_leaflet%>%
            addPolygons(data=BasinShp,weight=1,col = 'red',fillOpacity = 0.2,
                        highlight = highlightOptions(color='white',weight=1,
                                                     bringToFront = TRUE),
                        group = "Watershed")%>%
            addMarkers(lng = BasinShp$Longitude, lat = BasinShp$Latitude,
                       popup = popups,
                       label = paste0("BasinID = ",BasinShp$BasinID),
                       clusterOptions = markerClusterOptions())

        }

        if(!is.null(Tmean_avg)){
          print("Add temperature to basic map")
          ovlgrps<-c(ovlgrps,"Temperature")
          pal_tmp <- colorNumeric(c("cyan", "yellow", "red"), values(Tmean_avg),
                            na.color = "transparent")
          print(Tmean_avg)
          input_leaflet<-input_leaflet%>%
            addRasterImage(Tmean_avg, colors = pal_tmp, opacity = 0.8, group = "Temperature") %>%
            addLegend("bottomleft",pal = pal_tmp, values = values(Tmean_avg),
                      title = "Mean temp (°C)", group = "Temperature")
        }

        if(!is.null(Pre_avg)){
          print("Add precipitation to basic map")
          ovlgrps<-c(ovlgrps,"Precipitation")
        pal_prcp <- colorNumeric(c("azure", "cornflowerblue","darkblue"), values(Pre_avg),
                                na.color = "transparent")
        input_leaflet<-input_leaflet%>%
          addRasterImage(Pre_avg, colors = pal_prcp, opacity = 0.8, group = "Precipitation") %>%
          addLegend("bottomleft",pal = pal_prcp, values = values(Pre_avg),
                    title = "Annual precipation (mm/yr)", group = "Precipitation")
        }

        if(!is.null(dem)){
          print("Add dem to basic map")
          ovlgrps<-c(ovlgrps,"Elevation")
          pal_dem <- colorNumeric(c("green","yellow","red"), values(dem),
                                  na.color = "transparent")

          input_leaflet<-input_leaflet%>%
            addRasterImage(dem, colors = pal_dem, opacity = 0.8, group = "Elevation") %>%
            addLegend("bottomleft",pal = pal_dem, values = values(dem),
                    title = "Elevation (m)", group = "Elevation")
        }

        if(!is.null(Lai_avg)){
          print("Add lai to basic map")
          ovlgrps<-c(ovlgrps,"LAI")
          pal_lai <- colorNumeric(c("red","yellow","green"), values(Lai_avg),
                                  na.color = "transparent")

          input_leaflet<-input_leaflet%>%
            addRasterImage(Lai_avg, colors = pal_lai, opacity = 0.8, group = "LAI") %>%
            addLegend("bottomleft",pal = pal_lai, values = values(Lai_avg),
                      title = "Leaf area index (m2/m2)", group = "LAI")
        }

        if(!is.null(lc)){
          print("Add land cover to basic map")
          ovlgrps<-c(ovlgrps,"Land cover")
          pal_lc <- colorFactor("RdYlBu", levels = values(lc),values(lc),
                                  na.color = "transparent")

          input_leaflet<-input_leaflet%>%
            addRasterImage(lc, colors = pal_lc, opacity = 0.8, group = "Land cover") %>%
            addLegend("bottomleft",pal = pal_lc, values =  values(lc),
                      title = "Land cover",group = "Land cover")
        }

        # Layers control
        input_leaflet %>%
          addLayersControl(
            baseGroups = grps,
            overlayGroups = ovlgrps,
            options = layersControlOptions(collapsed = FALSE)
          )

      })

  })

  # Tab: Plot input data----
  ## Print: Print the processig infor ----
  output$printplottinginfo<-renderPrint({
    if(!is.null(inforprint$plotting)) writeLines(inforprint$plotting)
  })

  ## Plot: Plot selected input data ----
  observeEvent(input$plotdata,{
        inforprint$plotting<-"Plotting log: "

        data2plot<-input$daname2plot
        basinids<-as.integer(strsplit(input$plotBasinID,",")[[1]])
        plotvars<-strsplit(input$plotvar,",")[[1]]
        plotvarname<-strsplit(input$varnames,",")[[1]]
       # plotbasinname<-strsplit(input$Basinnames,",")[[1]]
        yr.start<-input$plotyrrange[1]; yr.end<-input$plotyrrange[2]
        plotmonths<-input$plotmonths

        print(paste0("Print data ....",data2plot))
        f_addinfo("plotting",paste0("Printing ",data2plot," data ..."))
        if(!data2plot %in% names(data_input)) {
          f_addinfo("plotting","It does have this data.")
          return()}

        if(sum(basinids %in% data_input[[data2plot]]$BasinID)<1) {

          f_addinfo("plotting",paste0("!!! Warning: No BasinID: ",paste(basinids,collapse=",")))
          return()
        }
        f_addinfo("plotting",paste0("Printed BasinID: ",paste(basinids,collapse=",")))
        df<-data_input[[data2plot]]%>%
            filter(BasinID %in% basinids)%>%
            mutate(BasinID=factor(paste0("BasinID = ",BasinID)))%>%
            filter(Year>=yr.start & Year<=yr.end)%>%
            filter(Month %in% plotmonths)
        if(nrow(df)<1) {
          f_addinfo("plotting","!!! Warning: No data in this time period!")
          f_addinfo("plotting",paste0("You can plot priod of ",paste(range(data_input[[data2plot]]$Year),collapse = " - ")))
          return()
          }else{
          df<-df%>%
              mutate(Date=as.Date(paste0(Year,"-",Month,"-","01")))
          }

        # test for LAI data
        if(data2plot=="LAI" & sum(plotvars %in% names(data_input[[data2plot]]))<1) {
          f_addinfo("plotting",paste0("!!! Warning: Could not find input variables: ",paste(plotvars,collapse=",")))
          f_addinfo("plotting",paste0("You can select these variables from LAI data: ",paste(names(data_input[[data2plot]])[-c(1:3)],collapse=",")))
          return()
        }else if (data2plot=="LAI" & sum(plotvars %in% names(data_input[[data2plot]]))>=1){

          if(length(plotvarname)!=length(plotvars)) f_addinfo("plotting","!!! Warning: You need to provide the same number of names as the variable.")
          var_indx<-which(!plotvars %in% names(data_input[[data2plot]]))
          if(sum(var_indx)>0) {
            f_addinfo("plotting",paste0("!!! Warning: These variables are not included: ",paste(plotvars[var_indx],collapse=",")))
            plotvars<-plotvars[-var_indx]
            plotvarname<-plotvarname[-var_indx]

          }
          f_addinfo("plotting",paste0("Printed variables: ",paste(plotvars,collapse=",")))

        }

    output$Plotinput <- renderPlot({
        input$plotdata
        if(data2plot == "Climate"){

          p1<-ggplot(df,aes(x=Date,y=Ppt_mm))+geom_line(col="blue")+geom_point(col="blue")+
            ggtitle(paste0("Monthly Precipitation (mm/month)"))+
            facet_wrap(BasinID~.,ncol=2)+
            scale_x_date(date_breaks = "1 year", date_labels = "%Y")+
            theme_ning(size.axis = 12,size.title = 14)
          p2<-ggplot(df,aes(x=Date,y=Tavg_C))+geom_line(col="forestgreen")+geom_point(col="forestgreen")+
            ggtitle(paste0("Monthly mean temperature (C)"))+
            facet_wrap(BasinID~.,ncol=2)+
            scale_x_date(date_breaks = "1 year", date_labels = "%Y")+
            theme_ning(size.axis = 12,size.title = 14)

          multiplot(p1,p2,cols = 1)

        }else if(data2plot == "LAI"){

          df<-df%>%
            dplyr::select(one_of(c("BasinID","Date"),plotvars))%>%
            gather(Lcs,LAI,plotvars)
          if(length(plotvarname)==length(plotvars)) df$Lcs<-factor(df$Lcs,levels =plotvars,labels = plotvarname)

          df%>%
            ggplot(aes(x=Date,y=LAI,col=Lcs,shape=Lcs))+geom_line()+geom_point(size=2.5)+
            ggtitle(paste0("Monthly Leaf area index (m2/m2)"))+
            facet_wrap(BasinID~.,ncol=2)+
            labs(col = "Land cover",shape="Land cover")+
            scale_x_date(date_breaks = "1 year", date_labels = "%Y")+
            theme_ning(size.axis = 12,size.title = 14)

        }

      })
  })
  # Tab: simulation----
  ## Print: Print the simulation infor ----
  output$printsimulatinginfo<-renderPrint({

    if(!is.null(inforprint$simulating)) writeLines(inforprint$simulating)
  })

  ## Action: Subset data by Station ID and period----
  observeEvent(input$runSimulation, {
    inforprint$simulating<-"Processing log:"
    Sim_dates[["Start"]]<<-as.Date(paste0(format(input$dateSimulation[1],"%Y"),"-01-01"))
    Sim_dates[["End"]]<<-as.Date(paste0(format(input$dateSimulation[2],"%Y"),"-12-01"))
    Sim_dates[["Seq_date"]]<<-seq.Date(Sim_dates[["Start"]],Sim_dates[["End"]],by = "month")

    StationID<-as.integer(input$StationID)
    if(StationID==0) {
      BasinID_sel <- data_input[["Cellinfo"]]$BasinID
      f_addinfo("simulating","No Station is provided, so that all HUCs are used in the simulation!")
    }else if(is.null(input$Input_routpar)){
      f_addinfo("simulating","!!! ERROR: Flow routing file is not provided!")
      return()
    }else{
      f_addinfo("simulating","Only upstreams of the station are selected for the simulation!")
       routpar<-read.csv(input$Input_routpar$datapath)
       BasinID_sel<-f_upstreamHUCs(BasinID =StationID ,routpar=routpar)
       BasinID_sel<-c(BasinID_sel,StationID)
      }

    NoBasins<-length(BasinID_sel)

    # process soilinfo
print("process soil")
    par_sacsma<-data_input[["Soilinfo"]]%>%
      as.data.frame()%>%
      filter(BasinID %in% BasinID_sel)%>%
      dplyr::select(-BasinID)%>%
      as.matrix()
    colnames(par_sacsma)<-toupper(colnames(par_sacsma))
    data_simulation[["par_sacsma"]]<<-par_sacsma
    f_addinfo("simulating","Finished subsetting soilinfo")

    # process climate data
    climate_sel<-data_input[["Climate"]]%>%
      filter(BasinID %in% BasinID_sel)%>%
      mutate(Date=as.Date(paste0(Year,"-",Month,"-","01")))%>%
      filter(Date>=Sim_dates[["Start"]] & Date<=Sim_dates[["End"]])
    if(nrow(climate_sel)<(12*NoBasins)){
      f_addinfo("simulating",paste0("!!! ERROR: Not enough climate data
      in the selected period, you can select Climate data
      for the period of ",paste(range(data_input[["Climate"]]$Year),collapse="-")))
      return()
      }
    # Climate time range
    Sim_dates[["Start_climate"]] <<- as.Date(paste0(min(data_input[["Climate"]]$Year),"/01/01"))
    Sim_dates[["End_climate"]] <<- as.Date(paste0(max(data_input[["Climate"]]$Year),"/12/01"))
    Sim_dates[["Seq_date_climate"]]<<-seq.Date(Sim_dates[["Start_climate"]],Sim_dates[["End_climate"]],by = "month")
    Sim_dates[["Sim_ind"]] <<-which(Sim_dates[["Seq_date_climate"]] %in% Sim_dates[["Seq_date"]])

    clim_prcp<-matrix(data_input[["Climate"]]$Ppt_mm,ncol=NoBasins)
    data_simulation[["clim_prcp"]]<<-clim_prcp
    clim_tavg<-matrix(data_input[["Climate"]]$Tavg_C,ncol=NoBasins)
    data_simulation[["clim_tavg"]]<<-clim_tavg

    f_addinfo("simulating","Finished subsetting Climate")

    ## Input LAI and ratio for each vegetation type
    huc_lc_ratio<-data_input[["Cellinfo"]]%>%
      filter(BasinID %in% BasinID_sel)%>%
      dplyr::select(starts_with("Lc_"))%>%
      as.matrix()
    # mutate(Date=as.Date(paste0(Year,"-",Month,"-","01")))%>%
    #   filter(Date>=SimStart & Date<=SimEnd)
    hru_lc_lai<-lapply(BasinID_sel, function(x) as.matrix(subset(data_input[["LAI"]],BasinID==x)[c(4:length(data_input[["LAI"]][1,]))]))
    # LAI time range
    Sim_dates[["Start_lai"]] <<- as.Date(paste0(min(data_input[["LAI"]]$Year),"/01/01") )
    Sim_dates[["End_lai"]] <<- as.Date(paste0(max(data_input[["LAI"]]$Year),"/12/01"))

    data_simulation[["huc_lc_ratio"]]<<-huc_lc_ratio
    data_simulation[["hru_lc_lai"]]<<-hru_lc_lai
     f_addinfo("simulating","Finished subsetting laiinfo")

     ## Hru info like Lat, long and elevation for each HUC
    hru_info<-data_input[["Cellinfo"]]%>%
      filter(BasinID %in% BasinID_sel)%>%
      dplyr::select(one_of("BasinID","Latitude","Area_m2","Elev_m","FlwLen_m"))
    if(!"Elev_m" %in% names(hru_info)) hru_info["Elev_m"]<-1000
    if(!"FlwLen_m" %in% names(hru_info)) hru_info["FlwLen_m"]<-1000
    data_simulation[["hru_info"]]<<-as.matrix(hru_info)
    f_addinfo("simulating","Finished processing cellinfo")

    # Load coefficients for SUN's ET and carbon calculation
    ET_coefs<-NULL
    WUE_coefs<-NULL
    if(! is.null(input$Input_ETmodel)) {
      ET_coefs<-read.csv(input$Input_ETmodel$datapath,nrows = length(huc_lc_ratio[1,]))
      names(ET_coefs)<-c("LC_ID","Intercept","P_coef","PET_coef","LAI_coef","P_PET_coef","P_LAI_coef","PET_LAI_coef","IGBP","LC_Name")
      data_simulation[["ET_coefs"]]<<-ET_coefs
      f_addinfo("simulating","User defined ET model is read.")
      }
    if(! is.null(input$Input_WUEmodel)) {
      WUE_coefs<-read.csv(input$Input_ETmodel$datapath,nrows = length(huc_lc_ratio[1,]))
      names(WUE_coefs)<-c("LC_ID","WUE","RECO_Intercept", "RECO_Slope","IGBP","LC_Name")
      data_simulation[["WUE_coefs"]]<<-WUE_coefs
      f_addinfo("simulating","User defined WUE model is read.")
    }
    ## PET parameters
    par_petHamon<- rep(1,NoBasins)
    print(par_petHamon)
    # Run the model
    print(" runing")
    f_addinfo("simulating","Runing simulation ...")
    flowR <<- dWaSSIC(sim.dates=Sim_dates, warmup = 2,mcores = 1,
                     par.sacsma = par_sacsma,par.petHamon=par_petHamon,par.routing =par_routing,
                     hru.info = hru_info,
                     clim.prcp = clim_prcp, clim.tavg = clim_tavg,
                     hru.lai=NULL,hru.lc.lai=hru_lc_lai,huc.lc.ratio=huc_lc_ratio,
                     WUE.coefs = WUE_coefs,ET.coefs = ET_coefs)

    print("finished runing")

})

  ## Action: Run model ----
  observeEvent(input$plotSimOut,{
    out_weight<-lapply(c(1:length(flowR$HUC)),function (x) flowR$HUC[[x]]*data_simulation[["hru_info"]][,"Area_m2"][x]/sum(data_simulation[["hru_info"]][,"Area_m2"]))
    flowR_Catchment<-as.data.frame(sapply(names(out_weight[[1]]), function(var) apply(sapply(out_weight, function(x) x[[var]]),1,sum)))

    Rain_mean<-apply(data_simulation[["clim_prcp"]], 1, mean,na.rm=T)
    df <- data.frame(date = Sim_dates$Seq_date,Srufflow = flowR$FLOW_SURF,
                     Baseflow=flowR$FLOW_BASE,Totflow=flowR$FLOW_SURF+flowR$FLOW_BASE,
                     AET=flowR_Catchment$totaet,Rain=Rain_mean[Sim_dates$Sim_ind],
                     nt=flowR_Catchment$tot)%>%
      filter(date>=as.Date(paste0(input$plotSimDaterange[1],"/01/01")) & date<=as.Date(paste0(input$plotSimDaterange[2],"/12/01")) )

    output$plotSimOut <- renderPlot({
      input$plotSimOut
      df%>%
        gather(Variable,Value,Srufflow:nt)%>%
        filter(Variable %in% c("Srufflow","Baseflow","Totflow","Rain"))%>%
        mutate(Variable=factor(Variable,levels = c("Srufflow","Baseflow","Totflow" , "AET","Rain","nt" ),labels = c("Surface flow","Base flow","Total flow","Actual evpotransipration","Rainfall","nn")))%>%
        ggplot(aes(x=date,y=Value,col=Variable))+
        labs(y="(mm/month)")+
        scale_x_date(date_breaks = "1 year", date_labels = "%Y")+
        geom_line()+geom_point()+
        theme_ning(size.axis = 12,size.title =14 )
    })
  })

  ## Action: Plot the basic result----
  observeEvent(input$plotresa,{

    output$plotresout<-renderPlot({
      req(input$plotvars)
      load("www/outp.RData")
      out<-merge(outp$Input,outp$Out)
      ind<-which(index(out)>= as.POSIXct(paste0(input$plotdaterange[1],"-01-01")) & index(out)<= as.POSIXct(paste0(input$plotdaterange[2],"-12-31")))
      plot(out[ind,input$plotvars],main = "",xlab = "Time")
      },
      height=1000,
      res = 100)
  })

  ## Action: select time period for calibration and validation----



  } # END
)
