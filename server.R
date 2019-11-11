#
# This is the server logic of a dWaSSI model's web application. You can run the
# application by clicking 'Run App' above.
#

library(shiny)
if(!exists("resultOutput")) resultOutput<-list()
# Define server logic required to draw a histogram
shinyServer(function(input, output,session) {
  options(shiny.maxRequestSize=100*1024^2)

  # Define the info printing
  inforprint<-reactiveValues(reading=NULL,
                            processing="Data process log:",
                            plotting="Plotting log: ",
                            simulating="Simulationg log: ",
                            simulplotting="Plotting log: ")
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
    if(is.null(c(input$Input_climate,input$Input_LAI,input$Input_cellinfo,input$Input_soilinfo))) {
      inforprint$reading<-"!!! Warining: You need to select at least one file!"
      return()
    }
    inforprint$reading<-"Reading input data ..."
    if (!is.null(input$Input_climate)) data_input[["Climate"]]<<-read.csv(input$Input_climate$datapath)
    if (!is.null(input$Input_LAI)) data_input[["LAI"]]<<-read.csv(input$Input_LAI$datapath)
    if (!is.null(input$Input_cellinfo)) data_input[["Cellinfo"]]<<-read.csv(input$Input_cellinfo$datapath)
    if (!is.null(input$Input_soilinfo)) data_input[["Soilinfo"]]<<-read.csv(input$Input_soilinfo$datapath)
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
       Basins_pro<-spTransform(Basins,CRS("+init=epsg:32648"))
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

    withProgress(message = 'Processing ', value = 0, {

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

      incProgress(0.1, detail = paste("Climate data ...!"))

      f_addinfo("processing","Processing climate data ...")
      print("Processsing Climate data ...")
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

    incProgress(0.1, detail = paste("Cellinfo ...!"))
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
      incProgress(0.2, detail = paste("Lai data ...!"))
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
      incProgress(0.5, detail = paste("Soilinfo ...!"))
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
  })


  ## PlotAction: Plot Basin and raster data  ----
      observeEvent(input$plotrasterdata,{

        Tmean_avg<-NULL
        Pre_avg<-NULL
        Lai_avg<-NULL
        dem<-NULL
        lc<-NULL

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


  ## SaveAction: save the processed input data  ----
  observeEvent(input$saveInputData,{

   for (var in  names(data_input)){

      write.csv(data_input[[var]],paste0("www/inputs/",var,".csv"),row.names=F)

   }

    f_addinfo("processing",
              paste0("The input data (",
                     paste(names(data_input),collapse = ","),
                     ") has been save to default foler 'www/inputs'"))

    })

  ## ExportAction: save the processed input data  ----
    output$downloadInputData <- downloadHandler(

      # This function returns a string which tells the client
      # browser what name to use when saving the file.
      filename = function() {

        paste("Inputs", "zip", sep=".")
      },

      # This function should write data to a file given to it by
      # the argument 'file'.
      content = function(fname) {

        fs=dir("www/inputs/","csv",full.names = T)
        withProgress(message = 'Compressing files ...', value = 0.5,{
        zip::zipr(zipfile=fname,files=fs)
        })

      },
      contentType="application/zip"
    )

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

        if(data2plot %in% c("Cellinfo","Soilinfo")){
          df<-data_input[[data2plot]]%>%
            filter(BasinID %in% basinids)
        }else{

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
        }
    output$printbasininfo<-renderPrint({

      if(data2plot %in% c("Cellinfo","Soilinfo")) df
    })

    output$Plotinput <- renderPlot({
        input$plotdata
        if(data2plot == "Climate"){

          if(input$plotannualinput){

            # plot annual result
            df<-df%>%
              dplyr::select(-Date,-Month)%>%
              group_by(BasinID,Year)%>%
              summarise(Ppt_mm=sum(Ppt_mm),Tavg_C=mean(Tavg_C))
            p1<-ggplot(df,aes(x=Year,y=Ppt_mm))+geom_line(col="blue")+geom_point(col="blue")+
              ggtitle(paste0("Annual Precipitation (mm/yr)"))+
              facet_wrap(BasinID~.,ncol=2)+
              scale_x_continuous(breaks = c(1982:2018))+
              theme_ning(size.axis = 12,size.title = 14)
            p2<-ggplot(df,aes(x=Year,y=Tavg_C))+geom_line(col="forestgreen")+geom_point(col="forestgreen")+
              ggtitle(paste0("Annual mean temperature (C)"))+
              facet_wrap(BasinID~.,ncol=2)+
              scale_x_continuous(breaks = c(1982:2018))+
              theme_ning(size.axis = 12,size.title = 14)

          }else{
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
          }


          multiplot(p1,p2,cols = 1)

        }else if(data2plot == "LAI"){

          df<-df%>%
            dplyr::select(one_of(c("BasinID","Date"),plotvars))%>%
            gather(Lcs,LAI,plotvars)
          if(length(plotvarname)==length(plotvars)) df$Lcs<-factor(df$Lcs,levels =plotvars,labels = plotvarname)

          if(input$plotannualinput){
            df%>%
              mutate(Year=as.integer(format(Date,"%Y")))%>%
              dplyr::select(-Date)%>%
              group_by(BasinID,Lcs,Year)%>%
              summarise(LAI=mean(LAI))%>%
              ggplot(aes(x=Year,y=LAI,col=Lcs,shape=Lcs))+geom_line()+geom_point(size=2.5)+
              ggtitle(paste0("Annual mean Leaf area index (m2/m2)"))+
              facet_wrap(BasinID~.,ncol=2)+
              labs(col = "Land cover",shape="Land cover")+
              scale_x_continuous(breaks = c(1982:2018))+
              theme_ning(size.axis = 12,size.title = 14)
            }else{

              df%>%
                ggplot(aes(x=Date,y=LAI,col=Lcs,shape=Lcs))+geom_line()+geom_point(size=2.5)+
                ggtitle(paste0("Monthly Leaf area index (m2/m2)"))+
                facet_wrap(BasinID~.,ncol=2)+
                labs(col = "Land cover",shape="Land cover")+
                scale_x_date(date_breaks = "1 year", date_labels = "%Y")+
                theme_ning(size.axis = 12,size.title = 14)

          }
        }
      })
  })
  # Tab: simulation----
  ## Print: Print the simulation infor ----
  output$printsimulatinginfo<-renderPrint({

    if(!is.null(inforprint$simulating)) writeLines(inforprint$simulating)
  })

  ## Action: Run model ----
  observeEvent(input$runSimulation, {
    inforprint$simulating<-"Processing log:"

    Sim_dates[["Start"]]<<-as.Date(format(input$dateSimulation[1],"%Y-%m-01"))
    Sim_dates[["End"]]<<-as.Date(format(input$dateSimulation[2],"%Y-%m-01"))
    # Only use year of the input
    # Sim_dates[["Start"]]<<-as.Date(paste0(format(input$dateSimulation[1],"%Y"),"-01-01"))
    # Sim_dates[["End"]]<<-as.Date(paste0(format(input$dateSimulation[2],"%Y"),"-12-01"))
    Sim_dates[["Seq_date"]]<<-seq.Date((Sim_dates[["Start"]]-years(warmup)),Sim_dates[["End"]],by = "month")

    climate_sel<-data_input[["Climate"]]%>%
      filter(BasinID==1)%>%
      mutate(Date=as.Date(paste0(Year,"-",Month,"-","01")))
    daterange_climate<-range(climate_sel$Date)
    # Check the input dataset
    if(daterange_climate[1]>Sim_dates[["Start"]] | daterange_climate[2]< Sim_dates[["End"]]){
      f_addinfo("simulating",paste0("!!! ERROR: You can select simulation data
      for the period of ",paste(range(climate_sel$Date),collapse="-")))
      return()
    }

    # Climate time range
    Sim_dates[["Start_climate"]] <<- daterange_climate[1]
    Sim_dates[["End_climate"]] <<- daterange_climate[2]

    # LAI time range
    lai_sel<-data_input[["LAI"]]%>%
      filter(BasinID==1)%>%
      mutate(Date=as.Date(paste0(Year,"-",Month,"-","01")))
    daterange_lai<-range(lai_sel$Date)
    Sim_dates[["Start_lai"]]<<-daterange_lai[1]
    Sim_dates[["End_lai"]]<<-daterange_lai[2]
    Lc_Nos<-length(lai_sel[1,])-4
    print(paste0("There are ",Lc_Nos," land cover types."))

    data_simulation<<-data_input
    # Select HUCs for simulation
    StationID<-as.integer(input$StationID)
    routpar<-NULL
    hru_flowlen<-NULL
    if(StationID==0) {
      BasinID_sel <- data_simulation[["Cellinfo"]]$BasinID
      f_addinfo("simulating","No Station is provided, so that all HUCs are used in the simulation!")
    }else if(is.null(input$Input_routpar)){
      f_addinfo("simulating","!!! Warning: Flow routing file is not provided!")

      BasinID_sel<-c(StationID)
      ## Subset each file
      data_simulation<<-lapply(data_simulation, function(x) subset(x,BasinID %in% BasinID_sel))

      f_addinfo("simulating",paste0("As no flow routing file is provided,
                                    only the selected BasinID = ",BasinID_sel, " is used for the simulation!"))

      #return()
    }else{
       routpar<-f_stream_level_pete(input$Input_routpar$datapath)
       BasinID_sel<-f_upstreamHUCs(BasinID =StationID ,routpar=routpar)
       BasinID_sel<-c(BasinID_sel,StationID)
       ## Subset each file
       data_simulation<<-lapply(data_simulation, function(x) subset(x,BasinID %in% BasinID_sel))

       f_addinfo("simulating",paste0("There are ",length(BasinID_sel),
                                    " HRUs in the upstreams of the station are selected for the simulation!"))
}

    # Update flow_len
    if("Flwlen_m" %in% names(data_simulation[["Cellinfo"]])) hru_flowlen<-data_simulation[["Cellinfo"]]$Flwlen_m

    # get the number of selected HUCs
    data_simulation[["BasinID"]]<<-BasinID_sel
    NoBasins<-length(BasinID_sel)

    # Select climate data for the HRU
     print("Processing climate data")
    ## Gapfill climate data for warming up
    if (min(Sim_dates[["Seq_date"]])<Sim_dates[["Start_climate"]]){

      climate_monthly<-data_simulation[["Climate"]]%>%
        group_by(BasinID,Month)%>%
        summarise_all(.funs = mean,na.rm=T)%>%
        dplyr::select(-Year)

      datesGap<-seq(min(Sim_dates[["Seq_date"]]),(Sim_dates[["Start_climate"]]-months(1)),by="month")

      climate_filled<-data.frame("BasinID"=rep(BasinID_sel,each=length(datesGap)),
                                 "Date"=datesGap)%>%
        mutate(Year=as.integer(format(Date,"%Y")),
               Month=as.integer(format(Date,"%m")))%>%
        merge(climate_monthly,by=c("BasinID","Month"))%>%
        dplyr::select(BasinID,Year,Month,Ppt_mm,Tavg_C)%>%
        rbind(data_simulation[["Climate"]])%>%
        arrange(BasinID,Year,Month)

      data_simulation[["Climate"]]<<-climate_filled
      Sim_dates[["Start_climate"]]<<-(Sim_dates[["Start"]]-years(warmup))

    }

    # Extract the sim period from the climate data for each HRU
    climate_date <- seq.Date(Sim_dates[["Start_climate"]], Sim_dates[["End_climate"]], by = "month")
    climate_ind <-which(climate_date %in% Sim_dates[["Seq_date"]])
    Sim_dates[["climate_index"]]<<-climate_ind
    print("Finished processing climate data")

    ## Gapfill LAI data
    print("Processing LAI")
    if (min(Sim_dates[["Seq_date"]])<Sim_dates[["Start_lai"]] | max(Sim_dates[["Seq_date"]])>Sim_dates[["End_lai"]]){

      lai_monthly<-data_simulation[["LAI"]]%>%
        group_by(BasinID,Month)%>%
        summarise_all(.funs = mean,na.rm=T)%>%
        dplyr::select(-Year)
      datesGap<-NULL
      if(min(Sim_dates[["Seq_date"]])<Sim_dates[["Start_lai"]]) datesGap<-seq(min(Sim_dates[["Seq_date"]]),(Sim_dates[["Start_lai"]]-months(1)),by="month")
      if(max(Sim_dates[["Seq_date"]])>Sim_dates[["End_lai"]]) datesGap<-c(datesGap,seq((Sim_dates[["End_lai"]]+months(1)),max(Sim_dates[["Seq_date"]]),by="month"))

      LAI_filled<-data.frame("BasinID"=rep(BasinID_sel,each=length(datesGap)),
                             "Date"=datesGap)%>%
        mutate(Year=as.integer(format(Date,"%Y")),
               Month=as.integer(format(Date,"%m")))%>%
        merge(lai_monthly,by=c("BasinID","Month"))%>%
        dplyr::select(BasinID,Year,Month,starts_with("Lc_"))%>%
        rbind(data_simulation[["LAI"]])%>%
        arrange(BasinID,Year,Month)

      data_simulation[["LAI"]]<<-LAI_filled
      Sim_dates[["Start_lai"]]<<-min(Sim_dates[["Seq_date"]],Sim_dates[["Start_lai"]])
      Sim_dates[["End_lai"]]<<-max(Sim_dates[["Seq_date"]], Sim_dates[["End_lai"]])

    }

    # Extract the sim period from the climate data for each HRU
    lai_date <- seq.Date(Sim_dates[["Start_lai"]], Sim_dates[["End_lai"]], by = "month")
    lai_ind <-which(lai_date %in% Sim_dates[["Seq_date"]])
    Sim_dates[["lai_index"]]<<-lai_ind
    print("Finished processing LAI")

    # date vectors
    jdate_index <- as.numeric(format(Sim_dates[["Seq_date"]], "%m"))
    jdate <-  c(15,46,76,107,137,168,198,229,259,290,321,351)
    Sim_dates[["jdate"]]<<-jdate[jdate_index]
    Sim_dates[["dateDays"]]<<-sapply(Sim_dates[["Seq_date"]],numberOfDays)

     # Load coefficients for SUN's ET and carbon calculation
    data_simulation[["ET_coefs"]]<<-NULL
    data_simulation[["WUE_coefs"]]<<-NULL
    print("Processing ET parameters")
    if(! is.null(input$Input_ETmodel)) {
      ET_coefs<-read.csv(input$Input_ETmodel$datapath,nrows = Lc_Nos)
      names(ET_coefs)<-c("LC_ID","Intercept","P_coef","PET_coef","LAI_coef","P_PET_coef","P_LAI_coef","PET_LAI_coef","IGBP","LC_Name")
      data_simulation[["ET_coefs"]]<<-ET_coefs
      f_addinfo("simulating","User defined ET model is read.")
    }
    if(! is.null(input$Input_WUEmodel)) {
      WUE_coefs<-read.csv(input$Input_ETmodel$datapath,nrows = Lc_Nos)
      names(WUE_coefs)<-c("LC_ID","WUE","RECO_Intercept", "RECO_Slope","IGBP","LC_Name")
      data_simulation[["WUE_coefs"]]<<-WUE_coefs
      f_addinfo("simulating","User defined WUE model is read.")
    }

    ## PET parameters
    par_petHamon<- rep(1,NoBasins)

    # Run the model
    print("runing the simulation")
    f_addinfo("simulating","Runing simulation ...")

    # get the real simulate period result
    if(mcores>1){
      result<-mclapply(BasinID_sel, WaSSI,datain=data_simulation,sim.dates=Sim_dates,mc.cores = mcores)
    }else{
      result<-lapply(BasinID_sel, WaSSI,datain=data_simulation,sim.dates=Sim_dates)
    }
    print("Finished runing the simulation")

    # process the result for HUC level
    lc_out<-lapply(result, "[[","lc_output")
    vars<-names(lc_out[[1]])
    lc_out<-lapply(vars, function(x) lapply(lc_out, "[[",x))
    names(lc_out)<-vars
    resultOutput[["lc_output"]]<<-lc_out
    print("Finished processing the Lc level output")

    # process the result for HUC level
    out<-lapply(result, "[[","output")
    Basin_out<-lapply(names(out[[1]]), function(x) sapply(out, "[[",x))
    names(Basin_out)<-names(out[[1]])
    resultOutput[["Basin_output"]]<<-Basin_out
    print("Finished processing the HUC level output")

    # routing based on catchment area
    hru_area<-data_simulation[["Cellinfo"]]$Area_m2
    print("Processing flow routing to the outlet HUC")
    if(input$Input_routingmethod=="NULL"){
      Station_out<-sapply(Basin_out, function(x) apply(x, 1, weighted.mean,hru_area) )
      Station_out<-round(Station_out,3)%>%
        as.data.frame()%>%
        mutate(Date=seq(Sim_dates$Start,Sim_dates$End,by="month"))
      resultOutput[["Station_output"]]<<-Station_out

      }else if(is.null(routpar) | is.null(hru_flowlen)){

      }else{
      # Channel flow routing from Lohmann routing model
      out2 <- mclapply(BasinID_sel, huc_routing,mc.cores = mcores)
      out3<-lapply(names(out2[[1]]), function(var) apply(sapply(out2, function(x) x[[var]]),1,sum))

    }
    f_addinfo("simulating","Finished simulation! ^_^")
    print("Finished!")

})

  ## Print: Print the simulation infor ----
  output$printsimplotinfo<-renderPrint({

    if(!is.null(inforprint$simulplotting)) writeLines(inforprint$simulplotting)
  })
  ## PlotAction: Plot model result----
  observeEvent(input$plotsimOut,{
    inforprint$simulplotting<-"Plotting log: "
    plot_dates<-c(input$plotSimDaterange[1],input$plotSimDaterange[2])

    output_monthly<- resultOutput[["Station_output"]]%>%
      mutate(Year=as.integer(format(Date,"%Y")),
             Month=as.integer(format(Date,"%m")))
    daterange_output<-range(output_monthly$Date)
    # Check the input dataset
    if(daterange_output[1]>plot_dates[1] | daterange_output[2]<plot_dates[2]){
      f_addinfo("simulplotting",paste0("!!! ERROR: You can only select the period of ",paste(range(daterange_output),collapse="-")))
      return()
    }

  df<-resultOutput[["Station_output"]]%>%
    filter(Date>=plot_dates[1] & Date<=plot_dates[2])

    output$SimOutplot <- renderPlot({
      input$plotsimOut
      #print(head(df))
      if(input$plotannualoutput){
        df%>%
          mutate(Year=as.integer(format(Date,"%Y")))%>%
          dplyr::select(-Date)%>%
          group_by(Year)%>%
          summarise_all(.funs = sum)%>%
          gather(Variable,Value,prcp:PET)%>%
          filter(Variable %in% c("flwTot","flwBase","flwSurface","aetTot","prcp"))%>%
          mutate(Variable=factor(Variable,levels = c("flwTot","flwBase","flwSurface","rain","prcp","aetTot"),
                                 labels = c("Total flow","Base flow","Surface flow","Effective rainfall","Precipitation","Actual evpotransipration")))%>%
          ggplot(aes(x=Year,y=Value,col=Variable))+
          labs(y="(mm/yr)")+
          scale_x_continuous(breaks = c(1980:2018))+
          geom_line()+geom_point()+
          theme_ning(size.axis = 12,size.title =14 )


      }else{
        df%>%
          gather(Variable,Value,prcp:PET)%>%
          filter(Variable %in% c("flwTot","flwBase","flwSurface","rain","prcp","aetTot"))%>%
          mutate(Variable=factor(Variable,levels = c("flwTot","flwBase","flwSurface","rain","prcp","aetTot"),
                                 labels = c("Total flow","Base flow","Surface flow","Effective rainfall","Precipitation","Actual evpotransipration")))%>%
          ggplot(aes(x=Date,y=Value,col=Variable))+
          labs(y="(mm/month)")+
          scale_x_date(date_breaks = "1 year", date_labels = "%Y")+
          geom_line()+geom_point()+
          theme_ning(size.axis = 12,size.title =14 )

      }

    })
})

  ## SaveAction: Save all output----
  observeEvent(input$savesimOut,{
    inforprint$simulplotting<-"Plotting log: "
    print("Processing the out to output format!")

    save(resultOutput,file="www/resultOutput.Rdata")

    f_addinfo("simulplotting","Processing the out to output format!")
    resultOutput[["Station_output"]]%>%
      dplyr::select(Date,prcp,rain,temp,LAI,PET,PET_hamon, aetTot,
                    flwTot, flwSurface, flwBase,everything())%>%
      write.csv(paste0("www/output/Output_outlet_monthly_TS.csv"),row.names=F)

    annual_output<-resultOutput[["Station_output"]]%>%
      mutate(Year=as.integer(format(Date,"%Y")))%>%
      dplyr::select(-Date)%>%
      group_by(Year)%>%
      summarise_all(.funs = sum,na.rm=T)%>%
      mutate(temp=temp/12,LAI=LAI/12,uztwc=uztwc/12,	uzfwc=uzfwc/12,
             lztwc=lztwc,	lzfpc=lzfpc/12,	lzfsc=lzfsc/12)%>%
      dplyr::select(Year,prcp,rain,temp,LAI,PET,PET_hamon, aetTot,
                    flwTot, flwSurface, flwBase,everything())
    resultOutput[["Output_station_annual"]]<<-annual_output
    annual_output%>%
      write.csv(paste0("www/output/Output_outlet_annual_TS.csv"),row.names=F)

    # Process output for each HUC
   f_ReshapebyBasinID<-function(da,BasinIDs,seq_dates){
     out<- da%>%
       do.call(what = rbind)%>%
       as.data.frame()
     names(out)<-BasinIDs
     out[["Date"]]=seq_dates
     out[["Variable"]]=rep(names(da),each=length(seq_dates))
     return(out)
   }
    BasinIDs<-data_simulation$BasinID
    seq_dates<-seq(Sim_dates$Start,Sim_dates$End,by="month")

    .a<-f_ReshapebyBasinID(da=resultOutput[["Basin_output"]],
                          BasinIDs=BasinIDs,
                          seq_dates=seq_dates)
    Output_BasinID<-.a%>%
        melt(id=c("Date","Variable"))%>%
        rename(BasinID=variable)%>%
        mutate(BasinID=as.integer(as.character(BasinID)))%>%
        dcast(BasinID+Date~Variable)%>%
        dplyr::select(BasinID,Date,prcp,rain,temp,LAI,PET,PET_hamon, aetTot,
                   flwTot, flwSurface, flwBase,everything())
    resultOutput[["Output_Basin_monthly"]]<<- Output_BasinID
    Output_BasinID%>%
        write.csv(paste0("www/output/Output_BasinID_monthly_TS.csv"),row.names=F)

    Output_BasinID_ann<- Output_BasinID%>%
      mutate(Year=as.integer(format(Date,"%Y")))%>%
      dplyr::select(-Date)%>%
      group_by(BasinID,Year)%>%
      summarise_all(.funs = sum,na.rm=T)%>%
      mutate(temp=temp/12,LAI=LAI/12,uztwc=uztwc/12,	uzfwc=uzfwc/12,
             lztwc=lztwc,	lzfpc=lzfpc/12,	lzfsc=lzfsc/12)%>%
      dplyr::select(BasinID,Year,prcp,rain,temp,LAI,PET,PET_hamon, aetTot,
                    flwTot, flwSurface, flwBase,everything())
    resultOutput[["Output_Basin_annual"]]<<- Output_BasinID_ann
    Output_BasinID_ann%>%
      write.csv(paste0("www/output/Output_BasinID_annual_TS.csv"),row.names=F)

    Output_BasinID_avg<-Output_BasinID_ann%>%
      dplyr::select(-Year)%>%
      group_by(BasinID)%>%
      summarise_all(.funs = mean,na.rm=T)
    resultOutput[["Output_Basin_avg"]]<<- Output_BasinID_avg
    Output_BasinID_avg%>%
      write.csv(paste0("www/output/Output_BasinID_avg.csv"),row.names=F)

    # Process output for each land cover
    f_ReshapebyLc<-function(da,Lcs,Output_BasinID){
      out<- da%>%
        do.call(what = rbind)%>%
        as.data.frame()
      names(out)<-Lcs
      out<-cbind(Output_BasinID[c("BasinID","Date")],out)
      return(out)
    }
    Lcs<-paste0("Lc_",c(1:9))

    for(var in names(resultOutput[["lc_output"]])){
      lc_mon_output<-f_ReshapebyLc(da=resultOutput[["lc_output"]][[var]],
                  Lcs=Lcs,
                  Output_BasinID=Output_BasinID)
      resultOutput[["Output_lc_monthly"]]<<- lc_mon_output
      lc_mon_output%>%
        write.csv(paste0("www/Output/Output_Lc_monthly_TS_",var,".csv"),row.names=F)

      dd_mon<-lc_mon_output%>%
        mutate(Year=as.integer(format(Date,"%Y")))%>%
        dplyr::select(-Date)
      if(var %in% c("aetTot","flwTot","flwSurface","flwBase","PET")){
        dd_ann<-dd_mon %>%
          group_by(BasinID,Year) %>%
          summarise_all(.funs = sum)
        dd_ann%>%
          write.csv(paste0("www/Output/Output_Lc_annual_TS_",var,".csv"),row.names=F)

        dd_ann %>%
          group_by(BasinID) %>%
          dplyr::select(-Year)%>%
          summarise_all(.funs = mean)%>%
          write.csv(paste0("www/Output/Output_Lc_avg_",var,".csv"),row.names=F)
      }else{
        dd_mon %>%
          group_by(BasinID,Year) %>%
          summarise_all(.funs = mean)%>%
          write.csv(paste0("www/Output/Output_Lc_annual_TS_",var,".csv"),row.names=F)

        dd_mon %>%
          group_by(BasinID) %>%
          dplyr::select(-Year)%>%
          summarise_all(.funs = mean)%>%
          write.csv(paste0("www/Output/Output_Lc_avg_",var,".csv"),row.names=F)
      }


    }
    print("Finished processing the output!")
    f_addinfo("simulplotting",
              paste0("The simulated result has been save to default foler 'www/output'"))

  })

  ## ExportAction: save the processed input data  ----
  output$downloadoutputData <- downloadHandler(

    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      f_addinfo("simulating","Files are downloaded! ^_^")
      paste("Outputs", "zip", sep=".")
    },

    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(fname) {

      fs=dir("www/output/","csv",full.names = T)

      withProgress(message = 'Compressing files ...', value = 0.5,{
        zip::zipr(zipfile=fname,files=fs)
      })

    },
    contentType="application/zip"
  )

    ## Action: select time period for calibration and validation----

  } # END
)


# WaSSI for each HRU and calculate adjust temp and pet values
WaSSI<-function(hru,datain,sim.dates){

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
