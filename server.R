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
                            plotting="Plotting log: ")
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

  ## Print: Print the processig infor ----
  output$printprocessinginfo<-renderPrint({
    if(!is.null(inforprint$processing)) writeLines(inforprint$processing)
  })

  ## Action: Read Basin shapefile ----
  if (!exists("BasinShp")) BasinShp<-NULL
  observeEvent(input$Input_basin,{
    f_addinfo("processing","Reading shapfile ...")
    Basins<- readOGR(input$Input_basin$datapath)
    Basins$BasinID<-c(1:length(Basins[,1]))

    # Add latitude and longitude infor to the Basin
    if(!"Latitude"  %in% names(Basins)) {

      # get the contral coordinates of each polygon
     if(is.na(Basins)) proj4string(Basins)<-CRS("+init=epsg:4326")
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

    if (is.null(BasinShp)) {
      f_addinfo("processing","You need to upload a shapefile including the basins.")
      return()
    }

    # process Climate input
    if(is.null(input$Input_temp_raster)| is.null(input$Input_precp_raster)){
      f_addinfo("processing","Warning: You need to provide the climate raster time series!")
    }else{
      f_addinfo("processing","Processing climate data ...")
      Tmean_catchment<-f_sta_shp_nc(input$Input_temp_raster$datapath,BasinShp,varname = "Tavg_C",start = 2000,zonal_field = "BasinID")
      Pre_catchment<-f_sta_shp_nc(input$Input_precp_raster$datapath,BasinShp,varname = "Ppt_mm",start = 2000,zonal_field = "BasinID")
      climate<-Pre_catchment%>%
        mutate(Tavg_C=Tmean_catchment$Tavg_C)%>%
        arrange(BasinID,Year,Month)%>%
        mutate(Ppt_mm=round(Ppt_mm,2))%>%
        mutate(Tavg_C=round(Tavg_C,2))%>%
        dplyr::select(BasinID,Year,Month,Ppt_mm,Tavg_C)
      data_input[["Climate"]]<<-climate
      f_addinfo("processing","Finished processing climate input!")
      print(paste0("Finished processsing Climate data ..."))
      #print(summary(climate))
    }

    # Process cellinfo
    if(is.null(input$Input_lc_raster)){
      f_addinfo("processing","Warning: There is no land cover input!")
    }else{
      f_addinfo("processing","processing Cellinfo data ...")
      print(paste0("processing Cellinfo ..."))
      cellinfo<-f_cellinfo(classfname=input$Input_lc_raster$datapath,
                           Basins=BasinShp,
                           byfield="BasinID",
                           demfname=input$Input_dem_raster$datapath)
      data_input[["Cellinfo"]]<<-cellinfo
      f_addinfo("processing","Finished processing cellinfo!")
    }
      # process LAI input
    if(is.null(input$Input_lc_raster) | is.null(input$Input_lai_raster) ){
      f_addinfo("processing","Warning: There is no land cover or LAI data!")
    }else{
      laiinfo<-f_landlai(lcfname=input$Input_lc_raster$datapath,
                          laifname=input$Input_lai_raster$datapath,
                          Basins=BasinShp,
                          byfield="BasinID",
                          yr.start=1982,
                          yr.end=2014)
      data_input[["LAI"]]<<-laiinfo
      f_addinfo("processing","Finished processing Lai data!")
    }
      # process SOIL input
    if(is.null(input$Input_soil_raster) | is.null(input$Input_lai_raster) ){
      f_addinfo("processing","Warning: There is no soil data!")
    }else{
      soilinfo<-f_soilinfo(soilfname=input$Input_soil_raster$datapath,
                          Basins=BasinShp)
      data_input[["Soilinfo"]]<<-soilinfo
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

        if(!is.null(input$Input_temp_raster$datapath)){
          f_addinfo("processing","Calculate average for temperature time series!")
          Tmean_br<-brick(input$Input_temp_raster$datapath)
          beginCluster()
          Tmean_avg <<- clusterR(Tmean_br, calc, args=list(fun=mean,na.rm=T))
          endCluster()
          }
        if(!is.null(input$Input_precp_raster$datapath)){
          Pre_br<-brick(input$Input_precp_raster$datapath)
          beginCluster()
          Pre_br <- clusterR(Pre_br, calc, args=list(fun=mean,na.rm=T))
          Pre_avg<<-Pre_br*12
          endCluster()
        }

        if(!is.null(input$Input_lai_raster$datapath)){
          Lai_br<-brick(input$Input_lai_raster$datapath)
          beginCluster()
          Lai_avg <<- clusterR(Lai_br, calc, args=list(fun=mean,na.rm=T))
          endCluster()
        }

      output$basinrastermap <- renderLeaflet({
        # Plot the BasinShp
        input_leaflet<-leaflet() %>%
            addTiles(group = "OSM (default)") %>%
            addProviderTiles('Esri.WorldImagery',group = "Esri.Imagery")%>%
            addFullscreenControl()
        grps<-c("OSM (default)", "Esri.Imagery")
        ovlgrps<-NULL

        if (!is.null(BasinShp)) {
          ovlgrps<-c(ovlgrps,"Watershed")
          input_leaflet<-input_leaflet%>%
            addPolygons(data=BasinShp,weight=1,col = 'red',fillOpacity = 0.2,
                        highlight = highlightOptions(color='white',weight=1,
                                                     bringToFront = TRUE),
                        group = "Watershed")%>%
            addMarkers(lng = BasinShp$Longitude, lat = BasinShp$Latitude,
                       popup = as.character(BasinShp$BasinID),
                       label = paste0("BasinID = ",BasinShp$BasinID))
        }

        if(!is.null(input$Input_temp_raster$datapath)){
          ovlgrps<-c(ovlgrps,"Temperature")
          pal_tmp <- colorNumeric(c("cyan", "yellow", "red"), values(Tmean_avg),
                            na.color = "transparent")
          print(Tmean_avg)
          input_leaflet<-input_leaflet%>%
            addRasterImage(Tmean_avg, colors = pal_tmp, opacity = 0.8, group = "Temperature") %>%
            addLegend("bottomleft",pal = pal_tmp, values = values(Tmean_avg),
                      title = "Mean temp (°C)", group = "Temperature")
        }

        if(!is.null(input$Input_precp_raster$datapath)){
          ovlgrps<-c(ovlgrps,"Precipitation")
        pal_prcp <- colorNumeric(c("azure", "cornflowerblue","darkblue"), values(Pre_avg),
                                na.color = "transparent")
        input_leaflet<-input_leaflet%>%
          addRasterImage(Pre_avg, colors = pal_prcp, opacity = 0.8, group = "Precipitation") %>%
          addLegend("bottomleft",pal = pal_prcp, values = values(Pre_avg),
                    title = "Annual precipation (mm/yr)", group = "Precipitation")
        }

        if(!is.null(input$Input_dem_raster$datapath)){
          ovlgrps<-c(ovlgrps,"Elevation")
          dem<-raster(input$Input_dem_raster$datapath)
          pal_dem <- colorNumeric(c("green","yellow","red"), values(dem),
                                  na.color = "transparent")

          input_leaflet<-input_leaflet%>%
            addRasterImage(dem, colors = pal_dem, opacity = 0.8, group = "Elevation") %>%
            addLegend("bottomleft",pal = pal_dem, values = values(dem),
                    title = "Elevation (m)", group = "Elevation")
        }

        if(!is.null(input$Input_lai_raster$datapath)){
          ovlgrps<-c(ovlgrps,"LAI")
          pal_lai <- colorNumeric(c("red","yellow","green"), values(Lai_avg),
                                  na.color = "transparent")

          input_leaflet<-input_leaflet%>%
            addRasterImage(Lai_avg, colors = pal_lai, opacity = 0.8, group = "LAI") %>%
            addLegend("bottomleft",pal = pal_lai, values = values(Lai_avg),
                      title = "Leaf area index (m2/m2)", group = "LAI")
        }

        if(!is.null(input$Input_lc_raster$datapath)){
          lc<-raster(input$Input_lc_raster$datapath)
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
  output$Plotinput <- renderPlot({

    if(input$plotdata>0){
      print(paste0("Print data ....",input$daname2plot))
      inforprint$info<-paste0("Print data ....",input$daname2plot)
      if(!input$daname2plot %in% names(data_input)) {
        inforprint$info<-"no data"
        return()}
      if(!input$plotvar %in% names(data_input[[input$daname2plot]])) {
        f_addinfo("plotting","The variable is not in the dataset")
        return()
      }
      df<-data_input[[input$daname2plot]]%>%
          filter(Year>=input$plotyrrange[1] & Year<=input$plotyrrange[2])%>%
          filter(Month %in% input$plotmonths)%>%
          mutate(Date=as.Date(paste0(Year,"-",Month,"-","01")))%>%
          dplyr::select(-Year,-Month)
      if(nrow(df)<1) {
        inforprint$info<-"no data"
        return()}

      if(input$daname2plot == "Climate"){

        names(df)<-c("HRU_ID","Precip","Tmean","Date")
        df$HRU_ID<-factor(df$HRU_ID)
        p1<-ggplot(df,aes(x=Date,y=Precip,col=HRU_ID))+geom_line()+geom_point()+
          ggtitle(paste0("Monthly Precipitation (mm/month)"))+
          facet_wrap(HRU_ID~.,ncol=3)+
          theme_ning(size.axis = 12,size.title = 14)
        p2<-ggplot(df,aes(x=Date,y=Tmean,col=HRU_ID))+geom_line()+geom_point()+
          ggtitle(paste0("Monthly mean temperature (C)"))+
          facet_wrap(HRU_ID~.,ncol=3)+
          theme_ning(size.axis = 12,size.title = 14)
        #print(p1)
        multiplot(p1,p2,cols = 1)
      }else if(input$daname2plot == "LAI"){
        names(df)<-c("HRU_ID","LAI","Date")
        df$HRU_ID<-factor(df$HRU_ID)
        ggplot(df,aes(x=Date,y=LAI,col=HRU_ID))+geom_line()+geom_point()+
          ggtitle(paste0("Monthly Leaf area index (m2/m2)"))+
          facet_wrap(HRU_ID~.,ncol=3)+
          theme_ning(size.axis = 12,size.title = 14)

      }

    }

  })

  # Tab: simulation----
  # Sidepanel: ----

  ## output$uploadfilename:Select a uploaded file for subsetting----
  output$uploadfilename <- renderUI({
    # If missing input, return to avoid error later in function
    if(is.null(input$files)){
      uploadedfiles<-"test.csv"
    }else{
      uploadedfiles <-dir("tmp/",'.csv')[-which(dir("tmp/",'.csv')=="test.csv")]
    }

    # Create the checkboxes and select them all by default
    selectInput("uploadfilename", "Select the upload filename for simulation",
                choices  = uploadedfiles,
                selected = "test.csv"
    )
  })


  ## Action: Subset data by Station ID and period----
  observeEvent(input$subdata, # subsetting input data
               output$subsetdata<-renderText({
                 print("Subseting data has started, please wait .....")
                 f_subset(input$StationID,input$subsetdaterange,input$uploadfilename)
                 })
  )

  ## Action: Print the input data info ----
  observeEvent(input$printinfo,
               output$inputfilesummary<-renderPrint({
                 if(!file.exists(paste0("www/",input$inputfiles,".csv"))) return(print(paste0("There is no ", input$inputfiles," input")))
                 df<-read.csv(paste0("www/",input$inputfiles,".csv"),sep=",",stringsAsFactors = F)
                 if("Timestamp" %in% names(df)) df$Timestamp<-as.Date(df$Timestamp,format="%Y-%m-%d")
                 print("Summary informaiton")
                 print(summary(df))
                 print("Header informaiton")
                 print(head(df))
                 print("Tailer informaiton")
                 print(tail(df))
               })

  )

  ## Mainpanel: ----
  ## Action: Merge climate and lai data ----
  observeEvent(input$mergeinput,
               output$mergeinputout<-renderPrint({
                 f_merge_input()
               })
  )
  ## Action: Calculate WUE based on vegetation ratio ----
  observeEvent(input$calwue,
               output$wues<-renderPrint({
                 if(!file.exists(paste0("www/LUCC.csv"))) return(print("There is no Landcover input"))
                 vegratio<-read.csv("www/LUCC.csv",sep=",",stringsAsFactors = F)
                 f_WaSSI_WUE(vegratio)
               })
  )

  ## Action: Calculate ET based on SUN ET function ----
  observeEvent(input$calsunET,
               output$ETs<-renderPrint({
                 if(!file.exists(paste0("www/LUCC.csv"))) return(print("There is no Landcover input"))
                 vegratio<-read.csv("www/LUCC.csv",sep=",",stringsAsFactors = F)
                 wues<-f_WaSSI_WUE(vegratio)
                 if(!file.exists(paste0("www/pars.csv"))) return(print("There is no pars input"))
                 print("Calculating monthly ET based on SUN's function")
                 pars<-read.csv("www/pars.csv",header = T)
                 pars<-c("LAT"=pars$LAT,"LONG"=pars$LONG,"ALT"=pars$ALT,wues)
                 print(pars)
                 load("www/data_TS.Rdata")
                 ET<-f_SUN_ET(data_TS,pars)
                 save(ET,file = "www/ET.RData")
               })
  )

  ## Action: Transfer monthly input to daily input ----
  observeEvent(input$mon2day,
               output$mon2dayout<-renderPrint({
                 load("www/ET.RData")
                 data_TS_daily<-month2daily(ET,fields = c("P","Q","PET","E"),ts = T)
                 save(data_TS_daily,file="www/data_TS_daily.RData")
               })
  )

  ## Action: Calculate snow ----
  observeEvent(input$calsnow,
               output$calsnowout<-renderPrint({
                 load("www/data_TS_daily.RData")
                 snow_out<-f_snow_calculation(data_daily)
                 save(snow_out,file="www/snow_out.RData")
               })
  )

  ## Action: Run model ----
  observeEvent(input$runmodel,
               output$ning2<-renderPrint({
                 print("Simulation has started, please wait .....")
                 load("www/data_TS_daily.RData")
                 soil<-as.vector(read.csv("www/SOIL.csv"))
                 soil_pars<-list("uztwm" = 1, "uzfwm" = 150, "uzk" = 0.1, "pctim" = 1e-06, "adimp" = 0,
                                 "zperc" = 28.8997, "rexp" = 5, "lztwm" = 205.652, "lzfsm" = 758.774,
                                 "lzfpm" = 1000, "lzsk" = 0.149213, "lzpk" = 0.00607691, "pfree" = 0.582714)
                 soil_pars[names(soil)]<-soil[names(soil)]

                 if(as.integer(format(range(index(data_TS_daily))[1],"%Y"))>input$simulatedaterange[1]) return(print(paste0("The start date should be later than ",range(index(data_TS_daily))[1])))
                 if(as.integer(format(range(index(data_TS_daily))[2],"%Y"))<input$simulatedaterange[2]) return(print(paste0("The end date should be earler than ",range(index(data_TS_daily))[2])))

                 outp<-f_WaSSI_Q(data_TS_daily,soil_pars,
                                 calibrate =F,
                                 y_s = input$simulatedaterange[1],
                                 y_e = input$simulatedaterange[2]
                                 )
                 save(outp,file="www/outp.RData")
                 })
  )

  ## Action: Plot the basic result----
  observeEvent(input$plotres,{

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
