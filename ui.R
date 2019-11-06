#
# This is the user-interface definition of a dWaSSI-C web application. You can
# run the application by clicking 'Run App' above.
#

library(shiny)

# Define UI for dWaSSI-C application
shinyUI(
  navbarPage(
    "dWaSSIC model",
    # Tab: About ----
    tabPanel("About",
             column(10, wellPanel(
               #tags$img(src = 'coweeta_logo.jpg',width=600),

              tags$h4("Description"),
              tags$html("This is the distributed R-based ", tags$b("WaSSI")," model. The original model code is written in Fortran and is publicly accessible. The R version is translated from the Fortran code developed by Peter Caldwell and R code of Sac-sma developed by Umit Taner."),
              tags$br(),
              tags$html("The web application for WaSSI allows users to define a custom simulation scenario,
                        view/download model inputs and outputs in tabular and graphical form for a location of interest,
                        and view/export model outputs spatially for a variety of time scales using an interactive map viewer.
                        Users may select their location in the map viewer,
                        select a specific HUC, or input a zip code to view model inputs and outputs"),
              tags$hr(),
              tags$h4("WaSSI"),
              tags$html("WaSSI is an integrated, process-based model that
                        can be used to project the effects of forest land cover change,
                        climate change, and water withdrawals on river flows, water supply stress,
                        and ecosystem productivity (i.e. carbon dynamics).
                        WaSSI operates on a monthly time step at the HUC-4 (8-digit HUC)
                        watershed scale (see more on HUCs) and across Mexico at the 0.5 degree scale.
                        For the conterminous U.S., the model can also be run at the HUC12
                        scale for water and carbon balances from 1960 to 2012.
                        As water yield and carbon sequestration are tightly coupled,
                        WaSSI can be used to evaluate trade-offs among management strategies for these ecosystem services."),
              tags$br(),
              tags$html("Further information is available at: USDA - Forest Service or the Guide"),
              tags$hr(),
               h4("Author:"),
               p("Dr Ning Liu.",p(),a("Email: LN1267@Gmail.com",href="mailto:LN1267@Gmail.com")),
               h4("Terms:"),
               p("This app is under the ",a("MIT",href="https://opensource.org/licenses/MIT"),"License.")

             ))
    ),
    # Tab: Upload data ----
    tabPanel("Upload input data",
             sidebarLayout(
               sidebarPanel(
                 # Application title
                 titlePanel("Select inputs data for dWaSSIC model"),
                 ## HTML: introduction of the input data
                 HTML("<p>The input data of dWaSSIC model includes:</p>
                        <ol>
                        <li>Watershed characters for each Hydrologic research unit (HRU);</li>
                        <li>Monthly climate time series for each HRU;</li>
                        <li>Monthly LAI time series for each land cover in each HRU;</li>
                        <li>Soil parameters for each HRU.</li>
                        </ol>"),
                 tags$hr(),
                 # Action: Read input data ----
                 actionButton("readdata","Read all input data!"),
                 HTML("<p>You can prepare those files as a 'csv/txt' format.</p>
                      "),
                 ## Choose each file----
                 fileInput("Input_cellinfo", "Cellinfo [BasinID, Latitude, Longitude, Elevation_m, Area_m2 and ratio for each lc]"),
                 fileInput("Input_climate", "Monthly Climate [BasinID, Year, Month, Precip_mm, Tmean_C]"),
                 fileInput("Input_LAI", "Monthly LAI [BasinID, Year, Month, and LAI for each lc]"),
                 fileInput("Input_soilinfo", "Soil info [BasinID, uztwm,uzfwm,uzk,zperc,rexp,lztwm,lzfsm,lzfpm,lzsk,lzpk,pfree]")

               ),

               # MainPanel for upload data----
               mainPanel(
                 ## Output: Summary of the selected plotting file----
                 h2("The summary of the input file"),
                 tags$h6("Process log: "),
                 verbatimTextOutput("printreadinginfo"),
                 verbatimTextOutput("printsummary")
               )
             )
    ),
    # Tab: Prepare input data ----
    tabPanel("Input preparation",
             sidebarLayout(
               sidebarPanel(
                 # Application title

                 h2("Upload watershed Shapefile and other raster data"),
                 HTML("<p>Please chose a '*gml' file to upload.</p>
                      "),
                 fileInput('Input_basin', 'Choose watershed shapefile', multiple=FALSE, accept="gml"),


                 # Action: Plot basin ----
                 actionButton("processrasters","process input!"),
                 # Action: Process input ----
                 actionButton("plotrasterdata","Plot input data!"),
                 actionButton("saveInputData","Save input data!"),
                 downloadButton("downloadInputData", label = "Export input data"),
                 #h2("Upload all original raster files"),
                  HTML("<p>Please chose a '*.tif' or '*.nc' monthly stacked precipitation and temperation file to upload.</p>
                      "),

                 fileInput('Input_dem_raster', 'Choose elevation raster file (Optional)', multiple=FALSE, accept="tif"),
                 checkboxInput("updateLc", "Update land cover", FALSE),
                 fileInput('Input_lc_raster', 'Choose land cover raster file', multiple=FALSE, accept="tif"),
                 textInput("yearStartClimate", "Start year of climate data","2000",width = "50%"),
                 checkboxInput("updateClimate", "Update", FALSE),
                 fileInput('Input_temp_raster', 'Choose monthly temperature raster file', multiple=FALSE, accept="tif"),
                 fileInput('Input_precp_raster', 'Choose monthly precipitation raster file', multiple=FALSE, accept="tif"),
                 textInput("yearStartLai", "Start year of Lai data","2000",width = "50%"),
                 checkboxInput("updateLai", "Update", FALSE),
                 # shinyFilesButton("Input_lai_fpath", "Choose a file" ,
                 #                  title = "Please select a file:", multiple = FALSE,
                 #                  buttonType = "default", class = NULL),
                 textInput("Input_lai_fpath", "Type in the lai data path","~"),
                 fileInput('Input_lai_raster', 'Choose monthly LAI raster file', multiple=FALSE, accept="tif"),
                 checkboxInput("updateSoil", "Update", FALSE),
                 fileInput('Input_soil_raster', 'Choose Soil raster file', multiple=FALSE, accept="tif"),
                 fileInput('Input_imp_raster', 'Choose impverious raster file [Optional]', multiple=FALSE, accept="tif")

               ),

               # MainPanel for upload data----
               mainPanel(
                 ## Output: Head of the selected plotting file----
                 # h2("The Map of the Watershed"),
                 # leafletOutput("basinmap"),


                 h2("Plot an interactice map and print summary of processed input."),
                 HTML("<p>It will take a while for the data process. Please wait ...</p>
                      "),
                 leafletOutput("basinrastermap"),
                 tags$h4("Processing log:"),
                 tags$h5("We recommand upload all raster are in the same projection and cellsize!"),
                 tags$h5("It will take more time to process,
                         if your land cover data's cellsize is much smaller than your LAI data."),
                 verbatimTextOutput("printprocessinginfo"),
                 verbatimTextOutput("txt_file"),
                 verbatimTextOutput("prntraster")


               )
             )
    ),
    # Tab: Plot data ----
    tabPanel("Plot input data",
             sidebarLayout(
               sidebarPanel(
                 # Application title

                 ## Plot the information----
                 h2("Select a dataset to plot"),
                 ## Input: Select a file to plot----
                 selectInput("daname2plot", "Name of the dataset.",
                             choices=c("Climate","LAI","Cellinfo","Soilinfo")),
                 ## Input: Select a variable ----
                 textInput("plotvar", "Type the variable name for ploting","Ppt_mm"),
                 textInput("varnames", "Type the printing name for those variables","Ppt_mm"),
                 textInput("plotBasinID", "Type the BasinID name for ploting","1"),
                 #textInput("Basinnames", "Type the printing name for those Basins","ID1"),
                 # Input: Year for plotting ---
                 sliderInput("plotyrrange", "Select the year range", 1970, 2017, value = c(2000, 2010)),
                 # Input: Month ----
                 checkboxGroupInput("plotmonths","Select the months to plot",c(1:12),selected=c(1:12),inline=T),

                 # Action: plot input data ----
                 actionButton("plotdata","Plot the selected input data!")
               ),

               # MainPanel for ploting data----
               mainPanel(

                 h2("Plot the selected input dataset!"),
                 tags$html("This is the processing log:"),
                 verbatimTextOutput("printplottinginfo"),
                 verbatimTextOutput("printbasininfo"),

                 plotOutput("Plotinput",height = 800)
               )
             )
    ),
    # Tab: Simulation-----
    tabPanel("Run the simulation",
             sidebarLayout(
               sidebarPanel(
                 # Application title
                 titlePanel("Select the outlet station's BasinID and timeperiod for simulation"),
                 ## HTML: Introduction ----

                 ## Input: Select a Station ----
                 textInput("StationID", "Type the BasinID of the hydrologic station!","0"),
                 fileInput('Input_routpar', 'Choose flow routing parameter file (Optional)', multiple=FALSE, accept="csv"),
                 selectInput("Input_routingmethod", "Select method for flow routing.",
                             choices=c("NULL","AreaWeighted","Method1","Method2","Method3")),
                 fileInput('Input_ETmodel', 'Choose your ET model file (Optional)', multiple=FALSE, accept="csv"),
                 fileInput('Input_WUEmodel', 'Choose your WUE model file (Optional)', multiple=FALSE, accept="csv"),

                 # Action: Subset data ----
                 #actionButton("subSimData","Subset the selected dataset"),

                 # Input: Simulation period
                 dateRangeInput("dateSimulation","Select the date range for simulation"),
                 #sliderInput("dateSimulation", "Year simulation", 1970, 2017, value = c(2000, 2010)),

                 # Input: Checkbox if calibate model ----
                 checkboxInput("Calibration", "Calibration", FALSE),
                 ## Action: run model
                 actionButton("runSimulation","Run Simulation!"),
                 actionButton("savesimOut","Save output!"),
                 downloadButton("downloadoutputData", label = "Export output data!"),
                 tags$hr(),
                 tags$h5("Simulation log:"),
                 verbatimTextOutput("printsimulatinginfo")
                  ),

                 # Mainpanel: of Simulation----
                 mainPanel(


                   h2("Plot water and carbon processes"),
                  # Input: Result plot variables ----
                  checkboxGroupInput("plotvars","Select the variables to plot",c("P","T","Q"),inline = T),

                  # Input: Year for plotting ---
                  #sliderInput("plotSimDaterange", "Year plot", 1970, 2017, value = c(2000, 2010)),
                  dateRangeInput("plotSimDaterange","Select the date range for plotting"),
                  textInput("plotSimuBasinID", "Type the BasinID for plotting!","0"),
                  # Action: Plot simulated result ---
                   actionButton("plotsimOut","Plot"),
                  verbatimTextOutput("printsimplotinfo"),
                    # Output: ploted simulated result
                   plotOutput("SimOutplot")
                 )
              )
        ),
    # Tab: Help-----
    tabPanel("Help",
             column(5, wellPanel(

               tags$h4("Instruction for using this app:"),
               tags$ol(
                 tags$li("Read the target format input data;"),
                 tags$li("Process data to the target format from the scratch;"),
                 tags$li("Plot the input data, if you want;"),
                 tags$li("Run the simulation and save the result;")
               ),
               tags$hr(),
               tags$h4("Framework of dWaSSI model:"),
               tags$img(src = 'Framework_WaSSIC.png',width=400)

             ))
    )
  )
)
