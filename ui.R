#
# This is the user-interface definition of a dWaSSI-C web application. You can
# run the application by clicking 'Run App' above.
#

library(shiny)

# Define UI for dWaSSI-C application
shinyUI(
  navbarPage(
    "dWaSSIC model",
    # Tab: Upload data ----
    tabPanel("Upload input data",
             sidebarLayout(
               sidebarPanel(
                 # Application title
                 titlePanel("Select inputs data for dWaSSIC model"),
                 ## HTML: introduction of the input data
                 HTML("<p>The input data of dWaSSIC model includes:</p>
                        <ol>
                        <li>Watershed infomation</li>
                        <li>Monthly climate time series for each HRU</li>
                        <li>Monthly LAI data</li>
                        <li>Landcover</li>
                        <li>Soil parameters</li>
                        <li>Shapefiles of your study area (Optional)</li>
                        </ol>"),
                 tags$hr(),
                 # Action: Read input data ----
                 actionButton("readdata","Read all input data!"),
                 HTML("<p>You can prepare those files as a 'csv/txt' format.</p>
                      "),
                 ## Choose each file----
                 fileInput("Input_HRUinfo", "Watershed info [HRUID, Latitude, Longitude, Elevation_m and Area_m2]"),
                 fileInput("Input_climate", "Monthly Climate [HRUID, Year, Month, Precip_mm, Tmean_C]"),
                 fileInput("Input_LAI", "Monthly LAI [HRUID, Year, Month, LAI]"),
                 fileInput("Input_landinfo", "Land cover info [HRUID, IGBP,Ratio]"),
                 fileInput("Input_soil", "Soil info [HRUID, uztwm,uzfwm,uzk,zperc,rexp,lztwm,lzfsm,lzfpm,lzsk,lzpk,pfree]")

               ),

               # MainPanel for upload data----
               mainPanel(
                 ## Output: Summary of the selected plotting file----
                 h2("The summary of the input file"),
                 verbatimTextOutput("printsummary")
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
                             choices=c("Climate","LAI")),
                 # Input: Year for plotting ---
                 sliderInput("plotyrrange", "Select the year range", 1970, 2017, value = c(2000, 2010)),
                 # Input: Month ----
                 checkboxGroupInput("plotmonths","Select the months to plot",c(1:12),selected=c(6:11),inline=T),

                 # Action: plot input data ----
                 actionButton("plotdata","Plot the selected input data!")
               ),

               # MainPanel for upload data----
               mainPanel(

                 ## Output: Plot of the selected two columns----
                 h2("Please select two columns to plot the data"),
                 verbatimTextOutput("distPloterror"),
                 plotOutput("Plotinput")
               )
             )
    ),
    # Tab: Simulation-----
    tabPanel("Run the simulation",
             sidebarLayout(
               sidebarPanel(
                 # Application title
                 titlePanel("Select the station and timeperiod for simulation"),
                 ## HTML: Introduction ----
                  HTML("<p>Please select a station for simulation:</p>
                      <ol>
                      <li>Climate data</li>
                      <li>LAI data</li>
                      <li>Landcover</li>
                      <li>Soil</li>
                      </ol>"),
                 tags$hr(),
                 ## Input: Select a Station ----
                 textInput("StationID", "Type the Station ID for simulation","1"
                 ),
                 ## UIoutput: select a file to subset ----
                 uiOutput("uploadfilename"),

                 # Input: Time range for subsetting ----
                 dateRangeInput("subsetdaterange","Select the date range for simulation"),

                 # Action: Subset data ----
                 actionButton("subdata","Run Subset"),
                 verbatimTextOutput("subsetdata"),

                 h2("Select a file to print the basic information"),
                 # Input: select an input data to print  ----
                 selectInput("inputfiles", "Select a file to print the information",c("Climate","LAI","SOIL","LUCC","Flow"),"Climate"
                 ),


                 tags$hr()
                  ),

                 # Mainpanel: of Simulation----
                 mainPanel(

                   h2("Print the input data informaiton"),
                   # Action: pint the selected inputdata----
                   actionButton("printinfo","Print input data"),
                   verbatimTextOutput("inputfilesummary"),

                   # Horizontal line ----
                   tags$hr(),
                   h2("Merge climate and LAI input data"),
                   # Action: Merge Climate and LAI data----
                   actionButton("mergeinput","Merge inputs"),
                   verbatimTextOutput("mergeinputout"),

                   tags$hr(),
                   h2("Calculate water use efficiency using vegetation ratio"),
                   # Action: Calculate WUE based on vegetation ratio----
                   actionButton("calwue","Calculate WUE"),
                   verbatimTextOutput("wues"),

                   h2("Calculate ET using Sun's function"),
                   # Action: Calculate monthly Sun ET----
                   actionButton("calsunET","Calculate ET"),
                   verbatimTextOutput("ETs"),

                   h2("Transfering monthly data to daily data"),
                   # Action: Convert monthly input to daily input ----
                   actionButton("mon2day","Convert monthly input to daily"),
                   verbatimTextOutput("mon2dayout"),

                   h2("Calculate snow and effective rainfall using daily data"),
                   # Action: Calculate snow with daily input ----
                   actionButton("calsnow","Calculate snow"),
                   verbatimTextOutput("calsnowout"),

                   tags$hr(),
                   h2("Simulate water and carbon processes"),
                   # Input: Checkbox if calibate model ----
                   checkboxInput("Calibration", "Calibration", FALSE),
                  # Input: Simulation period
                   sliderInput("simulatedaterange", "Year simulation", 1970, 2017, value = c(2000, 2010)),

                  ## Action: run model
                  actionButton("runmodel","Run Simulation"),

                   tags$hr(),
                   h2("Plot water and carbon processes"),
                  # Input: Result plot variables ----
                  checkboxGroupInput("plotvars","Select the variables to plot",c(input_vars,result_vars[1:8])),

                  # Input: Year for plotting ---
                  sliderInput("plotdaterange", "Year plot", 1970, 2017, value = c(2000, 2010)),

                  # Action: Plot simulated result ---
                   actionButton("plotres","Plot"),

                  # Output: ploted simulated result
                   plotOutput("plotresout")
                 )
              )
        ),
    # Tab: Prepare input data ----
    tabPanel("Input preparation",
             sidebarLayout(
               sidebarPanel(
                 # Application title

                 h2("Upload watershed Shapefile"),
                 HTML("<p>Please chose a '*gml' file to upload.</p>
                      "),
                 fileInput('Input_basin', 'Choose watershed shapefile', multiple=FALSE, accept="gml"),

                 h2("Upload all original raster files"),
                 HTML("<p>Please chose a '*.tif' or '*.nc' monthly stacked precipitation and temperation file to upload.</p>
                      "),
                 fileInput('Input_dem_raster', 'Choose elevation raster file', multiple=FALSE, accept="tif"),
                 fileInput('Input_temp_raster', 'Choose temperature raster file', multiple=FALSE, accept="tif"),
                 fileInput('Input_precp_raster', 'Choose precipitation raster file', multiple=FALSE, accept="tif"),
                 fileInput('Input_lc_raster', 'Choose land cover raster file', multiple=FALSE, accept="tif"),
                 fileInput('Input_imp_raster', 'Choose impverious raster file', multiple=FALSE, accept="tif"),
                 fileInput('Input_lai_raster', 'Choose LAI raster file', multiple=FALSE, accept="tif"),

                 # Action: Plot basin ----
                 actionButton("processrasters","process input"),
                 # Action: Process input ----
                 actionButton("plotrasterdata","Plot the input rasters!")

               ),

               # MainPanel for upload data----
               mainPanel(
                 ## Output: Head of the selected plotting file----
                 # h2("The Map of the Watershed"),
                 # leafletOutput("basinmap"),


                 h2("Plot mean climate data"),
                 HTML("<p>It will take a while for the data process. Please wait ...</p>
                      "),
                 verbatimTextOutput("prntraster"),
                 leafletOutput("basinrastermap")
               )
             )
      ),
    navbarMenu("Export result",
               tabPanel("Export water",

                        # Sidebar layout with input and output definitions ----
                        sidebarLayout(

                          # Sidebar panel for inputs ----
                          sidebarPanel(

                            # Input: Choose dataset ----
                            selectInput("water", "Choose a dataset:",
                                        choices = c("rock", "pressure", "cars")),

                            # Button
                            downloadButton("downloadData", "Download")

                          ),

                          # Main panel for displaying outputs ----
                          mainPanel(

                            tableOutput("table")

                          )
                        )

               ),
               tabPanel("Export Climate",

                        # Sidebar layout with input and output definitions ----
                        sidebarLayout(

                          # Sidebar panel for inputs ----
                          sidebarPanel(

                            # Input: Choose dataset ----
                            selectInput("dataset1", "Choose a dataset:",
                                        choices = c("rock", "pressure", "cars")
                            ),

                            # Button
                            downloadButton("downloadData1", "Download")

                          ),

                          # Main panel for displaying outputs ----
                          mainPanel(

                            tableOutput("table1")

                          )
                        )

               )
    )
  )
)
