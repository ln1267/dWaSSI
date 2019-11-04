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
                 fileInput("Input_cellinfo", "Cellinfo [HRUID, Latitude, Longitude, Elevation_m, Area_m2 and ratio for each lc]"),
                 fileInput("Input_climate", "Monthly Climate [HRUID, Year, Month, Precip_mm, Tmean_C]"),
                 fileInput("Input_LAI", "Monthly LAI [HRUID, Year, Month, and LAI for each lc]"),
                 fileInput("Input_soilinfo", "Soil info [HRUID, uztwm,uzfwm,uzk,zperc,rexp,lztwm,lzfsm,lzfpm,lzsk,lzpk,pfree]")

               ),

               # MainPanel for upload data----
               mainPanel(
                 ## Output: Summary of the selected plotting file----
                 h2("The summary of the input file"),
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
                 actionButton("plotrasterdata","Plot the input data!"),
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
                             choices=c("Climate","LAI")),
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
                 verbatimTextOutput("distPloterror"),

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
                 textInput("StationID", "Type the BasinID of the hydrologic station!","1"
                 ),

                 # # Input: Time range for subsetting ----
                 dateRangeInput("dateSimulation","Select the date range for simulation"),

                 # Action: Subset data ----
                 actionButton("subSimData","Subset the selected dataset"),

                 # Input: Simulation period
                 #sliderInput("dateSimulation", "Year simulation", 1970, 2017, value = c(2000, 2010)),
                 # Input: Checkbox if calibate model ----
                 checkboxInput("Calibration", "Calibration", FALSE),
                 ## Action: run model
                 actionButton("runSimulation","Run Simulation"),
                 tags$hr(),
                 tags$h5("Simulation log:"),
                 verbatimTextOutput("simulating")
                  ),

                 # Mainpanel: of Simulation----
                 mainPanel(


                   h2("Plot water and carbon processes"),
                  # Input: Result plot variables ----
                  checkboxGroupInput("plotvars","Select the variables to plot",c("P","T","Q"),inline = T),

                  # Input: Year for plotting ---
                  sliderInput("plotdaterange", "Year plot", 1970, 2017, value = c(2000, 2010)),

                  # Action: Plot simulated result ---
                   actionButton("plotres","Plot"),

                  # Output: ploted simulated result
                   plotOutput("plotresout")
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
