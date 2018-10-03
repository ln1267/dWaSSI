# dWaSSIC
This is a R version of dWaSSI-C model and some functions for processing spatial data
===================================================

***Disclaimer: currently under development***

## Installation

You can install the development version from GitHub with:

``` r
# load devtools package for installing github packages
  if("devtools" %in% installed.packages()){
    library(devtools)
  }else{
    install.packages("devtools")
    library(devtools)
  }

# Install dWaSSIC package from github
  if("dWaSSI" %in% installed.packages()){
    library(dWaSSI)
  }else{
    install_github("ln1267/dWaSSI")
    library(dWaSSI)
  }

# Load other required packages
f_lib_check(c("raster","reshape2","parallel","lubridate"))
```

## R Version

The original model code is written in Fortran and is publicly accessible. The R version is translated from the Fortran code developed by [Peter Caldwell](https://github.com/sungwookwi) and R code of Sac-sma developed by [Umit Taner](https://github.com/tanerumit/sacsmaR).

------------------------------------------------------------------------

### Description

WaSSI is an integrated, process-based model that can be used to project the effects of forest land cover change, climate change, and water withdrawals on river flows, water supply stress, and ecosystem productivity (i.e. carbon dynamics). WaSSI operates on a monthly time step at the HUC-4 (8-digit HUC) watershed scale (see more on HUCs) and across Mexico at the 0.5 degree scale. For the conterminous U.S., the model can also be run at the HUC12 scale for water and carbon balances from 1960 to 2012. As water yield and carbon sequestration are tightly coupled, WaSSI can be used to evaluate trade-offs among management strategies for these ecosystem services.

The web application for WaSSI allows users to define a custom simulation scenario, view/download model inputs and outputs in tabular and graphical form for a location of interest, and view/export model outputs spatially for a variety of time scales using an interactive map viewer. Users may select their location in the map viewer, select a specific HUC, or input a zip code to view model inputs and outputs



Further information is available at: [USDA - Forest  Service](https://www.fs.usda.gov/ccrc/tools/wassi) or the [Guide PDF]("WaSSIUserGuide_english_v1.1.pdf")

<p align="center">
<img src="data/Framework_WaSSIC.png">
</p>

------------------------------------------------------------------------

### The structure of original Fortran model ("src/")

Original fortran code is located in **"src"** folder:

GENERAL.F90 
	This is the main file for the whole program. It contants the main 3 nested do loops for grids, year and month. Besides, the year and month loops should be run orderedly, and they are relevant to previous loop. 
	It reads input data from 5 TXT files where colums are separated by ",".
	
Subroutines:  
	
WARMUP.F90  
	This subroutine is mainly used for reading data from input files and write some basic runnning information to ouput file.

PET.F90   
	It is used for simulating PAET(year,month) for the following simulation.(As PAET did not include grid information, We may be need to change it in the future.)
	There is no relevant in this matrix
	
WATERBAL.F90, WATERBAL_LC.F90 and WATERBAL_MON_LC.F90   
	This is the main file for data processing.
	It contants a day loop for the following simulation. The struction details can be found in the "Waterbal.png".( There is a bit differece between this framwork and the code, but the main calculation part is the same. The code does not have the K loop)" 
	
SUMMARY.F90  
	This is similar to output file, but it only have year loop.
	
OUTPUT.F90   
	This is used for summary the main results and write them to outputs files.
	It has two nested loops. (Year and month loop)

------------------------------------------------------------------------

### How to use Makefile to build model
To compile first make sure that the Makefile is updated for the machine that you are compiling on.

```bash
make
```
This should produce some objects with extensions '*.o' and *a.out* which is your executable.

To clean the project just type:
```bash
make clean 
```
and it should delete all except the source files. 

------------------------------------------------------------------------

### Included functions

The package consists of five functions:

``` r
"Wa"
```
