!**********************************************************************C
!     *** SUBROUTINE WATERBAL ***                                      C
!     SIMULATES MONTHLY WATER BALANCE USING 2 LAYER SOIL MOISTURE      C
!     ALGORITHM FROM NOAA NATIONAL WEATHER SERVICE SACRAMENTO SOIL     C
!     MOISTURE ACCOUNTING MODEL (SAC-SMA)                              C
!     Line 941 --- Area for each cell                                  C
!**********************************************************************C
      
SUBROUTINE WATERBAL_MON(I,J_S,M)
    Use Common_var
    implicit none       
! ----------------------------------------------------------------------------     
         
    INTEGER I,J,M,J_S
            
    REAL INFIL, EXCESS,SNOWW,RAINPPT,SNOWPPT
          
    REAL BF,SBF,SIF,SSUR,TWX
    
    REAL DEFR,DEL,HPL,ETLZTW,ETUZFW,ETUZTW,ET,RESIDET
    
    REAL GEP,RECO,NEE
    
    REAL PERC,PERCF,PERCM,PERCP,PERCS,PERCT
    
    REAL RATLP,RATLS,RATLZ,RATLZT,UZRAT,FRACP,LZDEF,SPERC
    
    REAL DINC,DLZP,DLZS,DUZ,MELTMAX,PINC,SMF,RAINTEMP,SNOWTEMP        
           
    INTEGER NINC,F,LC_N
       
! *****************************************************************************************************
!----Set the simulate ID for the start year
    J=J_S+IYSTART-1-NWARMUP
    
    LC_N=LADUSE(I)
! *****************************************************************************************************
! INITIAL  SOIL WATER CONDITION                
    IF (J_S .EQ. 1 .AND. M .EQ. 1) THEN
 
        LZFSC=0.0
        UZTWC = UZTWM(I)
        UZFWC = 0.0
        LZTWC = LZTWM(I)
        LZFSC = 0.75*LZFSM(I)
        LZFPC = 0.75*LZFPM(I)
    ENDIF 
         
! *****************************************************************************************************
! *****************************************************************************************************
! --ESTIMATE INITIAL SNOWPACK (DECEMBER OF  YEAR BEFORE FIRST YEAR OF SIMULATION) 
! --AS A FUNCTION OF LATITUDE (DECIMAL DEGREES)
! --INITIAL SNOWPACK = 0.0002*EXP(0.2644*LATUDE(I))
! --R2=0.42, N=1559
! --BASED ON MEASURED NSIDC 12/1999 SNOWPACK
    SNOWPACK = 0.0002*EXP(0.2644*LATUDE(I))
    
    RAINTEMP = 2.4
    
    SNOWTEMP = -6.8
    
    MELTMAX = 0.85
        
! ---- FOR TEMPERATURE LESS THAN SNOW TEMP, CALCULATE SNOWPPT

    IF (TEMP(I,J, M).LE.SNOWTEMP) THEN
        
        SNOWPPT = RAIN(I,J, M)
             
! ---- FOR TEMPERATURE GREATER THAN RAIN TEMP, CALCULATE RAINPPT

    ELSEIF (TEMP(I,J, M).GE.RAINTEMP) THEN
        
        RAINPPT = RAIN(I,J, M)
        
! ---- FOR TEMPERATURE BETWEEN RAINTEMP AND SNOWTEMP, CALCULATE SNOWPPT AND RAINPPT

    ELSE
        
        SNOWPPT = RAIN(I,J, M)*((RAINTEMP - TEMP(I,J, M))/(RAINTEMP - SNOWTEMP))
           
        RAINPPT = RAIN(I,J, M) - SNOWPPT
           
    ENDIF
! ---- ACCUMULATE SNOWPPT AS SNOWPACK
           
    SNOWPACK = SNOWPACK + SNOWPPT

! ---- CALCULATE SNOW MELT FRACTION BASED ON MAXIMUM MELT RATE (MELTMAX) AND MONTHLY TEMPERATURE

    SMF = ((TEMP(I,J, M)-SNOWTEMP)/(RAINTEMP - SNOWTEMP))*MELTMAX
        
    IF (SMF .GT. MELTMAX) THEN
        
        SMF = MELTMAX
           
    ENDIF
        
    IF (SMF .LT. 0.) THEN
        
    SMF = 0.
    
    ENDIF
      
! ---- CALCULATE AMOUNT OF SNOW MELTED (MM) TO CONTRIBUTE TO INFILTRATION & RUNOFF

! ---- IF SNOWPACK IS LESS THAN 10.0 MM, ASSUME IT WILL ALL MELT (MCCABE AND WOLOCK, 1999)
! ---- (GENERAL-CIRCULATION-MODEL SIMULATIONS OF FUTURE SNOWPACK IN THE WESTERN UNITED STATES)

    SNOWW = SNOWPACK * SMF
        
    SNOWPACK = SNOWPACK - SNOWW        

    IF (SNOWPACK .LT. 10.0) THEN
        
        SNOWW = SNOWW+SNOWPACK
        
        SNOWPACK = 0.0
          
    ENDIF      
      
! -- COMPUTE THE INFILTRATION FOR A GIVEN MONTH FOR EACH LAND USE

    INFIL=RAINPPT+SNOWW      
      
! *****************************************************************************************************
! *****************************************************************************************************
! *****************************************************************************************************

! *****************************************************************************************************
! -- SET ET, SURFACE RUNOFF, INTERFLOW, GEP TO ZERO IF TEMPERATURE IS LE -1.0
! -- BASEFLOW STILL OCCURS 
! --- INITIALIZE FLOWS
    SBF=0.0
    SSUR=0.0
    SIF=0.0
    SPERC=0.0 
    
! -- FOR OPEN WATER, ET=PAET=PET, INFILTRATION IN EXCESS OF ET GOES TO SURFACE RUNOFF

      ! IF(LC_N.EQ.8) THEN
                   
       ! ET = PAET(I,J,M)
                      
      ! SSUR=INFIL-ET
    
      ! ENDIF   

  
! *****************************************************************************************************
! --- COMPUTE AET GIVEN TOTAL WATER STORED IN UPPER SOIL LAYER STORAGES AND PAET CALCULATED IN PET.FOR
! --- ASSUME ET IS SUPPLIED ONLY FROM UPPER LAYER NO UPWARD FLUX FROM LOWER LAYER TO UPPER LAYER
! --- NOTE THAT SAC-SMA ALLOWS ET TO ALSO BE SUPPLIED UNRESTRICTED BY LZ TENSION WATER STORAGE
               
    ET = PAET(I,J,M)
                
! --- COMPUTE ET FROM UZ TENSION WATER STORAGE, RECALCULATE UZTWC, CALCULATE RESIDUAL ET DEMAND

    ETUZTW = ET * (UZTWC/UZTWM(I))
                   
    RESIDET = ET - ETUZTW
                   
    UZTWC = UZTWC - ETUZTW
                   
    ETUZFW = 0.0
                   
    IF (UZTWC.GE.0.0) GOTO 220
    
    ETUZTW = ETUZTW + UZTWC
                   
    UZTWC = 0.0
                   
    RESIDET = ET - ETUZTW
                   
! --- COMPUTE ET FROM UZ FREE WATER STORAGE, RECALCULATE UZFWC, CALCULATE RESIDUAL ET DEMAND                   
                   
    IF (UZFWC .GE. RESIDET) GO TO 221
                   
    ETUZFW = UZFWC
                   
    UZFWC = 0.0
                   
    RESIDET = RESIDET - ETUZFW
                   
    GO TO 225
                   
221 ETUZFW = RESIDET

    UZFWC = UZFWC - ETUZFW
                   
    RESIDET = 0.0
                   
! --- REDISTRIBUTE WATER BETWEEN UZ TENSION WATER AND FREE WATER STORAGES

220 IF((UZTWC/UZTWM(I)).GE.(UZFWC/UZFWM(I)))  GO TO 225

    UZRAT=(UZTWC+UZFWC)/(UZTWM(I)+UZFWM(I))
                      
    UZTWC = UZTWM(I) * UZRAT
                      
    UZFWC = UZFWM(I) * UZRAT
                        
225 IF (UZTWC .LT. 0.00001) UZTWC = 0.0

    IF (UZFWC .LT. 0.00001) UZFWC = 0.0
                   
! --- COMPUTE ET FROM LZ TENSION WATER STORAGE, RECALCULATE LZTWC, CALCULATE RESIDUAL ET DEMAND

    ETLZTW = RESIDET * (LZTWC / &
        (UZTWM(I) + LZTWM(I)))
                   
    LZTWC = LZTWC - ETLZTW
                   
    IF(LZTWC .GE. 0.0) GO TO 226
                   
    ETLZTW = ETLZTW + LZTWC
                   
    LZTWC = 0.0
                   
226 RATLZT = LZTWC / LZTWM(I)

    RATLZ = (LZTWC + LZFPC + LZFSC) /&
                   (LZTWM(I) + LZFPM(I) + LZFSM(I))
     
    IF (RATLZT .GE. RATLZ) GO TO 230
                  
    LZTWC = LZTWC + (RATLZ - RATLZT) * LZTWM(I)

    LZFSC = LZFSC - (RATLZ - RATLZT) * LZTWM(I)
                   
    IF(LZFSC .GE. 0.0) GO TO 230
                   
    LZFPC = LZFPC + LZFSC
                   
    LZFSC = 0.0
                   
230 IF (LZTWC .LT. 0.00001) LZTWC = 0.0

! --- CALCULATE TOTAL ET SUPPLIED BY UPPER AND LOWER LAYERS

    ET = ETUZTW + ETUZFW + ETLZTW
                   
    IF (ET .LE. 0)  ET=0.0 
    
! *****************************************************************************************************
! *****************************************************************************************************
! --- COMPUTE PERCOLATION INTO SOIL WATER STORAGES AND SURFACE RUNOFF

!    --- COMPUTE WATER IN EXCESS OF UZ TENSION WATER CAPACITY (TWX)

    TWX = INFIL + UZTWC - UZTWM(I)
           
    IF (TWX.GE.0.0) THEN
                 
!     --- IF INFIL EXCEEDS UZ TENSION WATER CAPACITY, SET UZ TENSION WATER STORAGE TO CAPACITY
        UZTWC = UZTWM(I)     
                
    ELSE
        
!     --- IF INFIL DOES NOT EXCEED UZ TENSION WATER CAPACITY, ALL INFIL GOES TO UZ TENSION WATER STORAGE
        UZTWC = UZTWC + INFIL
                     
        TWX=0.0
                     
    ENDIF

!-------------------------------------------------------------
!************************************************************
!------------------------------------------------------------------      
                  
! --- DETERMINE COMPUTATIONAL TIME INCREMENTS FOR BASIC TIME INTERVAL

    NINC=INT(1.0+0.2*(UZFWC+TWX))
   
!     NINC=NUMBER OF TIME INCREMENTS THAT THE TIME INTERVAL
!     IS DIVIDED INTO FOR FURTHER SOIL-MOISTURE ACCOUNTING.  
!     NO ONE INCREMENT WILL EXCEED 5.0 MILLIMETERS OF UZFWC+PAV

    DINC=(1.0/NINC)
!    DINC IS THE LENGTH OF EACH INCREMENT IN MONTHS

    PINC=TWX/NINC
!    PINC IS THE AMOUNT OF AVAILABLE MOISTURE FOR EACH INCREMENT

! --- COMPUTE FREE WATER DEPLETION COEFICIENTS FOR THE TIME INCREMENT BEING USED

    ! DUZ=UZK(I)*DINC
    ! DLZP=LZPK(I)*DINC
    ! DLZS=LZSK(I)*DINC
    
    DUZ  = 1 - (1 - UZK(I))**DINC
    DLZP = 1 - (1 - LZPK(I))**DINC
    DLZS = 1 - (1 - LZSK(I))**DINC 

!if (I==1.and.J==45) print*,"I=",I,'J=',J,'M=',M,"middle LZFSC=",LZFSC
          
! --- SET UP INCREMENTAL DO LOOP FOR THE TIME INTERVAL
      
    DO 450 F=1,NINC
      
      EXCESS=0.
      
! --  COMPUTE BASEFLOW AND KEEP TRACK OF TIME INTERVAL SUM.

      BF = LZFPC * DLZP
                
      LZFPC = LZFPC - BF
                
      IF (LZFPC .LE. 0.0001) THEN 
                
      BF = BF + LZFPC
                
      LZFPC = 0.0
                
      ENDIF

      SBF = SBF + BF
                 
      BF = LZFSC * DLZS

      LZFSC = LZFSC - BF
                
      IF (LZFSC .LE. 0.0001) THEN
                
      BF = BF + LZFSC
                
      LZFSC = 0.0
                
      ENDIF         
                 
      SBF = SBF + BF       
                 

! --- COMPUTE PERCOLATION TO LZ IF FREE WATER IS AVAILABLE IN UZ

      IF (PINC + UZFWC.GT.0.01) GOTO 460
          
      UZFWC = UZFWC + PINC
          
      GOTO 450

!     --- COMPUTE PERCOLATION DEMAND FROM LZ

460   PERCM = LZFPM(I) * DLZP + LZFSM(I) * DLZS
                
      PERC = PERCM * (UZFWC/UZFWM(I))
                
      DEFR=1.0-((LZTWC+LZFPC+LZFSC)/(LZTWM(I)+LZFPM(I)+LZFSM(I)))
                
      PERC = PERC * (1.0 + ZPERC(I) * (DEFR**REXP(I)))

!     --- COMPARE LZ PERCOLATION DEMAND TO UZ FREE WATER AVAILABLE AND COMPUTE ACTUAL PERCOLATION

      IF (PERC .GE. UZFWC) THEN
                
      PERC = UZFWC
                       
      ENDIF
                       
      UZFWC = UZFWC - PERC
                
!      --- CHECK TO SEE IF PERC EXCEEDS LZ TENSION AND FREE WATER DEFICIENCY, IF SO SET PERC TO LZ DEFICIENCY

      LZDEF = (LZTWC + LZFPC + LZFSC) - (LZTWM(I) + LZFPM(I) + LZFSM(I)) + PERC
                    
      IF (LZDEF .GT. 0.0) THEN
                    
      PERC = PERC - LZDEF
          
      UZFWC = UZFWC + LZDEF
                       
      ENDIF
                    
      SPERC = SPERC + PERC
!        SPERC IS THE TIME INTERVAL SUM OF PERC

! --- COMPUTE INTERFLOW FROM UZ

      DEL = UZFWC * DUZ
      
      IF (DEL.GT.UZFWC)THEN
      
      DEL=UZFWC
      
      UZFWC=0.
      
      ELSE
      
      UZFWC = UZFWC - DEL
      
      ENDIF
                 
      SIF = SIF + DEL

! --- DISRIBUTE PERCOLATED WATER INTO THE LZ STORAGES AND COMPUTE THE REMAINDER IN UZ FREE WATER STORAGE AND RESIDUAL AVAIL FOR RUNOFF
   
!     --- COMPUTE PERC WATER GOING INTO LZ TENSION WATER STORAGE AND COMPARE TO AVAILABLE STORAGE

      PERCT = PERC * (1.0 - PFREE(I))
                
      IF ((PERCT + LZTWC) .GT. LZTWM(I)) THEN
                
!     --- WHEN PERC IS GREATER THAN AVAILABLE TENSION WATER STORAGE, SET TENSION WATER STORAGE TO MAX, REMAINDER OF PERC GETS EVALUATED AGAINST FREE WATER STORAGE

      PERCF = PERCT + LZTWC - LZTWM(I)
                
      LZTWC = LZTWM(I)
                
      ELSE
                
!     --- WHEN PERC IS LESS THAN AVAILABLE TENSION WATER STORAGE, UPDATE TENSION WATER STORAGE

      LZTWC = LZTWC + PERCT
                
      PERCF = 0.0
                
      ENDIF
                
!     --- COMPUTE TOTAL PERC WATER GOING INTO LZ FREE WATER STORAGE

      PERCF = PERCF + PERC * PFREE(I)                

      IF(PERCF .EQ. 0.0) GOTO 470
                
!     --- COMPUTE RELATIVE SIZE OF LZ PRIMARY FREE WATER STORAGE COMPARED TO LZ TOTAL FREE WATER STORAGE

      HPL = LZFPM(I) / (LZFPM(I) + LZFSM(I))
                
!     --- COMPUTE LZ PRIMARY AND SECONDARY FREE WATER CONTENT TO CAPACITY RATIOS

      RATLP = LZFPC / LZFPM(I)
                
      RATLS = LZFSC / LZFSM(I)
                
!    --- COMPUTE FRACTIONS AND PERCENTAGES OF FREE WATER PERC TO GO TO LZ PRIMARY STORAGE

    FRACP = (HPL * 2.0 * (1.0 - RATLP)) / ((1.0 - RATLP) + (1.0 - RATLS))
                
    IF (FRACP .GT. 1.0) FRACP = 1.0

    PERCP = PERCF * FRACP
                
    PERCS = PERCF - PERCP
                
!    --- COMPUTE NEW PRIMARY AND SUPPLEMENTAL STORAGE

!    --- COMPUTE NEW SUPPLEMENTAL FREE WATER STORAGE

      LZFSC = LZFSC + PERCS

      IF(LZFSC .GT. LZFSM(I)) THEN
                
!    --- IF NEW SUPPLEMENTAL FREE WATER STORAGE EXCEEDS CAPACITY SET SUPPLEMENTAL STORAGE TO CAPACITY AND EXCESS GOES TO PRIMARY FREE WATER STORAGE

      PERCS = PERCS - LZFSC + LZFSM(I)
                
      LZFSC = LZFSM(I)
                          
      ENDIF
                
            
!    --- IF NEW LZ SUPPLEMENTAL FREE WATER STORAGE IS LESS THAN CAPACITY MOVE ON TO COMPUTE NEW PRIMARY FREE WATER STORAGE

      LZFPC = LZFPC + (PERCF - PERCS)
                
      IF (LZFPC .LE. LZFPM(I)) GOTO 470

!    --- IF LZ FREE PRIMARY WATER STORAGE EXCEEDS CAPACITY SET PRIMARY STORAGE TO CAPACITY AND EVALUATE EXCESS AGAINST LZ TENSION WATER STORAGE

      LZFPC = LZFPC + (PERCF - PERCS)

      IF (LZFPC .LE. LZFPM(I)) GOTO 470

!    --- IF LZ FREE PRIMARY WATER STORAGE EXCEEDS CAPACITY SET PRIMARY STORAGE TO CAPACITY AND EVALUATE EXCESS AGAINST LZ TENSION WATER STORAGE

      LZTWC = LZTWC + LZFPC - LZFPM(I)
                
      LZFPC = LZFPM(I)
      
! -- CHECK TO SEE IF LZTWC>LZTWM, IF SO, SET LZTWC=LZTWM AND SEND EXCESS FOR INCREMENT TO SURFACE RUNOFF

      IF(LZTWC.GT.LZTWM(I)) THEN
      
      EXCESS=LZTWC-LZTWM(I)
      
      LZTWC=LZTWM(I)
      
      ELSE
      
      EXCESS=0.
      
      ENDIF
                
!    ---DISTRIBUTE PINC BETWEEN UZFWC AND SURFACE RUNOFF

!    IF UZFWC CAN ABSORB EXCESS, NO SURFACE RUNOFF

470   IF(PINC+EXCESS.EQ.0.0) GOTO 450

      IF (PINC+UZFWC+EXCESS.GT.UZFWM(I)) GOTO 480
                       
      UZFWC=UZFWC+PINC+EXCESS
                          
      GOTO 450
                          
!    IF UZFWC CAN NOT ABSORB EXCESS, THE COMPUTE SURFACE RUNOFF
                       
480   SSUR=SSUR+PINC+UZFWC-UZFWM(I)+EXCESS
                          
      UZFWC=UZFWM(I)
      
450 CONTINUE

! if (I==1.and.J==45)       print*,BF,  SBF
 
! *****************************************************************************************************
! *****************************************************************************************************
! Calculate GEP based on ET and the equation
    IF(LC_N .LE. 0 ) THEN
        GEP=0.0
        NEE=0.0
        RECO=0.0
        GOTO 50201
    ENDIF
    
    GEP = wue_k(LC_N) * ET 

    IF (GEP .LE. 0 .OR. TEMP(I,J,M) .LE. -1.0)  GEP=0.0 

    RECO= reco_inter(LC_N) + reco_slope(LC_N) * GEP
        
    NEE=  RECO- GEP      

! *****************************************************************************************************
! --- CALCULATE THE FRACTION OF Monthly AET
50201    AET(I,J,M) = ET
               
! *****************************************************************************************************
! -- CALCULATE THE FRACTION OF Monthly SURFACE RUNOFF
    RUNOFF(I,J,M) = SSUR
    
    IF (RUNOFF(I,J,M) .LT. 0.)  RUNOFF(I,J,M)=0.    
    
! *****************************************************************************************************
! -- CALCULATE THE FRACTION OF Monthly SECONDARY BASEFLOW
    SECBF(I,J,M) = SBF
! *****************************************************************************************************
! -- CALCULATE THE FRACTION OF Monthly INTERFLOW
    INTF(I,J,M) = SIF
                   
    GEPM(I,J, M) = GEP
    RECOM(I,J,M)  = RECO
    NEEM(I,J,M) = NEE           
           
    SP(I,J,M) = SNOWPACK
    AVUZTWC(I,J,M) = UZTWC
    AVUZFWC(I,J,M) = UZFWC
    AVLZTWC(I,J,M) = LZTWC
    AVLZFPC(I,J,M) = LZFPC
    AVLZFSC(I,J,M) = LZFSC
    SMC(I,J,M) = UZTWC+UZFWC+LZTWC+LZFPC+LZFSC

! -- STREAMFLOW IN MILLION M3 FOR EACH HUC FOR MONTH M. HUCAREA IN SQ. METERS

    STRFLOW(I, J, M) = (RUNOFF(I,J,M) + SECBF(I,J,M) + INTF(I,J,M)) /1000./1000000.
 
!if (I==1 .and. J==45) print*,I,J,M,RUNOFF(I,J,M), SECBF(I,J,M), INTF(I,J,M), HUCAREA(I)
    RETURN
    
END
