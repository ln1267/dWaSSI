!**********************************************************************C
!                                                                      C
!     *** SUBROUTINE WATERBAL ***                                      C
!     SIMULATES MONTHLY WATER BALANCE USING 2 LAYER SOIL MOISTURE      C
!     ALGORITHM FROM NOAA NATIONAL WEATHER SERVICE SACRAMENTO SOIL     C
!     MOISTURE ACCOUNTING MODEL (SAC-SMA)                              C
!     Line 591 --- Carbon model                                        C
!     Line 941 --- Area for each cell                                  C
!**********************************************************************C
      
      SUBROUTINE WATERBAL_MON_LC(I,J_S,M)
        Use Common_var
        implicit none       
! ----------------------------------------------------------------------------     
         
      INTEGER I,J,M,K,J_S
      
      REAL AETTEMP, RUNOFFTEMP, PBFTEMP, SBFTEMP,IFTEMP, GEPTEMP,&
            RECOTEMP, NEETEMP,INTFTEMP
      
      REAL SNOWPACK, SNOWW, INFIL, EXCESS
      
      REAL TAREA
      
      REAL TAUZTWC, TAUZFWC, TALZTWC, TALZFPC, TALZFSC
     
      REAL TASM,AUZTWC,AUZFWC,ALZTWC,ALZFPC,ALZFSC,ASM
      REAL DINC,DLZP,DLZS,DUZ,MELTMAX,PINC,SMF,RAINTEMP,SNOWTEMP,RAINPPT,SNOWPPT        
           
      INTEGER NINC,F
      
    ! REAL :: RUNLAND(NGRID,NYEAR_S+NWARMUP,12,31,NLC)
    ! REAL :: ETLAND(NGRID,NYEAR_S+NWARMUP,12,31,NLC)
    ! REAL :: GEPLAND(NGRID,NYEAR_S+NWARMUP,12,31,NLC) 
           
! *****************************************************************************************************
!----Set the simulate ID for the start year
    J=J_S+IYSTART-1-NWARMUP

! *****************************************************************************************************
                  
        IF (J_S .EQ. 1 .AND. M .EQ. 1) THEN
 
        LZFSC_lc=0.0
        DO 50 K=1, NLC
                        
           UZTWC_lc(K) = UZTWM(I)
           UZFWC_lc(K) = 0.0
           LZTWC_lc(K) = LZTWM(I)
           LZFSC_lc(K) = 0.75*LZFSM(I)
           LZFPC_lc(K) = 0.75*LZFPM(I)
           SURFRO_lc(K)=0
50      CONTINUE
        
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
! --- INITIALIZE VARIABLES FOR START OF SIMULATION
             AETTEMP =0.
             RUNOFFTEMP = 0.
             PBFTEMP = 0.0
             SBFTEMP = 0.0
             IFTEMP = 0.0
             INTFTEMP=0.0
                          
            TAREA = 0.

            TASM = 0. 
            ASM = 0. 
            TAUZTWC=0.
            TAUZFWC=0.
            TALZTWC=0.
            TALZFPC=0.
            TALZFSC=0.           
              
            GEPTEMP=0.
            RECOTEMP=0.
            NEETEMP=0.

! *****************************************************************************************************
! *****************************************************************************************************
! -- LOOP THROUGH LAND COVERS IN THE HUC AND PERFORM WATER BALANCE COMPUTATIONS
! -- ASSUMES OPEN WATER LAND COVER IS NEG.
          
             DO 40 K=1, NLC
                  
! *****************************************************************************************************
! -- SET ET, SURFACE RUNOFF, INTERFLOW, GEP TO ZERO IF TEMPERATURE IS LE -1.0
! -- BASEFLOW STILL OCCURS 
! --- INITIALIZE FLOWS
                 SBF_lc(K)=0.0
                 SSUR_lc(K)=0.0
                 SIF_lc(K)=0.0
                 SPERC_lc(K)=0.0 
!if (I==1.and.J==45) print*,"I=",I,'J=',J,'M=',M,'K=',K,"start LZFSC_lc=",LZFSC_lc(K),"AET=",AET_lc(J, M, K)
! -- FOR OPEN WATER, ET=PAET=PET, INFILTRATION IN EXCESS OF ET GOES TO SURFACE RUNOFF

      ! IF(K.EQ.8) THEN
                   
       ! ET_lc(J, M, K) = AET_lc(J, M, K)
                      
      ! SSUR_lc(K)=INFIL-ET_lc(J, M, K)
    
      ! ENDIF   

  
! *****************************************************************************************************
! --- COMPUTE AET GIVEN TOTAL WATER STORED IN UPPER SOIL LAYER STORAGES AND PAET CALCULATED IN PET.FOR
! --- ASSUME ET IS SUPPLIED ONLY FROM UPPER LAYER NO UPWARD FLUX FROM LOWER LAYER TO UPPER LAYER
! --- NOTE THAT SAC-SMA ALLOWS ET TO ALSO BE SUPPLIED UNRESTRICTED BY LZ TENSION WATER STORAGE
               
                   ET_lc(J, M, K) = AET_lc(J, M, K)
                
! --- COMPUTE ET FROM UZ TENSION WATER STORAGE, RECALCULATE UZTWC, CALCULATE RESIDUAL ET DEMAND

                   ETUZTW_lc(J, M, K) = ET_lc(J, M, K) * (UZTWC_lc(K)/UZTWM(I))
                   
                   RESIDET_lc(J, M, K) = ET_lc(J, M, K) - ETUZTW_lc(J, M, K)
                   
                   UZTWC_lc(K) = UZTWC_lc(K) - ETUZTW_lc(J, M, K)
                   
                   ETUZFW_lc(J, M, K) = 0.0
                   
                   IF (UZTWC_lc(K).GE.0.0) GOTO 220
                   
                   ETUZTW_lc(J, M, K) = ETUZTW_lc(J, M, K) + UZTWC_lc(K)
                   
                   UZTWC_lc(K) = 0.0
                   
                   RESIDET_lc(J, M, K) = ET_lc(J, M, K) - ETUZTW_lc(J, M, K)
                   
! --- COMPUTE ET FROM UZ FREE WATER STORAGE, RECALCULATE UZFWC, CALCULATE RESIDUAL ET DEMAND                   
                   
                   IF (UZFWC_lc(K) .GE. RESIDET_lc(J, M, K)) GO TO 221
                   
                   ETUZFW_lc(J, M, K) = UZFWC_lc(K)
                   
                   UZFWC_lc(K) = 0.0
                   
                   RESIDET_lc(J, M, K) = RESIDET_lc(J, M, K) - ETUZFW_lc(J, M, K)
                   
                   GO TO 225
                   
221                ETUZFW_lc(J, M, K) = RESIDET_lc(J, M, K)

                   UZFWC_lc(K) = UZFWC_lc(K) - ETUZFW_lc(J, M, K)
                   
                   RESIDET_lc(J, M, K) = 0.0
                   
! --- REDISTRIBUTE WATER BETWEEN UZ TENSION WATER AND FREE WATER STORAGES

220                IF((UZTWC_lc(K)/UZTWM(I)).GE.(UZFWC_lc(K)/UZFWM(I)))  GO TO 225

                   UZRAT_lc(K)=(UZTWC_lc(K)+UZFWC_lc(K))/(UZTWM(I)+UZFWM(I))
                      
                   UZTWC_lc(K) = UZTWM(I) * UZRAT_lc(K)
                      
                   UZFWC_lc(K) = UZFWM(I) * UZRAT_lc(K)
                        
225                IF (UZTWC_lc(K) .LT. 0.00001) UZTWC_lc(K) = 0.0

                   IF (UZFWC_lc(K) .LT. 0.00001) UZFWC_lc(K) = 0.0
                   
                   
! --- COMPUTE ET FROM LZ TENSION WATER STORAGE, RECALCULATE LZTWC, CALCULATE RESIDUAL ET DEMAND

                   ETLZTW_lc(J, M, K) = RESIDET_lc(J, M, K) * (LZTWC_lc(K) / &
                  (UZTWM(I) + LZTWM(I)))
                   
                   LZTWC_lc(K) = LZTWC_lc(K) - ETLZTW_lc(J, M, K)
                   
                   IF(LZTWC_lc(K) .GE. 0.0) GO TO 226
                   
                   ETLZTW_lc(J, M, K) = ETLZTW_lc(J, M, K) + LZTWC_lc(K)
                   
                   LZTWC_lc(K) = 0.0
                   
226                RATLZT_lc(K) = LZTWC_lc(K) / LZTWM(I)

                   RATLZ_lc(K) = (LZTWC_lc(K) + LZFPC_lc(K) + LZFSC_lc(K)) /&
                   (LZTWM(I) + LZFPM(I) + LZFSM(I))
     
                   IF (RATLZT_lc(K) .GE. RATLZ_lc(K)) GO TO 230
                  
                   LZTWC_lc(K) = LZTWC_lc(K) + (RATLZ_lc(K) - RATLZT_lc(K)) * &
                  LZTWM(I)

                   
                   LZFSC_lc(K) = LZFSC_lc(K) - (RATLZ_lc(K) - RATLZT_lc(K)) * &
                  LZTWM(I)
                   
                   IF(LZFSC_lc(K) .GE. 0.0) GO TO 230
                   
                   LZFPC_lc(K) = LZFPC_lc(K) + LZFSC_lc(K)
                   
                   LZFSC_lc(K) = 0.0
                   
230                IF (LZTWC_lc(K) .LT. 0.00001) LZTWC_lc(K) = 0.0

! --- CALCULATE TOTAL ET SUPPLIED BY UPPER AND LOWER LAYERS

                   ET_lc(J, M, K) = ETUZTW_lc(J, M, K) + ETUZFW_lc(J, M, K) + &
                  ETLZTW_lc(J, M, K)
                   
                    IF (ET_lc(J,M,K) .LE. 0)   THEN
                        ET_lc(J,M,K)=0.0 
                    ENDIF     
! *****************************************************************************************************
! *****************************************************************************************************
! --- COMPUTE PERCOLATION INTO SOIL WATER STORAGES AND SURFACE RUNOFF

!     --- COMPUTE WATER IN EXCESS OF UZ TENSION WATER CAPACITY (TWX)

                  TWX_lc(K) = INFIL + UZTWC_lc(K) - UZTWM(I)
           
                  IF (TWX_lc(K).GE.0.0) THEN
                 
!     --- IF INFIL EXCEEDS UZ TENSION WATER CAPACITY, SET UZ TENSION WATER STORAGE TO CAPACITY

                     UZTWC_lc(K) = UZTWM(I)     
                
                  ELSE
        
!     --- IF INFIL DOES NOT EXCEED UZ TENSION WATER CAPACITY, ALL INFIL GOES TO UZ TENSION WATER STORAGE

                     UZTWC_lc(K) = UZTWC_lc(K) + INFIL
                     
                     TWX_lc(K)=0.0
                     
                  ENDIF

!-------------------------------------------------------------
!************************************************************
!------------------------------------------------------------------      
                  
! --- DETERMINE COMPUTATIONAL TIME INCREMENTS FOR BASIC TIME INTERVAL

      NINC=INT(1.0+0.2*(UZFWC_lc(K)+TWX_lc(K)))

      
!     NINC=NUMBER OF TIME INCREMENTS THAT THE TIME INTERVAL
!     IS DIVIDED INTO FOR FURTHER SOIL-MOISTURE ACCOUNTING.  
!     NO ONE INCREMENT WILL EXCEED 5.0 MILLIMETERS OF UZFWC+PAV

      DINC=(1.0/NINC)
!                 DINC IS THE LENGTH OF EACH INCREMENT IN MONTHS

      PINC=TWX_lc(K)/NINC
!                 PINC IS THE AMOUNT OF AVAILABLE MOISTURE FOR EACH INCREMENT


! --- COMPUTE FREE WATER DEPLETION COEFICIENTS FOR THE TIME INCREMENT BEING USED

      DUZ=UZK(I)*DINC
      
      DLZP=LZPK(I)*DINC
      
      DLZS=LZSK(I)*DINC
                  

! --- SET UP INCREMENTAL DO LOOP FOR THE TIME INTERVAL
!if (I==1.and.J==45) print*,"I=",I,'J=',J,'M=',M,'K=',K,"middle LZFSC_lc=",LZFSC_lc(K)
      DO 450 F=1,NINC
      
      EXCESS=0.
      
! --  COMPUTE BASEFLOW AND KEEP TRACK OF TIME INTERVAL SUM.

      BF_lc(K) = LZFPC_lc(K) * DLZP
                
      LZFPC_lc(K) = LZFPC_lc(K) - BF_lc(K)
                
      IF (LZFPC_lc(K) .LE. 0.0001) THEN 
                
      BF_lc(K) = BF_lc(K) + LZFPC_lc(K)
                
      LZFPC_lc(K) = 0.0
                
      ENDIF

      SBF_lc(K) = SBF_lc(K) + BF_lc(K)
                 
      BF_lc(K) = LZFSC_lc(K) * DLZS

      LZFSC_lc(K) = LZFSC_lc(K) - BF_lc(K)
                
      IF (LZFSC_lc(K) .LE. 0.0001) THEN
                
      BF_lc(K) = BF_lc(K) + LZFSC_lc(K)
                
      LZFSC_lc(K) = 0.0
                
      ENDIF         
                 
      SBF_lc(K) = SBF_lc(K) + BF_lc(K)       
                 

! --- COMPUTE PERCOLATION TO LZ IF FREE WATER IS AVAILABLE IN UZ

      IF (PINC + UZFWC_lc(K).GT.0.01) GOTO 460
          
      UZFWC_lc(K) = UZFWC_lc(K) + PINC
          
      GOTO 450

!     --- COMPUTE PERCOLATION DEMAND FROM LZ

460   PERCM_lc(K) = LZFPM(I) * DLZP + LZFSM(I) * DLZS
                
      PERC_lc(K) = PERCM_lc(K) * (UZFWC_lc(K)/UZFWM(I))
                
      DEFR_lc(K)=1.0-((LZTWC_lc(K)+LZFPC_lc(K)+LZFSC_lc(K))/&
     (LZTWM(I)+LZFPM(I)+LZFSM(I)))
                
      PERC_lc(K) = PERC_lc(K) * (1.0 + ZPERC(I) * (DEFR_lc(K)&
     **REXP(I)))

!     --- COMPARE LZ PERCOLATION DEMAND TO UZ FREE WATER AVAILABLE AND COMPUTE ACTUAL PERCOLATION

      IF (PERC_lc(K) .GE. UZFWC_lc(K)) THEN
                
      PERC_lc(K) = UZFWC_lc(K)
                       
      ENDIF
                       
      UZFWC_lc(K) = UZFWC_lc(K) - PERC_lc(K)
                
            
!      --- CHECK TO SEE IF PERC EXCEEDS LZ TENSION AND FREE WATER DEFICIENCY, IF SO SET PERC TO LZ DEFICIENCY


      LZDEF_lc(K) = (LZTWC_lc(K) + LZFPC_lc(K) + LZFSC_lc(K)) - &
     (LZTWM(I) + LZFPM(I) + LZFSM(I)) + PERC_lc(K)
                    
      IF (LZDEF_lc(K) .GT. 0.0) THEN
                    
      PERC_lc(K) = PERC_lc(K) - LZDEF_lc(K)
          
      UZFWC_lc(K) = UZFWC_lc(K) + LZDEF_lc(K)
                       
      ENDIF
                    
      SPERC_lc(K) = SPERC_lc(K) + PERC_lc(K)
!                 SPERC IS THE TIME INTERVAL SUM OF PERC

! --- COMPUTE INTERFLOW FROM UZ

      DEL_lc(K) = UZFWC_lc(K) * DUZ
      
      IF (DEL_lc(K).GT.UZFWC_lc(K))THEN
      
      DEL_lc(K)=UZFWC_lc(K)
      
      UZFWC_lc(K)=0.
      
      ELSE
      
      UZFWC_lc(K) = UZFWC_lc(K) - DEL_lc(K)
      
      ENDIF
                 
      SIF_lc(K) = SIF_lc(K) + DEL_lc(K)

! --- DISRIBUTE PERCOLATED WATER INTO THE LZ STORAGES AND COMPUTE THE REMAINDER IN UZ FREE WATER STORAGE AND RESIDUAL AVAIL FOR RUNOFF
   

!     --- COMPUTE PERC WATER GOING INTO LZ TENSION WATER STORAGE AND COMPARE TO AVAILABLE STORAGE

      PERCT_lc(K) = PERC_lc(K) * (1.0 - PFREE(I))
                
      IF ((PERCT_lc(K) + LZTWC_lc(K)) .GT. LZTWM(I)) THEN
                
!     --- WHEN PERC IS GREATER THAN AVAILABLE TENSION WATER STORAGE, SET TENSION WATER STORAGE TO MAX, REMAINDER OF PERC GETS EVALUATED AGAINST FREE WATER STORAGE

      PERCF_lc(K) = PERCT_lc(K) + LZTWC_lc(K) - LZTWM(I)
                
      LZTWC_lc(K) = LZTWM(I)
                
      ELSE
                
!     --- WHEN PERC IS LESS THAN AVAILABLE TENSION WATER STORAGE, UPDATE TENSION WATER STORAGE

      LZTWC_lc(K) = LZTWC_lc(K) + PERCT_lc(K)
                
      PERCF_lc(K) = 0.0
                
      ENDIF
                
!     --- COMPUTE TOTAL PERC WATER GOING INTO LZ FREE WATER STORAGE

      PERCF_lc(K) = PERCF_lc(K) + PERC_lc(K) * PFREE(I)                

      IF(PERCF_lc(K) .EQ. 0.0) GOTO 470
                
!     --- COMPUTE RELATIVE SIZE OF LZ PRIMARY FREE WATER STORAGE COMPARED TO LZ TOTAL FREE WATER STORAGE

      HPL_lc(K) = LZFPM(I) / (LZFPM(I) + LZFSM(I))
                
!     --- COMPUTE LZ PRIMARY AND SECONDARY FREE WATER CONTENT TO CAPACITY RATIOS

      RATLP_lc(K) = LZFPC_lc(K) / LZFPM(I)
                
      RATLS_lc(K) = LZFSC_lc(K) / LZFSM(I)
                
!     --- COMPUTE FRACTIONS AND PERCENTAGES OF FREE WATER PERC TO GO TO LZ PRIMARY STORAGE

      FRACP_lc(K) = (HPL_lc(K) * 2.0 * (1.0 - RATLP_lc(K))) &
      / ((1.0 - RATLP_lc(K)) + (1.0 - RATLS_lc(K)))
                
      IF (FRACP_lc(K) .GT. 1.0) FRACP_lc(K) = 1.0

      PERCP_lc(K) = PERCF_lc(K) * FRACP_lc(K)
                
      PERCS_lc(K) = PERCF_lc(K) - PERCP_lc(K)
                
!     --- COMPUTE NEW PRIMARY AND SUPPLEMENTAL STORAGE

!         --- COMPUTE NEW SUPPLEMENTAL FREE WATER STORAGE

      LZFSC_lc(K) = LZFSC_lc(K) + PERCS_lc(K)

      IF(LZFSC_lc(K) .GT. LZFSM(I)) THEN
                
!         --- IF NEW SUPPLEMENTAL FREE WATER STORAGE EXCEEDS CAPACITY SET SUPPLEMENTAL STORAGE TO CAPACITY AND EXCESS GOES TO PRIMARY FREE WATER STORAGE

      PERCS_lc(K) = PERCS_lc(K) - LZFSC_lc(K) + LZFSM(I)
                
      LZFSC_lc(K) = LZFSM(I)
                          
      ENDIF
                
            
!         --- IF NEW LZ SUPPLEMENTAL FREE WATER STORAGE IS LESS THAN CAPACITY MOVE ON TO COMPUTE NEW PRIMARY FREE WATER STORAGE


      LZFPC_lc(K) = LZFPC_lc(K) + (PERCF_lc(K) - PERCS_lc(K))

                
      IF (LZFPC_lc(K) .LE. LZFPM(I)) GOTO 470

!             --- IF LZ FREE PRIMARY WATER STORAGE EXCEEDS CAPACITY SET PRIMARY STORAGE TO CAPACITY AND EVALUATE EXCESS AGAINST LZ TENSION WATER STORAGE

      LZFPC_lc(K) = LZFPC_lc(K) + (PERCF_lc(K) - PERCS_lc(K))

                
      IF (LZFPC_lc(K) .LE. LZFPM(I)) GOTO 470

!             --- IF LZ FREE PRIMARY WATER STORAGE EXCEEDS CAPACITY SET PRIMARY STORAGE TO CAPACITY AND EVALUATE EXCESS AGAINST LZ TENSION WATER STORAGE

      LZTWC_lc(K) = LZTWC_lc(K) + LZFPC_lc(K) - LZFPM(I)
                
      LZFPC_lc(K) = LZFPM(I)
      
! -- CHECK TO SEE IF LZTWC>LZTWM, IF SO, SET LZTWC=LZTWM AND SEND EXCESS FOR INCREMENT TO SURFACE RUNOFF

      IF(LZTWC_lc(K).GT.LZTWM(I)) THEN
      
      EXCESS=LZTWC_lc(K)-LZTWM(I)
      
      LZTWC_lc(K)=LZTWM(I)
      
      ELSE
      
      EXCESS=0.
      
      ENDIF
                
!            ---DISTRIBUTE PINC BETWEEN UZFWC AND SURFACE RUNOFF

                       
!                 IF UZFWC CAN ABSORB EXCESS, NO SURFACE RUNOFF

470   IF(PINC+EXCESS.EQ.0.0) GOTO 450

      IF (PINC+UZFWC_lc(K)+EXCESS.GT.UZFWM(I)) GOTO 480
                       
      UZFWC_lc(K)=UZFWC_lc(K)+PINC+EXCESS
                          
      GOTO 450
                          
!                 IF UZFWC CAN NOT ABSORB EXCESS, THE COMPUTE SURFACE RUNOFF
                       
480   SSUR_lc(K)=SSUR_lc(K)+PINC+UZFWC_lc(K)-UZFWM(I)+EXCESS
                          
      UZFWC_lc(K)=UZFWM(I)
450   CONTINUE
   ! if (I==1.and.J==45)       print*,K,BF_lc(K),  SBF_lc(K)
 

! NOTE  the following is based on Law et al

! *****************************************************************************************************
! *****************************************************************************************************
    ! Calculate GEP based on ET and the equation
        GEP_lc(J,M,K) = wue_k(K) * ET_lc(J,M,K) 
        IF (GEP_lc(J,M,K) .LE. 0 .OR. TEMP(I,J,M) .LE. -1.0)   THEN
            GEP_lc(J,M,K)=0.0 
        ENDIF                 
                
        RECO_lc(J,M,K)= reco_inter(K) + reco_slope(K) * GEP_lc(J,M,K)
        
        NEE_lc(J,M,K)=  RECO_lc(J,M,K)- GEP_lc(J,M,K)      
        
! --- COMPUTE FRACTION OF EACH WATER BALANCE COMPONENT AND GEP FOR EACH LAND COVER

! *****************************************************************************************************
! -- CALCULATE THE FRACTION OF Monthly GEP

        GEPTEMP = GEPTEMP + GEP_lc(J,M,K)*LADUSE_lc(I,K)   
        RECOTEMP = RECOTEMP + RECO_lc(J,M,K)*LADUSE_lc(I,K)               
        NEETEMP = NEETEMP + NEE_lc(J,M,K)*LADUSE_lc(I,K)


! *****************************************************************************************************
! --- CALCULATE THE FRACTION OF Monthly AET
               
               AETTEMP = AETTEMP + ET_lc(J,M,K) * LADUSE_lc(I,K)
               
! *****************************************************************************************************
! -- CALCULATE THE FRACTION OF Monthly SURFACE RUNOFF
              
               RUNOFFTEMP = RUNOFFTEMP + SSUR_lc(K) * LADUSE_lc(I,K)
                                            
! *****************************************************************************************************
! -- CALCULATE THE FRACTION OF Monthly SECONDARY BASEFLOW
              
               SBFTEMP = SBFTEMP + SBF_lc(K) * LADUSE_lc(I,K)               

! *****************************************************************************************************
! -- CALCULATE THE FRACTION OF Monthly INTERFLOW
              
               IFTEMP = IFTEMP + SIF_lc(K) * LADUSE_lc(I,K)    
               
! *****************************************************************************************************

! ----CALCUALTE TOTAL RUNOFF FOR EACH LANDUSE, K (GE SUN OCT 19, 2010)             
            ! RUNLAND(I,J_S,M,K) = SURFRO_lc(K) + BF_lc(K) + &
          ! SBF_lc(K) + INF_lc(K)+SIF_lc(K
              
            ! ETLAND(I,J_S,M,K) =ET_lc(J,M,K)
            
            ! GEPLAND(I,J_S,M,K) =GEP_lc(J,M,K)
            
    ! print*,J_S,NYEAR_S+NWARMUP,K,  GEP_lc(J,M,K)
    
! *****************************************************************************************************
! -- ACCUMULATE THE LAND COVER AREA WEIGHTED AVERAGE SOIL MOISTURE FOR THE HUC          
              
              TAUZTWC = TAUZTWC + UZTWC_lc(K) * LADUSE_lc(I,K)
              
              TAUZFWC = TAUZFWC + UZFWC_lc(K) * LADUSE_lc(I,K)
              
              TALZTWC = TALZTWC + LZTWC_lc(K) * LADUSE_lc(I,K)
              
              TALZFPC = TALZFPC + LZFPC_lc(K) * LADUSE_lc(I,K)
              
              TALZFSC = TALZFSC + LZFSC_lc(K) * LADUSE_lc(I,K)
              
              TASM = TASM + (UZTWC_lc(K)+UZFWC_lc(K)+LZTWC_lc(K)+LZFPC_lc(K) &
             +LZFSC_lc(K)) * LADUSE_lc(I,K)
                   
              TAREA = TAREA + LADUSE_lc(I,K) 

! if (I==1.and.J==45) print*,"I=",I,'J=',J,'M=',M,'K=',K,"END LZFSC_lc=",LZFSC_lc(K)
40         CONTINUE ! end of land cover

           
! -- CALCULATE AVG SMC

              AUZTWC = TAUZTWC / TAREA
              AUZFWC = TAUZFWC / TAREA
              ALZTWC = TALZTWC / TAREA
              ALZFPC = TALZFPC / TAREA
              ALZFSC = TALZFSC / TAREA
              ASM = TASM/TAREA

           AET(I,J,M) = AETTEMP        
           RUNOFF(I,J,M) = RUNOFFTEMP
           SECBF(I,J,M) = SBFTEMP
           INTF(I,J,M) = IFTEMP
           SMC(I,J,M) = ASM          
           SP(I,J,M) = SNOWPACK
           AVUZTWC(I,J,M) = AUZTWC
           AVUZFWC(I,J,M) = AUZFWC
           AVLZTWC(I,J,M) = ALZTWC
           AVLZFPC(I,J,M) = ALZFPC
           AVLZFSC(I,J,M) = ALZFSC
                   
           IF (RUNOFF(I,J,M) .LT. 0.) THEN
           
           RUNOFF(I,J,M)=0.
           
           ENDIF            
     

! -- STREAMFLOW IN MILLION M3 FOR EACH HUC FOR MONTH M. HUCAREA IN SQ. METERS

          STRFLOW(I, J, M) = (RUNOFF(I,J,M) + SECBF(I,J,M) + INTF(I,J,M)) &
      * HUCAREA(I)/1000./1000000.
           
           GEPM(I,J, M) = GEPTEMP
           RECOM(I,J,M)  = RECOTEMP
           NEEM(I,J,M) = NEETEMP
!if (I==1 .and. J==45) print*,I,J,M,RUNOFF(I,J,M), SECBF(I,J,M), INTF(I,J,M), HUCAREA(I)
      RETURN
      END
