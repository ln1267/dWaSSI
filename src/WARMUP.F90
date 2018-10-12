!**********************************************************************C
!                                                                      C
!     *** SUBROUTINE RPSDF ***                                         C
!     Read basic input data from GENERAL.TXT and write to              !
!     BASICOUT.TXT, set up column headings in output files             C
!                                                                      !
!**********************************************************************!
SUBROUTINE RPSDF

    USE Common_var
      
    implicit none

    CHARACTER*4 HEADNG(20)
    
    INTEGER ISCENARIOS

! --- Read in data from GENERAL.TXT and write to BASICOUT.TXT

    READ(1,1000) HEADNG
      
1000 FORMAT(20A4)
 
    WRITE(77,2010) HEADNG
      
2010 FORMAT(' ',20A4/)

    READ(1,*) ISCENARIOS

    WRITE(77,*) ISCENARIOS

    READ(1,*) modelscale

    IF (modelscale ==0) THEN
        
        WRITE(77,*) "Model is set as catchment scale simulation using HUC"
        
        WRITE(*,*) "Model is set as catchment scale simulation using HUC"
        
    ELSE
        WRITE(77,*) "Model is set as grid scale simulation using cell information"
        
        WRITE(*,*) "Model is set as grid scale simulation using cell information"
        
    ENDIF

    READ(1,*) NGRID, NYEAR, BYEAR, NLC

    WRITE(77,2020) NGRID
2020 FORMAT(I5,'ACTIVE GRIDS'/)

    WRITE(77,2030) NYEAR, BYEAR

2030 FORMAT(I10,'YEARS TO BE SIMULATED AND',' FIRST YEAR =', I10)

    WRITE(77,2040) NLC

2040 FORMAT('NUMBER OF LAND COVER CATEGORIES: ',I10)

    READ (1,*) LAI_S_Y,LAI_E_Y
    
    WRITE(77,2051) LAI_S_Y,LAI_E_Y
    
2051  FORMAT('FOR LAI input data, the first year',I10, ' , END ',I10)

    READ (1,*) NWARMUP, YSTART, YEND

    IYSTART=YSTART-BYEAR+1

    IYEND=YEND-BYEAR+1

    NYEAR_S=YEND-YSTART+1

    WRITE(77,2050) NYEAR_S,YSTART,IYSTART,YEND, IYEND

    WRITE(*,2050) NYEAR_S,YSTART,IYSTART,YEND, IYEND
    
2050  FORMAT('FOR SIMULATION SUMMARY IN TOTAL',I10,' Years, YEAR TO START',I10,'(ID=',I10, ') , END ',I10,'ID=',I10)

    READ (1,*) FPERD

!--  reduction fraction of Leaf Area index for scenario analysis
    READ (1, *) FPERDLAI

    WRITE(77,2058) FPERD, FPERDLAI
    
2058  FORMAT('DEFOREST RATE%',F10.2, /'LAI REDUCTION%=', F10.2)

    READ (1,*) SNOWPACK

    WRITE(77,1047) SNOWPACK

1047  FORMAT('INTIAL SNOWPACK (MM) = :', F10.2)
    
!-----PRINT TITLE FOR MONTHLY OUTPUT FILE MONTHRUNOFF.TXT
    WRITE (78, 1050)
      
1050  FORMAT('CELL,YEAR,MONTH,PRECIP,TEMP,', &
      'SMC,SNWPK,PET,AET,Sun_ET,',&
      'RUNOFF,BASEFLOW,FLOWMCMMon')
      
!-----PRINT TITLE FOR MONTHLY OUTPUT FILE Soil Storage.TXT
    WRITE (900, 1055)
      
1055  FORMAT('CELL,YEAR,MONTH,UZTWC,UZFWC,',&
      'LZTWC,LZFPC,LZFSC')

!-----PRINT TITLE FOR ANNUAL OUTPUT ANNUALFLOW.TXT
    WRITE (79, 1060)
       
1060  FORMAT ( 'CELL,YEAR,RAIN,PET,',&
      'AET,Sun_ET,RUNOFF,RUN_Pratio,ET_Pratio,RUN_ETRatio,', &
       'SNWPCKMON,RFACTOR')

!-----PRINT TITLE FOR SUMMARY OUTPUT SUMMARRUNOFF.TXT
    WRITE (80, 1070)
      
!-----PRINT TITLE FOR SUMMARY OUTPUT SUMMARRUNOFF.TXT      
1070  FORMAT ('CELL,RAIN,PET,',&
      'AET,RUNOFF,RUNOFF/P,ET/P,(RUN+ET)/P',&
       'RFACTOR,Y_n')
       
!-----PRINT TITLE FOR SUMMARY OUTPUT SUMMARRUNOFF.TXT
    WRITE (910,1080)

1080    FORMAT ('WATERSHEDID,YEAR,LADUSEID,',&
      'HUCRUNOFF,FLOWVOL,LAND%,HUCAREA')

!-----PRINT TITLE FOR HUC FLOW HUCFLOW.TXT
    WRITE (920,1090)

1090    FORMAT ('WATERSHEDID,YEAR,CROPFLOW,',&
       'FORESTFLOW,GRASSFLOW,SHRUBSAVAFLOW,URBANWATERFLOW,TFLOW')

!-----PRINT TITLE FOR SUMMARY HUC CARBON HUCCARBON.TXT
    WRITE (400, 500)
500    FORMAT ('CELL,YEAR,MONTH,GEP(gC/m2/Month),Reco,NEE')

!-----PRINT TITLE FOR ANNUAL CARBON ANNUALCARBON.TXT
    WRITE (500, 600)
600    FORMAT ('CELL,YEAR,GEP(gC/m2/yr),Reco,NEE(gC/m2/yr),',&
              'AET(MM),PET(MM)')
              
!-----PRINT TITLE FOR SUMMARY OUTPUT SUMMARRUNOFF.TXT
    WRITE (600, 650)
    
650    FORMAT ('CELL,NO_YR,GEP(gC/m2/yr),Reco,NEE')

!-----PRINT TITLE FOR SUMMARY OUTPUT SUMMARRUNOFF.TXT
    WRITE (700, 700)

700     FORMAT ('CELL,YEAR,TREE,MAMMALS,BIRD, ', &
       'AMPHIB, REPTILES, VERTEB, AET, PET')
       
!-----PRINT TITLE FOR SUMMARY OUTPUT SUMMARRUNOFF.TXT
    WRITE (800,800)

800      FORMAT ('CELL,NO_YR, TREE, MAMMALS, BIRD,',&
        'AMPHIB, REPTILES, AHUCVERTEB')

      RETURN
      END


!**********************************************************************!
!                                                                      !
!     *** SUBROUTINE RPSWUE ***                                        !
!     Read in WUE paramters for each landuse data from WUE_input.TXT   !
!     landuse based on percent forest decrease (if desired), write     !
!     to BASICOUT.TXT                                                  !
!                                !
!**********************************************************************!
SUBROUTINE RPSWUE

    Use Common_var
    implicit none

    INTEGER I,K

    CHARACTER*1000 DUMY(30),LCname,IGBP

! --- Read and print land use data for each active cell IN THE BASIC.OUT FILE
    WRITE(77,20001)
    
20001  FORMAT(/'WUE paramters INFO FOR EACH SIMULATION CELL'/)

    READ (9,5001) DUMY
    READ (10,5001) DUMY
5001   FORMAT (1000A30)
    
    print*,"Landcover ID,WUE,RECO_inter,RECO_slope,Land cover name"
    
    DO 1011 K=1, NLC

        READ(9,*) I, wue_k(K), reco_inter(K) , reco_slope(K),IGBP,LCname
        READ(10,*) I, ET_Interc(K),P_coef(K),PET_coef(K),LAI_coef(K),&
                   P_PET_coef(K),P_LAI_coef(K),PET_LAI_coef(K),IGBP,LCname
        WRITE(77,11001) I, wue_k(K), reco_inter(K) , reco_slope(K),IGBP,LCname
        
        !print*,I, wue_k(K), reco_inter(K) , reco_slope(K)
    
        WRITE(*,11001) I, wue_k(K), reco_inter(K) , reco_slope(K),IGBP,LCname
        WRITE(*,11002) I, ET_Interc(K),P_coef(K),PET_coef(K),LAI_coef(K),&
                   P_PET_coef(K),P_LAI_coef(K),PET_LAI_coef(K),IGBP,LCname
11001  FORMAT(I5, 3F8.2,',',A30,',',A30)
11002  FORMAT(I5, 7F8.2,',',A30,',',A30)
1011    CONTINUE

    RETURN
END


!**********************************************************************!
!                                                                      !
!     *** SUBROUTINE RPSINT ***                                        !
!     Read in landuse data from CELLINFO.TXT AND SOIL PARAMETERS FROM SOILINFO.TXT!
!                                !
!**********************************************************************!
SUBROUTINE RPSINT

    Use Common_var
    
    implicit none

    INTEGER I,ID,K

    CHARACTER*1000 DUMY(30)

! --- Read and print land use data for each active cell IN THE BASIC.OUT FILE
    WRITE(77,2000)
      
2000 FORMAT(/'LANDUSE INFO FOR EACH SIMULATION CELL'/)

    READ (2,500) DUMY
      
500    FORMAT (1000A30)


! ----LANC = raw Landcover types

    WRITE (77,500) DUMY

    DO 10 I=1, NGRID
    
        IF (modelscale==0) THEN

            READ(2,*) ID, HUCNO(I),HUCAREA(I), LATUDE(I), LONGI(I),(LADUSE_lc(I,K),K=1, NLC)

        ELSE

            READ(2,*) ID, HUCNO(I), LATUDE(I), LONGI(I),LADUSE(I) !,HUCELE(I)

        ENDIF

10    CONTINUE

! --- Read and print SOIL PARAMETERS for each active cell IN THE BASIC.OUT FILE
!
      WRITE(77,2051)
      
2051  FORMAT(/'SOIL PARAMETERS FOR EACH SIMULATION CELL'/)

      READ (7,550) DUMY
      
550   FORMAT (30A8)

      WRITE (77,550) DUMY

      DO 15 I=1, NGRID

      READ(7,*) ID, HUCNO(I), UZTWM(I), UZFWM(I), UZK(I), ZPERC(I),&
     REXP(I), LZTWM(I), LZFSM(I), LZFPM(I), LZSK(I),&
     LZPK(I), PFREE(I)


15    CONTINUE

      RETURN
      END



!C**********************************************************************C
!C                                                                      C
!C     *** SUBROUTINE RPSLAI ***                                        C
!C     Input MONTHLY LAI, FILL IN GAPS FOR PERIODS WITH NO LAI DATA     C
!C                                                                      C
!C**********************************************************************C
      SUBROUTINE RPSLAI

      use Common_var
      
      implicit none
      
      INTEGER(kind=4) I
      
      INTEGER(kind=4) YEAR, J, M,Mon,K

      INTEGER(kind=4) IY_LAI_END,IY_LAI_START

      CHARACTER*100 TEMPHEAD3 (11)

    Print*, "For LAI, Start Year=",LAI_S_Y,"END Year=",LAI_E_Y,NEW_LINE('A')

!   SET THE YEAR ID OF INPUT LAI DATA-----
    IF (BYEAR .LT. LAI_S_Y ) then
           
        IY_LAI_START=LAI_S_Y-BYEAR+1
           
    ELSE
        
        IY_LAI_START=1
         
    ENDIF
       
    If (YEND .GT. LAI_E_Y) then
      
        IY_LAI_END=LAI_E_Y-BYEAR+1
        
    Else
       
        IY_LAI_END=LAI_E_Y-BYEAR+1
        
    Endif

    print*,"reading LAI",NEW_LINE('A')

! --- READ IN LAI DATA FROM LANDLAI.TXT

    IF (modelscale ==0) THEN ! read multiple landcover LAI

        DO 2011 I=1, NGRID

            DO 3011 J= IY_LAI_START,IY_LAI_END

                DO 4011 M=1, 12

                    ! this is for reading the header
                    IF (I .EQ. 1 .AND. J .EQ. IY_LAI_START .AND. M .EQ. 1) THEN

                       READ (8, 9021) TEMPHEAD3

         9021           FORMAT (100A11)

                    ENDIF


! --- read all landuse LAIs

                    READ(8,*) HUCNO(I),YEAR,Mon,(LAI_lc(I,J,M,K), K=1, NLC)

4011             CONTINUE

3011          CONTINUE

2011       CONTINUE

! --- ASSIGN YEAR LAI_S_Y LAI DATA TO YEARS BEFORE LAI_S_Y
        IF  ( BYEAR .LT. LAI_S_Y)  then
          DO 2021 I=1, NGRID

             DO 3021 J=1, IY_LAI_START-1

                DO 4021 M=1, 12

                    DO 5021 K=1, NLC

                        LAI_lc(I,J,M,K) = LAI_lc(I,IY_LAI_START,M,K)

5021                CONTINUE

4021             CONTINUE

3021          CONTINUE

2021        CONTINUE
!
        ENDIF

!C--- ASSIGN YEAR IY_LAI_END LAI DATA TO YEARS AFTER IY_LAI_END
        IF (IYEND .GT. IY_LAI_END) then
          DO 2031 I=1, NGRID

             DO 3031 J=IY_LAI_END+1, NYEAR

                DO 4031 M=1, 12

                    DO 5031 K=1, NLC

                        LAI_lc(I,J,M,K) = LAI_lc(I,IY_LAI_END,M,K)

5031                CONTINUE

4031             CONTINUE

3031          CONTINUE

2031        CONTINUE

        ENDIF

    ELSE

        DO 201 I=1, NGRID

            DO 301 J= IY_LAI_START,IY_LAI_END

                DO 401 M=1, 12

                    ! this is for reading the header
                    IF (I .EQ. 1 .AND. J .EQ. IY_LAI_START .AND. M .EQ. 1) THEN

                        READ (8, 902) TEMPHEAD3

902                       FORMAT (100A11)

                    ENDIF


! --- Read monthly data for each grid
                    READ(8,*) HUCNO(I),YEAR,Mon,LAI(I,J,M)

401             CONTINUE

301          CONTINUE

201       CONTINUE

    print*,"finished reading LAI",NEW_LINE('A')

! --- ASSIGN YEAR LAI_S_Y LAI DATA TO YEARS BEFORE LAI_S_Y
        IF  ( BYEAR .LT. LAI_S_Y)  then
          DO 202 I=1, NGRID

             DO 302 J=1, LAI_S_Y-1

                DO 402 M=1, 12

                LAI(I,J,M) = LAI(I,LAI_S_Y,M)

402             CONTINUE

302          CONTINUE

202        CONTINUE
!
        ENDIF
!
!C--- ASSIGN YEAR IY_LAI_END LAI DATA TO YEARS AFTER IY_LAI_END
      IF (IYEND .GT. IY_LAI_END) then
          DO 203 I=1, NGRID

             DO 303 J=IY_LAI_END+1, NYEAR

                DO 403 M=1, 12

                LAI(I,J,M) = LAI(I,IY_LAI_END,M)

403             CONTINUE

303          CONTINUE

203        CONTINUE

      ENDIF

    ENDIF
    RETURN
END

!C**********************************************************************C
!C                                                                      C
!C     *** SUBROUTINE RPSCLIMATE ***                                    C
!C     Input MONTHLY CLIMATE DATA, CALCULATE ANNUAL PPT                 C
!C                                                                      C
!C**********************************************************************C
      SUBROUTINE RPSCLIMATE

      USE Common_var
      implicit none

      INTEGER YEAR

      INTEGER I, J, M,Mon

      REAL,POINTER :: ANNPPT(:,:),SUMANPPT(:)

      CHARACTER*10 TEMPHEAD (10)

      ALLOCATE (ANNPPT(MAX_GRIDS,MAX_YEARS))
      ALLOCATE (SUMANPPT(MAX_GRIDS))
      
      ANNPPT =0.

      SUMANPPT = 0.

      AAPPT = 0.


      DO 5000 I=1,NGRID

         DO 5001 J=1,NYEAR

            DO 5002 M=1,12

               IF (I .EQ. 1 .AND. J .EQ. 1 .AND. M .EQ. 1) THEN

                  READ (4, 900) TEMPHEAD

                 WRITE (77, 900) TEMPHEAD

               ENDIF

900            FORMAT (10A10)

               READ(4,*) HUCNO(I), YEAR, Mon, RAIN(I,J,M), TEMP(I,J,M)

               ANNPPT(I, J) = ANNPPT(I, J) + RAIN(I,J,M)

5002        CONTINUE

            SUMANPPT(I) = SUMANPPT(I) + ANNPPT(I, J)

5001     CONTINUE

         AAPPT(I) = SUMANPPT(I)/NYEAR

5000  CONTINUE

      RETURN
      END
