!  ==============================================================================================================================\n
!  MODULE 	: readdim2
! 
!  CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
!  LICENCE      : IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        This module contains subroutines for reading the forcing file for the dim2_driver.
!!
!!
!!\n DESCRIPTION  : This module contains subroutines for reading the forcing file for the dim2_driver.
!!                  Following subroutines are public and called from dim2_driver :
!!                    - forcing_info : Open the forcing file and return information about the grid in the forcing file. 
!!                                     Prepare for a zoom if needed. 
!!                                     Initialization of parallelization related to the grid. 
!!                    - forcing_read : Return the forcing data for the current time step of the model. The forcing file will
!!                                     be read if it has not already been done for the current time-step in the forcing file.
!!                    - forcing_grid : Calculate the longitudes and latitudes of the model grid.
!! 
!! RECENT CHANGE(S): None 
!! 
!! REFERENCE(S) : None
!!   
!! SVN     :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_driver/readdim2.f90 $ 
!! $Date: 2017-10-26 15:35:04 +0200 (Thu, 26 Oct 2017) $
!! $Revision: 4719 $
!! \n
!_ ================================================================================================================================

MODULE readdim2

   USE ioipsl_para
   USE weather
   USE TIMER
   USE constantes
   USE time
   USE solar
   USE grid
   USE mod_orchidee_para

   IMPLICIT NONE

   PRIVATE
   PUBLIC  :: forcing_read, forcing_info, forcing_grid

   INTEGER, SAVE                            :: iim_full, jjm_full, llm_full, ttm_full
   INTEGER, SAVE                            :: iim_zoom, jjm_zoom
   INTEGER, SAVE                            :: iim_g_begin,jjm_g_begin,iim_g_end,jjm_g_end
   REAL, SAVE, ALLOCATABLE, DIMENSION(:,:)  :: data_full, lon_full, lat_full
   REAL, SAVE, ALLOCATABLE, DIMENSION(:)    :: lev_full
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: itau, i_index, j_index,j_index_g
   INTEGER, SAVE                            :: i_test, j_test
   INTEGER, SAVE                            :: printlev_loc   !! Local printlev
   LOGICAL, SAVE                            :: allow_weathergen, interpol, daily_interpol
   LOGICAL, SAVE, PUBLIC                    :: weathergen, is_watchout
   REAL, SAVE                               :: merid_res, zonal_res
   LOGICAL, SAVE                            :: have_zaxis=.FALSE.
!-
!- Heigh controls and data 
!- 
   LOGICAL, SAVE                            :: zfixed, zsigma, zhybrid, zlevels, zheight 
   LOGICAL, SAVE                            :: zsamelev_uv 
   REAL, SAVE                               :: zlev_fixed, zlevuv_fixed 
   REAL, SAVE                               :: zhybrid_a, zhybrid_b 
   REAL, SAVE                               :: zhybriduv_a, zhybriduv_b 

CONTAINS

!! ==============================================================================================================================\n
!! SUBROUTINE 	: forcing_info
!!
!>\BRIEF        Open the forcing file and return information about the grid in the forcing file. 
!!
!!\n DESCRIPTION : This subroutine will get all the information from the forcing file and prepare for the zoom if needed.
!!                 It returns to the caller the sizes of the data it will receive at the forcing_read call. 
!!                 This is important so that the caller can allocate the right space.
!!
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S)	: 
!!
!_ ================================================================================================================================
  SUBROUTINE forcing_info(filename, iim, jjm, llm, tm, date0, dt_force, force_id)

    IMPLICIT NONE

!! 0.1 Input variables
    CHARACTER(LEN=*), INTENT(in) :: filename     !! Name of the file to be opened

!! 0.2 Output variables
    INTEGER, INTENT(out)         :: iim          !! Size in x of the forcing data
    INTEGER, INTENT(out)         :: jjm          !! Size in y of the forcing data
    INTEGER, INTENT(out)         :: llm          !! Number of levels in the forcing data (not yet used)
    INTEGER, INTENT(out)         :: tm           !! Time dimension of the forcing
    REAL, INTENT(out)            :: date0        !! The date at which the forcing file starts (julian days)
    REAL, INTENT(out)            :: dt_force     !! Time-step of the forcing file in seconds
    INTEGER, INTENT(out)         :: force_id     !! Id of the forcing file

!! 0.3 Local variables
    CHARACTER(LEN=20)                  :: calendar_str
    CHARACTER(LEN=200)                 :: printstr       !! temporary character string to contain error message 
    REAL                               :: juld_1, juld_2
    REAL, ALLOCATABLE, DIMENSION(:,:)  :: fcontfrac
    REAL, ALLOCATABLE, DIMENSION(:,:)  :: qair
    LOGICAL                            :: contfrac_exists=.FALSE.
    INTEGER                            :: NbPoint
    INTEGER                            :: i_test,j_test
    INTEGER                            :: i,j,ind,ttm_part
    INTEGER, ALLOCATABLE, DIMENSION(:) :: index_l
    REAL, ALLOCATABLE, DIMENSION(:,:)  :: lon, lat
    REAL, ALLOCATABLE, DIMENSION(:)    :: lev, levuv

    !-
    CALL flininfo(filename,  iim_full, jjm_full, llm_full, ttm_full, force_id)
    !-
    IF ( printlev_loc>=3 ) WRITE(numout,*) 'forcing_info : Details from forcing file :', &
         iim_full, jjm_full, llm_full, ttm_full
    !-
    IF ( llm_full < 1 ) THEN
       have_zaxis = .FALSE.
    ELSE
       have_zaxis = .TRUE.
    ENDIF
    IF ( printlev_loc>=3 ) WRITE(numout,*) 'have_zaxis : ', llm_full, have_zaxis
    !-
    ttm_part = 2
    ALLOCATE(itau(ttm_part))
    ALLOCATE(data_full(iim_full, jjm_full),lon_full(iim_full, jjm_full),&
         & lat_full(iim_full, jjm_full))
    ALLOCATE(lev_full(llm_full))
    ALLOCATE(fcontfrac(iim_full,jjm_full))
    !-
    lev_full(:) = zero
    !-
    dt_force=zero
    CALL flinopen &
         &  (filename, .FALSE., iim_full, jjm_full, llm_full, lon_full, lat_full, &
         &   lev_full, ttm_part, itau, date0, dt_force, force_id)
    IF ( dt_force == zero ) THEN
       dt_force = itau(2) - itau(1) 
       itau(:) = itau(:) / dt_force
    ENDIF
    !  WRITE(numout,*) "forcing_info : Forcing time step out of flinopen : ",dt_force
    !-
    !- What are the alowed options for the temportal interpolation
    !-
    !Config Key   = ALLOW_WEATHERGEN
    !Config Desc  = Allow weather generator to create data
    !Config If    = [-]
    !Config Def   = n
    !Config Help  = This flag allows the forcing-reader to generate
    !Config         synthetic data if the data in the file is too sparse
    !Config         and the temporal resolution would not be enough to
    !Config         run the model.
    !Config Units = [FLAG]
    !-
    allow_weathergen = .FALSE.
    CALL getin_p('ALLOW_WEATHERGEN',allow_weathergen)
    !-
    !- The calendar was set by the forcing file. If no "calendar" attribute was
    !- found then it is assumed to be gregorian, 
    !MM => FALSE !! it is NOT assumed anything !
    !- else it is what ever is written in this attribute.
    !-
    CALL ioget_calendar(calendar_str)
    i=INDEX(calendar_str,ACHAR(0))
    IF ( i > 0 ) THEN
       calendar_str(i:20)=' '
    ENDIF
    !  WRITE(numout,*) "forcing_info : Calendar used : ",calendar_str
    IF ( calendar_str == 'XXXX' ) THEN
       IF (printlev_loc >= 1) WRITE(numout,*) "forcing_info : The calendar was not found in the forcing file."
       IF (allow_weathergen) THEN
          ! Then change the calendar
          CALL ioconf_calendar("noleap") 
       ELSE
          IF ( printlev_loc>=1 ) WRITE(numout,*) "forcing_info : We will force it to gregorian by default."
          CALL ioconf_calendar("gregorian") !! = 365.2425 ; "noleap" = 365.0; "360d"; "julian"=365.25
       ENDIF
    ENDIF
    IF (printlev_loc >= 1) WRITE(numout,*) "readdim2 : Calendar used by the model : ",calendar_str
    IF (ttm_full .GE. 2) THEN
       juld_1 = itau2date(itau(1), date0, dt_force)
       juld_2 = itau2date(itau(2), date0, dt_force)
    ELSE
       juld_1 = 0
       juld_2 = 0
       CALL ipslerr_p ( 3, 'forcing_info','What is that only one time step in the forcing file ?', &
            &         ' That can not be right.','verify forcing file.')
    ENDIF
    !-
    !- Initialize one_year / one_day
    CALL ioget_calendar (one_year, one_day)
    !-
    !- What is the distance between the two first states. From this we will deduce what is
    !- to be done.
    weathergen = .FALSE.
    interpol = .FALSE.
    daily_interpol = .FALSE.
    is_watchout = .FALSE.
    !-
    IF ( ABS(ABS(juld_2-juld_1)-30.) .LE. 2.) THEN
       IF ( allow_weathergen ) THEN
          weathergen = .TRUE.
          IF (printlev_loc >= 1) WRITE(numout,*) 'Using weather generator.' 
       ELSE
          CALL ipslerr_p ( 3, 'forcing_info', &
               &         'This seems to be a monthly file.', &
               &         'We should use a weather generator with this file.', &
               &         'This should be allowed in the run.def')
       ENDIF
    ELSEIF (( ABS(juld_1-juld_2) .LE. 1./4.) .OR. ( ABS(juld_1-juld_2) .EQ. 1.)) THEN
       interpol = .TRUE.
       IF (printlev_loc >= 1) WRITE(numout,*) 'We will interpolate between the forcing data time steps.' 
       IF ( ABS(juld_1-juld_2) .EQ. 1.) THEN
          daily_interpol = .TRUE.
       ENDIF
    ELSE
       ! Using the weather generator with data other than monthly ones probably
       ! needs some thinking.
       WRITE(numout,*) 'The time step is not suitable:',ABS(juld_1-juld_2),' days.'
       CALL ipslerr_p ( 3, 'forcing_info','The time step is not suitable.', &
            &         '','We cannot do anything with these forcing data.')
    ENDIF
    !-
    !- redefine the forcing time step if the weather generator is activated
    !-
    IF ( weathergen ) THEN
       !Config Key   = DT_WEATHGEN
       !Config Desc  = Calling frequency of weather generator
       !Config If    = ALLOW_WEATHERGEN
       !Config Def   = 1800.
       !Config Help  = Determines how often the weather generator
       !Config         is called (time step in s). Should be equal
       !Config         to or larger than Sechiba's time step (say,
       !Config         up to 6 times Sechiba's time step or so).
       !Config Units = [seconds]
       dt_force = 1800.
       CALL getin_p('DT_WEATHGEN',dt_force)
    ENDIF
    !-
    !- Define the zoom
    !-
    !Config Key   = LIMIT_WEST
    !Config Desc  = Western limit of region
    !Config If    = [-]
    !Config Def   = -180.
    !Config Help  = Western limit of the region we are 
    !Config         interested in. Between -180 and +180 degrees
    !Config         The model will use the smalest regions from
    !Config         region specified here and the one of the forcing file.
    !Config Units = [Degrees] 
    !- 
    limit_west = -180.
    CALL getin_p('LIMIT_WEST',limit_west)
    !-
    !Config Key   = LIMIT_EAST
    !Config Desc  = Eastern limit of region
    !Config If    = [-]
    !Config Def   = 180.
    !Config Help  = Eastern limit of the region we are
    !Config         interested in. Between -180 and +180 degrees
    !Config         The model will use the smalest regions from
    !Config         region specified here and the one of the forcing file.
    !Config Units = [Degrees] 
    !-
    limit_east = 180.
    CALL getin_p('LIMIT_EAST',limit_east)
    !-
    !Config Key   = LIMIT_NORTH
    !Config Desc  = Northern limit of region
    !Config If    = [-]
    !Config Def   = 90.
    !Config Help  = Northern limit of the region we are
    !Config         interested in. Between +90 and -90 degrees
    !Config         The model will use the smalest regions from
    !Config         region specified here and the one of the forcing file.
    !Config Units = [Degrees]
    !-
    limit_north = 90.
    CALL getin_p('LIMIT_NORTH',limit_north)
    !-
    !Config Key   = LIMIT_SOUTH
    !Config Desc  = Southern limit of region
    !Config If    = [-]
    !Config Def   = -90.
    !Config Help  = Southern limit of the region we are
    !Config         interested in. Between 90 and -90 degrees
    !Config         The model will use the smalest regions from
    !Config         region specified here and the one of the forcing file.
    !Config Units = [Degrees]
    !-
    limit_south = -90.
    CALL getin_p('LIMIT_SOUTH',limit_south)
    !- 
    !- Calculate domain size
    !-
    IF ( interpol ) THEN
       !-
       !- If we use temporal interpolation, then we cannot change the resolution (yet?)
       !-
       ALLOCATE(i_index(iim_full), j_index(jjm_full),j_index_g(jjm_full))
       IF (is_root_prc) THEN

          CALL domain_size (limit_west, limit_east, limit_north, limit_south,&
               &         iim_full, jjm_full, lon_full, lat_full, iim_zoom, jjm_zoom,&
               &         i_index, j_index_g)

          j_index(:)=j_index_g(:)

          ALLOCATE(qair(iim_full,jjm_full))
          CALL flinget_buffer (force_id,'Qair',iim_full, jjm_full, 1, ttm_full,  1, 1, data_full)
          CALL forcing_zoom(data_full, qair)

          CALL flinquery_var(force_id, 'contfrac', contfrac_exists)
          IF ( contfrac_exists ) THEN
             IF (printlev_loc >= 1) WRITE(numout,*) "contfrac exist in the forcing file."
             CALL flinget_buffer (force_id,'contfrac',iim_full, jjm_full, 1, ttm_full,  1, 1, data_full)
             CALL forcing_zoom(data_full, fcontfrac)
             IF (printlev_loc >= 2) WRITE(numout,*) "fcontfrac min/max :", &
                  MINVAL(fcontfrac(1:iim_zoom,1:jjm_zoom)),MAXVAL(fcontfrac(1:iim_zoom,1:jjm_zoom))
          ELSE
             fcontfrac(:,:)=1.
          ENDIF


          DO i=1,iim_zoom
             DO j=1,jjm_zoom
                IF ( fcontfrac(i,j) <= EPSILON(1.) ) THEN
                   qair(i,j) = 999999.
                ENDIF
             ENDDO
          ENDDO

          ALLOCATE(index_l(iim_zoom*jjm_zoom))
          !- In this point is returning the global NbPoint with the global index
          CALL forcing_landind(iim_zoom,jjm_zoom,qair,NbPoint,index_l,i_test,j_test)
          ! 
          ! Work out the vertical layers to be used 
          ! 
          CALL forcing_vertical_ioipsl(force_id) 
       ELSE
          ALLOCATE(index_l(1))
       ENDIF

       ! Initiate global grid and parallelism
       CALL bcast(iim_zoom)
       CALL bcast(jjm_zoom)
       CALL bcast(NbPoint)
       CALL grid_set_glo(iim_zoom,jjm_zoom,NbPoint)
       CALL grid_allocate_glo(4)
       
       ! Check consistency in the number of procs and the land points selected
       ! in order to prevent an exception 
       IF (NbPoint < mpi_size) THEN
           WRITE(printstr,*) 'The number of landpoints found (', NbPoint,') is less than the number of processors selected (', mpi_size,')'
           CALL ipslerr_p(3, 'forcing_info', 'Wrong parallelization options', &
                          TRIM(printstr), &
                          'Increase the window (EAST_BOUND, ...) or decrease the number of processors')
       ENDIF

       !
       !- global index index_g is the index_l of root proc 
       IF (is_root_prc) index_g(:)=index_l(1:NbPoint)

       DEALLOCATE(index_l)

       ! 
       ! Distribute to all processors the information on the forcing 
       ! 
       CALL bcast(index_g)
       CALL Init_orchidee_data_para_driver(nbp_glo,index_g)
       CALL init_ioipsl_para

       ! Initialize printlev_loc
       printlev_loc=get_printlev('readdim2')
       IF (printlev_loc >= 2) WRITE(numout,*) 'Standard PRINTLEV= ', printlev
       IF (printlev_loc >= 2) WRITE(numout,*) 'Local PRINTLEV_readdim2= ', printlev_loc

!     CALL Init_writeField_p

       CALL bcast(jjm_zoom)
       CALL bcast(i_index)
       CALL bcast(j_index_g)
       CALL bcast(zfixed) 
       CALL bcast(zsigma) 
       CALL bcast(zhybrid) 
       CALL bcast(zlevels) 
       CALL bcast(zheight) 
       CALL bcast(zsamelev_uv) 
       CALL bcast(zlev_fixed) 
       CALL bcast(zlevuv_fixed) 
       CALL bcast(zhybrid_a) 
       CALL bcast(zhybrid_b) 
       CALL bcast(zhybriduv_a)  
       CALL bcast(zhybriduv_b) 
       ind=0
       DO j=1,jjm_zoom
          IF ( (j >= jj_begin) .AND. (j <= jj_end) ) THEN
             ind=ind+1
             j_index(ind)=j_index_g(j)
          ENDIF
       ENDDO

       jjm_zoom=jj_nb
       iim_zoom=iim_g

       !-
       !- If we use the weather generator, then we read zonal and meridional resolutions
       !- This should be unified one day...
       !-
    ELSEIF ( weathergen ) THEN
       !-
       !Config Key   = MERID_RES
       !Config Desc  = North-South Resolution
       !Config Def   = 2.
       !Config If    = ALLOW_WEATHERGEN
       !Config Help  = North-South Resolution of the region we are
       !Config         interested in. 
       !Config Units = [Degrees]
       !-
       merid_res = 2.
       CALL getin_p('MERID_RES',merid_res)
       !-
       !Config Key   = ZONAL_RES
       !Config Desc  = East-West Resolution
       !Config Def   = 2.
       !Config If    = ALLOW_WEATHERGEN
       !Config Help  = East-West Resolution of the region we are
       !Config         interested in. In degrees
       !Config Units = [Degrees] 
       !-
       zonal_res = 2.
       CALL getin_p('ZONAL_RES',zonal_res)
       !-
       !- Number of time steps is meaningless in this case
       !-
       !    ttm_full = HUGE( ttm_full )
       !MM Number (realistic) of time steps for half hour dt
       ttm_full = NINT(one_year * 86400. / dt_force)
       !-
       IF (is_root_prc) THEN

          !- In this point is returning the global NbPoint with the global index
          CALL weathgen_domain_size (limit_west,limit_east,limit_north,limit_south, &
               zonal_res,merid_res,iim_zoom,jjm_zoom)
          ALLOCATE(index_l(iim_zoom*jjm_zoom))
       ENDIF
       CALL bcast(iim_zoom)
       CALL bcast(jjm_zoom)

       ALLOCATE(lon(iim_zoom,jjm_zoom))
       ALLOCATE(lat(iim_zoom,jjm_zoom))
       ALLOCATE(lev(llm_full),levuv(llm_full))
       
       ! We need lon and lat now for weathgen_init
       CALL forcing_grid (iim_zoom,jjm_zoom,llm_full,lon,lat,init_f=.TRUE.)
       CALL forcing_vertical_ioipsl(-1) 

       IF (is_root_prc) THEN
          CALL weathgen_init &
               &        (filename,dt_force,force_id,iim_zoom,jjm_zoom, &
               &         zonal_res,merid_res,lon,lat,index_l,NbPoint)
!!$,&
!!$               &         i_index,j_index_g)
       ELSE
          ALLOCATE(index_l(1))
          index_l(1)=1
       ENDIF

       CALL bcast(NbPoint)
       CALL grid_set_glo(iim_zoom,jjm_zoom,NbPoint)
       CALL grid_allocate_glo(4)

       !
       !- global index index_g is the index_l of root proc 
       IF (is_root_prc) index_g(:)=index_l(1:NbPoint)

       DEALLOCATE(index_l)

     CALL bcast(index_g)
     CALL Init_orchidee_data_para_driver(nbp_glo,index_g)
     CALL init_ioipsl_para
!     CALL Init_writeField_p
     CALL bcast(jjm_zoom)
!!$       CALL bcast(i_index)
!!$       CALL bcast(j_index_g)

!!$       ind=0
!!$       DO j=1,jjm_zoom
!!$          IF ( (j >= jj_begin) .AND. (j <= jj_end) ) THEN
!!$             ind=ind+1
!!$             j_index(ind)=j_index_g(j)
!!$          ENDIF
!!$       ENDDO

       jjm_zoom=jj_nb
       iim_zoom=iim_g
       !
       CALL weathgen_read_file(force_id,iim_zoom,jjm_zoom)

       !-
    ELSE
       !-
       CALL ipslerr_p(3,'forcing_info','Neither interpolation nor weather generator is specified.','','')
       !-
    ENDIF
    !-
    !- Transfer the right information to the caller
    !-
    iim = iim_zoom
    jjm = jjm_zoom
    llm = 1
    tm = ttm_full
    !-
    IF ( printlev_loc>=3 ) WRITE(numout,*) 'forcing_info : end : ', iim,jjm, llm,tm
    !-
  END SUBROUTINE forcing_info


!! ==============================================================================================================================\n
!! SUBROUTINE 	: forcing_read
!!
!>\BRIEF        Return forcing data for the current time step
!!
!!\n DESCRIPTION : Return the forcing data for the current time step of the model. The forcing file will
!!                 be read if it has not already been done for the current time-step in the forcing file. 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S)	: 
!!
!_ ================================================================================================================================
SUBROUTINE forcing_read &
  & (filename, rest_id, lrstread, lrstwrite, &
  &  itauin, istp, itau_split, split, nb_spread, lwdown_cons, swdown_cons, date0,   &
  &  dt_force, iim, jjm, lon, lat, zlev, zlevuv, ttm, &
  &  swdown, coszang, precip, snowf, tair, u, v, qair, pb, lwdown, &
  &  fcontfrac, fneighbours, fresolution, &
  &  SWnet, Eair, petAcoef, peqAcoef, petBcoef, peqBcoef, cdrag, ccanopy, &
  &  kindex, nbindex, force_id)

  IMPLICIT NONE

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
   CHARACTER(LEN=*), INTENT(IN) :: filename   !! name of the file to be opened
   INTEGER, INTENT(IN) :: force_id            !! FLINCOM file id. It is used to close the file at the end of the run.
   INTEGER, INTENT(IN) :: rest_id             !! ID of restart file
   LOGICAL, INTENT(IN) :: lrstread            !! read restart file?
   LOGICAL, INTENT(IN) :: lrstwrite           !! write restart file?
   INTEGER, INTENT(IN) :: itauin              !! time step for which we need the data
   INTEGER, INTENT(IN) :: istp                !! time step for restart file
   INTEGER, INTENT(IN) :: itau_split          !! Current step between 2 forcing times-step (it decides if it is time to read)
   INTEGER, INTENT(IN) :: split               !! The number of time steps between two time-steps of the forcing
   INTEGER, INTENT(IN) :: nb_spread           !! Over how many time steps do we spread the precipitation
   LOGICAL, INTENT(IN) :: lwdown_cons         !! Flag to conserve lwdown radiation from forcing
   LOGICAL, INTENT(IN) :: swdown_cons         !! Flag to conserve swdown radiation from forcing
   REAL, INTENT(IN)    :: date0               !! The date at which the forcing file starts (julian days)
   REAL, INTENT(IN)    :: dt_force            !! time-step of the forcing file in seconds
   INTEGER, INTENT(IN) :: iim                 !! Size of the grid in x
   INTEGER, INTENT(IN) :: jjm                 !! Size of the grid in y
   INTEGER, INTENT(IN) :: ttm                 !! number of time steps in all in the forcing file
   REAL,DIMENSION(iim,jjm), INTENT(IN) :: lon !! Longitudes
   REAL,DIMENSION(iim,jjm), INTENT(IN) :: lat !! Latitudes

    !! 0.2 Output variables
   REAL,DIMENSION(iim,jjm), INTENT(OUT) :: zlev        !! First Levels if it exists (ie if watchout file)
   REAL,DIMENSION(iim,jjm), INTENT(OUT) :: zlevuv      !! First Levels of the wind (equal precedent, if it exists)
   REAL,DIMENSION(iim,jjm), INTENT(OUT) :: swdown      !! Downward solar radiation (W/m^2)
   REAL,DIMENSION(iim,jjm), INTENT(OUT) :: coszang     !! Cosine of the solar zenith angle (unitless)
   REAL,DIMENSION(iim,jjm), INTENT(OUT) :: precip      !! Precipitation (Rainfall) (kg/m^2s) 
   REAL,DIMENSION(iim,jjm), INTENT(OUT) :: snowf       !! Snowfall (kg/m^2s)
   REAL,DIMENSION(iim,jjm), INTENT(OUT) :: tair        !! 1st level (2m ? in off-line) air temperature (K)
   REAL,DIMENSION(iim,jjm), INTENT(OUT) :: u           !! 1st level (2m/10m ? in off-line) (in theory !) wind speed (m/s)
   REAL,DIMENSION(iim,jjm), INTENT(OUT) :: v           !! 1st level (2m/10m ? in off-line) (in theory !) wind speed (m/s)
   REAL,DIMENSION(iim,jjm), INTENT(OUT) :: qair        !! 1st level (2m ? in off-line) humidity (kg/kg)
   REAL,DIMENSION(iim,jjm), INTENT(OUT) :: pb          !! Surface pressure (Pa)
   REAL,DIMENSION(iim,jjm), INTENT(OUT) :: lwdown      !! Downward long wave radiation (W/m^2)
   REAL,DIMENSION(iim,jjm), INTENT(OUT) :: fcontfrac   !! Continental fraction (no unit)
   REAL,DIMENSION(iim,jjm,2), INTENT(OUT) :: fresolution     !! resolution in x and y dimensions for each points
   INTEGER,DIMENSION(iim,jjm,8), INTENT(OUT) :: fneighbours  !! land neighbours

   !! From a WATCHOUT file :
   REAL,DIMENSION(iim,jjm), INTENT(OUT) :: SWnet       !! Net surface short-wave flux
   REAL,DIMENSION(iim,jjm), INTENT(OUT) :: Eair        !! Air potential energy
   REAL,DIMENSION(iim,jjm), INTENT(OUT) :: petAcoef    !! Coeficients A from the PBL resolution for T
   REAL,DIMENSION(iim,jjm), INTENT(OUT) :: peqAcoef    !! Coeficients A from the PBL resolution for q
   REAL,DIMENSION(iim,jjm), INTENT(OUT) :: petBcoef    !! Coeficients B from the PBL resolution for T
   REAL,DIMENSION(iim,jjm), INTENT(OUT) :: peqBcoef    !! Coeficients B from the PBL resolution for q
   REAL,DIMENSION(iim,jjm), INTENT(OUT) :: cdrag       !! Surface drag
   REAL,DIMENSION(iim,jjm), INTENT(OUT) :: ccanopy     !! CO2 concentration in the canopy

    !! 0.3 Modified variable
   INTEGER, INTENT(INOUT) :: nbindex                   !! Number of land points
   INTEGER,DIMENSION(iim*jjm), INTENT(INOUT) :: kindex !! Index of all land-points in the data (used for the gathering)

    !! 0.4 Local variables
   INTEGER :: ik,i,j

   IF ( interpol ) THEN

     CALL forcing_read_interpol &
         (filename, itauin, itau_split, split, nb_spread, lwdown_cons, swdown_cons, date0,   &
          dt_force, iim, jjm, lon, lat, zlev, zlevuv, ttm, &
          swdown, coszang, precip, snowf, tair, u, v, qair, pb, lwdown, &
          fcontfrac, fneighbours, fresolution, &
          SWnet, Eair, petAcoef, peqAcoef, petBcoef, peqBcoef, cdrag, ccanopy, &
          kindex, nbindex, force_id)

   ELSEIF ( weathergen ) THEN

      IF (lrstread) THEN
         fcontfrac(:,:) = 1.0
         IF (printlev_loc >= 2) WRITE(numout,*) 'contfrac : ', MINVAL(fcontfrac), MAXVAL(fcontfrac)
      ENDIF

      IF ( (itauin == 0).AND.(itau_split == 0) ) THEN
         CALL weathgen_main (istp, istp, filename, force_id, iim, jjm, 1, &
              rest_id, lrstread, lrstwrite, &
              limit_west, limit_east, limit_north, limit_south, &
              zonal_res, merid_res, lon, lat, date0, dt_force, &
              kindex, nbindex, &
              swdown, precip, snowf, tair, u, v, qair, pb, lwdown)
      ELSE
         CALL weathgen_main (itauin, istp, filename, force_id, iim, jjm, 1, &
              rest_id, lrstread, lrstwrite, &
              limit_west, limit_east, limit_north, limit_south, &
              zonal_res, merid_res, lon, lat, date0, dt_force, &
              kindex, nbindex, &
              swdown, precip, snowf, tair, u, v, qair, pb, lwdown)
      ENDIF

      IF ( (itauin == 0).AND.(itau_split == 0) ) THEN
         !---
         !--- Allocate grid stuff
         !---
         CALL grid_init ( nbindex, 4, "RegLonLat", "ForcingGrid" )
         !---
         !--- Compute
         !---
         CALL grid_stuff(nbp_glo, iim_g, jjm_g, lon_g, lat_g, index_g)
         !CALL grid_stuff (nbindex, iim, jjm, lon, lat, kindex)
         DO ik=1,nbindex
         
            j = ((kindex(ik)-1)/iim) + 1
            i = (kindex(ik) - (j-1)*iim)
            !-
            !- Store variable to help describe the grid
            !- once the points are gathered.
            !-
            fneighbours(i,j,:) = neighbours(ik,:)
            !
            fresolution(i,j,:) = resolution(ik,:)
         ENDDO
      ENDIF
   ELSE

      CALL ipslerr_p(3,'forcing_read','Neither interpolation nor weather generator is specified.','','')

   ENDIF

   IF (.NOT. is_watchout) THEN
      ! We have to compute Potential air energy
      WHERE(tair(:,:) < val_exp) 
         eair(:,:) = cp_air*tair(:,:)+cte_grav*zlev(:,:)
      ENDWHERE
   ENDIF

END SUBROUTINE forcing_read

!! ==============================================================================================================================\n
!! SUBROUTINE 	: forcing_read_interpol
!!
!>\BRIEF        
!!
!!\n DESCRIPTION : 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S)	: 
!!
!_ ================================================================================================================================
SUBROUTINE forcing_read_interpol &
  & (filename, itauin, itau_split, split, nb_spread, lwdown_cons, swdown_cons, date0,   &
  &  dt_force, iim, jjm, lon, lat, zlev, zlevuv, ttm, swdown, coszang, rainf, snowf, tair, &
  &  u, v, qair, pb, lwdown, &
  &  fcontfrac, fneighbours, fresolution, &
  &  SWnet, Eair, petAcoef, peqAcoef, petBcoef, peqBcoef, cdrag, ccanopy, &
  &  kindex, nbindex, force_id)
!---------------------------------------------------------------------
!- filename   : name of the file to be opened
!- itauin     : time step for which we need the data
!- itau_split : Where are we within the splitting
!-              of the time-steps of the forcing files
!-              (it decides IF we READ or not)
!- split      : The number of time steps we do
!-              between two time-steps of the forcing
!- nb_spread  : Over how many time steps do we spread the precipitation
!- lwdown_cons: flag that decides if lwdown radiation should be conserved.
!- swdown_cons: flag that decides if swdown radiation should be conserved.
!- date0      : The date at which the forcing file starts (julian days)
!- dt_force   : time-step of the forcing file in seconds
!- iim        : Size of the grid in x
!- jjm        : size of the grid in y
!- lon        : Longitudes
!- lat        : Latitudes
!- zlev       : First Levels if it exists (ie if watchout file)
!- zlevuv     : First Levels of the wind (equal precedent, if it exists)
!- ttm        : number of time steps in all in the forcing file
!- swdown     : Downward solar radiation (W/m^2)
!- coszang    : Cosine of the solar zenith angle (unitless)
!- rainf      : Rainfall (kg/m^2s)
!- snowf      : Snowfall (kg/m^2s)
!- tair       : 2m air temperature (K)
!- u and v    : 2m (in theory !) wind speed (m/s)
!- qair       : 2m humidity (kg/kg)
!- pb         : Surface pressure (Pa)
!- lwdown     : Downward long wave radiation (W/m^2)
!- fcontfrac  : Continental fraction (no unit)
!- fneighbours: land neighbours
!- fresolution: resolution in x and y dimensions for each points
!-
!- From a WATCHOUT file :
!- SWnet      : Net surface short-wave flux
!- Eair       : Air potential energy
!- petAcoef   : Coeficients A from the PBL resolution for T
!- peqAcoef   : Coeficients A from the PBL resolution for q
!- petBcoef   : Coeficients B from the PBL resolution for T
!- peqBcoef   : Coeficients B from the PBL resolution for q
!- cdrag      : Surface drag
!- ccanopy    : CO2 concentration in the canopy
!- 
!- kindex     : Index of all land-points in the data
!-              (used for the gathering)
!- nbindex    : Number of land points
!- force_id   : FLINCOM file id.
!-              It is used to close the file at the end of the run.
!---------------------------------------------------------------------
   IMPLICIT NONE
!-
   INTEGER,PARAMETER :: lm=1
!-
!- Input variables
!-
   CHARACTER(LEN=*) :: filename
   INTEGER :: itauin, itau_split, split, nb_spread
   LOGICAL, INTENT(IN) :: lwdown_cons, swdown_cons
   REAL    :: date0, dt_force
   INTEGER :: iim, jjm, ttm
   REAL,DIMENSION(:,:),INTENT(IN) :: lon, lat   !- LOCAL data array (size=iim,jjm)
   INTEGER, INTENT(IN) :: force_id
!-
!- Output variables
!-
   REAL,DIMENSION(:,:),INTENT(OUT) :: zlev, zlevuv, &  !- LOCAL data array (size=iim,jjm)
  &  swdown, coszang, rainf, snowf, tair, u, v, qair, pb, lwdown, &
  &  fcontfrac
   REAL,DIMENSION(:,:,:),INTENT(OUT) :: fresolution    !- LOCAL data array (size=iim,jjm,2)
   INTEGER,DIMENSION(:,:,:),INTENT(OUT) :: fneighbours !- LOCAL data array (size=iim,jjm,8)
   ! for watchout files
   REAL,DIMENSION(:,:),INTENT(OUT) :: &
  &  SWnet, Eair, petAcoef, peqAcoef, petBcoef, peqBcoef, cdrag, ccanopy
   INTEGER,DIMENSION(:),INTENT(INOUT) :: kindex    !- LOCAL index of the map
   INTEGER, INTENT(INOUT) :: nbindex
!-
!- Local variables
!-
   INTEGER, SAVE :: last_read=0
   INTEGER, SAVE :: itau_read, itau_read_nm1=0, itau_read_n=0
   REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: &
  &  zlev_nm1, zlevuv_nm1, swdown_nm1, rainf_nm1, snowf_nm1, tair_nm1, u_nm1, v_nm1, qair_nm1, & 
  &  pb_nm1, lwdown_nm1, & 
  &  zlev_n, zlevuv_n, swdown_n, rainf_n, snowf_n, tair_n, u_n, v_n, qair_n, &
  &  pb_n, lwdown_n, mean_coszang

   REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: &
  &  startday_n, startday_nm1, daylength_n, daylength_nm1, tmax_n, tmax_nm1, tmin_nm1, tmin_nm2, tmin_n,   &
  &  qsatta, qsattmin_n, qsattmin_nm1, qmin_n, qmin_nm1, qmax_n, qmax_nm1, qsa
   REAL,SAVE    :: hour

!  just for grid stuff if the forcing file is a watchout file
   REAL, ALLOCATABLE, DIMENSION(:,:) :: tmpdata
   ! variables to be read in watchout files
   REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: &
  &  SWnet_nm1, Eair_nm1, petAcoef_nm1, peqAcoef_nm1, petBcoef_nm1, peqBcoef_nm1, cdrag_nm1, ccanopy_nm1, &
  &  SWnet_n, Eair_n, petAcoef_n, peqAcoef_n, petBcoef_n, peqBcoef_n, cdrag_n, ccanopy_n
   REAL, SAVE :: julian_for ! Date of the forcing to be read
   REAL :: julian, ss, rw
!jur, 
   REAL, SAVE :: julian0    ! First day of this year
   INTEGER :: yy, mm, dd, is, i, j, ik
   REAL(r_std), DIMENSION(2) :: min_resol, max_resol
!  if Wind_N and Wind_E are in the file (and not just Wind)
   LOGICAL, SAVE :: wind_N_exists=.FALSE.
   LOGICAL       :: wind_E_exists=.FALSE.
   LOGICAL, SAVE :: contfrac_exists=.FALSE.
   LOGICAL, SAVE :: neighbours_exists=.FALSE.
   LOGICAL, SAVE :: initialized = .FALSE.
!  to bypass FPE problem on integer convertion with missing_value heigher than precision
   INTEGER, PARAMETER                         :: undef_int = 999999999
!---------------------------------------------------------------------

   itau_read = MOD((itauin-1),ttm)+1

   IF (printlev_loc >= 5) THEN
      WRITE(numout,*) &
           " FORCING READ : itauin, itau_read, itau_split : ",&
           itauin, itau_read, itau_split
   ENDIF

!-
!- This part initializes the reading of the forcing. As you can see
!- we only go through here if both time steps are zero.
!- 
   IF ( (itau_read == 0).AND.(itau_split == 0) ) THEN
!-
!- Tests on forcing file type
     CALL flinquery_var(force_id, 'Wind_N', wind_N_exists)
     IF ( wind_N_exists ) THEN
        CALL flinquery_var(force_id, 'Wind_E', wind_E_exists)
        IF ( .NOT. wind_E_exists ) THEN
           CALL ipslerr_p(3,'forcing_read_interpol', &
   &             'Variable Wind_E does not exist in forcing file', &
   &             'But variable Wind_N exists.','Please, rename Wind_N to Wind;') 
        ENDIF
     ENDIF
     CALL flinquery_var(force_id, 'levels', is_watchout)
     IF ( is_watchout ) THEN
        WRITE(numout,*) "Read a Watchout File."
     ENDIF
     CALL flinquery_var(force_id, 'contfrac', contfrac_exists)
!-
     IF (printlev_loc >= 5) WRITE(numout,*) 'ALLOCATE all the memory needed'
!-
     ALLOCATE &
    &  (swdown_nm1(iim,jjm), rainf_nm1(iim,jjm), snowf_nm1(iim,jjm), &
    &   tair_nm1(iim,jjm), u_nm1(iim,jjm), v_nm1(iim,jjm), qair_nm1(iim,jjm), &
    &   pb_nm1(iim,jjm), lwdown_nm1(iim,jjm))
     ALLOCATE &
    &  (swdown_n(iim,jjm), rainf_n(iim,jjm), snowf_n(iim,jjm), &
    &   tair_n(iim,jjm), u_n(iim,jjm), v_n(iim,jjm), qair_n(iim,jjm), &
    &   pb_n(iim,jjm), lwdown_n(iim,jjm))

     IF(daily_interpol) THEN
        ALLOCATE &
       &  (startday_n(iim,jjm), startday_nm1(iim,jjm), daylength_n(iim,jjm), &
       &   daylength_nm1(iim,jjm), tmax_n(iim,jjm), tmax_nm1(iim,jjm), tmin_n(iim,jjm), &
       &   tmin_nm1(iim,jjm), tmin_nm2(iim,jjm), qsatta(iim,jjm), qsattmin_n(iim,jjm), qsattmin_nm1(iim,jjm), &
       &   qmin_n(iim,jjm), qmin_nm1(iim,jjm), qmax_n(iim,jjm), qmax_nm1(iim,jjm), qsa(iim,jjm) )
     ENDIF


     ALLOCATE &
    &  (zlev_nm1(iim,jjm), zlev_n(iim,jjm), zlevuv_nm1(iim,jjm), zlevuv_n(iim,jjm), &
    &   SWnet_nm1(iim,jjm), Eair_nm1(iim,jjm), cdrag_nm1(iim,jjm), ccanopy_nm1(iim,jjm), &
    &   petAcoef_nm1(iim,jjm), peqAcoef_nm1(iim,jjm), petBcoef_nm1(iim,jjm), peqBcoef_nm1(iim,jjm), &
    &   SWnet_n(iim,jjm), Eair_n(iim,jjm), cdrag_n(iim,jjm), ccanopy_n(iim,jjm), &
    &   petAcoef_n(iim,jjm), peqAcoef_n(iim,jjm), petBcoef_n(iim,jjm), peqBcoef_n(iim,jjm))
     ALLOCATE &
    &  (mean_coszang(iim,jjm))
!-
     IF (printlev_loc >= 5) WRITE(numout,*) 'Memory ALLOCATED'
!-
!- We need for the driver surface air temperature and humidity before the
!- the loop starts. Thus we read it here.
!-     
     CALL forcing_just_read (iim, jjm, zlev, zlevuv, ttm, 1, 1, &
          &  swdown, rainf, snowf, tair, &
          &  u, v, qair, pb, lwdown, &
          &  SWnet, Eair, petAcoef, peqAcoef, petBcoef, peqBcoef, cdrag, ccanopy, &
          &  force_id, wind_N_exists)
!----

!-- First in case it's not a watchout file
     IF ( .NOT. is_watchout ) THEN
        IF ( contfrac_exists ) THEN
           IF (printlev_loc >= 1) WRITE(numout,*) "contfrac exist in the forcing file."
           CALL flinget_buffer (force_id,'contfrac',iim_full, jjm_full, llm_full, ttm,  1, 1, data_full)
           CALL forcing_zoom(data_full, fcontfrac)
           IF (printlev_loc >= 2) WRITE(numout,*) "fcontfrac min/max :", &
                MINVAL(fcontfrac(1:iim_zoom,jjm_zoom)),MAXVAL(fcontfrac(1:iim_zoom,jjm_zoom))
           !
           ! We need to make sure that when we gather the points we pick all
           ! the points where contfrac is above 0. Thus we prepare tair for
           ! subroutine forcing_landind
           !
           DO i=1,iim
              DO j=1,jjm
                 IF ( j==1 .AND. i<ii_begin) fcontfrac(i,j)=0.     ! bande de recouvrement du scatter2D 
                 IF ( j==jjm .AND. i>ii_end) fcontfrac(i,j)=0.     ! => on mets les donn�es qu'on veut pas du noeud � missing_value
                 IF ( fcontfrac(i,j) <= EPSILON(1.) ) THEN
                    tair(i,j) = 999999.
                 ENDIF
              ENDDO
           ENDDO
        ELSE
           fcontfrac(:,:) = 1.0
        ENDIF
        !---
        !--- Create the index table
        !---
        !--- This job return a LOCAL kindex
        CALL forcing_landind(iim, jjm, tair, nbindex, kindex, i_test, j_test)
#ifdef CPP_PARA
        ! We keep previous function forcing_landind, only to get a valid (i_test,j_test) 
        ! Force nbindex points for parallel computation
        nbindex=nbp_loc
        CALL scatter(index_g,kindex(1:nbindex))
        kindex(1:nbindex)=kindex(1:nbindex)-(jj_begin-1)*iim_g
#endif
        ik=MAX(nbindex/2,1)
        j_test = (((kindex(ik)-1)/iim) + 1)
        i_test = (kindex(ik) - (j_test-1)*iim)
        IF (printlev_loc >= 5) THEN
           WRITE(numout,*) 'New test point is : ', i_test, j_test
        ENDIF
        !---
        !--- Allocate grid stuff
        !---
        CALL grid_init ( nbindex, 4, "RegLonLat", "ForcingGrid" )
        !---
        !--- All grid variables
        !---
        CALL grid_stuff(nbp_glo, iim_g, jjm_g, lon_g, lat_g, index_g)
        DO ik=1,nbindex
           !
           j = ((kindex(ik)-1)/iim) + 1
           i = (kindex(ik) - (j-1)*iim)
           !-
           !- Store variable to help describe the grid
           !- once the points are gathered.
              !-
           fneighbours(i,j,:) = neighbours(ik,:)
           !
           fresolution(i,j,:) = resolution(ik,:)
        ENDDO
     ELSE
!-- Second, in case it is a watchout file
        ALLOCATE (tmpdata(iim,jjm))
        tmpdata(:,:) = 0.0
!--
        IF ( .NOT. contfrac_exists ) THEN
           CALL ipslerr_p (3,'forcing_read_interpol', &
 &          'Could get contfrac variable in a watchout file :',TRIM(filename), &
 &          '(Problem with file ?)')
        ENDIF
        CALL flinget_buffer (force_id,'contfrac',iim_full, jjm_full, llm_full, ttm,  1, 1, data_full)
        CALL forcing_zoom(data_full, fcontfrac)
        !
        ! We need to make sure that when we gather the points we pick all
        ! the points where contfrac is above 0. Thus we prepare tair for
        ! subroutine forcing_landind
        !
        DO i=1,iim
           DO j=1,jjm
	      IF ( j==1 .AND. i<ii_begin) fcontfrac(i,j)=0.
	      IF ( j==jjm .AND. i>ii_end) fcontfrac(i,j)=0.
              IF ( fcontfrac(i,j) <= EPSILON(1.) ) THEN
                 tair(i,j) = 999999.
              ENDIF
           ENDDO
        ENDDO
        !---
        !--- Create the index table
        !---
        !--- This job return a LOCAL kindex
        CALL forcing_landind(iim, jjm, tair, nbindex, kindex, i_test, j_test)
#ifdef CPP_PARA
        ! We keep previous function forcing_landind, only to get a valid (i_test,j_test) 
        ! Force nbindex points for parallel computation
        nbindex=nbp_loc
        CALL scatter(index_g,kindex)
        kindex(:)=kindex(:)-offset
!        kindex(:)=kindex(:)-(jj_begin-1)*iim_g
#endif
        ik=MAX(nbindex/2,1)
        j_test = (((kindex(ik)-1)/iim) + 1)
        i_test = (kindex(ik) - (j_test-1)*iim)
        IF (printlev_loc >= 5) THEN
           WRITE(numout,*) 'New test point is : ', i_test, j_test
        ENDIF
        !---
        !--- Allocate grid stuff
        !---
        CALL grid_init ( nbindex, 4, "RegLonLat", "ForcingGrid" )
        neighbours(:,:) = -1
        resolution(:,:) = 0.
        min_resol(:) = 1.e6
        max_resol(:) = -1.
        !---
        !--- All grid variables
        !---
        !-
        !- Get variables to help describe the grid
        CALL flinquery_var(force_id, 'neighboursNN', neighbours_exists)
        IF ( .NOT. neighbours_exists ) THEN
           CALL ipslerr_p (3,'forcing_read_interpol', &
 &          'Could get neighbours in a watchout file :',TRIM(filename), &
 &          '(Problem with file ?)')
        ELSE
           IF (printlev_loc >= 2) WRITE(numout,*) "Watchout file contains neighbours and resolutions."
        ENDIF
        !
        fneighbours(:,:,:) = undef_int
        !
        !- once the points are gathered.
        CALL flinget_buffer (force_id,'neighboursNN',iim_full, jjm_full, llm_full, ttm,  1, 1, data_full)
        CALL forcing_zoom(data_full, tmpdata)
        WHERE(tmpdata(:,:) < undef_int)
           fneighbours(:,:,1) = NINT(tmpdata(:,:))
        ENDWHERE
        !
        CALL flinget_buffer (force_id,'neighboursNE',iim_full, jjm_full, llm_full, ttm,  1, 1, data_full)
        CALL forcing_zoom(data_full, tmpdata)
        WHERE(tmpdata(:,:) < undef_int)
           fneighbours(:,:,2) = NINT(tmpdata(:,:))
        ENDWHERE
        !
        CALL flinget_buffer (force_id,'neighboursEE',iim_full, jjm_full, llm_full, ttm,  1, 1, data_full)
        CALL forcing_zoom(data_full, tmpdata)
        WHERE(tmpdata(:,:) < undef_int)
           fneighbours(:,:,3) = NINT(tmpdata(:,:))
        ENDWHERE
        !
        CALL flinget_buffer (force_id,'neighboursSE',iim_full, jjm_full, llm_full, ttm,  1, 1, data_full)
        CALL forcing_zoom(data_full, tmpdata)
        WHERE(tmpdata(:,:) < undef_int)
           fneighbours(:,:,4) = NINT(tmpdata(:,:))
        ENDWHERE
        !
        CALL flinget_buffer (force_id,'neighboursSS',iim_full, jjm_full, llm_full, ttm,  1, 1, data_full)
        CALL forcing_zoom(data_full, tmpdata)
        WHERE(tmpdata(:,:) < undef_int)
           fneighbours(:,:,5) = NINT(tmpdata(:,:))
        ENDWHERE
        !
        CALL flinget_buffer (force_id,'neighboursSW',iim_full, jjm_full, llm_full, ttm,  1, 1, data_full)
        CALL forcing_zoom(data_full, tmpdata)
        WHERE(tmpdata(:,:) < undef_int)
           fneighbours(:,:,6) = NINT(tmpdata(:,:))
        ENDWHERE
        !
        CALL flinget_buffer (force_id,'neighboursWW',iim_full, jjm_full, llm_full, ttm,  1, 1, data_full)
        CALL forcing_zoom(data_full, tmpdata)
        WHERE(tmpdata(:,:) < undef_int)
           fneighbours(:,:,7) = NINT(tmpdata(:,:))
        ENDWHERE
        !
        CALL flinget_buffer (force_id,'neighboursNW',iim_full, jjm_full, llm_full, ttm,  1, 1, data_full)
        CALL forcing_zoom(data_full, tmpdata)
        WHERE(tmpdata(:,:) < undef_int)
           fneighbours(:,:,8) = NINT(tmpdata(:,:))
        ENDWHERE
        !
        ! now, resolution of the grid
        CALL flinget_buffer (force_id,'resolutionX',iim_full, jjm_full, llm_full, ttm,  1, 1, data_full)
        CALL forcing_zoom(data_full, tmpdata)
        fresolution(:,:,1) = tmpdata(:,:)
        !
        CALL flinget_buffer (force_id,'resolutionY',iim_full, jjm_full, llm_full, ttm,  1, 1, data_full)
        CALL forcing_zoom(data_full, tmpdata)
        fresolution(:,:,2) = tmpdata(:,:)
        !
        DO ik=1,nbindex
           !
           j = ((kindex(ik)-1)/iim) + 1
           i = (kindex(ik) - (j-1)*iim)
           !-
           !- Store variable to help describe the grid
           !- once the points are gathered.
           !-
           neighbours(ik,:) = fneighbours(i,j,:)  
           !
           resolution(ik,:) = fresolution(i,j,:)
           !
        
        ENDDO
        CALL gather(neighbours,neighbours_g)
        CALL gather(resolution,resolution_g)
        min_resol(1) = MINVAL(resolution(:,1))
        min_resol(2) = MAXVAL(resolution(:,2))
        max_resol(1) = MAXVAL(resolution(:,1))
        max_resol(2) = MAXVAL(resolution(:,2))
        !
        area(:) = resolution(:,1)*resolution(:,2)
        CALL gather(area,area_g)
!--
        DEALLOCATE (tmpdata)
     ENDIF
     IF (printlev_loc >= 2) WRITE(numout,*) 'contfrac : ', MINVAL(fcontfrac), MAXVAL(fcontfrac)
!---
   ENDIF
!---
   IF (printlev_loc >= 5) THEN
      WRITE(numout,*) &
           & 'The dates : ',itau_read,itau_split,itau_read_nm1,itau_read_n
   ENDIF
!---
!--- Here we do the work in case only interpolation is needed.
!---
   IF ( initialized .AND. interpol ) THEN
!---
      IF ( daily_interpol ) THEN

         IF (split > 1) THEN
            IF ( itau_split <= (split/2.) ) THEN
               rw = REAL(itau_split+split/2.)/split
            ELSE 
               rw = REAL(itau_split-split/2.)/split
            ENDIF
         ELSE
            rw = 1.
         ENDIF

         IF ((last_read == 0) .OR. ( rw==(1./split)) ) THEN
   !---
   !-----   Start or Restart
            IF (last_read == 0) THEN
               ! Case of a restart or a shift in the forcing file.
               IF (itau_read > 1) THEN
                  itau_read_nm1=itau_read-1
                  CALL forcing_just_read (iim, jjm, zlev_nm1, zlevuv_nm1, ttm, itau_read_nm1, itau_read_nm1, &
                       &  swdown_nm1, rainf_nm1, snowf_nm1, tmin_nm1, &
                       &  u_nm1, v_nm1, qair_nm1, pb_nm1, lwdown_nm1, &
                       &  SWnet_nm1, Eair_nm1, petAcoef_nm1, peqAcoef_nm1, petBcoef_nm1, peqBcoef_nm1, cdrag_nm1, ccanopy_nm1, &
                       &  force_id, wind_N_exists)
                  CALL forcing_just_read_tmax (iim, jjm, ttm, itau_read_nm1, itau_read_nm1, tmax_nm1, force_id )
               ! Case of a simple start.
               ELSE 
                  itau_read_nm1 = un
                  IF (printlev_loc >= 2) WRITE(numout,*) "we will use the forcing of the first day to initialize "
                  CALL forcing_just_read (iim, jjm, zlev_nm1, zlevuv_nm1, ttm, itau_read_nm1, itau_read_nm1, &
                       &  swdown_nm1, rainf_nm1, snowf_nm1, tmin_nm1, &
                       &  u_nm1, v_nm1, qair_nm1, pb_nm1, lwdown_nm1, &
                       &  SWnet_nm1, Eair_nm1, petAcoef_nm1, peqAcoef_nm1, petBcoef_nm1, peqBcoef_nm1, cdrag_nm1, ccanopy_nm1, &
                       &  force_id, wind_N_exists)
                  CALL forcing_just_read_tmax (iim, jjm, ttm, itau_read_nm1, itau_read_nm1, tmax_nm1, force_id )
               ENDIF
               tmin_nm2(:,:)=tmin_nm1(:,:)
               IF ( dt_force .GT. 3600. ) THEN
                  mean_coszang(:,:) = 0.0
                  daylength_n(:,:) = 0.
                  DO is=1,split
                     !MM we compute mean SWdown between t and t+Dt then I take t+Dt/2. 
                     julian = julian_for+((is-0.5)/split)*dt_force/one_day
      !!$               julian = julian_for+(FLOAT(is)/split)*dt_force/one_day
                     CALL solarang (julian, julian0, iim, jjm, lon*0, lat, coszang)
                     mean_coszang(:,:) = mean_coszang(:,:)+coszang(:,:)
                     WHERE( coszang(:,:) > 0. ) 
                        daylength_n(:,:)=daylength_n(:,:)+1./split*24
                     ENDWHERE
                  ENDDO
                  mean_coszang(:,:) = mean_coszang(:,:)/split
                  daylength_nm1(:,:)=daylength_n(:,:)
      !            WRITE(*,*) "mean_coszang =",MAXVAL(mean_coszang)
               ENDIF
            ELSE
   !-----   Normal mode : copy old step
               swdown_nm1(:,:) = swdown_n(:,:)
               rainf_nm1(:,:) = rainf_n(:,:)
               snowf_nm1(:,:)  = snowf_n(:,:) 
               tair_nm1(:,:)   = tair_n(:,:)
               u_nm1(:,:)      = u_n(:,:)
               v_nm1(:,:)      = v_n(:,:)
               qair_nm1(:,:)   = qair_n(:,:)
               pb_nm1(:,:)     = pb_n(:,:)
               lwdown_nm1(:,:) = lwdown_n(:,:)
               tmin_nm2(:,:)   = tmin_nm1(:,:)
               tmin_nm1(:,:)   = tmin_n(:,:)
               tmax_nm1(:,:)   = tmax_n(:,:)

               IF (is_watchout) THEN
                  zlev_nm1(:,:)   = zlev_n(:,:)
                  zlevuv_nm1(:,:) = zlevuv_n(:,:)
                  ! Net surface short-wave flux
                  SWnet_nm1(:,:) = SWnet_n(:,:)
                  ! Air potential energy
                  Eair_nm1(:,:)   = Eair_n(:,:)
                  ! Coeficients A from the PBL resolution for T
                  petAcoef_nm1(:,:) = petAcoef_n(:,:)
                  ! Coeficients A from the PBL resolution for q
                  peqAcoef_nm1(:,:) = peqAcoef_n(:,:)
                  ! Coeficients B from the PBL resolution for T
                  petBcoef_nm1(:,:) = petBcoef_n(:,:)
                  ! Coeficients B from the PBL resolution for q
                  peqBcoef_nm1(:,:) = peqBcoef_n(:,:)
                  ! Surface drag
                  cdrag_nm1(:,:) = cdrag_n(:,:)
                  ! CO2 concentration in the canopy
                  ccanopy_nm1(:,:) = ccanopy_n(:,:)
               ENDIF
               itau_read_nm1 = itau_read_n
            ENDIF
   !-----
   !-----
            IF(last_read==0)THEN
               itau_read_n = itau_read
            ELSE
               itau_read_n = itau_read+1
            ENDIF

            IF (itau_read_n > ttm) THEN
               WRITE(numout,*) 'WARNING --WARNING --WARNING --WARNING '
               WRITE(numout,*) &
                    &  'WARNING : We are going back to the start of the file'
               itau_read_n =1
            ENDIF
            IF (printlev_loc >= 5) THEN
               WRITE(numout,*) &
                    & 'The dates 2 : ',itau_read,itau_split,itau_read_nm1,itau_read_n
            ENDIF
   !-----
   !----- Get a reduced julian day !
   !----- This is needed because we lack the precision on 32 bit machines.
   !-----
            IF ( dt_force .GT. 3600. ) THEN
               julian_for = itau2date(itau_read-1, date0, dt_force)
               CALL ju2ymds (julian_for, yy, mm, dd, ss)
   
               ! first day of this year
               CALL ymds2ju (yy,1,1,0.0, julian0)
   !-----
               IF (printlev_loc >= 5) THEN
                  WRITE(numout,*) 'Forcing for Julian day ',julian_for,'is read'
                  WRITE(numout,*) 'Date for this day ',yy,' / ',mm,' / ',dd,"  ",ss
               ENDIF
            ENDIF
   !-----
            CALL forcing_just_read (iim, jjm, zlev_n, zlevuv_n, ttm, itau_read_n, itau_read_n, &
                 &  swdown_n, rainf_n, snowf_n, tmin_n, &
                 &  u_n, v_n, qair_n, pb_n, lwdown_n, &
                 &  SWnet_n, Eair_n, petAcoef_n, peqAcoef_n, petBcoef_n, peqBcoef_n, cdrag_n, ccanopy_n, &
                 &  force_id, wind_N_exists)
            CALL forcing_just_read_tmax (iim, jjm, ttm, itau_read_n, itau_read_n, tmax_n, force_id )

   !---
            last_read = itau_read_n
   !-----
   !----- Compute mean solar angle for the comming period
   !-----
            IF (printlev_loc >= 5) WRITE(numout,*) 'Going into  solarang', split, one_day
   !-----

   !-----
         ENDIF
   !---
         IF ( itau_split == 1. ) THEN
            IF ( dt_force .GT. 3600. ) THEN
               mean_coszang(:,:) = 0.0
               daylength_nm1(:,:)=daylength_n(:,:)
               daylength_n(:,:) = 0.
               DO is=1,split
                  !MM we compute mean SWdown between t and t+Dt then I take t+Dt/2. 
                  julian = julian_for+((is-0.5)/split)*dt_force/one_day
   !!$               julian = julian_for+(FLOAT(is)/split)*dt_force/one_day
                  CALL solarang (julian, julian0, iim, jjm, lon*0, lat, coszang)
                  mean_coszang(:,:) = mean_coszang(:,:)+coszang(:,:)
                  WHERE( coszang(:,:) > 0. ) 
                     daylength_n(:,:)=daylength_n(:,:)+1./split*24
                  ENDWHERE
               ENDDO
               mean_coszang(:,:) = mean_coszang(:,:)/split
   !            WRITE(*,*) "mean_coszang =",MAXVAL(mean_coszang)
            ENDIF
         ENDIF
 
   !--- Do the interpolation
         IF (printlev_loc >= 5) WRITE(numout,*) 'Doing the interpolation between time steps'
   !---

         IF (printlev_loc >= 5) WRITE(numout,*) 'Coeff of interpollation : ',rw
   !---

         pb(:,:) = (pb_n(:,:)-pb_nm1(:,:))*rw + pb_nm1(:,:)
         u(:,:)  = (u_n(:,:)-u_nm1(:,:))*rw + u_nm1(:,:)
         v(:,:)  = (v_n(:,:)-v_nm1(:,:))*rw + v_nm1(:,:)

   !--- Take care of the height of the vertical levels 
         zlev(:,:) = (zlev_n(:,:)-zlev_nm1(:,:))*rw + zlev_nm1(:,:) 
         zlevuv(:,:) = (zlevuv_n(:,:)-zlevuv_nm1(:,:))*rw + zlevuv_nm1(:,:) 

         hour=REAL(itau_split)/split*24
         startday_n(:,:)=12.-daylength_n(:,:)/2.
         startday_nm1(:,:)=12.-daylength_nm1(:,:)/2.

         WHERE ( ( hour >= startday_n(:,:) ) .AND. ( hour > 12) .AND. ( hour <= 14) )
            tair(:,:)=(tmax_nm1(:,:)-tmin_nm1(:,:))/2 * ( sin(pi/(14-startday_n(:,:))*(hour-0.5* &
           &    (14.-startday_n(:,:))-startday_n(:,:))) )+ (tmax_nm1(:,:)+tmin_nm1(:,:))/2.
         ELSEWHERE( ( hour >= startday_n(:,:) ) .AND. ( hour <= 12) )
            tair(:,:)=(tmax_n(:,:)-tmin_n(:,:))/2 * ( sin(pi/(14-startday_n(:,:))*(hour-0.5* &
           &    (14.-startday_n(:,:))-startday_n(:,:))) )+ (tmax_n(:,:)+tmin_n(:,:))/2.
         ELSEWHERE ( hour < startday_n(:,:) )
            tair(:,:)=(tmax_nm1(:,:)-tmin_n(:,:))/2.*sin(pi/(24.-14.+startday_nm1(:,:) )* &
           &    (hour + 24.+0.5*(24.-14.+startday_nm1(:,:) )-14.))+(tmax_nm1(:,:)+tmin_n(:,:))/2.
         ELSEWHERE
            tair(:,:)=(tmax_nm1(:,:)-tmin_n(:,:))/2.*sin(pi/(24.-14.+startday_n(:,:))*(hour+0.5* &
           &    (24.-14.+startday_n(:,:))-14.))+(tmax_nm1(:,:)+tmin_n(:,:))/2.
         ENDWHERE

         CALL weathgen_qsat_2d (iim,jjm,tmin_n,pb,qsattmin_n)
         CALL weathgen_qsat_2d (iim,jjm,tmin_nm1,pb,qsattmin_nm1)
         CALL weathgen_qsat_2d (iim,jjm,tair,pb,qsatta)

         !---
         qmin_nm1(:,:) = MIN(qair_nm1(:,:),0.99*qsattmin_nm1(:,:))
         qmin_n(:,:) = MIN(qair_n(:,:),0.99*qsattmin_n(:,:))
         qmax_nm1(:,:) = (qair_nm1(:,:)-qmin_nm1(:,:)) + qair_nm1(:,:)
         qmax_n(:,:) = (qair_n(:,:)-qmin_n(:,:)) + qair_n(:,:)

         qsa(:,:)  = 0.99*qsatta(:,:)


         WHERE ( ( hour >= startday_n(:,:) ) .AND. ( hour > 12) .AND. ( hour <= 14) )
            qair(:,:)=MIN(qsa(:,:),(qmax_nm1(:,:)-qmin_nm1(:,:))/2 * ( sin(pi/(14-startday_n(:,:))*(hour-0.5* &
           &    (14.-startday_n(:,:))-startday_n(:,:))) )+ (qmax_nm1(:,:)+qmin_nm1(:,:))/2.)
         ELSEWHERE( ( hour >= startday_n(:,:) ) .AND. ( hour <= 12) )
            qair(:,:)=MIN(qsa(:,:),(qmax_n(:,:)-qmin_n(:,:))/2 * ( sin(pi/(14-startday_n(:,:))*(hour-0.5* &
           &    (14.-startday_n(:,:))-startday_n(:,:))) )+ (qmax_n(:,:)+qmin_n(:,:))/2.)
         ELSEWHERE ( hour < startday_n(:,:) )
            qair(:,:)=MIN(qsa(:,:),(qmax_nm1(:,:)-qmin_n(:,:))/2.*sin(pi/(24.-14.+startday_nm1(:,:) )* &
           &    (hour + 24.+0.5*(24.-14.+startday_nm1(:,:) )-14.))+(qmax_nm1(:,:)+qmin_n(:,:))/2.)
         ELSEWHERE
            qair(:,:)=MIN(qsa(:,:),(qmax_nm1(:,:)-qmin_n(:,:))/2.*sin(pi/(24.-14.+startday_n(:,:))*(hour+0.5* &
           &    (24.-14.+startday_n(:,:))-14.))+(qmax_nm1(:,:)+qmin_n(:,:))/2.)
         ENDWHERE

         IF (is_watchout) THEN
            SWnet(:,:) = (SWnet_n(:,:)-SWnet_nm1(:,:))*rw + SWnet_nm1(:,:)
            Eair(:,:) = (Eair_n(:,:)-Eair_nm1(:,:))*rw + Eair_nm1(:,:)
            petAcoef(:,:) = (petAcoef_n(:,:)-petAcoef_nm1(:,:))*rw + petAcoef_nm1(:,:)
            peqAcoef(:,:) = (peqAcoef_n(:,:)-peqAcoef_nm1(:,:))*rw + peqAcoef_nm1(:,:)
            petBcoef(:,:) = (petBcoef_n(:,:)-petBcoef_nm1(:,:))*rw + petBcoef_nm1(:,:)
            peqBcoef(:,:) = (peqBcoef_n(:,:)-peqBcoef_nm1(:,:))*rw + peqBcoef_nm1(:,:)
            cdrag(:,:) = (cdrag_n(:,:)-cdrag_nm1(:,:))*rw + cdrag_nm1(:,:)
            ccanopy(:,:) = (ccanopy_n(:,:)-ccanopy_nm1(:,:))*rw + ccanopy_nm1(:,:)
         ENDIF
   !---
   !--- Here we need to allow for an option
   !--- where radiative energy is conserved
   !---
         IF ( lwdown_cons ) THEN
            lwdown(:,:) = lwdown_n(:,:)
         ELSE
            lwdown(:,:) = (lwdown_n(:,:)-lwdown_nm1(:,:))*rw + lwdown_nm1(:,:)
         ENDIF
   !---
   !--- For the solar radiation we decompose the mean value
   !--- using the zenith angle of the sun, conservative approach under 2000W/m2
   !----
         IF (printlev_loc >= 5) WRITE(numout,*) 'Ready to deal with the solar radiation'
   !----
         ! We compute mean SWdown between t and t+Dt then we take t+Dt/2. 
         julian = julian_for + (itau_split-0.5)/split*dt_force/one_day
!!$         julian = julian_for + rw*dt_force/one_day
         IF (printlev_loc >= 5) THEN
            WRITE(numout,'(a,f20.10,2I3)') &
                 &  'JULIAN BEFORE SOLARANG : ',julian,itau_split,split
         ENDIF
 
         CALL solarang(julian, julian0, iim, jjm, lon*0, lat, coszang)
 
         WHERE ((mean_coszang(:,:) > 0.) .AND. (hour <= 12 ))
            swdown(:,:) = swdown_n(:,:) *coszang(:,:)/mean_coszang(:,:)
         ELSEWHERE ((mean_coszang(:,:) > 0.) .AND. (hour > 12 ))
            swdown(:,:) = swdown_nm1(:,:) *coszang(:,:)/mean_coszang(:,:)
         ELSEWHERE
            swdown(:,:) = 0.0
         END WHERE
 
         WHERE (swdown(:,:) > 2000. )
            swdown(:,:) = 2000.
         END WHERE
         
         IF (printlev_loc >= 5) THEN
            WRITE(numout,*) '__ Forcing read at ',itau_split,' :',i_test, j_test
            WRITE(numout,*) 'SWdown  : ',swdown_nm1(i_test, j_test), &
                 &           ' < ', swdown(i_test, j_test), ' < ', swdown_n(i_test, j_test)
            IF (is_watchout) THEN
               WRITE(numout,*) 'SWnet  : ',swnet_nm1(i_test, j_test), &
                    &           ' < ', swnet(i_test, j_test), ' < ', swnet_n(i_test, j_test)
               WRITE(numout,*) 'levels  :',zlev_nm1(i_test, j_test), &
                    &           ' < ', zlev(i_test, j_test), ' < ', zlev_n(i_test, j_test)
               WRITE(numout,*) 'EAIR  :',Eair_nm1(i_test, j_test), &
                    &           ' < ', eair(i_test, j_test), ' < ', Eair_n(i_test, j_test)
            ENDIF
            WRITE(numout,*) 'TAIR  :',tair_nm1(i_test, j_test), &
                 &           ' < ', tair(i_test, j_test), ' < ', tair_n(i_test, j_test)
            WRITE(numout,*) 'QAIR  :',qair_nm1(i_test, j_test), &
                 &           ' < ', qair(i_test, j_test), ' < ', qair_n(i_test, j_test)
            WRITE(numout,*) 'U  :',u_nm1(i_test, j_test), &
                 &           ' < ', u(i_test, j_test), ' < ', u_n(i_test, j_test)
            WRITE(numout,*) 'V  :',v_nm1(i_test, j_test), &
                 &           ' < ', v(i_test, j_test), ' < ', v_n(i_test, j_test)
         ENDIF
   !---
   !--- For precip we suppose that the rain
   !--- is the sum over the next 6 hours
   !---
         WHERE ((itau_split <= nb_spread).AND.(hour<=12).AND.(tair(:,:)>=273.15)) 
            rainf(:,:) = rainf_n(:,:) *(split/REAL(nb_spread))
            snowf(:,:) = 0.0
         ELSEWHERE ((itau_split <= nb_spread).AND.(hour<=12).AND.(tair(:,:)<273.15)) 
            snowf(:,:) = rainf_n(:,:) *(split/REAL(nb_spread))
            rainf(:,:) = 0.0
         ELSEWHERE ((itau_split <= nb_spread).AND.(hour>12).AND.(tair(:,:)>=273.15)) 
            rainf(:,:) = rainf_nm1(:,:) *(split/REAL(nb_spread))
            snowf(:,:) = 0.0
         ELSEWHERE ((itau_split <= nb_spread).AND.(hour>12).AND.(tair(:,:)<273.15)) 
            snowf(:,:) = rainf_nm1(:,:) *(split/REAL(nb_spread))
            rainf(:,:) = 0.0
         ELSEWHERE
            snowf(:,:) = 0.0
            rainf(:,:) = 0.0
         ENDWHERE

         IF (printlev_loc >= 5) THEN
            WRITE(numout,*) '__ Forcing read at ',itau_split,' :'
            WRITE(numout,*) 'Rainf  : ',rainf_nm1(i_test, j_test), &
                 &           ' < ', rainf(i_test, j_test), ' < ', rainf_n(i_test, j_test)
            WRITE(numout,*) 'Snowf  : ',snowf_nm1(i_test, j_test), &
                 &           ' < ', snowf(i_test, j_test), ' < ', snowf_n(i_test, j_test)
         ENDIF
   !---

      ELSE ! If not daily_interpol
         
         IF (itau_read /= last_read) THEN
   !---
   !-----   Start or Restart
            IF (itau_read_n == 0) THEN
               ! Case of a restart or a shift in the forcing file.
               IF (itau_read > 1) THEN
                  itau_read_nm1=itau_read-1
                  CALL forcing_just_read (iim, jjm, zlev_nm1, zlevuv_nm1, ttm, itau_read_nm1, itau_read_nm1, &
                       &  swdown_nm1, rainf_nm1, snowf_nm1, tair_nm1, &
                       &  u_nm1, v_nm1, qair_nm1, pb_nm1, lwdown_nm1, &
                       &  SWnet_nm1, Eair_nm1, petAcoef_nm1, peqAcoef_nm1, petBcoef_nm1, peqBcoef_nm1, cdrag_nm1, ccanopy_nm1, &
                       &  force_id, wind_N_exists)
               ! Case of a simple start.
               ELSE IF (dt_force*ttm > one_day-1. ) THEN
                  ! if the forcing file contains at least 24 hours, 
                  ! we will use the last forcing step of the first day
                  ! as initiale condition to prevent first shift off reading.
                  itau_read_nm1 = NINT (one_day/dt_force)
                  IF (printlev_loc >= 1) WRITE(numout,*) "The forcing file contains 24 hours :",dt_force*ttm,one_day-1.
                  IF (printlev_loc >= 1) WRITE(numout,*) "We will use the last forcing step of the first day : itau_read_nm1 ",&
                       itau_read_nm1
                  CALL forcing_just_read (iim, jjm, zlev_nm1, zlevuv_nm1, ttm, itau_read_nm1, itau_read_nm1, &
                       &  swdown_nm1, rainf_nm1, snowf_nm1, tair_nm1, &
                       &  u_nm1, v_nm1, qair_nm1, pb_nm1, lwdown_nm1, &
                       &  SWnet_nm1, Eair_nm1, petAcoef_nm1, peqAcoef_nm1, petBcoef_nm1, peqBcoef_nm1, cdrag_nm1, ccanopy_nm1, &
                       &  force_id, wind_N_exists)
               ELSE
                  ! if the forcing file contains less than 24 hours, 
                  ! just say error !
                  CALL ipslerr_p(3,'forcing_read_interpol', &
      &             'The forcing file contains less than 24 hours !', &
      &             'We can''t intialize interpolation with such a file.','') 
               ENDIF
            ELSE
   !-----   Normal mode : copy old step
               swdown_nm1(:,:) = swdown_n(:,:)
               rainf_nm1(:,:) = rainf_n(:,:)
               snowf_nm1(:,:)  = snowf_n(:,:) 
               tair_nm1(:,:)   = tair_n(:,:)
               u_nm1(:,:)      = u_n(:,:)
               v_nm1(:,:)      = v_n(:,:)
               qair_nm1(:,:)   = qair_n(:,:)
               pb_nm1(:,:)     = pb_n(:,:)
               lwdown_nm1(:,:) = lwdown_n(:,:)
               IF (is_watchout) THEN
                  zlev_nm1(:,:)   = zlev_n(:,:)
                  ! Net surface short-wave flux
                  SWnet_nm1(:,:) = SWnet_n(:,:)
                  ! Air potential energy
                  Eair_nm1(:,:)   = Eair_n(:,:)
                  ! Coeficients A from the PBL resolution for T
                  petAcoef_nm1(:,:) = petAcoef_n(:,:)
                  ! Coeficients A from the PBL resolution for q
                  peqAcoef_nm1(:,:) = peqAcoef_n(:,:)
                  ! Coeficients B from the PBL resolution for T
                  petBcoef_nm1(:,:) = petBcoef_n(:,:)
                  ! Coeficients B from the PBL resolution for q
                  peqBcoef_nm1(:,:) = peqBcoef_n(:,:)
                  ! Surface drag
                  cdrag_nm1(:,:) = cdrag_n(:,:)
                  ! CO2 concentration in the canopy
                  ccanopy_nm1(:,:) = ccanopy_n(:,:)
               ENDIF
               itau_read_nm1 = itau_read_n
            ENDIF
   !-----
            itau_read_n = itau_read
            IF (itau_read_n > ttm) THEN
               WRITE(numout,*) 'WARNING --WARNING --WARNING --WARNING '
               WRITE(numout,*) &
                    &  'WARNING : We are going back to the start of the file'
               itau_read_n =1
            ENDIF
            IF (printlev_loc >= 5) THEN
               WRITE(numout,*) &
                    & 'The dates 2 : ',itau_read,itau_split,itau_read_nm1,itau_read_n
            ENDIF
   !-----
   !----- Get a reduced julian day !
   !----- This is needed because we lack the precision on 32 bit machines.
   !-----
            IF ( dt_force .GT. 3600. ) THEN
               julian_for = itau2date(itau_read-1, date0, dt_force)
               CALL ju2ymds (julian_for, yy, mm, dd, ss)
   
               ! first day of this year
               CALL ymds2ju (yy,1,1,0.0, julian0)
   !-----
               IF (printlev_loc >= 5) THEN
                  WRITE(numout,*) 'Forcing for Julian day ',julian_for,'is read'
                  WRITE(numout,*) 'Date for this day ',yy,' / ',mm,' / ',dd,"  ",ss
               ENDIF
            ENDIF
   !-----
            CALL forcing_just_read (iim, jjm, zlev_n, zlevuv_n, ttm, itau_read_n, itau_read_n, &
                 &  swdown_n, rainf_n, snowf_n, tair_n, &
                 &  u_n, v_n, qair_n, pb_n, lwdown_n, &
                 &  SWnet_n, Eair_n, petAcoef_n, peqAcoef_n, petBcoef_n, peqBcoef_n, cdrag_n, ccanopy_n, &
                 &  force_id, wind_N_exists)
   !---
            last_read = itau_read_n
   !-----
   !----- Compute mean solar angle for the comming period
   !-----
            IF (printlev_loc >= 5) WRITE(numout,*) 'Going into  solarang', split, one_day
   !-----
            IF ( dt_force .GT. 3600. ) THEN
               mean_coszang(:,:) = 0.0
               DO is=1,split
                  !MM we compute mean SWdown between t and t+Dt then I take t+Dt/2. 
                  julian = julian_for+((is-0.5)/split)*dt_force/one_day
   !!$               julian = julian_for+(FLOAT(is)/split)*dt_force/one_day
                  CALL solarang (julian, julian0, iim, jjm, lon, lat, coszang)
                  mean_coszang(:,:) = mean_coszang(:,:)+coszang(:,:)
               ENDDO
               mean_coszang(:,:) = mean_coszang(:,:)/split
   !            WRITE(*,*) "mean_coszang =",MAXVAL(mean_coszang)
            ENDIF
   !-----
         ENDIF
   !---
   !--- Do the interpolation
         IF (printlev_loc >= 5) WRITE(numout,*) 'Doing the interpolation between time steps'
   !---
         IF (split > 1) THEN
            rw = REAL(itau_split)/split
         ELSE
            rw = 1.
         ENDIF
         IF (printlev_loc >= 5) WRITE(numout,*) 'Coeff of interpollation : ',rw
   !---
         qair(:,:) = (qair_n(:,:)-qair_nm1(:,:))*rw + qair_nm1(:,:)
         tair(:,:) = (tair_n(:,:)-tair_nm1(:,:))*rw + tair_nm1(:,:)
         pb(:,:) = (pb_n(:,:)-pb_nm1(:,:))*rw + pb_nm1(:,:)
         u(:,:)  = (u_n(:,:)-u_nm1(:,:))*rw + u_nm1(:,:)
         v(:,:)  = (v_n(:,:)-v_nm1(:,:))*rw + v_nm1(:,:)
         IF (is_watchout) THEN
            zlev(:,:) = (zlev_n(:,:)-zlev_nm1(:,:))*rw + zlev_nm1(:,:)
            zlevuv(:,:) = zlev(:,:)
            SWnet(:,:) = (SWnet_n(:,:)-SWnet_nm1(:,:))*rw + SWnet_nm1(:,:)
            Eair(:,:) = (Eair_n(:,:)-Eair_nm1(:,:))*rw + Eair_nm1(:,:)
            petAcoef(:,:) = (petAcoef_n(:,:)-petAcoef_nm1(:,:))*rw + petAcoef_nm1(:,:)
            peqAcoef(:,:) = (peqAcoef_n(:,:)-peqAcoef_nm1(:,:))*rw + peqAcoef_nm1(:,:)
            petBcoef(:,:) = (petBcoef_n(:,:)-petBcoef_nm1(:,:))*rw + petBcoef_nm1(:,:)
            peqBcoef(:,:) = (peqBcoef_n(:,:)-peqBcoef_nm1(:,:))*rw + peqBcoef_nm1(:,:)
            cdrag(:,:) = (cdrag_n(:,:)-cdrag_nm1(:,:))*rw + cdrag_nm1(:,:)
            ccanopy(:,:) = (ccanopy_n(:,:)-ccanopy_nm1(:,:))*rw + ccanopy_nm1(:,:)
         ENDIF
   !---
   !--- Here we need to allow for an option
   !--- where radiative energy is conserved
   !---
         IF ( lwdown_cons ) THEN
            lwdown(:,:) = lwdown_n(:,:)
         ELSE
            lwdown(:,:) = (lwdown_n(:,:)-lwdown_nm1(:,:))*rw + lwdown_nm1(:,:)
         ENDIF
   !---
   !--- For the solar radiation we decompose the mean value
   !--- using the zenith angle of the sun if the time step in the forcing data is 
   !---- more than an hour. Else we use the standard linera interpolation
   !----
         IF (printlev_loc >= 5) WRITE(numout,*) 'Ready to deal with the solar radiation'
   !----
         IF ( dt_force .GT. 3600. ) THEN

            ! In this case solar radiation will be conserved
  
            ! We compute mean SWdown between t and t+Dt then we take t+Dt/2. 
            julian = julian_for + (itau_split-0.5)/split*dt_force/one_day
   !!$         julian = julian_for + rw*dt_force/one_day
            IF (printlev_loc >= 5) THEN
               WRITE(numout,'(a,f20.10,2I3)') &
                    &  'JULIAN BEFORE SOLARANG : ',julian,itau_split,split
            ENDIF
   !---
            CALL solarang(julian, julian0, iim, jjm, lon, lat, coszang)
   !---
            WHERE (mean_coszang(:,:) > 0.)
               swdown(:,:) = swdown_n(:,:) *coszang(:,:)/mean_coszang(:,:)
            ELSEWHERE
               swdown(:,:) = 0.0
            END WHERE
   !---
            WHERE (swdown(:,:) > 2000. )
               swdown(:,:) = 2000.
            END WHERE
   !---
         ELSE ! If dt_force < 3600 (1h)

            IF ( swdown_cons ) THEN
               ! Conserve swdown radiation
               swdown(:,:) = swdown_n(:,:)
            ELSE
               swdown(:,:) = (swdown_n(:,:)-swdown_nm1(:,:))*rw + swdown_nm1(:,:)
            ENDIF
   !---
         ENDIF
   !---
         IF (printlev_loc >= 5) THEN
            WRITE(numout,*) '__ Forcing read at ',itau_split,' :',i_test, j_test
            WRITE(numout,*) 'SWdown  : ',swdown_nm1(i_test, j_test), &
                 &           ' < ', swdown(i_test, j_test), ' < ', swdown_n(i_test, j_test)
            IF (is_watchout) THEN
               WRITE(numout,*) 'SWnet  : ',swnet_nm1(i_test, j_test), &
                    &           ' < ', swnet(i_test, j_test), ' < ', swnet_n(i_test, j_test)
               WRITE(numout,*) 'levels  :',zlev_nm1(i_test, j_test), &
                    &           ' < ', zlev(i_test, j_test), ' < ', zlev_n(i_test, j_test)
               WRITE(numout,*) 'EAIR  :',Eair_nm1(i_test, j_test), &
                    &           ' < ', eair(i_test, j_test), ' < ', Eair_n(i_test, j_test)
            ENDIF
            WRITE(numout,*) 'TAIR  :',tair_nm1(i_test, j_test), &
                 &           ' < ', tair(i_test, j_test), ' < ', tair_n(i_test, j_test)
            WRITE(numout,*) 'QAIR  :',qair_nm1(i_test, j_test), &
                 &           ' < ', qair(i_test, j_test), ' < ', qair_n(i_test, j_test)
            WRITE(numout,*) 'U  :',u_nm1(i_test, j_test), &
                 &           ' < ', u(i_test, j_test), ' < ', u_n(i_test, j_test)
            WRITE(numout,*) 'V  :',v_nm1(i_test, j_test), &
                 &           ' < ', v(i_test, j_test), ' < ', v_n(i_test, j_test)
         ENDIF
   !---
   !--- For precip we suppose that the rain
   !--- is the sum over the next 6 hours
   !---
         IF (itau_split <= nb_spread) THEN
            rainf(:,:) = rainf_n(:,:)*(split/REAL(nb_spread))
            snowf(:,:) = snowf_n(:,:)*(split/REAL(nb_spread))
         ELSE
            rainf(:,:) = 0.0
            snowf(:,:) = 0.0
         ENDIF
         IF (printlev_loc >= 5) THEN
            WRITE(numout,*) '__ Forcing read at ',itau_split,' :'
            WRITE(numout,*) 'Rainf  : ',rainf_nm1(i_test, j_test), &
                 &           ' < ', rainf(i_test, j_test), ' < ', rainf_n(i_test, j_test)
            WRITE(numout,*) 'Snowf  : ',snowf_nm1(i_test, j_test), &
                 &           ' < ', snowf(i_test, j_test), ' < ', snowf_n(i_test, j_test)
         ENDIF
   !---
      ENDIF ! (daily_interpol)
   ENDIF
!---
!--- Here we might put the call to the weather generator ... one day.
!--- Pour le moment, le branchement entre interpolation et generateur de temps 
!--- est fait au-dessus.
!---
!-   IF ( initialized .AND. weathergen ) THEN
!-      ....
!-   ENDIF
!---
!--- At this point the code should be initialized. If not we have a problem !
!---
   IF ( (itau_read == 0).AND.(itau_split == 0) ) THEN
!---
      initialized = .TRUE.
!---
   ELSE
      IF ( .NOT. initialized ) THEN
         WRITE(numout,*) 'Why is the code forcing_read not initialized ?'
         WRITE(numout,*) 'Have you called it with both time-steps set to zero ?'
         CALL ipslerr_p(3,'forcing_read_interpol','Pb in initialization','','')
      ENDIF
   ENDIF

END SUBROUTINE forcing_read_interpol


!! ==============================================================================================================================\n
!! SUBROUTINE 	: forcing_just_read
!!
!>\BRIEF        
!!
!!\n DESCRIPTION : 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S)	: 
!!
!_ ================================================================================================================================
SUBROUTINE forcing_just_read &
  & (iim, jjm, zlev, zlev_uv, ttm, itb, ite, &
  &  swdown, rainf, snowf, tair, &
  &  u, v, qair, pb, lwdown, &
  &  SWnet, Eair, petAcoef, peqAcoef, petBcoef, peqBcoef, cdrag, ccanopy, &
  &  force_id, wind_N_exists)
!---------------------------------------------------------------------
!- iim        : Size of the grid in x
!- jjm        : size of the grid in y
!- zlev       : height of the varibales T and Q 
!- zlev_uv    : height of the varibales U and V 
!- ttm        : number of time steps in all in the forcing file
!- itb, ite   : index of respectively begin and end of read for each variable
!- swdown     : Downward solar radiation (W/m^2)
!- rainf      : Rainfall (kg/m^2s)
!- snowf      : Snowfall (kg/m^2s)
!- tair       : 2m air temperature (K)
!- u and v    : 2m (in theory !) wind speed (m/s)
!- qair       : 2m humidity (kg/kg)
!- pb         : Surface pressure (Pa)
!- lwdown     : Downward long wave radiation (W/m^2)
!-
!- From a WATCHOUT file :
!- SWnet      : Net surface short-wave flux
!- Eair       : Air potential energy
!- petAcoef   : Coeficients A from the PBL resolution for T
!- peqAcoef   : Coeficients A from the PBL resolution for q
!- petBcoef   : Coeficients B from the PBL resolution for T
!- peqBcoef   : Coeficients B from the PBL resolution for q
!- cdrag      : Surface drag
!- ccanopy    : CO2 concentration in the canopy
!- force_id   : FLINCOM file id.
!-              It is used to close the file at the end of the run.
!- wind_N_exists : if Wind_N and Wind_E are in the file (and not just Wind)
!---------------------------------------------------------------------
   IMPLICIT NONE
!-
   INTEGER, INTENT(in) :: iim, jjm, ttm
   INTEGER, INTENT(in) :: itb, ite
   REAL, DIMENSION(iim,jjm), INTENT(out) ::  zlev, zlev_uv, &
  &  swdown, rainf, snowf, tair, u, v, qair, pb, lwdown
   ! for watchout files
   REAL, DIMENSION(iim,jjm), INTENT(out) :: &
  &  SWnet, Eair, petAcoef, peqAcoef, petBcoef, peqBcoef, cdrag, ccanopy
   INTEGER, INTENT(in) :: force_id
!  if Wind_N and Wind_E are in the file (and not just Wind)
   LOGICAL, INTENT(in) :: wind_N_exists
   INTEGER :: i, j 
   REAL :: rau 

!-
!---------------------------------------------------------------------
   IF ( daily_interpol ) THEN
      CALL flinget_buffer (force_id,'Tmin'  , iim_full, jjm_full, llm_full, ttm, itb, ite, data_full)
      CALL forcing_zoom(data_full, tair)
      CALL flinget_buffer (force_id,'precip' , iim_full, jjm_full, llm_full, ttm, itb, ite, data_full)
      CALL forcing_zoom(data_full, rainf)
   ELSE
      CALL flinget_buffer (force_id,'Tair'  , iim_full, jjm_full, llm_full, ttm, itb, ite, data_full)
      CALL forcing_zoom(data_full, tair)
      CALL flinget_buffer (force_id,'Snowf' , iim_full, jjm_full, llm_full, ttm, itb, ite, data_full)
      CALL forcing_zoom(data_full, snowf)
      CALL flinget_buffer (force_id,'Rainf' , iim_full, jjm_full, llm_full, ttm, itb, ite, data_full)
      CALL forcing_zoom(data_full, rainf)
   ENDIF


   CALL flinget_buffer (force_id,'SWdown', iim_full, jjm_full, llm_full, ttm, itb, ite, data_full)
   CALL forcing_zoom(data_full, swdown)
   CALL flinget_buffer (force_id,'LWdown', iim_full, jjm_full, llm_full, ttm, itb, ite, data_full)
   CALL forcing_zoom(data_full, lwdown)

   CALL flinget_buffer (force_id,'PSurf' , iim_full, jjm_full, llm_full, ttm, itb, ite, data_full)
   CALL forcing_zoom(data_full, pb)
   CALL flinget_buffer (force_id,'Qair'  , iim_full, jjm_full, llm_full, ttm, itb, ite, data_full)
   CALL forcing_zoom(data_full, qair)
!---
   IF ( wind_N_exists ) THEN
      CALL flinget_buffer (force_id,'Wind_N', iim_full, jjm_full, llm_full, ttm, itb, ite, data_full)
      CALL forcing_zoom(data_full, u)
      CALL flinget_buffer (force_id,'Wind_E', iim_full, jjm_full, llm_full, ttm, itb, ite, data_full)
      CALL forcing_zoom(data_full, v)
   ELSE
      CALL flinget_buffer (force_id,'Wind',   iim_full, jjm_full, llm_full, ttm, itb, ite, data_full)
      CALL forcing_zoom(data_full, u)
      v=0.0
   ENDIF

!- 
!- Deal with the height of the atmospheric forcing varibles 
!- 
!---- 
   IF ( zheight ) THEN 
      zlev(:,:) = zlev_fixed 
   ELSE IF ( zsigma .OR. zhybrid ) THEN 
      DO i=1,iim 
         DO j=1,jjm 
            IF ( tair(i,j) < val_exp ) THEN 
               rau = pb(i,j)/(cte_molr*tair(i,j)) 
 	  	 
               zlev(i,j) =  (pb(i,j) - (zhybrid_a + zhybrid_b*pb(i,j)))/(rau * cte_grav) 
            ELSE 
               zlev(i,j) = 0.0 
            ENDIF
         ENDDO
      ENDDO
   ELSE IF ( zlevels ) THEN 
      CALL flinget_buffer (force_id,'Levels', iim_full, jjm_full, llm_full, ttm, itb, ite, data_full) 
      CALL forcing_zoom(data_full, zlev) 
   ELSE 
      CALL ipslerr(3, 'forcing_just_read','No case for the vertical levels was specified.', & 
           &         'We cannot determine the height for T and Q.','stop readdim2') 
   ENDIF
   
   IF ( zsamelev_uv ) THEN 
      zlev_uv(:,:) = zlev(:,:) 
   ELSE 
      IF ( zheight ) THEN 
         zlev_uv(:,:) = zlevuv_fixed 
      ELSE IF ( zsigma .OR. zhybrid ) THEN 
         DO i=1,iim 
            DO j=1,jjm 
               IF ( tair(i,j) < val_exp ) THEN 
                  rau = pb(i,j)/(cte_molr*tair(i,j)) 
                  
                  zlev_uv(i,j) =  (pb(i,j) - (zhybriduv_a + zhybriduv_b*pb(i,j)))/(rau * cte_grav) 
               ELSE 
                  zlev_uv(i,j) = 0.0 
               ENDIF
            ENDDO
         ENDDO
      ELSE IF ( zlevels ) THEN 
         CALL flinget_buffer (force_id,'Levels_uv', iim_full, jjm_full, llm_full, ttm, itb, ite, data_full) 
         CALL forcing_zoom(data_full, zlev_uv) 
      ELSE 
         CALL ipslerr(3, 'forcing_just_read','No case for the vertical levels was specified.', & 
              &         'We cannot determine the height for U and V.','stop readdim2') 
      ENDIF
   ENDIF
   !----
   IF ( is_watchout ) THEN
      CALL flinget_buffer (force_id,'levels', iim_full, jjm_full, llm_full, ttm, itb, ite, data_full)
      CALL forcing_zoom(data_full, zlev)
      ! 
      ! If we are in WATHCOUT it means T,Q are at the same height as U,V 
      ! 
      zlev_uv(:,:) = zlev(:,:) 
      ! Net surface short-wave flux
      CALL flinget_buffer (force_id,'SWnet', iim_full, jjm_full, llm_full, ttm, itb, ite, data_full)
      CALL forcing_zoom(data_full, SWnet)
      ! Air potential energy
      CALL flinget_buffer (force_id,'Eair', iim_full, jjm_full, llm_full, ttm, itb, ite, data_full)
      CALL forcing_zoom(data_full, Eair)
      ! Coeficients A from the PBL resolution for T
      CALL flinget_buffer (force_id,'petAcoef', iim_full, jjm_full, llm_full, ttm, itb, ite, data_full)
      CALL forcing_zoom(data_full, petAcoef)
      ! Coeficients A from the PBL resolution for q
      CALL flinget_buffer (force_id,'peqAcoef', iim_full, jjm_full, llm_full, ttm, itb, ite, data_full)
      CALL forcing_zoom(data_full, peqAcoef)
      ! Coeficients B from the PBL resolution for T
      CALL flinget_buffer (force_id,'petBcoef', iim_full, jjm_full, llm_full, ttm, itb, ite, data_full)
      CALL forcing_zoom(data_full, petBcoef)
      ! Coeficients B from the PBL resolution for q
      CALL flinget_buffer (force_id,'peqBcoef', iim_full, jjm_full, llm_full, ttm, itb, ite, data_full)
      CALL forcing_zoom(data_full, peqBcoef)
      ! Surface drag
      CALL flinget_buffer (force_id,'cdrag', iim_full, jjm_full, llm_full, ttm, itb, ite, data_full)
      CALL forcing_zoom(data_full, cdrag)
      ! CO2 concentration in the canopy
      CALL flinget_buffer (force_id,'ccanopy', iim_full, jjm_full, llm_full, ttm, itb, ite, data_full)
      CALL forcing_zoom(data_full, ccanopy)
   ENDIF
!
!----
     IF (printlev_loc >= 5) WRITE(numout,*) 'Variables have been extracted between ',itb,&
          ' and ',ite,' iterations of the forcing file.'
!-------------------------
END SUBROUTINE forcing_just_read


!! ==============================================================================================================================\n
!! SUBROUTINE 	: forcing_just_read_tmax
!!
!>\BRIEF        
!!
!!\n DESCRIPTION : 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S)	: 
!!
!_ ================================================================================================================================
SUBROUTINE forcing_just_read_tmax &
  & (iim, jjm, ttm, itb, ite, tmax, force_id )
!---------------------------------------------------------------------
!- iim        : Size of the grid in x
!- jjm        : size of the grid in y
!- ttm        : number of time steps in all in the forcing file
!- itb, ite   : index of respectively begin and end of read for each variable
!- tmax       : 2m air temperature (K)
!- force_id   : FLINCOM file id.
!-              It is used to close the file at the end of the run.
!---------------------------------------------------------------------
   IMPLICIT NONE
!-
   INTEGER, INTENT(in) :: iim, jjm, ttm
   INTEGER, INTENT(in) :: itb, ite
   REAL, DIMENSION(iim,jjm), INTENT(out) ::  tmax
   INTEGER, INTENT(in) :: force_id
!-
!---------------------------------------------------------------------
   CALL flinget_buffer (force_id,'Tmax'  , iim_full, jjm_full, llm_full, ttm, itb, ite, data_full)
   CALL forcing_zoom(data_full, tmax)
!-------------------------
END SUBROUTINE forcing_just_read_tmax


!! ==============================================================================================================================\n
!! SUBROUTINE 	: forcing_landind
!!
!>\BRIEF        
!!
!!\n DESCRIPTION : This subroutine finds the indices of the land points over which the land
!!                 surface scheme is going to run.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S)	: 
!!
!_ ================================================================================================================================
SUBROUTINE forcing_landind(iim, jjm, tair, nbindex, kindex, i_test, j_test)
!---
!--- 
!---
  IMPLICIT NONE
!-
!- ARGUMENTS
!-
  INTEGER, INTENT(IN) :: iim, jjm
  REAL, INTENT(IN)    :: tair(iim,jjm)
  INTEGER, INTENT(OUT) :: i_test, j_test, nbindex
  INTEGER, INTENT(OUT) :: kindex(iim*jjm)
!-
!- LOCAL
  INTEGER :: i, j, ig
!-
!-
  ig = 0
  i_test = 0
  j_test = 0
!---
  IF (MINVAL(tair(:,:)) < 100.) THEN
!----- In this case the 2m temperature is in Celsius
     DO j=1,jjm
        DO i=1,iim
           IF (tair(i,j) < 100.) THEN
              ig = ig+1
              kindex(ig) = (j-1)*iim+i
              !
              !  Here we find at random a land-point on which we can do
              !  some printouts for test.
              !
              IF (ig .GT. (iim*jjm)/2 .AND. i_test .LT. 1) THEN
                 i_test = i
                 j_test = j
                 IF (printlev_loc >= 5) THEN
                    WRITE(numout,*) 'The test point chosen for output is : ', i_test, j_test
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
     ENDDO
  ELSE 
!----- 2m temperature is in Kelvin
     DO j=1,jjm
        DO i=1,iim
           IF (tair(i,j) < 500.) THEN
              ig = ig+1
              kindex(ig) = (j-1)*iim+i
              !
              !  Here we find at random a land-point on which we can do
              !  some printouts for test.
              !
              IF (ig .GT. (iim*jjm)/2 .AND. i_test .LT. 1) THEN
                 i_test = i
                 j_test = j
                 IF (printlev_loc >= 5) THEN
                    WRITE(numout,*) 'The test point chosen for output is : ', i_test, j_test
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
     ENDDO
  ENDIF

  nbindex = ig

END SUBROUTINE forcing_landind


!! ==============================================================================================================================\n
!! SUBROUTINE 	: forcing_grid
!!
!>\BRIEF        
!!
!!\n DESCRIPTION : This subroutine calculates the longitudes and latitudes of the model grid.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S)	: 
!!
!_ ================================================================================================================================
SUBROUTINE forcing_grid(iim,jjm,llm,lon,lat,init_f)

  IMPLICIT NONE

  INTEGER, INTENT(in)                   :: iim, jjm, llm
  LOGICAL, INTENT(in)                   :: init_f
  REAL, DIMENSION(iim,jjm), INTENT(out) :: lon, lat
!-
  INTEGER :: i,j
!-
!- Should be unified one day
!-
  IF ( printlev_loc>=3 ) WRITE(numout,*) 'forcing_grid : options : ', weathergen, interpol
!-
  IF ( weathergen ) THEN
     IF (init_f) THEN
        DO i = 1, iim
           lon(i,:) = limit_west + merid_res/2. + &
                FLOAT(i-1)*(limit_east-limit_west)/FLOAT(iim)
        ENDDO
        !-
        DO j = 1, jjm
           lat(:,j) = limit_north - zonal_res/2. - &
                FLOAT(j-1)*(limit_north-limit_south)/FLOAT(jjm)
        ENDDO
     ELSE
        IF (is_root_prc) THEN
           DO i = 1, iim_g
              lon_g(i,:) = limit_west + merid_res/2. + &
                   FLOAT(i-1)*(limit_east-limit_west)/FLOAT(iim_g)
           ENDDO
           !-
           DO j = 1, jjm_g
              lat_g(:,j) = limit_north - zonal_res/2. - &
                   FLOAT(j-1)*(limit_north-limit_south)/FLOAT(jjm_g)
           ENDDO
        ENDIF
        CALL bcast(lon_g)
        CALL bcast(lat_g)
        lon=lon_g(:,jj_para_begin(mpi_rank):jj_para_end(mpi_rank))
        lat=lat_g(:,jj_para_begin(mpi_rank):jj_para_end(mpi_rank))
     ENDIF
!-
  ELSEIF ( interpol ) THEN
!-
    CALL forcing_zoom(lon_full, lon)
    IF ( printlev_loc>=3 ) WRITE(numout,*) 'forcing_grid : out of zoom on lon'
    CALL forcing_zoom(lat_full, lat)
    IF ( printlev_loc>=3 ) WRITE(numout,*) 'forcing_grid : out of zoom on lat'
 
 ELSE 
    CALL ipslerr_p(3,'forcing_grid','Neither interpolation nor weather generator is specified.','','')
 ENDIF
 
 IF ( printlev_loc>=3 ) WRITE(numout,*) 'forcing_grid : ended'
 
END SUBROUTINE forcing_grid


!! ==============================================================================================================================\n
!! SUBROUTINE 	: forcing_zoom
!!
!>\BRIEF        This subroutine extracts the region of data over which we wish to run the model.
!!
!!\n DESCRIPTION : This subroutine extracts the region of data over which we wish to run the model.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S)	: 
!!
!_ ================================================================================================================================
SUBROUTINE forcing_zoom(x_f, x_z)

  IMPLICIT NONE
!-
  REAL, DIMENSION(iim_full, jjm_full), INTENT(IN) :: x_f
  REAL, DIMENSION(iim_zoom, jjm_zoom), INTENT(OUT) :: x_z

  INTEGER :: i,j
!-
  DO i=1,iim_zoom
     DO j=1,jjm_zoom
        x_z(i,j) = x_f(i_index(i),j_index(j))
     ENDDO
  ENDDO
!-
END SUBROUTINE forcing_zoom


!! ==============================================================================================================================\n
!! SUBROUTINE 	: forcing_vertical_ioipsl
!!
!>\BRIEF        
!!
!!\n DESCRIPTION : This subroutine explores the forcing file to decide what information is available
!!                 on the vertical coordinate.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S)	: 
!!
!_ ================================================================================================================================
SUBROUTINE forcing_vertical_ioipsl(force_id)

  INTEGER, INTENT(IN) :: force_id

  LOGICAL :: var_exists, vara_exists, varb_exists, varuv_exists
  LOGICAL :: foundvar
  LOGICAL :: levlegacy

  !
  !- Set all the defaults
  !
  zfixed=.FALSE.
  zsigma=.FALSE.
  zhybrid=.FALSE.
  zlevels=.FALSE.
  zheight=.FALSE.
  zsamelev_uv = .TRUE.
  levlegacy = .FALSE.
  !
  foundvar = .FALSE.
  !
  !- We have a forcing file to explore so let us see if we find any of the conventions
  !- which allow us to find the height of T,Q,U and V.
  !
  IF ( force_id > 0 ) THEN
     !
     ! Case for sigma levels
     !
     IF ( .NOT. foundvar ) THEN
        CALL flinquery_var(force_id, 'Sigma', var_exists)
        CALL flinquery_var(force_id, 'Sigma_uv', varuv_exists)
        IF ( var_exists ) THEN
           foundvar = .TRUE.
           zsigma = .TRUE.
           IF ( varuv_exists ) zsamelev_uv = .FALSE.
        ENDIF
     ENDIF
     !
     ! Case for Hybrid levels
     !
     IF ( .NOT. foundvar ) THEN
        CALL flinquery_var(force_id, 'HybSigA', vara_exists)
        IF ( vara_exists ) THEN
           CALL flinquery_var(force_id, 'HybSigB', varb_exists)
           IF ( varb_exists ) THEN
              var_exists=.TRUE.
           ELSE
              CALL ipslerr ( 3, 'forcing_vertical_ioipsl','Missing the B coefficient for', &
                   &         'Hybrid vertical levels for T and Q','stop readdim2')
           ENDIF
        ENDIF
        CALL flinquery_var(force_id, 'HybSigA_uv', vara_exists)
        IF ( vara_exists ) THEN
           CALL flinquery_var(force_id, 'HybSigB_uv', varb_exists)
           IF ( varb_exists ) THEN
              varuv_exists=.TRUE.
           ELSE
              CALL ipslerr ( 3, 'forcing_vertical_ioipsl','Missing the B coefficient for', &
                   &         'Hybrid vertical levels for U and V','stop readdim2')
           ENDIF
        ENDIF
        IF ( var_exists ) THEN
           foundvar = .TRUE.
           zhybrid = .TRUE.
           IF ( varuv_exists ) zsamelev_uv = .FALSE.
        ENDIF
     ENDIF
     !
     ! Case for levels (i.e. a 2d time varying field with the height in meters)
     !
     IF ( .NOT. foundvar ) THEN
        CALL flinquery_var(force_id, 'Levels', var_exists)
        CALL flinquery_var(force_id, 'Levels_uv', varuv_exists)
        IF ( var_exists ) THEN
           foundvar = .TRUE.
           zlevels = .TRUE.
           IF ( varuv_exists ) zsamelev_uv = .FALSE.
        ENDIF
     ENDIF
     !
     ! Case where a fixed height is provided in meters
     !
     IF ( .NOT. foundvar ) THEN
        CALL flinquery_var(force_id, 'Height_Lev1', var_exists)
        CALL flinquery_var(force_id, 'Height_Levuv', varuv_exists)
       IF ( var_exists ) THEN
           foundvar = .TRUE.
           zheight = .TRUE.
           IF ( varuv_exists ) zsamelev_uv = .FALSE.
        ENDIF
     ENDIF
     !
     ! Case where a fixed height is provided in meters in the lev variable
     !
     IF ( .NOT. foundvar ) THEN
        CALL flinquery_var(force_id, 'lev', var_exists)
        IF ( var_exists ) THEN
           foundvar = .TRUE.
           zheight = .TRUE.
           levlegacy = .TRUE.
        ENDIF
     ENDIF
     !
  ENDIF
  !
  ! We found forcing variables so we need to extract the values if we are dealing with constant values (i.e. all
  ! except the case zlevels
  !
  IF ( foundvar .AND. .NOT. zlevels ) THEN
     !
     IF ( zheight ) THEN
        !
        ! Constant height
        !
        IF ( levlegacy ) THEN
           CALL flinget (force_id,'lev',1, 1, 1, 1,  1, 1, zlev_fixed)
        ELSE
           CALL flinget (force_id,'Height_Lev1',1, 1, 1, 1,  1, 1, zlev_fixed)
           IF ( .NOT. zsamelev_uv ) THEN
              CALL flinget (force_id,'Height_Levuv',1, 1, 1, 1,  1, 1, zlevuv_fixed)
           ENDIF
        ENDIF
        IF (printlev_loc >= 1) WRITE(numout,*) "forcing_vertical_ioipsl : case ZLEV : Read from forcing file :", &
             zlev_fixed, zlevuv_fixed
        !
     ELSE IF ( zsigma .OR. zhybrid ) THEN
        !
        ! Sigma or hybrid levels
        !
        IF ( zsigma ) THEN
           CALL flinget (force_id,'Sigma',1, 1, 1, 1,  1, 1, zhybrid_b)
           zhybrid_a = zero
           IF ( .NOT. zsamelev_uv ) THEN
              CALL flinget (force_id,'Sigma_uv',1, 1, 1, 1,  1, 1, zhybriduv_b)
              zhybriduv_a = zero
           ENDIF
        ELSE
           CALL flinget (force_id,'HybSigB',1, 1, 1, 1,  1, 1, zhybrid_b)
           CALL flinget (force_id,'HybSigA',1, 1, 1, 1,  1, 1, zhybrid_a)
           IF ( .NOT. zsamelev_uv ) THEN
              CALL flinget (force_id,'HybSigB_uv',1, 1, 1, 1,  1, 1, zhybriduv_b)
              CALL flinget (force_id,'HybSigA_uv',1, 1, 1, 1,  1, 1, zhybriduv_a)
           ENDIF
        ENDIF
        IF (printlev_loc >= 1) WRITE(numout,*) "forcing_vertical_ioipsl : case Pressure coordinates : "
        IF (printlev_loc >= 1) WRITE(numout,*) "Read from forcing file :", zhybrid_b, zhybrid_a, zhybriduv_b, zhybriduv_a
     ELSE
        !
        ! Why are we here ???
        !
        CALL ipslerr ( 3, 'forcing_vertical_ioipsl','What is the option used to describe the height of', &
             &         'the atmospheric forcing ?','Please check your forcing file.')
     ENDIF
  ENDIF
  !
  !- We have no forcing file to explore or we did not find anything. So revert back to the run.def and
  !- read what has been specified by the user.
  !
  IF ( force_id < 0 .OR. .NOT. foundvar ) THEN
     !
     !-
     !Config Key   = HEIGHT_LEV1
     !Config Desc  = Height at which T and Q are given
     !Config Def   = 2.0
     !Config If    = offline mode
     !Config Help  = The atmospheric variables (temperature and specific
     !Config         humidity) are measured at a specific level.
     !Config         The height of this level is needed to compute
     !Config         correctly the turbulent transfer coefficients.
     !Config         Look at the description of the forcing
     !Config         DATA for the correct value.
     !Config Units = [m]
     !-
     zlev_fixed = 2.0
     CALL getin('HEIGHT_LEV1', zlev_fixed)

     !-
     !Config Key  = HEIGHT_LEVW
     !Config Desc = Height at which the wind is given
     !Config Def  = 10.0
     !Config If   = offline mode
     !Config Help = The height at which wind is needed to compute
     !Config        correctly the turbulent transfer coefficients.
     !Config Units= [m]
     !-
     zlevuv_fixed = 10.0
     CALL getin('HEIGHT_LEVW', zlevuv_fixed)

     zheight = .TRUE.

     IF ( ABS(zlevuv_fixed-zlev_fixed) > EPSILON(zlev_fixed)) THEN
        zsamelev_uv = .FALSE.
     ELSE
        zsamelev_uv = .TRUE.
     ENDIF

     CALL ipslerr ( 2, 'forcing_vertical_ioipsl','The height of the atmospheric forcing variables', &
          &         'was not found in the netCDF file.','Thus the values in run.def were used ... or their defaults.')
  ENDIF

END SUBROUTINE forcing_vertical_ioipsl


!! ==============================================================================================================================\n
!! SUBROUTINE 	: domain_size
!!
!>\BRIEF        
!!
!!\n DESCRIPTION : 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S)	: 
!!
!_ ================================================================================================================================
SUBROUTINE domain_size (limit_west, limit_east, limit_north, limit_south, &
     &  iim_f, jjm_f, lon, lat, iim, jjm, iind, jind)
  
  IMPLICIT NONE
  !
  ! ARGUMENTS
  !
  REAL, INTENT(inout)  :: limit_west,limit_east,limit_north,limit_south
  INTEGER, INTENT(in)  :: iim_f, jjm_f
  REAL, INTENT(in)     :: lon(iim_f, jjm_f), lat(iim_f, jjm_f)
  INTEGER, INTENT(out) :: iim,jjm
  INTEGER, INTENT(out) :: iind(iim_f), jind(jjm_f)
  !
  ! LOCAL
  !
  INTEGER :: i, j
  REAL :: lolo
  LOGICAL :: over_dateline = .FALSE.
  !
  !
  IF ( ( ABS(limit_east) .GT. 180. ) .OR. &
       ( ABS(limit_west) .GT. 180. ) ) THEN
     WRITE(numout,*) 'Limites Ouest, Est: ',limit_west,limit_east
     CALL ipslerr_p (3,'domain_size', &
 &        'Longitudes problem.','In run.def file :', &
 &        'limit_east > 180. or limit_west > 180.')
  ENDIF
  !
  IF ( limit_west .GT. limit_east ) over_dateline = .TRUE.
  !
  IF ( ( limit_south .LT. -90. ) .OR. &
       ( limit_north .GT. 90. ) .OR. &
       ( limit_south .GE. limit_north ) ) THEN
     WRITE(numout,*) 'Limites Nord, Sud: ',limit_north,limit_south
     CALL ipslerr_p (3,'domain_size', &
 &        'Latitudes problem.','In run.def file :', &
 &        'limit_south < -90. or limit_north > 90. or limit_south >= limit_north')
  ENDIF
  !
  ! Here we assume that the grid of the forcing data is regular. Else we would have
  ! to do more work to find the index table.
  !
  iim = 0
  DO i=1,iim_f
     !
     lolo =  lon(i,1)
     IF ( lon(i,1) .GT. 180. ) lolo =  lon(i,1) - 360.
     IF ( lon(i,1) .LT. -180. ) lolo =  lon(i,1) + 360.
     !
     IF (lon(i,1) < limit_west) iim_g_begin = i+1
     IF (lon(i,1) < limit_east) iim_g_end = i
     !
     IF ( over_dateline ) THEN
        IF ( lolo .LE. limit_west .OR. lolo .GE. limit_east ) THEN
           iim = iim + 1
           iind(iim) = i
        ENDIF
     ELSE
        IF ( lolo .GE. limit_west .AND. lolo .LE. limit_east ) THEN
           iim = iim + 1
           iind(iim) = i
        ENDIF
     ENDIF
     !
  ENDDO
  !
  jjm = 0
  DO j=1,jjm_f
     IF (lat(1,j) > limit_north) jjm_g_begin = j+1
     IF (lat(1,j) > limit_south) jjm_g_end = j
     !
     IF ( lat(1,j) .GE. limit_south .AND. lat(1,j) .LE. limit_north) THEN
        jjm = jjm + 1
        jind(jjm) = j
     ENDIF
  ENDDO

  IF (printlev_loc >= 1) WRITE(numout,*) 'Domain zoom size: iim, jjm = ', iim, jjm

  END SUBROUTINE domain_size


!! ==============================================================================================================================\n
!! SUBROUTINE 	: flinget_buffer
!!
!>\BRIEF        
!!
!!\n DESCRIPTION :  This subroutine is a wrap of flinget/IOIPSL. The arguments are the same. 
!!                  flinget_buffer will call flinget and buffer the forcing data localy in this subroutine. 
!!                  According to the variable NBUFF set in run.def, several time steps can be read at the same time 
!!                  from the forcing file. If NBUFF=0, the full forcing file is read. 
!!                  The output, data_full, from this subroutine is always only one time step of corresponding to itb.
!!                  itb must be equal to ite. 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S)	: 
!!
!_ ================================================================================================================================
  SUBROUTINE flinget_buffer(force_id, varname, iim_full, jjm_full, llm_full, ttm, itb, ite, data_full)
    
    !! Input arguments 
    INTEGER, INTENT(in)          :: force_id                      !! Id for forcing file
    CHARACTER(len=*), INTENT(in) :: varname                       !! Name of current variable to be read
    INTEGER, INTENT(in)          :: iim_full, jjm_full, llm_full  !! Horizontal and vertical domaine
    INTEGER, INTENT(in)          :: ttm                           !! Full lenght of forcing file
    INTEGER, INTENT(in)          :: itb, ite                      !! Time step to be read from forcing file. itb must be equal to ite

    !! Output argument
    REAL, DIMENSION(iim_full, jjm_full), INTENT(out) :: data_full !! Data for time step itb. 

    !! Define specific type to buffer data together with name and index
    TYPE buffer_type
       CHARACTER(len=20) :: name                   !! Name of variable in forcing file
       INTEGER :: istart                           !! Start index of current buffered data
       INTEGER :: iend                             !! End index of current buffered data
       REAL, ALLOCATABLE, DIMENSION(:,:,:) :: data !! Data read from forcing file for intervall [istart,iend]
    END TYPE buffer_type
    
    !! Local variables 
    INTEGER, PARAMETER :: maxvar=20                          !! Max number of variables to be buffered
    TYPE(buffer_type), DIMENSION(maxvar),SAVE :: data_buffer !! Containing all variables and the current buffered data
    INTEGER, SAVE   :: nbuff                                 !! Number of time steps to be buffered
    INTEGER, SAVE   :: lastindex=0                           !! Current number of variables stored in data_buffer
    INTEGER, SAVE   :: ttm0                                  !! Time lenght of forcing file, stored for test purpose
    LOGICAL, SAVE   :: first=.TRUE.                          !! First call to this subroutine
    INTEGER         :: index                                 !! Index in data_buffer for current variable
    INTEGER         :: i, ierr                               !! Loop and error variables
    INTEGER         :: tsend                                 !! Ending timestep
    INTEGER         :: new_buf_sz                            !! New buffer size
    INTEGER         :: nbuff_new                             !! Temporary variable for nbuff used in the end of the forcing file
    
    !! 1. Initialization
    IF (first) THEN
       data_buffer(:)%name='undef'
       ! Read NBUFF from run.def
       ! Note that getin_p is not used because this subroutine might be called only by master process

       !Config Key  = NBUFF
       !Config Desc = Number of time steps of data to buffer between each reading of the forcing file
       !Config If   = OFF_LINE
       !Config Help = The full simulation time length will be read if NBUFF equal 0. 
       !Config        NBUFF > 1 can be used for smaller regions or site simulations only. 
       !Config Def  = 15
       !Config Units= -

       nbuff=15
       CALL getin('NBUFF', nbuff)

       IF (nbuff == 0 .OR. nbuff >ttm) THEN
          ! Set nbuff as the full forcing file lenght
          nbuff=ttm
       ELSE IF (nbuff < 0) THEN
          ! Negativ nbuff not possible
          CALL ipslerr_p(3,'flinget_buffer','NBUFF must be a positiv number','Set NBUFF=0 for full simulation lenght','')
       END IF
       IF (printlev_loc >= 1) WRITE(numout,*)'flinget_buffer: NBUFF=',nbuff,' number of time step will be buffered'
       IF (printlev_loc >= 1) WRITE(numout,*)'flinget_buffer: Choose a lower value for NBUFF if problem with memory'
       
       ! Save dimensions to check following timesteps
       ! ttm is the full lenght of forcing file
       ttm0=ttm
       
       first=.FALSE.
    END IF

    !! 2. Coeherence tests on input arguments
    IF (ttm /= ttm0) THEN
       WRITE(numout,*)'Problem with ttm=',ttm,' ttm0=',ttm0
       CALL ipslerr_p(3,'flinget_buffer','Error with ttm and ttm0','','')
    END IF
    IF (itb /= ite) THEN
       WRITE(numout,*) 'There is a problem. Why is itb not equal ite ?'
       WRITE(numout,*) 'itb=',itb,' ite=',ite,' varname=',varname
       CALL ipslerr_p(3,'flinget_buffer','ite not equal itb','','')
    END IF


    !! 3. Find index for current variable
    index=0
    DO i=1, maxvar
       IF ( trim(varname) == data_buffer(i)%name ) THEN
          index=i
          CYCLE
       END IF
    END DO
    
    !! 4. Initialize and allocate if necesary the current variable
    IF ( index == 0 ) THEN
       ! The variable was not found
       ! This must be the first time for current variable
       index=lastindex+1
       lastindex=index
       IF (index > maxvar) CALL ipslerr_p(3,'flinget_buffer','to many variables','maxvar is too small','')
       
       ! Initialize the data_buffer for this index
       data_buffer(index)%name=trim(varname)
       ALLOCATE(data_buffer(index)%data(iim_full,jjm_full,nbuff),stat=ierr)
       IF (ierr /= 0) CALL ipslerr_p(3,'flinget_buffer','pb alloc data_buffer%data','for variable=',varname)
       data_buffer(index)%istart=0
       data_buffer(index)%iend=0
    END IF
    
    
    !! 5. Call flinget if current time step (itb) is outside the buffered intervall
    IF (( itb > data_buffer(index)%iend ) .OR. ( itb < data_buffer(index)%istart )) THEN
       ! itb is not in the time slice previously read or it is the first time to read
       ! Reading of forcing file will now be done
       ! First recalculate index to be read
       data_buffer(index)%istart = itb

       tsend = itb + nbuff - 1
       ! when NBUFF is not divisible by total forcing timesteps it requires some extra management 
       ! This avoid requesting more data at the end of the simulation 
       IF (tsend > ttm) tsend = ttm 
       new_buf_sz = tsend-itb+1

       ! Resize data buffer
       IF (SIZE(data_buffer(index)%data, DIM=3) .NE. new_buf_sz ) THEN
         DEALLOCATE(data_buffer(index)%data)
         ALLOCATE (data_buffer(index)%data(iim_full,jjm_full, new_buf_sz), stat=ierr )
         IF (ierr /= 0) CALL ipslerr_p(3,'flinget_buffer','pb alloc data_buffer%data','for variable=',varname)
       ENDIF
       data_buffer(index)%iend   = tsend 
       
       ! Check and correct if data_buffer(index)%iend is exceeding file size
       IF (data_buffer(index)%iend > ttm) THEN
          ! iend is exceeding the limit of the file. Change iend to the last time step in the file. 
          data_buffer(index)%iend = ttm

          ! Calculate a new smaller nbuff
          nbuff_new = ttm - itb + 1

          ! Resize data buffer
          DEALLOCATE(data_buffer(index)%data)
          ALLOCATE(data_buffer(index)%data(iim_full,jjm_full, nbuff_new), stat=ierr )
          IF (ierr /= 0) CALL ipslerr_p(3,'flinget_buffer','pb realloc data_buffer%data with new nbuff','','')
       END IF

!       WRITE(numout,*) 'Now do flinget for ',varname,', itb=',itb,', istart=',&
!            data_buffer(index)%istart,', iend=',data_buffer(index)%iend
       CALL flinget (force_id,varname, iim_full, jjm_full, llm_full, ttm, data_buffer(index)%istart, &
            data_buffer(index)%iend, data_buffer(index)%data(:,:,:))
    END IF
    
    !! 6. Initialize the output variable with data from buffered variable
    ! Find index for the time step corrsponding to itb in the time slice previously read from forcing file
    i=itb-data_buffer(index)%istart+1
    ! Initialize output variable
    data_full(:,:) = data_buffer(index)%data(:,:,i)
    
    
  END SUBROUTINE flinget_buffer

END MODULE readdim2
