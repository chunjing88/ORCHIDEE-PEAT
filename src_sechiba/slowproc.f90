! =================================================================================================================================
! MODULE       : slowproc
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         Groups the subroutines that: (1) initialize all variables used in 
!! slowproc_main, (2) prepare the restart file for the next simulation, (3) Update the 
!! vegetation cover if needed, and (4) handle all slow processes if the carbon
!! cycle is activated (call STOMATE) or update the vegetation properties (LAI and 
!! fractional cover) in the case of a run with only SECHIBA.
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	:
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_sechiba/slowproc.f90 $
!! $Date: 2018-05-30 15:40:04 +0200 (Wed, 30 May 2018) $
!! $Revision: 5268 $
!! \n
!_ ================================================================================================================================

MODULE slowproc

  USE defprec
  USE constantes 
  USE constantes_soil
  USE constantes_soil_var !! crop irrigation needs, xuhui
!!!qcj++ peatland
  USE constantes_var
  USE pft_parameters
  USE ioipsl
  USE xios_orchidee
  USE ioipsl_para
  USE sechiba_io_p
  USE interpol_help
  USE stomate
  USE stomate_data
  USE grid
  USE time, ONLY : dt_sechiba, dt_stomate, one_day, FirstTsYear, LastTsDay, LastTsYear, FirstTsMonth
  USE time, ONLY : year_start, month_start, day_start, sec_start
  USE time, ONLY : month_end, day_end
  USE mod_orchidee_para

  IMPLICIT NONE

  ! Private & public routines

  PRIVATE
  PUBLIC slowproc_main, slowproc_clear, slowproc_initialize, slowproc_finalize, slowproc_change_frac, slowproc_veget

  !
  ! variables used inside slowproc module : declaration and initialisation
  !
  REAL(r_std), SAVE                                  :: slope_default = 0.1
!$OMP THREADPRIVATE(slope_default)
  INTEGER(i_std) , SAVE                              :: veget_update        !! update frequency in years for landuse (nb of years)
!$OMP THREADPRIVATE(veget_update)
  !spitfire
  REAL(r_std), SAVE                                :: m_lightn_default = 0.02
!$OMP THREADPRIVATE(m_lightn_default)
  LOGICAL, SAVE                                   :: read_popdens
!$OMP THREADPRIVATE(read_popdens)
  LOGICAL, SAVE                                   :: read_humign
!$OMP THREADPRIVATE(read_humign)
  REAL(r_std), SAVE                                :: popdens_default = 1.
!$OMP THREADPRIVATE(popdens_default)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)     :: m_lightn
!$OMP THREADPRIVATE(m_lightn)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)     :: glccNetLCC            !! the land-cover-change (LCC) matrix in case a gross LCC is 
!$OMP THREADPRIVATE(glccNetLCC)
                                                                               !! used.
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)     :: glccSecondShift       
!$OMP THREADPRIVATE(glccSecondShift)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)     :: glccPrimaryShift
!$OMP THREADPRIVATE(glccPrimaryShift)
                                                                               
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)     :: harvest_matrix        
!$OMP THREADPRIVATE(harvest_matrix)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)     :: harvest_biomass        
!$OMP THREADPRIVATE(harvest_biomass)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)     :: bound_spa
!$OMP THREADPRIVATE(bound_spa)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)     :: newvegfrac            !! fraction of different MTCs that is used to guide on how to
                                                                               !! allocate the newly created MTC in gross LUC.
!$OMP THREADPRIVATE(newvegfrac)
                                                                               
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)     :: proxy_anidens
!$OMP THREADPRIVATE(proxy_anidens)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)     :: popd
!$OMP THREADPRIVATE(popd)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)     :: humign
!$OMP THREADPRIVATE(humign)
  REAL(r_std), SAVE                                :: m_ba_default = 0.
!$OMP THREADPRIVATE(m_ba_default)
  LOGICAL, SAVE                                   :: read_observed_ba
!$OMP THREADPRIVATE(read_observed_ba)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)   :: m_observed_ba
!$OMP THREADPRIVATE(m_observed_ba)

  REAL(r_std), SAVE                                :: m_cf_coarse_default = 0.
!$OMP THREADPRIVATE(m_cf_coarse_default)
  LOGICAL, SAVE                                   :: read_cf_coarse
!$OMP THREADPRIVATE(read_cf_coarse)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)   :: m_cf_coarse
!$OMP THREADPRIVATE(m_cf_coarse)

  REAL(r_std), SAVE                                :: m_cf_fine_default = 0.
!$OMP THREADPRIVATE(m_cf_fine_default)
  LOGICAL, SAVE                                   :: read_cf_fine
!$OMP THREADPRIVATE(read_cf_fine)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)   :: m_cf_fine
!$OMP THREADPRIVATE(m_cf_fine)

  REAL(r_std), SAVE                                :: m_ratio_default = 0.
!$OMP THREADPRIVATE(m_ratio_default)
  LOGICAL, SAVE                                   :: read_ratio
!$OMP THREADPRIVATE(read_ratio)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)   :: m_ratio
!$OMP THREADPRIVATE(m_ratio)

  REAL(r_std), SAVE                                :: m_ratio_flag_default = 0.
!$OMP THREADPRIVATE(m_ratio_flag_default)
  LOGICAL, SAVE                                   :: read_ratio_flag
!$OMP THREADPRIVATE(read_ratio_flag)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)   :: m_ratio_flag
!$OMP THREADPRIVATE(m_ratio_flag)

  !endspit
  !
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: clayfraction            !! Clayfraction (0-1, unitless)
!$OMP THREADPRIVATE(clayfraction)
  INTEGER, SAVE                                      :: printlev_loc        !! Local printlev in slowproc module 
!$OMP THREADPRIVATE(printlev_loc)  
REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: laimap              !! LAI map when the LAI is prescribed and not calculated by STOMATE
!$OMP THREADPRIVATE(laimap)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: soilclass_default
!$OMP THREADPRIVATE(soilclass_default)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: veget_max_new       !! New year fraction of vegetation type (0-1, unitless)
!$OMP THREADPRIVATE(veget_max_new)                                              
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: frac_nobio_new      !! New year fraction of ice+lakes+cities+... (0-1, unitless)
!$OMP THREADPRIVATE(frac_nobio_new)                                             
 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: frac_nobio_lastyear !! last year fraction of ice+lakes+cities+... (0-1, unitless)
!$OMP THREADPRIVATE(frac_nobio_lastyear)
  INTEGER(i_std), SAVE                               :: lcanop              !! canopy levels used for LAI
!$OMP THREADPRIVATE(lcanop)
  INTEGER(i_std) , SAVE                              :: veget_year          !! year for vegetation update
!$OMP THREADPRIVATE(veget_year)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: vegetnew_firstday          !! next year fraction of vegetation type (0-1, unitless)
!$OMP THREADPRIVATE(vegetnew_firstday)
!gmjc 15Feb2016 avoid grazing wet
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: fc_grazing
!$OMP THREADPRIVATE(fc_grazing)
!end gmjc
!!!qcj++ peatland
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: veget_cov_max_adjusted
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: peatPET_lastyear
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: peatPET_thisyear
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: precipitation_lastsummer
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: precipitation_thissummer
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: summerp_long
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: summerpet_long
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: peatC
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: peatC_ok
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: growth_day     !!! days with t2m_month greater than 0
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: GSL            !!! growing season length
  REAL(r_std), SAVE                                  :: PWT_lim   !!!Water balance threshold for peatland expanding
  REAL(r_std), SAVE                                  :: PC_lim    !!!Carbon threshold for peatland expanding
  REAL(r_std), SAVE                                  :: sat_gsl
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)      :: fpeat_map

!  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: t2m_lastyear
!  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: t2m_thisyear
CONTAINS

!! ================================================================================================================================
!! SUBROUTINE 	: slowproc_initialize
!!
!>\BRIEF         Initialize slowproc module and call initialization of stomate module
!!
!! DESCRIPTION : Allocate module variables, read from restart file or initialize with default values
!!               Call initialization of stomate module.
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_initialize (kjit,         kjpij,          kjpindex,       date0,          &
                                IndexLand,      indexveg,       lalo,           neighbours,     &
                                resolution,     contfrac,       soiltile,       reinf_slope,    &
                                t2m,                                                            &
                                deadleaf_cover, assim_param,    lai,            frac_age,       &
                                height,         veget,          frac_nobio,     njsc,           &
                                veget_max,      totfrac_nobio,  qsintmax,       rest_id,        &
                                rest_id_stom,   hist_id_stom,   tot_bare_soil,                  &
                                hist_id_stom_IPCC, co2_flux,    fco2_lu,        temp_growth,    &
                                soilc_total,   thawed_humidity, depth_organic_soil, heat_Zimov, &
                                f_rot_sech, altmax,fpeat,wtpmax,wtp_year) !!!qcj++ peatland
  

!! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                          :: kjit                !! Time step number
    INTEGER(i_std), INTENT(in)                          :: kjpij               !! Total size of the un-compressed grid
    INTEGER(i_std),INTENT(in)                           :: kjpindex            !! Domain size - terrestrial pixels only
    REAL(r_std),INTENT (in)                             :: date0               !! Initial date of what ???
    INTEGER(i_std),INTENT (in)                          :: rest_id             !! Restart file identifier
    INTEGER(i_std),INTENT (in)                          :: rest_id_stom        !! STOMATE's _Restart_ file identifier
    INTEGER(i_std),INTENT (in)                          :: hist_id_stom        !! STOMATE's _history_ file identifier
    INTEGER(i_std),INTENT(in)                           :: hist_id_stom_IPCC   !! STOMATE's IPCC _history_ file identifier
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: IndexLand           !! Indices of the points on the land map
    INTEGER(i_std),DIMENSION (kjpindex*nvm), INTENT (in):: indexveg            !! Indices of the points on the vegetation (3D map ???) 
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)     :: lalo                !! Geogr. coordinates (latitude,longitude) (degrees)
    INTEGER(i_std), DIMENSION (kjpindex,NbNeighb), INTENT(in):: neighbours     !! neighbouring grid points if land.
    REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)     :: resolution          !! size in x an y of the grid (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: contfrac            !! Fraction of continent in the grid (0-1, unitless)
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)        :: t2m                 !! 2 m air temperature (K)
    LOGICAL,DIMENSION(kjpindex), INTENT(in)             :: f_rot_sech          !! whether a grid is under rotation
    
!! 0.2 Output variables 
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(out)     :: co2_flux       !! CO2 flux per average ground area (gC m^{-2} dt_stomate^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: fco2_lu        !! CO2 flux from land-use (without forest management) (gC m^{-2} dt_stomate^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: temp_growth    !! Growth temperature (�C) - Is equal to t2m_month 
    INTEGER(i_std), DIMENSION(kjpindex), INTENT(out)       :: njsc           !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std), DIMENSION(kjpindex,ndeep,nvm),   INTENT (out)   :: heat_Zimov !! heating associated with decomposition
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)     :: lai            !! Leaf area index (m^2 m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)     :: height         !! height of vegetation (m)
    REAL(r_std),DIMENSION (kjpindex,nvm,nleafages), INTENT(out):: frac_age   !! Age efficacity from STOMATE for isoprene
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)     :: veget          !! Fraction of vegetation type including none biological fraction (unitless)
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (out)  :: frac_nobio     !! Fraction of ice, lakes, cities etc. in the mesh
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)     :: veget_max      !! Maximum fraction of vegetation type including none biological fraction (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)     :: altmax
!!!qcj++ peatland
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)          ::  fpeat
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)     ::  wtpmax
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)     ::  wtp_year

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: tot_bare_soil  !! Total evaporating bare soil fraction 
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: totfrac_nobio  !! Total fraction of ice+lakes+cities etc. in the mesh
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT(out)    :: soiltile       !! Fraction of each soil tile (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)          :: reinf_slope    !! slope coef for reinfiltration
    REAL(r_std),DIMENSION (kjpindex,nvm,npco2),INTENT (out):: assim_param    !! min+max+opt temperatures & vmax for photosynthesis (K, \mumol m^{-2} s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: deadleaf_cover !! Fraction of soil covered by dead leaves (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)     :: qsintmax       !! Maximum water storage on vegetation from interception (mm)

!! 0.3 Modified variables
    REAL(r_std), DIMENSION(kjpindex,ndeep,nvm), INTENT (inout) :: soilc_total  !! total soil carbon for use in thermal
    REAL(r_std), DIMENSION(kjpindex), INTENT (inout)         :: thawed_humidity!! specified humidity of thawed soil
    REAL(r_std), DIMENSION(kjpindex), INTENT (inout)         :: depth_organic_soil !! how deep is the organic soil?

!_ ================================================================================================================================

    !! 1. Perform the allocation of all variables, define some files and some flags. 
    !     Restart file read for Sechiba.
    CALL slowproc_init (kjit, kjpindex, IndexLand, lalo, neighbours, resolution, contfrac, &
         rest_id, lai, frac_age, veget, frac_nobio, totfrac_nobio, soiltile, reinf_slope, &
         veget_max, tot_bare_soil, njsc, &
         height, lcanop, veget_update, veget_year, f_rot_sech)
    

    !! 2. Define Time step in days for stomate
    dt_days = dt_stomate / one_day
    

    !! 3. check time step coherence between slow processes and fast processes
    IF ( dt_stomate .LT. dt_sechiba ) THEN
       WRITE(numout,*) 'slow_processes: time step smaller than forcing time step, dt_sechiba=',dt_sechiba,' dt_stomate=',dt_stomate
       CALL ipslerr_p(3,'slowproc_initialize','Coherence problem between dt_stomate and dt_sechiba',&
            'Time step smaller than forcing time step','')
    ENDIF
    
    !! 4. Call stomate to initialize all variables manadged in stomate,
    IF ( ok_stomate ) THEN

       CALL stomate_initialize &
            (kjit,           kjpij,                  kjpindex,                        &
             rest_id_stom,   hist_id_stom,           hist_id_stom_IPCC,               &
             indexLand,      lalo,                   neighbours,   resolution,        &
             contfrac,       totfrac_nobio,          clayfraction, t2m,               &
             lai,            veget,                  veget_max,                       &
             co2_flux,       fco2_lu,                                                 &
             deadleaf_cover, assim_param,            thawed_humidity, depth_organic_soil, &
             soilc_total,    heat_Zimov,             temp_growth,altmax,fpeat,wtpmax,wtp_year) !!qcj++ peatland
    ENDIF
    
    !! 5. Specific run without the carbon cycle (STOMATE not called): 
    !!     Need to initialize some variables that will be used in SECHIBA:
    !!     height, deadleaf_cover, assim_param, qsintmax.
    IF (.NOT. ok_stomate ) THEN
       CALL slowproc_derivvar (kjpindex, veget, lai, &
            qsintmax, deadleaf_cover, assim_param, height, temp_growth)
    ELSE
       qsintmax(:,:) = qsintcst * veget(:,:) * lai(:,:)
       qsintmax(:,1) = zero
    ENDIF
    
  END SUBROUTINE slowproc_initialize


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_main
!!
!>\BRIEF         Main routine that manage variable initialisation (slowproc_init), 
!! prepare the restart file with the slowproc variables, update the time variables 
!! for slow processes, and possibly update the vegetation cover, before calling 
!! STOMATE in the case of the carbon cycle activated or just update LAI (and possibly
!! the vegetation cover) for simulation with only SECHIBA   
!!
!!
!! DESCRIPTION  : (definitions, functional, design, flags): The subroutine manages 
!! diverses tasks:
!! (1) Initializing all variables of slowproc (first call)
!! (2) Preparation of the restart file for the next simulation with all prognostic variables
!! (3) Compute and update time variable for slow processes
!! (4) Update the vegetation cover if there is some land use change (only every years)
!! (5) Call STOMATE for the runs with the carbone cycle activated (ok_stomate) and compute the respiration
!!     and the net primary production
!! (6) Compute the LAI and possibly update the vegetation cover for run without STOMATE 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S):  ::co2_flux, ::fco2_lu, ::lai, ::height, ::veget, ::frac_nobio,  
!! ::veget_max, ::totfrac_nobio, ::soiltype, ::assim_param, ::deadleaf_cover, ::qsintmax,
!! and resp_maint, resp_hetero, resp_growth, npp that are calculated and stored
!! in stomate is activated.  
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : 
! \latexonly 
! \includegraphics(scale=0.5){SlowprocMainFlow.eps} !PP to be finalize!!)
! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_main (kjit, kjpij, kjpindex, date0, &
       IndexLand, indexveg, lalo, neighbours, resolution, contfrac, soiltile, &
       t2m, temp_sol, stempdiag, &
       humrel, shumdiag, litterhumdiag, precip_rain, precip_snow, &
       !spitfire
       wspeed, &
       !endspit
       gpp, &
       deadleaf_cover, &
       assim_param, &
       lai, frac_age, height, veget, frac_nobio, veget_max, totfrac_nobio, qsintmax, &
       rest_id, hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC, &
       co2_flux, fco2_lu, temp_growth,&
       swdown, evapot_corr, & ! crops, xuhui
       tdeep, hsdeep_long, snow, heat_Zimov, pb, &
       sfluxCH4_deep, sfluxCO2_deep, &
       thawed_humidity, depth_organic_soil, zz_deep, zz_coef_deep, &
       soilc_total,snowdz,snowrho, &
       tot_bare_soil, f_rot_sech, rot_cmd, &
!gmjc
       tmc_topgrass, humcste_use,altmax, &
! end gmjc
!!!qcj++ peatland
       wtp,fwet_new,fpeat,mc_peat_above,peat_PET,liqwt_ratio,shumdiag_peat,&
       wtpmax,veget_cov_max_new_peatdgvm,&
       wtp_year,wettile_dgvm,tile_name_dgvm)
  
!! INTERFACE DESCRIPTION

!! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                          :: kjit                !! Time step number
    INTEGER(i_std), INTENT(in)                          :: kjpij               !! Total size of the un-compressed grid
    INTEGER(i_std),INTENT(in)                           :: kjpindex            !! Domain size - terrestrial pixels only
    REAL(r_std),INTENT (in)                             :: date0               !! Initial date of what ???
    INTEGER(i_std),INTENT (in)                          :: rest_id,hist_id     !! _Restart_ file and _history_ file identifier
    INTEGER(i_std),INTENT (in)                          :: hist2_id            !! _history_ file 2 identifier
    INTEGER(i_std),INTENT (in)                          :: rest_id_stom        !! STOMATE's _Restart_ file identifier
    INTEGER(i_std),INTENT (in)                          :: hist_id_stom        !! STOMATE's _history_ file identifier
    INTEGER(i_std),INTENT(in)                           :: hist_id_stom_IPCC   !! STOMATE's IPCC _history_ file identifier
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: IndexLand           !! Indices of the points on the land map
    INTEGER(i_std),DIMENSION (kjpindex*nvm), INTENT (in):: indexveg            !! Indices of the points on the vegetation (3D map ???) 
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)     :: lalo                !! Geogr. coordinates (latitude,longitude) (degrees)
    INTEGER(i_std), DIMENSION (kjpindex,NbNeighb), INTENT(in)  :: neighbours   !! neighbouring grid points if land
    REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)     :: resolution          !! size in x an y of the grid (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: contfrac            !! Fraction of continent in the grid (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in)  :: humrel              !! Relative humidity ("moisture stress") (0-1, unitless)
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)        :: t2m                 !! 2 m air temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: temp_sol            !! Surface temperature (K)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)  :: stempdiag           !! Soil temperature (K)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)  :: shumdiag            !! Relative soil moisture (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: litterhumdiag       !! Litter humidity  (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: precip_rain         !! Rain precipitation (mm dt_stomate^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: precip_snow         !! Snow precipitation (mm dt_stomate^{-1})
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)    :: gpp                 !! GPP of total ground area (gC m^{-2} time step^{-1}). 
                                                                               !! Calculated in sechiba, account for vegetation cover and 
                                                                               !! effective time step to obtain gpp_d
!!!qcj++ peatland
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)               :: fwet_new
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)               :: liqwt_ratio
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in)          :: wtp
    REAL(r_std), DIMENSION (kjpindex), INTENT (inout)            :: fpeat
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in)               :: mc_peat_above
!    LOGICAL                                                      :: FirstTsMonth
!    REAL(r_std)                                                  :: sum_veg    ! sum of vegets
!    REAL(r_std)                 :: sum_nobio                   ! sum of nobios
!    REAL(r_std)                                                  :: sumvpeat_old, sumvpeat ! last an new sum of peatland vegetation
!    REAL(r_std)                                                  :: rapport      ! (S-B) / (S-A)
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)               :: peat_PET 
    LOGICAL,SAVE                                                 :: Judge_PC
    REAL(r_std), DIMENSION(kjpindex,nslm,nvm),   INTENT (in)    :: shumdiag_peat
!    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(inout)     :: veget_max_adjusted !! Maximum fraction of vegetation adjusted with peatland
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(inout)        :: veget_cov_max_new_peatdgvm
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (inout)      :: wtpmax
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (inout)      :: wtp_year
    LOGICAL,DIMENSION (nstm), INTENT (in)                 :: wettile_dgvm
    CHARACTER(LEN=25), DIMENSION(nstm),INTENT (in)        :: tile_name_dgvm


!!!!! crops
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: swdown            !!downward shortwave radiation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: evapot_corr      !!potential evaportranspiration (mm)
    LOGICAL, DIMENSION(kjpindex), INTENT(out)           :: f_rot_sech
    INTEGER(i_std), DIMENSION(kjpindex,rot_cmd_max), INTENT(out) :: rot_cmd
!!!!! crops, xuhui 
    REAL(r_std), DIMENSION(kjpindex,ndeep,nvm),   INTENT (in)    :: tdeep      !! deep temperature profile
    REAL(r_std), DIMENSION(kjpindex,ndeep,nvm),   INTENT (in)    :: hsdeep_long!! deep long term soil humidity profile
    REAL(r_std), DIMENSION(kjpindex),         INTENT (in)        :: snow       !! Snow mass [Kg/m^2]
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)               :: pb         !! Lowest level pressure
    REAL(r_std), DIMENSION(ndeep),   INTENT (in)                 :: zz_deep    !! deep vertical profile
    REAL(r_std), DIMENSION(ndeep),   INTENT (in)                 :: zz_coef_deep!! deep vertical profile   
    REAL(r_std), DIMENSION(kjpindex,nsnow),INTENT(in)            :: snowdz     !! snow depth for each layer
    REAL(r_std), DIMENSION(kjpindex,nsnow),INTENT(in)            :: snowrho    !! snow density for each layer
    !spitfire
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: wspeed              !!Wind speed (m/s)
    !endspit  
!gmjc top 5 layer grassland soil moisture for grazing
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: tmc_topgrass
!end gmjc 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(in)       :: humcste_use
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(inout)       :: altmax
!! 0.2 Output variables 
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(out)  :: co2_flux            !! CO2 flux per average ground area (gC m^{-2} dt_stomate^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: fco2_lu             !! CO2 flux from land-use (without forest management) (gC m^{-2} dt_stomate^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: temp_growth         !! Growth temperature (�C) - Is equal to t2m_month 
    REAL(r_std), DIMENSION(kjpindex,ndeep,nvm),   INTENT (out)   :: heat_Zimov !! heating associated with decomposition
    REAL(r_std), DIMENSION(kjpindex),     INTENT (out)           :: sfluxCH4_deep      !! surface flux of CH4 to atmosphere from permafrost
    REAL(r_std), DIMENSION(kjpindex),     INTENT (out)           :: sfluxCO2_deep      !! surface flux of CO2 to atmosphere from permafrost
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)      :: tot_bare_soil       !! Total evaporating bare soil fraction in the mesh
    INTEGER(i_std), DIMENSION (kjpindex)                   :: sat_duration
!! 0.3 Modified variables
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)     :: lai            !! Leaf area index (m^2 m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)     :: height         !! height of vegetation (m)
    REAL(r_std),DIMENSION (kjpindex,nvm,nleafages), INTENT(inout):: frac_age   !! Age efficacity from STOMATE for isoprene
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)     :: veget          !! Fraction of vegetation type including none biological fractionin the mesh (unitless)
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (inout)  :: frac_nobio     !! Fraction of ice, lakes, cities etc. in the mesh
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)     :: veget_max      !! Maximum fraction of vegetation type in the mesh (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)         :: totfrac_nobio  !! Total fraction of ice+lakes+cities etc. in the mesh
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT(inout)    :: soiltile       !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm,npco2),INTENT (inout):: assim_param    !! min+max+opt temperatures & vmax for photosynthesis (K, \mumol m^{-2} s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)         :: deadleaf_cover !! Fraction of soil covered by dead leaves (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)     :: qsintmax       !! Maximum water storage on vegetation from interception (mm)
    REAL(r_std), DIMENSION(kjpindex),   INTENT (inout)           :: thawed_humidity    !! specified humidity of thawed soil
    REAL(r_std), DIMENSION(kjpindex),   INTENT (inout)           :: depth_organic_soil !! how deep is the organic soil?
    REAL(r_std), DIMENSION(kjpindex,ndeep,nvm),   INTENT (inout) :: soilc_total        !! total soil carbon for use in thermal

!! 0.4 Local variables
    INTEGER(i_std)                                     :: i, j, jv, ji            !! indices 
    REAL(r_std), DIMENSION(kjpindex,nvm)               :: resp_maint           !! Maitanance component of autotrophic respiration in (gC m^{-2} dt_stomate^{-1})
    REAL(r_std), DIMENSION(kjpindex,nvm)               :: resp_hetero          !! heterotrophic resp. (gC/(m**2 of total ground)/time step)
    REAL(r_std), DIMENSION(kjpindex,nvm)               :: resp_growth          !! Growth component of autotrophic respiration in gC m^{-2} dt_stomate^{-1})
    REAL(r_std), DIMENSION(kjpindex,nvm)               :: npp                  !! Net Ecosystem Exchange (gC/(m**2 of total ground)/time step)
    REAL(r_std),DIMENSION (kjpindex,nvm)               :: veget_nextyear       !! Temporary variable for new veget_max
    REAL(r_std),DIMENSION (kjpindex,nnobio)            :: frac_nobio_nextyear  !! Temporary variable for new frac_nobio
    REAL(r_std),DIMENSION (kjpindex)                   :: totfrac_nobio_lastyear !! Total fraction for the previous year
    REAL(r_std),DIMENSION (kjpindex)                   :: totfrac_nobio_new    !! Total fraction for the next year
    !spitfire
    REAL(r_std),DIMENSION (kjpindex)                    :: lightn            !!lightning (flashes/km2/day)
    REAL(r_std),DIMENSION (kjpindex)                    :: animal_density 
    REAL(r_std),DIMENSION (kjpindex)                    :: observed_ba       !!observed burned area (ha)
    INTEGER(i_std)                                      :: yy, mm, dd
    REAL(r_std)                                         :: ss
    REAL(r_std),DIMENSION (kjpindex)                    :: cf_coarse       !!observed burned area (ha)
    REAL(r_std),DIMENSION (kjpindex)                    :: cf_fine       !!observed burned area (ha)
    REAL(r_std),DIMENSION (kjpindex)                    :: ratio       !!observed burned area (ha)
    REAL(r_std),DIMENSION (kjpindex)                    :: ratio_flag       !!observed burned area (ha)
    REAL(r_std)                                         :: cropshare_old, cropshare_new
    !endspit
    INTEGER(i_std)                                     :: ivma,ivm,jvm         !! Indices
    INTEGER(i_std) , SAVE                              :: veget_year_tmp           !! year for landuse

    REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE      :: manage  !! temporary matrix for rotation input
    INTEGER(i_std)                                :: yrlen
    CHARACTER(LEN=30)                             :: strManage, strVar
    LOGICAL                                       :: is_update_rotation_cycle  = .FALSE. !! enable update of rotation cycle for CROP module

    REAL(r_std)                                         :: in_julian
!_ ================================================================================================================================

    !! 1. Compute and update all variables linked to the date and time
    IF (printlev_loc>=5) WRITE(numout,*) 'Entering slowproc_main' !  year_start, month_start, day_start, sec_start='

!         year_start, month_start,day_start,sec_start   
    
!    IF ( sec_start == dt_sechiba .AND. month==1 .AND. day==1 ) THEN
       ! The current time step is the first sechiba timestep of a new year
!       IF (printlev_loc>=4) WRITE(numout,*) "This is a new day and a new year: month, day, sec_start=", month, day, sec_start
!       FirstTsYear=.TRUE.
!    ELSE
!       FirstTsYear=.FALSE.
!    END IF
    
!    IF ( sec_start == 0 ) THEN
       ! The current time step is the last sechiba time step on a day
!       LastTsDay=.TRUE.
       
!       IF ( month == 1 .AND. day == 1 ) THEN
          ! The current time step is the last sechiba time step on a year
          ! JG : note that month=1, day=1, sec=O is the last day of the year. 
          !      This is due to a problem before the first call to slowproc. 
          !      slowproc_main enters the first time on the 2nd time step (1800s)
!          LastTsYear = .TRUE.
!          IF (printlev_loc>=4) WRITE(numout,*) "This is the last sechiba time step of a year, LastTsYear is activated"
!       ELSE
!          LastTsYear = .FALSE.
!       END IF
!    ELSE
!       LastTsDay = .FALSE.
!       LastTsYear = .FALSE.
!    END IF

! the above time calculation has been moved to time.f90
! but we need it for peatland calculation by qcj
    !!!qcj ++ peatland
    IF (FirstTsYear) THEN
           growth_day(:)=0
    END IF
    !!! qcj++ peatland, only when peat C density exceeds a certain amount, can the
    !peatland expands
    Judge_PC = .FALSE.
    CALL getin_p('JUDGE_PC',Judge_PC)

    IF (sec_start == dt_sechiba .AND. month_start==1 .AND. day_start==2) THEN
       DO ji=1, kjpindex
          IF ( .NOT. dynpeat_PC) THEN
            peatC_ok(ji) = 1.
          ELSE
            IF (Judge_PC) THEN
              IF (peatC(ji) .GE. PC_lim ) THEN
                 peatC_ok(ji) = 1.
              ELSE
                 peatC_ok(ji) = 0.
              ENDIF
            ENDIF
          ENDIF
       ENDDO
    ENDIF
!!! qcj ++ peatland
    IF ( sec_start == 0 ) THEN  !!!the last time step of the day
      DO ji = 1,kjpindex
         IF (temp_growth(ji) .GT. 5.) THEN
            growth_day(ji)= growth_day(ji)+1
         ELSE
            growth_day(ji)= growth_day(ji)+0
         ENDIF
      ENDDO
    ENDIF

!!!qcj++ peatland
    IF (dyn_peat) THEN
!!!Water balance, annual evapotranspiration*
       IF (dyn_peat .AND. dynpeat_PWT )  THEN
          IF (month_start <= 9 .AND. month_start >= 5) THEN  
              precipitation_thissummer(:) = precipitation_thissummer(:)+precip_rain(:)+precip_snow(:)
              peatPET_thisyear(:)=peatPET_thisyear(:)+peat_PET(:)
          ENDIF
       ENDIF


       IF (dyn_peat .AND. dynpeat_PWT .AND. LastTsYear ) THEN
          peatPET_lastyear(:)=peatPET_thisyear(:)
          peatPET_thisyear(:)=zero
          precipitation_lastsummer(:)=precipitation_thissummer(:)
          precipitation_thissummer(:)= zero
          GSL(:) = growth_day(:)/365.
          GSL(:)=GSL(:)*sat_gsl
       ENDIF
       sat_duration(:)= INT(30*12*GSL(:))  !!!Length of growing season
       sat_duration(:) = MIN (360,MAX(1,360-sat_duration(:)))
    ENDIF
!!!update fraction of peatland vegetion every month, peatland vegetation
!fraction = fpeat
    IF (dyn_peat) THEN
       IF (FirstTsMonth) THEN
          update_peatfrac= .TRUE.
       ENDIF
    ENDIF

    IF (dyn_peat .AND. dynpeat_PWT .AND. LastTsYear) THEN
!!!longterm precipitation and PET, 30 years mean
       IF (day_start==365) THEN
          summerp_long(:)= precipitation_lastsummer(:)
          summerpet_long(:)= peatPET_lastyear(:)
       ENDIF
       summerp_long(:)=( summerp_long(:) * ( 30. - 1. ) + &
                        precipitation_lastsummer(:) * 1. ) / 30.
       summerpet_long(:)=( summerpet_long(:) * ( 30. - 1. ) + &
                         peatPET_lastyear(:) * 1. ) / 30.
    ENDIF




    !! 2. Activate slow processes if it is the end of the day
    IF ( LastTsDay ) THEN
       ! 3.2.2 Activate slow processes in the end of the day
       do_slow = .TRUE.
       
       ! 3.2.3 Count the number of days 
       days_since_beg = days_since_beg + 1
       IF (printlev_loc>=4) WRITE(numout,*) "New days_since_beg : ",days_since_beg
    ELSE
       do_slow = .FALSE.
    ENDIF

    !! 3. Update the vegetation if it is time to do so.
    !!    This is done at the first sechiba time step on a new year and only every "veget_update" years. 
    !!    veget_update correspond to a number of years between each vegetation updates.
    !!    Nothing is done if veget_update=0.
    !!    Update of the vegetation map can not be done if map_pft_format=false.

    IF ( (map_pft_format) .AND. (veget_update .GT. 0) ) THEN

       IF ( days_since_beg == 1 ) THEN !!! This is in case the simulation did not start from the beginning of the year (noted xuhui)
          !
          veget_year_tmp = veget_year + 1

          ! Update of the vegetation cover with Land Use only if 
          ! the current year match the requested condition (a multiple of
          ! "veget_update")
          IF ( ( MOD(veget_year_tmp - veget_year_orig, veget_update) == 0) .AND. (.NOT. ok_dgvm_peat) ) THEN  !!qcj++ peatland
          
             WRITE(numout,*)  'We read the new vegetmax map for year =' , veget_year_tmp
             
             ! Call the routine to update the vegetation (output is vegetnew_firstday)
             CALL slowproc_readvegetmax(kjpindex, lalo, neighbours, resolution, contfrac, &
               &               veget_max, veget_nextyear, frac_nobio_nextyear, veget_year, .FALSE.)
             !!!! veget_nextyear is saved for used next year, but since it is
             !now input veget_map in the FirstTsYear, it is useless... (noted xuhui)
            !print *,'veget_nextyear after reading the new map firsttime zz',veget_nextyear
             !
          ENDIF
          !
       ENDIF
    ENDIF 

    ! Update vegetation and fraction and save old values
    frac_nobio_lastyear = frac_nobio

    IF ( map_pft_format .AND. (veget_update > 0) .AND. FirstTsYear ) THEN
       veget_year = veget_year + 1

       ! Update of the vegetation cover with Land Use only if 
       ! the current year match the requested condition (a multiple of "veget_update")
       IF ( MOD(veget_year - veget_year_orig, veget_update) == 0 ) THEN
          IF (printlev_loc>=1) WRITE(numout,*)  'We are updating the vegetation map for year =' , veget_year
          
          ! Read the new the vegetation from file. Output is veget_nextyear and frac_nobio_nextyear.
          IF ( .NOT. ok_dgvm .AND. (dyn_peat .OR. agri_peat)) THEN
              CALL slowproc_readvegetmax(kjpindex, lalo, neighbours, resolution, contfrac, &
               veget_max, veget_nextyear, frac_nobio_nextyear, veget_year, .TRUE.)
          ELSE
              CALL slowproc_readvegetmax(kjpindex, lalo, neighbours, resolution, contfrac, &
               veget_max, veget_nextyear, frac_nobio_nextyear, veget_year, .FALSE.)  !.FALSE. 
          ENDIF      
          IF (.NOT. use_age_class) THEN

            IF (ok_rotate) THEN 
            ! when rotation is activated, the conversion among different croplands will
            ! not follow the vegetation map but follow the rotation commands.
            ! Vegetation map is used as change of cropland share in the grid.
            ! This change is shared proportionally to contemporary croplands.
            ! xuhui 20160503
                DO ji = 1,kjpindex
                    cropshare_old = SUM(veget_max(ji,:),MASK=ok_LAIdev(:))
                    cropshare_new = SUM(veget_nextyear(ji,:),MASK=ok_LAIdev(:))
                    IF (printlev>=4) THEN
                        WRITE(numout,*) 'ji, cropshare_old, cropshare_new',ji, cropshare_old, cropshare_new
                    ENDIF
                    IF (cropshare_old .LT. min_sechiba) THEN
                    ! special case 1, no croplands previously existed
                        IF (printlev>=4) WRITE(numout,*) 'new croplands included, no change to veget_nextyear'
                    ELSEIF (cropshare_new .LT. min_sechiba) THEN
                    !special case 2, all croplands are removed
                        IF (printlev>=4) WRITE(numout,*) 'all croplands are killed, a case has not been tested'
                        !xuhui: need to consider how could it be well treated
                    ELSE
                        DO jv = 2,nvm
                            IF (ok_LAIdev(jv)) THEN
                                veget_nextyear(ji,jv) = veget_max(ji,jv) * cropshare_new / cropshare_old
                            ENDIF
                        ENDDO
                    ENDIF
                ENDDO
            ENDIF
            veget_max_new       = veget_nextyear
            frac_nobio_new      = frac_nobio_nextyear
          ! [chaoyue]
          ! the new veget_max will be calculated in the gross land use change
          ! module when use_age_class is True. Here we will just make the 
          ! veget_max_new and frac_nobio_new being the old values.

!!            !!!! xuhui: so the veget_max no longer updates here, so re-consider
!!            !how to maintain crop fraction....
!!            !!! should maintain veget_max_new so that when slowproc_change_frac
!!            !applies, it works fine
!!              IF (ok_rotate) THEN
!!                  DO jv = 1,nvm
!!                    !!!! proportionally keeping the bare soil and croplands as it
!!                    !was previously in order to maintain the rotation cycle
!!                  ENDDO
!!              ENDIF
          ELSE !!!! if use_age_class, readvegetmax is no longer used
            veget_max_new      = veget_max
            frac_nobio_new     = frac_nobio_lastyear
          ENDIF
          
          ! Verification and correction on veget_max, calculation of veget and soiltile.
          ! [chaoyue] this call in the trunk equivalent is lacking, it's done 
          ! in sechiba.f90 by calling slowproc_change_frac
          CALL slowproc_veget (kjpindex, f_rot_sech, lai, frac_nobio, totfrac_nobio, veget_max, veget, soiltile)
          
          ! Set the flag do_now_stomate_lcchange to activate stomate_lcchange.
          ! This flag will be kept to true until stomate_lcchange has been done. 
          ! The variable totfrac_nobio_lastyear will only be used in stomate when this flag is activated
          do_now_stomate_lcchange=.TRUE.
          IF ( .NOT. ok_stomate ) THEN
             ! Special case if stomate is not activated : set the variable done_stomate_lcchange=true 
             ! so that the subroutine slowproc_change_frac will be called in the end of sechiba_main.
             done_stomate_lcchange=.TRUE.
          END IF
       
    !print *,'veget_max after reading the new map zz',veget_max
    !print *,'frac_nobio_lastyear after reading the new map zz',frac_nobio_lastyear
    !print *,'totfrac_nobio after reading the new map zz',totfrac_nobio
       ENDIF
    ENDIF

    !WRITE(numout,*),'do_now_stomate_lcchange in slowproc,',do_now_stomate_lcchange
    !WRITE(numout,*),'veget_year in slowproc,',veget_year
    IF ( (rotation_update .GT. 0) .AND. FirstTsYear ) THEN 
    ! rotation_update is necessarily  zero, if not ok_rotate
        IF ( MOD(veget_year - veget_year_orig, rotation_update) == 0 ) THEN ! update rotation cycle
            IF (printlev_loc>=1) THEN
                WRITE(numout,*) 'xuhui: updating rotation system at year ', veget_year
                WRITE(numout,*) 'rotation_update, ', rotation_update
            ENDIF
            is_update_rotation_cycle  = .TRUE. ! to stomate_main
        ENDIF ! start rotation update
    ENDIF ! FirstTsOfYear

       IF (dyn_peat .AND. update_peatfrac) THEN
          DO ji = 1,kjpindex
!update peat fraction, to be fpeat
             DO jv = 1,nvm
                IF (is_peat(jv)) THEN
                   veget_max_new(ji,jv)= zero
                !!!need positive summer water balance for peatland expanding
                   IF (dynpeat_PWT) THEN
                     ! meant2m(ji) = t2m_lastyear(ji)/365*dt_sechiba/one_day-273.15
                      !!!peatland expand
                      
                      IF (fpeat(ji) .GT. veget_max(ji,jv)) THEN
                         IF (fpeat(ji) .GT. min_vegfrac) THEN   !ini_peat=1.E-6_r_std  !min_stomate
                         !!! positive summer water balance is needed for peatland expansion
                           IF (precipitation_lastsummer(ji)-peatPET_lastyear(ji) .GE. PWT_lim)   THEN  !-60, 0.
                         !! longterm summer water balance  
                         !  IF (summerp_long (ji) - summerpet_long(ji) .GE. PWT_lim) THEN
                              IF (veget_max(ji,jv) .LE. zero) THEN
                                 !!! first initiation of peatland is not limited by carbon  
                                 !!! so that we can have a seed of peatland development  
                                 veget_max_new(ji,jv)=fpeat(ji)

                              ELSE
                                 !!! after first initiation, an specific amount of carbon is needed before expanding 
                                 IF (peatC_ok(ji) .GT. 0.) THEN
                                     veget_max_new(ji,jv)=fpeat(ji)
                                 ELSE
                                     veget_max_new(ji,jv)=veget_max(ji,jv)
                                 ENDIF
                              ENDIF
                           ELSE
                             veget_max_new(ji,jv)=veget_max(ji,jv)
                           ENDIF
                         ELSE
                      !!!if fpeat near zero, pealtand won't initiate
                            veget_max_new(ji,jv)=veget_max(ji,jv)
                         ENDIF                            
                      ELSE
                      !!!peatland shrink, no limitation
                          IF (fpeat(ji) .GT. min_stomate) THEN
                             veget_max_new(ji,jv)=fpeat(ji)
                          ELSE
                             veget_max_new(ji,jv)= zero
                          ENDIF
                      ENDIF

                   ELSE
                !!! water balance and peat C are not important
                      IF (fpeat(ji) .GT. veget_max(ji,jv)) THEN
!!!peatland expand
                         IF (fpeat(ji) .GT. ini_peat) THEN !ini_peat=1.E-6_r_std  !min_stomate
                             veget_max_new(ji,jv)=fpeat(ji)
                         ELSE
                             veget_max_new(ji,jv)=veget_max(ji,jv)
                         ENDIF
                      ELSE
!!!peatland contract
                         IF (fpeat(ji) .GT. min_stomate) THEN
                             veget_max_new(ji,jv)=fpeat(ji)
                         ELSE
                             veget_max_new(ji,jv)= zero
                         ENDIF
                      ENDIF
                   ENDIF 
                ENDIF
             ENDDO !nvm
          ENDDO !kjpindex
       ENDIF  !update_peatfrac

    ! 5 call STOMATE, either because we want to keep track of
    !   long-term variables (WATCHOUT case) or just because STOMATE is
    !   activated

    !spitfire
     in_julian = itau2date(kjit, date0, dt_sechiba)
     CALL ju2ymds(in_julian, yy, mm, dd, ss) 
     lightn(:)=m_lightn(:,mm)
     animal_density(:)=proxy_anidens(:,mm)
     observed_ba(:)=m_observed_ba(:,mm)
     cf_coarse(:)=m_cf_coarse(:,mm)
     cf_fine(:)=m_cf_fine(:,mm)
     ratio(:)=m_ratio(:,mm)
     ratio_flag(:)=m_ratio_flag(:,mm)
    !endspit

    IF (FirstTsYear .AND. ok_dgvm_peat) THEN
       CALL slowproc_read_wetland(kjpindex,lalo,resolution,wettile_dgvm,tile_name_dgvm)
    ENDIF

    IF ( ok_stomate ) THEN
       ! Caluclate totfrac_nobio_lastyear
       totfrac_nobio_lastyear(:) = zero
       DO jv = 1, nnobio
          totfrac_nobio_lastyear(:) = totfrac_nobio_lastyear(:) + frac_nobio_lastyear(:,jv)
       ENDDO

       ! Caluclate totfrac_nobio_new only for the case when the land use map has been read previously
       IF (do_now_stomate_lcchange) THEN                               
          IF (.NOT. use_age_class) THEN         
             totfrac_nobio_new(:) = zero                                           
             DO jv = 1, nnobio                                                     
                totfrac_nobio_new(:) = totfrac_nobio_new(:) + frac_nobio_new(:,jv) 
             ENDDO                                                                 
          ELSE
             totfrac_nobio_new = totfrac_nobio_lastyear
          ENDIF
       ELSE                                                                     
          totfrac_nobio_new(:) = zero                                           
       END IF                                                                   

       !! 4.1 Call stomate main routine that will call all c-cycle routines       !
       CALL stomate_main (kjit, kjpij, kjpindex, &
            IndexLand, lalo, neighbours, resolution, contfrac, totfrac_nobio_lastyear, clayfraction, &
            t2m, temp_sol, stempdiag, &
            humrel, shumdiag, litterhumdiag, precip_rain, precip_snow, &
            !spitfire
            wspeed, lightn, popd, read_observed_ba, observed_ba, humign, & 
            read_cf_fine,cf_fine,read_cf_coarse,cf_coarse,read_ratio_flag,ratio_flag,read_ratio,ratio,&
            !endspit
            gpp, &
            deadleaf_cover, &
            assim_param, &
            lai, frac_age, height, veget, veget_max, &
            veget_max_new,vegetnew_firstday, totfrac_nobio_new, &
            glccNetLCC,glccSecondShift,glccPrimaryShift, &
            harvest_matrix, harvest_biomass,bound_spa, newvegfrac, &
            hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC, &
            co2_flux, fco2_lu, resp_maint,resp_hetero,resp_growth,temp_growth,&
            swdown, evapot_corr, is_update_rotation_cycle, &   !!! xuhui added for crops
            tdeep, hsdeep_long, snow, heat_Zimov, pb, &
            sfluxCH4_deep, sfluxCO2_deep, &
            thawed_humidity, depth_organic_soil, zz_deep, &
            zz_coef_deep, soilc_total,snowdz,snowrho, &
            LastTsYear, f_rot_sech, rot_cmd, &
!gmjc top 5 layer grassland soil moisture for grazing
            tmc_topgrass,fc_grazing,humcste_use,altmax, &
!end gmjc
!!!qcj++ peatland
            wtp,fwet_new,fpeat,shumdiag_peat,mc_peat_above, &
            peatC,sat_duration,liqwt_ratio,veget_cov_max_adjusted,&
            wtpmax,fpeat_map,veget_cov_max_new_peatdgvm,&
            wtp_year,wettile_dgvm)


       ! [chaoyue] in case of use_age_class and gross land use change, veget_max
       ! will be updated in the land use change module rather than using veget_max_new,
       ! thus we should pass the updated veget_max into veget_max_new
       IF (done_stomate_lcchange .AND. use_age_class) THEN
         veget_max_new = veget_max

         IF (min_vegfrac .GT. min_stomate) THEN
            WRITE(numout,*) "the vlaue of min_vegfrac is too big: ",min_vegfrac
            WRITE(numout,*) "When calling slowproc_veget, a too small vlaue of"
            WRITE(numout,*) "min_vegfrac will leads to readjustment of veget_max"
            WRITE(numout,*) "and will break the conservation of veget_max in land use change"
            STOP
         ENDIF
       ENDIF
       !! 4.2 Output the respiration terms and the net primary
       !!     production (NPP) that are calculated in STOMATE

       ! 4.2.1 Output the 3 respiration terms
       CALL xios_orchidee_send_field("GSL",GSL)
       CALL xios_orchidee_send_field("peatPET_lastyear",peatPET_lastyear)
       CALL xios_orchidee_send_field("precipitation_lastsummer",precipitation_lastsummer)
       CALL xios_orchidee_send_field('summerp_long',summerp_long)
       CALL xios_orchidee_send_field('summerpet_long',summerpet_long)

       CALL xios_orchidee_send_field("maint_resp",resp_maint/dt_sechiba)
       CALL xios_orchidee_send_field("hetero_resp",resp_hetero/dt_sechiba)
       CALL xios_orchidee_send_field("growth_resp",resp_growth/dt_sechiba)
       
       CALL histwrite_p(hist_id, 'maint_resp', kjit, resp_maint, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'hetero_resp', kjit, resp_hetero, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'growth_resp', kjit, resp_growth, kjpindex*nvm, indexveg)
       
       ! 4.2.2 Compute the net primary production as the diff from
       ! Gross primary productin and the growth and maintenance
       ! respirations
       npp(:,1)=zero
       DO j = 2,nvm
          npp(:,j) = gpp(:,j) - resp_growth(:,j) - resp_maint(:,j)
       ENDDO
       
       CALL xios_orchidee_send_field("npp",npp/dt_sechiba)
       
       CALL histwrite_p(hist_id, 'npp', kjit, npp, kjpindex*nvm, indexveg)
       
       IF ( hist2_id > 0 ) THEN
          CALL histwrite_p(hist2_id, 'maint_resp', kjit, resp_maint, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'hetero_resp', kjit, resp_hetero, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'growth_resp', kjit, resp_growth, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'npp', kjit, npp, kjpindex*nvm, indexveg)
       ENDIF
     
    ELSE
       !! ok_stomate is not activated
       !! Define the CO2 flux from the grid point to zero (no carbone cycle)
       co2_flux(:,:) = zero
    ENDIF

 
    !! 5. Do daily processes if necessary
    !!
    IF ( do_slow ) THEN

       !!  5.1 Calculate the LAI if STOMATE is not activated
       IF ( .NOT. ok_stomate ) THEN
          CALL slowproc_lai (kjpindex, lcanop,stempdiag, &
               lalo,resolution,lai,laimap)
          
          frac_age(:,:,1) = un
          frac_age(:,:,2) = zero
          frac_age(:,:,3) = zero
          frac_age(:,:,4) = zero
       ENDIF

!       !! 5.2.0 crop rotation, if rotation started
!       DO ji = 1,kjpindex
!          DO jv = 2,nvm
!              IF (ok_LAIdev(jv) .AND. f_rot_sechiba(ji,jv)) THEN !! whether to rotate
!               
!              ENDIF
!          ENDDO 
!       ENDDO
!        
!       !! end rotation, xuhui

       !! 5.2 Update veget
       CALL slowproc_veget (kjpindex, f_rot_sech, lai, frac_nobio, totfrac_nobio, veget_max, veget, soiltile)

       !! 5.3 updates qsintmax and other derived variables
       IF ( .NOT. ok_stomate ) THEN
          CALL slowproc_derivvar (kjpindex, veget, lai, &
               qsintmax, deadleaf_cover, assim_param, height, temp_growth)
       ELSE
          qsintmax(:,:) = qsintcst * veget(:,:) * lai(:,:)
          qsintmax(:,1) = zero
       ENDIF
    END IF

    !! 6. Calculate tot_bare_soil needed in hydrol, diffuco and condveg (fraction in the mesh)
    tot_bare_soil(:) = veget_max(:,1)
    DO jv = 2, nvm
       DO ji =1, kjpindex
          tot_bare_soil(ji) = tot_bare_soil(ji) + (veget_max(ji,jv) - veget(ji,jv))
       ENDDO
    END DO
    

    !! 7. Do some basic tests on the surface fractions updated above, only if
    !!    slowproc_veget has been done (do_slow). No change of the variables. 
    IF (do_slow) THEN
        CALL slowproc_checkveget(kjpindex, frac_nobio, veget_max, veget, tot_bare_soil, soiltile)
    END IF  

    !! 8. Write output fields
    CALL xios_orchidee_send_field("tot_bare_soil",tot_bare_soil)
    
    IF ( .NOT. almaoutput) THEN
       CALL histwrite_p(hist_id, 'tot_bare_soil', kjit, tot_bare_soil, kjpindex, IndexLand)
    END IF


    IF (printlev_loc>=3) WRITE (numout,*) ' slowproc_main done '

  END SUBROUTINE slowproc_main


!! ================================================================================================================================
!! SUBROUTINE 	: slowproc_finalize
!!
!>\BRIEF         Write to restart file variables for slowproc module and call finalization of stomate module
!!
!! DESCRIPTION : 
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_finalize (kjit,       kjpindex,  rest_id,  IndexLand,  &
                                njsc,       lai,       height,   veget,      &
                                frac_nobio, veget_max, reinf_slope,          & 
                                zz_deep, zz_coef_deep, thawed_humidity, depth_organic_soil, &
                                assim_param, frac_age, altmax, fpeat, wtpmax, wtp_year) !!!qcj++

!! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                           :: kjit           !! Time step number
    INTEGER(i_std),INTENT(in)                            :: kjpindex       !! Domain size - terrestrial pixels only
    INTEGER(i_std),INTENT (in)                           :: rest_id        !! Restart file identifier
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)     :: IndexLand      !! Indices of the points on the land map
    INTEGER(i_std), DIMENSION(kjpindex), INTENT(in)      :: njsc           !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: lai            !! Leaf area index (m^2 m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: height         !! height of vegetation (m)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: veget          !! Fraction of vegetation type including none biological fraction (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)    :: altmax
!!!qcj++ peatland
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)        :: fpeat
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)    :: wtpmax 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)    :: wtp_year

    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in) :: frac_nobio     !! Fraction of ice, lakes, cities etc. in the mesh
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: veget_max      !! Maximum fraction of vegetation type including none biological fraction (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)         :: reinf_slope    !! slope coef for reinfiltration
    REAL(r_std), DIMENSION(ndeep),   INTENT (in)            :: zz_deep        !! deep vertical profile
    REAL(r_std), DIMENSION(ndeep),   INTENT (in)            :: zz_coef_deep   !! deep vertical profile   
    REAL(r_std),DIMENSION (kjpindex,nvm,npco2),INTENT (in):: assim_param   !! min+max+opt temperatures & vmax for photosynthesis (K, \mumol m^{-2} s^{-1})
    REAL(r_std),DIMENSION (kjpindex,nvm,nleafages), INTENT(in):: frac_age  !! Age efficacity from STOMATE for isoprene
    REAL(r_std), DIMENSION(kjpindex), INTENT (inout)        :: thawed_humidity!! specified humidity of thawed soil
    REAL(r_std), DIMENSION(kjpindex), INTENT (inout)        :: depth_organic_soil !! how deep is the organic soil?

!! 0.4 Local variables
    REAL(r_std)                                          :: tmp_day(1)     !! temporary variable for I/O
    INTEGER                                              :: jf             !! Indice
    CHARACTER(LEN=4)                                     :: laistring      !! Temporary character string
    CHARACTER(LEN=80)                                    :: var_name       !! To store variables names for I/O
!_ ================================================================================================================================

    IF (printlev_loc>=3) WRITE (numout,*) 'Write restart file with SLOWPROC variables '

    ! 2.1 Write a series of variables controled by slowproc: day
    ! counter, vegetation fraction, max vegetation fraction, LAI
    ! variable from stomate, fraction of bare soil, soiltype
    ! fraction, clay fraction, height of vegetation, map of LAI
    
    CALL restput_p (rest_id, 'veget', nbp_glo, nvm, 1, kjit, veget, 'scatter',  nbp_glo, index_g)

    CALL restput_p (rest_id, 'veget_max', nbp_glo, nvm, 1, kjit, veget_max, 'scatter',  nbp_glo, index_g)

    CALL restput_p (rest_id, 'lai', nbp_glo, nvm, 1, kjit, lai, 'scatter',  nbp_glo, index_g)

    CALL restput_p (rest_id, 'frac_nobio', nbp_glo, nnobio, 1, kjit, frac_nobio, 'scatter',  nbp_glo, index_g)

    CALL restput_p (rest_id, 'frac_age', nbp_glo, nvm, nleafages, kjit, frac_age, 'scatter',  nbp_glo, index_g)

    ! Add the soil_classif as suffix for the variable name of njsc when it is stored in the restart file. 
    IF (soil_classif == 'zobler') THEN
       var_name= 'njsc_zobler'
    ELSE IF (soil_classif == 'usda') THEN
       var_name= 'njsc_usda'
    END IF
    CALL restput_p (rest_id, var_name, nbp_glo, 1, 1, kjit, REAL(njsc, r_std), 'scatter',  nbp_glo, index_g)
    
    IF ( hydrol_cwrr ) THEN
       CALL restput_p (rest_id, 'reinf_slope', nbp_glo, 1, 1, kjit, reinf_slope, 'scatter',  nbp_glo, index_g)
    END IF
       
    CALL restput_p (rest_id, 'clay_frac', nbp_glo, 1, 1, kjit, clayfraction, 'scatter',  nbp_glo, index_g)
    !
    ! The height of the vegetation could in principle be recalculated at the beginning of the run.
    ! However, this is very tedious, as many special cases have to be taken into account. This variable
    ! is therefore saved in the restart file.
    CALL restput_p (rest_id, 'height', nbp_glo, nvm, 1, kjit, height, 'scatter',  nbp_glo, index_g)
    !
    ! Specific case where the LAI is read and not calculated by STOMATE: need to be saved
    IF (read_lai) THEN     
       CALL restput_p (rest_id, 'laimap', nbp_glo, nvm, 12, kjit, laimap)
    ENDIF
    !
    ! If there is some land use change, write the year for the land use ??? 
    IF (map_pft_format) THEN
       CALL restput_p (rest_id, 'veget_year', kjit, veget_year)
    ENDIF
!!!qcj++ peatland
    CALL restput_p (rest_id, 'peatPET_last', nbp_glo, 1, 1, kjit, peatPET_lastyear, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'growth_day', nbp_glo, 1, 1, kjit,growth_day, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'GSL', nbp_glo, 1, 1, kjit,GSL, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'peatPET_this', nbp_glo, 1, 1, kjit, peatPET_thisyear, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'precipitation_last', nbp_glo, 1, 1, kjit,precipitation_lastsummer, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'precipitation_this', nbp_glo, 1, 1, kjit,precipitation_thissummer, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'summerpet_longterm', nbp_glo, 1, 1, kjit, summerpet_long, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'summerp_longterm', nbp_glo, 1, 1, kjit,summerp_long, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'peatC', nbp_glo, 1, 1, kjit,peatC,'scatter', nbp_glo, index_g)
    CALL restput_p (rest_id, 'peatC_ok', nbp_glo, 1, 1, kjit, peatC_ok, 'scatter', nbp_glo, index_g)
    CALL restput_p (rest_id, 'fpeat_map',nbp_glo, nstm, 1, kjit,fpeat_map,'scatter',nbp_glo, index_g)
!    CALL restput_p (rest_id, 't2m_last', nbp_glo, 1, 1, kjit,t2m_lastyear, 'scatter',  nbp_glo, index_g)
!    CALL restput_p (rest_id, 't2m_this', nbp_glo, 1, 1, kjit,t2m_thisyear, 'scatter',  nbp_glo, index_g)

    ! 2.2 Write restart variables managed by STOMATE
    IF ( ok_stomate ) THEN
       CALL stomate_finalize (kjit,    kjpindex,     indexLand,       clayfraction, & 
                              zz_deep, zz_coef_deep, thawed_humidity, depth_organic_soil, &
                              assim_param,altmax,fpeat,wtpmax,wtp_year)  !!!qcj++ peatland

    ENDIF
    
  END SUBROUTINE slowproc_finalize


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_init
!!
!>\BRIEF         Initialisation of all variables linked to SLOWPROC
!!
!! DESCRIPTION  : (definitions, functional, design, flags): The subroutine manages 
!! diverses tasks:
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::lcanop, ::veget_update, ::veget_year,
!! ::lai, ::veget, ::frac_nobio, ::totfrac_nobio, ::veget_max, ::height, ::soiltype
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_init (kjit, kjpindex, IndexLand, lalo, neighbours, resolution, contfrac, &
       rest_id, lai, frac_age, veget, frac_nobio, totfrac_nobio, soiltile, reinf_slope, &
       veget_max, tot_bare_soil, njsc, &
       height, lcanop, veget_update, veget_year, f_rot_sech)
    
    !! INTERFACE DESCRIPTION

    !! 0.1 Input variables
    INTEGER(i_std), INTENT (in)                           :: kjit           !! Time step number
    INTEGER(i_std), INTENT (in)                           :: kjpindex       !! Domain size - Terrestrial pixels only 
    INTEGER(i_std), INTENT (in)                           :: rest_id        !! Restart file identifier
    
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)      :: IndexLand      !! Indices of the land points on the map
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)       :: lalo           !! Geogr. coordinates (latitude,longitude) (degrees)
    INTEGER(i_std), DIMENSION (kjpindex,NbNeighb), INTENT(in):: neighbours  !! Vector of neighbours for each grid point
                                                                            !! (1=North and then clockwise)
    REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)       :: resolution     !! size in x and y of the grid (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)         :: contfrac       !! Fraction of continent in the grid (unitless)
    
    !! 0.2 Output variables
    INTEGER(i_std), INTENT(out)                           :: lcanop         !! Number of Canopy level used to compute LAI
    INTEGER(i_std), INTENT(out)                           :: veget_update   !! update frequency in timesteps (years) for landuse
    INTEGER(i_std), INTENT(out)                           :: veget_year     !! first year for landuse   (year or index ???)
    
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)    :: lai            !! Leaf Area index (m^2 / m^2)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)    :: veget          !! Fraction of vegetation type in the mesh (unitless)
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (out) :: frac_nobio     !! Fraction of ice,lakes,cities, ... in the mesh (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)        :: totfrac_nobio  !! Total fraction of ice+lakes+cities+... in the mesh (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)    :: veget_max      !! Max fraction of vegetation type in the mesh (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)        :: tot_bare_soil  !! Total evaporating bare soil fraction in the mesh
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)    :: height         !! Height of vegetation or surface in genral ??? (m)
    REAL(r_std),DIMENSION (kjpindex,nvm,nleafages), INTENT (out):: frac_age !! Age efficacity from STOMATE for isoprene
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT(out)   :: soiltile       !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)        :: reinf_slope    !! slope coef for reinfiltration
    INTEGER(i_std), DIMENSION(kjpindex), INTENT(out)      :: njsc           !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    LOGICAL, DIMENSION(kjpindex), INTENT(in)              :: f_rot_sech     !! whether a grid is under rotation
    
    !! 0.3 Local variables 
    REAL(r_std)                                           :: tmp_veget_year(1) !! temporary variable
    REAL(r_std)                                           :: zcanop            !! ???? soil depth taken for canopy
    INTEGER(i_std)                                        :: vtmp(1)           !! temporary variable
    REAL(r_std), DIMENSION(nslm)                          :: zsoil             !! soil depths at diagnostic levels
    CHARACTER(LEN=4)                                      :: laistring         !! Temporary character string
    INTEGER(i_std)                                        :: l, jf             !! Indices
    CHARACTER(LEN=80)                                     :: var_name          !! To store variables names for I/O
    INTEGER(i_std)                                        :: ji, jv, ier,jst   !! Indices 
    LOGICAL                                               :: get_slope
    REAL(r_std)                                           :: frac_nobio1       !! temporary variable for frac_nobio(see above)
    REAL(r_std), DIMENSION(kjpindex)                      :: tmp_real
    REAL(r_std), DIMENSION(kjpindex,nslm)                 :: stempdiag2_bid    !! matrix to store stempdiag_bid
    REAL(r_std), DIMENSION (kjpindex,nscm)                :: soilclass         !! Fractions of each soil textural class in the grid cell (0-1, unitless)
    CHARACTER(LEN=30), SAVE                               :: veget_str         !! update frequency for landuse
    !$OMP THREADPRIVATE(veget_str)
    REAL(r_std), DIMENSION(kjpindex)                      :: frac_crop_tot     !! Total fraction occupied by crops (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm)                  :: veget_nextyear    !! Temporary variable for new veget_max
    REAL(r_std),DIMENSION (kjpindex,nnobio)               :: frac_nobio_nextyear!! Temporary variable for new frac_nobio
    LOGICAL                                               :: found_restart     !! found_restart=true if all 3 variables veget_max, veget and 
                                                                               !! frac_nobio are read from restart file

  LOGICAL, SAVE                                         :: read_veg_map_fr_restfile = .FALSE. 
  !spitfire
  CHARACTER(LEN=80)                                     :: data_filename
  !endspit
    INTEGER(i_std)                                     :: ivma,ivm,jvm         !! Indices
  REAL(r_std),DIMENSION (nbp_glo,nvm)                   :: veget_max_g       !! Fraction of vegetation type at global scale
  REAL(r_std),DIMENSION (nbp_glo,nnobio)                :: frac_nobio_g      !! Fraction of ice, lakes, cities etc. in the mesh (global)
  REAL(r_std),DIMENSION (kjpindex,nvmap)                :: veget_max_map     !! Fraction of vegetation for MTCs, used to hold the vegetation
                                                                             !! fractions for the first year of spin-up read from a MTC land cover map.
  REAL(r_std),DIMENSION (kjpindex,nvmap)                :: veget_ny_map      !! Fraction of vegetation for MTCs, it's a dummy variable
                                                                             !! used to initialization of veget_max_map for the first year
                                                                             !! run from scratch.
!gmjc soil field capacity
    INTEGER(i_std)                                        :: temp_njsc
!end gmjc
!_ ================================================================================================================================

    ! Initialize local printlev
    printlev_loc=get_printlev('slowproc')
    IF (printlev_loc>=3) WRITE (numout,*) "In slowproc_init"
    
    
    !! 1. Allocation 

    ALLOCATE (clayfraction(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable clayfraction','','')
    clayfraction(:)=undef_sechiba

    ! Initialisation of the fraction of the different vegetation: Start with 100% of bare soil
    ALLOCATE (soilclass_default(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable soilclass_default','','')
    soilclass_default(:)=undef_sechiba

    ! Allocation of the fraction of non biospheric areas 
    ALLOCATE(frac_nobio_lastyear(kjpindex, nnobio), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable frac_nobio_lastyear','','')
    frac_nobio_lastyear(:,:) = zero

    ALLOCATE(veget_max_new(kjpindex, nvm), STAT=ier)                            
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable veget_max_new','','')

    ! Allocation of the fraction of non biospheric areas 
    ALLOCATE(frac_nobio_new(kjpindex, nnobio), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable frac_nobio_new','','')
                                                                               
    ! Allocate laimap
    IF (read_lai)THEN
       ALLOCATE (laimap(kjpindex,nvm,12),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable laimap','','')
    ELSE
       ALLOCATE (laimap(1,1,1), stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable laimap(1,1,1)','','')
    ENDIF 
!gmjc allocate fc_grazing
    ALLOCATE(fc_grazing(kjpindex), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable fc_grazing','','')
!end gmjc
    ! Allocation of next year vegetation fraction in case of land use change
    ier=-1
    ALLOCATE(vegetnew_firstday(kjpindex, nvm), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable vegetnew_firstday','','')
    vegetnew_firstday(:,1) = un
    vegetnew_firstday(:,2:nvm) = zero

!!!qcj++ peatland
    ALLOCATE(veget_cov_max_adjusted(kjpindex, nvm), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable veget_cov_max_adjusted','','')

    ALLOCATE(peatPET_lastyear(kjpindex), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable peatPET_lastyear','','')

    ALLOCATE(peatPET_thisyear(kjpindex), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable peatPET_thisyear','','')

    ALLOCATE(precipitation_lastsummer(kjpindex), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable precipitation_lastsummer','','')

    ALLOCATE(precipitation_thissummer(kjpindex), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable precipitation_thissummer','','')

    ALLOCATE(summerp_long(kjpindex), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable summerp_long','','')

    ALLOCATE(summerpet_long(kjpindex), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable summerpet_long','','')

    ALLOCATE(peatC(kjpindex), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable peatC','','')

    ALLOCATE(peatC_ok(kjpindex), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable peatC_ok','','')

    ALLOCATE(growth_day(kjpindex), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable growth_day','','')

    ALLOCATE(GSL(kjpindex), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable GSL','','')

    ALLOCATE(fpeat_map(kjpindex,nstm), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable fpeat_map','','')

 !   ALLOCATE(t2m_lastyear(kjpindex), STAT=ier)
 !   IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable t2m_lastyear','','')

  !  ALLOCATE(t2m_thisyear(kjpindex), STAT=ier)
  !  IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable t2m_thisyear','','')

    !spitfire 
    ier=-1
    ALLOCATE (m_lightn(kjpindex,12),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable m_lightn','','')
    m_lightn(:,:) = zero
    !
    ier=-1
    ALLOCATE (glccSecondShift(kjpindex,12),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable glccSecondShift','','')
    glccSecondShift(:,:) = zero
    !
    ier=-1
    ALLOCATE (glccPrimaryShift(kjpindex,12),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable glccPrimaryShift','','')
    glccPrimaryShift(:,:) = zero
    !
    ier=-1
    ALLOCATE (glccNetLCC(kjpindex,12),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable glccNetLCC','','')
    glccNetLCC(:,:) = zero
    !
    ier=-1
    ALLOCATE (harvest_matrix(kjpindex,12),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable harvest_matrix','','')
    harvest_matrix(:,:) = zero
    !
    ier=-1
    ALLOCATE (harvest_biomass(kjpindex,12),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable harvest_biomass','','')
    harvest_biomass(:,:) = zero
    !
    ier=-1
    ALLOCATE (bound_spa(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable bound_spa','','')
    bound_spa(:,:) = zero
    !
    ier=-1
    ALLOCATE (newvegfrac(kjpindex,nvmap),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable newvegfrac','','')
    newvegfrac(:,:) = zero
    !
    ier=-1
    ALLOCATE (proxy_anidens(kjpindex,12),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable proxy_anidens','','')
    proxy_anidens(:,:) = zero
    !
    ier=-1
    ALLOCATE (m_observed_ba(kjpindex,12),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable m_observed_ba','','')
    m_observed_ba(:,:) = zero
    !
    ier=-1
    ALLOCATE (m_cf_coarse(kjpindex,12),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable m_cf_coarse','','')
    m_cf_coarse(:,:) = zero
    !
    ier=-1
    ALLOCATE (m_cf_fine(kjpindex,12),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable m_cf_fine','','')
    m_cf_fine(:,:)= zero
    !
    ier=-1
    ALLOCATE (m_ratio(kjpindex,12),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable m_ratio','','')
    m_ratio(:,:) =zero
    !
    ier=-1
    ALLOCATE (m_ratio_flag(kjpindex,12),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable m_ratio_flag','','')
    m_ratio_flag(:,:) = zero
    !
    ier=-1
    ALLOCATE (popd(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable popd','','')
    popd(:) = zero
    !
    ier=-1
    ALLOCATE (humign(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable humign','','')
    humign(:) = zero
    !endspit
    
    !! 2. Read variables from restart file

    found_restart=.TRUE.
    var_name= 'veget'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Vegetation fraction')
    CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., veget, "gather", nbp_glo, index_g)
    IF ( ALL( veget(:,:) .EQ. val_exp ) ) found_restart=.FALSE.

    var_name= 'veget_max'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Maximum vegetation fraction')
    CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., veget_max, "gather", nbp_glo, index_g)
    IF ( ALL( veget_max(:,:) .EQ. val_exp ) ) found_restart=.FALSE.
    !chaoprint
    !print *,'veget_max after reading restart',veget_max

    ! Get frac_nobio from the restart file 
    var_name= 'frac_nobio'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Special soil type fraction')
    CALL restget_p (rest_id, var_name, nbp_glo, nnobio, 1, kjit, .TRUE., frac_nobio, "gather", nbp_glo, index_g)
    IF ( ALL( frac_nobio(:,:) .EQ. val_exp ) ) found_restart=.FALSE.
!
    IF (map_pft_format .AND. .NOT. impveg) THEN
          IF (veget_reinit) THEN
             ! Do not take the value read from restart file
             veget_year=veget_year_orig
          ELSE 
            var_name= 'veget_year'
            CALL ioconf_setatt_p('UNITS', '-')
            CALL ioconf_setatt_p('LONG_NAME','Last year get in Land Use file.')
            CALL restget_p (rest_id, 'veget_year', kjit, .TRUE., val_exp, veget_year)
            IF (veget_year == val_exp) THEN
               ! veget_year was not found in restart file
               veget_year=veget_year_orig
            ENDIF
          ENDIF
!       ENDIF
!       CALL bcast(veget_year)

       !
       !Config Key   = VEGET_UPDATE
       !Config Desc  = Update vegetation frequency
       !Config If    = MAP_PFT_FORMAT
       !Config Def   = 0Y
       !Config Help  = The veget datas will be update each this time step.
       !Config Units = [years]
       !
       veget_update=0
       WRITE(veget_str,'(a)') '0Y'
       CALL getin_p('VEGET_UPDATE', veget_str)
       l=INDEX(TRIM(veget_str),'Y')
       READ(veget_str(1:(l-1)),"(I2.2)") veget_update
       IF (printlev_loc >= 2) WRITE(numout,*) "Update frequency for land use in years :",veget_update

       ! Coherence test
       IF (veget_update > 0 .AND. ok_dgvm .AND. .NOT. agriculture) THEN
          CALL ipslerr_p(3,'slowproc_init',&
               'The combination DGVM=TRUE, AGRICULTURE=FALSE and VEGET_UPDATE>0 is not possible', &
               'Set VEGET_UPDATE=0Y in run.def','')
       END IF
    ELSE
       ! map_pft_format=FALSE or impveg=TRUE: there can not be any land use change, veget_update must be =0
       ! Read VEGET_UPDATE from run.def and exit if it is different from 0Y
       veget_update=0
       WRITE(veget_str,'(a)') '0Y'
       CALL getin_p('VEGET_UPDATE', veget_str)
       l=INDEX(TRIM(veget_str),'Y')
       READ(veget_str(1:(l-1)),"(I2.2)") veget_update
       IF (veget_update /= 0) THEN
          WRITE(numout,*) 'veget_update=',veget_update,' is not coeherent with map_pft_format=',map_pft_format,' or impveg=',impveg
          CALL ipslerr_p(3,'slowproc_init','Incoherent values between impveg, map_pft_format and veget_update', &
               'veget_update must be equal to 0 if map_pft_format=false or if impveg=true','')
       END IF

    ENDIF

    IF (printlev_loc>=3) WRITE (numout,*) 'slowproc_init : End of Land Use configuration'

    IF ( hydrol_cwrr ) THEN
       var_name= 'reinf_slope'
       CALL ioconf_setatt_p('UNITS', '-')
       CALL ioconf_setatt_p('LONG_NAME','Slope coef for reinfiltration')
       CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., reinf_slope, "gather", nbp_glo, index_g)
    END IF
    
    ! Below we define the soil texture of the grid-cells
    ! Add the soil_classif as suffix for the variable name of njsc when it is stored in the restart file. 
    IF (soil_classif == 'zobler') THEN
       var_name= 'njsc_zobler'
    ELSE IF (soil_classif == 'usda') THEN
       var_name= 'njsc_usda'
    ELSE
       CALL ipslerr_p(3,'slowproc_init','Non supported soil type classification','','')
    END IF

    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Index of soil type')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., tmp_real, "gather", nbp_glo, index_g)
    IF ( ALL( tmp_real(:) .EQ. val_exp) ) THEN
       njsc (:) = undef_int
    ELSE
       njsc = NINT(tmp_real)
    ENDIF
    
    var_name= 'clay_frac'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Fraction of clay in each mesh')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., clayfraction, "gather", nbp_glo, index_g)

    IF (printlev_loc>=3) WRITE (numout,*) 'slowproc_init : End CWRR configuration'
    !
    var_name= 'lai'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Leaf area index')
    CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., lai, "gather", nbp_glo, index_g)

    ! The height of the vegetation could in principle be recalculated at the beginning of the run.
    ! However, this is very tedious, as many special cases have to be taken into account. This variable
    ! is therefore saved in the restart file.
    var_name= 'height'
    CALL ioconf_setatt_p('UNITS', 'm')
    CALL ioconf_setatt_p('LONG_NAME','Height of vegetation')
    CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., height, "gather", nbp_glo, index_g)
 
    IF (read_lai)THEN
       var_name= 'laimap'
       CALL ioconf_setatt_p('UNITS', '-')
       CALL ioconf_setatt_p('LONG_NAME','Leaf area index read')
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, 12, kjit, .TRUE., laimap)
    ENDIF

    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Fraction of leaves in leaf age class ')
    CALL restget_p (rest_id, 'frac_age', nbp_glo, nvm, nleafages, kjit, .TRUE.,frac_age, "gather", nbp_glo, index_g)

!!!qcj++ peatland
    peatPET_lastyear(:) = val_exp
    CALL ioconf_setatt_p('UNITS', 'mm/year')
    CALL ioconf_setatt_p('LONG_NAME','peatPET_lastyear')
    CALL restget_p (rest_id, 'peatPET_last', nbp_glo, 1, 1, kjit, .TRUE.,peatPET_lastyear, "gather", nbp_glo, index_g)
    IF ( ALL( peatPET_lastyear(:) .EQ. val_exp ) ) THEN
       peatPET_lastyear(:) = undef_sechiba   !!! to make peatPET_lastyear > precip_lastyear (precip_crit) in the first year, so that peatland will not expand
    ENDIF
    
    precipitation_lastsummer(:) = val_exp
    CALL ioconf_setatt_p('UNITS', 'mm')
    CALL ioconf_setatt_p('LONG_NAME','precipitation_lastsummer')
    CALL restget_p (rest_id, 'precipitation_last', nbp_glo, 1, 1, kjit, .TRUE.,precipitation_lastsummer, "gather", nbp_glo, index_g)
    IF ( ALL( precipitation_lastsummer(:) .EQ. val_exp ) ) THEN
       precipitation_lastsummer(:) = precip_crit  
    ENDIF

    precipitation_thissummer(:) =val_exp
    CALL ioconf_setatt_p('UNITS', 'mm')
    CALL ioconf_setatt_p('LONG_NAME','precipitation_thissummer')
    CALL restget_p (rest_id, 'precipitation_this', nbp_glo, 1, 1, kjit, .TRUE.,precipitation_thissummer, "gather", nbp_glo, index_g)
    IF ( ALL( precipitation_thissummer(:) .EQ. val_exp ) ) THEN
       precipitation_thissummer(:) = zero
    ENDIF
    
    peatPET_thisyear(:) = val_exp
    CALL ioconf_setatt_p('UNITS', 'mm/year')
    CALL ioconf_setatt_p('LONG_NAME','peatPET_thisyear')
    CALL restget_p (rest_id, 'peatPET_this', nbp_glo, 1, 1, kjit, .TRUE.,peatPET_thisyear, "gather", nbp_glo, index_g)
    IF ( ALL( peatPET_thisyear(:) .EQ. val_exp ) ) THEN
       peatPET_thisyear(:) = zero  
    ENDIF

    growth_day(:) = val_exp
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Days with monthly mean air T greater than 0')   
    CALL restget_p (rest_id, 'growth_day', nbp_glo, 1, 1, kjit, .TRUE.,growth_day, "gather", nbp_glo, index_g)
    IF ( ALL( growth_day(:) .EQ. val_exp ) ) THEN
       growth_day(:) = 0
    ENDIF

    GSL(:) = val_exp
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','growing season length')
    CALL restget_p (rest_id, 'GSL', nbp_glo, 1, 1, kjit, .TRUE., GSL, "gather", nbp_glo, index_g)
    IF ( ALL( GSL(:) .EQ. val_exp ) ) THEN
       GSL(:) = 1.
    ENDIF

    summerp_long(:) = val_exp
    CALL ioconf_setatt_p('UNITS', 'mm/year')
    CALL ioconf_setatt_p('LONG_NAME','summerp_long')
    CALL restget_p (rest_id, 'summerp_longterm', nbp_glo, 1, 1, kjit,.TRUE.,summerp_long, "gather", nbp_glo, index_g)
    IF ( ALL( summerp_long(:) .EQ. val_exp ) ) THEN 
       summerp_long(:) = zero
    ENDIF

    summerpet_long(:) = val_exp
    CALL ioconf_setatt_p('UNITS', 'mm/year')
    CALL ioconf_setatt_p('LONG_NAME','summerpet_long')
    CALL restget_p (rest_id, 'summerpet_longterm', nbp_glo, 1, 1, kjit, .TRUE.,summerpet_long, "gather", nbp_glo, index_g)
    IF ( ALL( summerpet_long(:) .EQ. val_exp ) ) THEN
       summerpet_long(:) = undef_sechiba !!! to make sure peatland will not expand in the first year
    ENDIF

    peatC(:) = val_exp
    CALL ioconf_setatt_p('UNITS', 'kg/m2')
    CALL ioconf_setatt_p('LONG_NAME','peatC')
    CALL restget_p (rest_id, 'peatC', nbp_glo, 1, 1, kjit, .TRUE.,peatC,"gather", nbp_glo, index_g)
    IF ( ALL( peatC(:) .EQ. val_exp ) ) THEN
      peatC(:) = 0  !!!Carbon is not enough for peat to form
    ENDIF

    peatC_ok(:) = val_exp
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','peatC_ok')
    CALL restget_p (rest_id, 'peatC_ok', nbp_glo, 1, 1, kjit, .TRUE.,peatC_ok, "gather", nbp_glo, index_g)
    IF ( ALL( peatC_ok(:) .EQ. val_exp ) ) THEN
       peatC_ok(:) = 0.  !!!Carbon is not enough for peat to form
    ENDIF

    fpeat_map(:,:) = val_exp
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','peatland fraction from vegetation map')
    CALL restget_p (rest_id, 'fpeat_map', nbp_glo, nstm, 1, kjit, .TRUE.,fpeat_map,"gather", nbp_glo, index_g)
    IF ( ALL( fpeat_map(:,:) .EQ. val_exp ) ) THEN
      fpeat_map(:,:) = 0  
    ENDIF

!    t2m_lastyear(:) = val_exp
!    CALL ioconf_setatt_p('UNITS', 'K')
!    CALL ioconf_setatt_p('LONG_NAME','t2m_lastyear')
!    CALL restget_p (rest_id, 't2m_last', nbp_glo, 1, 1, kjit, .TRUE.,t2m_lastyear, "gather", nbp_glo, index_g)
!    IF ( ALL( t2m_lastyear(:) .EQ. val_exp ) ) THEN
!       t2m_lastyear(:) = undef_sechiba   !!! to make meant2m > 15 C, so that peatland will not expand
!    ENDIF

!    t2m_thisyear(:) = val_exp
!    CALL ioconf_setatt_p('UNITS', 'K')
!    CALL ioconf_setatt_p('LONG_NAME','t2m_thisyear')
!    CALL restget_p (rest_id, 't2m_this', nbp_glo, 1, 1, kjit, .TRUE.,t2m_thisyear, "gather", nbp_glo, index_g)
!    IF ( ALL( t2m_thisyear(:) .EQ. val_exp ) ) THEN
!       t2m_thisyear(:) = zero
!    ENDIF

    !! 3. Some other initializations

    !Config Key   = SECHIBA_ZCANOP
    !Config Desc  = Soil level used for canopy development (if STOMATE disactivated)
    !Config If    = OK_SECHIBA and .NOT. OK_STOMATE  
    !Config Def   = 0.5
    !Config Help  = The temperature at this soil depth is used to determine the LAI when
    !Config         STOMATE is not activated.
    !Config Units = [m]
    zcanop = 0.5_r_std
    CALL setvar_p (zcanop, val_exp, 'SECHIBA_ZCANOP', 0.5_r_std)

    ! depth at center of the levels
    zsoil(1) = diaglev(1) / 2.
    DO l = 2, nslm
       zsoil(l) = ( diaglev(l) + diaglev(l-1) ) / 2.
    ENDDO

    ! index of this level
    vtmp = MINLOC ( ABS ( zcanop - zsoil(:) ) )
    lcanop = vtmp(1)

    !
    !  Interception reservoir coefficient
    !
    !Config Key   = SECHIBA_QSINT 
    !Config Desc  = Interception reservoir coefficient
    !Config If    = OK_SECHIBA 
    !Config Def   = 0.1
    !Config Help  = Transforms leaf area index into size of interception reservoir
    !Config         for slowproc_derivvar or stomate
    !Config Units = [m]
    CALL getin_p('SECHIBA_QSINT', qsintcst)
    IF (printlev >= 2) WRITE(numout, *)' SECHIBA_QSINT, qsintcst = ', qsintcst


!qcj++ peatland
    PWT_lim = 60. !!mm/year
    CALL getin_p("PWT_LIM", PWT_lim)
    PC_lim = 50. !!kg/m2
    CALL getin_p("PC_LIM", PC_lim)
    sat_gsl = 1.
    CALL getin_p("SAT_GSL", sat_gsl)
    !! 4. Initialization of variables not found in restart file

    IF ( impveg ) THEN

       !! 4.1.a Case impveg=true: Initialization of variables by reading run.def
       !!       The routine setvar_p will only initialize the variable if it was not found in restart file. 
       !!       We are on a point and thus we can read the information from the run.def
       
       !Config Key   = SECHIBA_VEGMAX
       !Config Desc  = Maximum vegetation distribution within the mesh (0-dim mode)
       !Config If    = IMPOSE_VEG
       !Config Def   = 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0
       !Config Help  = The fraction of vegetation is read from the restart file. If
       !Config         it is not found there we will use the values provided here.
       !Config Units = [-]
       CALL setvar_p (veget_max, val_exp, 'SECHIBA_VEGMAX', veget_ori_fixed_test_1)

       !Config Key   = SECHIBA_FRAC_NOBIO
       !Config Desc  = Fraction of other surface types within the mesh (0-dim mode)
       !Config If    = IMPOSE_VEG
       !Config Def   = 0.0
       !Config Help  = The fraction of ice, lakes, etc. is read from the restart file. If
       !Config         it is not found there we will use the values provided here.
       !Config         For the moment, there is only ice.
       !Config Units = [-]
       frac_nobio1 = frac_nobio(1,1)
       CALL setvar_p (frac_nobio1, val_exp, 'SECHIBA_FRAC_NOBIO', frac_nobio_fixed_test_1)
       frac_nobio(:,:) = frac_nobio1
       
       !Config Key   = SECHIBA_LAI
       !Config Desc  = LAI for all vegetation types (0-dim mode)
       !Config Def   = 0., 8., 8., 4., 4.5, 4.5, 4., 4.5, 4., 2., 2., 2., 2.
       !Config If    = IMPOSE_VEG and .NOT. STOMATE
       !Config Help  = The maximum LAI used in the 0dim mode. The values should be found
       !Config         in the restart file. The new values of LAI will be computed anyway
       !Config         at the end of the current day. The need for this variable is caused
       !Config         by the fact that the model may stop during a day and thus we have not
       !Config         yet been through the routines which compute the new surface conditions.
       !Config Units = [-]
       ! MICT: lai cannot be modified when ok_stomate is enabled. This breaks consistency between
       ! other variables through biomass
       IF (.NOT. ok_stomate) CALL setvar_p (lai, val_exp, 'SECHIBA_LAI', llaimax)

       IF (.NOT. found_restart) THEN
          ! Call slowproc_veget to correct veget_max and to calculate veget and soiltiles
          CALL slowproc_veget (kjpindex, f_rot_sech, lai, frac_nobio, totfrac_nobio, veget_max, veget, soiltile)
       END IF

       IF (impsoilt) THEN

          ! If njsc is not in restart file, then initialize soilclass from values 
          ! from run.def file and recalculate njsc
          IF ( ALL(njsc(:) .EQ. undef_int )) THEN
             !Config Key   = SOIL_FRACTIONS
             !Config Desc  = Fraction of the 3 soil types (0-dim mode)
             !Config Def   = undef_sechiba
             !Config If    = IMPOSE_VEG and IMPOSE_SOILT
             !Config Help  = Determines the fraction for the 3 soil types
             !Config         in the mesh in the following order : sand loam and clay.
             !Config Units = [-]
          
             soilclass(1,:) = soilclass_default(:)
             CALL getin_p('SOIL_FRACTIONS',soilclass(1,:))
             ! Assign for each grid-cell the % of the different textural classes (up to 12 if 'usda')
             DO ji=2,kjpindex
                ! here we read, for the prescribed grid-cell, the % occupied by each of the soil texture classes 
                soilclass(ji,:) = soilclass(1,:)
             ENDDO

             ! Simplify an heterogeneous grid-cell into an homogeneous one with the dominant texture
             njsc(:) = 0
             DO ji = 1, kjpindex
                ! here we reduce to the dominant texture class
                njsc(ji) = MAXLOC(soilclass(ji,:),1)
             ENDDO
          END IF

          !Config Key   = CLAY_FRACTION
          !Config Desc  = Fraction of the clay fraction (0-dim mode)
          !Config Def   = 0.2
          !Config If    = IMPOSE_VEG and IMPOSE_SOIL
          !Config Help  = Determines the fraction of clay in the grid box.
          !Config Units = [-] 
          
          ! If clayfraction was not in restart file it will be read fro run.def file instead of deduced 
          ! based on fractions of each textural class
          CALL setvar_p (clayfraction, val_exp, 'CLAY_FRACTION', clayfraction_default)
       ELSE
          ! Case impsoilt=false and impveg=true
          IF ( MINVAL(clayfraction) .EQ. MAXVAL(clayfraction) .AND. MAXVAL(clayfraction) .EQ. val_exp .OR. &
               MINVAL(njsc) .EQ. MAXVAL(njsc) .AND. MAXVAL(njsc) .EQ. undef_int ) THEN
             
             CALL slowproc_soilt(kjpindex, lalo, neighbours, resolution, contfrac, soilclass, clayfraction)
             njsc(:) = 0
             DO ji = 1, kjpindex
                njsc(ji) = MAXLOC(soilclass(ji,:),1)
             ENDDO
          ENDIF
       ENDIF

       !Config Key   = REINF_SLOPE
       !Config Desc  = Slope coef for reinfiltration 
       !Config Def   = 0.1
       !Config If    = IMPOSE_VEG
       !Config Help  = Determines the reinfiltration ratio in the grid box due to flat areas
       !Config Units = [-]
       !
       slope_default=0.0
       CALL setvar_p (reinf_slope, val_exp, 'SLOPE', slope_default)

       !Config Key   = SLOWPROC_HEIGHT
       !Config Desc  = Height for all vegetation types 
       !Config Def   = 0., 30., 30., 20., 20., 20., 15., 15., 15., .5, .6, 1.0, 1.0
       !Config If    = OK_SECHIBA
       !Config Help  = The height used in the 0dim mode. The values should be found
       !Config         in the restart file. The new values of height will be computed anyway
       !Config         at the end of the current day. The need for this variable is caused
       !Config         by the fact that the model may stop during a day and thus we have not
       !Config         yet been through the routines which compute the new surface conditions.
       !Config Units = [m]
       CALL setvar_p (height, val_exp, 'SLOWPROC_HEIGHT', height_presc)


    ELSE IF ( .NOT. found_restart ) THEN
        !!! rotation is not concerned if no restarting of veget_max

       !! 4.1.b Case impveg=false and no restart files: Initialization by reading vegetation map
       
       ! Initialize veget_max and frac_nobio
       IF ( map_pft_format ) THEN
         IF(use_age_class) THEN
           veget_ny_map(:,:)=zero
           ! Note that veget_ny_map is only a dummy variable and not used at all,
           ! because veget_ny_map takes the position of veget_last, which is used
           ! only to retain the agriculture fraction when land cover is partially
           ! updated.
           ! `veget_max_map` is the variable that holds for us the MTC fraction 
           ! information from the land cover map.
           CALL slowproc_readvegetmax(kjpindex, lalo, neighbours, resolution, contfrac, &
                &               veget_ny_map, veget_max_map, frac_nobio_nextyear, veget_year, .TRUE.)
           ! Now we need to change this into the map with age classes. 
           ! Since this is an initialization, we will just assign the whole PFT to the youngest age class.
           veget_max(:,:)=zero
           DO ivma=1,nvmap
              veget_max(:,start_index(ivma))=veget_max_map(:,ivma)
           ENDDO
           frac_nobio          = frac_nobio_nextyear          
   
         !we're not using age groups so it's just "normal PFT map initialization"
         ELSE
           ! Case without restart file and map_pft_format=true
           IF (printlev_loc>=3) WRITE(numout,*) 'Before call slowproc_readvegetmax in initialization phase without restart files'
           IF (printlev_loc>=3) WRITE(numout,*) 'veget_year=', veget_year
           
           ! Call the routine to update the vegetation (output is veget_nextyear)
           CALL slowproc_readvegetmax(kjpindex, lalo, neighbours, resolution, contfrac, &
                veget_max, veget_nextyear, frac_nobio_nextyear, veget_year, .TRUE.)
           IF (printlev_loc>=4) WRITE (numout,*) 'After slowproc_readvegetmax in initialization phase'
          
           ! Update vegetation with values read from the file
           veget_max           = veget_nextyear
           frac_nobio          = frac_nobio_nextyear          
         END IF

       ELSE
          ! The interpolation of vegetation has changed.
          CALL getin_p('read_veg_map_fr_restfile', read_veg_map_fr_restfile)
          IF ( .NOT. read_veg_map_fr_restfile ) THEN
              ! slowproc interpol :
              CALL slowproc_interpol_g(nbp_glo, lalo_g, neighbours_g, resolution_g, contfrac_g, veget_max_g, frac_nobio_g)
          ELSE
              CALL slowproc_read_veg_restfile(nbp_glo, lalo_g, neighbours_g, resolution_g, contfrac_g, &
                   veget_max_g, frac_nobio_g)
          ENDIF

           CALL scatter(veget_max_g, veget_max)
           CALL scatter(frac_nobio_g, frac_nobio)

           ! map_pft_format=FALSE: Read and interpolate Olson type map
           CALL slowproc_interpol(kjpindex, lalo, neighbours, resolution, contfrac, veget_max, frac_nobio)
       END IF
       
       !! Reset totaly or partialy veget_max if using DGVM
       IF ( ok_dgvm  ) THEN
          ! If we are dealing with dynamic vegetation then all natural PFTs should be set to veget_max = 0
          ! In case no agriculture is desired, agriculture PFTS should be set to 0 as well
          IF (agriculture) THEN
             DO jv = 2, nvm
! dgvmjc consider pasture as not natural
                IF (natural(jv) .AND. .NOT. pasture(jv)) THEN
                   veget_max(:,jv)=zero
                ENDIF
             ENDDO
             
             ! Calculate the fraction of crop for each point.
             ! Sum only on the indexes corresponding to the non_natural pfts
             frac_crop_tot(:) = zero
             DO jv = 2, nvm
! dgvmjc consider pasture as not natural
                IF(.NOT. natural(jv) .OR. pasture(jv)) THEN
                   DO ji = 1, kjpindex
                      frac_crop_tot(ji) = frac_crop_tot(ji) + veget_max(ji,jv)
                   ENDDO
                ENDIF
             ENDDO
            
             ! Calculate the fraction of bare soil
             DO ji = 1, kjpindex
                veget_max(ji,1) = un - frac_crop_tot(ji) - SUM(frac_nobio(ji,:))
             ENDDO
          ELSE
             veget_max(:,:) = zero
             DO ji = 1, kjpindex
                veget_max(ji,1) = un  - SUM(frac_nobio(ji,:))
             ENDDO
          ENDIF ! end agriculture
       ENDIF ! end ok_dgvm
       

       ! correct veget_max and to calculate veget and soiltiles

       CALL slowproc_veget (kjpindex, f_rot_sech, lai, frac_nobio, totfrac_nobio, veget_max, veget, soiltile)


    END IF ! end impveg

    !! 4.2 Continue initializing variables not found in restart file. Case for both impveg=true and false.

    ! Initialize laimap for the case read_lai if not found in restart file
    IF (read_lai) THEN
       IF ( ALL( laimap(:,:,:) .EQ. val_exp) ) THEN
          ! Interpolation of LAI
          CALL slowproc_interlai (kjpindex, lalo, resolution,  neighbours, contfrac, laimap)
       ENDIF
    ENDIF
    
    ! Initialize lai if not found in restart file and not already initialized using impveg
    IF ( MINVAL(lai) .EQ. MAXVAL(lai) .AND. MAXVAL(lai) .EQ. val_exp) THEN
       IF (read_lai) THEN
          stempdiag2_bid(1:kjpindex,1:nslm) = stempdiag_bid
          CALL slowproc_lai (kjpindex, lcanop, stempdiag2_bid, &
               lalo,resolution,lai,laimap)
       ELSE
          ! If we start from scratch, we set lai to zero for consistency with stomate
          lai(:,:) = zero
       ENDIF
       
       frac_age(:,:,1) = un
       frac_age(:,:,2) = zero
       frac_age(:,:,3) = zero
       frac_age(:,:,4) = zero
    ENDIF
    
    ! Initialize heigth if not found in restart file and not already initialized using impveg
    IF ( MINVAL(height) .EQ. MAXVAL(height) .AND. MAXVAL(height) .EQ. val_exp) THEN
       ! Impose height
       DO jv = 1, nvm
          height(:,jv) = height_presc(jv)
       ENDDO
    ENDIF
    
    ! Initialize clayfraction and njsc if not found in restart file and not already initialized using impveg
    IF ( MINVAL(clayfraction) .EQ. MAXVAL(clayfraction) .AND. MAXVAL(clayfraction) .EQ. val_exp .OR. &
         MINVAL(njsc) .EQ. MAXVAL(njsc) .AND. MAXVAL(njsc) .EQ. undef_int ) THEN
       
       IF (printlev_loc>=4) WRITE (numout,*) 'clayfraction or njcs were not in restart file, call slowproc_soilt'
       CALL slowproc_soilt(kjpindex, lalo, neighbours, resolution, contfrac, soilclass, clayfraction)
       IF (printlev_loc>=4) WRITE (numout,*) 'After slowproc_soilt'
       
       njsc(:) = 0
       DO ji = 1, kjpindex
          njsc(ji) = MAXLOC(soilclass(ji,:),1)
       ENDDO
    ENDIF

!gmjc 15Feb2016 avoid grazing wet
    ! NOTE only for usda 12 classes soil map
    WRITE (numout,*) 'get fc_grazing gmjc'
    SELECTCASE(soil_classif)
    CASE('none')
       DO ji = 1, kjpindex
         temp_njsc = njsc(ji)
         fc_grazing(ji) = mcs_fao(temp_njsc)
       ENDDO    
    CASE('zobler')
       DO ji = 1, kjpindex
         temp_njsc = njsc(ji)
         fc_grazing(ji) = mcs_fao(temp_njsc)
       ENDDO 
    CASE("usda")
       DO ji = 1, kjpindex
         temp_njsc = njsc(ji)
         fc_grazing(ji) = mcs_usda(temp_njsc) 
!         fc_grazing(ji) = 1.0
       ENDDO
    CASE DEFAULT
       WRITE(*,*) 'A non supported soil type classification has been chosen'
       CALL ipslerr_p(3,'slowproc_soilt','non supported soil type classification','','')
    ENDSELECT
!end gmjc
    !Config Key   = GET_SLOPE
    !Config Desc  = Read slopes from file and do the interpolation
    !Config Def   = n
    !Config If    =
    !Config Help  = Needed for reading the slopesfile and doing the interpolation. This will be
    !               used by the re-infiltration parametrization
    !Config Units = [FLAG]
    get_slope = .FALSE.
    CALL getin_p('GET_SLOPE',get_slope)
    
    IF ( hydrol_cwrr ) THEN
       IF ( MINVAL(reinf_slope) .EQ. MAXVAL(reinf_slope) .AND. MAXVAL(reinf_slope) .EQ. val_exp .OR. get_slope) THEN
          IF (printlev_loc>=4) WRITE (numout,*) 'reinf_slope was not in restart file. Now call slowproc_slope'
          
          CALL slowproc_slope(kjpindex, lalo, neighbours, resolution, contfrac, reinf_slope)
          IF (printlev_loc>=4) WRITE (numout,*) 'After slowproc_slope'
          
       ENDIF
    END IF

    
    !! 5. Some calculations always done, with and without restart files
       
    ! The variables veget, veget_max and frac_nobio were all read from restart file or initialized above.
    ! Calculate now totfrac_nobio and soiltiles using these variables.
    
    ! Calculate totfrac_nobio
    totfrac_nobio(:) = zero
    DO jv = 1, nnobio
       totfrac_nobio(:) = totfrac_nobio(:) + frac_nobio(:,jv)
    ENDDO
    
    ! Calculate soiltile. This variable do not need to be in the restart file.
    ! The sum of all soiltiles makes one, and corresponds to the bio fraction
    ! of the grid cell (called vegtot in hydrol)
    soiltile(:,:) = zero
    DO jv = 1, nvm
       jst = pref_soil_veg(jv)
       DO ji = 1, kjpindex
          soiltile(ji,jst) = soiltile(ji,jst) + veget_max(ji,jv)
       ENDDO
    ENDDO
    DO ji = 1, kjpindex 
       IF (totfrac_nobio(ji) .LT. (1-min_sechiba)) THEN
          soiltile(ji,:)=soiltile(ji,:)/(1-totfrac_nobio(ji))
       ENDIF
    ENDDO
    
    ! Always calculate tot_bare_soil
    ! Fraction of bare soil in the mesh (bio+nobio) 
    tot_bare_soil(:) = veget_max(:,1)
    DO jv = 2, nvm
       DO ji =1, kjpindex
          tot_bare_soil(ji) = tot_bare_soil(ji) + (veget_max(ji,jv) - veget(ji,jv))
       ENDDO
    END DO
    
    !! 6. Verify consistency between different fractions. No change of the variables.
    IF (ok_stomate .AND. .NOT.disable_fire) THEN
        !spitfire
        !    
        !Config  Key  = LIGHTNING
        !Config  Desc = Read the ligntning map
        !Config Def  = ?
        !Config If   = NOT FIRE_DISABLE
        !Config Help = reads a 12 month lightning map which will
        !Config        then be interpolated to daily values as needed.
        !    
        data_filename='lightn.nc'
        CALL setvar_p (m_lightn, val_exp, 'm_lightn', m_lightn_default)
        CALL getin_p('LIGHTNING_FILE',data_filename)
        CALL slowproc_read_data(kjpindex, lalo, resolution, m_lightn,data_filename,'lightn')

        !    
        !Config  Key  = LCC_MATRIX
        !Config  Desc = Read the ligntning map
        !Config Def  = ?
        !Config Help = reads a 12 month lightning map which will
        !Config        then be interpolated to daily values as needed.
        !    

        CALL setvar_p (glccNetLCC, val_exp, 'glccNetLCC', 0.)
        CALL setvar_p (glccSecondShift, val_exp, 'glccSecondShift', 0.)
        CALL setvar_p (glccPrimaryShift, val_exp, 'glccPrimaryShift', 0.)
        CALL setvar_p (harvest_matrix, val_exp, 'harvest_matrix', 0.)
        CALL setvar_p (harvest_biomass, val_exp, 'harvest_biomass', 0.)
        IF ( (use_age_class) .AND. (veget_update .GT. 0) ) THEN
          data_filename = 'GLUC_NET_LCC_FILE.nc'
          CALL getin_p('GLUC_NET_LCC_FILE',data_filename)
          CALL slowproc_read_data(kjpindex, lalo, resolution, glccNetLCC, data_filename, 'matrix')

          data_filename = 'GLUC_SHIFT_SEC_FILE.nc'
          CALL getin_p('GLUC_SHIFT_SEC_FILE',data_filename)
          CALL slowproc_read_data(kjpindex, lalo, resolution, glccSecondShift, data_filename, 'matrix')

          data_filename = 'GLUC_SHIFT_PRI_FILE.nc'
          CALL getin_p('GLUC_SHIFT_PRI_FILE',data_filename)
          CALL slowproc_read_data(kjpindex, lalo, resolution, glccPrimaryShift, data_filename, 'matrix')

          data_filename = 'GLUC_NewVegFrac_File.nc'
          IF (gluc_newfrac_guide) THEN
            CALL getin_p('GLUC_NewVegFrac_File',data_filename)
            CALL slowproc_read_data(kjpindex, lalo, resolution, newvegfrac, data_filename, 'value')
          ENDIF

          IF (allow_forestry_harvest) THEN
            data_filename = 'GLUC_FORESTRY_HARVEST_FILE.nc'
            CALL getin_p('GLUC_FORESTRY_HARVEST_FILE',data_filename)
            CALL slowproc_read_data(kjpindex, lalo, resolution, harvest_matrix, data_filename, 'matrix')

            data_filename = 'GLUC_HARVEST_BIOMASS_FILE.nc'
            IF (gluc_use_harvest_biomass) THEN
              CALL getin_p('GLUC_HARVEST_BIOMASS_FILE',data_filename)
              ! harvest_matrix should have the 1st dim as industrial wood havest, 2nd as fuel wood, and 3rd
              ! dim as the fuel wood fraction. The unit should be GgC.
              CALL slowproc_read_data(kjpindex, lalo, resolution, harvest_biomass, data_filename, 'matrix')
              harvest_biomass(:,1:2) = harvest_biomass(:,1:2) * 1e9 !change to gC
            ENDIF
          ENDIF

          IF (gluc_allow_trans_bioe) THEN
            data_filename = 'GLUC_TRANS_BIOE.nc'
            !CALL getin_p('GLUC_TRANS_BIOE1_FILE',data_filename)
            !CALL slowproc_read_data(kjpindex, lalo, resolution, trans_bioe1_matrix, data_filename, 'matrix')
          ENDIF

        ENDIF

        IF ( use_age_class .AND. use_bound_spa ) THEN
          CALL getin_p('GLUC_AGE_THRESHOLD_FILE',data_filename)
          CALL slowproc_read_data(kjpindex, lalo, resolution, bound_spa, data_filename, 'value')
        ENDIF

        !Config  Key  = proxy_anidens 
        !Config  Desc = Read the general animal_density map
        !Config Def  = ?
        !Config Help = reads a 12 month lightning map which will
        !Config        then be interpolated to daily values as needed.
        !    


        !Config  Key  = OBSERVED_BA_FLAG 
        !Config  Desc = Read the observed burned_area flag
        !Config Def  = ?
        !    
        read_observed_ba = .FALSE.
        CALL getin_p('READ_OBSERVED_BA',read_observed_ba)
        WRITE(numout,*) 'flag for READ_OBSERVED_BA ', read_observed_ba

        !    
        !Config  Key  = OBSERVED_BURNED_AREA 
        !Config  Desc = Read the population density map
        !Config Def  = ?
        !Config Help = reads a one-year monthly burned area map with 12 as the value
        !of time dimension 
        !    
        IF(read_observed_ba) THEN
          CALL setvar_p (m_observed_ba, val_exp, 'm_observed_ba', m_ba_default)
          CALL getin_p('BA_FILE',data_filename)
          CALL slowproc_read_data(kjpindex, lalo, resolution, m_observed_ba,data_filename,'ba')
        ENDIF

        ! 1.
        !Config  Key  = OBSERVED_BA_FLAG 
        !Config  Desc = Read the observed burned_area flag
        !Config Def  = ?
        !    
        read_cf_coarse = .FALSE.
        CALL getin_p('READ_CF_COARSE',read_cf_coarse)
        WRITE(numout,*) 'flag for READ_CF_COARSE ', read_cf_coarse

        !    
        !Config  Key  = OBSERVED_BURNED_AREA 
        !Config  Desc = Read the population density map
        !Config Def  = ?
        !Config Help = reads a one-year monthly burned area map with 12 as the value
        !of time dimension 
        !    
        IF(read_cf_coarse) THEN
          CALL setvar_p (m_cf_coarse, val_exp, 'm_cf_coarse', m_cf_coarse_default)
          CALL getin_p('CF_COARSE_FILE',data_filename)
          CALL slowproc_read_data(kjpindex, lalo, resolution, m_cf_coarse,data_filename,'cf')
        ENDIF

        ! 2.
        !Config  Key  = OBSERVED_BA_FLAG 
        !Config  Desc = Read the observed burned_area flag
        !Config Def  = ?
        !    
        read_cf_fine = .FALSE.
        CALL getin_p('READ_CF_FINE',read_cf_fine)
        WRITE(numout,*) 'flag for READ_CF_FINE ', read_cf_fine

        !    
        !Config  Key  = OBSERVED_BURNED_AREA 
        !Config  Desc = Read the population density map
        !Config Def  = ?
        !Config Help = reads a one-year monthly burned area map with 12 as the value
        !of time dimension 
        !    
        IF(read_cf_fine) THEN
          CALL setvar_p (m_cf_fine, val_exp, 'm_cf_fine', m_cf_fine_default)
          CALL getin_p('CF_FINE_FILE',data_filename)
          CALL slowproc_read_data(kjpindex, lalo, resolution, m_cf_fine,data_filename,'cf')
        ENDIF

        ! 3.
        !Config  Key  = OBSERVED_BA_FLAG 
        !Config  Desc = Read the observed burned_area flag
        !Config Def  = ?
        !    
        read_ratio = .TRUE.
        CALL getin_p('READ_RATIO',read_ratio)
        WRITE(numout,*) 'flag for READ_RATIO ', read_ratio

        !    
        !Config  Key  = OBSERVED_BURNED_AREA 
        !Config  Desc = Read the population density map
        !Config Def  = ?
        !Config Help = reads a one-year monthly burned area map with 12 as the value
        !of time dimension 
        !    
        IF(read_ratio) THEN
          CALL setvar_p (m_ratio, val_exp, 'm_ratio', m_ratio_default)
          CALL getin_p('RATIO_FILE',data_filename)
          CALL slowproc_read_data(kjpindex, lalo, resolution, m_ratio,data_filename,'ratio')
        ENDIF

        ! 4.
        !Config  Key  = OBSERVED_BA_FLAG 
        !Config  Desc = Read the observed burned_area flag
        !Config Def  = ?
        !    
        read_ratio_flag = .TRUE.
        CALL getin_p('READ_RATIO_FLAG',read_ratio_flag)
        WRITE(numout,*) 'flag for READ_RATIO_FLAG ', read_ratio_flag

        !    
        !Config  Key  = OBSERVED_BURNED_AREA 
        !Config  Desc = Read the population density map
        !Config Def  = ?
        !Config Help = reads a one-year monthly burned area map with 12 as the value
        !of time dimension 
        !    
        IF(read_ratio_flag) THEN
          CALL setvar_p (m_ratio_flag, val_exp, 'm_ratio_flag', m_ratio_flag_default)
          CALL getin_p('RATIO_FLAG_FILE',data_filename)
          CALL slowproc_read_data(kjpindex, lalo, resolution, m_ratio_flag,data_filename,'invalid_flag')
        ENDIF

        !Config  Key  = POPDENS_FLAG 
        !Config  Desc = Read the popdens flag
        !Config Def  = n
        !    
        read_popdens = .FALSE.
        CALL getin_p('READ_POPDENS',read_popdens)
        WRITE(numout,*) 'flag for READ_POPDENS ',read_popdens

        !    
        !Config  Key  = POPDENS 
        !Config  Desc = Read the population density map
        !Config Def  = ?
        !Config Help = reads a yearly map 
        !    
        popd=0.
        IF(read_popdens) THEN
          CALL setvar_p (popd, val_exp, 'popdens', popdens_default)
          CALL getin_p('POPDENS_FILE',data_filename)
          CALL slowproc_read_annual(kjpindex, lalo, resolution, popd,data_filename,'popdens')
        ENDIF

        !Config  Key  = HUMIGN_FLAG 
        !Config  Desc = Read the human ignition parameter flag
        !Config Def  = n
        !    
        read_humign = .FALSE.
        CALL getin_p('READ_HUMIGN',read_humign)
        WRITE(numout,*) 'flag for READ_HUMIGN:',read_humign

        !    
        !Config  Key  = POPDENS 
        !Config  Desc = Read the population density map
        !Config Def  = ?
        !Config Help = reads a yearly map 
        !    
        humign=0.22
        IF(read_humign) THEN
          CALL setvar_p (humign, val_exp, 'HUMIGN_FILE', 0.22)
          CALL getin_p('HUMIGN_FILE',data_filename)
          CALL slowproc_read_annual(kjpindex, lalo, resolution, humign,data_filename,'humign')
        ENDIF
        !endspit
    ENDIF !! NOT disable_fire

    IF (printlev_loc>=3) WRITE (numout,*) ' slowproc_init done '


  END SUBROUTINE slowproc_init

!! ================================================================================================================================
!! SUBROUTINE   : slowproc_clear
!!
!>\BRIEF          Clear all variables related to slowproc and stomate modules  
!!
!_ ================================================================================================================================

  SUBROUTINE slowproc_clear 

  ! 1 clear all the variables defined as common for the routines in slowproc 

    IF (ALLOCATED (clayfraction)) DEALLOCATE (clayfraction)
    IF (ALLOCATED (laimap)) DEALLOCATE (laimap)
    IF (ALLOCATED (frac_nobio_lastyear)) DEALLOCATE (frac_nobio_lastyear)
    IF (ALLOCATED (vegetnew_firstday)) DEALLOCATE (vegetnew_firstday)
    IF (ALLOCATED (veget_max_new)) DEALLOCATE (veget_max_new)
    IF (ALLOCATED (frac_nobio_new)) DEALLOCATE (frac_nobio_new)
    IF ( ALLOCATED (soilclass_default)) DEALLOCATE (soilclass_default)
    !spitfire
    IF (ALLOCATED(m_lightn)) DEALLOCATE (m_lightn)
    IF (ALLOCATED(glccNetLCC)) DEALLOCATE (glccNetLCC)
    IF (ALLOCATED(glccSecondShift)) DEALLOCATE (glccSecondShift)
    IF (ALLOCATED(glccPrimaryShift)) DEALLOCATE (glccPrimaryShift)
    IF (ALLOCATED(harvest_matrix)) DEALLOCATE (harvest_matrix)
    IF (ALLOCATED(harvest_biomass)) DEALLOCATE (harvest_biomass)
    IF (ALLOCATED(bound_spa)) DEALLOCATE (bound_spa)
    IF (ALLOCATED(newvegfrac)) DEALLOCATE (newvegfrac)
    IF (ALLOCATED(proxy_anidens)) DEALLOCATE (proxy_anidens)
    IF (ALLOCATED(popd)) DEALLOCATE (popd)
    IF (ALLOCATED(humign)) DEALLOCATE (humign)
    IF (ALLOCATED(m_observed_ba)) DEALLOCATE (m_observed_ba)
    IF (ALLOCATED(m_cf_coarse)) DEALLOCATE (m_cf_coarse)
    IF (ALLOCATED(m_cf_fine)) DEALLOCATE (m_cf_fine)
    IF (ALLOCATED(m_ratio)) DEALLOCATE (m_ratio)
    IF (ALLOCATED(m_ratio_flag)) DEALLOCATE (m_ratio_flag)
    !endspit
!gmjc fc_grazing
    IF (ALLOCATED (fc_grazing)) DEALLOCATE (fc_grazing)
!end gmjc
!!!qcj++ peatland
    IF (ALLOCATED (veget_cov_max_adjusted)) DEALLOCATE (veget_cov_max_adjusted)
    IF (ALLOCATED (peatPET_lastyear)) DEALLOCATE (peatPET_lastyear)
    IF (ALLOCATED (peatPET_thisyear)) DEALLOCATE (peatPET_thisyear)
    IF (ALLOCATED (summerp_long)) DEALLOCATE (summerp_long)
    IF (ALLOCATED (summerpet_long)) DEALLOCATE (summerpet_long)
    IF (ALLOCATED (growth_day)) DEALLOCATE (growth_day)
    IF (ALLOCATED (GSL)) DEALLOCATE (GSL)
!    IF (ALLOCATED (t2m_lastyear)) DEALLOCATE (t2m_lastyear)
!    IF (ALLOCATED (t2m_thisyear)) DEALLOCATE (t2m_thisyear)
    IF (ALLOCATED (precipitation_lastsummer)) DEALLOCATE (precipitation_lastsummer)
    IF (ALLOCATED (precipitation_thissummer)) DEALLOCATE (precipitation_thissummer)
    IF (ALLOCATED (peatC)) DEALLOCATE (peatC)
    IF (ALLOCATED (peatC_ok)) DEALLOCATE (peatC_ok)
    IF (ALLOCATED (fpeat_map)) DEALLOCATE (fpeat_map)

 ! 2. Clear all the variables in stomate 

    CALL stomate_clear 
    !
  END SUBROUTINE slowproc_clear

!! ================================================================================================================================
!! SUBROUTINE   : slowproc_derivvar
!!
!>\BRIEF         Initializes variables related to the
!! parameters to be assimilated, the maximum water on vegetation, the vegetation height, 
!! and the fraction of soil covered by dead leaves and the vegetation height 
!!
!! DESCRIPTION  : (definitions, functional, design, flags):
!! (1) Initialization of the variables relevant for the assimilation parameters  
!! (2) Intialization of the fraction of soil covered by dead leaves
!! (3) Initialization of the Vegetation height per PFT
!! (3) Initialization the maximum water on vegetation for interception with a particular treatement of the PFT no.1
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::qsintmax, ::deadleaf_cover, ::assim_param, ::height  
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_derivvar (kjpindex, veget, lai, &
       qsintmax, deadleaf_cover, assim_param, height, temp_growth)

    !! INTERFACE DESCRIPTION

    !! 0.1 Input scalar and fields 
    INTEGER(i_std), INTENT (in)                                :: kjpindex       !! Domain size - terrestrial pixels only
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)          :: veget          !! Fraction of pixel covered by PFT. Fraction accounts for none-biological land covers (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)          :: lai            !! PFT leaf area index (m^{2} m^{-2})

    !! 0.2. Output scalar and fields 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)          :: qsintmax       !! Maximum water on vegetation for interception(mm)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)              :: deadleaf_cover !! fraction of soil covered by dead leaves (unitless)
    REAL(r_std), DIMENSION (kjpindex,nvm,npco2), INTENT (out)   :: assim_param    !! min+max+opt temperatures & vmax for photosynthesis (K, \mumol m^{-2} s^{-1})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)          :: height         !! height of the vegetation or surface in general ??? (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)              :: temp_growth    !! growth temperature (�C)  
    !
    !! 0.3 Local declaration
    INTEGER(i_std)                                              :: jv             !! Local indices
!_ ================================================================================================================================

    !
    ! 1. Initialize (why here ??) the variables revelant for the assimilation parameters
    !
    DO jv = 1, nvm
       assim_param(:,jv,ivcmax) = vcmax_fix(jv)
    ENDDO

    !
    ! 2. Intialize the fraction of soil covered by dead leaves 
    !
    deadleaf_cover(:) = zero

    !
    ! 3. Initialize the Vegetation height per PFT
    !
    DO jv = 1, nvm
       height(:,jv) = height_presc(jv)
    ENDDO
    !
    ! 4. Initialize the maximum water on vegetation for interception
    !
    qsintmax(:,:) = qsintcst * veget(:,:) * lai(:,:)

    ! Added by Nathalie - July 2006
    !  Initialize the case of the PFT no.1 to zero 
    qsintmax(:,1) = zero

    temp_growth(:)=25.

  END SUBROUTINE slowproc_derivvar


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_mean
!!
!>\BRIEF          Accumulates field_in over a period of dt_tot.
!! Has to be called at every time step (dt). 
!! Mean value is calculated if ldmean=.TRUE.
!! field_mean must be initialized outside of this routine! 
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!! (1) AcumAcuumlm 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::field_main
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_mean (npts, n_dim2, dt_tot, dt, ldmean, field_in, field_mean)

    !
    !! 0 declarations

    !! 0.1 input scalar and variables 
    INTEGER(i_std), INTENT(in)                           :: npts     !! Domain size- terrestrial pixels only 
    INTEGER(i_std), INTENT(in)                           :: n_dim2   !! Number of PFTs 
    REAL(r_std), INTENT(in)                              :: dt_tot   !! Time step of stomate (in days). The period over which the accumulation or the mean is computed 
    REAL(r_std), INTENT(in)                              :: dt       !! Time step in days 
    LOGICAL, INTENT(in)                                  :: ldmean   !! Flag to calculate the mean after the accumulation ???
    REAL(r_std), DIMENSION(npts,n_dim2), INTENT(in)      :: field_in !! Daily field 

    !! 0.3 Modified field; The computed sum or mean field over dt_tot time period depending on the flag ldmean 
    REAL(r_std), DIMENSION(npts,n_dim2), INTENT(inout)   :: field_mean !! Accumulated field at dt_tot time period or mean field over dt_tot 
 

!_ ================================================================================================================================

    !
    ! 1. Accumulation the field over dt_tot period 
    !
    field_mean(:,:) = field_mean(:,:) + field_in(:,:) * dt

    !
    ! 2. If the flag ldmean set, the mean field is computed over dt_tot period  
    !
    IF (ldmean) THEN
       field_mean(:,:) = field_mean(:,:) / dt_tot
    ENDIF

  END SUBROUTINE slowproc_mean


  
!! ================================================================================================================================
!! SUBROUTINE   : slowproc_long
!!
!>\BRIEF        Calculates a temporally smoothed field (field_long) from
!! instantaneous input fields.Time constant tau determines the strength of the smoothing.
!! For tau -> infinity??, field_long becomes the true mean value of field_inst
!! (but  the spinup becomes infinietly long, too).
!! field_long must be initialized outside of this routine! 
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!! (1) Testing the time coherence betwen the time step dt and the time tau over which
!! the rescaled of the mean is performed   
!!  (2) Computing the rescaled mean over tau period 
!! MAIN OUTPUT VARIABLE(S): field_long  
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::field_long
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_long (npts, n_dim2, dt, tau, field_inst, field_long)

    !
    ! 0 declarations
    !

    ! 0.1 input scalar and fields 

    INTEGER(i_std), INTENT(in)                                 :: npts        !! Domain size- terrestrial pixels only
    INTEGER(i_std), INTENT(in)                                 :: n_dim2      !! Second dimension of the fields, which represents the number of PFTs
    REAL(r_std), INTENT(in)                                    :: dt          !! Time step in days   
    REAL(r_std), INTENT(in)                                    :: tau         !! Integration time constant (has to have same unit as dt!)  
    REAL(r_std), DIMENSION(npts,n_dim2), INTENT(in)            :: field_inst  !! Instantaneous field 


    ! 0.2 modified field

    ! Long-term field
    REAL(r_std), DIMENSION(npts,n_dim2), INTENT(inout)         :: field_long  !! Mean value of the instantaneous field rescaled at tau time period 

!_ ================================================================================================================================

    !
    ! 1 test coherence of the time 

    IF ( ( tau .LT. dt ) .OR. ( dt .LE. zero ) .OR. ( tau .LE. zero ) ) THEN
       WRITE(numout,*) 'slowproc_long: Problem with time steps'
       WRITE(numout,*) 'dt=',dt
       WRITE(numout,*) 'tau=',tau
    ENDIF

    !
    ! 2 integration of the field over tau 

    field_long(:,:) = ( field_inst(:,:)*dt + field_long(:,:)*(tau-dt) ) / tau

  END SUBROUTINE slowproc_long


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_veget
!!
!>\BRIEF        Set small fractions to zero and normalize to keep the sum equal 1. Calucate veget and soiltile.
!!
!! DESCRIPTION  : Set small fractions to zero and normalize to keep the sum equal 1. Calucate veget and soiltile.
!! (1) Set veget_max and frac_nobio for fraction smaller than min_vegfrac.
!! (2) Reset some variables in stomate for small fractions
!! (3) Calculate veget
!! (5) Calculate totfrac_nobio
!! (6) Calculate soiltile
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: frac_nobio, totfrac_nobio, veget_max, veget, soiltile
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_veget (kjpindex, f_rot_sech, lai, frac_nobio, totfrac_nobio, veget_max, veget, soiltile)
    !
    ! 0. Declarations
    !
    ! 0.1 Input variables 
    INTEGER(i_std), INTENT(in)                             :: kjpindex    !! Domain size - terrestrial pixels only
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)       :: lai         !! PFT leaf area index (m^{2} m^{-2})

    ! 0.2 Modified variables 
    REAL(r_std), DIMENSION(kjpindex,nnobio), INTENT(inout) :: frac_nobio  !! Fraction of the mesh which is covered by ice, lakes, ...
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(inout)    :: veget_max   !! Maximum fraction of vegetation type including none biological fraction (unitless)

    ! 0.3 Output variables 
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)      :: veget       !! Fraction of pixel covered by PFT. Fraction accounts for none-biological land covers (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: totfrac_nobio
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT(out)    :: soiltile     !! Fraction of each soil tile within vegtot (0-1, unitless)
    LOGICAL,DIMENSION(kjpindex), INTENT(in)         :: f_rot_sech        !! whether a grid point is under rotation

    ! 0.4 Local scalar and varaiables 
    INTEGER(i_std)                                         :: ji, jv, jst !! indices 
    REAL(r_std)                                            :: SUMveg      

!_ ================================================================================================================================
    IF (printlev_loc > 8) WRITE(numout,*) 'Entering slowproc_veget'

    ! 0. Normalize fractions of frac_nobio and veget_max smaller than min_vegfrac
    !    This is due to precision issues (float 64). It might lead to negative values. 
    !    At some point, a floating point exception.
    DO ji = 1, kjpindex
       IF ( SUM(frac_nobio(ji,:)) .LT. min_vegfrac ) THEN
          frac_nobio(ji,:) = zero
       ENDIF
    
       IF (.NOT. ok_dgvm) THEN
          IF (.NOT. ok_dgvm_peat) THEN
             DO jv = 1, nvm
               IF ( veget_max(ji,jv) .LT. min_vegfrac ) THEN
                 veget_max(ji,jv) = zero
               ENDIF
             ENDDO
          ELSE
!!!qcj++peatland
             DO jv = 1, nvm
               IF ( (veget_max(ji,jv) .LT. min_vegfrac) .AND. ( .NOT. is_peat(jv)) ) THEN
                 veget_max(ji,jv) = zero
               ENDIF
             ENDDO
          ENDIF
       END IF
 
       !! Normalize to keep the sum equal 1.
       SUMveg = SUM(frac_nobio(ji,:))+SUM(veget_max(ji,:))
       frac_nobio(ji,:) = frac_nobio(ji,:)/SUMveg
       veget_max(ji,:) = veget_max(ji,:)/SUMveg
    ENDDO

    !! 3. Calculate veget
    !!    If lai of a vegetation type (jv > 1) is small, increase soil part
    !!    stomate-like calculation
    DO ji = 1, kjpindex
       veget(ji,1)=veget_max(ji,1)
       DO jv = 2, nvm
          veget(ji,jv) = veget_max(ji,jv) * ( un - exp( - lai(ji,jv) * ext_coeff_vegetfrac(jv) ) )
       ENDDO
    ENDDO


    !! 4. Calculate totfrac_nobio
    totfrac_nobio(:) = zero
    DO jv = 1, nnobio
       totfrac_nobio(:) = totfrac_nobio(:) + frac_nobio(:,jv)
    ENDDO
    

    !! 5. Calculate soiltiles
    !! Soiltiles are only used in hydrol, but we fix them in here because some time it might depend
    !! on a changing vegetation (but then some adaptation should be made to hydrol) and be also used
    !! in the other modules to perform separated energy balances
    ! The sum of all soiltiles makes one, and corresponds to the bio fraction
    ! of the grid cell (called vegtot in hydrol)   
    soiltile(:,:) = zero
    DO jv = 1, nvm
       jst = pref_soil_veg(jv)
       DO ji = 1, kjpindex
          soiltile(ji,jst) = soiltile(ji,jst) + veget_max(ji,jv)
       ENDDO
    ENDDO
    DO ji = 1, kjpindex 
       IF (totfrac_nobio(ji) .LT. (1-min_sechiba)) THEN
          soiltile(ji,:)=soiltile(ji,:)/(1.-totfrac_nobio(ji))
       ENDIF
    ENDDO   

  END SUBROUTINE slowproc_veget
 
 
!! ================================================================================================================================
!! SUBROUTINE   : slowproc_lai
!!
!>\BRIEF        Do the interpolation of lai for the PFTs in case the laimap is not read   
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!! (1) Interplation by using the mean value of laimin and laimax for the PFTs    
!! (2) Interpolation between laimax and laimin values by using the temporal
!!  variations 
!! (3) If problem occurs during the interpolation, the routine stops 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::lai
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_lai (kjpindex,lcanop,stempdiag,lalo,resolution,lai,laimap)
    !
    ! 0. Declarations
    !
    !! 0.1 Input variables 
    INTEGER(i_std), INTENT(in)                          :: kjpindex   !! Domain size - terrestrial pixels only
    INTEGER(i_std), INTENT(in)                          :: lcanop     !! soil level used for LAI
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)  :: stempdiag  !! Soil temperature (K) ???
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)     :: lalo       !! Geogr. coordinates (latitude,longitude) (degrees)
    REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)     :: resolution !! Size in x an y of the grid (m) - surface area of the gridbox
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)           :: laimap     !! map of lai read 

    !! 0.2 Output
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)   :: lai        !! PFT leaf area index (m^{2} m^{-2})LAI

    !! 0.4 Local
    INTEGER(i_std)                                      :: ji,jv      !! Local indices 
!_ ================================================================================================================================

    !
    IF  ( .NOT. read_lai ) THEN
    
       lai(: ,1) = zero
       ! On boucle sur 2,nvm au lieu de 1,nvm
       DO jv = 2,nvm
          SELECT CASE (type_of_lai(jv))
             
          CASE ("mean ")
             !
             ! 1. do the interpolation between laimax and laimin
             !
             lai(:,jv) = undemi * (llaimax(jv) + llaimin(jv))
             !
          CASE ("inter")
             !
             ! 2. do the interpolation between laimax and laimin
             !
             DO ji = 1,kjpindex
                lai(ji,jv) = llaimin(jv) + tempfunc(stempdiag(ji,lcanop)) * (llaimax(jv) - llaimin(jv))
             ENDDO
             !
          CASE default
             !
             ! 3. Problem
             !
             WRITE (numout,*) 'This kind of lai choice is not possible. '// &
                  ' We stop with type_of_lai ',jv,' = ', type_of_lai(jv) 
             CALL ipslerr_p(3,'slowproc_lai','Bad value for type_of_lai','read_lai=false','')
          END SELECT
          
       ENDDO
       !
    ELSE
       lai(: ,1) = zero
       ! On boucle sur 2,nvm au lieu de 1,nvm
       DO jv = 2,nvm

          SELECT CASE (type_of_lai(jv))
             
          CASE ("mean ")
             !
             ! 1. force MAXVAL of laimap on lai on this PFT
             !
             DO ji = 1,kjpindex
                lai(ji,jv) = MAXVAL(laimap(ji,jv,:))
             ENDDO
             !
          CASE ("inter")
             !
             ! 2. do the interpolation between laimax and laimin
             !
             !
             ! If January
             !
             IF (month_end .EQ. 1 ) THEN
                IF (day_end .LE. 15) THEN
                   lai(:,jv) = laimap(:,jv,12)*(1-(day_end+15)/30.) + laimap(:,jv,1)*((day_end+15)/30.)
                ELSE
                   lai(:,jv) = laimap(:,jv,1)*(1-(day_end-15)/30.) + laimap(:,jv,2)*((day_end-15)/30.)
                ENDIF
                !
                ! If December
                !
             ELSE IF (month_end .EQ. 12) THEN
                IF (day_end .LE. 15) THEN
                   lai(:,jv) = laimap(:,jv,11)*(1-(day_end+15)/30.) + laimap(:,jv,12)*((day_end+15)/30.)
                ELSE
                   lai(:,jv) = laimap(:,jv,12)*(1-(day_end-15)/30.) + laimap(:,jv,1)*((day_end-15)/30.)
                ENDIF
          !
          ! ELSE
          !
             ELSE
                IF (day_end .LE. 15) THEN
                   lai(:,jv) = laimap(:,jv,month_end-1)*(1-(day_end+15)/30.) + laimap(:,jv,month_end)*((day_end+15)/30.)
                ELSE
                   lai(:,jv) = laimap(:,jv,month_end)*(1-(day_end-15)/30.) + laimap(:,jv,month_end+1)*((day_end-15)/30.)
                ENDIF
             ENDIF
             !
          CASE default
             !
             ! 3. Problem
             !
             WRITE (numout,*) 'This kind of lai choice is not possible. '// &
                  ' We stop with type_of_lai ',jv,' = ', type_of_lai(jv) 
             CALL ipslerr_p(3,'slowproc_lai','Bad value for type_of_lai','read_lai=true','')
          END SELECT
          
       ENDDO
    ENDIF

  END SUBROUTINE slowproc_lai

!! ================================================================================================================================
!! SUBROUTINE   : slowproc_interlai
!!
!>\BRIEF         Interpolate the LAI map to the grid of the model 
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::laimap
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_interlai(nbpt, lalo, resolution, neighbours, contfrac, laimap)

    USE interpweight

    IMPLICIT NONE

    !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)          :: nbpt                  !! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)             :: lalo(nbpt,2)          !! Vector of latitude and longitudes 
                                                                 !! (beware of the order = 1 : latitude, 2 : longitude)
    REAL(r_std), INTENT(in)             :: resolution(nbpt,2)    !! The size in km of each grid-box in X and Y
    INTEGER(i_std), INTENT(in)          :: neighbours(nbpt,NbNeighb)!! Vector of neighbours for each grid point
                                                                 !! (1=North and then clockwise)
    REAL(r_std), INTENT(in)             :: contfrac(nbpt)        !! Fraction of land in each grid box.
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), INTENT(out)    ::  laimap(nbpt,nvm,12)          !! lai read variable and re-dimensioned
    !
    !  0.3 LOCAL
    !
    CHARACTER(LEN=80) :: filename                               !! name of the LAI map read
    INTEGER(i_std) :: ib, ip, jp, it, jv
    REAL(r_std) :: lmax, lmin, ldelta
    LOGICAL ::           renormelize_lai  ! flag to force LAI renormelization
    INTEGER                  :: ier

    REAL(r_std), DIMENSION(nbpt)                         :: alaimap          !! availability of the lai interpolation 
    INTEGER, DIMENSION(4)                                :: invardims
    REAL(r_std), DIMENSION(:,:,:), ALLOCATABLE           :: lairefrac        !! lai fractions re-dimensioned
    REAL(r_std), DIMENSION(:), ALLOCATABLE               :: vmin, vmax       !! min/max values to use for the 
                                                                             !!   renormalization
    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat names in input file
    REAL(r_std), DIMENSION(nvm)                          :: variabletypevals !! Values for all the types of the variable
                                                                             !!   (variabletypevals(1) = -un, not used)
    CHARACTER(LEN=50)                                    :: fractype         !! method of calculation of fraction
                                                                             !!   'XYKindTime': Input values are kinds 
                                                                             !!     of something with a temporal 
                                                                             !!     evolution on the dx*dy matrix'
    LOGICAL                                              :: nonegative       !! whether negative values should be removed
    CHARACTER(LEN=50)                                    :: maskingtype      !! Type of masking
                                                                             !!   'nomask': no-mask is applied
                                                                             !!   'mbelow': take values below maskvals(1)
                                                                             !!   'mabove': take values above maskvals(1)
                                                                             !!   'msumrange': take values within 2 ranges;
                                                                             !!      maskvals(2) <= SUM(vals(k)) <= maskvals(1)
                                                                             !!      maskvals(1) < SUM(vals(k)) <= maskvals(3)
                                                                             !!        (normalized by maskvals(3))
                                                                             !!   'var': mask values are taken from a 
                                                                             !!     variable inside the file (>0)
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to 
                                                                             !!   `maskingtype') 
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask 
!_ ================================================================================================================================

    !
    !Config Key   = LAI_FILE
    !Config Desc  = Name of file from which the vegetation map is to be read
    !Config If    = LAI_MAP
    !Config Def   = lai2D.nc
    !Config Help  = The name of the file to be opened to read the LAI
    !Config         map is to be given here. Usualy SECHIBA runs with a 5kmx5km
    !Config         map which is derived from a Nicolas VIOVY one. 
    !Config Units = [FILE]
    !
    filename = 'lai2D.nc'
    CALL getin_p('LAI_FILE',filename)

    variablename = 'LAI'

    IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_interlai: Read and interpolate " &
         // TRIM(filename) //" for variable " //TRIM(variablename)

    ! invardims: shape of variable in input file to interpolate
    invardims = interpweight_get_var4dims_file(filename, variablename)
    ! Check coherence of dimensions read from the file
    IF (invardims(4) /= 12)  CALL ipslerr_p(3,'slowproc_interlai','Wrong dimension of time dimension in input file for lai','','')
    IF (invardims(3) /= nvm) CALL ipslerr_p(3,'slowproc_interlai','Wrong dimension of PFT dimension in input file for lai','','')

    ALLOCATE(vmin(nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_interlai','Problem in allocation of variable vmin','','')

    ALLOCATE(vmax(nvm), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_interlai','Problem in allocation of variable vmax','','')

    ALLOCATE(lairefrac(nbpt,nvm,invardims(4)), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_interlai','Problem in allocation of variable lairefrac','','')

! Assigning values to vmin, vmax
    vmin = un
    vmax = nvm*un

    variabletypevals = -un

    !! Variables for interpweight
    ! Type of calculation of cell fractions
    fractype = 'default'
    ! Name of the longitude and latitude in the input file
    lonname = 'longitude'
    latname = 'latitude'
    ! Should negative values be set to zero from input file?
    nonegative = .TRUE.
    ! Type of mask to apply to the input data (see header for more details)
    maskingtype = 'mbelow'
    ! Values to use for the masking
    maskvals = (/ 20., undef_sechiba, undef_sechiba /)
    ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
    namemaskvar = ''

    CALL interpweight_4D(nbpt, nvm, variabletypevals, lalo, resolution, neighbours,        &
      contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
      maskvals, namemaskvar, nvm, invardims(4), -1, fractype,                            &
      -1., -1., lairefrac, alaimap)

    IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_interlai after interpweight_4D'

    !
    !
    !Config Key   = RENORM_LAI
    !Config Desc  = flag to force LAI renormelization
    !Config If    = LAI_MAP
    !Config Def   = n
    !Config Help  = If true, the laimap will be renormalize between llaimin and llaimax parameters.
    !Config Units = [FLAG]
    !
    renormelize_lai = .FALSE.
    CALL getin_p('RENORM_LAI',renormelize_lai)

    !
    laimap(:,:,:) = zero
    !
    IF (printlev_loc >= 5) THEN
      WRITE(numout,*)'  slowproc_interlai before starting loop nbpt:', nbpt
    END IF 

    ! Assigning the right values and giving a value where information was not found
    DO ib=1,nbpt
      IF (alaimap(ib) < 0.) THEN
        DO jv=1,nvm
          laimap(ib,jv,:) = (llaimax(jv)+llaimin(jv))/deux
        ENDDO
      ELSE
        DO jv=1, nvm
          DO it=1, invardims(4)
            laimap(ib,jv,it) = lairefrac(ib,jv,it)
          ENDDO
        ENDDO
      END IF
    ENDDO
    !
    ! Normelize the read LAI by the values SECHIBA is used to
    !
    IF ( renormelize_lai ) THEN
       DO ib=1,nbpt
          DO jv=1, nvm
             lmax = MAXVAL(laimap(ib,jv,:))
             lmin = MINVAL(laimap(ib,jv,:))
             ldelta = lmax-lmin
             IF ( ldelta < min_sechiba) THEN
                ! LAI constante ... keep it constant
                laimap(ib,jv,:) = (laimap(ib,jv,:)-lmin)+(llaimax(jv)+llaimin(jv))/deux
             ELSE
                laimap(ib,jv,:) = (laimap(ib,jv,:)-lmin)/(lmax-lmin)*(llaimax(jv)-llaimin(jv))+llaimin(jv)
             ENDIF
          ENDDO
       ENDDO
    ENDIF

    ! Write diagnostics
    CALL xios_orchidee_send_field("alaimap",alaimap)
    
    IF (printlev_loc >= 3) WRITE(numout,*) '  slowproc_interlai ended'

  END SUBROUTINE slowproc_interlai

!! ================================================================================================================================
!! SUBROUTINE   : slowproc_readvegetmax
!!
!>\BRIEF          Interpolate a vegetation map (by pft)
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_readvegetmax(nbpt, lalo, neighbours,  resolution, contfrac, & 
       veget_last, veget_next, frac_nobio_next, veget_year, init)

    USE interpweight

    IMPLICIT NONE

    !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)                             :: nbpt            !! Number of points for which the data needs 
                                                                              !! to be interpolated
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: lalo            !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(nbpt,NbNeighb), INTENT(in)   :: neighbours      !! Vector of neighbours for each grid point
                                                                              !! (1=North and then clockwise)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: resolution      !! The size in km of each grid-box in X and Y
    REAL(r_std), DIMENSION(nbpt), INTENT(in)               :: contfrac        !! Fraction of continent in the grid
    !
    REAL(r_std), DIMENSION(nbpt,nvm), INTENT(in)           :: veget_last      !! old max vegetfrac
    INTEGER(i_std), INTENT(in)         :: veget_year            !! first year for landuse (0 == NO TIME AXIS)
    LOGICAL, INTENT(in)                :: init                  !! initialisation
                                                                !! In case of dgvm == FALSE, whatever its value, 
                                                                !! all PFT fractions will be updated.
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), DIMENSION(nbpt,nvmap), INTENT(out)          :: veget_next       !! new max vegetfrac
    REAL(r_std), DIMENSION(nbpt,nnobio), INTENT(out)       :: frac_nobio_next  !! new fraction of the mesh which is 
                                                                               !! covered by ice, lakes, ...
    
    !
    !  0.3 LOCAL
    !
    !
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: ib, inobio, jv
    REAL(r_std) :: sumf, err, norm
    !
    ! for DGVM case :
    REAL(r_std)                 :: sum_veg                     ! sum of vegets
    REAL(r_std)                 :: sum_nobio                   ! sum of nobios
    REAL(r_std)                 :: sumvAnthro_old, sumvAnthro  ! last an new sum of antrhopic vegets
    REAL(r_std)                 :: rapport                     ! (S-B) / (S-A)
    LOGICAL                     :: partial_update              ! if TRUE, partialy update PFT (only anthropic ones) 
                                                               ! e.g. in case of DGVM and not init (optional parameter)
    REAL(r_std), DIMENSION(nbpt,nvm)                     :: vegetrefrac      !! veget fractions re-dimensioned
    REAL(r_std), DIMENSION(nbpt)                         :: aveget           !! Availability of the soilcol interpolation
    REAL(r_std), DIMENSION(nvm)                          :: vmin, vmax       !! min/max values to use for the renormalization
    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat names in input file
    REAL(r_std), DIMENSION(nvm)                          :: variabletypevals !! Values for all the types of the variable
                                                                             !!   (variabletypevals(1) = -un, not used)
    CHARACTER(LEN=50)                                    :: fractype         !! method of calculation of fraction
                                                                             !!   'XYKindTime': Input values are kinds 
                                                                             !!     of something with a temporal 
                                                                             !!     evolution on the dx*dy matrix'
    LOGICAL                                              :: nonegative       !! whether negative values should be removed
    CHARACTER(LEN=50)                                    :: maskingtype      !! Type of masking
                                                                             !!   'nomask': no-mask is applied
                                                                             !!   'mbelow': take values below maskvals(1)
                                                                             !!   'mabove': take values above maskvals(1)
                                                                             !!   'msumrange': take values within 2 ranges;
                                                                             !!      maskvals(2) <= SUM(vals(k)) <= maskvals(1)
                                                                             !!      maskvals(1) < SUM(vals(k)) <= maskvals(3)
                                                                             !!        (normalized by maskvals(3))
                                                                             !!   'var': mask values are taken from a 
                                                                             !!     variable inside the file (>0)
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to 
                                                                             !!   `maskingtype') 
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask 
    CHARACTER(LEN=250)                                   :: msg

!_ ================================================================================================================================

    IF (printlev_loc >= 5) PRINT *,'  In slowproc_readvegetmax'

    !
    !Config Key   = VEGETATION_FILE
    !Config Desc  = Name of file from which the vegetation map is to be read
    !Config If    = MAP_PFT_FORMAT
    !Config Def   = PFTmap.nc
    !Config Help  = The name of the file to be opened to read a vegetation
    !Config         map (in pft) is to be given here. 
    !Config Units = [FILE]
    !
    filename = 'PFTmap.nc'
    CALL getin_p('VEGETATION_FILE',filename)
    variablename = 'maxvegetfrac'

    IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_readvegetmax: Read and interpolate " &
         // TRIM(filename) // " for variable " // TRIM(variablename)

    ! Assigning values to vmin, vmax
    vmin = 1
    ! chaoyue: be careful of this when do merging with trunk. This has to be
    ! persistently diverged from trunk because trunk does not have, neither uses
    ! the concept of age calsses. Here nvmap is the number of metaclasses.
    vmax = nvmap*1._r_std

    variabletypevals = -un

    !! Variables for interpweight
    ! Type of calculation of cell fractions
    fractype = 'default'
    ! Name of the longitude and latitude in the input file
    lonname = 'lon'
    latname = 'lat'
    ! Should negative values be set to zero from input file?
    nonegative = .FALSE.
    ! Type of mask to apply to the input data (see header for more details)
    maskingtype = 'msumrange'
    ! Values to use for the masking
    maskvals = (/ 1.-1.e-7, 0., 2. /)
    ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
    namemaskvar = ''

    ! persistent divergence with trunk. nvmap rather than nvm should be used.
    ! see the comment above.
    CALL interpweight_3D(nbpt, nvmap, variabletypevals, lalo, resolution, neighbours,        &
      contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
      maskvals, namemaskvar, nvmap, 0, veget_year, fractype,                                 &
      -1., -1., vegetrefrac, aveget)
    IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_readvegetmax after interpeeight_3D'
    
    !
    ! Compute the logical for partial (only anthropic) PTFs update
    IF (ok_dgvm .AND. .NOT. init) THEN
       partial_update= .TRUE.
    ELSE
       partial_update=.FALSE.
    END IF

    IF (printlev_loc >= 5) THEN
      WRITE(numout,*)'  slowproc_readvegetmax before updating loop nbpt:', nbpt
    END IF

    IF ( .NOT. partial_update ) THEN
       veget_next(:,:)=zero
       
       IF (printlev_loc >=3 .AND. ANY(aveget < min_sechiba)) THEN
          WRITE(numout,*) 'Some grid cells on the model grid did not have any points on the source grid.'
          IF (init) THEN
             WRITE(numout,*) 'Initialization with full fraction of bare soil are done for the below grid cells.'
          ELSE
             WRITE(numout,*) 'Old values are kept for the below grid cells.'
          ENDIF
          WRITE(numout,*) 'List of grid cells (ib, lat, lon):'
       END IF
 
      DO ib = 1, nbpt
          ! vegetrefrac is already normalized to sum equal one for each grid cell
          veget_next(ib,:) = vegetrefrac(ib,:)

          IF (aveget(ib) < min_sechiba) THEN
             IF (printlev_loc >=3) WRITE(numout,*) ib,lalo(ib,1),lalo(ib,2)
             IF (init) THEN
                veget_next(ib,1) = un
                veget_next(ib,2:nvmap) = zero
             ELSE
                veget_next(ib,:) = veget_last(ib,:)
             ENDIF
          ENDIF
       ENDDO
    ! `partial_update` is TRUE, the natural PFT fraction will be copied from 
    ! `veget_last` to `veget_next`, only the anthropogenic PFT are updated
    ! from the input veget_max map. This is used only when DGVM is activated.
    !print *,"second time within slowproc_readvegetmax"
    !print *,"veget_last",veget_last
    !print *,"veget_next",veget_next
    !print *,"partial_update",partial_update
    !print *,"vegmap",vegmap
    ELSE
       ! Partial update
       DO ib = 1, nbpt
          IF (aveget(ib) > min_sechiba) THEN
             ! For the case with properly interpolated grid cells (aveget>0)

             ! last veget for this point
             sum_veg=SUM(veget_last(ib,:))
             !
             ! If the DGVM is activated, only anthropic PFTs are utpdated, the others are copied from previous time-step 
             veget_next(ib,:) = veget_last(ib,:)
             
             DO jv = 2, nvmap
                IF ( .NOT. natural(jv) .OR. pasture(jv)) THEN       
                   veget_next(ib,jv) = vegetrefrac(ib,jv)
                ENDIF
             ENDDO

             sumvAnthro_old = zero
             sumvAnthro     = zero
             DO jv = 2, nvmap
                IF ( .NOT. natural(jv) .OR. pasture(jv)) THEN
                   sumvAnthro = sumvAnthro + veget_next(ib,jv)
                   sumvAnthro_old = sumvAnthro_old + veget_last(ib,jv)
                ENDIF
             ENDDO

             IF ( sumvAnthro_old < sumvAnthro ) THEN
                ! Increase of non natural vegetations (increase of agriculture)
                ! The proportion of natural PFT's must be preserved
                ! ie the sum of vegets is preserved
                !    and natural PFT / (sum of veget - sum of antropic veget)
                !    is preserved. 
                rapport = ( sum_veg - sumvAnthro ) / ( sum_veg - sumvAnthro_old )
                DO jv = 1, nvmap
                   IF ( natural(jv) .AND. .NOT. pasture(jv)) THEN
                      veget_next(ib,jv) = veget_last(ib,jv) * rapport
                   ENDIF
                ENDDO
             ELSE
                ! Increase of natural vegetations (decrease of agriculture)
                ! The decrease of agriculture is replaced by bare soil. The DGVM will
                ! re-introduce natural PFT's.
                DO jv = 1, nvmap
                   IF ( natural(jv) .AND. .NOT. pasture(jv)) THEN
                      veget_next(ib,jv) = veget_last(ib,jv)
                   ENDIF
                ENDDO
                veget_next(ib,1) = veget_next(ib,1) + sumvAnthro_old - sumvAnthro
             ENDIF

             ! test
             IF ( ABS( SUM(veget_next(ib,:)) - sum_veg ) > 10*EPSILON(un) ) THEN
                WRITE(numout,*) 'slowproc_readvegetmax _______'
                msg = "  No conservation of sum of veget for point "
                WRITE(numout,*) TRIM(msg), ib, ",(", lalo(ib,1),",", lalo(ib,2), ")" 
                WRITE(numout,*) "  last sum of veget ", sum_veg, " new sum of veget ",                &
                  SUM(veget_next(ib,:)), " error : ", SUM(veget_next(ib,:))-sum_veg
                WRITE(numout,*) "  Anthropic modifications : last ",sumvAnthro_old," new ",sumvAnthro     
                CALL ipslerr_p (3,'slowproc_readvegetmax',                                            &
                     &          'No conservation of sum of veget_next',                               &
                     &          "The sum of veget_next is different after reading Land Use map.",     &
                     &          '(verify the dgvm case model.)')
             ENDIF
          ELSE
             ! For the case when there was a propblem with the interpolation, aveget < min_sechiba
             WRITE(numout,*) 'slowproc_readvegetmax _______'
             WRITE(numout,*) "  No land point in the map for point ", ib, ",(", lalo(ib,1), ",",      &
               lalo(ib,2),")" 
             CALL ipslerr_p (2,'slowproc_readvegetmax',                                               &
                  &          'Problem with vegetation file for Land Use.',                            &
                  &          "No land point in the map for point",                                    & 
                  &          '(verify your land use file.)')
             veget_next(ib,:) = veget_last(ib,:)
          ENDIF
          
       ENDDO
    ENDIF

    IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_readvegetmax after updating'
    !
    frac_nobio_next (:,:) = un
    !
    ! Work only for one nnobio !! (ie ice)
    DO inobio=1,nnobio
       DO jv=1,nvmap
          DO ib = 1, nbpt
             frac_nobio_next(ib,inobio) = frac_nobio_next(ib,inobio) - veget_next(ib,jv)
          ENDDO
       ENDDO
    ENDDO

    DO ib = 1, nbpt
       sum_veg = SUM(veget_next(ib,:))
       sum_nobio = SUM(frac_nobio_next(ib,:))
       IF (sum_nobio < 0.) THEN
          frac_nobio_next(ib,:) = zero
          veget_next(ib,1) = veget_next(ib,1) + sum_nobio
          sum_veg = SUM(veget_next(ib,:))
       ENDIF
       sumf = sum_veg + sum_nobio
       IF (sumf > min_sechiba) THEN
          veget_next(ib,:) = veget_next(ib,:) / sumf
          frac_nobio_next(ib,:) = frac_nobio_next(ib,:) / sumf
          norm=SUM(veget_next(ib,:))+SUM(frac_nobio_next(ib,:))
          err=norm-un
          IF (printlev_loc >=5) WRITE(numout,*) "  slowproc_readvegetmax: ib ",ib,                    &
            " SUM(veget_next(ib,:)+frac_nobio_next(ib,:))-un, sumf",err,sumf
          IF (abs(err) > -EPSILON(un)) THEN
             IF ( SUM(frac_nobio_next(ib,:)) > min_sechiba ) THEN
                frac_nobio_next(ib,1) = frac_nobio_next(ib,1) - err
             ELSE
                veget_next(ib,1) = veget_next(ib,1) - err
             ENDIF
             norm=SUM(veget_next(ib,:))+SUM(frac_nobio_next(ib,:))
             err=norm-un
             IF (printlev_loc >=5) WRITE(numout,*) "  slowproc_readvegetmax: ib ", ib,                &
               " SUM(veget_next(ib,:)+frac_nobio_next(ib,:))-un",err
             IF (abs(err) > EPSILON(un)) THEN
                WRITE(numout,*) '  slowproc_readvegetmax _______'
                WRITE(numout,*) "update : Problem with point ",ib,",(",lalo(ib,1),",",lalo(ib,2),")" 
                WRITE(numout,*) "         err(sum-1.) = ",abs(err)
                CALL ipslerr_p (2,'slowproc_readvegetmax', &
                     &          'Problem with sum vegetation + sum fracnobio for Land Use.',          &
                     &          "sum not equal to 1.", &
                     &          '(verify your land use file.)')
                aveget(ib) = -0.6
             ENDIF
          ENDIF
       ELSE
          ! sumf < min_sechiba
          WRITE(numout,*) '  slowproc_readvegetmax _______'
          WRITE(numout,*)"    No vegetation nor frac_nobio for point ", ib, ",(", lalo(ib,1), ",",    &
            lalo(ib,2),")" 
          WRITE(numout,*)"    Replaced by bare_soil !! "
          veget_next(ib,1) = un
          veget_next(ib,2:nvmap) = zero
          frac_nobio_next(ib,:) = zero
       ENDIF
    ENDDO

    ! Write diagnostics
    CALL xios_orchidee_send_field("aveget",aveget)

    IF (printlev_loc >= 3) WRITE(numout,*) '  slowproc_readvegetmax ended'
    
  END SUBROUTINE slowproc_readvegetmax

!! ================================================================================================================================
!! SUBROUTINE   : slowproc_interpol
!!
!>\BRIEF         Interpolate the IGBP vegetation map to the grid of the model
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_interpol(nbpt, lalo, neighbours, resolution, contfrac, veget, frac_nobio)

    USE interpweight

    IMPLICIT NONE

    !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)          :: nbpt                  !! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)             :: lalo(nbpt,2)          !! Vector of latitude and longitudes (beware of the order!)
    INTEGER(i_std), INTENT(in)          :: neighbours(nbpt,NbNeighb)!! Vector of neighbours for each grid point
                                                                    !! (1=North and then clockwise)
    REAL(r_std), INTENT(in)              :: resolution(nbpt,2)   !! The size in km of each grid-box in X and Y
    REAL(r_std),DIMENSION (nbpt), INTENT (in) :: contfrac        !! Fraction of continent in the grid
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), INTENT(out)    ::  veget(nbpt,nvm)              !! Vegetation fractions
    REAL(r_std), INTENT(out)    ::  frac_nobio(nbpt,nnobio)      !! Fraction of the mesh which is covered by ice, lakes, ...
    !
    !  0.3 LOCAL
    !
    INTEGER(i_std), PARAMETER  :: nolson = 94                    !! Number of Olson classes
    REAL(r_std)                :: resollon, resollat             !! resolution of the longitudes and the latitudes 
                                                                 !!   in the input data which it is in a Goode compressed projection
    !
    !
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: iml, jml, lml, tml, fid, ib, vid
    REAL(r_std), DIMENSION(1)                            :: lev
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: lat_ful, lon_ful, vegmap
    REAL(r_std) :: vegcorr(nolson,nvm)
    REAL(r_std) :: nobiocorr(nolson,nnobio)
    REAL(r_std) :: sumf
    INTEGER(i_std) :: jv, inear
    INTEGER                  :: ALLOC_ERR
    INTEGER                                              :: Ndimslonlat      !! Number of dimensions of lon/lat
    CHARACTER(LEN=1)                                     :: dimlLS           
    INTEGER                                              :: dim1Dlonlat      !! Length of 1D longitudes, latitudes
    INTEGER, DIMENSION(2)                                :: invardims2D
    REAL(r_std), DIMENSION(nbpt,nolson)                  :: vegetrefrac      !! vegegt fractions re-dimensioned
    REAL(r_std), DIMENSION(nbpt)                         :: aveget5k         !! Availability of the interpolation
    REAL(r_std), ALLOCATABLE, DIMENSION(:)               :: aveget5k_glob    !! Availability of the interpolation
    REAL(r_std)                                          :: vmin, vmax       !! min/max values to use for the 
                                                                             !!   renormalization
    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
    CHARACTER(LEN=80)                                    :: lonname, latname !! lonm lat names in input file
    REAL(r_std), DIMENSION(nolson)                       :: variabletypevals !! Values for all the types of the variable
                                                                             !!   (variabletypevals(1) = -un, not used)
    CHARACTER(LEN=50)                                    :: fractype         !! method of calculation of fraction
                                                                             !!   'XYKindTime': Input values are kinds 
                                                                             !!     of something with a temporal 
                                                                             !!     evolution on the dx*dy matrix'
    LOGICAL                                              :: nonegative       !! whether negative values should be removed
    CHARACTER(LEN=50)                                    :: maskingtype      !! Type of masking
                                                                             !!   'nomask': no-mask is applied
                                                                             !!   'mbelow': take values below maskvals(1)
                                                                             !!   'mabove': take values above maskvals(1)
                                                                             !!   'msumrange': take values within 2 ranges;
                                                                             !!      maskvals(2) <= SUM(vals(k)) <= maskvals(1)
                                                                             !!      maskvals(1) < SUM(vals(k)) <= maskvals(3)
                                                                             !!        (normalized by maskvals(3))
                                                                             !!   'var': mask values are taken from a 
                                                                             !!     variable inside the file (>0)
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to 
                                                                             !!   `maskingtype') 
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask 
    LOGICAL                                              :: foundnegvals     !! whether negative aveget5k values 
                                                                             !!   where found
    CHARACTER(LEN=250)                                   :: msg
    INTEGER                                              :: rcode

!_ ================================================================================================================================
    variablename = 'vegetation_map'

    CALL get_vegcorr (nolson,vegcorr,nobiocorr)
    !Config Key   = VEGETATION_FILE
    !Config Desc  = Name of file from which the vegetation map is to be read
    !Config If    = NOT(IMPOSE_VEG) and NOT(MAP_PFT_FORMAT)
    !Config Def   = carteveg5km.nc
    !Config Help  = The name of the file to be opened to read the vegetation
    !Config         map is to be given here. Usualy SECHIBA runs with a 5kmx5km
    !Config         map which is derived from the IGBP one. We assume that we have
    !Config         a classification in 87 types. This is Olson modified by Viovy.
    !Config Units = [FILE]
    !
    filename = 'carteveg5km.nc'
    CALL getin_p('VEGETATION_FILE',filename)
    
    IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_interpol: Read and interpolate " &
         // TRIM(filename) // " for variable " // TRIM(variablename)

! Assigning values to vmin, vmax
    vmin = un
    vmax = nolson*un
    !
    !
    ALLOC_ERR=-1
    variabletypevals = -un

    !! Variables for interpweight
    ! Type of calculation of cell fractions
    fractype = 'default'
    ! Name of the longitude and latitude in the input file
    lonname = 'longitude'
    latname = 'latitude'
    ! Should negative values be set to zero from input file?
    nonegative = .FALSE.
    ! Type of mask to apply to the input data (see header for more details)
    maskingtype = 'mabove'
    ! Values to use for the masking
    maskvals = (/ min_sechiba, undef_sechiba, undef_sechiba /)
    ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
    namemaskvar = ''
    ! Meridional and zonal resolutions of the input data [m]
    resollon = 5000.*un
    resollat = 5000.*un

    CALL interpweight_1D(nbpt, nolson, variabletypevals, lalo, resolution, neighbours,        &
      contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
      maskvals, namemaskvar, 0, 0, -1, fractype,                                                      &
      resollon, resollat, vegetrefrac, aveget5k)
    IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_interpol after interpweight_1D'


    !
    ! Some assumptions on the vegetation file. This information should be
    ! be computed or read from the file. 
    ! It is the resolution in meters of the grid of the vegetation file.
    !
    !
    ! Now we know how many points of which Olson type from the fine grid fall
    ! into each box of the (coarse) model grid: n_origveg(nbpt,nolson)
    !
    ! vegetrefrac is already normalized in subroutine interpweight_1D
    !
    ! now finally calculate coarse vegetation map
    ! Find which model vegetation corresponds to each Olson type 
    !
    veget(:,:) = zero
    frac_nobio(:,:) = zero
    
    DO vid = 1, nolson
       DO jv = 1, nvm
          veget(:,jv) = veget(:,jv) + vegetrefrac(:,vid) * vegcorr(vid,jv)
       ENDDO
    
       DO jv = 1, nnobio
          frac_nobio(:,jv) = frac_nobio(:,jv) + vegetrefrac(:,vid) * nobiocorr(vid,jv)
       ENDDO
    ENDDO
 
    IF (printlev_loc >= 5) THEN
      WRITE(numout,*)'  slowproc_interpol before starting loop nbpt:', nbpt
    END IF 

    ! Getting input longitude and latitude matrices for looking nearest cell
    ! Looking on the global grid if there are points without interpolated values
    IF (is_root_prc) THEN
      ALLOCATE(aveget5k_glob(iim_g*jjm_g))
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_interpol','Problem in allocation of variable aveget5k_glo','','')
    ELSE
      ALLOCATE (aveget5k_glob(1))
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_interpol','Problem in allocation of variable aveget5k_glo','','')
    ENDIF
    CALL gather(aveget5k,aveget5k_glob)

    IF (is_root_prc) THEN
      foundnegvals = ANY(aveget5k_glob .lt. zero)
    END IF
    CAll bcast(foundnegvals)

    IF (foundnegvals) THEN
      ! lon, lat matrices of the input data have to be recupered...
      WRITE(numout,*) '  Looking for nearest point on the 5 km map'
      IF (is_root_prc) THEN 
        CALL flininfo(filename, dim1Dlonlat, jml, lml, tml, fid) 
        !Ndimslonlat = interpweight_get_varNdims_file(filename, TRIM(lonname))
      END IF
      !CALL bcast(Ndimslonlat)
      CALL bcast(dim1Dlonlat)
      Ndimslonlat = 1
      IF (Ndimslonlat ==1) THEN
        ALLOCATE(lon_ful(dim1Dlonlat), STAT=ALLOC_ERR)
        IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_interpol','Problem in allocation of variable lon_ful','','')
        ALLOCATE(lat_ful(dim1Dlonlat), STAT=ALLOC_ERR)
        IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_interpol','Problem in allocation of variable lat_ful','','')
        ALLOCATE(vegmap(dim1Dlonlat), STAT=ALLOC_ERR)
        IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_interpol','Problem in allocation of variable vegmap','','')
        IF (is_root_prc) THEN 
          CALL flinget(fid, TRIM(lonname), dim1Dlonlat, 0, 0, 0, 1, 1, lon_ful)
          CALL flinget(fid, TRIM(latname), dim1Dlonlat, 0, 0, 0, 1, 1, lat_ful)
          CALL flinget(fid, TRIM(variablename), dim1Dlonlat, 0, 0, 0, 1, 1, vegmap)
          CALL flinclo(fid)
        END IF
      ELSE
        WRITE(dimlLS,'(A1)')Ndimslonlat
        msg = "Problem in rank of '" // TRIM(lonname) // "': " // dimlLS // " Not ready !!"
        CALL ipslerr_p(3,'slowproc_interpol',TRIM(msg),'','')
      END IF
    END IF

    CALL bcast(lon_ful)
    CALL bcast(lat_ful)
    CALL bcast(vegmap)

    DEALLOCATE(aveget5k_glob)
   
    !
    !   Clean up the point of the map
    !
    DO ib = 1, nbpt
       !
       !  Let us see if all points found something in the 5km map !
       !
       IF ( aveget5k(ib) .EQ. -1 ) THEN
          !
          ! Now we need to handle some exceptions
          !
          IF ( lalo(ib,1) .LT. -56.0) THEN
             ! Antartica
             frac_nobio(ib,:) = zero
             frac_nobio(ib,iice) = un
             veget(ib,:) = zero
!             aveget5k(ib) = -1.2
          ELSE IF ( lalo(ib,1) .GT. 70.0) THEN
             ! Artica
             frac_nobio(ib,:) = zero
             frac_nobio(ib,iice) = un
             veget(ib,:) = zero
!             aveget5k(ib) = -1.2
          ELSE IF ( lalo(ib,1) .GT. 55.0 .AND. lalo(ib,2) .GT. -65.0 .AND. lalo(ib,2) .LT. -20.0) THEN
             ! Greenland
             frac_nobio(ib,:) = zero
             frac_nobio(ib,iice) = un
             veget(ib,:) = zero
!             aveget5k(ib) = -1.2
          ELSE
             WRITE(numout,*) '  slowproc_interpol _______'
             WRITE(numout,*) '  PROBLEM, no point in the 5km map found for this grid box',ib
             WRITE(numout,*) '  Longitude range : ', lalo(ib,2)
             WRITE(numout,*) '  Latitude range : ', lalo(ib,1)
             
             CALL slowproc_nearest (dim1Dlonlat, lon_ful, lat_ful, &
                  lalo(ib,2), lalo(ib,1), inear)
             WRITE(numout,*) '  Coordinates of the nearest point:', &
                  lon_ful(inear),lat_ful(inear)
             
             DO jv = 1, nvm
                veget(ib,jv) = vegcorr(NINT(vegmap(inear)),jv)
             ENDDO
             
             DO jv = 1, nnobio
                frac_nobio(ib,jv) = nobiocorr(NINT(vegmap(inear)),jv)
             ENDDO
          ENDIF
       ENDIF
       !
       !
       !  Limit the smallest vegetation fraction to 0.5%
       !
       DO vid = 1, nvm
          IF ( veget(ib,vid) .LT. min_vegfrac ) THEN
             veget(ib,vid) = zero
          ENDIF
       ENDDO
  
       sumf = SUM(frac_nobio(ib,:))+SUM(veget(ib,:))
       frac_nobio(ib,:) = frac_nobio(ib,:)/sumf
       veget(ib,:) = veget(ib,:)/sumf
    ENDDO
    
    IF (ALLOCATED(vegmap)) DEALLOCATE(vegmap)
    IF (ALLOCATED(lon_ful)) DEALLOCATE(lon_ful)
    IF (ALLOCATED(lat_ful)) DEALLOCATE(lat_ful)

    ! Write diagnostics
    CALL xios_orchidee_send_field("aveget5k",aveget5k)
    
    IF (printlev_loc >= 3) WRITE(numout,*) '  slowproc_interpol ended'

  END SUBROUTINE slowproc_interpol

!! ================================================================================================================================
!! SUBROUTINE   : slowproc_interpol_g
!!
!>\BRIEF         Interpolate the IGBP vegetation map to the grid of the model
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::veget, ::frac_nobio
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_interpol_g(nbpt, lalo, neighbours, resolution, contfrac, veget, frac_nobio )
    !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)           :: nbpt                  !! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)              :: lalo(nbpt,2)          !! Vector of latitude and longitudes 
                                                                  !! (beware of the order : 1=latitude ; 2=longitude)
    INTEGER(i_std), INTENT(in)           :: neighbours(nbpt,8)    !! Vector of neighbours for each grid point 
                                                                  !! (1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW)
    REAL(r_std), INTENT(in)              :: resolution(nbpt,2)    !! The size in km of each grid-box in X and Y
    REAL(r_std),DIMENSION (nbpt), INTENT (in) :: contfrac         !! Fraction of continent in the grid
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), INTENT(out)    ::  veget(nbpt,nvm)               !! Vegetation fractions
    REAL(r_std), INTENT(out)    ::  frac_nobio(nbpt,nnobio)       !! Fraction of the mesh which is covered by ice, lakes, ...
    !
    LOGICAL ::           ok_interpol                              !! optionnal return of aggregate_vec
    !
    !  0.3 LOCAL
    !
    INTEGER(i_std), PARAMETER                       :: nolson = 94      !! Number of Olson classes
    !
    !
    CHARACTER(LEN=80) :: filename                                       !!vegetation map filename
    INTEGER(i_std) :: iml, jml, lml, tml, fid, ib, ip, vid              
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: lat_ful, lon_ful, vegmap  !! for 5km vegetation map 
                                                                        !! latitude vector, longitude vector, and 
                                                                        !! value of Olson's classes for each location
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: sub_area                !! the area of the fine grid in the model grid ??? 
                                                                        !! cf src_global/interpol_help.f90, line 377, called "areaoverlap"
    INTEGER(i_std),ALLOCATABLE, DIMENSION(:,:) :: sub_index             !! the indexes from the grid boxes from the data that go 
                                                                        !! into the model's boxes 
                                                                        !! cf src_global/interpol_help.f90,line 300, called "ip"
    REAL(r_std), DIMENSION(nbpt,nolson) :: n_origveg                    !! number of points of each Olson type from the fine grid 
                                                                        !! in each box of the (coarse) model grid 
    REAL(r_std), DIMENSION(nbpt) :: n_found                             !! total number of different Olson types found in each 
                                                                        !! box of the (coarse) model grid
    REAL(r_std), DIMENSION(nbpt,nolson) :: frac_origveg                 !! fraction of each Olson type in each box of the (coarse) model grid
    REAL(r_std) :: vegcorr(nolson,nvm)                                  !! correspondance table between Olson and the following SECHIBA Classes.
                                                                        !!   vegcorr(i,:)+nobiocorr(i,:) = 1.  for all i 
                                                                        !! see each class in src_parameters/constantes_veg.f90

    REAL(r_std) :: nobiocorr(nolson,nnobio)                             !! non-biospheric surface typesi
    CHARACTER(LEN=40) :: callsign                                       !! Allows to specify which variable is beeing treated
    REAL(r_std) :: sumf, resol_lon, resol_lat                           !! sumf = sum veget + sum nobio
                                                                        !! resol_lon, resol_lat  reolution in meters of the grid of the vegetation file
    INTEGER(i_std) :: idi, jv, inear, nbvmax                            !! idi : counter for nbvmax, see below   
                                                                        !! jv : counter for nvm, number of PFT
                                                                        !! inear : location of the point of vegmap, which is the closest from the modelled point
                                                                        !! nbvmax : number of maximum vegetation map points in the GCM grid 
    INTEGER(i_std) :: nix, njx
    !
    INTEGER                  :: ALLOC_ERR                               !! location of the eventual missing value in vegmap

!_ ================================================================================================================================
    !
    n_origveg(:,:) = zero
    n_found(:) = zero
    !
    CALL get_vegcorr (nolson,vegcorr,nobiocorr)
    !
    !Config Key   = VEGETATION_FILE
    !Config Desc  = Name of file from which the vegetation map is to be read
    !Config If    = NOT(IMPOSE_VEG) and NOT(MAP_PFT_FORMAT)
    !Config Def   = carteveg5km.nc
    !Config Help  = The name of the file to be opened to read the vegetation
    !Config         map is to be given here. Usualy SECHIBA runs with a 5kmx5km
    !Config         map which is derived from the IGBP one. We assume that we have
    !Config         a classification in 87 types. This is Olson modified by Viovy.
    !Config Units = [FILE]
    !
    filename = 'carteveg5km.nc'
    CALL getin_p('VEGETATION_FILE',filename)  ! GETIN_P !!
    !
    CALL flininfo(filename, iml, jml, lml, tml, fid)   
    !
    ! see IOIPSL/src/flincom.f90, line 665
    ! fid      : File ID
    !- iml      | These 4 variables give the size of the variables
    !- jml      | to be read. It will be verified that the variables
    !- lml      | fits in there.
    !- tml     
    ! iml, jml : horizontal size of the grid, lml = vertical size
    ! tml : size of time axis

    ! TL : pourquoi 2 variables pour la taille horizontale ? cf
    ! IOIPSL/src/flincom.f90 , line 160 

    ALLOC_ERR=-1
    ALLOCATE(lat_ful(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lat_ful : ",ALLOC_ERR
      STOP 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(lon_ful(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lon_ful : ",ALLOC_ERR
      STOP 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(vegmap(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of vegmap : ",ALLOC_ERR
      STOP 
    ENDIF
    !
    WRITE(numout,*) 'Reading the OLSON type vegetation file'
    !
    CALL flinget(fid, 'longitude', iml, jml, lml, tml, 1, 1, lon_ful)
    CALL flinget(fid, 'latitude', iml, jml, lml, tml, 1, 1, lat_ful)
    CALL flinget(fid, 'vegetation_map', iml, jml, lml, tml, 1, 1, vegmap)
    !
    WRITE(numout,*) 'File name : ', filename
    WRITE(numout,*) 'Min and max vegetation numbers : ', MINVAL(vegmap), MAXVAL(vegmap)
    !
    CALL flinclo(fid)
    !
    IF (MAXVAL(vegmap) .LT. nolson) THEN
       WRITE(numout,*) 'WARNING -- WARNING'
       WRITE(numout,*) 'The vegetation map has too few vegetation types.'
       WRITE(numout,*) 'If you are lucky it will work but please check'
    ELSE IF ( MAXVAL(vegmap) .GT. nolson) THEN
       WRITE(numout,*) 'More vegetation types in file than the code can'
       WRITE(numout,*) 'deal with.: ',  MAXVAL(vegmap),  nolson
       STOP 'slowproc_interpol'
    ENDIF
    !
    ! Some assumptions on the vegetation file. This information should be
    ! be computed or read from the file. 
    ! It is the reolution in meters of the grid of the vegetation file.
    !
    
    !TL : CODE EN DUR ????? 
    resol_lon = 5000.
    resol_lat = 5000.
    !
    !
    ! The number of maximum vegetation map points in the GCM grid is estimated.
    ! Some margin is taken.
    !
    nix=INT(MAXVAL(resolution_g(:,1)*2)/resol_lon)+1
    njx=INT(MAXVAL(resolution_g(:,2)*2)/resol_lon)+1
    nbvmax = nix*njx
    !
    ! No need to broadcast as this routine is only called on root_proc
    !
    callsign="Vegetation map"
    !
    ok_interpol = .FALSE.
    DO WHILE ( .NOT. ok_interpol )
       WRITE(numout,*) "Projection arrays for ",callsign," : "
       WRITE(numout,*) "nbvmax = ",nbvmax
       !
       ALLOC_ERR=-1
       ALLOCATE(sub_index(nbpt, nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) THEN
          WRITE(numout,*) "ERROR IN ALLOCATION of sub_index : ",ALLOC_ERR
          STOP 
       ENDIF
       sub_index(:,:)=0
       ALLOC_ERR=-1
       ALLOCATE(sub_area(nbpt, nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) THEN
          WRITE(numout,*) "ERROR IN ALLOCATION of sub_area : ",ALLOC_ERR
          STOP 
       ENDIF
       sub_area(:,:)=zero
       !
       CALL aggregate_p (nbpt, lalo, neighbours, resolution, contfrac, &
            &                iml, lon_ful, lat_ful, resol_lon, resol_lat, callsign, &
            &                nbvmax, sub_index, sub_area, ok_interpol)
       !
       ! Defined as aggregate_2d or aggregate_vec in src_global/interpol_help.f90, depending
       ! on the dimensions (2D region or vector)i. 
       ! This routing will get for each point of the coarse grid the
       ! indexes of the finer grid and the area of overlap.
       ! This routine is designed for a fine grid which is regular in lat/lon.

       IF ( .NOT. ok_interpol ) THEN
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
          !
          nbvmax = nbvmax * 2
       ELSE
          !
          DO ib = 1, nbpt
             DO idi=1, nbvmax
                ! Leave the do loop if all sub areas are treated, sub_area <= 0
                IF ( sub_area(ib,idi) <= zero ) EXIT

                ip = sub_index(ib,idi)
                n_origveg(ib,NINT(vegmap(ip))) = n_origveg(ib,NINT(vegmap(ip))) + sub_area(ib,idi)
                n_found(ib) =  n_found(ib) + sub_area(ib,idi)
             ENDDO
          ENDDO
          !
       ENDIF
    ENDDO
    !
    ! Now we know how many points of which Olson type from the fine grid fall
    ! into each box of the (coarse) model grid: n_origveg(nbpt,nolson)
    !
    !
    ! determine fraction of Olson vegetation type in each box of the coarse grid
    !
    DO vid = 1, nolson
       WHERE ( n_found(:) .GT. 0 ) 
          frac_origveg(:,vid) = n_origveg(:,vid) / n_found(:)
       ELSEWHERE
          frac_origveg(:,vid) = zero
       ENDWHERE
    ENDDO
    !
    ! now finally calculate coarse vegetation map
    ! Find which model vegetation corresponds to each Olson type 
    !
    veget(:,:) = zero
    frac_nobio(:,:) = zero
    !
    DO vid = 1, nolson
       !
       DO jv = 1, nvm
          veget(:,jv) = veget(:,jv) + frac_origveg(:,vid) * vegcorr(vid,jv)
       ENDDO
       !
       DO jv = 1, nnobio
          frac_nobio(:,jv) = frac_nobio(:,jv) + frac_origveg(:,vid) * nobiocorr(vid,jv)
       ENDDO
       !
    ENDDO
    !
    WRITE(numout,*) 'slowproc_interpol : Interpolation Done'
    !
    !   Clean up the point of the map
    !
    DO ib = 1, nbpt
       !
       !  Let us see if all points found something in the 5km map !
       !
       IF ( n_found(ib) .EQ. 0 ) THEN
          !
          ! Now we need to handle some exceptions
          !
          IF ( lalo(ib,1) .LT. -56.0) THEN
             ! Antartica
             frac_nobio(ib,:) = zero
             frac_nobio(ib,iice) = un
             veget(ib,:) = zero
             !
          ELSE IF ( lalo(ib,1) .GT. 70.0) THEN
             ! Artica
             frac_nobio(ib,:) = zero
             frac_nobio(ib,iice) = un
             veget(ib,:) = zero
             !
          ELSE IF ( lalo(ib,1) .GT. 55.0 .AND. lalo(ib,2) .GT. -65.0 .AND. lalo(ib,2) .LT. -20.0) THEN
             ! Greenland
             frac_nobio(ib,:) = zero
             frac_nobio(ib,iice) = un
             veget(ib,:) = zero
             !
          ELSE
             !
             WRITE(numout,*) 'PROBLEM, no point in the 5km map found for this grid box',ib
             WRITE(numout,*) 'Longitude range : ', lalo(ib,2)
             WRITE(numout,*) 'Latitude range : ', lalo(ib,1)
             !
             WRITE(numout,*) 'Looking for nearest point on the 5 km map'
             CALL slowproc_nearest (iml, lon_ful, lat_ful, &
                  lalo(ib,2), lalo(ib,1), inear)
             WRITE(numout,*) 'Coordinates of the nearest point:', &
                  lon_ful(inear),lat_ful(inear)
             !
             DO jv = 1, nvm
                veget(ib,jv) = vegcorr(NINT(vegmap(inear)),jv)
             ENDDO
             !
             DO jv = 1, nnobio
                frac_nobio(ib,jv) = nobiocorr(NINT(vegmap(inear)),jv)
             ENDDO
             !
          ENDIF
          !
       ENDIF
       !
       !
       !  Limit the smallest vegetation fraction to 0.5%
       !
       DO vid = 1, nvm
          IF ( veget(ib,vid) .LT. min_vegfrac ) THEN  ! min_vegfrac=0.001 in constantes_veg.f90
             veget(ib,vid) = zero
          ENDIF
       ENDDO
       !
       sumf = SUM(frac_nobio(ib,:))+SUM(veget(ib,:))
       frac_nobio(ib,:) = frac_nobio(ib,:)/sumf
       veget(ib,:) = veget(ib,:)/sumf
       !
       !       
    ENDDO
    !
    DEALLOCATE(vegmap)
    DEALLOCATE(lat_ful, lon_ful)
    DEALLOCATE(sub_index)
    DEALLOCATE(sub_area)

    !
    RETURN
    !
  END SUBROUTINE slowproc_interpol_g


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_nearest
!!
!>\BRIEF         looks for nearest grid point on the fine map
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::inear
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_nearest(iml, lon5, lat5, lonmod, latmod, inear)

    !! INTERFACE DESCRIPTION
    
    !! 0.1 input variables

    INTEGER(i_std), INTENT(in)                   :: iml             !! size of the vector
    REAL(r_std), DIMENSION(iml), INTENT(in)      :: lon5, lat5      !! longitude and latitude vector, for the 5km vegmap
    REAL(r_std), INTENT(in)                      :: lonmod, latmod  !! longitude  and latitude modelled

    !! 0.2 output variables
    
    INTEGER(i_std), INTENT(out)                  :: inear           !! location of the grid point from the 5km vegmap grid
                                                                    !! closest from the modelled grid point

    !! 0.4 Local variables

    REAL(r_std)                                  :: pa, p
    REAL(r_std)                                  :: coscolat, sincolat
    REAL(r_std)                                  :: cospa, sinpa
    REAL(r_std), ALLOCATABLE, DIMENSION(:)       :: cosang
    INTEGER(i_std)                               :: i
    INTEGER(i_std), DIMENSION(1)                 :: ineartab
    INTEGER                                      :: ALLOC_ERR

!_ ================================================================================================================================

    ALLOCATE(cosang(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) CALL ipslerr_p(3,'slowproc_nearest','Error in allocation for cosang','','')

    pa = pi/2.0 - latmod*pi/180.0 ! dist. between north pole and the point a 
                                                      !! COLATITUDE, in radian
    cospa = COS(pa)
    sinpa = SIN(pa)

    DO i = 1, iml

       sincolat = SIN( pi/2.0 - lat5(i)*pi/180.0 ) !! sinus of the colatitude
       coscolat = COS( pi/2.0 - lat5(i)*pi/180.0 ) !! cosinus of the colatitude

       p = (lonmod-lon5(i))*pi/180.0 !! angle between a & b (between their meridian)in radians

       !! dist(i) = ACOS( cospa*coscolat + sinpa*sincolat*COS(p))
       cosang(i) = cospa*coscolat + sinpa*sincolat*COS(p) !! TL : cosang is maximum when angle is at minimal value  
!! orthodromic distance between 2 points : cosang = cosinus (arc(AB)/R), with
!R = Earth radius, then max(cosang) = max(cos(arc(AB)/R)), reached when arc(AB)/R is minimal, when
! arc(AB) is minimal, thus when point B (corresponding grid point from LAI MAP) is the nearest from
! modelled A point
    ENDDO

    ineartab = MAXLOC( cosang(:) )
    inear = ineartab(1)

    DEALLOCATE(cosang)
  END SUBROUTINE slowproc_nearest

!! ================================================================================================================================
!! SUBROUTINE   : slowproc_soilt
!!
!>\BRIEF         Interpolate the Zobler soil type map
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::soiltype, ::clayfraction
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE slowproc_soilt(nbpt, lalo, neighbours, resolution, contfrac, soilclass, clayfraction)

    USE interpweight

    IMPLICIT NONE
    !
    !
    !   This subroutine should read the Zobler map and interpolate to the model grid. The method
    !   is to get fraction of the three main soiltypes for each grid box.
    !   The soil fraction are going to be put into the array soiltype in the following order :
    !   coarse, medium and fine.
    !
    !
    !!  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)    :: nbpt                   !! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)       :: lalo(nbpt,2)           !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), INTENT(in)    :: neighbours(nbpt,NbNeighb)!! Vector of neighbours for each grid point
                                                              !! (1=North and then clockwise)
    REAL(r_std), INTENT(in)       :: resolution(nbpt,2)     !! The size in km of each grid-box in X and Y
    REAL(r_std), INTENT(in)       :: contfrac(nbpt)         !! Fraction of land in each grid box.
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), INTENT(out)      :: soilclass(nbpt, nscm)  !! Soil type map to be created from the Zobler map
    REAL(r_std), INTENT(out)      :: clayfraction(nbpt)     !! The fraction of clay as used by STOMATE
    !
    !
    !  0.3 LOCAL
    !
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: ib, ilf, nbexp, i
    INTEGER(i_std) :: fopt                                  !! Nb of pts from the texture map within one ORCHIDEE grid-cell
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: soiltext, soiltext2
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)  :: sub_area
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:)  :: sub_index
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:) :: solt       !! Texture the different points from the input texture map 
                                                            !! in one ORCHIDEE grid cell (unitless)
    !
    ! Number of texture classes in Zobler
    !
    INTEGER(i_std), PARAMETER :: nzobler = 7
    REAL(r_std),ALLOCATABLE   :: textfrac_table(:,:)
    !   
    INTEGER                  :: ALLOC_ERR
    INTEGER                                              :: ntextinfile      !! number of soil textures in the in the file
    REAL(r_std), DIMENSION(:,:), ALLOCATABLE             :: textrefrac       !! text fractions re-dimensioned
    REAL(r_std), DIMENSION(nbpt)                         :: atext            !! Availability of the text interpolation
    REAL(r_std)                                          :: vmin, vmax       !! min/max values to use for the 

    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat name in input file
    REAL(r_std), DIMENSION(:), ALLOCATABLE               :: variabletypevals !! Values for all the types of the variable
                                                                             !!   (variabletypevals(1) = -un, not used)
    CHARACTER(LEN=50)                                    :: fractype         !! method of calculation of fraction
                                                                             !!   'XYKindTime': Input values are kinds 
                                                                             !!     of something with a temporal 
                                                                             !!     evolution on the dx*dy matrix'
    LOGICAL                                              :: nonegative       !! whether negative values should be removed
    CHARACTER(LEN=50)                                    :: maskingtype      !! Type of masking
                                                                             !!   'nomask': no-mask is applied
                                                                             !!   'mbelow': take values below maskvals(1)
                                                                             !!   'mabove': take values above maskvals(1)
                                                                             !!   'msumrange': take values within 2 ranges;
                                                                             !!      maskvals(2) <= SUM(vals(k)) <= maskvals(1)
                                                                             !!      maskvals(1) < SUM(vals(k)) <= maskvals(3)
                                                                             !!        (normalized by maskvals(3))
                                                                             !!   'var': mask values are taken from a 
                                                                             !!     variable inside the file (>0)
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to 
                                                                             !!   `maskingtype') 
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask 
    INTEGER(i_std), DIMENSION(:), ALLOCATABLE            :: vecpos
    REAL(r_std)                                          :: sgn              !! sum of fractions excluding glaciers and ocean
!_ ================================================================================================================================

    IF (printlev_loc>=3) WRITE (numout,*) 'slowproc_soilt'
    !
    !  Needs to be a configurable variable
    !
    !
    !Config Key   = SOILCLASS_FILE
    !Config Desc  = Name of file from which soil types are read
    !Config Def   = soils_param.nc
    !Config If    = NOT(IMPOSE_VEG)
    !Config Help  = The name of the file to be opened to read the soil types. 
    !Config         The data from this file is then interpolated to the grid of
    !Config         of the model. The aim is to get fractions for sand loam and
    !Config         clay in each grid box. This information is used for soil hydrology
    !Config         and respiration.
    !Config Units = [FILE]
    !
    filename = 'soils_param.nc'
    CALL getin_p('SOILCLASS_FILE',filename)

    variablename = 'soiltext'

    !! Variables for interpweight
    ! Type of calculation of cell fractions
    fractype = 'default'
    ! Name of the longitude and latitude in the input file
    lonname = 'nav_lon'
    latname = 'nav_lat'

    IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_soilt: Read and interpolate " &
         // TRIM(filename) // " for variable " // TRIM(variablename)

    IF ( TRIM(soil_classif) /= 'none' ) THEN

       ! Define a variable for the number of soil textures in the input file
       SELECTCASE(soil_classif)
       CASE('zobler')
          ntextinfile=nzobler
       CASE('usda')
          ntextinfile=nscm
       CASE DEFAULT
          WRITE(numout,*) 'slowproc_soilt:'
          WRITE(numout,*) '  A non supported soil type classification has been chosen'
          CALL ipslerr_p(3,'slowproc_soilt','non supported soil type classification','','')
       ENDSELECT

       ALLOCATE(textrefrac(nbpt,ntextinfile), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_soilt','Problem in allocation of variable textrefrac',&
         '','')

       ! Assigning values to vmin, vmax
       vmin = un
       vmax = ntextinfile*un
       
       ALLOCATE(variabletypevals(ntextinfile), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_soilt','Problem in allocation of variabletypevals','','')
       variabletypevals = -un
       
       !! Variables for interpweight
       ! Should negative values be set to zero from input file?
       nonegative = .FALSE.
       ! Type of mask to apply to the input data (see header for more details)
       maskingtype = 'mabove'
       ! Values to use for the masking
       maskvals = (/ min_sechiba, undef_sechiba, undef_sechiba /)
       ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') ( not used)
       namemaskvar = ''
       
       CALL interpweight_2D(nbpt, ntextinfile, variabletypevals, lalo, resolution, neighbours,        &
          contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,    & 
          maskvals, namemaskvar, 0, 0, -1, fractype, -1., -1., textrefrac, atext)

       ALLOCATE(vecpos(ntextinfile), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_soilt','Problem in allocation of variable vecpos','','')
       ALLOCATE(solt(ntextinfile), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_soilt','Problem in allocation of variable solt','','')
       
       IF (printlev_loc >= 5) THEN
          WRITE(numout,*)'  slowproc_soilt after interpweight_2D'
          WRITE(numout,*)'  slowproc_soilt before starting loop nbpt:', nbpt
          WRITE(numout,*)"  slowproc_soilt starting classification '" // TRIM(soil_classif) // "'..."
       END IF
    ELSE
      IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_soilt using default values all points are propertly ' // &
        'interpolated atext = 1. everywhere!'
      atext = 1.
    END IF

    nbexp = 0
    SELECTCASE(soil_classif)
    CASE('none')
       ALLOCATE(textfrac_table(nscm,ntext), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) CALL ipslerr_p(3,'slowproc_soilt','Error in allocation for textfrac_table','','')
       DO ib=1, nbpt
          soilclass(ib,:) = soilclass_default_fao
          clayfraction(ib) = clayfraction_default
       ENDDO
    CASE('zobler')
       !
       soilclass_default=soilclass_default_fao ! FAO means here 3 final texture classes
       !
       IF (printlev_loc>=2) WRITE(numout,*) "Using a soilclass map with Zobler classification"
       !
       ALLOCATE(textfrac_table(nzobler,ntext), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) CALL ipslerr_p(3,'slowproc_soilt','Error in allocation for textfrac_table','','')
       CALL get_soilcorr_zobler (nzobler, textfrac_table)
       !
       !
       IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_soilt after getting table of textures'
       DO ib =1, nbpt
          soilclass(ib,:) = zero
          clayfraction(ib) = zero
          !
          ! vecpos: List of positions where textures were not zero
          ! vecpos(1): number of not null textures found
          vecpos = interpweight_ValVecR(textrefrac(ib,:),nzobler,zero,'neq')
          fopt = vecpos(1)

          IF ( fopt .EQ. 0 ) THEN
             ! No points were found for current grid box, use default values
             nbexp = nbexp + 1
             soilclass(ib,:) = soilclass_default(:)
             clayfraction(ib) = clayfraction_default
          ELSE
             IF (fopt == nzobler) THEN
                ! All textures are not zero
                solt=(/(i,i=1,nzobler)/)
             ELSE
               DO ilf = 1,fopt
                 solt(ilf) = vecpos(ilf+1)
               END DO
             END IF
             !
             !   Compute the fraction of each textural class
             !
             sgn = 0.
             DO ilf = 1,fopt
                !
                ! Here we make the correspondance between the 7 zobler textures and the 3 textures in ORCHIDEE
                ! and soilclass correspond to surfaces covered by the 3 textures of ORCHIDEE (coase,medium,fine)
                ! For type 6 = glacier, default values are set and it is also taken into account during the normalization 
                ! of the fractions (done in interpweight_2D)
                ! Note that type 0 corresponds to ocean but it is already removed using the mask above.
                !
                IF ( (solt(ilf) .LE. nzobler) .AND. (solt(ilf) .GT. 0) .AND. & 
                     (solt(ilf) .NE. 6) ) THEN
                   SELECT CASE(solt(ilf))
                     CASE(1)
                        soilclass(ib,1) = soilclass(ib,1) + textrefrac(ib,solt(ilf))
                     CASE(2)
                        soilclass(ib,2) = soilclass(ib,2) + textrefrac(ib,solt(ilf))
                     CASE(3)
                        soilclass(ib,2) = soilclass(ib,2) + textrefrac(ib,solt(ilf))
                     CASE(4)
                        soilclass(ib,2) = soilclass(ib,2) + textrefrac(ib,solt(ilf))
                     CASE(5)
                        soilclass(ib,3) = soilclass(ib,3) + textrefrac(ib,solt(ilf))
                     CASE(7)
                        soilclass(ib,2) = soilclass(ib,2) + textrefrac(ib,solt(ilf))
                     CASE DEFAULT
                        WRITE(numout,*) 'We should not be here, an impossible case appeared'
                        CALL ipslerr_p(3,'slowproc_soilt','Bad value for solt','','')
                   END SELECT
                   ! clayfraction is the sum of the % of clay (as a mineral of small granulometry, and not as a texture)
                   ! over the zobler pixels composing the ORCHIDEE grid-cell
                   clayfraction(ib) = clayfraction(ib) + &
                        & textfrac_table(solt(ilf),3) * textrefrac(ib,solt(ilf))
                   ! Sum the fractions which are not glaciers nor ocean
                   sgn = sgn + textrefrac(ib,solt(ilf))
                ELSE
                   IF (solt(ilf) .GT. nzobler) THEN
                      WRITE(numout,*) 'The file contains a soil color class which is incompatible with this program'
                      CALL ipslerr_p(3,'slowproc_soilt','Problem soil color class incompatible','','')
                   ENDIF
                ENDIF
             ENDDO

             IF ( sgn .LT. min_sechiba) THEN
                ! Set default values if grid cells were only covered by glaciers or ocean 
                ! or if now information on the source grid was found.
                nbexp = nbexp + 1
                soilclass(ib,:) = soilclass_default(:)
                clayfraction(ib) = clayfraction_default
             ELSE
                ! Normalize using the fraction of surface not including glaciers and ocean
                soilclass(ib,:) = soilclass(ib,:)/sgn
                clayfraction(ib) = clayfraction(ib)/sgn
             ENDIF
          ENDIF
       ENDDO
   
    !
    CASE("fao")
       !
       soilclass_default=soilclass_default_fao
       !
       WRITE(numout,*) "Using a soilclass map with fao classification"
       !
       ALLOCATE(textfrac_table(nscm,ntext))
       !
       CALL get_soilcorr_zobler (nscm, textfrac_table)
       !
       DO ib =1, nbpt
          !
          ! GO through the point we have found
          !
          !
          fopt = COUNT(sub_area(ib,:) > zero)
          !
          !    Check that we found some points
          !
          soilclass(ib,:) = 0.0
          clayfraction(ib) = 0.0
          !
          IF ( fopt .EQ. 0) THEN
             nbexp = nbexp + 1
             soilclass(ib,:) = soilclass_default(:)
             clayfraction(ib) = clayfraction_default
          ELSE
             !
             DO ilf = 1,fopt
                solt(ilf) = soiltext(sub_index(ib,ilf,1),sub_index(ib,ilf,2))
             ENDDO
             !
             !
             !   Compute the average bare soil albedo parameters
             !
             sgn = zero
             !
             DO ilf = 1,fopt
                !
                ! 
                !
                IF ( (solt(ilf) .LE. nscm) .AND. (solt(ilf) .GT. 0) ) THEN
                   soilclass(ib,solt(ilf)) = soilclass(ib,solt(ilf)) + sub_area(ib,ilf)
                   clayfraction(ib) = clayfraction(ib) + textfrac_table(solt(ilf),3) * sub_area(ib,ilf)
                   sgn = sgn + sub_area(ib,ilf)
                ELSE
                   IF (solt(ilf) .GT. nscm) THEN
                      WRITE(*,*) 'The file contains a soil color class which is incompatible with this program'
                      STOP 'slowproc_soilt'
                   ENDIF
                ENDIF
                !
             ENDDO
             !
             ! Normalize the surface
             !
             IF ( sgn .LT. min_sechiba) THEN
                nbexp = nbexp + 1
                soilclass(ib,:) = soilclass_default(:)
                clayfraction(ib) = clayfraction_default
             ELSE
                soilclass(ib,:) = soilclass(ib,:)/sgn
                clayfraction(ib) = clayfraction(ib)/sgn
             ENDIF
             !
          ENDIF
          !
       ENDDO
   
    ! The "USDA" case reads a map of the 12 USDA texture classes, 
    ! such as to assign the corresponding soil properties
    CASE("usda")
       IF (printlev_loc>=2) WRITE(numout,*) "Using a soilclass map with usda classification"

       soilclass_default=soilclass_default_usda

       ALLOCATE(textfrac_table(nscm,ntext), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) CALL ipslerr_p(3,'slowproc_soilt','Error in allocation for textfrac_table','','')

       CALL get_soilcorr_usda (nscm, textfrac_table)

       IF (printlev_loc>=4) WRITE (numout,*) 'slowproc_soilt: After get_soilcorr_usda'
       !
       DO ib =1, nbpt
          !
          ! GO through the point we have found
          !
          !
          ! Provide which textures were found
          ! vecpos: List of positions where textures were not zero
          !   vecpos(1): number of not null textures found
          vecpos = interpweight_ValVecR(textrefrac(ib,:),ntextinfile,zero,'neq')
          fopt = vecpos(1)
         
          !
          !    Check that we found some points
          !
          soilclass(ib,:) = 0.0
          clayfraction(ib) = 0.0
          
          IF ( fopt .EQ. 0) THEN
             ! No points were found for current grid box, use default values
             IF (printlev_loc>=3) WRITE(numout,*)'slowproc_soilt: no soil class in input file found for point=', ib
             nbexp = nbexp + 1
             soilclass(ib,:) = soilclass_default
             clayfraction(ib) = clayfraction_default
          ELSE
             IF (fopt == nscm) THEN
                ! All textures are not zero
                solt(:) = (/(i,i=1,nscm)/)
             ELSE
               DO ilf = 1,fopt
                 solt(ilf) = vecpos(ilf+1) 
               END DO
             END IF
             !
             !
             !   Compute the fraction of each textural class  
             !
             !
             DO ilf = 1,fopt
                IF ( (solt(ilf) .LE. nscm) .AND. (solt(ilf) .GT. 0) ) THEN
                   soilclass(ib,solt(ilf)) = textrefrac(ib,solt(ilf))
                   clayfraction(ib) = clayfraction(ib) + textfrac_table(solt(ilf),3) *                &
                        textrefrac(ib,solt(ilf))
                ELSE
                   IF (solt(ilf) .GT. nscm) THEN
                      WRITE(*,*) 'The file contains a soil color class which is incompatible with this program'
                      CALL ipslerr_p(3,'slowproc_soilt','Problem soil color class incompatible 2','','')
                   ENDIF
                ENDIF
                !
             ENDDO

             ! Set default values if the surface in source file is too small
             IF ( atext(ib) .LT. min_sechiba) THEN
                nbexp = nbexp + 1
                soilclass(ib,:) = soilclass_default(:)
                clayfraction(ib) = clayfraction_default
             ENDIF
          ENDIF

       ENDDO
       
       IF (printlev_loc>=4) WRITE (numout,*) '  slowproc_soilt: End case usda'
       
    CASE DEFAULT
       WRITE(numout,*) 'slowproc_soilt _______'
       WRITE(numout,*) '  A non supported soil type classification has been chosen'
       CALL ipslerr_p(3,'slowproc_soilt','non supported soil type classification','','')
    ENDSELECT
    IF (printlev_loc >= 5 ) WRITE(numout,*)'  slowproc_soilt end of type classification'

    IF ( nbexp .GT. 0 ) THEN
       WRITE(numout,*) 'slowproc_soilt:'
       WRITE(numout,*) '  The interpolation of the bare soil albedo had ', nbexp
       WRITE(numout,*) '  points without data. This are either coastal points or ice covered land.'
       WRITE(numout,*) '  The problem was solved by using the default soil types.'
    ENDIF

    IF (ALLOCATED(variabletypevals)) DEALLOCATE (variabletypevals)
    IF (ALLOCATED(textrefrac)) DEALLOCATE (textrefrac)
    IF (ALLOCATED(solt)) DEALLOCATE (solt)
    IF (ALLOCATED(textfrac_table)) DEALLOCATE (textfrac_table)

    ! Write diagnostics
    CALL xios_orchidee_send_field("atext",atext)
    
    IF (printlev_loc >= 3) WRITE(numout,*) '  slowproc_soilt ended'

  END SUBROUTINE slowproc_soilt
 
!! ================================================================================================================================
!! SUBROUTINE   : slowproc_slope
!!
!>\BRIEF         Calculate mean slope coef in each  model grid box from the slope map
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::reinf_slope
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_slope(nbpt, lalo, neighbours, resolution, contfrac, reinf_slope)

    USE interpweight

    IMPLICIT NONE

    !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)          :: nbpt                  ! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)              :: lalo(nbpt,2)          ! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), INTENT(in)          :: neighbours(nbpt,NbNeighb)! Vector of neighbours for each grid point
                                                                    ! (1=North and then clockwise)
    REAL(r_std), INTENT(in)              :: resolution(nbpt,2)    ! The size in km of each grid-box in X and Y
    REAL(r_std), INTENT (in)             :: contfrac(nbpt)         !! Fraction of continent in the grid
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), INTENT(out)    ::  reinf_slope(nbpt)                   ! slope coef 
    !
    !  0.3 LOCAL
    !
    !
    REAL(r_std)  :: slope_noreinf                 ! Slope above which runoff is maximum
    CHARACTER(LEN=80) :: filename
    REAL(r_std)                                          :: vmin, vmax       !! min/max values to use for the 
                                                                             !!   renormalization
    REAL(r_std), DIMENSION(nbpt)                         :: aslope           !! slope availability 

    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat name in the input file
    CHARACTER(LEN=50)                                    :: fractype         !! method of calculation of fraction
                                                                             !!   'XYKindTime': Input values are kinds 
                                                                             !!     of something with a temporal 
                                                                             !!     evolution on the dx*dy matrix'
    LOGICAL                                              :: nonegative       !! whether negative values should be removed
    CHARACTER(LEN=50)                                    :: maskingtype      !! Type of masking
                                                                             !!   'nomask': no-mask is applied
                                                                             !!   'mbelow': take values below maskvals(1)
                                                                             !!   'mabove': take values above maskvals(1)
                                                                             !!   'msumrange': take values within 2 ranges;
                                                                             !!      maskvals(2) <= SUM(vals(k)) <= maskvals(1)
                                                                             !!      maskvals(1) < SUM(vals(k)) <= maskvals(3)
                                                                             !!        (normalized by maskvals(3))
                                                                             !!   'var': mask values are taken from a 
                                                                             !!     variable inside the file  (>0)
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to 
                                                                             !!   `maskingtype') 
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask 

!_ ================================================================================================================================
    
    !
    !Config Key   = SLOPE_NOREINF
    !Config Desc  = See slope_noreinf above
    !Config If    = 
    !Config Def   = 0.5
    !Config Help  = The slope above which there is no reinfiltration
    !Config Units = [-]
    !
    slope_noreinf = 0.5
    !
    CALL getin_p('SLOPE_NOREINF',slope_noreinf)
    !
    !Config Key   = TOPOGRAPHY_FILE
    !Config Desc  = Name of file from which the topography map is to be read
    !Config If    = 
    !Config Def   = cartepente2d_15min.nc
    !Config Help  = The name of the file to be opened to read the orography
    !Config         map is to be given here. Usualy SECHIBA runs with a 2'
    !Config         map which is derived from the NGDC one. 
    !Config Units = [FILE]
    !
    filename = 'cartepente2d_15min.nc'
    CALL getin_p('TOPOGRAPHY_FILE',filename)

    variablename = 'pente'
    IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_slope: Read and interpolate " &
         // TRIM(filename) // " for variable " // TRIM(variablename)

    ! For this case there are not types/categories. We have 'only' a continuos field
    ! Assigning values to vmin, vmax
    vmin = 0.
    vmax = 9999.

    !! Variables for interpweight
    ! Type of calculation of cell fractions
    fractype = 'slopecalc'
    ! Name of the longitude and latitude in the input file
    lonname = 'longitude'
    latname = 'latitude'
    ! Should negative values be set to zero from input file?
    nonegative = .FALSE.
    ! Type of mask to apply to the input data (see header for more details)
    maskingtype = 'mabove'
    ! Values to use for the masking
    maskvals = (/ min_sechiba, undef_sechiba, undef_sechiba /)
    ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
    namemaskvar = ''

    CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,                                &
      contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
      maskvals, namemaskvar, -1, fractype, slope_default, slope_noreinf,                              &
      reinf_slope, aslope)
    IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_slope after interpweight_2Dcont'

    ! Write diagnostics
    CALL xios_orchidee_send_field("aslope",aslope)

    IF (printlev_loc >= 3) WRITE(numout,*) '  slowproc_slope ended'

  END SUBROUTINE slowproc_slope

!! ================================================================================================================================
!! SUBROUTINE 	: get_vegcorr
!!
!>\BRIEF         The "get_vegcorr" routine defines the table of correspondence
!!               between the 94 Olson vegetation types and the 13 Plant Functional Types known 
!!               by SECHIBA and STOMATE. Used by slowproc for the old interpolation.
!!
!!\DESCRIPTION : get_vegcorr is needed if you use the old_map carteveg5km.nc. \n
!!               Usually SECHIBA can run with a 5kmx5km map which is derived from the IGBP one. \n
!!               We assume that we have a classification in 94 types. This is Olson one modified by Nicolas Viovy.\n
!!               ORCHIDEE has to convert the Olson vegetation types into PFTs for the run (interpolation step).\n
!!               Each Olson matches to a combination of fractions of one or several PFTs.\n
!!               This routine uses the same process for the non-biospheric map (not used).\n
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::vegcorr, ::nobiocorr.
!!
!! REFERENCE(S)	: 
!! - Olson, J.S., J.A. Watts, and L.J. Allison., 1983. 
!! "Carbon in Live Vegetation of Major World Ecosystems." 
!! Report ORNL-5862. Oak Ridge National Laboratory, Oak Ridge, Tennessee. 
!! - Olson, J.S., J.A. Watts, and L.J. Allison., 1985. 
!! "Major World Ecosystem Complexes Ranked by Carbon in Live Vegetation: A Database." 
!! NDP-017. Carbon Dioxide Information Center, Oak Ridge National Laboratory, Oak Ridge, Tennessee.
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================

  SUBROUTINE get_vegcorr (nolson,vegcorr,nobiocorr)

    IMPLICIT NONE

    !! 0. Variables and parameters declaration
    
    INTEGER(i_std),PARAMETER :: nolson94 = 94                       !! Number of Olson vegetation types (unitless)
    INTEGER(i_std),PARAMETER :: nvm13 = 13                          !! Number of PFTS of ORCHIDEE (unitless) 

    !! 0.1 Input variables

    INTEGER(i_std),INTENT(in) :: nolson                             !! Number of Olson vegetation types (unitless)
    
    !! 0.2 Output variables 

    REAL(r_std),DIMENSION(nolson,nvm),INTENT(out) :: vegcorr        !! Correspondence array between Olson types and PFTS 
                                                                    !! (0-1, unitless)
    REAL(r_std),DIMENSION(nolson,nnobio),INTENT(out) :: nobiocorr   !! Correspondence array between non-vegetation types and nobio 
                                                                    !! types (lake,etc..) (0-1, unitless)

    !! 0.4 Local variable
    
    INTEGER(i_std) :: ib                                            !! Indice (unitless)
    
 !_ ================================================================================================================================

    !-
    ! 0. Check consistency
    !-
    IF (nolson /= nolson94) THEN
       WRITE(numout,*) nolson,nolson94
       CALL ipslerr_p(3,'get_vegcorr', '', '',&
            &                 'wrong number of OLSON vegetation types.') ! Fatal error
    ENDIF !(nolson /= nolson94)
    
    IF (nvm /= nvm13) THEN
       WRITE(numout,*) nvm,nvm13
       CALL ipslerr_p(3,'get_vegcorr', '', '',&
            &                 'wrong number of SECHIBA vegetation types.') ! Fatal error
    ENDIF !(nvm /= nvm13)

    ! The carteveg5km cannot be used if the PFTs are not in the standard order
    DO ib = 1,nvm
       IF (pft_to_mtc(ib) /= ib ) THEN
          CALL ipslerr_p(3,'get_vegcorr','You have redefined the order of the 13 PFTS', & 
               &          'You can not use carteveg5km', 'Use the standard configuration of PFTS' )
       ENDIF
    ENDDO

    !-
    ! 1 set the indices of non-biospheric surface types to 0.
    !-
    nobiocorr(:,:) = zero
    !-
    ! 2 Here we construct the correspondance table
    !   between Olson and the following SECHIBA Classes.
    !   vegcorr(i,:)+nobiocorr(i,:) = 1.  for all i.
    !-
    ! The modified OLSON types found in file carteveg5km.nc
    ! created by Nicolas Viovy :
    !  1 Urban
    vegcorr( 1,:) = &
         & (/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    !  2 Cool low sparse grassland
    vegcorr( 2,:) = &
         & (/0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0/)
    !  3 Cold conifer forest
    vegcorr( 3,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    !  4 Cold deciduous conifer forest
    vegcorr( 4,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0/)
    !  5 Cool Deciduous broadleaf forest
    vegcorr( 5,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    !  6 Cool evergreen broadleaf forests
    vegcorr( 6,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    !  7 Cool tall grasses and shrubs
    vegcorr( 7,:) = &
         & (/0.1, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0/)
    !  8 Warm C3 tall grasses and shrubs
    vegcorr( 8,:) = &
         & (/0.1, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0/)
    !  9 Warm C4 tall grases and shrubs
    vegcorr( 9,:) = &
         & (/0.1, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0/)
    ! 10 Bare desert
    vegcorr(10,:) = &
         & (/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    ! 11 Cold upland tundra
    vegcorr(11,:) = &
         & (/0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0/)
    ! 12 Cool irrigated grassland
    vegcorr(12,:) = &
         & (/0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 0.0/)
    ! 13 Semi desert
    vegcorr(13,:) = &
         & (/0.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0/)
    ! 14 Glacier ice
    vegcorr(14,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    nobiocorr(14,iice) = 1.
    ! 15 Warm wooded wet swamp
    vegcorr(15,:) = &
         & (/0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0/)
    ! 16 Inland water
    vegcorr(16,:) = &
         & (/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    ! 17 sea water
    vegcorr(17,:) = &
         & (/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    ! 18 cool shrub evergreen
    vegcorr(18,:) = &
         & (/0.1, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0/)
    ! 19 cold shrub deciduous
    vegcorr(19,:) = &
         & (/0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.6, 0.0, 0.0, 0.0/)
    ! 20 Cold evergreen forest and fields
    vegcorr(20,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0/)
    ! 21 cool rain forest
    vegcorr(21,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0/)
    ! 22 cold conifer boreal forest
    vegcorr(22,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0/)
    ! 23 cool conifer forest
    vegcorr(23,:) = &
         & (/0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0/)
    ! 24 warm mixed forest
    vegcorr(24,:) = &
         & (/0.0, 0.4, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0/)
    ! 25 cool mixed forest
    vegcorr(25,:) = &
         & (/0.0, 0.0, 0.0, 0.4, 0.0, 0.4, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0/)
    ! 26 cool broadleaf forest
    vegcorr(26,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0/)
    ! 27 cool deciduous broadleaf forest
    vegcorr(27,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.3, 0.5, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0/)
    ! 28 warm montane tropical forest
    vegcorr(28,:) = &
         & (/0.0, 0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0/)
    ! 29 warm seasonal tropical forest
    vegcorr(29,:) = &
         & (/0.0, 0.5, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0/)
    ! 30 cool crops and towns
    vegcorr(30,:) = &
         & (/0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0/)
    ! 31 warm crops and towns
    vegcorr(31,:) = &
         & (/0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8/)
    ! 32 cool crops and towns
    vegcorr(32,:) = &
         & (/0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0/)
    ! 33 warm dry tropical woods
    vegcorr(33,:) = &
         & (/0.2, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0/)
    ! 34 warm tropical rain forest
    vegcorr(34,:) = &
         & (/0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    ! 35 warm tropical degraded forest
    vegcorr(35,:) = &
         & (/0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0/)
    ! 36 warm corn and beans cropland
    vegcorr(36,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0/)
    ! 37 cool corn and bean cropland
    vegcorr(37,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0/)
    ! 38 warm rice paddy and field
    vegcorr(38,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0/)
    ! 39 hot irrigated cropland
    vegcorr(39,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0/)
    ! 40 cool irrigated cropland
    vegcorr(40,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0/)
    ! 41 cold irrigated cropland
    vegcorr(41,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0/)
    ! 42 cool grasses and shrubs
    vegcorr(42,:) = &
         & (/0.1, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.7, 0.0, 0.0, 0.0/)
    ! 43 hot and mild grasses and shrubs
    vegcorr(43,:) = &
         & (/0.2, 0.0, 0.1, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0/)
    ! 44 cold grassland
    vegcorr(44,:) = &
         & (/0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 0.0/)
    ! 45 Savanna (woods) C3
    vegcorr(45,:) = &
         & (/0.1, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7, 0.0, 0.0, 0.0/)
    ! 46 Savanna woods C4
    vegcorr(46,:) = &
         & (/0.1, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7, 0.0, 0.0/)
    ! 47 Mire, bog, fen
    vegcorr(47,:) = &
         & (/0.1, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7, 0.0, 0.0, 0.0/)
    ! 48 Warm marsh wetland
    vegcorr(48,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0/)
    ! 49 cold marsh wetland
    vegcorr(49,:) = &
         & (/0.0, 0.0, 0.0, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0/)
    ! 50 mediteraean scrub
    vegcorr(50,:) = &
         & (/0.1, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0/)
    ! 51 Cool dry woody scrub
    vegcorr(51,:) = &
         & (/0.3, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0/)
    ! 52 Warm dry evergreen woods
    vegcorr(52,:) = &
         & (/0.1, 0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    ! 53 Volcanic rocks
    vegcorr(53,:) = &
         & (/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    ! 54 sand desert
    vegcorr(54,:) = &
         & (/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    ! 55 warm semi desert shrubs
    vegcorr(55,:) = &
         & (/0.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0/)
    ! 56 cool semi desert shrubs
    vegcorr(56,:) = &
         & (/0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.0, 0.0, 0.0/)
    ! 57 semi desert sage
    vegcorr(57,:) = &
         & (/0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0/)
    ! 58 Barren tundra
    vegcorr(58,:) = &
         & (/0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.0, 0.0, 0.0/)
    ! 59 cool southern hemisphere mixed forest
    vegcorr(59,:) = &
         & (/0.1, 0.0, 0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0/)
    ! 60 cool fields and woods
    vegcorr(60,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0/)
    ! 61 warm forest and filed
    vegcorr(61,:) = &
         & (/0.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6/)
    ! 62 cool forest and field
    vegcorr(62,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0/)
    ! 63 warm C3 fields and woody savanna
    vegcorr(63,:) = &
         & (/0.1, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0/)
    ! 64 warm C4 fields and woody savanna
    vegcorr(64,:) = &
         & (/0.1, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6/)
    ! 65 cool fields and woody savanna
    vegcorr(65,:) = &
         & (/0.0, 0.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0/)
    ! 66 warm succulent and thorn scrub
    vegcorr(66,:) = &
         & (/0.1, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0/)
    ! 67 cold small leaf mixed woods
    vegcorr(67,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.3, 0.0, 0.5, 0.0, 0.0, 0.0/)
    ! 68 cold deciduous and mixed boreal fores
    vegcorr(68,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0/)
    ! 69 cold narrow conifers
    vegcorr(69,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0/)
    ! 70 cold wooded tundra
    vegcorr(70,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.7, 0.0, 0.0, 0.0/)
    ! 71 cold heath scrub
    vegcorr(71,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.7, 0.0, 0.0, 0.0/)
    ! 72 Polar and alpine desert
    vegcorr(72,:) = &
         & (/0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0/)
    ! 73 warm Mangrove
    vegcorr(73,:) = &
         & (/0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    ! 74 cool crop and water mixtures
    vegcorr(74,:) = &
         & (/0.1, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0/)
    ! 75 cool southern hemisphere mixed forest
    vegcorr(75,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.4, 0.4, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0/)
    ! 76 cool moist eucalyptus
    vegcorr(76,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0/)
    ! 77 warm rain green tropical forest
    vegcorr(77,:) = &
         & (/0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    ! 78 warm C3 woody savanna
    vegcorr(78,:) = &
         & (/0.0, 0.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0/)
    ! 79 warm C4 woody savanna
    vegcorr(79,:) = &
         & (/0.0, 0.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0/)
    ! 80 cool woody savanna
    vegcorr(80,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.0, 0.6, 0.0, 0.0, 0.0/)
    ! 81 cold woody savanna
    vegcorr(81,:) = &
         & (/0.0, 0.0, 0.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0/)
    ! 82 warm broadleaf crops
    vegcorr(82,:) = &
         & (/0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9, 0.0/)
    ! 83 warm C3 grass crops
    vegcorr(83,:) = &
         & (/0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9, 0.0/)
    ! 84 warm C4 grass crops
    vegcorr(84,:) = &
         & (/0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9/)
    ! 85 cool grass crops
    vegcorr(85,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0/)
    ! 86 warm C3 crops grass,shrubs
    vegcorr(86,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0/)
    ! 87 cool crops,grass,shrubs
    vegcorr(87,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.0/)
    ! 88 warm evergreen tree crop
    vegcorr(88,:) = &
         & (/0.0, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2/)
    ! 89 cool evergreen tree crop
    vegcorr(89,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0/)
    ! 90 cold evergreen tree crop
    vegcorr(90,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0/)
    ! 91 warm deciduous tree crop
    vegcorr(91,:) = &
         & (/0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2/)
    ! 92 cool deciduous tree crop
    vegcorr(92,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0/)
    ! 93 cold deciduous tree crop
    vegcorr(93,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.2, 0.0/)
    ! 94 wet sclerophylic forest
    vegcorr(94,:) = &
         & (/0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    !-
    ! 3 Check the mapping for the Olson types which are going into the
    !   the veget and nobio array.
    !-
    DO ib=1,nolson
       !
       IF ( ABS(SUM(vegcorr(ib,:))+SUM(nobiocorr(ib,:))-1.0) &
            &       > EPSILON(1.0)) THEN
          WRITE(numout,*) 'Wrong correspondance for Olson type :', ib
          CALL ipslerr_p(3,'get_vegcorr', '', '',&
               &                 'Wrong correspondance for Olson type.') ! Fatal error
       ENDIF
       !
    ENDDO ! Loop over the # Olson type


  END SUBROUTINE get_vegcorr

!! ================================================================================================================================
!! SUBROUTINE 	: get_soilcorr_zobler
!!
!>\BRIEF         The "get_soilcorr_zobler" routine defines the table of correspondence
!!               between the Zobler types and the three texture types known by SECHIBA and STOMATE :
!!               silt, sand and clay. 
!!
!! DESCRIPTION : get_soilcorr_zobler is needed if you use soils_param.nc .\n
!!               The data from this file is then interpolated to the grid of the model. \n
!!               The aim is to get fractions for sand loam and clay in each grid box.\n
!!               This information is used for soil hydrology and respiration.
!!
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S) : ::texfrac_table
!!
!! REFERENCE(S)	: 
!! - Zobler L., 1986, A World Soil File for global climate modelling. NASA Technical memorandum 87802. NASA 
!!   Goddard Institute for Space Studies, New York, U.S.A.
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================

  SUBROUTINE get_soilcorr_zobler (nzobler,textfrac_table)

    IMPLICIT NONE

    !! 0. Variables and parameters declaration
    
    INTEGER(i_std),PARAMETER :: nbtypes_zobler = 7                    !! Number of Zobler types (unitless)

    !! 0.1  Input variables
    
    INTEGER(i_std),INTENT(in) :: nzobler                              !! Size of the array (unitless)
    
    !! 0.2 Output variables 
    
    REAL(r_std),DIMENSION(nzobler,ntext),INTENT(out) :: textfrac_table !! Table of correspondence between soil texture class
                                                                       !! and granulometric composition (0-1, unitless)
    
    !! 0.4 Local variables
    
    INTEGER(i_std) :: ib                                              !! Indice (unitless)
    
!_ ================================================================================================================================

    !-
    ! 0. Check consistency
    !-  
    IF (nzobler /= nbtypes_zobler) THEN 
       CALL ipslerr_p(3,'get_soilcorr_zobler', 'nzobler /= nbtypes_zobler',&
          &   'We do not have the correct number of classes', &
          &                 ' in the code for the file.')  ! Fatal error
    ENDIF

    !-
    ! 1. Textural fraction for : silt        sand         clay
    !-
    textfrac_table(1,:) = (/ 0.12, 0.82, 0.06 /)
    textfrac_table(2,:) = (/ 0.32, 0.58, 0.10 /)
    textfrac_table(3,:) = (/ 0.39, 0.43, 0.18 /)
    textfrac_table(4,:) = (/ 0.15, 0.58, 0.27 /)
    textfrac_table(5,:) = (/ 0.34, 0.32, 0.34 /)
    textfrac_table(6,:) = (/ 0.00, 1.00, 0.00 /)
    textfrac_table(7,:) = (/ 0.39, 0.43, 0.18 /)


    !-
    ! 2. Check the mapping for the Zobler types which are going into the ORCHIDEE textures classes 
    !-
    DO ib=1,nzobler ! Loop over # classes soil
       
       IF (ABS(SUM(textfrac_table(ib,:))-1.0) > EPSILON(1.0)) THEN ! The sum of the textural fractions should not exceed 1 !
          WRITE(numout,*) &
               &     'Error in the correspondence table', &
               &     ' sum is not equal to 1 in', ib
          WRITE(numout,*) textfrac_table(ib,:)
          CALL ipslerr_p(3,'get_soilcorr_zobler', 'SUM(textfrac_table(ib,:)) /= 1.0',&
               &                 '', 'Error in the correspondence table') ! Fatal error
       ENDIF
       
    ENDDO ! Loop over # classes soil

    
  END SUBROUTINE get_soilcorr_zobler

!! ================================================================================================================================
!! SUBROUTINE 	: get_soilcorr_usda
!!
!>\BRIEF         The "get_soilcorr_usda" routine defines the table of correspondence
!!               between the 12 USDA textural classes and their granulometric composition, 
!!               as % of silt, sand and clay. This is used to further defien clayfraction.
!!
!! DESCRIPTION : get_soilcorr is needed if you use soils_param.nc .\n
!!               The data from this file is then interpolated to the grid of the model. \n
!!               The aim is to get fractions for sand loam and clay in each grid box.\n
!!               This information is used for soil hydrology and respiration.
!!               The default map in this case is derived from Reynolds et al 2000, \n
!!               at the 1/12deg resolution, with indices that are consistent with the \n
!!               textures tabulated below
!!
!! RECENT CHANGE(S): Created by A. Ducharne on July 02, 2014
!!
!! MAIN OUTPUT VARIABLE(S) : ::texfrac_table
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================

  SUBROUTINE get_soilcorr_usda (nusda,textfrac_table)

    IMPLICIT NONE

    !! 0. Variables and parameters declaration
    
    !! 0.1  Input variables
    
    INTEGER(i_std),INTENT(in) :: nusda                               !! Size of the array (unitless)
    
    !! 0.2 Output variables 
    
    REAL(r_std),DIMENSION(nusda,ntext),INTENT(out) :: textfrac_table !! Table of correspondence between soil texture class
                                                                     !! and granulometric composition (0-1, unitless)
    
    !! 0.4 Local variables

    INTEGER(i_std),PARAMETER :: nbtypes_usda = 12                    !! Number of USDA texture classes (unitless)
    INTEGER(i_std) :: n                                              !! Index (unitless)
    
!_ ================================================================================================================================

    !-
    ! 0. Check consistency
    !-  
    IF (nusda /= nbtypes_usda) THEN 
       CALL ipslerr_p(3,'get_soilcorr_usda', 'nusda /= nbtypes_usda',&
          &   'We do not have the correct number of classes', &
          &                 ' in the code for the file.')  ! Fatal error
    ENDIF

    !! Parameters for soil type distribution :
    !! Sand, Loamy Sand, Sandy Loam, Silt Loam, Silt, Loam, Sandy Clay Loam, Silty Clay Loam, Clay Loam, Sandy Clay, Silty Clay, Clay
    ! The order comes from constantes_soil.f90
    ! The corresponding granulometric composition comes from Carsel & Parrish, 1988

    !-
    ! 1. Textural fractions for : sand, clay
    !-
    textfrac_table(1,2:3)  = (/ 0.93, 0.03 /) ! Sand
    textfrac_table(2,2:3)  = (/ 0.81, 0.06 /) ! Loamy Sand
    textfrac_table(3,2:3)  = (/ 0.63, 0.11 /) ! Sandy Loam
    textfrac_table(4,2:3)  = (/ 0.17, 0.19 /) ! Silt Loam
    textfrac_table(5,2:3)  = (/ 0.06, 0.10 /) ! Silt
    textfrac_table(6,2:3)  = (/ 0.40, 0.20 /) ! Loam
    textfrac_table(7,2:3)  = (/ 0.54, 0.27 /) ! Sandy Clay Loam
    textfrac_table(8,2:3)  = (/ 0.08, 0.33 /) ! Silty Clay Loam
    textfrac_table(9,2:3)  = (/ 0.30, 0.33 /) ! Clay Loam
    textfrac_table(10,2:3) = (/ 0.48, 0.41 /) ! Sandy Clay
    textfrac_table(11,2:3) = (/ 0.06, 0.46 /) ! Silty Clay
    textfrac_table(12,2:3) = (/ 0.15, 0.55 /) ! Clay

    ! Fraction of silt

    DO n=1,nusda
       textfrac_table(n,1) = 1. - textfrac_table(n,2) - textfrac_table(n,3)
    END DO
       
  END SUBROUTINE get_soilcorr_usda

!! ================================================================================================================================
!! FUNCTION 	: tempfunc
!!
!>\BRIEF        ! This function interpolates value between ztempmin and ztempmax
!! used for lai detection. 
!!
!! DESCRIPTION   : This subroutine calculates a scalar between 0 and 1 with the following equation :\n
!!                 \latexonly
!!                 \input{constantes_veg_tempfunc.tex}
!!                 \endlatexonly
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : tempfunc_result
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  FUNCTION tempfunc (temp_in) RESULT (tempfunc_result)


    !! 0. Variables and parameters declaration

    REAL(r_std),PARAMETER    :: ztempmin=273._r_std   !! Temperature for laimin (K)
    REAL(r_std),PARAMETER    :: ztempmax=293._r_std   !! Temperature for laimax (K)
    REAL(r_std)              :: zfacteur              !! Interpolation factor   (K^{-2})

    !! 0.1 Input variables

    REAL(r_std),INTENT(in)   :: temp_in               !! Temperature (K)

    !! 0.2 Result

    REAL(r_std)              :: tempfunc_result       !! (unitless)
    
!_ ================================================================================================================================

    !! 1. Define a coefficient
    zfacteur = un/(ztempmax-ztempmin)**2
    
    !! 2. Computes tempfunc
    IF     (temp_in > ztempmax) THEN
       tempfunc_result = un
    ELSEIF (temp_in < ztempmin) THEN
       tempfunc_result = zero
    ELSE
       tempfunc_result = un-zfacteur*(ztempmax-temp_in)**2
    ENDIF !(temp_in > ztempmax)


  END FUNCTION tempfunc

!!
!================================================================================================================================
!! SUBROUTINE   : slowproc_read_veg_restfile
!!
!>\BRIEF         read a vegetation map from a restart file, and interpolate if
!necessary to current model grid
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: none
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_
!================================================================================================================================

  SUBROUTINE slowproc_read_veg_restfile(nbpt, lalo, neighbours, resolution, contfrac, veget, frac_nobio )
    !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)          :: nbpt                  ! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)              :: lalo(nbpt,2)          ! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), INTENT(in)          :: neighbours(nbpt,8)    ! Vector of neighbours for each grid point 
    ! (1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW)
    REAL(r_std), INTENT(in)              :: resolution(nbpt,2)    ! The size in km of each grid-box in X and Y
    REAL(r_std),DIMENSION (nbpt), INTENT (in) :: contfrac         !! Fraction of continent in the grid
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), INTENT(out)    ::  veget(nbpt,nvm)         ! Vegetation fractions
    REAL(r_std), INTENT(out)    ::  frac_nobio(nbpt,nnobio) ! Fraction of the mesh which is covered by ice, lakes, ...
    !
    !  0.3 LOCAL
    !
    !
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: iml, jml, lml, tml, fid, ib, jp, it, jj, jv, im, jm
    INTEGER(i_std) :: il, ils, ip, ix, iy, imin, jmin, ier
    REAL(r_std) :: dlon, dlonmin, dlat, dlatmin
    REAL(r_std), PARAMETER :: maxmargin = 5.
    !
    !  0.4 allocatable
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:) :: veget_max_in
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) ::  frac_nobio_in
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: xx, yy
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: x, y

    ! get filename
    CALL getin_p('vegetation_map_restartfile_name', filename)

    ! load the sechiba restart file, and read in values of veget_max and
    ! frac_nobio
    CALL flininfo(filename, iml, jml, lml, tml, fid)
    !! for now, assume hard-coded  96 x 72 gridpoints because iml and jml
    !readings are wrong
    im = 96
    jm = 72
    ALLOCATE(xx(im, jm))
    ALLOCATE(yy(im,jm))
    ALLOCATE (x(im),y(jm))
    ALLOCATE(veget_max_in(im,jm,nvm))
    ALLOCATE(frac_nobio_in(im,jm))
    WRITE(*,*) 'cdk debug: '
    WRITE(*,*) 'filename: ', filename
    WRITE(*,*) 'im, jm, lml, tml, fid: ', im, jm, lml, tml, fid

    CALL flinget(fid, 'nav_lon', im, jm, 1, 1, 1, 1, xx)
    WRITE(*,*) 'cdk debug read nav_lon '
    CALL flinget(fid, 'nav_lat', im, jm, 1, 1, 1, 1, yy)
    WRITE(*,*) 'cdk debug read nav_lat '
    CALL flinget(fid, 'veget_max', im, jm, nvm, 1, 1, 1, veget_max_in)  !! these need to be double-precision
    WRITE(*,*) 'cdk debug read veget_max '
    CALL flinget(fid, 'frac_nobio', im, jm, 1, 1, 1, 1, frac_nobio_in)  !! these need to be double-precision
    WRITE(*,*) 'cdk debug read frac_nobio '

    x(:) = xx(:,1)
    y(:) = yy(1,:)

    ! prendre la valeur la plus proche
    DO ip = 1, nbpt
       dlonmin = HUGE(1.)
       DO ix = 1,im
          dlon = MIN( ABS(lalo(ip,2)-x(ix)), ABS(lalo(ip,2)+360.-x(ix)), ABS(lalo(ip,2)-360.-x(ix)) )
          IF ( dlon .LT. dlonmin ) THEN
             imin = ix
             dlonmin = dlon
          ENDIF
       ENDDO
       dlatmin = HUGE(1.)
       DO iy = 1,jm
          dlat = ABS(lalo(ip,1)-y(iy))
          IF ( dlat .LT. dlatmin ) THEN
             jmin = iy
             dlatmin = dlat
          ENDIF
       ENDDO
       WRITE(*,*) 'cdk debug found: ip,  dlonmin, dlatmin: ', ip, dlonmin,dlatmin
       ! imin, jmin, x(imin), y(jmin):,
       ! if nothing is close, then set as all nobio
       IF ( ( dlonmin .LE. maxmargin ) .AND. ( dlatmin .LE. maxmargin ) ) THEN
          veget(ip,:) = veget_max_in(imin,jmin,:)
          frac_nobio(ip,:) = frac_nobio_in(imin,jmin)
       ELSE
          veget(ip,:) = zero
          frac_nobio(ip,:) = un
       ENDIF
       WRITE(*,*) 'cdk debug: ip,  veget, frac_nobio: ', veget(ip,:),frac_nobio(ip,:)

    ENDDO

    DEALLOCATE(xx)
    DEALLOCATE(yy)
    DEALLOCATE(x)
    DEALLOCATE(y)
    DEALLOCATE(veget_max_in)
    DEALLOCATE(frac_nobio_in)
    !
    RETURN
    !
  END SUBROUTINE slowproc_read_veg_restfile


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_checkveget
!!
!>\BRIEF         To verify the consistency of the various fractions defined within the grid box after having been
!!               been updated by STOMATE or the standard procedures.
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: none
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!
  SUBROUTINE slowproc_checkveget(nbpt, frac_nobio, veget_max, veget, tot_bare_soil, soiltile)

    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)                      :: nbpt       ! Number of points for which the data needs to be interpolated
    REAL(r_std),DIMENSION (nbpt,nnobio), INTENT(in) :: frac_nobio ! Fraction of ice,lakes,cities, ... (unitless)
    REAL(r_std),DIMENSION (nbpt,nvm), INTENT(in)    :: veget_max  ! Maximum fraction of vegetation type including none biological fraction (unitless) 
    REAL(r_std),DIMENSION (nbpt,nvm), INTENT(in)    :: veget      ! Vegetation fractions
    REAL(r_std),DIMENSION (nbpt), INTENT(in)        :: tot_bare_soil ! Total evaporating bare soil fraction within the mesh
    REAL(r_std),DIMENSION (nbpt,nstm), INTENT(in)   :: soiltile   ! Fraction of soil tiles in the gridbox (unitless)

    !  0.3 LOCAL
    !
    INTEGER(i_std) :: ji, jn, jv
    REAL(r_std)  :: epsilocal  !! A very small value
    REAL(r_std)  :: totfrac
    CHARACTER(len=80) :: str1, str2
    
!_ ================================================================================================================================
    
    !
    ! There is some margin added as the computing errors might bring us above EPSILON(un)
    !
    epsilocal = EPSILON(un)*1000.
    
    !! 1.0 Verify that none of the fractions are smaller than min_vegfrac, without beeing zero.
    !!
    DO ji=1,nbpt
       DO jn=1,nnobio
          IF ( frac_nobio(ji,jn) > epsilocal .AND. frac_nobio(ji,jn) < min_vegfrac ) THEN
             WRITE(str1,'("Occurs on grid box", I8," and nobio type ",I3 )') ji, jn
             WRITE(str2,'("The small value obtained is ", E14.4)') frac_nobio(ji,jn)
             CALL ipslerr_p (3,'slowproc_checkveget', &
                  "frac_nobio is larger than zero but smaller than min_vegfrac.", str1, str2)
          ENDIF
       ENDDO
    END DO
    
    IF ( (.NOT. ok_dgvm) .AND. (.NOT. ok_dgvm_peat) ) THEN  !!qcj++ peatland      
       DO ji=1,nbpt
          DO jv=1,nvm
             IF ( veget_max(ji,jv) > epsilocal .AND. veget_max(ji,jv) < min_vegfrac ) THEN
                WRITE(str1,'("Occurs on grid box", I8," and nobio type ",I3 )') ji, jv
                WRITE(str2,'("The small value obtained is ", E14.4)') veget_max(ji,jv)
                CALL ipslerr_p (3,'slowproc_checkveget', &
                     "veget_max is larger than zero but smaller than min_vegfrac.", str1, str2)
             ENDIF
          ENDDO
       ENDDO
    END IF
    
    !! 2.0 verify that with all the fractions we cover the entire grid box  
    !!
    DO ji=1,nbpt
       totfrac = zero
       DO jn=1,nnobio
          totfrac = totfrac + frac_nobio(ji,jn)
       ENDDO
       DO jv=1,nvm
          totfrac = totfrac + veget_max(ji,jv)
       ENDDO
       IF ( ABS(totfrac - un) > epsilocal) THEN
             WRITE(str1,'("This occurs on grid box", I8)') ji
             WRITE(str2,'("The sum over all fraction and error are ", E14.4, E14.4)') totfrac, ABS(totfrac - un)
             CALL ipslerr_p (3,'slowproc_checkveget', &
                   "veget_max + frac_nobio is not equal to 1.", str1, str2)
             WRITE(*,*) "EPSILON =", epsilocal 
       ENDIF
    ENDDO
    
    !! 3.0 Verify that veget is smaller or equal to veget_max
    !!
    DO ji=1,nbpt
       DO jv=1,nvm
          IF ( jv == ibare_sechiba ) THEN
             IF ( ABS(veget(ji,jv) - veget_max(ji,jv)) > epsilocal ) THEN
                WRITE(str1,'("This occurs on grid box", I8)') ji
                WRITE(str2,'("The difference is ", E14.4)') veget(ji,jv) - veget_max(ji,jv)
                CALL ipslerr_p (3,'slowproc_checkveget', &
                     "veget is not equal to veget_max on bare soil.", str1, str2)
             ENDIF
          ELSE
             !IF ( veget(ji,jv) > veget_max(ji,jv) ) THEN
             IF ( veget(ji,jv) - veget_max(ji,jv) > epsilocal ) THEN !!!qcj++ peatland
                WRITE(str1,'("This occurs on grid box", I8)') ji
                WRITE(str2,'("The values for veget and veget_max :", F8.4, F8.4)') veget(ji,jv), veget_max(ji,jv)
                CALL ipslerr_p (3,'slowproc_checkveget', &
                     "veget is greater than veget_max.", str1, str2)
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    
    !! 4.0 Test tot_bare_soil in relation to the other variables
    !!
    DO ji=1,nbpt
       totfrac = zero
       DO jv=1,nvm
          totfrac = totfrac + (veget_max(ji,jv) - veget(ji,jv))
       ENDDO
       ! add the bare soil fraction to totfrac
       totfrac = totfrac + veget(ji,ibare_sechiba)
       ! do the test
       IF ( ABS(totfrac - tot_bare_soil(ji)) > epsilocal ) THEN
          WRITE(str1,'("This occurs on grid box", I8)') ji
          WRITE(str2,'("The values for tot_bare_soil, tot frac and error :", F8.4, F8.4, E14.4)') &
               &  tot_bare_soil(ji), totfrac, ABS(totfrac - tot_bare_soil(ji))
          CALL ipslerr_p (3,'slowproc_checkveget', &
               "tot_bare_soil does not correspond to the total bare soil fraction.", str1, str2)
       ENDIF
    ENDDO
    
    !! 5.0 Test that soiltile has the right sum
    !!
    DO ji=1,nbpt
       totfrac = SUM(soiltile(ji,:))
       IF ( ABS(totfrac - un) > epsilocal ) THEN
          WRITE(numout,*) "soiltile does not sum-up to one. This occurs on grid box", ji
          WRITE(numout,*) "The soiltile for ji are :", soiltile(ji,:)
          CALL ipslerr_p (2,'slowproc_checkveget', &
               "soiltile does not sum-up to one.", "", "")
       ENDIF
    ENDDO
    
  END SUBROUTINE slowproc_checkveget
    
!! ================================================================================================================================
!! SUBROUTINE   : slowproc_change_frac
!!
!>\BRIEF        Update the vegetation fractions
!!
!! DESCRIPTION  : Update the vegetation fractions. This subroutine is called in the same time step as lcchange in stomatelpj has
!!                has been done. This subroutine is called after the diagnostics have been written in sechiba_main.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: veget_max, veget, frac_nobio, totfrac_nobio, tot_bare_soil, soiltile
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
   
  SUBROUTINE slowproc_change_frac(kjpindex, f_rot_sech, lai, &
                                  veget_max, veget, frac_nobio, totfrac_nobio, tot_bare_soil, soiltile,veget_cov_max_new_peatdgvm)!, veget_max_adjusted)
    !
    ! 0. Declarations
    !
    ! 0.1 Input variables 
    INTEGER(i_std), INTENT(in)                           :: kjpindex      !! Domain size - terrestrial pixels only
    LOGICAL,DIMENSION(kjpindex),INTENT(in)               :: f_rot_sech    !! whether a grid is under rotation
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(in)     :: lai           !! Leaf area index (m^2 m^{-2})    
!    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(in)     :: veget_max_adjusted !! Maximum fraction of vegetation adjusted with peatland 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(in)     :: veget_cov_max_new_peatdgvm
    ! 0.2 Output variables 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(out)    :: veget_max      !! Maximum fraction of vegetation type in the mesh (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(out)    :: veget          !! Fraction of vegetation type in the mesh (unitless)
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(out) :: frac_nobio     !! Fraction of ice, lakes, cities etc. in the mesh
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)        :: totfrac_nobio  !! Total fraction of ice+lakes+cities etc. in the mesh
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)       :: tot_bare_soil  !! Total evaporating bare soil fraction in the mesh
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT(out)  :: soiltile       !! Fraction of each soil tile within vegtot (0-1, unitless)
    
    ! 0.3 Local variables
    INTEGER(i_std)                                       :: ji, jv         !! Loop index
    
       
    !! Update vegetation fractions with the values coming from the vegetation file read in slowproc_readvegetmax.
    !! Partial update has been taken into account for the case with DGVM and AGRICULTURE in slowproc_readvegetmax.
    

    IF (agri_peat) THEN
       veget_max = veget_cov_max_adjusted
       frac_nobio = frac_nobio_new
    ELSE
       frac_nobio = frac_nobio_new
       IF (ok_dgvm_peat) THEN
         veget_max = veget_cov_max_new_peatdgvm
       ELSE
         veget_max  = veget_max_new
       ENDIF
    ENDIF
 
    !! Verification and correction on veget_max, calculation of veget and soiltile.
    CALL slowproc_veget (kjpindex, f_rot_sech, lai, frac_nobio, totfrac_nobio, veget_max, veget, soiltile)

    !! Calculate tot_bare_soil needed in hydrol, diffuco and condveg (fraction of bare soil in the mesh)
    tot_bare_soil(:) = veget_max(:,1)
    DO jv = 2, nvm
       DO ji =1, kjpindex
          tot_bare_soil(ji) = tot_bare_soil(ji) + (veget_max(ji,jv) - veget(ji,jv))
       ENDDO
    END DO

    !! Do some basic tests on the surface fractions updated above 
    CALL slowproc_checkveget(kjpindex, frac_nobio, veget_max, veget, tot_bare_soil, soiltile)
     
  END SUBROUTINE slowproc_change_frac 

  !spitfire
  SUBROUTINE slowproc_read_data(nbpt, lalo,  resolution, proxydata,data_filename,field_name)

    !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)          :: nbpt                  ! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)             :: lalo(nbpt,2)          ! Vector of latitude and longitudes (beware of the order !)
    REAL(r_std), INTENT(in)             :: resolution(nbpt,2)    ! The size in km of each grid-box in X and Y
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), INTENT(out)    ::  proxydata(:,:)         ! lightn read variable and re-dimensioned
    !
    !  0.3 LOCAL
    !
    REAL(r_std), PARAMETER                          :: R_Earth = 6378000., min_sechiba=1.E-8
    !
    !
    CHARACTER(LEN=*),INTENT(in) :: data_filename
    CHARACTER(LEN=*),INTENT(in) :: field_name
    INTEGER(i_std) :: iml, jml, ijml, i, j, ik, lml, tml, fid, ib, jb,ip, jp, vid, ai,iki,jkj
    INTEGER(i_std) :: nb_coord,nb_dim,nb_var,nb_gat
    LOGICAL :: l_ex
    INTEGER,DIMENSION(1)                        :: l_d_w
    REAL(r_std) :: lev(1), date, dt, coslat, pi
    INTEGER(i_std) :: itau(1)
    
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)                     ::  mask_lu
    REAL(r_std), ALLOCATABLE, DIMENSION(:)                       :: lat_lu, lon_lu, mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:)                       :: lat_ful, lon_ful
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)                     :: lightn_orig
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)                   :: lightn_lu !LOOP: as 4D needed?
    REAL(r_std), ALLOCATABLE, DIMENSION(:)                         :: lon_up, lon_low, lat_up, lat_low
    INTEGER, DIMENSION(nbpt) :: n_origlightn
    INTEGER, DIMENSION(nbpt) :: n_found
    
    CHARACTER(LEN=80) :: meter
    REAL(r_std) :: prog, sumf
    LOGICAL :: found
    INTEGER :: idi,jdi, ilast, jlast, jj, ii, jv, inear, iprog
    REAL(r_std) :: domaine_lon_min, domaine_lon_max, domaine_lat_min, domaine_lat_max
    !
    pi = 4. * ATAN(1.)
    !
    !Config Key  = Lightn_FILE
    !Config Desc = Name of file from which the lightn climatology is to be read
    !Config If   = !lightn
    !Config Def  = lightn_climatology_otd.nc
    !Config Help = The name of the file to be opened to read the lightnin flash rate
    !Config        map is to be given here. Resolution is 1.0x1.0deg.

    !
    CALL flininfo(data_filename, iml, jml, lml, tml, fid)

    !
    !
    ALLOCATE(lon_lu(iml))
    ALLOCATE(lat_lu(jml))
    ALLOCATE(lightn_lu(iml,jml,tml))
    ALLOCATE(mask_lu(iml,jml))
    !
    !
    WRITE(numout,*) 'input filename : ', data_filename
    CALL flinget(fid, 'lon', iml, 0, 0, 0, 1, 1, lon_lu)
    CALL flinget(fid, 'lat', jml, 0, 0, 0, 1, 1, lat_lu)
    CALL flinget(fid, field_name, iml, jml, 0, tml, 1, tml, lightn_lu)
    !
    WRITE(numout,*) 'cordinate information: ', lon_lu(1), lon_lu(iml),lat_lu(1), lat_lu(jml),tml

    !
    !
    ijml=iml*jml
    ALLOCATE(lon_ful(ijml))
    ALLOCATE(lat_ful(ijml))
    ALLOCATE(lightn_orig(ijml,tml))
    ALLOCATE(mask(ijml))


    DO i=1,iml
      DO j=1,jml
        iki=(j-1)*iml+i
        lon_ful(iki)=lon_lu(i)
        lat_ful(iki)=lat_lu(j)
        lightn_orig(iki,:)=lightn_lu(i,j,:)
        IF (lightn_lu(i,j,1).gt.-9000.)        &
           mask(iki) = 1.0
      ENDDO
    ENDDO

    !
    ALLOCATE(lon_up(nbpt)) 
    ALLOCATE(lon_low(nbpt))
    ALLOCATE(lat_up(nbpt))
    ALLOCATE(lat_low(nbpt))
    !
    DO ib =1, nbpt
      !
      !  We find the 4 limits of the grid-box. As we transform the resolution of the model
      !  into longitudes and latitudes we do not have the problem of periodicity.
      !  coslat is a help variable here !
      !
      coslat = MAX(COS(lalo(ib,1) * pi/180. ), 0.001 )*pi/180. * R_Earth
      !
      lon_up(ib) = lalo(ib,2) + resolution(ib,1)/(2.0*coslat) 
      lon_low(ib) = lalo(ib,2) - resolution(ib,1)/(2.0*coslat) 
      !
      coslat = pi/180. * R_Earth
      !
      lat_up(ib) = lalo(ib,1) + resolution(ib,2)/(2.0*coslat) 
      lat_low(ib) = lalo(ib,1) - resolution(ib,2)/(2.0*coslat) 
      !
      !
      !
    ENDDO
    !
    !  Get the limits of the integration domaine so that we can speed up the calculations
    !
    domaine_lon_min = MINVAL(lon_low)
    domaine_lon_max = MAXVAL(lon_up)
    domaine_lat_min = MINVAL(lat_low)
    domaine_lat_max = MAXVAL(lat_up)
    !
    ! Ensure that the fine grid covers the whole domain
    WHERE ( lon_ful(:) .LT. domaine_lon_min )
        lon_ful(:) = lon_ful(:) + 360.
    ENDWHERE
    !
    WHERE ( lon_ful(:) .GT. domaine_lon_max )
        lon_ful(:) = lon_ful(:) - 360.
    ENDWHERE
    !
    WRITE(numout,*) 'Interpolating the --',field_name,'-- map :'
    WRITE(numout,'(2a40)')'0%--------------------------------------', &
       & '------------------------------------100%'
    !
    ilast = 1
    n_origlightn(:) = 0.
    proxydata(:,:) = 0.    
    !
    DO ip=1,ijml
      !
      !   Give a progress meter
      !
      iprog = NINT(float(ip)/float(ijml)*79.) - NINT(float(ip-1)/float(ijml)*79.)
      IF ( iprog .NE. 0 ) THEN
          WRITE(numout,'(a1,$)') 'y'
      ENDIF
      !
      !  Only start looking for its place in the smaler grid if we are within the domaine
      !  That should speed up things !
      !
      IF ( ( lon_ful(ip) .GE. domaine_lon_min ) .AND. &
         ( lon_ful(ip) .LE. domaine_lon_max ) .AND. &
         ( lat_ful(ip) .GE. domaine_lat_min ) .AND. &
         ( lat_ful(ip) .LE. domaine_lat_max )        ) THEN
          !
          ! look for point on GCM grid which this point on fine grid belongs to.
          ! First look at the point on the model grid where we arrived just before. There is 
          ! a good chance that neighbouring points on the fine grid fall into the same model
          ! grid box.
          !
          IF ( ( lon_ful(ip) .GE. lon_low(ilast) ) .AND. &
             ( lon_ful(ip) .LT. lon_up(ilast) ) .AND. &
             ( lat_ful(ip) .GE. lat_low(ilast) ) .AND. &
             ( lat_ful(ip) .LT. lat_up(ilast) )         ) THEN
              !
              ! We were lucky
              !
              IF (mask(ip) .GT. 0) THEN 
                  n_origlightn(ilast) =  n_origlightn(ilast) + 1  
                  DO j=1,tml
                    proxydata(ilast,j) = proxydata(ilast,j) + lightn_orig(ip,j)
                  ENDDO
              ENDIF
              !
          ELSE
              !
              ! Otherwise, look everywhere.
              ! Begin close to last grid point.
              !
              found = .FALSE. 
              idi = 1
              !
              DO WHILE ( (idi .LT. nbpt) .AND. ( .NOT. found ) )

                !
                ! forward and backward
                !
                DO ii = 1,2
                  !
                  IF ( ii .EQ. 1 ) THEN
                      ib = ilast - idi
                  ELSE
                      ib = ilast + idi
                  ENDIF
                  !
                  IF ( ( ib .GE. 1 ) .AND. ( ib .LE. nbpt ) ) THEN 
                      IF ( ( lon_ful(ip) .GE. lon_low(ib) ) .AND. &
                         ( lon_ful(ip) .LT. lon_up(ib) ) .AND. &
                         ( lat_ful(ip) .GE. lat_low(ib) ) .AND. &
                         ( lat_ful(ip) .LT. lat_up(ib) )         ) THEN
                          !
                          IF (mask(ip) .gt. 0) THEN
                              !                            DO i=1,nvm
                              DO j=1,tml
                                proxydata(ib,j) = proxydata(ib,j) + lightn_orig(ip,j) 
                              ENDDO
                              !                            ENDDO
                              n_origlightn(ib) =  n_origlightn(ib) + 1
                          ENDIF
                          ilast = ib
                          found = .TRUE.
                          !
                      ENDIF
                  ENDIF
                  !
                ENDDO
                !
                idi = idi + 1
                !
              ENDDO
              !
          ENDIF ! lucky/not lucky
          !
      ENDIF     ! in the domain
    ENDDO

    WRITE(numout,*) ''

    ! determine fraction of lightning points in each box of the coarse grid
    DO ip=1,nbpt
      IF ( n_origlightn(ip) .GT. 0 ) THEN
          proxydata(ip,:) = proxydata(ip,:)/REAL(n_origlightn(ip),r_std)
      ELSE
          !
          WRITE(numout,*) 'PROBLEM, no point in the ba map found for this grid box'
          WRITE(numout,*) 'Longitude range : ', lon_low(ip), lon_up(ip)
          WRITE(numout,*) 'Latitude range : ', lat_low(ip), lat_up(ip)
          !
          WRITE(numout,*) 'Looking for nearest point on the ba map file'
          CALL slowproc_nearest (ijml, lon_ful, lat_ful, &
             lalo(ip,2), lalo(ip,1), inear)
          WRITE(numout,*) 'Coordinates of the nearest point, ',inear,' :', &
             lon_ful(inear),lat_ful(inear)
          proxydata(ip,:) = lightn_orig(inear,:)
      ENDIF
    ENDDO
    ! 
    IF (printlev_loc>=1) WRITE(numout,*) '--',field_name,'--: Interpolation Done'
    WRITE(numout,*) ''

    !
    DEALLOCATE(lon_up)
    DEALLOCATE(lon_low)
    DEALLOCATE(lat_up)
    DEALLOCATE(lat_low)
    DEALLOCATE(lat_ful)
    DEALLOCATE(lon_ful)
    DEALLOCATE(lat_lu)
    DEALLOCATE(lon_lu)
    DEALLOCATE(lightn_lu)
    DEALLOCATE(lightn_orig)
    DEALLOCATE(mask)
    !
    RETURN
    !
  END SUBROUTINE slowproc_read_data

  !--LOOP
  SUBROUTINE slowproc_read_annual(nbpt, lalo,  resolution, popd,data_filename,field_name)

    !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)          :: nbpt                  ! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)             :: lalo(nbpt,2)          ! Vector of latitude and longitudes (beware of the order !)
    REAL(r_std), INTENT(in)             :: resolution(nbpt,2)    ! The size in km of each grid-box in X and Y
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), INTENT(out)    ::  popd(nbpt)         ! population density  read variable and re-dimensioned
    !
    !  0.3 LOCAL
    !
    REAL(r_std), PARAMETER                          :: R_Earth = 6378000., min_sechiba=1.E-8
    !
    !
    CHARACTER(LEN=*),INTENT(in) :: data_filename
    CHARACTER(LEN=*),INTENT(in) :: field_name
    INTEGER(i_std) :: iml, jml, ijml, i, j, ik, lml, tml, fid, ib, jb,ip, jp, vid, ai,iki,jkj
    INTEGER(i_std) :: nb_coord,nb_dim,nb_var,nb_gat
    LOGICAL :: l_ex
    INTEGER,DIMENSION(1)                        :: l_d_w
    REAL(r_std) :: lev(1), date, dt, coslat, pi
    INTEGER(i_std) :: itau(1)

    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)                     ::  mask_lu
    REAL(r_std), ALLOCATABLE, DIMENSION(:)                       :: lat_lu, lon_lu, mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:)                       :: lat_ful, lon_ful
    REAL(r_std), ALLOCATABLE, DIMENSION(:)                     :: popden_orig
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)                   :: popden_lu !LOOP: as 4D needed?
    REAL(r_std), ALLOCATABLE, DIMENSION(:)                         :: lon_up, lon_low, lat_up, lat_low
    INTEGER, DIMENSION(nbpt) :: n_origpopden
    INTEGER, DIMENSION(nbpt) :: n_found

    CHARACTER(LEN=80) :: meter
    REAL(r_std) :: prog, sumf
    LOGICAL :: found
    INTEGER :: idi,jdi, ilast, jlast, jj, ii, jv, inear, iprog
    REAL(r_std) :: domaine_lon_min, domaine_lon_max, domaine_lat_min, domaine_lat_max
    !
    pi = 4. * ATAN(1.)

    !
    !Config Key  = POPDENS_FILE
    !Config Desc = Name of file from which the popden climatology is to be read
    !Config If   = !popden
    !Config Def  = popdens.nc
    !Config Help = The name of the file to be opened to read the population density
    !Config        map is to be given here. Resolution is 1.0x1.0deg.

    !
    CALL flininfo(data_filename, iml, jml, lml, tml, fid)
    
    !
    !
    ALLOCATE(lon_lu(iml))
    ALLOCATE(lat_lu(jml))
    ALLOCATE(popden_lu(iml,jml))
    ALLOCATE(mask_lu(iml,jml))
    !
    !
    WRITE(numout,*) 'input filename : ', data_filename
    CALL flinget(fid, 'lon', iml, 0, 0, 0, 1, 1, lon_lu)
    CALL flinget(fid, 'lat', jml, 0, 0, 0, 1, 1, lat_lu)
    CALL flinget(fid, field_name, iml, jml, 0, tml, 1, 1, popden_lu)

    !
    WRITE(numout,*) 'cordinate information: ', lon_lu(1), lon_lu(iml),lat_lu(1), lat_lu(jml),tml

    !
    !
    ijml=iml*jml
    ALLOCATE(lon_ful(ijml))
    ALLOCATE(lat_ful(ijml))
    ALLOCATE(popden_orig(ijml))
    ALLOCATE(mask(ijml))


    DO i=1,iml
       DO j=1,jml
          iki=(j-1)*iml+i
          lon_ful(iki)=lon_lu(i)
          lat_ful(iki)=lat_lu(j)
          popden_orig(iki)=popden_lu(i,j)
       IF (popden_lu(i,j).gt.-9000.)        &
          mask(iki) = 1.0
       ENDDO
    ENDDO

    !
    WHERE  ( popden_orig(:) .LT. 0 )
       popden_orig(:) = 0.
    ENDWHERE
    !
    !
    ALLOCATE(lon_up(nbpt)) 
    ALLOCATE(lon_low(nbpt))
    ALLOCATE(lat_up(nbpt))
    ALLOCATE(lat_low(nbpt))
    !
    DO ib =1, nbpt
       !
       !  We find the 4 limits of the grid-box. As we transform the resolution of the model
       !  into longitudes and latitudes we do not have the problem of periodicity.
       !  coslat is a help variable here !
       !
       coslat = MAX(COS(lalo(ib,1) * pi/180. ), 0.001 )*pi/180. * R_Earth
       !
       lon_up(ib) = lalo(ib,2) + resolution(ib,1)/(2.0*coslat) 
       lon_low(ib) = lalo(ib,2) - resolution(ib,1)/(2.0*coslat) 
       !
       coslat = pi/180. * R_Earth
       !
       lat_up(ib) = lalo(ib,1) + resolution(ib,2)/(2.0*coslat) 
       lat_low(ib) = lalo(ib,1) - resolution(ib,2)/(2.0*coslat) 
       !
       !
       !
    ENDDO
    !
    !  Get the limits of the integration domaine so that we can speed up the calculations
    !
    domaine_lon_min = MINVAL(lon_low)
    domaine_lon_max = MAXVAL(lon_up)
    domaine_lat_min = MINVAL(lat_low)
    domaine_lat_max = MAXVAL(lat_up)
    !
    !
    ! Ensure that the fine grid covers the whole domain
    WHERE ( lon_ful(:) .LT. domaine_lon_min )
      lon_ful(:) = lon_ful(:) + 360.
    ENDWHERE
    !
    WHERE ( lon_ful(:) .GT. domaine_lon_max )
      lon_ful(:) = lon_ful(:) - 360.
    ENDWHERE
    !
    WRITE(numout,*) 'Interpolating the --',field_name,'-- map :'
    WRITE(numout,'(2a40)')'0%--------------------------------------', &
                   & '------------------------------------100%'
    !
    ilast = 1
    n_origpopden(:) = 0
    popd(:) = 0.    
    !
    DO ip=1,ijml
       !
       !   Give a progress meter
       !
       iprog = NINT(float(ip)/float(ijml)*79.) - NINT(float(ip-1)/float(ijml)*79.)
       IF ( iprog .NE. 0 ) THEN
          WRITE(numout,'(a1,$)') 'y'
       ENDIF
       !
       !  Only start looking for its place in the smaler grid if we are within the domaine
       !  That should speed up things !
       !
       IF ( ( lon_ful(ip) .GE. domaine_lon_min ) .AND. &
            ( lon_ful(ip) .LE. domaine_lon_max ) .AND. &
            ( lat_ful(ip) .GE. domaine_lat_min ) .AND. &
            ( lat_ful(ip) .LE. domaine_lat_max )        ) THEN
          !
          ! look for point on GCM grid which this point on fine grid belongs to.
          ! First look at the point on the model grid where we arrived just before. There is 
          ! a good chance that neighbouring points on the fine grid fall into the same model
          ! grid box.
          !
          IF ( ( lon_ful(ip) .GE. lon_low(ilast) ) .AND. &
               ( lon_ful(ip) .LT. lon_up(ilast) ) .AND. &
               ( lat_ful(ip) .GE. lat_low(ilast) ) .AND. &
               ( lat_ful(ip) .LT. lat_up(ilast) )         ) THEN
             !
             ! We were lucky
             !
             IF (mask(ip) .GT. 0) THEN 
               n_origpopden(ilast) =  n_origpopden(ilast) + 1  
               popd(ilast) = popd(ilast) + popden_orig(ip)
             ENDIF
             !
          ELSE
             !
             ! Otherwise, look everywhere.
             ! Begin close to last grid point.
             !
             found = .FALSE. 
             idi = 1
             !
             DO WHILE ( (idi .LT. nbpt) .AND. ( .NOT. found ) )

                !
                ! forward and backward
                !
                DO ii = 1,2
                   !
                   IF ( ii .EQ. 1 ) THEN
                      ib = ilast - idi
                   ELSE
                      ib = ilast + idi
                   ENDIF
                   !
                   IF ( ( ib .GE. 1 ) .AND. ( ib .LE. nbpt ) ) THEN 
                      IF ( ( lon_ful(ip) .GE. lon_low(ib) ) .AND. &
                           ( lon_ful(ip) .LT. lon_up(ib) ) .AND. &
                           ( lat_ful(ip) .GE. lat_low(ib) ) .AND. &
                           ( lat_ful(ip) .LT. lat_up(ib) )         ) THEN
                         !
                         IF (mask(ip) .gt. 0) THEN
                            popd(ib) = popd(ib) + popden_orig(ip) 
                            n_origpopden(ib) =  n_origpopden(ib) + 1
                         ENDIF
                         ilast = ib
                         found = .TRUE.
                         !
                      ENDIF
                   ENDIF
                   !
                ENDDO
                !
                idi = idi + 1
                !
             ENDDO
             !
          ENDIF ! lucky/not lucky
          !
       ENDIF     ! in the domain
    ENDDO

    WRITE(numout,*) ''

    ! determine fraction of popdening points in each box of the coarse grid
    DO ip=1,nbpt
       IF ( n_origpopden(ip) .GT. 0 ) THEN
             popd(ip) = popd(ip)/REAL(n_origpopden(ip),r_std)
       ELSE
             WRITE(numout,*) 'PROBLEM, no point in the popdens map found for this grid box'
             WRITE(numout,*) 'Longitude range : ', lon_low(ip), lon_up(ip)
             WRITE(numout,*) 'Latitude range : ', lat_low(ip), lat_up(ip)
             !
             WRITE(numout,*) 'Looking for nearest point on the popdens map file'
             CALL slowproc_nearest (ijml, lon_ful, lat_ful, &
                  lalo(ip,2), lalo(ip,1), inear)
             WRITE(numout,*) 'Coordinates of the nearest point, ',inear,' :', &
                  lon_ful(inear),lat_ful(inear)
             !
                popd(ip) = popden_orig(inear)
       ENDIF
    ENDDO
    ! 
    IF (printlev_loc>=1) WRITE(numout,*) '--',field_name,'--: Interpolation Done'
    WRITE(numout,*) ''
   
    !
    DEALLOCATE(lon_up)
    DEALLOCATE(lon_low)
    DEALLOCATE(lat_up)
    DEALLOCATE(lat_low)
    DEALLOCATE(lat_ful)
    DEALLOCATE(lon_ful)
    DEALLOCATE(lat_lu)
    DEALLOCATE(lon_lu)
    DEALLOCATE(popden_lu)
    DEALLOCATE(popden_orig)
    DEALLOCATE(mask)
    !
    RETURN
    !
  END SUBROUTINE slowproc_read_annual

  !endspit

!================================================================================================================================
!! SUBROUTINE   : slowproc_read_wetland
!!
!>\BRIEF          
!!
!! DESCRIPTION  : Read file of wetland map
!!                
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): Distribution of bogs, fens and mineral wetland
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_
!================================================================================================================================
  SUBROUTINE slowproc_read_wetland(nbpt, lalo,resolution,wettile_dgvm,tile_name_dgvm)

    INTEGER(i_std), INTENT(in)                    :: nbpt                  !! Number of points for which the data needs to be interpolated (unitless)             
    REAL(r_std), INTENT(in)                       :: lalo(nbpt,2)          !! Vector of latitude and longitudes (degree)        
    REAL(r_std), INTENT(in)                       :: resolution(nbpt,2)    !! The size of each grid cell in X and Y (km)
    LOGICAL,DIMENSION (nstm), INTENT (in)                 :: wettile_dgvm
    CHARACTER(LEN=25), DIMENSION(nstm),INTENT (in)        :: tile_name_dgvm
    REAL(r_std), DIMENSION(nbpt)                  :: fpeat_map_fen
    REAL(r_std), DIMENSION(nbpt)                  :: fpeat_map_bog
    REAL(r_std), DIMENSION(nbpt)                  :: fpeat_map_mineral

    REAL(r_std), PARAMETER                        :: R_Earth = 6378000.
    REAL(r_std), ALLOCATABLE, DIMENSION (:,:)     ::bog_map
    REAL(r_std), ALLOCATABLE, DIMENSION (:,:)     ::fen_map
    REAL(r_std), ALLOCATABLE, DIMENSION (:,:)     ::mineral_map
    CHARACTER(LEN=80)                             :: filenamenew
    REAL(r_std), ALLOCATABLE, DIMENSION(:)              ::bog_map_orig
    REAL(r_std), ALLOCATABLE, DIMENSION(:)              ::fen_map_orig
    REAL(r_std), ALLOCATABLE, DIMENSION(:)              ::mineral_map_orig
    INTEGER                                       :: ALLOC_ERR !!Help varialbe to count allocation error
    INTEGER(i_std) :: iml, jml, ijml, i, j, ik, lml, tml, fid, ib, ip,  iki, jst
    INTEGER(i_std) :: nwet,iFen,iBog,iMineral,iWetland
    REAL(r_std) :: coslat, pi
    REAL(r_std), ALLOCATABLE, DIMENSION(:)                       :: lat_lu,lon_lu
    REAL(r_std), ALLOCATABLE, DIMENSION(:)                       :: lat_ful,lon_ful
    REAL(r_std), ALLOCATABLE, DIMENSION(:)                         :: lon_up,lon_low, lat_up, lat_low
    INTEGER, DIMENSION(nbpt) :: n_grid

    LOGICAL :: found
    INTEGER :: idi, ilast, ii, inear, iprog
    REAL(r_std) :: domaine_lon_min, domaine_lon_max, domaine_lat_min,domaine_lat_max
    REAL(r_std)                                  :: pa,cospa, sinpa
    REAL(r_std), ALLOCATABLE, DIMENSION(:)       :: cosang
!_
!================================================================================================================================

  filenamenew = 'wetland_type_file.nc'
  CALL getin_p('WETLAND_TYPE_FILE',filenamenew)

  pi = 4. * ATAN(1.)

  IF (is_root_prc) CALL flininfo(filenamenew,iml, jml, lml, tml, fid)
  CALL bcast(iml)
  CALL bcast(jml)
  CALL bcast(lml)
  CALL bcast(tml)

  ALLOCATE(lon_lu(iml), STAT=ALLOC_ERR)
  IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_read_wetland','Problem in allocation of variable lon_lu','','')
  ALLOCATE(lat_lu(jml), STAT=ALLOC_ERR)
  IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_read_wetland','Problem in allocation of variable lat_lu','','')

  ALLOCATE(bog_map(iml,jml), STAT=ALLOC_ERR)
  IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_read_wetland','Pb in allocation for bog_map','','')
  ALLOCATE(fen_map(iml,jml), STAT=ALLOC_ERR)
  IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_read_wetland','Pb in allocation for fen_map','','')
  ALLOCATE(mineral_map(iml,jml), STAT=ALLOC_ERR)
  IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_read_wetland','Pb in allocation for mineral_map','','')



  IF (is_root_prc) THEN
     CALL flinget(fid, 'longitude', iml, 0, 0, 0, 1, 1, lon_lu)
     CALL flinget(fid, 'latitude', jml, 0, 0, 0, 1, 1, lat_lu)
     !CALL flinget(fid, 'Fens', iml, jml, lml, tml, 1, 1, fen_map)
     !CALL flinget(fid, 'Bogs', iml, jml, lml, tml, 1, 1, bog_map)
     !CALL flinget(fid, 'MineralWet', iml, jml, lml, tml, 1, 1, mineral_map)
     CALL flinget(fid, 'Fens', iml, jml, 0, 0, 1, 1, fen_map)
     CALL flinget(fid, 'Bogs', iml, jml, 0, 0, 1, 1, bog_map)
     CALL flinget(fid, 'MineralWet', iml, jml, 0, 0, 1, 1, mineral_map)
     CALL flinclo(fid)
  ENDIF

  CALL bcast(lon_lu)
  CALL bcast(lat_lu)
  CALL bcast(fen_map)
  CALL bcast(bog_map)
  CALL bcast(mineral_map)

  ijml=iml*jml
  ALLOCATE(lon_ful(ijml))
  ALLOCATE(lat_ful(ijml))
  ALLOCATE(fen_map_orig(ijml))
  ALLOCATE(bog_map_orig(ijml))
  ALLOCATE(mineral_map_orig(ijml))

  DO i=1,iml
    DO j=1,jml
      iki=(j-1)*iml+i
      lon_ful(iki)=lon_lu(i)
      lat_ful(iki)=lat_lu(j)
      fen_map_orig(iki)=fen_map(i,j)
      bog_map_orig(iki)=bog_map(i,j)
      mineral_map_orig(iki)=mineral_map(i,j)
    ENDDO
  ENDDO

  ALLOCATE(lon_up(nbpt))
  ALLOCATE(lon_low(nbpt))
  ALLOCATE(lat_up(nbpt))
  ALLOCATE(lat_low(nbpt))

  DO ib =1, nbpt
  !  We find the 4 limits of the grid-box. As we transform the resolution of
  !  the model into longitudes and latitudes we do not have the problem of
  !  periodicity. coslat is a help variable here !

     coslat = MAX(COS(lalo(ib,1) * pi/180. ), 0.001 )*pi/180. * R_Earth
      !
     lon_up(ib) = lalo(ib,2) + resolution(ib,1)/(2.0*coslat)
     lon_low(ib) = lalo(ib,2) - resolution(ib,1)/(2.0*coslat)
      !
     coslat = pi/180. * R_Earth
      !
     lat_up(ib) = lalo(ib,1) + resolution(ib,2)/(2.0*coslat)
     lat_low(ib) = lalo(ib,1) - resolution(ib,2)/(2.0*coslat)
  ENDDO
  !  Get the limits of the integration domaine so that we can speed up the calculations
    domaine_lon_min = MINVAL(lon_low)
    domaine_lon_max = MAXVAL(lon_up)
    domaine_lat_min = MINVAL(lat_low)
    domaine_lat_max = MAXVAL(lat_up)
    !
    ! Ensure that the fine grid covers the whole domain
    WHERE ( lon_ful(:) .LT. domaine_lon_min )
        lon_ful(:) = lon_ful(:) + 360.
    ENDWHERE
    !
    WHERE ( lon_ful(:) .GT. domaine_lon_max )
        lon_ful(:) = lon_ful(:) - 360.
    ENDWHERE
  ilast = 1
  n_grid(:) = 0.
  fpeat_map(:,:)=zero
  fpeat_map_fen(:)=zero
  fpeat_map_bog(:)=zero
  fpeat_map_mineral(:)=zero
  DO ip=1,ijml
      !   Give a progress meter
    iprog = NINT(float(ip)/float(ijml)*79.) - NINT(float(ip-1)/float(ijml)*79.)
      !  Only start looking for its place in the smaler grid if we are within
      !  the domaine
      !  That should speed up things !
    IF ( ( lon_ful(ip) .GE. domaine_lon_min ) .AND. &
       ( lon_ful(ip) .LE. domaine_lon_max ) .AND. &
       ( lat_ful(ip) .GE. domaine_lat_min ) .AND. &
       ( lat_ful(ip) .LE. domaine_lat_max )        ) THEN
      ! look for point on GCM grid which this point on fine grid belongs to.
      ! First look at the point on the model grid where we arrived just before. There is 
      ! a good chance that neighbouring points on the fine grid fall into the same model
      ! grid box.
       IF ( ( lon_ful(ip) .GE. lon_low(ilast) ) .AND. &
          ( lon_ful(ip) .LT. lon_up(ilast) ) .AND. &
          ( lat_ful(ip) .GE. lat_low(ilast) ) .AND. &
          ( lat_ful(ip) .LT. lat_up(ilast) )         ) THEN

             n_grid(ilast) =  n_grid(ilast) + 1
             fpeat_map_fen(ilast)=fpeat_map_fen(ilast)+fen_map_orig(ip)
             fpeat_map_bog(ilast)=fpeat_map_bog(ilast)+bog_map_orig(ip)
             fpeat_map_mineral(ilast)=fpeat_map_mineral(ilast)+mineral_map_orig(ip)
       ELSE   
              ! Otherwise, look everywhere.             
              ! Begin close to last grid point.
          found = .FALSE.
          idi = 1
          DO WHILE ( (idi .LT. nbpt) .AND. ( .NOT. found ) )
          
             ! forward and backward
            DO ii = 1,2
               IF ( ii .EQ. 1 ) THEN
                  ib = ilast - idi
               ELSE
                  ib = ilast + idi
               ENDIF
               IF ( ( ib .GE. 1 ) .AND. ( ib .LE. nbpt ) ) THEN
                  IF ( ( lon_ful(ip) .GE. lon_low(ib) ) .AND. &
                     ( lon_ful(ip) .LT. lon_up(ib) ) .AND. &
                     ( lat_ful(ip) .GE. lat_low(ib) ) .AND. &
                     ( lat_ful(ip) .LT. lat_up(ib) )         ) THEN

                        n_grid(ib) =  n_grid(ib) + 1
                        fpeat_map_fen(ib)=fpeat_map_fen(ib)+ fen_map_orig(ip)
                        fpeat_map_bog(ib)=fpeat_map_bog(ib)+ bog_map_orig(ip)
                        fpeat_map_mineral(ib)=fpeat_map_mineral(ib)+ mineral_map_orig(ip)
                     ilast = ib
                     found = .TRUE.
                  ENDIF
               ENDIF
            ENDDO
            idi = idi + 1
          ENDDO
       ENDIF ! lucky/not lucky
    ENDIF     ! in the domain
  ENDDO

  DO ip=1,nbpt
    IF ( n_grid(ip) .GT. 0 ) THEN
       fpeat_map_fen(ip) = fpeat_map_fen(ip)/REAL(n_grid(ip),r_std)
       fpeat_map_bog(ip) = fpeat_map_bog(ip)/REAL(n_grid(ip),r_std)
       fpeat_map_mineral(ip) = fpeat_map_mineral(ip)/REAL(n_grid(ip),r_std)
    ELSE
       WRITE(numout,*) 'PROBLEM, no point in the WETLAND_TYPE_FILE found for this gridbox'
       WRITE(numout,*) 'Location : ', lalo(ip,2), lalo(ip,1)
        fpeat_map_fen(ip)= zero
        fpeat_map_bog(ip)= zero
        fpeat_map_mineral(ip)= zero
    ENDIF
  ENDDO

  nwet=0
  iFen=0
  iBog=0
  iMineral=0
  DO jst=1,nstm
     IF (wettile_dgvm(jst)) THEN
        nwet=nwet+1 
     ENDIF
     IF (tile_name_dgvm(jst)=='Fen' .OR. tile_name_dgvm(jst)=='Peatland') THEN
        iFen=jst
     ENDIF
     IF (tile_name_dgvm(jst)=='Bog') THEN
        iBog=jst
     ENDIF
     IF (tile_name_dgvm(jst)=='Mineral') THEN
        iMineral=jst
     ENDIF
  ENDDO

  iWetland=iFen+iBog+iMineral
  IF (iWetland .LE. 0) THEN
    CALL ipslerr_p(3,'slowproc_read_wetland','When the peatland module is activated, there must be at least one baresoil soiltile','','')
  ENDIF

  IF (nwet==1) THEN
    WRITE (numout,*) 'NOTE, nwet=1, Peatland without distinct types'
  ELSEIF (nwet==2) THEN 
    WRITE (numout,*) 'NOTE, nwet=2, Check which two wetland types among Fen, Bog, Peatland(Fen+Bog), and Mineral are being simulated'
  ELSEIF (nwet==3) THEN
    WRITE (numout,*) 'NOTE, nwet=3, Fen, Bog, Mineral'
  ELSE
    CALL ipslerr_p(3,'slowproc_read_wetland','Not enough wetland types when calculating fpeat_map','','')
  ENDIF 

  DO ip=1,nbpt
    IF (nwet==1) THEN
      fpeat_map(ip,iFen)=fpeat_map_fen(ip)+fpeat_map_bog(ip)
    ELSEIF (nwet==2) THEN 
      DO jst=1,nstm
        IF (tile_name_dgvm(jst)=='Fen') THEN
           fpeat_map(ip,jst)=fpeat_map_fen(ip)
        ELSEIF (tile_name_dgvm(jst)=='Peatland') THEN
           fpeat_map(ip,jst)=fpeat_map_fen(ip)+fpeat_map_bog(ip)
        ELSEIF (tile_name_dgvm(jst)=='Bog') THEN
           fpeat_map(ip,jst)=fpeat_map_bog(ip)
        ELSEIF (tile_name_dgvm(jst)=='Mineral') THEN
           fpeat_map(ip,jst)=fpeat_map_mineral(ip)
        ENDIF
      ENDDO    
    ELSEIF (nwet==3) THEN
       fpeat_map(ip,iFen)=fpeat_map_fen(ip)
       fpeat_map(ip,iBog)=fpeat_map_bog(ip)
       fpeat_map(ip,iMineral)=fpeat_map_mineral(ip)
    ELSE
      CALL ipslerr_p(3,'slowproc_read_wetland','Not enough wetland types when calculating fpeat_map','','')
    ENDIF
  ENDDO

  DEALLOCATE(lat_lu)
  DEALLOCATE(lon_lu)
  DEALLOCATE (bog_map)
  DEALLOCATE (fen_map)
  DEALLOCATE (mineral_map)
  DEALLOCATE(lat_ful)
  DEALLOCATE(lon_ful)
  DEALLOCATE (bog_map_orig)
  DEALLOCATE (fen_map_orig)
  DEALLOCATE (mineral_map_orig)
  DEALLOCATE(lon_up)
  DEALLOCATE(lon_low)
  DEALLOCATE(lat_up)
  DEALLOCATE(lat_low)

  END SUBROUTINE slowproc_read_wetland


END MODULE slowproc
