! =================================================================================================================================
! MODULE        : stomate_season
!
! CONTACT       : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE       : IPSL (2006). This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       This module calculates long-term meteorological parameters from daily temperatures
!! and precipitations (essentially for phenology). 
!!      
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_stomate/stomate_season.f90 $ 
!! $Date: 2017-10-26 15:35:04 +0200 (Thu, 26 Oct 2017) $
!! $Revision: 4719 $
!! \n
!_ =================================================================================================================================

MODULE stomate_season

  ! modules used:
  USE xios_orchidee
  USE ioipsl_para
  USE stomate_data
  USE constantes
  USE constantes_soil
  USE pft_parameters 
  USE solar
  USE grid
  USE time, ONLY : one_year, julian_diff, LastTsYear

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC season,season_clear

  LOGICAL, SAVE              :: firstcall_season = .TRUE.  !! first call (true/false)
!$OMP THREADPRIVATE(firstcall_season)

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE   : season_clear
!!
!>\BRIEF          Flag setting 
!!
!! DESCRIPTION  : This subroutine sets flags ::firstcall_season, to .TRUE., and therefore activates   
!!                section 1.1 of the ::season subroutine which writes messages to the output. \n
!!                This subroutine is called at the end of the subroutine ::stomate_clear, in the 
!!                module ::stomate.
!!
!! RECENT CHANGE(S):None
!!
!! MAIN OUTPUT VARIABLE(S): ::firstcall_season
!!
!! REFERENCE(S)  : None 
!!
!! FLOWCHART     : None
!! \n             
!_ =================================================================================================================================

  SUBROUTINE season_clear
    firstcall_season=.TRUE.
  END SUBROUTINE season_clear



!! ================================================================================================================================
!! SUBROUTINE   : season
!!
!>\BRIEF          This subroutine calculates many of the long-term biometeorological variables
!!                needed in the phenology subroutines and in the calculation of long-term vegetation
!!                dynamics in the LPJ DGVM. 
!!
!! DESCRIPTION    This subroutine is called by the module ::stomate before LPJ, and mainly deals  
!!                with the calculation of long-term meteorological variables, carbon fluxes and 
!!                vegetation-related variables that are used to calculate vegetation dynamics 
!!                in the stomate modules relating to the phenology and to the longer-term 
!!                changes in vegetation type and fractional cover in the LPJ DGVM modules. \n
!!                In sections 2 to 5, longer-term meteorological variables are calculated. The
!!                long-term moisture availabilities are used in the leaf onset and senescence
!!                phenological models ::stomate_phenology and ::stomate_turnover that require
!!                a moisture condition. The long term temperatures are also required for phenology
!!                but in addition they are used in calculations of C flux and the presence and
!!                establishment of vegetation patterns on a longer timescale in the LPJ DGVM modules.
!!                Finally the monthly soil humidity/relative soil moisture is used in the C
!!                allocation module ::stomate_alloc. \n
!!                Sections 12 to 14 also calculate long-term variables of C fluxes, including NPP,
!!                turnover and GPP. The first two are used to calculate long-term vegetation 
!!                dynamics and land cover change in the LPJ DVGM modules. The weekly GPP is used
!!                to determine the dormancy onset and time-length, as described below. \n
!!                The long-term variables described above are are used in certain vegetation
!!                dynamics processes in order to maintain consistency with the basic hypotheses of the 
!!                parameterisations of LPJ, which operates on a one year time step (Krinner et al., 2005).
!!                In order to reduce the computer memory requirements, short-term variables (e.g. daily
!!                temperatures) are not stored for averaging over a longer period. Instead
!!                the long-term variables (e.g. monthly temperature) are calculated at each time step
!!                using a linear relaxation method, following the equation:
!!                \latexonly
!!                \input{season_lin_relax_eqn1.tex}
!!                \endlatexonly
!!                \n 
!!                The long-term variables are therefore updated on a daily basis. This method allows for 
!!                smooth temporal variation of the long-term variables which are used to calculate 
!!                vegetation dynamics (Krinner et al., 2005). \n
!!                Sections 6 to 11 calculate the variables required for to determine leaf onset in
!!                the module ::stomate_phenology. 
!!                These include :
!!                - the dormance onset and time-length, when GPP is below a certain threshold \n
!!                - the growing degree days (GDD), which is the sum of daily temperatures greater than
!!                -5 degrees C, either since the onset of the dormancy period, or since midwinter \n
!!                - the number of chilling days, which is the number of days with a daily temperature
!!                lower than a PFT-dependent threshold since the beginning of the dormancy period \n
!!                - the number of growing days, which is the the number of days with a temperature
!!                greater than -5 degrees C since the onset of the dormancy period \n
!!                - the time since the minimum moisture availability in the dormancy period. \n
!!                These variables are used to determine the conditions needed for the start of the
!!                leaf growing season. The specific models to which they correspond are given below. \n
!!                Sections 15 to 20 are used to update the maximum/minimum or the sum of various 
!!                meteorological and biological variables during the year that are required for
!!                calculating either leaf onset or longer-term vegetation dynamics the following year. \n
!!                At the end of the year, these variables are updated from "thisyear" to "lastyear", 
!!                in Section 21 of this subroutine, for use the following year. \n
!!                Finally the probably amount of herbivore consumption is calculated in Section 22,
!!                following McNaughton et al. (1989).
!!
!! RECENT CHANGE(S): None
!!                
!! MAIN OUTPUT VARIABLE(S): :: herbivores,
!!                        :: maxmoiavail_lastyear, :: maxmoiavail_thisyear, 
!!                        :: minmoiavail_lastyear, :: minmoiavail_thisyear, 
!!                        :: maxgppweek_lastyear, :: maxgppweek_thisyear, 
!!                        :: gdd0_lastyear, :: gdd0_thisyear, 
!!                        :: precip_lastyear, :: precip_thisyear, 
!!                        :: lm_lastyearmax, :: lm_thisyearmax, 
!!                        :: maxfpc_lastyear, :: maxfpc_thisyear, 
!!                        :: moiavail_month, :: moiavail_week, 
!!                        :: t2m_longterm, :: t2m_month, :: t2m_week, 
!!                        :: tsoil_month, :: soilhum_month, 
!!                        :: npp_longterm, :: turnover_longterm, :: gpp_week, 
!!                        :: gdd_m5_dormance, :: gdd_midwinter, 
!!                        :: ncd_dormance, :: ngd_minus5, :: time_lowgpp, 
!!                        :: time_hum_min, :: hum_min_dormance 
!!      
!! REFERENCES   : 
!! - Krinner, G., N. Viovy, N. de Noblet-Ducoudre, J. Ogee, J. Polcher, P. 
!! Friedlingstein, P. Ciais, S. Sitch and I.C. Prentice (2005), A dynamic global
!! vegetation model for studies of the coupled atmosphere-biosphere system, Global
!! Biogeochemical Cycles, 19, doi:10.1029/2003GB002199.
!! - McNaughton, S.J., M. Oesterheld, D.A. Frank and K.J. Williams (1989), 
!! Ecosystem-level patterns of primary productivity and herbivory in terrestrial
!! habitats, Nature, 341, 142-144.
!!                       
!! FLOWCHART    : 
!! \latexonly
!! \includegraphics[scale = 1]{season_flowchart_part1.png}
!! \includegraphics[scale = 1]{season_flowchart_part2.png}
!! \includegraphics[scale = 1]{season_flowchart_part3.png}
!! \endlatexonly
!! \n   
!_ =================================================================================================================================

  SUBROUTINE season (npts, dt, veget_cov, veget_cov_max, &
       moiavail_daily, t2m_daily, &
       tsurf_daily, & !pss:+-
       tsoil_daily, soilhum_daily, lalo,  &
       precip_daily, npp_daily, biomass, turnover_daily, gpp_daily, when_growthinit, &
       maxmoiavail_lastyear, maxmoiavail_thisyear, &
       minmoiavail_lastyear, minmoiavail_thisyear, &
       maxgppweek_lastyear, maxgppweek_thisyear, &
       gdd0_lastyear, gdd0_thisyear, &
       precip_lastyear, precip_thisyear, &
       lm_lastyearmax, lm_thisyearmax, &
       maxfpc_lastyear, maxfpc_thisyear, &
       moiavail_month, moiavail_week, t2m_longterm, tau_longterm, t2m_month, t2m_week, &
       tsurf_year, & !pss:+-
       tsoil_month, soilhum_month, &
       npp_longterm, turnover_longterm, gpp_week, &
       gdd_m5_dormance, gdd_midwinter, ncd_dormance, ngd_minus5, &
       time_hum_min, hum_min_dormance, gdd_init_date , & 
       gdd_from_growthinit, herbivores, Tseason, Tseason_length, Tseason_tmp, &
       Tmin_spring_time, t2m_min_daily, begin_leaves, onset_date, &!)
!gmjc
       t2m_14, sla_calc, &
!end gmjc
       npp0_cumul,fpeat_map,wtp_daily,wtp_month,wtp_year) !!!qcj++ peatland

    !
    !! 0. Variable and parameter declaration
    !

    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                             :: npts              !! Domain size - number of grid cells (unitless)
    REAL(r_std), INTENT(in)                                :: dt                !! time step in days (dt_days)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)           :: veget_cov         !! coverage fraction of a PFT. Here: fraction of 
                                                                                !! total ground. (0-1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)           :: veget_cov_max     !! "maximal" coverage fraction of a PFT (for LAI -> 
                                                                                !! infinity) (0-1, unitless)Here: fraction of 
                                                                                !! total ground. 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)           :: moiavail_daily    !! Daily moisture availability (0-1, unitless)
    REAL(r_std), DIMENSION(npts), INTENT(in)               :: t2m_daily         !! Daily 2 meter temperature (K)
    REAL(r_std), DIMENSION(npts,nslm), INTENT(in)          :: tsoil_daily       !! Daily soil temperature (K)
    REAL(r_std), DIMENSION(npts,nslm), INTENT(in)          :: soilhum_daily     !! Daily soil humidity (0-1, unitless)
    REAL(r_std), DIMENSION(npts,2), INTENT(in)             :: lalo              !!  array of lat/lon
    REAL(r_std), DIMENSION(npts), INTENT(in)               :: precip_daily      !! Daily mean precipitation @tex ($mm day^{-1}$) 
                                                                                !! @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)           :: npp_daily         !! daily net primary productivity @tex ($gC m^{-2} 
                                                                                !! day^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(in) :: biomass     !! biomass @tex ($gC m^{-2} of ground$) @endtex
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(in) :: turnover_daily  !! Turnover rates @tex ($gC m^{-2} day^{-1}$) 
                                                                                !! @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)           :: gpp_daily         !! daily gross primary productivity  
                                                                                !! (Here: @tex $gC m^{-2} of total ground 
                                                                                !! day^{-1}$) @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)           :: when_growthinit   !! how many days ago was the beginning of the 
                                                                                !! growing season (days) 
    REAL(r_std), DIMENSION(npts), INTENT(in)           :: tsurf_daily !!pss:+ daily surface temperature

    LOGICAL, DIMENSION(npts,nvm), INTENT(in)               :: begin_leaves      !! signal to start putting leaves on (true/false)

    REAL(r_std), DIMENSION(npts), INTENT(in)               :: t2m_min_daily     !! Daily minimum 2-meter temperature (K)
 

    !
    !! 0.2 Output variables 
    ! (diagnostic)
    !
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)          :: herbivores        !! time constant of probability of a leaf to be 
                                                                                !! eaten by a herbivore (days) 

    !
    !! 0.3 Modified variables
    !
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: maxmoiavail_lastyear      !! last year's maximum moisture 
                                                                                        !! availability (0-1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: maxmoiavail_thisyear      !! this year's maximum moisture 
                                                                                        !! availability (0-1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: minmoiavail_lastyear      !! last year's minimum moisture 
                                                                                        !! availability (0-1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: minmoiavail_thisyear      !! this year's minimum moisture 
                                                                                        !! availability (0-1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: maxgppweek_lastyear       !! last year's maximum weekly GPP 
                                                                                        !! @tex ($gC m^{-2} week^{-1}$) @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: maxgppweek_thisyear       !! this year's maximum weekly GPP 
                                                                                        !! @tex ($gC m^{-2} week^{-1}$) @endtex 
    REAL(r_std), DIMENSION(npts), INTENT(inout)            :: gdd0_lastyear             !! last year's annual GDD0 (C)
    REAL(r_std), DIMENSION(npts), INTENT(inout)            :: gdd0_thisyear             !! this year's annual GDD0 (C)
    REAL(r_std), DIMENSION(npts), INTENT(inout)            :: precip_lastyear           !! last year's annual precipitation 
                                                                                        !! @tex ($mm year^{-1}$) @endtex 
    REAL(r_std), DIMENSION(npts), INTENT(inout)            :: precip_thisyear           !! this year's annual precipitation 
                                                                                        !! @tex ($mm year^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: lm_lastyearmax            !! last year's maximum leaf mass, for each 
                                                                                        !! PFT @tex ($gC m^{-2}$) @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: lm_thisyearmax            !! this year's maximum leaf mass, for each 
                                                                                        !! PFT @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: maxfpc_lastyear           !! last year's maximum fpc for each PFT, on 
                                                                                        !! ground (0-1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: maxfpc_thisyear           !! this year's maximum fpc for each PFT, on 
                                                                                        !! ground (0-1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: moiavail_month            !! "monthly" moisture availability 
                                                                                        !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: moiavail_week             !! "weekly" moisture availability 
                                                                                        !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts), INTENT(inout)            :: t2m_longterm              !! "long term" 2-meter temperatures (K)
    REAL(r_std), INTENT(inout)                             :: tau_longterm     
    REAL(r_std), DIMENSION(npts), INTENT(inout)            :: t2m_month                 !! "monthly" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(inout)            :: t2m_week                  !! "weekly" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(inout)            :: Tseason                   !! "seasonal" 2-meter temperatures (K), 
                                                                                        !! used to constrain boreal treeline
    REAL(r_std), DIMENSION(npts), INTENT(inout)            :: Tseason_tmp               !! temporary variable to calculate Tseason 
    REAL(r_std), DIMENSION(npts), INTENT(inout)            :: Tseason_length            !! temporary variable to calculate Tseason: 
                                                                                        !! number of days when t2m_week is higher than 0 degree 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: Tmin_spring_time          !! Number of days after begin_leaves (leaf onset)

    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: onset_date                !! Date in the year at when the leaves started to grow(begin_leaves)

    REAL(r_std), DIMENSION(npts,nslm), INTENT(inout)       :: tsoil_month               !! "monthly" soil temperatures (K)
    REAL(r_std), DIMENSION(npts,nslm), INTENT(inout)       :: soilhum_month             !! "monthly" soil humidity (0-1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: npp_longterm              !! "long term" net primary productivity 
                                                                                        !! @tex ($gC m^{-2} year^{-1}$) @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: turnover_longterm !! "long term" turnover rate 
                                                                                        !! @tex ($gC m^{-2} year^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: gpp_week                  !! "weekly" GPP @tex ($gC m^{-2} day^{-1}$) 
                                                                                        !! @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: gdd_m5_dormance           !! growing degree days above threshold -5 
                                                                                        !! deg. C (C) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: gdd_midwinter             !! growing degree days since midwinter (C)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: ncd_dormance              !! number of chilling days since leaves 
                                                                                        !! were lost (days) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: ngd_minus5                !! number of growing days (days)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: time_hum_min              !! time elapsed since strongest moisture 
                                                                                        !! availability (days) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: hum_min_dormance          !! minimum moisture during dormance 
                                                                                        !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts,2), INTENT(inout)          :: gdd_init_date             !! inital date for gdd count
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: gdd_from_growthinit       !! growing degree days, threshold 0 deg. C 
                                                                                        !! since beginning of season
    !pss:+
    REAL(r_std), DIMENSION(npts), INTENT(inout)            :: tsurf_year                !! "annual" surface temperature (K)
    !pss-
!gmjc
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)           :: sla_calc
    ! "14days" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(inout)             :: t2m_14
!end gmjc
!!!qcj++ peatland
    REAL(r_std), DIMENSION(npts,nstm),INTENT(in)                 :: fpeat_map
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)         :: npp0_cumul
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)            :: wtp_daily
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)         :: wtp_month
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)         :: wtp_year

    !
    !! 0.4 Local variables
    !
    INTEGER(i_std)                                          :: j,jst                        !! indices (unitless)
    REAL(r_std)                                             :: ncd_max                  !! maximum ncd (to avoid floating point 
                                                                                        !! underflows) (days) 
    REAL(r_std), DIMENSION(npts)                            :: sumfpc_nat               !! sum of natural fpcs (0-1, unitless)
                                                                                        !! [DISPENSABLE] 
    REAL(r_std), DIMENSION(npts,nvm)                        :: weighttot                !! weight of biomass @tex ($gC m^{-2}$) 
                                                                                        !! @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                        :: nlflong_nat              !! natural long-term leaf NPP 
                                                                                        !! @tex ($gC m^{-2} year^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm)                        :: green_age                !! residence time of green tissue (years)
    REAL(r_std), DIMENSION(npts)                            :: consumption              !! herbivore consumption 
                                                                                        !! @tex ($gC m^{-2} day^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts)                            :: fracnat                  !! fraction of each gridcell occupied by 
                                                                                        !! natural vegetation (0-1, unitless)
    REAL(r_std), DIMENSION(npts)                            :: solad                    !! maximal radation during current day 
                                                                                        !! (clear sky condition)
    REAL(r_std), DIMENSION(npts)                            :: solai                    !! maximal radation during current day
                                                                                        !! (clear sky condition)
    REAL(r_std), DIMENSION(npts)                            :: cloud                    !! cloud fraction

!_ =================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering season'

    !
    !! 1. Initializations
    !! 1.1 Calculate ::ncd_max - the maximum possible NCD (number of chilling days) as:
    !!     \latexonly
    !!     \input{season_ncdmax_eqn2.tex}
    !!     \endlatexonly
    !!     \n
    !!     where one_year is 1 year in seconds (defined in ::constantes).
    !
    ncd_max = ncd_max_year * one_year

    IF ( firstcall_season ) THEN

       !
       !! 1.2 first call - output message giving the setting of the ::gppfrac_dormance,
       !!     ::hvc1, ::hvc2 and ::leaf_frac (as a percentage) parameters which are 
       !!     defined at the beginning of this subroutine. Also outputs the value of 
       !!     ::ncd_max.
       !

       IF ( printlev>=2 ) THEN

          WRITE(numout,*) 'season: '

          WRITE(numout,*) '   > maximum GPP/GGP_max ratio for dormance (::gppfrac_dormance) :',gppfrac_dormance !*

          WRITE(numout,*) '   > maximum possible ncd(days) (::ncd_max) : ',ncd_max !*

          WRITE(numout,*) '   > herbivore consumption C (gC/m2/d) as a function of NPP (gC/m2/d): (::hvc1) (::hvc2)' !*
          WRITE(numout,*) '     C=',hvc1,' * NPP^',hvc2
          WRITE(numout,*) '   > for herbivores, suppose that (::leaf_frac_hvc)',leaf_frac_hvc*100., &
               '% of NPP is allocated to leaves'

       ENDIF

       !
       !! 1.3 Check whether longer-term meteorological and phenology-related variables are not initialized
       !!     to less than zero. If the following meteorological variables are less than ::min_stomate which is set 
       !!     to 1.E-8 in ::grid, then they are set to the current daily value. If the phenology-
       !!     related variables are less than ::min_stomate, they are set to "undef". 
       !!     Warning messages are output in each case, and are also output if "long-term" C fluxes GPP,
       !!     NPP and turnover are less than ::min_stomate.
       !

       !! 1.3.1 moisture availabilities (i.e. relative soil moisture, including the root soil profile):

       !! 1.3.1.1 "monthly" (::moiavail_month)
       IF ( ABS( SUM( moiavail_month(:,2:nvm) ) ) .LT. min_stomate ) THEN

          ! in this case, set them it today's moisture availability
          WRITE(numout,*) 'Warning! We have to initialize the ''monthly'' moisture availabilities. '
          moiavail_month(:,:) = moiavail_daily(:,:)

       ENDIF

       !! 1.3.1.2 "weekly" (::moiavail_week)

       IF ( ABS( SUM( moiavail_week(:,2:nvm) ) ) .LT. min_stomate ) THEN

          ! in this case, set them it today's moisture availability
          WRITE(numout,*) 'Warning! We have to initialize the ''weekly'' moisture availabilities. '
          moiavail_week(:,:) = moiavail_daily(:,:)

       ENDIF

       !! 1.3.2 2-meter temperatures:

       !!  "monthly" (::t2m_month)
       IF ( ABS( SUM( t2m_month(:) ) ) .LT. min_stomate ) THEN

          ! in this case, set them to today's temperature
          WRITE(numout,*) 'Warning! We have to initialize the ''monthly'' 2m temperatures.'
          t2m_month(:) = t2m_daily(:)

       ENDIF

       ! DZ modify: "seasonal"

       IF ( ABS( SUM( Tseason(:) ) ) .LT. min_stomate ) THEN

          ! in this case, set them to today's temperature
          WRITE(numout,*) 'Warning! We have to initialize the ''monthly'' 2m temperatures.'
          Tseason(:) = t2m_daily(:)

       ENDIF

       !!  "weekly" (::2m_week)
       IF ( ABS( SUM( t2m_week(:) ) ) .LT. min_stomate ) THEN

          ! in this case, set them to today's temperature
          WRITE(numout,*) 'Warning! We have to initialize the ''weekly'' 2m temperatures.'
          t2m_week(:) = t2m_daily(:)

       ENDIF

!qcj++ peatland
       IF ( ABS( SUM( wtp_month(:,:) ) ) .LT. min_stomate ) THEN
          ! in this case, set them to today's temperature
          WRITE(numout,*) 'Warning! We have to initialize the ''monthly'' water table depth.'
          wtp_month(:,:) = wtp_daily(:,:)
       ENDIF

       IF ( ABS( SUM( wtp_year(:,:) ) ) .LT. min_stomate ) THEN
          ! in this case, set them to today's temperature
          WRITE(numout,*) 'Warning! We have to initialize the ''yearly'' water table depth.'
          wtp_year(:,:) = wtp_daily(:,:)
       ENDIF

!gmjc
       !! 1.3.2.2 "14days" (::t2m_14)

       IF ( ABS( SUM( t2m_14(:) ) ) .LT. min_stomate ) THEN

          ! in this case, set them to today's temperature
          WRITE(numout,*) 'Warning! We have to initialize the ''14days'' 2m temperatures.'
          t2m_14(:) = t2m_daily(:)

       ENDIF
!end gmjc

       !pss:+ initial tsurf_year
       IF ( ABS( SUM( tsurf_year(:) ) ) .LT. min_stomate ) THEN
          ! in this case, set them to today's temperature
          WRITE(numout,*) 'Warning! We have to initialize the ''annual'' surface temperature.'
          tsurf_year(:) = tsurf_daily(:)
       ENDIF
       !pss:-

       !! 1.3.3 "monthly" soil temperatures (::tsoil_month):
       !MM PAS PARALLELISE!!
       IF ( ABS( SUM( tsoil_month(:,:) ) ) .LT. min_stomate ) THEN

          ! in this case, set them to today's temperature
          WRITE(numout,*) 'Warning!'// &
               ' We have to initialize the ''monthly'' soil temperatures.'
          tsoil_month(:,:) = tsoil_daily(:,:)

       ENDIF

       !! 1.3.4 "monthly" soil humidity (::soilhum_month):
       IF ( ABS( SUM( soilhum_month(:,:) ) ) .LT. min_stomate ) THEN

          ! in this case, set them to today's humidity
          WRITE(numout,*) 'Warning!'// &
               ' We have to initialize the ''monthly'' soil humidity.'
          soilhum_month(:,:) = soilhum_daily(:,:)

       ENDIF

       !! 1.3.5 growing degree days, threshold -5 deg C (::gdd_m5_dormance):
       IF ( ABS( SUM( gdd_m5_dormance(:,2:nvm) ) ) .LT. min_stomate ) THEN
          WRITE(numout,*) 'Warning! Growing degree days (-5 deg) are initialized to ''undef''.'
          gdd_m5_dormance(:,:) = undef
       ENDIF

       !! 1.3.6 growing degree days since midwinter (::gdd_midwinter):
       IF ( ABS( SUM( gdd_midwinter(:,2:nvm) ) ) .LT. min_stomate ) THEN
          WRITE(numout,*) 'Warning! Growing degree days since midwinter' // &
               ' are initialized to ''undef''.'
          gdd_midwinter(:,:) = undef
       ENDIF

       !! 1.3.7 number of chilling days since leaves were lost (::ncd_dormance):
       IF ( ABS( SUM( ncd_dormance(:,2:nvm) ) ) .LT. min_stomate ) THEN
          WRITE(numout,*) 'Warning! Number of chilling days is initialized to ''undef''.'
          ncd_dormance(:,:) = undef
       ENDIF

       !! 1.3.8 number of growing days, threshold -5 deg C (::ngd_minus5):
       IF ( ABS( SUM( ngd_minus5(:,2:nvm) ) ) .LT. min_stomate ) THEN
          WRITE(numout,*) 'Warning! Number of growing days (-5 deg) is initialized to 0.'
       ENDIF

       !! 1.3.9 "long term" npp (::npp_longterm):
       IF ( ABS( SUM( npp_longterm(:,2:nvm) ) ) .LT. min_stomate ) THEN
          WRITE(numout,*) 'Warning! Long term NPP is initialized to 0.'
       ENDIF

       !! 1.3.10 "long term" turnover (::turnover_longterm):
       IF ( ABS( SUM( turnover_longterm(:,2:nvm,:,:) ) ) .LT. min_stomate ) THEN
          WRITE(numout,*) 'Warning! Long term turnover is initialized to 0.'
       ENDIF

       !! 1.3.11 "weekly" GPP (::gpp_week):
       IF ( ABS( SUM( gpp_week(:,2:nvm) ) ) .LT. min_stomate ) THEN
          WRITE(numout,*) 'Warning! Weekly GPP is initialized to 0.'
       ENDIF

       !! 1.3.12 minimum moisture availabilities (::minmoiavail_thisyear)
       !!        In this case, if the minimum moisture
       !!        is less than ::min_stomate, the variable is set to ::large_value,
       !!        which is defined in ::stomate_constants as 1.E33_r_std:
       IF ( ABS( SUM( minmoiavail_thisyear(:,2:nvm) ) ) .LT. min_stomate ) THEN

          ! in this case, set them to a very high value
          WRITE(numout,*) 'Warning! We have to initialize this year''s minimum '// &
               'moisture availabilities.'
          minmoiavail_thisyear(:,:) = large_value

       ENDIF

!!!qcj++ peatland
        !! 1.3.13 cumulative, number of days with npp<=0
        IF ( ABS( SUM( npp0_cumul(:,2:nvm) ) ) .LT. min_stomate ) THEN
           WRITE(numout,*) 'Warning! The cumulative days with null/negative NPP is initialized to 0.'
           npp0_cumul(:,:)=zero
        ENDIF

       !
       !! 1.4 reset firstcall_season flag
       !

       firstcall_season = .FALSE.

    ENDIF

    ! determine min yearly clear sky radiation (solstice) as beginning for gdd count
    cloud(:) = zero
    CALL downward_solar_flux (npts, lalo(:,1),julian_diff,12.,cloud,1,solad,solai)

    WHERE (solad(:) .LT. gdd_init_date(:,2))
       gdd_init_date(:,1)= julian_diff
       gdd_init_date(:,2)= solad(:)
    ENDWHERE


    !
    !! NOTE: Sections 2. to 5. compute slowly-varying, "long-term" (i.e. weekly/monthly)
    !! input variables using a linear relaxation method, following the equation:
    !! \latexonly
    !! \input{season_lin_relax_eqn1.tex}
    !! \endlatexonly
    !! \n 
    !! as described in the introduction to this subroutine (see Krinner et al., 2005).  
    !! The time constant in the above equation is given for each of the variables
    !! described.
    !

    !
    !! 2. Moisture availability (relative soil moisture, including the root profile) :
    !!    The time constants (as in the above equation) for these calculations are 
    !!    given by the parameters ::tau_hum_month (for the monthly
    !!    calculation) or ::tau_hum_week (for the weekly),
    !!    which are set in ::stomate_data to be 20 and 7 days respectively.
    !!    If the moisture availability is less than zero, it is set to zero. 
    !!    These variables are mostly used in the phenological leaf onset and senescence
    !!    models in the modules ::stomate_phenology and ::stomate_turnover, 
    !!    respectively. They are used in models which require a moisture limitation
    !!    condition. These are the 'hum', 'moi', 'humgdd' and 'moigdd' onset models,
    !!    and the 'mixed' and 'dry' (weekly only) senescence models. In addition,
    !!    the weekly moisture availability is used to calculate the limitation
    !!    on the fraction of NPP allocated to different compartments in the module
    !!    ::stomate_alloc. 
    !

    !
    !! 2.1 "monthly" (::moiavail_month)
    !

    moiavail_month = ( moiavail_month * ( tau_hum_month - dt ) + &
         moiavail_daily * dt ) / tau_hum_month

    DO j = 2,nvm  ! Loop over # PFTs
       WHERE ( ABS(moiavail_month(:,j)) .LT. EPSILON(zero) )
          moiavail_month(:,j) = zero
       ENDWHERE
    ENDDO

    !
    !! 2.2 "weekly" (::moiavail_week)
    !

    moiavail_week = ( moiavail_week * ( tau_hum_week - dt ) + &
         moiavail_daily * dt ) / tau_hum_week

    DO j = 2,nvm ! Loop over # PFTs
       WHERE ( ABS(moiavail_week(:,j)) .LT. EPSILON(zero) ) 
          moiavail_week(:,j) = zero
       ENDWHERE
    ENDDO

    !
    !! 3. 2-meter temperatures
    !!    The time constants for the "long-term", "monthly" and "weekly" 2-meter
    !!    temperatures are given by the parameters ::tau_longterm_max, 
    !!    ::tau_t2m_month, and ::tau_t2m_week,
    !!    which are set in ::stomate_data to be 3 * one year (in seconds, as 
    !!    described above) and 20 and 7 days respectively.
    !!    If the temperature is less than zero, it is set to zero.
    !!    In addition the "long term" temperature is written to the history file. \n
    !!    These variables are used in many different modules of STOMATE. 
    !!    The longterm t2m is limied to the intervall [tlong_ref_min, tlong_ref_max]
    !!    using the equation:
    !!      \latexonly
    !!      \input{season_t2m_long_ref_eqn3.tex}
    !!      \endlatexonly
    !!      \n
    !!    The "monthly" and "weekly" temperature is used in the modules 
    !!    ::stomate_phenology and ::stomate_turnover for the onset and senescence 
    !!    phenological models which require a temperature condition. In addition
    !!    the "monthly" temperature is used in ::lpj_constraints to determine 
    !!    the presence and regeneration of the vegetation and in the module
    !!    ::stomate_assimtemp to calculate photosynthesis temperatures.
        
    ! Update tau_longterm
    tau_longterm = MIN(tau_longterm+dt,tau_longterm_max)
    
    ! Recalculate the reference temperature using the old reference temperature and the current temperature
    t2m_longterm(:) = ( t2m_longterm(:) * ( tau_longterm - dt ) + &
         t2m_daily(:) * dt ) / tau_longterm

    ! The longterm reference is not allowed to go outside the interval [tlong_ref_min, tlong_ref_max]
    t2m_longterm(:) = MAX( tlong_ref_min, MIN( tlong_ref_max, t2m_longterm(:) ) )

    CALL xios_orchidee_send_field("t2m_longterm",t2m_longterm)

    CALL histwrite_p (hist_id_stomate, 'T2M_LONGTERM', itime, &
         t2m_longterm, npts, hori_index)

    !
    !! 3.3 "monthly" (::t2m_month)
    !

    t2m_month = ( t2m_month * ( tau_t2m_month - dt ) + &
         t2m_daily * dt ) / tau_t2m_month

    WHERE ( ABS(t2m_month(:)) .LT. EPSILON(zero) )
       t2m_month(:) = zero
    ENDWHERE

    ! Store the day in onset_date when the leaves started grow. julian_diff is the day in the year [1-366].
    WHERE ( begin_leaves(:,:) )
       onset_date(:,:)=julian_diff
    ELSEWHERE
       onset_date(:,:)=0
    ENDWHERE
    
    ! Calculate "seasonal" temperature
    WHERE ( t2m_week (:) .GT. ZeroCelsius )
       Tseason_tmp(:) = Tseason_tmp(:) + t2m_week(:)
       Tseason_length(:)=Tseason_length(:) + dt
    ENDWHERE

    ! calculate Tmin in spring (after onset)
    DO j=1, nvm
       IF (leaf_tab(j)==1 .AND. pheno_type(j)==2) THEN
          ! leaf_tab=broadleaf and pheno_typ=summergreen
          ! Only treat broadleag and summergreen pfts
          
          WHERE ( Tmin_spring_time(:,j)>0 .AND. (Tmin_spring_time(:,j)<spring_days_max) )
             Tmin_spring_time(:,j)=Tmin_spring_time(:,j)+1
          ELSEWHERE
             Tmin_spring_time(:,j)=0
          ENDWHERE
          
          WHERE ( begin_leaves(:,j) )
             Tmin_spring_time(:,j)=1
          ENDWHERE
       END IF
    END DO

    !
    ! Store the day in onset_date when the leaves started grow. julian_diff is the day in the year [1-366].
    WHERE ( begin_leaves(:,:) )
       onset_date(:,:)=julian_diff
    ELSEWHERE 
       onset_date(:,:)=0
    ENDWHERE
    !
    !
    ! DZ modify: calculate "seasonal" temperature
    !

!qcj debug: duplicate calculation, double counting 
   ! WHERE ( t2m_week (:) .GT. ZeroCelsius )
   !    Tseason_tmp(:) = Tseason_tmp(:) + t2m_week(:)
   !    Tseason_length(:)=Tseason_length(:) + dt
   ! ENDWHERE

    WHERE ( ABS(Tseason_tmp(:)) .LT. EPSILON(zero) )
       Tseason_tmp(:) = zero
    ENDWHERE
    WHERE ( ABS(Tseason_length(:) ) .LT. EPSILON(zero) )
       Tseason_length(:) = zero
    ENDWHERE
    
    !
    ! calculate Tmin in spring (after onset)
    !
!qcj debug: duplicate calculation, double counting 
   ! DO j = 2,nvm
   !   IF (pft_to_mtc(j) == 8 .OR. pft_to_mtc(j) == 6) THEN
   !     WHERE ( Tmin_spring_time(:,j)>0 .AND. (Tmin_spring_time(:,j)<40) )
   !        Tmin_spring_time(:,j)=Tmin_spring_time(:,j)+1
   !     ELSEWHERE 
   !        Tmin_spring_time(:,j)=0
   !     ENDWHERE

   !     WHERE ( begin_leaves(:,j) )
   !        Tmin_spring_time(:,j)=1
   !     ENDWHERE

   !   ENDIF
   ! ENDDO

    !
    !! 3.4 "weekly" (::t2m_week)
    !

    t2m_week = ( t2m_week * ( tau_t2m_week - dt ) + &
         t2m_daily * dt ) / tau_t2m_week

    WHERE ( ABS(t2m_week(:)) .LT. EPSILON(zero) )
       t2m_week(:) = zero
    ENDWHERE
!gmjc
    !
    ! 3.5 "14days"
    !

    t2m_14 = ( t2m_14 * ( tau_t2m_14 - dt ) + &
                 t2m_daily * dt ) / tau_t2m_14

    WHERE ( ABS(t2m_14(:)) .LT. EPSILON(zero) )
       t2m_14(:) = t2m_daily(:)
    ENDWHERE
!end gmjc
    !
    !! 4. "monthly" soil temperatures (::tsoil_month)
    !!    The time constant is given by the parameter ::tau_tsoil_month, 
    !!    which is set in ::stomate_data to be 20 days.
    !!    If the monthly soil temperature is less than zero, it is set to zero. \n
    !!    This variable is used in the ::stomate_allocation module.
    !

    tsoil_month = ( tsoil_month * ( tau_tsoil_month - dt ) + &
         tsoil_daily(:,:) * dt ) / tau_tsoil_month

    WHERE ( ABS(tsoil_month(:,:)) .LT. EPSILON(zero) )
       tsoil_month(:,:) = zero
    ENDWHERE

    !
    !! 5. ''monthly'' soil humidity (or relative soil moisture - ::soilhum_month) 
    !!    The time constant is given by the parameter ::tau_soilhum_month, 
    !!    which is set in ::stomate_data to be 20 days.
    !!    If the monthly soil humidity is less than zero, it is set to zero.
    !

    soilhum_month = ( soilhum_month * ( tau_soilhum_month - dt ) + &
         soilhum_daily * dt ) / tau_soilhum_month

    WHERE ( ABS(soilhum_month(:,:)) .LT. EPSILON(zero) )
       soilhum_month(:,:) = zero
    ENDWHERE

!qcj++ peatland !'monthly' water table position
    DO j=1,nvm
      wtp_month(:,j) = ( wtp_month(:,j) * ( tau_soilhum_month - dt ) + &
         wtp_daily(:,j) * dt ) / tau_soilhum_month

      WHERE ( ABS(wtp_month(:,j)) .LT. EPSILON(zero) )
         wtp_month(:,j) = zero
      ENDWHERE
    ENDDO

    DO j=1,nvm
      wtp_year(:,j) = ( wtp_year(:,j) * ( tau_wtp_year - dt ) + &
         wtp_daily(:,j) * dt ) / tau_wtp_year

      WHERE ( ABS(wtp_year(:,j)) .LT. EPSILON(zero) )
         wtp_year(:,j) = zero
      ENDWHERE
    ENDDO

    CALL xios_orchidee_send_field("WTP_YEAR",wtp_year)

    !
    !! 6. Calculate the dormance time-length (::time_lowgpp).
    !!    The dormancy time-length is increased by the stomate time step when   
    !!    the weekly GPP is below one of two thresholds, 
    !!    OR if the growing season started more than 2 years ago AND the amount of biomass
    !!    in the carbohydrate reserve is 4 times as much as the leaf biomass.
    !!    i.e. the plant is declared dormant if it is accumulating carbohydrates but not 
    !!    using them, so that the beginning of the growing season can be detected .
    !     This last condition was added by Nicolas Viovy.
    !!    - NV: special case (3rd condition)
    !!    Otherwise, set ::time_lowgpp to zero. \n
    !!    The weekly GPP must be below either the parameter ::min_gpp_allowed, which is set at 
    !!    the beginning of this subroutine to be 0.3gC/m^2/year, 
    !!    OR it must be less than last year's maximum weekly GPP multiplied by the
    !!    maximal ratio of GPP to maximum GPP (::gppfrac_dormance), which is set to 0.2 at the 
    !!    beginning of this subroutine. \n
    !!    This variable is used in the ::stomate_phenology module to allow the detection of the
    !!    beginning of the vegetation growing season for different onset models.
    !!    Each PFT is looped over.
    !
    !NVMODIF
!!$    DO j = 2,nvm ! Loop over # PFTs
!!$       WHERE ( ( gpp_week(:,j) .LT. min_gpp_allowed ) .OR. & 
!!$            ( gpp_week(:,j) .LT. gppfrac_dormance * maxgppweek_lastyear(:,j) ) .OR. &
!!$            ( ( when_growthinit(:,j) .GT. 2.*one_year ) .AND. &
!!$            ( biomass(:,j,icarbres,icarbon) .GT. biomass(:,j,ileaf,icarbon)*4. ) ) )
!!$          !       WHERE ( ( gpp_week(:,j) .EQ. zero ) .OR. & 
!!$          !            ( gpp_week(:,j) .LT. gppfrac_dormance * maxgppweek_lastyear(:,j) ) .OR. &
!!$          !            ( ( when_growthinit(:,j) .GT. 2.*one_year ) .AND. &
!!$          !            ( biomass(:,j,icarbres,icarbon) .GT. biomass(:,j,ileaf,icarbon)*4. ) ) )
!!$       
!!$          time_lowgpp(:,j) = time_lowgpp(:,j) + dt
!!$          
!!$       ELSEWHERE
!!$          
!!$          time_lowgpp(:,j) = zero
!!$
!!$       ENDWHERE
!!$    ENDDO

    !
    !! 7. Calculate the growing degree days (GDD - ::gdd_m5_dormance).
    !!    This variable is the GDD sum of the temperatures higher than -5 degrees C.
    !!    It is inititalised to 0 at the beginning of the dormancy period (i.e.
    !!    when ::time_lowgpp>0), and is set to "undef" when there is GPP 
    !!    (::time_lowgpp=0), as it is not used during the growing season. \n 
    !!    ::gdd_m5_dormance is further scaled by (::tau_gdd -dt / ::tau_gdd), 
    !!    - ::tau_gdd is set to 40. in the module ::stomate_constants. \n
    !     Nicolas Viovy - not sure why this is...
    !!    This variable is used in the ::stomate_phenology module for the 
    !!    'humgdd' and 'moigdd' leaf onset models. \n
    !!    Each PFT is looped over but ::gdd_m5_dormance is only calculated for 
    !!    those PFTs for which a critical GDD is defined, i.e. those which are 
    !!    assigned to the 'humgdd' or 'moigdd' models. \n
    !!    Finally if GDD sum is less than zero, then it is set to zero. 
    !

    DO j = 2,nvm ! Loop over # PFTs

       IF (.NOT. natural(j)) THEN
          ! reset counter: start of the growing season
          WHERE ((when_growthinit(:,j) .EQ. zero))
	     gdd_from_growthinit(:,j) = zero
	  ENDWHERE
          ! increase gdd counter
          WHERE ( t2m_daily(:) .GT. (ZeroCelsius)  )
             gdd_from_growthinit(:,j) = gdd_from_growthinit(:,j) + &
                                 dt * ( t2m_daily(:) - (ZeroCelsius) )
          ENDWHERE
       ELSE
          gdd_from_growthinit(:,j) = undef
       ENDIF

       ! only for PFTs for which critical gdd is defined
       ! gdd_m5_dormance is set to 0 at the end of the growing season. It is set to undef
       ! at the beginning of the growing season.

       IF ( ALL(pheno_gdd_crit(j,:) .NE. undef) ) THEN

          !
          !! 7.1 set to zero if undefined and there is no GPP
          !
          IF (printlev>=4) THEN
            WRITE(numout,*) 'j', j
            WRITE(numout,*) 'gdd_init_date(:,1)', gdd_init_date(:,1)
            WRITE(numout,*) 'julian_diff',julian_diff
            WRITE(numout,*) 'gdd_m5_dormance(:,j)', gdd_m5_dormance(:,j)
          ENDIF

          WHERE (gdd_init_date(:,1) .EQ. julian_diff)

             gdd_m5_dormance(:,j) = zero

          ENDWHERE

          !
          !! 7.2 set to undef if there is GPP
          !

          WHERE ( when_growthinit(:,j) .EQ. zero )

             gdd_m5_dormance(:,j) = undef

          ENDWHERE

          !
          !! 7.3 normal update as described above where ::gdd_m5_dormance is defined (not set to "undef").
          !

          WHERE ( ( t2m_daily(:) .GT. ( ZeroCelsius - gdd_threshold ) ) .AND. &
               ( gdd_m5_dormance(:,j) .NE. undef )           )
             gdd_m5_dormance(:,j) = gdd_m5_dormance(:,j) + &
                  dt * ( t2m_daily(:) - ( ZeroCelsius - gdd_threshold ) )
          ENDWHERE

          WHERE ( gdd_m5_dormance(:,j) .NE. undef ) 
             gdd_m5_dormance(:,j) = gdd_m5_dormance(:,j) * &
                  ( tau_gdd - dt ) / tau_gdd
          ENDWHERE
          IF (printlev>=4) THEN 
            WRITE(numout,*) 'j',j
            WRITE(numout,*) 'when_growthinit(:,j)',when_growthinit(:,j)
            WRITE(numout,*) 'gdd_m5_dormance(:,j)',gdd_m5_dormance(:,j)
          ENDIF

       ENDIF

    ENDDO

    !
    !! 7.4 Set to zero if GDD is less than zero.

    DO j = 2,nvm ! Loop over # PFTs
       WHERE ( ABS(gdd_m5_dormance(:,j)) .LT. EPSILON(zero) )
          gdd_m5_dormance(:,j) = zero
       ENDWHERE
    ENDDO

    !
    !! 8. Calculate the growing degree days (GDD) since midwinter (::gdd_midwinter)
    !!    This variable represents the GDD sum of temperatures higher than a PFT-dependent
    !!    threshold (::ncdgdd_temp), since midwinter.
    !!    Midwinter is detected if the monthly temperature (::t2m_month) is lower than the weekly 
    !!    temperature (::t2m_week) AND if the monthly temperature is lower than the long-term 
    !!    temperature (t2m_longterm). These variables were calculated earlier in this subroutine. \n
    !!    ::gdd_midwinter is initialised to 0.0 when midwinter is detected, and from then on 
    !!    increased with each temperature greater than ::ncdgdd_temp, which is defined
    !!    in the module ::stomate_constants. \n
    !!    ::gdd_midwinter is set to "undef" when midsummer is detected, following the opposite
    !!    conditions to those used to define midwinter. \n
    !!    The variable is used in the ::stomate_phenology module for the leaf onset model 'ncdgdd'.
    !!    Each PFT is looped over but the ::gdd_midwinter is only calculated for those
    !!    PFTs for which a critical 'ncdgdd' temperature is defined, i.e. those which are 
    !!    assigned to the 'ncdgdd' model. \n
    !

    DO j = 2,nvm ! Loop over # PFTs

       ! only for PFTs for which ncdgdd_crittemp is defined

       IF ( ncdgdd_temp(j) .NE. undef ) THEN

          !
          !! 8.1 set to 0 if undef and if we detect "midwinter"
          !

!!$          WHERE ( ( gdd_midwinter(:,j) .EQ. undef ) .AND. &
!!$               ( t2m_month(:) .LT. t2m_week(:) ) .AND. &
!!$               ( t2m_month(:) .LT. t2m_longterm(:) )    )
          WHERE (gdd_init_date(:,1) .EQ. julian_diff)

             gdd_midwinter(:,j) = zero

          ENDWHERE

          !
          !! 8.2 set to undef if we detect "midsummer"
          !

          WHERE ( ( t2m_month(:) .GT. t2m_week(:) ) .AND. &
               ( t2m_month(:) .GT. t2m_longterm(:) )    )

             gdd_midwinter(:,j) = undef

          ENDWHERE

          !
          !! 8.3 normal update as described above
          !

          WHERE ( ( gdd_midwinter(:,j) .NE. undef ) .AND. &
               ( t2m_daily(:) .GT. ncdgdd_temp(j)+ZeroCelsius ) )

             gdd_midwinter(:,j) = &
                  gdd_midwinter(:,j) + &
                  dt * ( t2m_daily(:) - ( ncdgdd_temp(j)+ZeroCelsius ) )

          ENDWHERE

       ENDIF

    ENDDO

    !
    !! 9. Calculate the number of chilling days (NCD) since leaves were lost (::ncd_dormance).
    !!    This variable is initialised to 0 at the beginning of the dormancy period (::time_lowgpp>0)
    !!    and increased by the stomate time step when the daily temperature is lower than the 
    !!    PFT-dependent threshold ::ncdgdd_temp, which is defined for each PFT 
    !!    in a table (::ncdgdd_temp_tab) in the module ::stomate_data. \n
    !!    It is set to "undef" when there is GPP (::time_lowgpp=0) as it is not needed during
    !!    the growing season. \n
    !!    The variable is used in the ::stomate_phenology module for the leaf onset model 'ncdgdd'.
    !!    Each PFT is looped over but the ::ncd_dormance is only calculated for those
    !!    PFTs for which a critical 'ncdgdd' temperature is defined, i.e. those which are 
    !!    assigned to the 'ncdgdd' model.
    !

    DO j = 2,nvm ! Loop over # PFTs

       IF ( ncdgdd_temp(j) .NE. undef ) THEN

          !
          !! 9.1 set to zero if undefined and there is no GPP
          !

          WHERE (gdd_init_date(:,1) .EQ. julian_diff)

             ncd_dormance(:,j) = zero

          ENDWHERE

          !
          !! 9.2 set to undef if there is GPP
          !

          WHERE ( when_growthinit(:,j) .EQ. zero )

             ncd_dormance(:,j) = undef

          ENDWHERE

          !
          !! 9.3 normal update, as described above, where ::ncd_dormance is defined
          !

          WHERE ( ( ncd_dormance(:,j) .NE. undef ) .AND. &
               ( t2m_daily(:) .LE. ncdgdd_temp(j)+ZeroCelsius ) )

             ncd_dormance(:,j) = MIN( ncd_dormance(:,j) + dt, ncd_max )

          ENDWHERE

       ENDIF

    ENDDO

    !
    !! 10. Calculate the number of growing days (NGD) since leaves were lost (::ngd_minus5).
    !!     This variable is initialised to 0 at the beginning of the dormancy period (::time_lowgpp>0)
    !!     and increased by the stomate time step when the daily temperature is higher than the threshold
    !!     -5 degrees C. \n    
    !!     ::ngd_minus5 is further scaled by (::tau_ngd -dt / ::tau_ngd), 
    !!     - ::tau_ngd is set to 50. in the module ::stomate_constants. \n
    !      Nicolas Viovy - not sure why this is...
    !!     The variable is used in the ::stomate_phenology module for the leaf onset model 'ngd'.
    !!     Each PFT is looped over.
    !!     If the NGD is less than zero, it is set to zero.
    !

    DO j = 2,nvm ! Loop over # PFTs

       !
       !! 10.1 Where there is GPP (i.e. ::time_lowgpp=0), set NGD to 0. 
       !!      This means that we only take into account NGDs when the leaves are off
       !

       WHERE (gdd_init_date(:,1) .EQ. julian_diff)
          ngd_minus5(:,j) = zero
       ENDWHERE

       !
       !! 10.2 normal update, as described above.
       !

       WHERE ( t2m_daily(:) .GT. (ZeroCelsius - gdd_threshold) )
          ngd_minus5(:,j) = ngd_minus5(:,j) + dt
       ENDWHERE

       ngd_minus5(:,j) = ngd_minus5(:,j) * ( tau_ngd - dt ) / tau_ngd

    ENDDO

    DO j = 2,nvm ! Loop over # PFTs
       WHERE ( ABS(ngd_minus5(:,j)) .LT. EPSILON(zero) )
          ngd_minus5(:,j) = zero
       ENDWHERE
    ENDDO

    !
    !! 11. Calculate the minimum humidity/relative soil moisture since dormance began (::hum_min_dormance) 
    !!     and the time elapsed since this minimum (::time_hum_min). \n
    !!     The minimum moisture availability occurs when the monthly moisture availability, which is updated 
    !!     daily earlier in this subroutine as previously described, is at a minimum during the dormancy period.
    !!     Therefore the ::hum_min_dormance is initialised to the monthly moisture availability 
    !!     at the beginning of the dormancy period (i.e. if it was previously set to undefined and the
    !!     ::time_lowgpp>0) AND then whenever the monthly moisture availability is less than it has previously
    !!     been. \n
    !!     Consequently, the time counter (::time_hum_min) is initialised to 0 at the beginning of the dormancy  
    !!     period (i.e. if it was previously set to undefined and the ::time_lowgpp>0) AND when the minimum 
    !!     moisture availability is reached, and is increased by the stomate time step every day throughout
    !!     the dormancy period. \n
    !!     ::time_hum_min is used in the ::stomate_phenology module for the leaf onset models 'moi' and 'moigdd'.
    !!     Each PFT is looped over but the two variables are only calculated for those
    !!     PFTs for which the critical parameter ::hum_min_time is defined, i.e.   
    !!     those which are assigned to the 'moi' or 'moigdd' models.
    !

    DO j = 2,nvm ! Loop over # PFTs

       IF ( hum_min_time(j) .NE. undef ) THEN

          !
          !! 11.1 initialize if undefined and there is no GPP
          !

          WHERE (when_growthinit(:,j) .EQ. zero)

             time_hum_min(:,j) = zero
             hum_min_dormance(:,j) = moiavail_month(:,j)

          ENDWHERE

          !
          !! 11.3 normal update, as described above, where ::time_hum_min and ::hum_min_dormance are defined
          !

          !! 11.3.1 increase time counter by the stomate time step

          WHERE ( hum_min_dormance(:,j) .NE. undef )

             time_hum_min(:,j) = time_hum_min(:,j) + dt

          ENDWHERE

          !! 11.3.2 set time counter to zero if minimum is reached

          WHERE (( hum_min_dormance(:,j) .NE. undef ) .AND. &
               ( moiavail_month(:,j) .LE. hum_min_dormance(:,j) ) )

             hum_min_dormance(:,j) = moiavail_month(:,j)
             time_hum_min(:,j) = zero

          ENDWHERE

       ENDIF

    ENDDO

    !
    !! NOTE: Sections 12. to 14. compute slowly-varying, "long-term" (i.e. weekly/monthly)
    !! C fluxes (NPP, turnover, GPP) using the same linear relaxation method as in 
    !! Sections 2. and 5. and described in the introduction to this section, as given 
    !! by the following equation:
    !! \latexonly
    !! \input{season_lin_relax_eqn1.tex}
    !! \endlatexonly
    !! \n 
    !! The following variables are calculated using the above equation, and the time constant 
    !! is given for each. 
    !

    !
    !! 12. Update the "long term" NPP. (::npp_daily in gC/m^2/day, ::npp_longterm in gC/m^2/year.)
    !!     The time constant is given by the parameter ::tau_longterm, 
    !!     which is set in ::stomate_data to be 3 * one year (in seconds, as 
    !!     described above). If the ::npp_longterm is less than zero then set it to zero. \n
    !!     ::npp_longterm is used in ::stomate_lpj in the calculation of long-term
    !!     vegetation dynamics and in ::stomate_lcchange in the calculation of 
    !!     land cover change. It is also used to calculate diagnose the hebivory activity in 
    !!     Section 22 of this subroutine.
    !

    npp_longterm = ( npp_longterm * ( tau_longterm - dt ) + &
         (npp_daily*one_year) * dt                          ) / &
         tau_longterm

    DO j = 2,nvm ! Loop over # PFTs
       WHERE ( ABS(npp_longterm(:,j)) .LT. EPSILON(zero) )
          npp_longterm(:,j) = zero
       ENDWHERE
    ENDDO

!!!qcj++ peatland
    !!  Following Druel et al. 2017
    !!  Update the counter of number of day where npp is negative or =0
    !!  If this number is too long, in stomate_turnover.f90 the  turnover go increase by stress !
    !!  For the moment if one time npp > 0 one day, npp0_cumul go to be reset... maybe that need more than 1 day
    !!  For the moment it is only for moss/lichen !!

    DO j = 2,nvm ! Loop over # PFTs
       IF ( is_mosspeat(j) ) THEN
          WHERE ( npp_daily(:,j) .LT. min_stomate )
              npp0_cumul(:,j) = npp0_cumul(:,j) + 1
          ELSEWHERE
              npp0_cumul(:,j) = 0
          ENDWHERE
       ENDIF
    ENDDO

    !
    !! 13. Update the "long term" turnover rates (in gC/m^2/year).
    !!     The time constant is given by the parameter ::tau_longterm, 
    !!     which is set in ::stomate_data to be 3 * one year (in seconds, as 
    !!     described above). If the ::turnover_longterm is less than zero then set it to zero.\n
    !!     ::turnover_longterm is used in ::stomate_lpj and :: lpg_gap in the calculation 
    !!     of long-term vegetation dynamics.
    !

    turnover_longterm(:,:,:,:) = ( turnover_longterm(:,:,:,:) * ( tau_longterm - dt ) + &
         (turnover_daily(:,:,:,:)*one_year) * dt                          ) / &
         tau_longterm

    DO j = 2,nvm ! Loop over # PFTs
       WHERE ( ABS(turnover_longterm(:,j,:,:)) .LT. EPSILON(zero) )
          turnover_longterm(:,j,:,:) = zero
       ENDWHERE
    ENDDO

    !
    !! 14. Update the "weekly" GPP (where there is vegetation), otherwise set to zero.
    !!     The time constant is given by the parameter ::tau_gpp_week, 
    !!     which is set in ::stomate_data to be 7 days. If the ::gpp_week is 
    !!     less than zero then set it to zero. \n
    !!     ::gpp_week is used to update the annual maximum weekly GPP (::maxgppweek_thisyear)
    !!     in Section 16 of this subroutine, which is then used to update the variable
    !!     ::maxgppweek_lastyear in Section 21 of this subroutine. Both ::gpp_week and 
    !!     ::maxgppweek_lastyear are used in Section 6 of this subroutine to calculate
    !!     the onset and time-length of the dormancy period. 
    !      Note: Used to be weekly GPP divided by veget_cov_max, i.e. per ground covered, but not anymore.
    !
    WHERE ( veget_cov_max .GT. zero )

       gpp_week = ( gpp_week * ( tau_gpp_week - dt ) + &
            gpp_daily * dt ) / tau_gpp_week

    ELSEWHERE
       gpp_week = zero

    ENDWHERE

    DO j = 2,nvm ! Loop over # PFTs
       WHERE ( ABS(gpp_week(:,j)) .LT. EPSILON(zero) )
          gpp_week(:,j) = zero
       ENDWHERE
    ENDDO

    !
    !! 15. Update the maximum and minimum moisture availabilities (::maxmoiavail_thisyear 
    !!     and ::minmoiavail_thisyear). If the daily moisture availability, ::moiavail_daily,
    !!     which was calculated earlier in this subroutine, is greater or less than the current
    !!     value of the maximum or minimum moisture availability, then set the value for the maximum
    !!     or minimum to the current daily value respectively. \n
    !!     ::maxmoiavail_thisyear and ::minmoiavail_thisyear are used to update the variables 
    !!     ::maxmoiavail_lastyear and ::minmoiavail_lastyear in Section 21 of this subroutine, 
    !!     which are used in the module ::stomate_phenology for the leaf onset models 'hum' and 
    !!     'humgdd', and in the module ::stomate_turnover for the leaf senescence models 'dry' 
    !!     and 'mixed'.
    !

    WHERE ( moiavail_daily .GT. maxmoiavail_thisyear )
       maxmoiavail_thisyear = moiavail_daily
    ENDWHERE

    WHERE ( moiavail_daily .LT. minmoiavail_thisyear )
       minmoiavail_thisyear = moiavail_daily
    ENDWHERE

    !
    !! 16. Update the annual maximum weekly GPP (::maxgppweek_thisyear), if the weekly GPP
    !!     is greater than the current value of the annual maximum weekly GPP.
    !!     The use of this variable is described in Section 14.
    !

    WHERE ( gpp_week .GT. maxgppweek_thisyear )
       maxgppweek_thisyear = gpp_week
    ENDWHERE

    !
    !! 17. Update the annual GDD0 by adding the current daily 2-meter temperature 
    !!     (multiplied by the stomate time step, which is one day), if it is greater
    !!     than zero degrees C. \n
    !!     This variable is mostly used in the module ::lpj_establish.
    !

    WHERE ( t2m_daily .GT. ZeroCelsius )
       gdd0_thisyear = gdd0_thisyear + dt * ( t2m_daily - ZeroCelsius )
    ENDWHERE

    !
    !! 18. Update the annual precipitation by adding the current daily precipitation
    !!     amount (multiplied by the stomate time step, which is usually one day). \n
    !

    precip_thisyear = precip_thisyear + dt * precip_daily

    
!
    !! 19. Update the annual maximum leaf mass for each PFT (::lm_thisyearmax) and the maximum fractional
    !!     plant cover if the LPJ DGVM is activated.       
    !
    
    !
    !! 19.1 If the LPJ DGVM is activated first the fraction of natural vegetation (::fracnat), i.e. non-
    !!      agricultural vegetation (PFTs 2-11) is calculated for each PFT. Each PFT is looped over.
    !

    IF (ok_dgvm ) THEN

       fracnat(:) = un
       DO j = 2,nvm ! Loop over # PFTs
          IF ( .NOT. natural(j) .OR. pasture(j)) THEN
             fracnat(:) = fracnat(:) - veget_cov_max(:,j)
          ENDIF
       ENDDO

    ENDIF

    !
    !! 19.2 If LPJ and STOMATE are activated, first the maximum fractional plant cover needs to be updated, 
    !!      and then this year's leaf biomass. Each PFT is looped over.
    !!      Both are updated according to the linear relaxation method described above. 
    !!      \latexonly
    !!      \input{season_lin_relax_eqn1.tex}
    !!      \endlatexonly
    !!      \n
    !!      The time constant for this process is set as one year (in seconds) divided by the parameter
    !!      ::leaflife_tab, which gives the leaf lifetime in years and is set for each PFT in the module 
    !!      ::stomate_constants. Each PFT is looped over. \n
    !

    IF ( ok_stomate ) THEN
       IF(ok_dgvm ) THEN
          DO j=2,nvm ! Loop over # PFTs

             !
             !! 19.2.1 Calculate maximum fractional plant cover (::maxfpc_lastyear).
             !!        If natural vegetation is present in the grid cell, and the leaf
             !!        biomass is greater than three-quarters of last year's leaf biomass, the maximum fractional plant
             !!        cover for last year is updated. \n 
             !!        The short-term variable (Xs in the above equation) that is being used to update the long-term 
             !!        maximum fractional plant cover is the fractional cover of natural vegetation, specified as 
             !!        ::veget_cov/::fracnat. Last year's value is then set to be this year's value.
             !

             IF ( natural(j) .AND. ok_dgvm .AND. .NOT. pasture(j)) THEN

                WHERE ( fracnat(:) .GT. min_stomate .AND. biomass(:,j,ileaf,icarbon).GT. lm_lastyearmax(:,j)*0.75 )
                   maxfpc_lastyear(:,j) = ( maxfpc_lastyear(:,j) * ( one_year/leaflife_tab(j)- dt ) + &
                        veget_cov(:,j) / fracnat(:) * dt ) / (one_year/leaflife_tab(j))
                ENDWHERE
                maxfpc_thisyear(:,j) = maxfpc_lastyear(:,j) ! just to initialise value
                
             ENDIF

!NV : correct initialization
!!$             WHERE(biomass(:,j,ileaf,icarbon).GT. lm_lastyearmax(:,j)*0.75)
!!$                lm_lastyearmax(:,j) = ( lm_lastyearmax(:,j) * ( one_year/leaflife_tab(j)- dt ) + &
!!$                     biomass(:,j,ileaf,icarbon) * dt ) / (one_year/leaflife_tab(j))
!!$             ENDWHERE
!!$             lm_thisyearmax(:,j)=lm_lastyearmax(:,j) ! just to initialise value

             
             !
             !! 19.2.2 Update this year's leaf biomass (::lm_thisyearmax).
             !!        The short-term variable (Xs in the above equation) that is being used to update the long-term 
             !!        this year's leaf biomass is the leaf biomass pool (::biomass(i,j,ileaf).
             !
             WHERE (lm_thisyearmax(:,j) .GT. min_stomate)
                WHERE(biomass(:,j,ileaf,icarbon).GT. lm_thisyearmax(:,j)*0.75)
                   lm_thisyearmax(:,j) = ( lm_thisyearmax(:,j) * ( one_year/leaflife_tab(j)- dt ) + &
                        biomass(:,j,ileaf,icarbon) * dt ) / (one_year/leaflife_tab(j))
                ENDWHERE
             ELSEWHERE
                lm_thisyearmax(:,j) =biomass(:,j,ileaf,icarbon)
             ENDWHERE

          ENDDO

       ELSE

!!!qcj++ peatland
         IF (ok_dgvm_peat) THEN 

            DO j=2,nvm 
               jst=pref_soil_veg(j)
               IF ( natural(j) .AND. (.NOT. pasture(j)) .AND. is_peat(j) ) THEN

                  WHERE ( fpeat_map(:,jst) .GT. min_stomate .AND. biomass(:,j,ileaf,icarbon).GT. lm_lastyearmax(:,j)*0.75 )            
                     maxfpc_lastyear(:,j) = ( maxfpc_lastyear(:,j) * (one_year/leaflife_tab(j)- dt ) + &
                         veget_cov(:,j) / fpeat_map(:,jst) * dt ) /(one_year/leaflife_tab(j))
                  ENDWHERE
                  maxfpc_thisyear(:,j) = maxfpc_lastyear(:,j)
 
                  WHERE (lm_thisyearmax(:,j) .GT. min_stomate)
                     WHERE(biomass(:,j,ileaf,icarbon).GT. lm_thisyearmax(:,j)*0.75)
                         lm_thisyearmax(:,j) = ( lm_thisyearmax(:,j) * ( one_year/leaflife_tab(j)- dt ) + &
                               biomass(:,j,ileaf,icarbon) * dt ) / (one_year/leaflife_tab(j))
                     ENDWHERE
                  ELSEWHERE
                     lm_thisyearmax(:,j) =biomass(:,j,ileaf,icarbon)
                  ENDWHERE

               ELSE
                  WHERE ( biomass(:,j,ileaf,icarbon) .GT. lm_thisyearmax(:,j) )
                       lm_thisyearmax(:,j) = biomass(:,j,ileaf,icarbon)
                  ENDWHERE
               ENDIF

            ENDDO 

         ELSE          
          !
          !! 19.3 If LPJ DGVM is not activated but STOMATE is, the maximum leaf mass is set to be the same  
          !!      as the leaf biomass (::biomass(i,j,ileaf), without a change in the maximum fractional plant cover. 
          !
          DO j = 2,nvm ! Loop over # PFTs
             WHERE ( biomass(:,j,ileaf,icarbon) .GT. lm_thisyearmax(:,j) )
                lm_thisyearmax(:,j) = biomass(:,j,ileaf,icarbon)
             ENDWHERE
          ENDDO

         ENDIF

       ENDIF !(ok_dgvm)
    ELSE

          !
          !! 19.4 If STOMATE is not activated, the maximum leaf biomass is set to be the maximum possible 
          !!      LAI of the PFT.
          !

       DO j = 2,nvm ! Loop over # PFTs
!JCMODIF
!          lm_thisyearmax(:,j) = lai_max(j)  / sla(j)
          lm_thisyearmax(:,j) = lai_max(j)  / sla_calc(:,j)
!ENDJCMODIF
       ENDDO

    ENDIF  !(ok_stomate)

    !
    !! 20. Update the annual maximum fractional plant cover for each PFT if the current fractional cover 
    !!     (::veget_cov), is larger than the current maximum value.
    !!     ::veget_cov is defined as fraction of total ground. Therefore, maxfpc_thisyear has the same unit.
    !

    WHERE ( veget_cov(:,:) .GT. maxfpc_thisyear(:,:) )
       maxfpc_thisyear(:,:) = veget_cov(:,:)
    ENDWHERE

    !
    !! 21. At the end of the every year, last year's maximum and minimum moisture availability,
    !!     annual GDD0, annual precipitation, annual max weekly GPP, and maximum leaf mass are 
    !!     updated with the value calculated for the current year, and then their value reset.
    !!     Either to zero for the maximum variables, or to ::large_value (set to be 1.E33 in ::
    !!     the module ::stomate_constants) for the minimum variable ::minmoiavail_thisyear.
    !

    IF ( LastTsYear ) THEN

       !
       !! 21.1 Update last year's values. \n
       !!      The variables ::maxmoiavail_lastyear, ::minmoiavail_lastyear and ::maxgppweek_lastyear  
       !!      are updated using the relaxation method:
       !!      \latexonly
       !!      \input{season_lin_relax_eqn1.tex}
       !!      \endlatexonly
       !!      \n
       !!      where Xs is this year's value, dt (Delta-t) is set to 1 (year), and the time constant (tau) !??
       !!      is set by the parameter ::tau_climatology, which is set to 20 years at the beginning 
       !!      of this subroutine. \n
       !!      The other variables (::gdd0_lastyear, ::precip_lastyear, ::lm_lastyearmax and 
       !!      ::maxfpc_lastyear) are just replaced with this year's value.
       !
       !NVMODIF
       maxmoiavail_lastyear(:,:) = (maxmoiavail_lastyear(:,:)*(tau_climatology-1)+ maxmoiavail_thisyear(:,:))/tau_climatology
       minmoiavail_lastyear(:,:) = (minmoiavail_lastyear(:,:)*(tau_climatology-1)+ minmoiavail_thisyear(:,:))/tau_climatology
       maxgppweek_lastyear(:,:) =( maxgppweek_lastyear(:,:)*(tau_climatology-1)+ maxgppweek_thisyear(:,:))/tau_climatology
       !       maxmoiavail_lastyear(:,:) = maxmoiavail_thisyear(:,:)
       !       minmoiavail_lastyear(:,:) = minmoiavail_thisyear(:,:)
       !       maxgppweek_lastyear(:,:) = maxgppweek_thisyear(:,:)
       
       gdd0_lastyear(:) = gdd0_thisyear(:)

       precip_lastyear(:) = precip_thisyear(:)

       lm_lastyearmax(:,:) = lm_thisyearmax(:,:)

       maxfpc_lastyear(:,:) = maxfpc_thisyear(:,:)

       Tseason(:) = zero
       WHERE ( Tseason_length(:) .GT. min_sechiba )
          Tseason(:) = Tseason_tmp(:) / Tseason_length(:)
       ENDWHERE
       !
       !! 21.2 Reset new values for the "this year" variables. \n
       !!      The maximum variables are set to zero and the minimum variable ::minmoiavail_thisyear
       !!      is set to ::large_value (set to be 1.E33 in the module ::stomate_constants).
       !

       maxmoiavail_thisyear(:,:) = zero
       minmoiavail_thisyear(:,:) = large_value

       maxgppweek_thisyear(:,:) = zero

       gdd0_thisyear(:) = zero

       precip_thisyear(:) = zero

       lm_thisyearmax(:,:) = zero

       maxfpc_thisyear(:,:) = zero

       Tseason_tmp(:) = zero
       Tseason_length(:) =zero

       Tmin_spring_time(:,:)=zero
       onset_date(:,:)=zero
       !
       ! 21.3 Special treatment for maxfpc. !?? 
       !! 21.3 Set the maximum fractional plant cover for non-natural vegetation 
       !!      (i.e. agricultural C3 and C4 - PFT 12 and 13) vegetation for last year to be zero.
       !

       !
       ! 21.3.1 Only take into account natural PFTs
       !

       DO j = 2,nvm ! Loop over # PFTs
          IF ( .NOT. natural(j) .OR. pasture(j)) THEN
             maxfpc_lastyear(:,j) = zero
          ENDIF
       ENDDO

       ! 21.3.2 In Stomate, veget_cov is defined as a fraction of ground, not as a fraction 
       !        of total ground. maxfpc_lastyear will be compared to veget_cov in lpj_light.
       !        Therefore, we have to transform maxfpc_lastyear.


       ! 21.3.3 The sum of the maxfpc_lastyear for natural PFT must not exceed fpc_crit (=.95).
       !        However, it can slightly exceed this value as not all PFTs reach their maximum 
       !        fpc at the same time. Therefore, if sum(maxfpc_lastyear) for the natural PFTs
       !        exceeds fpc_crit, we scale the values of maxfpc_lastyear so that the sum is
       !        fpc_crit.

!!$       ! calculate the sum of maxfpc_lastyear
!!$       sumfpc_nat(:) = zero
!!$       DO j = 2,nvm ! Loop over # PFTs
!!$          sumfpc_nat(:) = sumfpc_nat(:) + maxfpc_lastyear(:,j)
!!$       ENDDO
!!$
!!$       ! scale so that the new sum is fpc_crit
!!$       DO j = 2,nvm ! Loop over # PFTs 
!!$          WHERE ( sumfpc_nat(:) .GT. fpc_crit )
!!$             maxfpc_lastyear(:,j) = maxfpc_lastyear(:,j) * (fpc_crit/sumfpc_nat(:))
!!$          ENDWHERE
!!$       ENDDO

    ENDIF  ! LastTsYear

    !
    !! 22. Diagnose herbivore activity (::herbivores), determined as probability for a leaf to be
    !!     eaten in a day (follows McNaughton et al., 1989). \n
    !!     The amount of herbivore activity is used in the modules ::lpj_establish  
    !!     and ::stomate_turnover.
    !

    !
    !! 22.1 First calculate the mean long-term leaf NPP in grid box, mean residence
    !!      time (years) of green tissue to give the biomass available for herbivore consumption.
    !!      Each PFT is looped over, though the calculation is only made for
    !!      natural vegetation (PFTs 2-11).
    !

    nlflong_nat(:,:) = zero
    weighttot(:,:) = zero
    green_age(:,:) = zero
    !
    DO j = 2,nvm ! Loop over # PFTs
       !
       IF ( natural(j) ) THEN
          !
          !! 22.1.1 Calculate the total weight of the leaves (::weighttot) as last year's leaf biomass
          !
          weighttot(:,j) = lm_lastyearmax(:,j)
          
          !
          !! 22.1.2 Calculate the mean long-term leaf NPP as the long-term NPP calculated in Section 12 of 
          !!        this subroutine weighted by the leaf fraction (::leaf_frac), which is defined to be 0.33
          !!        at the beginning of this subroutine.
          !
          nlflong_nat(:,j) = npp_longterm(:,j) * leaf_frac_hvc
          
          !
          !! 22.1.3 Calculate the mean residence time of the green tissue (::green_age) 
          !!        This is calculated as the sum of 6 months 
          !!        for natural seasonal vegetation (i.e. PFTs 3, 6, 8-11), and 2 years for evergreen 
          !!        (PFTs 2, 4, 5, 8), multiplied by last year's leaf biomass for each PFT, divided by the 
          !!        total weight of the leaves for all PFTs.
          !         This is a crude approximation.  !!?? By whom?
          !!        The difference between seasonal and evergreen vegetation is determined by the parameter
          !!        ::pheno_model, which specifies the onset model of the PFT.
 
          !
          IF ( pheno_model(j) .EQ. 'none' ) THEN
             green_age(:,j) = green_age_ever * lm_lastyearmax(:,j)
          ELSE
             green_age(:,j) = green_age_dec * lm_lastyearmax(:,j)
          ENDIF
          !
       ENDIF
       !
    ENDDO
    !
    WHERE ( weighttot(:,:) .GT. min_sechiba )
       green_age(:,:) = green_age(:,:) / weighttot(:,:)
    ELSEWHERE
       green_age(:,:) = un
    ENDWHERE

    !
    !! 22.2 McNaughton et al. (1989) give herbivore consumption as a function of mean, long-term leaf NPP.
    !!      as it gives an estimate of the edible biomass. The consumption of biomass by herbivores and the 
    !!      resultant herbivore activity are calculated following the equations:
    !!      \latexonly
    !!      \input{season_consumption_eqn4.tex}
    !!      \endlatexonly
    !!      \n 
    !!      and 
    !!      \latexonly
    !!      \input{season_herbivore_eqn5.tex}
    !!      \endlatexonly
    !!      \n 
    !       

    DO j = 2,nvm ! Loop over # PFTs
       !
       IF ( natural(j) ) THEN
          !
          WHERE ( nlflong_nat(:,j) .GT. zero )
             consumption(:) = hvc1 * nlflong_nat(:,j) ** hvc2
             herbivores(:,j) = one_year * green_age(:,j) * nlflong_nat(:,j) / consumption(:)
          ELSEWHERE
             herbivores(:,j) = 100000.
          ENDWHERE
          !
       ELSE
          !
          herbivores(:,j) = 100000.
          !
       ENDIF
       !
    ENDDO
    herbivores(:,ibare_sechiba) = zero

    IF (printlev>=4) WRITE(numout,*) 'Leaving season'

  END SUBROUTINE season

END MODULE stomate_season
