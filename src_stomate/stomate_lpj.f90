! ================================================================================================================================
! MODULE       : stomate_lpj
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Main entry point for daily processes in STOMATE and LPJ (phenology, 
!! allocation, npp_calc, kill, turn, light, establish, crown, cover, lcchange)
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S) : None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_stomate/stomate_lpj.f90 $
!! $Date: 2018-06-12 14:54:59 +0200 (Tue, 12 Jun 2018) $
!! $Revision: 5290 $
!! \n
!_ ================================================================================================================================

MODULE stomate_lpj

  ! modules used:

  USE ioipsl_para
  USE xios_orchidee
  USE grid
  USE stomate_data
  USE constantes
  USE constantes_soil
  USE pft_parameters
  USE lpj_constraints
  USE lpj_pftinout
  USE lpj_kill
  USE lpj_crown
  USE lpj_fire
  USE lpj_spitfire
  USE lpj_gap
  USE lpj_light
  USE lpj_establish
  USE lpj_cover
  USE stomate_prescribe
  USE stomate_phenology
  USE stomate_alloc
  USE stomate_npp
  USE stomate_turnover
  USE stomate_litter
  USE stomate_soilcarbon
  USE stomate_vmax
  USE stomate_lcchange
  USE stomate_gluc_common
  USE stomate_glcchange_SinAgeC
  USE stomate_glcchange_MulAgeC
  USE stomate_fharvest_SinAgeC
  USE stomate_fharvest_MulAgeC
  USE stomate_glcc_bioe1
  USE stomate_lai
  USE stomate_gluc_constants
!gmjc
  USE grassland_management
!end gmjc
  USE stomate_check

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC StomateLpj,StomateLpj_clear

  LOGICAL, SAVE                         :: firstcall = .TRUE.             !! first call
!$OMP THREADPRIVATE(firstcall)
!gmjc
  ! flag that enable grazing
  LOGICAL, SAVE :: enable_grazing
!$OMP THREADPRIVATE(enable_grazing)
!end gmjc

! check mass balance for stomate 
  LOGICAL, SAVE :: STO_CHECK_MASSBALANCE = .FALSE. 
!$OMP THREADPRIVATE(STO_CHECK_MASSBALANCE)

CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : StomateLpj_clear
!!
!>\BRIEF        Re-initialisation of variable
!!
!! DESCRIPTION  : This subroutine reinitializes variables. To be used if we want to relaunch 
!! ORCHIDEE but the routine is not used in current version.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE StomateLpj_clear

    CALL prescribe_clear
    CALL phenology_clear
    CALL npp_calc_clear
    CALL turn_clear
    CALL soilcarbon_clear
    CALL constraints_clear
    CALL establish_clear
    CALL fire_clear
    CALL spitfire_clear
    CALL gap_clear
    CALL light_clear
    CALL pftinout_clear
    CALL alloc_clear
    
    CALL grassmanag_clear

  END SUBROUTINE StomateLpj_clear


!! ================================================================================================================================
!! SUBROUTINE   : StomateLPJ
!!
!>\BRIEF        Main entry point for daily processes in STOMATE and LPJ, structures the call sequence 
!!              to the different processes such as dispersion, establishment, competition and mortality of PFT's.
!! 
!! DESCRIPTION  : This routine is the main entry point to all processes calculated on a 
!! daily time step. Is mainly devoted to call the different STOMATE and LPJ routines 
!! depending of the ok_dgvm (is dynamic veg used) and lpj_constant_mortality (is background mortality used).
!! It also prepares the cumulative 
!! fluxes or pools (e.g TOTAL_M TOTAL_BM_LITTER etc...)
!!
!! This routine makes frequent use of "weekly", "monthly" and "long term" variables. Quotion is used because
!! by default "weekly" denotes 7 days, by default "monthly" denotes 20 days and by default "Long term" denotes
!! 3 years. dtslow refers to 24 hours (1 day).
!!
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): All variables related to stomate and required for LPJ dynamic vegetation mode.
!!
!! REFERENCE(S) : 
!! - Krinner, G., N. Viovy, N. de Noblet-Ducoudré, J. Ogeé, J. Polcher, P. Friedlingstein, P. Ciais, S. Sitch, 
!! and I. C. Prentice. 2005. A dynamic global vegetation model for studies of the coupled atmosphere-biosphere 
!! system. Global Biogeochemical Cycles 19:GB1015, doi:1010.1029/2003GB002199.
!! - Sitch, S., B. Smith, I. C. Prentice, A. Arneth, A. Bondeau, W. Cramer, J. O. Kaplan, S. Levis, W. Lucht, 
!! M. T. Sykes, K. Thonicke, and S. Venevsky. 2003. Evaluation of ecosystem dynamics, plant geography and 
!! terrestrial carbon cycling in the LPJ dynamic global vegetation model. Global Change Biology 9:161-185.
!!
!! FLOWCHART    : Update with existing flowchart from N Viovy (Jan 19, 2012)
!! \n
!_ ================================================================================================================================
 
  SUBROUTINE StomateLpj (npts, dt_days, &
       lalo, neighbours, resolution, contfrac, &
       clay, herbivores, &
       tsurf_daily, tsoil_daily, t2m_daily, t2m_min_daily, &
       !spitfire
       t2m_max_daily, precip_daily, wspeed_daily, lightn, popd, a_nd, &
       read_observed_ba, observed_ba, &
       read_cf_fine,cf_fine,read_cf_coarse,cf_coarse,read_ratio_flag,ratio_flag,read_ratio,ratio,date,&
       !endspit
       litterhum_daily, soilhum_daily, &
       maxmoiavail_lastyear, minmoiavail_lastyear, &
       gdd0_lastyear, precip_lastyear, &
       moiavail_month, moiavail_week, t2m_longterm, t2m_month, t2m_week, &
       tsoil_month, soilhum_month, &
       gdd_m5_dormance, gdd_from_growthinit, gdd_midwinter, ncd_dormance, ngd_minus5, &
       turnover_longterm, gpp_daily, gpp_week, &
       time_hum_min, hum_min_dormance, maxfpc_lastyear, resp_maint_part, &
       PFTpresent, age, fireindex, firelitter, &
       leaf_age, leaf_frac, biomass, ind, adapted, regenerate, &
       senescence, when_growthinit, &
       litterpart, litter, litter_avail, litter_not_avail, litter_avail_frac, &
       dead_leaves, carbon,carbon_surf, lignin_struc, &
       !spitfire
       ni_acc,fire_numday,fuel_1hr,fuel_10hr,fuel_100hr,fuel_1000hr, &
       lcc,bafrac_deforest_accu,emideforest_litter_accu,emideforest_biomass_accu,&              
       deforest_litter_remain,deforest_biomass_remain,&
       def_fuel_1hr_remain,def_fuel_10hr_remain,&         
       def_fuel_100hr_remain,def_fuel_1000hr_remain,&   
       !endspit
       veget_cov_max, veget_cov_max_new, npp_longterm, lm_lastyearmax, &
       veget_lastlight, everywhere, need_adjacent, RIP_time, &
       lai, npp_daily, turnover_daily, turnover_time,&
       control_moist, control_temp, soilcarbon_input, &
       co2_to_bm, co2_fire, resp_hetero, resp_maint, resp_growth, &
       height, deadleaf_cover, vcmax, &
       bm_to_litter, &
       prod10,prod100,flux10, flux100, &
       vegetnew_firstday, glccNetLCC, &
       glccSecondShift,glccPrimaryShift, &
       harvest_matrix, harvest_biomass, bound_spa, newvegfrac, glcc_pft, &
       convflux,cflux_prod10,cflux_prod100, &
       harvest_above, carb_mass_total, &
       fpc_max, Tseason, Tseason_length, Tseason_tmp, &
       Tmin_spring_time, begin_leaves, onset_date, &
       MatrixA,&
!!!! crop variables
       pdlai, slai, & 
       ! for crop allocation
       in_cycle, deltai, dltaisen, ssla, pgrain, deltgrain, reprac, &
       nger, nlev, ndrp,  nlax, &
       c_reserve, c_leafb, nmat, nrec, N_limfert, tday_counter, &
!!!! end crop variables, xuhui
       zz_coef_deep, deepC_a, deepC_s, deepC_p, & 
       tsurf_year, & !pss:-
!gmjc
       wshtotsum, sr_ugb, compt_ugb, nb_ani, grazed_frac, &
       import_yield, sla_age1, t2m_14, sla_calc, snowfall_daily, day_of_year, &
       when_growthinit_cut, nb_grazingdays, &
       EndOfYear, &
       moiavail_daily,tmc_topgrass_daily,fc_grazing,snowmass_daily,&
       after_snow, after_wet, wet1day, wet2day, &
!end gmjc
!!!qcj++ peatland
       tcarbon_acro,tcarbon_cato, carbon_acro, carbon_cato,height_acro,height_cato, &
       resp_acro_oxic_d, resp_acro_anoxic_d,resp_cato_d,acro_to_cato_d,litter_to_acro_d,&
       carbon_save,deepC_a_save,deepC_s_save,deepC_p_save,delta_fsave, &
       liqwt_max_lastyear,veget_cov_max_adjusted,npp0_cumul,fpeat_map, &
       veget_cov_max_new_peatdgvm, wtp_month, wtpmax_month,fwet_month, &
       wettile_dgvm)
    
  !! 0. Variable and parameter declaration

    !! 0.1 input

    INTEGER(i_std), INTENT(in)                                 :: npts                 !! Domain size (unitless)
    REAL(r_std), INTENT(in)                                    :: dt_days              !! Time step of Stomate (days)
    INTEGER(i_std), DIMENSION(npts,NbNeighb), INTENT(in)       :: neighbours           !! Indices of the 8 neighbours of each grid 
                                                                                       !! point [1=North and then clockwise] 
    REAL(r_std), DIMENSION(npts,2), INTENT(in)                 :: resolution           !! Resolution at each grid point (m)  
                                                                                       !! [1=E-W, 2=N-S] 
    REAL(r_std),DIMENSION(npts,2),INTENT(in)                   :: lalo                 !! Geographical coordinates (latitude,longitude)
                                                                                       !! for pixels (degrees)
    REAL(r_std),DIMENSION (npts), INTENT (in)   :: contfrac          !! Fraction of continent in the grid cell (unitless)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: clay                 !! Clay fraction (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: herbivores           !! Time constant of probability of a leaf to 
                                                                                       !! be eaten by a herbivore (days) 
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: tsurf_daily          !! Daily surface temperatures (K)
    REAL(r_std), DIMENSION(npts,nslm), INTENT(in)              :: tsoil_daily          !! Daily soil temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: t2m_daily            !! Daily 2 meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: t2m_min_daily        !! Daily minimum 2 meter temperatures (K)
    !spitfire
    INTEGER(i_std),INTENT(in)                            :: date               !! Date (days) 
    ! daily maximum 2 meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                       :: t2m_max_daily
    ! daily precip (mm/day)
    REAL(r_std), DIMENSION(npts), INTENT(in)                        :: precip_daily
    ! Wind speed 
    REAL(r_std), DIMENSION(npts), INTENT(in)                       :: wspeed_daily
    ! Lightning flash rate
    REAL(r_std), DIMENSION(npts), INTENT(in)                       :: lightn
    ! Population density rate
    REAL(r_std), DIMENSION(npts), INTENT(inout)                    :: popd !popd declared and allocated and input in slowproc.f90
    ! Flag for read in observed burned area
    LOGICAL, INTENT (in)                                           :: read_observed_ba
    ! Observed burned area
    REAL(r_std),DIMENSION (npts), INTENT (in)                      :: observed_ba

    ! Flag for read in observed burned area
    LOGICAL, INTENT (in)                                   :: read_cf_coarse
    ! Observed burned area
    REAL(r_std),DIMENSION (npts), INTENT (in)       :: cf_coarse
    ! Flag for read in observed burned area
    LOGICAL, INTENT (in)                                   :: read_cf_fine
    ! Observed burned area
    REAL(r_std),DIMENSION (npts), INTENT (in)       :: cf_fine
    ! Flag for read in observed burned area
    LOGICAL, INTENT (in)                                   :: read_ratio
    ! Observed burned area
    REAL(r_std),DIMENSION (npts), INTENT (in)       :: ratio
    ! Flag for read in observed burned area
    LOGICAL, INTENT (in)                                   :: read_ratio_flag
    ! Observed burned area
    REAL(r_std),DIMENSION (npts), INTENT (in)       :: ratio_flag

    !endspit
!!!qcj++ peatland
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements) :: biomass_old
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements) :: bm_to_litter_old
    REAL(r_std), DIMENSION(npts,nvm)                  :: co2_to_bm_old
    REAL(r_std), DIMENSION(npts,nvm)                  :: flux_co2_to_bm
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements) :: flux_bm_to_litter
    REAL(r_std), DIMENSION(npts,nvm)                  ::    Kwt

    REAL(r_std),DIMENSION(npts,nvm), INTENT(in)   :: carbon_acro
    REAL(r_std),DIMENSION(npts,nvm), INTENT(in)   :: carbon_cato
    REAL(r_std),DIMENSION(npts), INTENT(in)     :: height_acro
    REAL(r_std),DIMENSION(npts), INTENT(in)     :: height_cato
    REAL(r_std),DIMENSION(npts), INTENT(in)     :: tcarbon_acro
    REAL(r_std),DIMENSION(npts), INTENT(in)     :: tcarbon_cato
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)    ::resp_acro_oxic_d !!respiration of acrotelm( oxic ),daily
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)    ::resp_acro_anoxic_d !!respiration of acrotelm( anoxic ),daily
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)    ::resp_cato_d !!respiration of catotelm,daily
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)    ::litter_to_acro_d
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)    ::acro_to_cato_d
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)           ::veget_cov_max_adjusted
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)    :: veget_cov_max_new_peatdgvm
    REAL(r_std), DIMENSION(npts,ncarb,nvm), INTENT(inout)     :: carbon_save
    REAL(r_std), DIMENSION(npts,ndeep), INTENT(inout)         :: deepC_a_save
    REAL(r_std), DIMENSION(npts,ndeep), INTENT(inout)         :: deepC_s_save
    REAL(r_std), DIMENSION(npts,ndeep), INTENT(inout)         :: deepC_p_save
    REAL(r_std), DIMENSION(npts), INTENT(inout)               :: delta_fsave
    REAL(r_std), DIMENSION(npts), INTENT(in)                  :: liqwt_max_lastyear
  !  REAL(r_std),DIMENSION(npts,nvm,nparts,nelements)          :: biomass_remove
    REAL(r_std), DIMENSION(npts,nvm)                          :: harvest_bio
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: npp0_cumul
    REAL(r_std), DIMENSION(npts,nstm),INTENT(in)                   :: fpeat_map
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: wtp_month 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: wtpmax_month
    REAL(r_std), DIMENSION(npts),INTENT(in)                   :: fwet_month
    LOGICAL,DIMENSION (nstm), INTENT (in)                 :: wettile_dgvm


    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: litterhum_daily      !! Daily litter humidity (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts,nslm), INTENT(in)              :: soilhum_daily        !! Daily soil humidity (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: maxmoiavail_lastyear !! Last year's maximum moisture availability 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: minmoiavail_lastyear !! Last year's minimum moisture availability 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: gdd0_lastyear        !! Last year's GDD0 (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: precip_lastyear      !! Lastyear's precipitation 
                                                                                       !! @tex $(mm year^{-1})$ @endtex
                                                                                       !! to determine if establishment possible
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: moiavail_month       !! "Monthly" moisture availability (0 to 1, 
                                                                                       !! unitless) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: moiavail_week        !! "Weekly" moisture availability 
                                                                                       !! (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: t2m_longterm         !! "Long term" 2 meter reference 
                                                                                       !! temperatures (K) 
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: t2m_month            !! "Monthly" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: t2m_week             !! "Weekly" 2-meter temperatures (K)
    ! "seasonal" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: Tseason
    ! temporary variable to calculate Tseason
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: Tseason_length
    ! temporary variable to calculate Tseason
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: Tseason_tmp

    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: Tmin_spring_time     !! Number of days after begin_leaves (leaf onset) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: onset_date           !! Date in the year at when the leaves started to grow(begin_leaves)

    !pss:+
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: tsurf_year           ! annual surface temperatures (K)
    !pss:-
    REAL(r_std), DIMENSION(npts,nslm), INTENT(in)              :: tsoil_month          !! "Monthly" soil temperatures (K)
    REAL(r_std), DIMENSION(npts,nslm), INTENT(in)              :: soilhum_month        !! "Monthly" soil humidity
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)               :: gdd_m5_dormance      !! Growing degree days (K), threshold -5 deg 
                                                                                       !! C (for phenology) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)               :: gdd_from_growthinit  !! growing degree days, since growthinit for crops
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: gdd_midwinter        !! Growing degree days (K), since midwinter 
                                                                                       !! (for phenology) - this is written to the history files 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)               :: ncd_dormance         !! Number of chilling days (days), since 
                                                                                       !! leaves were lost (for phenology) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)               :: ngd_minus5           !! Number of growing days (days), threshold 
                                                                                       !! -5 deg C (for phenology) 
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(in)        :: turnover_longterm    !! "Long term" turnover rate  
                                                                                       !! @tex $(gC m^{-2} year^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: gpp_daily            !! Daily gross primary productivity  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)       :: gpp_week                !! Mean weekly gross primary productivity 
                                                                                !! @tex $(gC m^{-2} day^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)               :: time_hum_min         !! Time elapsed since strongest moisture 
                                                                                       !! availability (days) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)             :: hum_min_dormance    !! minimum moisture during dormance 
                                                                              !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: maxfpc_lastyear      !! Last year's maximum foliage projected
                                                                                       !! coverage for each natural PFT,
                                                                                       !! @tex $(m^2 m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(in)        :: resp_maint_part      !! Maintenance respiration of different 
                                                                                       !! plant parts  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: fpc_max              !! "Maximal" coverage fraction of a PFT (LAI 
                                                                                       !! -> infinity) on ground  
                                                                                       !! @tex $(m^2 m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(ndeep),   INTENT (in)               :: zz_coef_deep         !! deep vertical profile
!!!!! crop variables
    ! FOR CROP---STICS
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)                  :: pdlai
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)               :: slai
   
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                 :: in_cycle
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: deltai
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: dltaisen
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)              :: ssla
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: pgrain
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: deltgrain
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: reprac
    INTEGER(i_std), DIMENSION(npts,nvm), INTENT(in)              :: nger
    INTEGER(i_std), DIMENSION(npts,nvm), INTENT(in)              :: nlev
    INTEGER(i_std), DIMENSION(npts,nvm), INTENT(in)              :: ndrp
    INTEGER(i_std), DIMENSION(npts,nvm), INTENT(in)              :: nmat
    INTEGER(i_std), DIMENSION(npts,nvm), INTENT(in)              :: nlax

!    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: N_limfert !!
!    defined already in GM
    INTEGER(i_std), INTENT(in)                                :: tday_counter

!!!!! end crop variables, xuhui

!gmjc
LOGICAL, INTENT(in)                                        :: EndOfYear
!end gmjc
  !! 0.2 Output variables
    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: npp_daily            !! Net primary productivity 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(out) :: turnover_daily       !! Turnover rates 
                                                                                       !! @tex $(gC m^{-2} year^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: co2_to_bm            !! CO2 taken up from atmosphere when 
                                                                                       !! introducing a new PFT (introduced for 
                                                                                       !! carbon balance closure) 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: co2_fire             !! Carbon emitted into the atmosphere by 
                                                                                       !! fire (living and dead biomass)  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: resp_hetero          !! Heterotrophic respiration
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: resp_maint           !! Maintenance respiration  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: resp_growth          !! Growth respiration  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    
    REAL(r_std), DIMENSION(npts), INTENT(inout)                :: deadleaf_cover       !! Fraction of soil covered by dead leaves 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: vcmax                !! Maximum rate of carboxylation 
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(out):: bm_to_litter      !! Conversion of biomass to litter 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    LOGICAL, DIMENSION(npts,nvm), INTENT(out)                  :: begin_leaves         !! signal to start putting leaves on (true/false)

!!!!! crop variables
    ! for crop c pools
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)             :: c_reserve
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)             :: c_leafb

    ! for crop turnover
    INTEGER(i_std), DIMENSION(npts,nvm), INTENT(in)             :: nrec   !harvest date
!!!!! end crop variables, xuhui 

    !! 0.3 Modified variables
    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: height               !! Height of vegetation (m) 
    REAL(r_std), DIMENSION(npts,nlevs), INTENT(inout)          :: control_moist        !! Moisture control of heterotrophic 
                                                                                       !! respiration (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,nlevs), INTENT(inout)          :: control_temp         !! Temperature control of heterotrophic 
                                                                                       !! respiration, above and below 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,ncarb,nvm), INTENT(inout)      :: soilcarbon_input     !! Quantity of carbon going into carbon 
                                                                                       !! pools from litter decomposition  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: lai                  !! Leaf area index OF AN INDIVIDUAL PLANT,
										       !! where a PFT contains n indentical plants
										       !! i.e., using the mean individual approach 
                                                                                       !! @tex $(m^2 m^{-2})$ @endtex 
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)                :: PFTpresent           !! Tab indicating which PFTs are present in 
                                                                                       !! each pixel 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: age                  !! Age (years)    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: fireindex            !! Probability of fire (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: firelitter           !! Longer term litter above the ground that 
                                                                                       !! can be burned, @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout)  :: leaf_age             !! Leaf age (days)
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout)  :: leaf_frac            !! Fraction of leaves in leaf age class, 
                                                                                       !! (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: biomass        !! Biomass @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: ind                  !! Density of individuals 
                                                                                       !! @tex $(m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: adapted              !! Adaptation of PFT (killed if too cold) 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: regenerate           !! "Fitness": Winter sufficiently cold for 
                                                                                       !! PFT regeneration ? (0 to 1, unitless) 
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)                :: senescence           !! Flag for setting senescence stage (only 
                                                                                       !! for deciduous trees) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: when_growthinit      !! How many days ago was the beginning of 
                                                                                       !! the growing season (days) 
    REAL(r_std), DIMENSION(npts,nvm,nlitt), INTENT(inout)      :: litterpart           !! Fraction of litter above the ground 
                                                                                       !! belonging to different PFTs
                                                                                       !! (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts,nlitt,nvm,nlevs,nelements), INTENT(inout):: litter     !! Metabolic and structural litter, above 
                                                                                       !! and below ground 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
!gmjc for grazing litter
    REAL(r_std), DIMENSION(npts,nlitt,nvm), INTENT(out):: litter_avail
    REAL(r_std), DIMENSION(npts,nlitt,nvm) , INTENT(out):: litter_not_avail
    REAL(r_std), DIMENSION(npts,nlitt,nvm), INTENT(in):: litter_avail_frac
!end gmjc
    REAL(r_std), DIMENSION(npts,nvm,nlitt), INTENT(inout)      :: dead_leaves          !! Dead leaves on ground, per PFT, metabolic 
                                                                                       !! and structural,  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,ncarb,nvm), INTENT(inout)      :: carbon               !! Carbon pool: active, slow, or passive, 
                                                                                       !! @tex $(gC m^{-2})$ @endtex  
    REAL(r_std), DIMENSION(npts,ncarb,nvm), INTENT(inout)      :: carbon_surf
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nlevs), INTENT(inout)      :: lignin_struc         !! Ratio of Lignin/Carbon in structural 
                                                                                       !! litter, above and below ground,  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: veget_cov_max        !! "Maximal" coverage fraction of a PFT (LAI 
                                                                                       !! -> infinity) on ground 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: veget_cov_max_new       !! "Maximal" coverage fraction of a PFT  
                                                                                       !! (LAI-> infinity) on ground (unitless) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: npp_longterm         !! "Long term" mean yearly primary 
                                                                                       !! productivity 
                                                                                       !! @tex $(gC m^{-2} year^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: lm_lastyearmax       !! Last year's maximum leaf mass, for each 
                                                                                       !! PFT @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: veget_lastlight      !! Vegetation fractions (on ground) after 
                                                                                       !! last light competition  
                                                                                       !! @tex $(m^2 m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: everywhere           !! Is the PFT everywhere in the grid box or 
                                                                                       !! very localized (after its introduction) 
                                                                                       !! (unitless) 
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)                :: need_adjacent        !! In order for this PFT to be introduced, 
                                                                                       !! does it have to be present in an 
                                                                                       !! adjacent grid box? 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: RIP_time             !! How much time ago was the PFT eliminated 
                                                                                       !! for the last time (y) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: turnover_time        !! Turnover_time of leaves for grasses 
                                                                                       !! (days)
    REAL(r_std), DIMENSION(npts,nvm),INTENT(in)             :: vegetnew_firstday        !! New "maximal" coverage fraction of a PFT 
                                                                                       !! (LAI -> infinity) (unitless) 
    REAL(r_std),DIMENSION(npts,12),INTENT(inout)               :: glccNetLCC      
    REAL(r_std),DIMENSION(npts,12),INTENT(inout)               :: glccSecondShift 
    REAL(r_std),DIMENSION(npts,12),INTENT(inout)               :: glccPrimaryShift  
    REAL(r_std), DIMENSION(npts,12),INTENT(inout)              :: harvest_matrix       !! The gross land use change matrix in case 
                                                                                       !! of gross land cover change is simulated.
    REAL(r_std), DIMENSION(npts,12),INTENT(in)              :: harvest_biomass       !! The gross land use change matrix in case 
                                                                                       !! of gross land cover change is simulated.
    REAL(r_std), DIMENSION(npts,nvm),INTENT(inout)              :: bound_spa           !! The gross land use change matrix in case 
                                                                                       !! of gross land cover change is simulated.
    REAL(r_std), DIMENSION(npts,nvmap),INTENT(in)               :: newvegfrac           !! The gross land use change matrix in case 
                                                                                       !! of gross land cover change is simulated.
    REAL(r_std),DIMENSION(npts,0:10,nwp), INTENT(inout)            :: prod10               !! Products remaining in the 10
                                                                                       !! year-turnover pool after the annual 
                                                                                       !! release for each compartment (10
                                                                                       !! + 1 : input from year of land cover 
                                                                                       !! change) @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts,0:100,nwp), INTENT(inout)           :: prod100              !! Products remaining in the 100 
                                                                                       !! year-turnover pool after the annual 
                                                                                       !! release for each compartment (100 
                                                                                       !! + 1 : input from year of land cover 
                                                                                       !! change) @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts,10,nwp), INTENT(inout)              :: flux10               !! Annual release from the 10
                                                                                       !! year-turnover pool compartments  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts,100,nwp), INTENT(inout)             :: flux100              !! Annual release from the 100 
                                                                                       !! year-turnover pool compartments  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts,nwp), INTENT(inout)                 :: convflux             !! Release during first year following land 
                                                                                       !! cover change @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts,nwp), INTENT(inout)                 :: cflux_prod10         !! Total annual release from the 10 
                                                                                       !! year-turnover pool 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts,nwp), INTENT(inout)                 :: cflux_prod100        !! Total annual release from the 100 
                                                                                       !! year-turnover pool 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts), INTENT(inout)                :: harvest_above        !! Harvest above ground biomass for 
                                                                                       !! agriculture @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts), INTENT(inout)                :: carb_mass_total      !! Carbon Mass total (soil, litter, veg) 
                                                                                       !! @tex $(gC m^{-2})$ @endtex  
    REAL(r_std), DIMENSION(npts,nvm,nbpools,nbpools), INTENT(inout) :: MatrixA         !! Matrix containing the fluxes  
                                                                                       !! between the carbon pools
                                                                                       !! per sechiba time step 
                                                                                       !! @tex $(gC.m^2.day^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,ndeep,nvm), INTENT(inout)         :: deepC_a           !! permafrost soil carbon (g/m**3) active
    REAL(r_std), DIMENSION(npts,ndeep,nvm), INTENT(inout)         :: deepC_s           !! permafrost soil carbon (g/m**3) slow
    REAL(r_std), DIMENSION(npts,ndeep,nvm), INTENT(inout)         :: deepC_p           !! permafrost soil carbon (g/m**3) passive
!gmjc

!glcc
    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: glcc_pft       !! a temporary variable to hold the fraction of ipft->ivma, i.e., from
                                                                              !! PFT_{ipft} to the youngest age class of MTC_{ivma}
!spitfire
    ! Nesterov index accumulated
    REAL(r_std), DIMENSION(npts), INTENT(inout)                           :: ni_acc
    REAL(r_std), DIMENSION(npts), INTENT(inout)                           :: fire_numday
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)                       :: bafrac_deforest_accu 
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements), INTENT(inout)       :: emideforest_litter_accu 
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout)      :: emideforest_biomass_accu 
    REAL(r_std), DIMENSION(npts,nlitt,nvm,nlevs,nelements), INTENT(inout) :: deforest_litter_remain   !! Vegetmax-weighted remaining litter on the ground for 
                                                                                                      !! deforestation region.
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout)      :: deforest_biomass_remain  !! Vegetmax-weighted remaining biomass on the ground for 
                                                                                                      !! deforestation region.
    ! fuel classes (1, 10, 100, 1000 hours)
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements), INTENT(inout)                 :: fuel_1hr
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements), INTENT(inout)                 :: fuel_10hr
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements), INTENT(inout)                 :: fuel_100hr
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements), INTENT(inout)                 :: fuel_1000hr
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements), INTENT(inout)                 :: def_fuel_1hr_remain
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements), INTENT(inout)                 :: def_fuel_10hr_remain
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements), INTENT(inout)                 :: def_fuel_100hr_remain
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements), INTENT(inout)                 :: def_fuel_1000hr_remain
    REAL(r_std), DIMENSION(npts)                                         :: d_area_burnt
    REAL(r_std), DIMENSION(npts)                                         :: d_numfire
    REAL(r_std), DIMENSION(npts)                                         :: fc_crown
    ! parameter for potential human-caused ignitions, ignitions ind^{-1}day{-1}, used in lpj_spitfire.f90
    REAL(r_std), DIMENSION(npts), INTENT(in)                             :: a_nd
!endspit
!gmjc
   ! snowfall_daily (mm/d?)
    REAL(r_std), DIMENSION(npts), INTENT(in)         :: snowfall_daily
   ! snowmass_daily (kg/m2)
    REAL(r_std), DIMENSION(npts), INTENT(in)         :: snowmass_daily
    ! "14days" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)         ::  t2m_14
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  sla_calc
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  wshtotsum
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  sr_ugb
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  compt_ugb
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  nb_ani
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  grazed_frac
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  import_yield
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  sla_age1
    INTEGER(i_std), INTENT(in)                       ::  day_of_year
    REAL(r_std), DIMENSION(:,:), INTENT(inout)       ::  N_limfert
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  when_growthinit_cut
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  nb_grazingdays
!    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  resp_hetero_litter_d
!    REAL(r_std), DIMENSION(npts,ncarb,nvm), INTENT(inout)  :: resp_hetero_soil_d
! top 5 layer grassland soil moisture for grazing
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)           :: moiavail_daily
!! Daily moisture availability (0-1, unitless)
    REAL(r_std),DIMENSION (npts), INTENT(in)       :: tmc_topgrass_daily
    REAL(r_std),DIMENSION (npts), INTENT(in)       :: fc_grazing
    REAL(r_std),DIMENSION (npts), INTENT(inout)    :: after_snow
    REAL(r_std),DIMENSION (npts), INTENT(inout)    :: after_wet
    REAL(r_std),DIMENSION (npts), INTENT(inout)    :: wet1day
    REAL(r_std),DIMENSION (npts), INTENT(inout)    :: wet2day
!end gmjc
    !! 0.4 Local variables

    REAL(r_std), DIMENSION(npts,nvm,nelements)                  :: tot_bm_to_litter    !! Total conversion of biomass to litter 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nelements)                  :: tot_live_biomass    !! Total living biomass  
                                                                                       !! @tex $(gC m{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements)           :: bm_alloc            !! Biomass increase, i.e. NPP per plant part 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nelements)                  :: tot_turnover        !! Total turnover rate  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                            :: tot_litter_soil_carb!! Total soil and litter carbon  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                            :: tot_litter_carb     !! Total litter carbon 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                            :: tot_soil_carb       !! Total soil carbon  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts)                                :: carb_mass_variation !! Carbon Mass variation  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                            :: cn_ind              !! Crown area of individuals 
                                                                                       !! @tex $(m^{2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                            :: woodmass_ind        !! Woodmass of individuals (gC) 
    REAL(r_std), DIMENSION(npts,nvm,nparts)                     :: f_alloc             !! Fraction that goes into plant part 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,nstm)                           :: avail_tree          !! Space availability for trees 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,nstm)                           :: avail_grass         !! Space availability for grasses 
                                                                                       !! (0 to 1, unitless) 
    INTEGER                                                     :: i,j,ivm,ivma,ipts,ilev,ilitt
    REAL(r_std),DIMENSION(npts)                                 :: prod10_total        !! Total products remaining in the pool 
                                                                                       !! after the annual release 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts)                                 :: prod100_total       !! Total products remaining in the pool 
                                                                                       !! after the annual release 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts)                                 :: cflux_prod_total    !! Total flux from conflux and the 10/100 
                                                                                       !! year-turnover pool 
                                                                                       !! @tex $(gC m^{-2} year^{-1})$ @endtex 
    REAL(r_std),DIMENSION(npts)                                 :: prod10_harvest_total!! Total products remaining in the pool 
                                                                                       !! after the annual release 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts)                                 :: prod100_harvest_total!! Total products remaining in the pool 
                                                                                       !! after the annual release 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts)                                 :: cflux_prod_harvest_total!! Total flux from conflux and the 10/100 
                                                                                       !! year-turnover pool 
                                                                                       !! @tex $(gC m^{-2} year^{-1})$ @endtex 
    REAL(r_std),DIMENSION(npts,nvm)                             :: veget_cov_max_tmp   !! "Maximal" coverage fraction of a PFT  
    REAL(r_std), DIMENSION(npts)                                :: area_land_m2        !! Land area within gridcel excluding water body [m2]
                                                                                       !! (LAI-> infinity) on ground (unitless) 
!!!!! crop
    REAL(r_std), DIMENSION(npts,nvm)                            :: crop_export         !! Cropland export (harvest & straw)
!!!!! end crop, xuhui
    REAL(r_std), DIMENSION(npts,nvm)                            :: mortality           !! Fraction of individual dying this time 
                                                                                       !! step (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts)                                :: vartmp              !! Temporary variable used to add history
    REAL(r_std), DIMENSION(npts,nvm)                            :: histvar             !! History variables
    REAL(r_std), DIMENSION(npts,nvm)                            :: lcc                 !! land cover change, i.e., loss of each PFT,
                                                                                       !! positive value indicates loss.
    REAL(r_std),DIMENSION(npts,nvm)                             :: deflitsup_total
    REAL(r_std),DIMENSION(npts,nvm)                             :: defbiosup_total
    !glcc
    REAL(r_std), DIMENSION(npts,nvm,nvmap)                      :: glcc_pftmtc    !! a temporary variable to hold the fractions each PFT is going to lose
    REAL(r_std), DIMENSION(npts,12)                             :: glccReal       !! The "real" glcc matrix that we apply in the model
                                                                                  !! after considering the consistency between presribed
                                                                                  !! glcc matrix and existing vegetation fractions.
    REAL(r_std), DIMENSION(npts,12)                              :: IncreDeficit   !! "Increment" deficits, negative values mean that 
                                                                                  !! there are not enough fractions in the source PFTs
                                                                                  !! /vegetations to target PFTs/vegetations. I.e., these
                                                                                  !! fraction transfers are presribed in LCC matrix but
                                                                                  !! not realized.
    REAL(r_std), DIMENSION(npts,12)                             :: glccZero       !! "Increment" deficits, negative values mean that 
    REAL(r_std), DIMENSION(npts)                       :: Deficit_pf2yf_final     !! 
    REAL(r_std), DIMENSION(npts)                       :: Deficit_sf2yf_final     !! 
    REAL(r_std), DIMENSION(npts)                       :: pf2yf_compen_sf2yf      !! 
    REAL(r_std), DIMENSION(npts)                       :: sf2yf_compen_pf2yf      !!
    REAL(r_std), DIMENSION(npts)                       :: fuelfrac      !!
    REAL(r_std)                                                 :: valtmp

!gmjc lcchange of managed grassland
    ! "maximal" coverage fraction of a PFT (LAI -> infinity) on ground
    INTEGER(i_std)                       :: ier
    LOGICAL                               :: l_error =.FALSE.
      ! variables to get closure carbon cycle for nbp
    REAL(r_std), DIMENSION(npts,nvm)                            :: harvest_gm
    REAL(r_std), DIMENSION(npts,nvm)                            :: ranimal_gm
    REAL(r_std), DIMENSION(npts,nvm)                            :: ch4_pft_gm
    REAL(r_std), DIMENSION(npts)                                :: ch4_gm
    REAL(r_std), DIMENSION(npts,nvm)                            :: cinput_gm
    REAL(r_std), DIMENSION(npts)                                :: co2_gm
    REAL(r_std),DIMENSION(npts,nvm)                             :: veget_cov_max_gm
    REAL(r_std), DIMENSION(npts)                                :: veget_exist_gm
    REAL(r_std),DIMENSION(npts,nvm)                             :: n2o_pft_gm
    REAL(r_std), DIMENSION(npts)                                :: n2o_gm
!end gmjc



    REAL(r_std), DIMENSION(npts,7) :: outflux_sta,outflux_end
    REAL(r_std), DIMENSION(npts,3) :: influx_sta,influx_end
    REAL(r_std), DIMENSION(npts,4) :: pool_sta,pool_end
    INTEGER(i_std) :: ind_biomass,ind_litter,ind_soil,ind_prod,ind_co2tobm,ind_gpp,ind_npp,&
                     ind_bm2lit,ind_resph,ind_respm,ind_respg,ind_convf,ind_cflux,ind_fire,&
                     l,m

    REAL(r_std), DIMENSION(npts,nvm,nelements)     :: mass_before
    REAL(r_std), DIMENSION(npts,nvm,nelements)     :: mass_change   !! positive
                               ! sign as inrease in mass
    REAL(r_std), DIMENSION(npts,nvm,nelements)     :: mass_after        !! Temporary variable 

    REAL(r_std), DIMENSION(npts,nelements)     :: mass_before_2rd
    REAL(r_std), DIMENSION(npts,nelements)     :: mass_change_2rd   !! positive
                               ! sign as inrease in mass
    REAL(r_std), DIMENSION(npts,nelements)     :: mass_after_2rd    !! Temporary variable 
    REAL(r_std), DIMENSION(npts,nelements)     :: mass_balance_2rd    !! Temporary variable 
    INTEGER :: f2g=1, f2p=2, f2c=3
    INTEGER :: g2f=4, g2p=5, g2c=6, p2f=7, p2g=8, p2c=9, c2f=10, c2g=11, c2p=12
!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering stomate_lpj'

  !! 1. Initializations

    lcc(:,:) = zero    
    glccZero = zero
    area_land_m2(:) = (resolution(:,1)*resolution(:,2)*contfrac(:))  !unit:m2
    glccReal(:,:) = zero
    glcc_pftmtc(:,:,:) = zero

    mass_before = zero
    mass_change = zero
    mass_after = zero

!gmjc
    IF (firstcall) THEN

        firstcall = .FALSE.

        !Config  Key  = GRM_ENABLE_GRAZING
        !Config  Desc = grazing allowed
        !Config  Def  = n
        !Config  Help = flag for choose if you want animals or not.
        !
        enable_grazing = .FALSE.
        CALL getin_p('GRM_ENABLE_GRAZING',enable_grazing)
        WRITE (numout,*) 'enable_grazing',enable_grazing
        WRITE (numout,*) 'manage',is_grassland_manag
        WRITE (numout,*) 'cut',is_grassland_cut
        WRITE (numout,*) 'grazed',is_grassland_grazed

        !Config  Key  = STO_CHECK_MASSBALANCE
        !Config  Desc = Check for mass balance in stomate_lpj
        !Config  Def  = n
        !Config  Help = flag to enable mass balance in stomate_lpj
        !
        STO_CHECK_MASSBALANCE = .FALSE.
        CALL getin_p('STO_CHECK_MASSBALANCE',STO_CHECK_MASSBALANCE)

        emideforest_biomass_accu(:,:,:,:) = zero
        emideforest_litter_accu(:,:,:,:) = zero
        bafrac_deforest_accu(:,:) = zero


      IF (do_now_stomate_lcchange) THEN
        IF (use_age_class) THEN
          IF (SingleAgeClass) THEN
          !  CALL glcc_SinAgeC_firstday(npts,veget_cov_max,newvegfrac,harvest_matrix, &
          !              glccSecondShift,glccPrimaryShift,glccNetLCC,&
          !              glccReal,glcc_pft,glcc_pftmtc,IncreDeficit, &
          !              Deficit_pf2yf_final, Deficit_sf2yf_final,   &
          !              pf2yf_compen_sf2yf, sf2yf_compen_pf2yf)

          ELSE
          !  CALL glcc_MulAgeC_firstday(npts,veget_cov_max,newvegfrac,harvest_matrix, &
          !              glccSecondShift,glccPrimaryShift,glccNetLCC,&
          !              glccReal,glcc_pft,glcc_pftmtc,IncreDeficit, &
          !              Deficit_pf2yf_final, Deficit_sf2yf_final,   &
          !              pf2yf_compen_sf2yf, sf2yf_compen_pf2yf)
          ENDIF
          ! We put only the conversion of tree->Notree as deforestation
          DO ipts = 1,npts
            DO ivm = 1,nvm
              DO ivma = 1,nvmap
                IF (is_tree(ivm) .AND. .NOT. is_tree(start_index(ivma))) THEN
                  lcc(ipts,ivm) = lcc(ipts,ivm) + glcc_pftmtc(ipts,ivm,ivma)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ELSE ! (.NOT. use_age_class), i.e., net land cover change.
          ! note here veget_cov_max is the last-year veget_cov_max; vegetnew_firstday
          ! is the veget_cov_max of next year, we need lcc to be >0 where forest
          ! area decreased, i.e., when vegetnew_firstday < veget_cov_max
          lcc(:,:) = veget_cov_max(:,:) - vegetnew_firstday(:,:)
        ENDIF

        IF (.NOT. allow_deforest_fire) lcc(:,:) = zero

        ! we creat the proxy that's needed for deforestation fire simulation.
        DO ipts = 1,npts
          DO ivm = 1,nvm
            deforest_litter_remain(ipts,:,ivm,:,:) = litter(ipts,:,ivm,:,:)*lcc(ipts,ivm)
            deforest_biomass_remain(ipts,ivm,:,:) = biomass(ipts,ivm,:,:)*lcc(ipts,ivm)
            def_fuel_1hr_remain(ipts,ivm,:,:) = fuel_1hr(ipts,ivm,:,:)*lcc(ipts,ivm)
            def_fuel_10hr_remain(ipts,ivm,:,:) = fuel_10hr(ipts,ivm,:,:)*lcc(ipts,ivm)
            def_fuel_100hr_remain(ipts,ivm,:,:) = fuel_100hr(ipts,ivm,:,:)*lcc(ipts,ivm)
            def_fuel_1000hr_remain(ipts,ivm,:,:) = fuel_1000hr(ipts,ivm,:,:)*lcc(ipts,ivm)
          ENDDO
        ENDDO
      ENDIF


    END IF !firstcall
    
    !! 1.1 Initialize variables to zero
    co2_to_bm(:,:) = zero
    co2_fire(:,:) = zero
    npp_daily(:,:) = zero
    resp_maint(:,:) = zero
    resp_growth(:,:) = zero
    harvest_above(:) = zero
    bm_to_litter(:,:,:,:) = zero
    cn_ind(:,:) = zero
    woodmass_ind(:,:) = zero
    turnover_daily(:,:,:,:) = zero
    crop_export(:,:) = zero
    deflitsup_total(:,:) = zero
    defbiosup_total(:,:) = zero
!gmjc
    !! Initialize gm variables for nbp to zero
    harvest_gm(:,:) = zero
    ranimal_gm(:,:) = zero
    ch4_pft_gm(:,:) = zero
    cinput_gm(:,:) = zero
    co2_gm(:) = zero
    ch4_gm(:) = zero
    n2o_gm(:) = zero
    n2o_pft_gm(:,:) = zero
    veget_cov_max_gm(:,:) = zero
    veget_exist_gm(:) = zero
!end gmjc

    ! GLUC
    IncreDeficit = zero    
    
    !! 1.2  Initialize variables to veget_cov_max
    veget_cov_max_tmp(:,:) = veget_cov_max(:,:)

    !! 1.3 Calculate some vegetation characteristics
    !! 1.3.1 Calculate some vegetation characteristics 
    !        Calculate cn_ind (individual crown mass) and individual height from
    !        state variables if running DGVM or dynamic mortality in static cover mode
    !??        Explain (maybe in the header once) why you mulitply with veget_cov_max in the DGVM
    !??        and why you don't multiply with veget_cov_max in stomate.
    IF ( ok_dgvm .OR. .NOT.lpj_gap_const_mort) THEN
       IF(ok_dgvm) THEN
          WHERE (ind(:,:).GT.min_stomate)
             woodmass_ind(:,:) = &
                  ((biomass(:,:,isapabove,icarbon)+biomass(:,:,isapbelow,icarbon) &
                  +biomass(:,:,iheartabove,icarbon)+biomass(:,:,iheartbelow,icarbon)) & 
                  *veget_cov_max(:,:))/ind(:,:)
          ENDWHERE
       ELSE
          WHERE (ind(:,:).GT.min_stomate)
             woodmass_ind(:,:) = &
                  (biomass(:,:,isapabove,icarbon)+biomass(:,:,isapbelow,icarbon) &
                  +biomass(:,:,iheartabove,icarbon)+biomass(:,:,iheartbelow,icarbon))/ind(:,:)
          ENDWHERE
       ENDIF

       CALL crown (npts,  PFTpresent, &
            ind, biomass, woodmass_ind, &
            veget_cov_max, cn_ind, height)
    ENDIF

!!!qcj++ peatland
    IF ( ok_dgvm_peat .AND. (.NOT. ok_dgvm) ) THEN
       DO j = 2,nvm
         IF (is_peat(j)) THEN
           WHERE (ind(:,j).GT.min_stomate)
             woodmass_ind(:,j) = &
                  ((biomass(:,j,isapabove,icarbon)+biomass(:,j,isapbelow,icarbon) &
                  +biomass(:,j,iheartabove,icarbon)+biomass(:,j,iheartbelow,icarbon)) &
                  *veget_cov_max(:,j))/ind(:,j)
           ENDWHERE
         ENDIF
       ENDDO

       CALL crown (npts,  PFTpresent, &
            ind, biomass, woodmass_ind, &
            veget_cov_max, cn_ind, height)
    ENDIF

    !! 1.3.2 Prescribe characteristics if the vegetation is not dynamic
    !        IF the DGVM is not activated, the density of individuals and their crown
    !        areas don't matter, but they should be defined for the case we switch on
    !        the DGVM afterwards. At the first call, if the DGVM is not activated, 
    !        impose a minimum biomass for prescribed PFTs and declare them present.


    !DSG mass conservation ========================================
    mass_before(:,:,:) = SUM(biomass(:,:,:,:),DIM=3)
 
    IF(STO_CHECK_MASSBALANCE) THEN
     DO j=2,nvm
       IF (is_peat(j)) THEN
         ! WRITE (numout,*) 'QCJ check before prescribe, PFT',j,'biomass',SUM(biomass(:,j,:,icarbon),DIM=2),'bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2),'co2_to_bm',co2_to_bm(:,j)
          biomass_old(:,j,:,:)=biomass(:,j,:,:)
          bm_to_litter_old(:,j,:,:)=bm_to_litter(:,j,:,:)
          co2_to_bm_old(:,j)=co2_to_bm(:,j)
       ENDIF
     ENDDO
    ENDIF

    
    CALL prescribe (npts, &
         veget_cov_max, dt_days, PFTpresent, everywhere, when_growthinit, &
         biomass, leaf_frac, ind, cn_ind, co2_to_bm)
  
    IF(STO_CHECK_MASSBALANCE) THEN
     DO i=1, npts
      DO j=2,nvm
       IF (is_peat(j)) THEN
          flux_co2_to_bm(:,j)=co2_to_bm(:,j)-co2_to_bm_old(:,j)
         ! WRITE (numout,*) 'QCJ check after prescribe, PFT',j,'biomass',SUM(biomass(:,j,:,icarbon),DIM=2),'bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2),'co2_to_bm',co2_to_bm(:,j)
         ! WRITE (numout,*) 'QCJ check after prescribe, PFT',j,'flux_co2_to_bm',flux_co2_to_bm(:,j),'delta biomass',SUM(biomass(:,j,:,icarbon),DIM=2)-SUM(biomass_old(:,j,:,icarbon),DIM=2)
          IF (ABS(SUM(biomass(i,j,:,icarbon))-flux_co2_to_bm(i,j)-SUM(biomass_old(i,j,:,icarbon))) .GT. min_stomate) THEN
             WRITE (numout,*) 'QCJ mass cons error, in prescribe'
          ENDIF
       ENDIF
      ENDDO
     ENDDO
    ENDIF

    IF(STO_CHECK_MASSBALANCE) THEN
       !DSG mass conservation ============================================
       CALL  stomate_check_mass_values(npts,biomass(:,:,:,:),'lpj: after prescribe')
       WRITE (numout,*) 'No mass cons test; prescribe violates mass per design' !:DSG right?
       !DSG mass_change(:,:,icarbon)     = zero
       !DSG CALL stomate_check_cons_mass(npts, nvm, nelements,     &  ! dimensions
       !DSG        mass_before(:,:,:),               &  ! mass before
       !DSG        SUM(biomass(:,:,:,:),DIM=3),      &  ! mass after
       !DSG        mass_change(:,:,:)                &  ! net of fluxes
       !DSG        )
    ENDIF

  !! 2. Climatic constraints for PFT presence and regenerativeness

    !   Call this even when DGVM is not activated so that "adapted" and "regenerate"
    !   are kept up to date for the moment when the DGVM is activated.
    CALL constraints (npts, dt_days, &
         t2m_month, t2m_min_daily,when_growthinit, &
         adapted, regenerate, Tseason, &
!qcj++ peatland
         wtp_month,Kwt)

    
  !! 3. Determine introduction and elimination of PTS based on climate criteria
 
    IF ( ok_dgvm .OR. ok_dgvm_peat) THEN !!!qcj++ peatland

       !DSG mass conservation ========================================
       mass_before(:,:,:) = SUM(biomass(:,:,:,:),DIM=3)

       IF(STO_CHECK_MASSBALANCE) THEN
        DO j=2,nvm
          IF (is_peat(j)) THEN
            ! WRITE (numout,*) 'QCJ check before pftinout, PFT',j,'biomass',SUM(biomass(:,j,:,icarbon),DIM=2),'bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2),'co2_to_bm',co2_to_bm(:,j)
            biomass_old(:,j,:,:)=biomass(:,j,:,:)
            bm_to_litter_old(:,j,:,:)=bm_to_litter(:,j,:,:)
            co2_to_bm_old(:,j)=co2_to_bm(:,j)
          ENDIF
        ENDDO
       ENDIF

       !! 3.1 Calculate introduction and elimination
       CALL pftinout (npts, dt_days, adapted, regenerate, &
            neighbours, veget_cov_max, &
            biomass, ind, cn_ind, age, leaf_frac, npp_longterm, lm_lastyearmax, senescence, &
            PFTpresent, everywhere, when_growthinit, need_adjacent, RIP_time, &
            co2_to_bm, &
            avail_tree, avail_grass, &
!gmjc
            sla_calc, &
!end gmjc 
            fpeat_map,fwet_month,wettile_dgvm) !!!qcj++ peatland

       IF (STO_CHECK_MASSBALANCE) THEN
          CALL  stomate_check_mass_values(npts,biomass(:,:,:,:), 'lpj: after pftinout')
          !DSG mass conservation ============================================
          WRITE (numout,*) 'MASS CONSERVATON: pftinout is designed to violate against mass conservation'
          WRITE (numout,*) 'but it should be not detectable according to comments'
          mass_change(:,:,icarbon)     = zero
          mass_after = SUM(biomass(:,:,:,:),DIM=3)

          CALL stomate_check_cons_mass(lalo, mass_before(:,:,:),              &  ! mass before
                 mass_after,                      &  ! mass after
                 mass_change(:,:,:),              &  ! net of fluxes
                 'lpj: after pftinout' )
       ENDIF

       IF(STO_CHECK_MASSBALANCE) THEN
         DO i=1, npts
          DO j=2,nvm
           IF (is_peat(j)) THEN
             flux_co2_to_bm(:,j)=co2_to_bm(:,j)-co2_to_bm_old(:,j)
             !  WRITE (numout,*) 'QCJ check after pftinout, PFT',j,'biomass',SUM(biomass(:,j,:,icarbon),DIM=2),'bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2),'co2_to_bm',co2_to_bm(:,j)
             !  WRITE (numout,*) 'QCJ check after pftinout, PFT',j,'flux_co2_to_bm',flux_co2_to_bm(:,j),'delta biomass',SUM(biomass(:,j,:,icarbon),DIM=2)-SUM(biomass_old(:,j,:,icarbon),DIM=2)
             IF ( ABS(SUM(biomass(i,j,:,icarbon))-flux_co2_to_bm(i,j)-SUM(biomass_old(i,j,:,icarbon))) .GT. min_stomate) THEN
               WRITE (numout,*) 'QCJ mass cons error, in pftinout'
             ENDIF
           ENDIF
          ENDDO
         ENDDO
       ENDIF

       !! 3.2 Reset attributes for eliminated PFTs.
       !     This also kills PFTs that had 0 leafmass during the last year. The message
       !     "... after pftinout" is misleading in this case.


       !DSG mass conservation ========================================
       mass_before(:,:,:) = SUM(biomass(:,:,:,:),DIM=3)

       IF(STO_CHECK_MASSBALANCE) THEN
        DO j=2,nvm
          IF (is_peat(j)) THEN
           !  WRITE (numout,*) 'QCJ check before kill pftinout, PFT',j,'biomass',SUM(biomass(:,j,:,icarbon),DIM=2),'bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2),'co2_to_bm',co2_to_bm(:,j)
           biomass_old(:,j,:,:)=biomass(:,j,:,:)
           bm_to_litter_old(:,j,:,:)=bm_to_litter(:,j,:,:)
           co2_to_bm_old(:,j)=co2_to_bm(:,j)
          ENDIF
        ENDDO
       ENDIF

       CALL kill (npts, 'pftinout  ', lm_lastyearmax, &
            ind, PFTpresent, cn_ind, biomass, senescence, RIP_time, &
            lai, age, leaf_age, leaf_frac, npp_longterm, &
            when_growthinit, everywhere, veget_cov_max, bm_to_litter)

       IF(STO_CHECK_MASSBALANCE) THEN
          CALL  stomate_check_mass_values(npts,biomass(:,:,:,:), 'lpj: after kill')
          !DSG mass conservation ============================================
          mass_change(:,:,icarbon)     = SUM(bm_to_litter(:,:,:,icarbon),DIM=3)
          mass_after = SUM(biomass(:,:,:,:),DIM=3)
          CALL stomate_check_cons_mass(lalo, mass_before(:,:,:),               &  ! mass before
                 mass_after,                       &  ! mass after
                 mass_change(:,:,:),               &  ! net of fluxes
                 'lpj: after kill' )
       ENDIF
   
       IF(STO_CHECK_MASSBALANCE) THEN
        DO i=1, npts
         DO j=2,nvm
          IF (is_peat(j)) THEN
           flux_bm_to_litter(:,j,:,:)=bm_to_litter(:,j,:,:)-bm_to_litter_old(:,j,:,:)
           !  WRITE (numout,*) 'QCJ check after kill pftinout, PFT',j,'biomass',SUM(biomass(:,j,:,icarbon),DIM=2),'bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2),'co2_to_bm',co2_to_bm(:,j)
           !  WRITE (numout,*) 'QCJ check after kill pftinout, PFT',j,'flux_bm_to_litter',SUM(flux_bm_to_litter(:,j,:,icarbon),DIM=2),'delta biomass',SUM(biomass(:,j,:,icarbon),DIM=2)-SUM(biomass_old(:,j,:,icarbon),DIM=2),'delta bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2)-SUM(bm_to_litter_old(:,j,:,icarbon),DIM=2)
           IF (ABS(SUM(biomass(i,j,:,icarbon))+SUM(flux_bm_to_litter(i,j,:,icarbon))-SUM(biomass_old(i,j,:,icarbon))) .GT. min_stomate) THEN
             WRITE (numout,*) 'QCJ mass cons error, in kill after pftinout'
           ENDIF
          ENDIF
         ENDDO
        ENDDO
       ENDIF

    ENDIF !( ok_dgvm .OR. ok_dgvm_peat)

    IF ( ok_dgvm ) THEN     
       !! 3.3 Calculate woodmass of individual tree
       IF(ok_dgvm) THEN
          WHERE ((ind(:,:).GT.min_stomate))
             woodmass_ind(:,:) = &
                  ((biomass(:,:,isapabove,icarbon)+biomass(:,:,isapbelow,icarbon) &
                  +biomass(:,:,iheartabove,icarbon)+biomass(:,:,iheartbelow,icarbon)) & 
                  *veget_cov_max(:,:))/ind(:,:)
          ENDWHERE
       ELSE
          WHERE (ind(:,:).GT.min_stomate)
             woodmass_ind(:,:) = &
                  (biomass(:,:,isapabove,icarbon)+biomass(:,:,isapbelow,icarbon) &
                  +biomass(:,:,iheartabove,icarbon)+biomass(:,:,iheartbelow,icarbon))/ind(:,:)
          ENDWHERE
       ENDIF

       ! Calculate crown area and diameter for all PFTs (including the newly established)
       CALL crown (npts, PFTpresent, &
            ind, biomass, woodmass_ind, &
            veget_cov_max, cn_ind, height)
    ENDIF
    
!!!qcj++ peatland
    IF ( ok_dgvm_peat .AND. (.NOT. ok_dgvm) ) THEN
       DO j = 2,nvm
         IF (is_peat(j)) THEN
           WHERE (ind(:,j).GT.min_stomate)
             woodmass_ind(:,j) = &
                  ((biomass(:,j,isapabove,icarbon)+biomass(:,j,isapbelow,icarbon) &
                  +biomass(:,j,iheartabove,icarbon)+biomass(:,j,iheartbelow,icarbon)) &
                  *veget_cov_max(:,j))/ind(:,j)
           ENDWHERE
         ENDIF
       ENDDO

       ! Calculate crown area and diameter for all PFTs (including the newly established)
       CALL crown (npts, PFTpresent, &
            ind, biomass, woodmass_ind, &
            veget_cov_max, cn_ind, height)
    ENDIF

  !! 4. Phenology

    !! 4.1 Write values to history file
    !      Current values for ::when_growthinit 
    CALL xios_orchidee_send_field("WHEN_GROWTHINIT",when_growthinit)

    CALL histwrite_p (hist_id_stomate, 'WHEN_GROWTHINIT', itime, when_growthinit, npts*nvm, horipft_index)

    ! Set and write values for ::PFTpresent
    WHERE(PFTpresent)
       histvar=un
    ELSEWHERE
       histvar=zero
    ENDWHERE

    CALL xios_orchidee_send_field("PFTPRESENT",histvar)

    CALL histwrite_p (hist_id_stomate, 'PFTPRESENT', itime, histvar, npts*nvm, horipft_index)

    ! Set and write values for gdd_midwinter
    WHERE(gdd_midwinter.EQ.undef)
       histvar=val_exp
    ELSEWHERE
       histvar=gdd_midwinter
    ENDWHERE

    CALL xios_orchidee_send_field("GDD_MIDWINTER",histvar)

    CALL histwrite_p (hist_id_stomate, 'GDD_MIDWINTER', itime, histvar, npts*nvm, horipft_index)

    ! Set and write values for gdd_m5_dormance
    WHERE(gdd_m5_dormance.EQ.undef)
       histvar=val_exp
    ELSEWHERE
       histvar=gdd_m5_dormance
    ENDWHERE
    
    CALL xios_orchidee_send_field('GDD_M5_DORMANCE',histvar)
    CALL histwrite_p (hist_id_stomate, 'GDD_M5_DORMANCE', itime, histvar, npts*nvm, horipft_index)

    ! Set and write values for ncd_dormance
    WHERE(ncd_dormance.EQ.undef)
       histvar=val_exp
    ELSEWHERE
       histvar=ncd_dormance
    ENDWHERE

    CALL xios_orchidee_send_field("NCD_DORMANCE",histvar)

    CALL histwrite_p (hist_id_stomate, 'NCD_DORMANCE', itime, histvar, npts*nvm, horipft_index)

    !DSG mass conservation ========================================
    mass_before(:,:,:) = SUM(biomass(:,:,:,:),DIM=3)

    !! 4.2 Calculate phenology
    IF(STO_CHECK_MASSBALANCE) THEN
     DO j=2,nvm
       IF (is_peat(j)) THEN
      !    WRITE (numout,*) 'QCJ check before phenology, PFT',j,'biomass',SUM(biomass(:,j,:,icarbon),DIM=2),'bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2),'co2_to_bm',co2_to_bm(:,j)
          biomass_old(:,j,:,:)=biomass(:,j,:,:)
          bm_to_litter_old(:,j,:,:)=bm_to_litter(:,j,:,:)
          co2_to_bm_old(:,j)=co2_to_bm(:,j)
       ENDIF
     ENDDO
    ENDIF

    CALL phenology (npts, dt_days, PFTpresent, &
         veget_cov_max, &
         t2m_longterm, t2m_month, t2m_week, gpp_daily, &
         maxmoiavail_lastyear, minmoiavail_lastyear, &
         moiavail_month, moiavail_week, &
         gdd_m5_dormance, gdd_midwinter, ncd_dormance, ngd_minus5, &
         senescence, time_hum_min, &
         biomass, leaf_frac, leaf_age, &
         when_growthinit, co2_to_bm, &
         pdlai, slai, deltai, ssla, &
         begin_leaves, &!)
!gmjc
         sla_calc)
!end gmjc

    IF(STO_CHECK_MASSBALANCE) THEN
     DO i=1, npts
      DO j=2,nvm
       !IF (is_peat(j)) THEN
        !WRITE (numout,*) 'QCJ check after phenology, PFT',j,'biomass',SUM(biomass(:,j,:,icarbon),DIM=2),'bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2),'co2_to_bm',co2_to_bm(:,j)
       !ENDIF  
       IF (is_peat(j)) THEN
          IF (ABS(SUM(biomass(i,j,:,icarbon))-SUM(biomass_old(i,j,:,icarbon))) .GT. min_stomate) THEN
             WRITE (numout,*) 'QCJ mass cons error, in phenology, biomass'
          ENDIF
          IF (ABS(SUM(bm_to_litter(i,j,:,icarbon))-SUM(bm_to_litter_old(i,j,:,icarbon))) .GT. min_stomate) THEN
             WRITE (numout,*) 'QCJ mass cons error, in phenology, bm_to_litter'
          ENDIF
          IF (ABS(co2_to_bm(i,j)-co2_to_bm_old(i,j)) .GT. min_stomate) THEN
             WRITE (numout,*) 'QCJ mass cons error, in phenology, co2_to_bm'
          ENDIF
       ENDIF
      ENDDO
     ENDDO
    ENDIF

    IF (STO_CHECK_MASSBALANCE) THEN
      
       CALL  stomate_check_mass_values(npts,biomass(:,:,:,:), 'lpj: after phenology_prognostic')
       WRITE (numout,*)'no mass conservation check after phenology: as the routine is designed to violate mass'
       !DSG as we sometimes have the case 'There is leaf or root carbon that
       !should not be here, something could have gone' see stomate_phenology.f90

       !DSG mass conservation ========================================
       mass_before(:,:,:)           = SUM(biomass(:,:,:,:),DIM=3)
    ENDIF
    
  !! 5. Allocate C to different plant parts
    !WRITE(numout,*) 'slai before lpj_alloc: ',slai(1,12:14)

    IF(STO_CHECK_MASSBALANCE) THEN
     DO j=2,nvm
       IF (is_peat(j)) THEN
        !  WRITE (numout,*) 'QCJ check before alloc, PFT',j,'biomass',SUM(biomass(:,j,:,icarbon),DIM=2),'bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2),'co2_to_bm',co2_to_bm(:,j)
          biomass_old(:,j,:,:)=biomass(:,j,:,:)
          bm_to_litter_old(:,j,:,:)=bm_to_litter(:,j,:,:)
          co2_to_bm_old(:,j)=co2_to_bm(:,j)
       ENDIF
     ENDDO
    ENDIF

    CALL alloc (npts, dt_days, &
         lai, veget_cov_max, senescence, when_growthinit, &
         moiavail_week, tsoil_month, soilhum_month, &
         biomass, age, leaf_age, leaf_frac, f_alloc, &!)
         deltai, ssla, & !added for crop by xuhui
!gmjc
         sla_calc, when_growthinit_cut)
!end gmjc

    IF(STO_CHECK_MASSBALANCE) THEN
     DO i=1, npts
      DO j=2,nvm
       IF (is_peat(j)) THEN
      !    WRITE (numout,*) 'QCJ check after alloc, PFT',j,'biomass',SUM(biomass(:,j,:,icarbon),DIM=2),'bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2),'co2_to_bm',co2_to_bm(:,j)
          IF (ABS(SUM(biomass(i,j,:,icarbon))-SUM(biomass_old(i,j,:,icarbon))) .GT. min_stomate) THEN
             WRITE (numout,*) 'QCJ mass cons error, in alloc, biomass'
          ENDIF
          IF (ABS(SUM(bm_to_litter(i,j,:,icarbon))-SUM(bm_to_litter_old(i,j,:,icarbon))) .GT. min_stomate) THEN
             WRITE (numout,*) 'QCJ mass cons error, in alloc, bm_to_litter'
          ENDIF
          IF (ABS(co2_to_bm(i,j)-co2_to_bm_old(i,j)) .GT. min_stomate) THEN
             WRITE (numout,*) 'QCJ mass cons error, in alloc, co2_to_bm'
          ENDIF
       ENDIF
      ENDDO
     ENDDO
    ENDIF

    !! 5.1. Recalculate lai
    !!      This should be done whenever biomass is modified
    CALL setlai(biomass, sla_calc, slai)

  !! 6. NPP, maintenance and growth respiration
    !! 6.1 Calculate NPP and respiration terms
    CALL npp_calc (npts, dt_days, &
         PFTpresent, &
         t2m_daily, tsoil_daily, lai, &
         gpp_daily, f_alloc, bm_alloc, resp_maint_part,&
         biomass, leaf_age, leaf_frac, age, &
         resp_maint, resp_growth, npp_daily, &!)
         ! for crop bm_alloc
!!! crop variables
         in_cycle, deltai, dltaisen, ssla, pgrain, deltgrain, reprac, &
         nger, nlev, ndrp, nlax, nmat, nrec, &
         c_reserve, c_leafb, slai, tday_counter, veget_cov_max, &
!!! end crop, xuhui
!gmjc
         sla_calc, sla_age1,N_limfert)
!end gmjc

    IF(STO_CHECK_MASSBALANCE) THEN
     DO j=2,nvm
       IF (is_peat(j)) THEN
        !  WRITE (numout,*) 'QCJ check before kill npp, PFT',j,'biomass',SUM(biomass(:,j,:,icarbon),DIM=2),'bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2),'co2_to_bm',co2_to_bm(:,j)
          biomass_old(:,j,:,:)=biomass(:,j,:,:)
          bm_to_litter_old(:,j,:,:)=bm_to_litter(:,j,:,:)
          co2_to_bm_old(:,j)=co2_to_bm(:,j)
       ENDIF
     ENDDO
    ENDIF
  
    !! 6.2 Kill slow growing PFTs in DGVM or STOMATE with constant mortality
    IF ( ok_dgvm .OR. .NOT.lpj_gap_const_mort .OR. ok_dgvm_peat) THEN !!!qcj++ peatland
       CALL kill (npts, 'npp       ', lm_lastyearmax,  &
            ind, PFTpresent, cn_ind, biomass, senescence, RIP_time, &
            lai, age, leaf_age, leaf_frac, npp_longterm, &
            when_growthinit, everywhere, veget_cov_max, bm_to_litter)
    ENDIF

       !! 6.2.1 Update wood biomass
       !        For the DGVM
    IF ( ok_dgvm .OR. .NOT.lpj_gap_const_mort) THEN 
       IF(ok_dgvm) THEN
          WHERE (ind(:,:).GT.min_stomate)
             woodmass_ind(:,:) = &
                  ((biomass(:,:,isapabove,icarbon) + biomass(:,:,isapbelow,icarbon) &
                  + biomass(:,:,iheartabove,icarbon) + biomass(:,:,iheartbelow,icarbon)) & 
                  *veget_cov_max(:,:))/ind(:,:)
          ENDWHERE

       ! For all pixels with individuals
       ELSE
          WHERE (ind(:,:).GT.min_stomate)
             woodmass_ind(:,:) = &
                  (biomass(:,:,isapabove,icarbon) + biomass(:,:,isapbelow,icarbon) &
                  + biomass(:,:,iheartabove,icarbon) + biomass(:,:,iheartbelow,icarbon))/ind(:,:)
          ENDWHERE
       ENDIF ! ok_dgvm

       !! 6.2.2 New crown area and maximum vegetation cover after growth
       CALL crown (npts, PFTpresent, &
            ind, biomass, woodmass_ind,&
            veget_cov_max, cn_ind, height)

    ENDIF ! ok_dgvm

    IF(STO_CHECK_MASSBALANCE) THEN
     DO i=1, npts
      DO j=2,nvm
       IF (is_peat(j)) THEN
          flux_bm_to_litter(i,j,:,:)=bm_to_litter(i,j,:,:)-bm_to_litter_old(i,j,:,:)
        !  WRITE (numout,*) 'QCJ check after kill npp, PFT',j,'biomass',SUM(biomass(:,j,:,icarbon),DIM=2),'bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2),'co2_to_bm',co2_to_bm(:,j)
        !  WRITE (numout,*) 'QCJ check after kill npp, PFT',j,'flux_bm_to_litter',SUM(flux_bm_to_litter(i,j,:,icarbon)),'delta biomass',SUM(biomass(i,j,:,icarbon))-SUM(biomass_old(:,j,:,icarbon),DIM=2),'delta bm_to_litter',SUM(bm_to_litter(i,j,:,icarbon))-SUM(bm_to_litter_old(:,j,:,icarbon),DIM=2)
          IF (ABS(SUM(biomass(i,j,:,icarbon))+SUM(flux_bm_to_litter(i,j,:,icarbon))-SUM(biomass_old(i,j,:,icarbon))) .GT. min_stomate) THEN
             WRITE (numout,*) 'QCJ mass cons error, in kill after npp'
          ENDIF
       ENDIF
       ENDDO
      ENDDO
    ENDIF
!!!
    IF (ok_dgvm_peat .AND. (.NOT. ok_dgvm) ) THEN
       DO j = 2,nvm
         IF (is_peat(j)) THEN
           WHERE (ind(:,j).GT.min_stomate)
             woodmass_ind(:,j) = &
                  ((biomass(:,j,isapabove,icarbon)+biomass(:,j,isapbelow,icarbon) &
                  +biomass(:,j,iheartabove,icarbon)+biomass(:,j,iheartbelow,icarbon)) &
                  *veget_cov_max(:,j))/ind(:,j)
           ENDWHERE
         ENDIF
       ENDDO

       !! 6.2.2 New crown area and maximum vegetation cover after growth
       CALL crown (npts, PFTpresent, &
            ind, biomass, woodmass_ind,&
            veget_cov_max, cn_ind, height)
    ENDIF

  !! 7. fire
    !! 7.1. Burn PFTs
    !CALL fire (npts, dt_days, litterpart, &
    !     litterhum_daily, t2m_daily, lignin_struc, veget_cov_max, &
    !     fireindex, firelitter, biomass, ind, &
    !     litter, dead_leaves, bm_to_litter, &
    !     co2_fire, MatrixA)

!gmjc update available and not available litter for grazing litter
!spitfire

        !disable_fire and allow_deforest_fire are defined as constants in src_parameters/constantes_var.f90
        !disable_fire to DISABLE fire when being TRUE
        !allow_deforest_fire to activate deforestation fire module when being TRUE
        IF(.NOT.disable_fire) THEN
           CALL spitfire(npts, dt_days, veget_cov_max,resolution,contfrac,   &
                PFTpresent,t2m_min_daily,t2m_max_daily,                  &    
                precip_daily,wspeed_daily,soilhum_daily(:,1),            &    
                lightn,litter(:,:,:,:,icarbon),ni_acc,fire_numday,       &
                fuel_1hr(:,:,:,icarbon),fuel_10hr(:,:,:,icarbon),fuel_100hr(:,:,:,icarbon), &    
                fuel_1000hr(:,:,:,icarbon),ind,biomass(:,:,:,icarbon),popd,a_nd,height,     &    
                read_observed_ba, observed_ba, read_cf_fine,cf_fine,                        &
                read_cf_coarse,cf_coarse,read_ratio_flag,                                   &
                ratio_flag,read_ratio,ratio,date,                                           &
                bm_to_litter(:,:,:,icarbon),co2_fire,                                       &
                lcc,bafrac_deforest_accu,emideforest_litter_accu(:,:,:,icarbon),emideforest_biomass_accu(:,:,:,icarbon),&
                deforest_litter_remain(:,:,:,:,icarbon),deforest_biomass_remain(:,:,:,icarbon),&
                def_fuel_1hr_remain(:,:,:,icarbon),def_fuel_10hr_remain(:,:,:,icarbon),&
                def_fuel_100hr_remain(:,:,:,icarbon),def_fuel_1000hr_remain(:,:,:,icarbon))
        ELSE
                ni_acc=0. 
                fire_numday=0.
                 
        ENDIF
!endspit
! after fire burning
  litter_avail(:,:,:) = litter(:,:,:,iabove,icarbon) * &
            litter_avail_frac(:,:,:)
  litter_not_avail(:,:,:) = litter(:,:,:,iabove,icarbon) * &
            (1.0 - litter_avail_frac(:,:,:))
!end gmjc
    !! 7.2 Kill PFTs in DGVM
    IF ( ok_dgvm .OR. ok_dgvm_peat) THEN

       ! reset attributes for eliminated PFTs
       CALL kill (npts, 'fire      ', lm_lastyearmax, &
            ind, PFTpresent, cn_ind, biomass, senescence, RIP_time, &
            lai, age, leaf_age, leaf_frac, npp_longterm, &
            when_growthinit, everywhere, veget_cov_max, bm_to_litter)

    ENDIF ! ok_dgvm

    IF (STO_CHECK_MASSBALANCE) THEN
       CALL stomate_check_mass_values(npts,biomass(:,:,:,:), 'lpj: after fire')
    ENDIF
  !! 8. Tree mortality
    !DSG mass conservation ========================================
    mass_before(:,:,:) = SUM(biomass(:,:,:,:),DIM=3) + SUM(bm_to_litter(:,:,:,:),DIM=3)

    IF(STO_CHECK_MASSBALANCE) THEN
     DO j=2,nvm
       IF (is_peat(j)) THEN
        !  WRITE (numout,*) 'QCJ check before gap, PFT',j,'biomass',SUM(biomass(:,j,:,icarbon),DIM=2),'bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2),'co2_to_bm',co2_to_bm(:,j)
          biomass_old(:,j,:,:)=biomass(:,j,:,:)
          bm_to_litter_old(:,j,:,:)=bm_to_litter(:,j,:,:)
          co2_to_bm_old(:,j)=co2_to_bm(:,j)
       ENDIF
     ENDDO
    ENDIF

    ! Does not depend on age, therefore does not change crown area.
    CALL gap (npts, dt_days, &
         npp_longterm, turnover_longterm, lm_lastyearmax, &
         PFTpresent, biomass, ind, bm_to_litter, mortality, t2m_min_daily, Tmin_spring_time, &!)
!gmjc
         sla_calc,Kwt,t2m_daily)
!end gmjc
    IF(STO_CHECK_MASSBALANCE) THEN
     DO i=1, npts
      DO j=2,nvm
       IF (is_peat(j)) THEN
          flux_bm_to_litter(:,j,:,:)=bm_to_litter(:,j,:,:)-bm_to_litter_old(:,j,:,:)
        !  WRITE (numout,*) 'QCJ check after gap, PFT',j,'biomass',SUM(biomass(:,j,:,icarbon),DIM=2),'bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2),'co2_to_bm',co2_to_bm(:,j)
        !  WRITE (numout,*) 'QCJ check after gap, PFT',j,'flux_bm_to_litter',SUM(flux_bm_to_litter(:,j,:,icarbon),DIM=2),'delta biomass',SUM(biomass(:,j,:,icarbon),DIM=2)-SUM(biomass_old(:,j,:,icarbon),DIM=2),'delta bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2)-SUM(bm_to_litter_old(:,j,:,icarbon),DIM=2)
          IF (ABS(SUM(biomass(i,j,:,icarbon))+SUM(flux_bm_to_litter(i,j,:,icarbon))-SUM(biomass_old(i,j,:,icarbon))) .GT. min_stomate) THEN
             WRITE (numout,*) 'QCJ mass cons error, in gap'
          ENDIF
       ENDIF
      ENDDO
     ENDDO
    ENDIF

    IF(STO_CHECK_MASSBALANCE) THEN
     DO j=2,nvm
       IF (is_peat(j)) THEN
        !  WRITE (numout,*) 'QCJ check before kill gap, PFT',j,'biomass',SUM(biomass(:,j,:,icarbon),DIM=2),'bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2),'co2_to_bm',co2_to_bm(:,j)
          biomass_old(:,j,:,:)=biomass(:,j,:,:)
          bm_to_litter_old(:,j,:,:)=bm_to_litter(:,j,:,:)
          co2_to_bm_old(:,j)=co2_to_bm(:,j)
       ENDIF
     ENDDO
    ENDIF

    IF ( ok_dgvm .OR. ok_dgvm_peat ) THEN

       ! reset attributes for eliminated PFTs
       CALL kill (npts, 'gap       ', lm_lastyearmax, &
            ind, PFTpresent, cn_ind, biomass, senescence, RIP_time, &
            lai, age, leaf_age, leaf_frac, npp_longterm, &
            when_growthinit, everywhere, veget_cov_max, bm_to_litter)

    ENDIF

    IF(STO_CHECK_MASSBALANCE) THEN
     DO i=1, npts
      DO j=2,nvm
       IF (is_peat(j)) THEN
          flux_bm_to_litter(:,j,:,:)=bm_to_litter(:,j,:,:)-bm_to_litter_old(:,j,:,:)
      !    WRITE (numout,*) 'QCJ check after kill gap, PFT',j,'biomass',SUM(biomass(:,j,:,icarbon),DIM=2),'bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2),'co2_to_bm',co2_to_bm(:,j)
      !    WRITE (numout,*) 'QCJ check after kill gap, PFT',j,'flux_bm_to_litter',SUM(flux_bm_to_litter(:,j,:,icarbon),DIM=2),'delta biomass',SUM(biomass(:,j,:,icarbon),DIM=2)-SUM(biomass_old(:,j,:,icarbon),DIM=2),'delta bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2)-SUM(bm_to_litter_old(:,j,:,icarbon),DIM=2)
          IF (ABS(SUM(biomass(i,j,:,icarbon))+SUM(flux_bm_to_litter(i,j,:,icarbon))-SUM(biomass_old(i,j,:,icarbon))) .GT. min_stomate) THEN
             WRITE (numout,*) 'QCJ mass cons error, in kill after gap'
          ENDIF
       ENDIF
      ENDDO
     ENDDO
    ENDIF

    IF(STO_CHECK_MASSBALANCE) THEN
       !DSG debug: ===================================================
       CALL stomate_check_mass_values(npts,biomass(:,:,:,:), 'lpj: after gap n kill')

       !DSG mass conservation ============================================
       mass_change(:,:,icarbon)     = -SUM(bm_to_litter(:,:,:,icarbon),DIM=3)

       mass_after = SUM(biomass(:,:,:,:),DIM=3)
       CALL stomate_check_cons_mass(lalo, mass_before(:,:,:),               &  ! mass before
             mass_after,                       &  ! mass after
             mass_change(:,:,:),               &  ! net of fluxes
             'lpj: after gap n kill')
     ENDIF


!! =====================================
!! DSG: mass conservation issue: START

  !! 10. Leaf senescence, new lai and other turnover processes
    !DSG mass conservation ========================================
    mass_before(:,:,:) = SUM(biomass(:,:,:,:),DIM=3)  + SUM(bm_to_litter(:,:,:,:),DIM=3)
    IF(STO_CHECK_MASSBALANCE) THEN
     DO j=2,nvm
       IF (is_peat(j)) THEN
       !   WRITE (numout,*) 'QCJ check before turn, PFT',j,'biomass',SUM(biomass(:,j,:,icarbon),DIM=2),'bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2),'co2_to_bm',co2_to_bm(:,j)
       !   WRITE (numout,*) 'QCJ check before turn, PFT,',j,'turnover_daily',SUM(turnover_daily(:,j,:,icarbon),DIM=2)
          biomass_old(:,j,:,:)=biomass(:,j,:,:)
       ENDIF
     ENDDO
    ENDIF

    CALL turn (npts, dt_days, PFTpresent, &
         herbivores, &
         maxmoiavail_lastyear, minmoiavail_lastyear, &
         moiavail_week,  moiavail_month,t2m_longterm, t2m_month, t2m_week, veget_cov_max, &
         gdd_from_growthinit, leaf_age, leaf_frac, age, lai, biomass, &
         turnover_daily, senescence,turnover_time, &!)
         nrec,crop_export, &
!gmjc
         sla_calc,&
         npp0_cumul) !!!qcj++ peatland
!end gmjc

    IF(STO_CHECK_MASSBALANCE) THEN
     DO i=1, npts
      DO j=2,nvm
       IF (is_peat(j)) THEN
       !   WRITE (numout,*) 'QCJ check after turn, PFT',j,'biomass',SUM(biomass(:,j,:,icarbon),DIM=2),'bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2),'co2_to_bm',co2_to_bm(:,j)
       !   WRITE (numout,*) 'QCJ check after turn, PFT',j,'turnover_daily',SUM(turnover_daily(:,j,:,icarbon),DIM=2),'delta bioamss',SUM(biomass(i,j,:,icarbon))-SUM(biomass_old(i,j,:,icarbon))
          IF ( ABS(SUM(biomass(i,j,:,icarbon))+SUM(turnover_daily(i,j,:,icarbon))-SUM(biomass_old(i,j,:,icarbon))) .GT. min_stomate) THEN
             WRITE (numout,*) 'QCJ mass cons error, in turn'
          ENDIF
       ENDIF
      ENDDO
     ENDDO
    ENDIF
    !! 10. Light competition
    
    !! If not using constant mortality then kill with light competition
!    IF ( ok_dgvm .OR. .NOT.(lpj_gap_const_mort) ) THEN
    IF ( ok_dgvm .OR. ok_dgvm_peat) THEN !!!qcj++ peatland

       IF(STO_CHECK_MASSBALANCE) THEN
        DO j=2,nvm
         IF (is_peat(j)) THEN
          !  WRITE (numout,*) 'QCJ check before light, PFT',j,'biomass',SUM(biomass(:,j,:,icarbon),DIM=2),'bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2),'co2_to_bm',co2_to_bm(:,j)
          biomass_old(:,j,:,:)=biomass(:,j,:,:)
          bm_to_litter_old(:,j,:,:)=bm_to_litter(:,j,:,:)
          co2_to_bm_old(:,j)=co2_to_bm(:,j)
         ENDIF
        ENDDO
       ENDIF

       !! 10.1 Light competition
       CALL light (npts, dt_days, &
            veget_cov_max, fpc_max, PFTpresent, cn_ind, lai, maxfpc_lastyear, &
            lm_lastyearmax, ind, biomass, veget_lastlight, bm_to_litter, mortality, &!)
!gmjc
         sla_calc, &
!end gmjc
         fpeat_map,wettile_dgvm) !!!qcj++ peatland

       IF(STO_CHECK_MASSBALANCE) THEN
        DO i=1, npts
         DO j=2,nvm
          IF (is_peat(j)) THEN
            flux_bm_to_litter(:,j,:,:)=bm_to_litter(:,j,:,:)-bm_to_litter_old(:,j,:,:)
            !  WRITE (numout,*) 'QCJ check after light,PFT',j,'biomass',SUM(biomass(:,j,:,icarbon),DIM=2),'bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2),'co2_to_bm',co2_to_bm(:,j)
            !  WRITE (numout,*) 'QCJ check after light,PFT',j,'flux_bm_to_litter',SUM(flux_bm_to_litter(:,j,:,icarbon),DIM=2),'delta biomass',SUM(biomass(:,j,:,icarbon),DIM=2)-SUM(biomass_old(:,j,:,icarbon),DIM=2),'delta bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2)-SUM(bm_to_litter_old(:,j,:,icarbon),DIM=2)
            IF (ABS(SUM(biomass(i,j,:,icarbon))+SUM(flux_bm_to_litter(i,j,:,icarbon))-SUM(biomass_old(i,j,:,icarbon))) .GT. min_stomate) THEN
             WRITE (numout,*) 'QCJ mass cons error, in light'
            ENDIF
          ENDIF
         ENDDO
        ENDDO
       ENDIF

       IF(STO_CHECK_MASSBALANCE) THEN
        DO j=2,nvm
         IF (is_peat(j)) THEN
          !  WRITE (numout,*) 'QCJ check before kill light, PFT',j,'biomass',SUM(biomass(:,j,:,icarbon),DIM=2),'bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2),'co2_to_bm',co2_to_bm(:,j)
          biomass_old(:,j,:,:)=biomass(:,j,:,:)
          bm_to_litter_old(:,j,:,:)=bm_to_litter(:,j,:,:)
          co2_to_bm_old(:,j)=co2_to_bm(:,j)
         ENDIF
        ENDDO
       ENDIF 
   
       !! 10.2 Reset attributes for eliminated PFTs
       CALL kill (npts, 'light     ', lm_lastyearmax, &
            ind, PFTpresent, cn_ind, biomass, senescence, RIP_time, &
            lai, age, leaf_age, leaf_frac, npp_longterm, &
            when_growthinit, everywhere, veget_cov_max, bm_to_litter)

       IF(STO_CHECK_MASSBALANCE) THEN
        DO i=1, npts
         DO j=2,nvm
          IF (is_peat(j)) THEN
           flux_bm_to_litter(:,j,:,:)=bm_to_litter(:,j,:,:)-bm_to_litter_old(:,j,:,:)
           !    WRITE (numout,*) 'QCJ check after kill light,PFT',j,'biomass',SUM(biomass(:,j,:,icarbon),DIM=2),'bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2),'co2_to_bm',co2_to_bm(:,j)
           !    WRITE (numout,*) 'QCJ check after kill light,PFT',j,'flux_bm_to_litter',SUM(flux_bm_to_litter(:,j,:,icarbon),DIM=2),'delta biomass',SUM(biomass(:,j,:,icarbon),DIM=2)-SUM(biomass_old(:,j,:,icarbon),DIM=2),'delta bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2)-SUM(bm_to_litter_old(:,j,:,icarbon),DIM=2)
           IF (ABS(SUM(biomass(i,j,:,icarbon))+SUM(flux_bm_to_litter(i,j,:,icarbon))-SUM(biomass_old(i,j,:,icarbon))) .GT. min_stomate) THEN
             WRITE (numout,*) 'QCJ mass cons error, in kill after light'
           ENDIF
          ENDIF
         ENDDO
        ENDDO
       ENDIF

    ENDIF !ok_dgvm .OR. ok_dgvm_peat

    IF(STO_CHECK_MASSBALANCE) THEN
       !DSG debug: ===================================================
       CALL  stomate_check_mass_values(npts,biomass(:,:,:,:), 'lpj: after turn')

       !DSG mass conservation ============================================
       mass_change(:,:,icarbon)     = -SUM(turnover_daily(:,:,:,icarbon),DIM=3) &
                                      -SUM(bm_to_litter(:,:,:,icarbon),DIM=3)  ! in case ok_dvgm=true

       mass_after = SUM(biomass(:,:,:,:),DIM=3)
       CALL stomate_check_cons_mass(lalo, mass_before(:,:,:),               &  ! mass before
             mass_after,                       &  ! mass afte
             mass_change(:,:,:),               &  ! net of fluxes
             'lpj: after turn')
    ENDIF

!! DSG: mass conservation issue: END
!! =====================================
    
  !! 11. Establishment of saplings
    
    IF ( ok_dgvm .OR. .NOT. lpj_gap_const_mort .OR. ok_dgvm_peat) THEN !!!qcj++ peatland
       !DSG mass conservation ========================================
       mass_before(:,:,:) = SUM(biomass(:,:,:,:),DIM=3) + SUM(bm_to_litter(:,:,:,:),DIM=3)

       IF(STO_CHECK_MASSBALANCE) THEN
        DO j=2,nvm
         IF (is_peat(j)) THEN
          !   WRITE (numout,*) 'QCJ check before establish, PFT',j,'biomass',SUM(biomass(:,j,:,icarbon),DIM=2),'bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2),'co2_to_bm',co2_to_bm(:,j)
          biomass_old(:,j,:,:)=biomass(:,j,:,:)
          bm_to_litter_old(:,j,:,:)=bm_to_litter(:,j,:,:)
          co2_to_bm_old(:,j)=co2_to_bm(:,j)
         ENDIF
        ENDDO
       ENDIF

       !! 11.1 Establish new plants
       CALL establish (npts, lalo, dt_days, PFTpresent, regenerate, &
            neighbours, resolution, need_adjacent, herbivores, &
            precip_lastyear, gdd0_lastyear, lm_lastyearmax, &
            cn_ind, lai, avail_tree, avail_grass, npp_longterm, &
            leaf_age, leaf_frac, &
            ind, biomass, age, everywhere, co2_to_bm, veget_cov_max, woodmass_ind, &
            mortality, bm_to_litter, &
!gmjc
            sla_calc, &
!end gmjc
            fpeat_map,wettile_dgvm) !!!qcj++ peatland

       !! 11.2 Calculate new crown area (and maximum vegetation cover)
       CALL crown (npts, PFTpresent, &
            ind, biomass, woodmass_ind, &
            veget_cov_max, cn_ind, height)

       IF(STO_CHECK_MASSBALANCE) THEN
        DO i=1, npts
         DO j=2,nvm
          IF (is_peat(j)) THEN
           flux_co2_to_bm(:,j)=co2_to_bm(:,j)-co2_to_bm_old(:,j)
           flux_bm_to_litter(:,j,:,:)=bm_to_litter(:,j,:,:)-bm_to_litter_old(:,j,:,:)
           !   WRITE (numout,*) 'QCJ check after establish, PFT',j,'biomass',SUM(biomass(:,j,:,icarbon),DIM=2),'bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2),'co2_to_bm',co2_to_bm(:,j)
           !   WRITE (numout,*) 'QCJ check after establish, PFT',j,'flux_bm_to_litter',SUM(flux_bm_to_litter(:,j,:,icarbon),DIM=2),'delta biomass',SUM(biomass(:,j,:,icarbon),DIM=2)-SUM(biomass_old(:,j,:,icarbon),DIM=2),'delta bm_to_litter',SUM(bm_to_litter(:,j,:,icarbon),DIM=2)-SUM(bm_to_litter_old(:,j,:,icarbon),DIM=2)
           IF (ABS(SUM(biomass(i,j,:,icarbon))+SUM(flux_bm_to_litter(i,j,:,icarbon))-flux_co2_to_bm(i,j)-SUM(biomass_old(i,j,:,icarbon))) .GT. min_stomate) THEN
             WRITE (numout,*) 'QCJ mass cons error, in establish'
           ENDIF
          ENDIF
         ENDDO
        ENDDO
       ENDIF

       IF(STO_CHECK_MASSBALANCE) THEN
          !DSG debug: ===================================================
          CALL  stomate_check_mass_values(npts,biomass(:,:,:,:), 'lpj: after establish')

          !DSG mass conservation ============================================
          mass_change(:,:,icarbon)     = -SUM(bm_to_litter(:,:,:,icarbon),DIM=3)

          mass_after =  SUM(biomass,DIM=3)
          CALL stomate_check_cons_mass(lalo, mass_before, &  ! mass before
                mass_after,                               &     ! mass after
                mass_change,                              &  ! net of fluxes
                'lpj: after establish')
       ENDIF
    ENDIF

!gmjc Grassland_management
    !
    ! 13 calculate grazing by animals or cutting for forage
    !
    IF (enable_grazing) THEN
        CALL main_grassland_management(&
           npts, lalo, neighbours, resolution, contfrac, &
           dt_days        , &
           day_of_year    , &
           t2m_daily      , &
           t2m_min_daily  , &
           t2m_14         , &
           tsurf_daily    , &
           snowfall_daily , &
           biomass        , &
           bm_to_litter   , &
           litter         , &
           litter_avail   , &
           litter_not_avail , &
           !spitfire
           fuel_1hr(:,:,:,icarbon), &
           fuel_10hr(:,:,:,icarbon), &
           fuel_100hr(:,:,:,icarbon), &
           fuel_1000hr(:,:,:,icarbon), &
           !end spitfire
           .TRUE.         , &
           EndOfYear    , & 
           when_growthinit_cut, nb_grazingdays, &
           lai,sla_calc,leaf_age,leaf_frac, &
           wshtotsum,sr_ugb,compt_ugb, &
           nb_ani,grazed_frac,import_yield,N_limfert, &
           moiavail_daily,tmc_topgrass_daily,fc_grazing, snowmass_daily, &
           after_snow, after_wet, wet1day, wet2day, &
           harvest_gm, ranimal_gm, ch4_pft_gm, cinput_gm, n2o_pft_gm)
    ENDIF
!end gmjc
  !! 12. Calculate final LAI and vegetation cover

!!!qcj++ peatland
   IF (update_peatfrac) THEN
       CALL lpj_cover_peat(npts,lalo, cn_ind, ind, biomass, &
            veget_cov_max_new, veget_cov_max, veget_cov_max_tmp, &
            litter, litter_avail, litter_not_avail, carbon, &
            fuel_1hr, fuel_10hr, fuel_100hr, fuel_1000hr, &
            turnover_daily, bm_to_litter, &
            co2_to_bm, co2_fire, resp_hetero, resp_maint, resp_growth, gpp_daily, &
            deepC_a, deepC_s, deepC_p, &
            dt_days, age, PFTpresent, senescence, when_growthinit,&
            everywhere, leaf_frac, lm_lastyearmax, npp_longterm,&
            carbon_save,deepC_a_save,deepC_s_save,deepC_p_save,delta_fsave,liqwt_max_lastyear)
      update_peatfrac = .FALSE.
      done_update_peatfrac = .TRUE.
   ELSE   

    ![chaoyue] veget_cov_max_tmp is used as veget_cov_max_old in cover SUBROUTINE
    CALL cover (npts, cn_ind, ind, biomass, &
         veget_cov_max, veget_cov_max_tmp, lai, &
         litter, litter_avail, litter_not_avail, carbon, & 
         fuel_1hr, fuel_10hr, fuel_100hr, fuel_1000hr, &
         turnover_daily, bm_to_litter, &
         co2_to_bm, co2_fire, resp_hetero, resp_maint, resp_growth, gpp_daily, &
         deepC_a, deepC_s,deepC_p, &
         fpeat_map,veget_cov_max_new,lalo,wettile_dgvm,date) !!!qcj++ peatland

   ENDIF

  !! 13. Update litter pools to account for harvest
 
    ! the whole litter stuff:
    !    litter update, lignin content, PFT parts, litter decay, 
    !    litter heterotrophic respiration, dead leaf soil cover.
    !    No vertical discretisation in the soil for litter decay.\n
    ! added by shilong for harvest
    IF(harvest_agri) THEN  !!! DO NOT activate harves_agri if ok_LAIdev
       CALL harvest(npts, dt_days, veget_cov_max, &
            bm_to_litter, turnover_daily, &
            harvest_above, harvest_bio)
    ENDIF


    IF (do_now_stomate_lcchange) THEN

    ! We initialize these variables to zero. If LUC modules are called, their
    ! values will be changed. Otherwise they stay as zero as is expected.
    ! JC comment, we can not initialize it every day, since value should be kept
    ! for the whole year for history write.
      convflux(:,:)         = zero
      cflux_prod10(:,:)     = zero
      cflux_prod100(:,:)    = zero
      ! Initialize LUC specific variables
      prod10(:,0,:)         = zero
      prod100(:,0,:)        = zero   

      IF (use_age_class) THEN

        ! Build the constants specific for gross land use change
        CALL stomate_gluc_constants_init

        WRITE(numout,*) 'Buiding indecies for age classes:'
        WRITE(numout,*) 'Buiding indecies for tress:'
        WRITE(numout,*) indall_tree
        WRITE(numout,*) indagec_tree
        WRITE(numout,*) indold_tree
        WRITE(numout,*) 'Buiding indecies for natural grass:'
        WRITE(numout,*) indall_grass
        WRITE(numout,*) indagec_grass
        WRITE(numout,*) indold_grass
        WRITE(numout,*) 'Buiding indecies for pasture:'
        WRITE(numout,*) indall_pasture
        WRITE(numout,*) indagec_pasture
        WRITE(numout,*) indold_pasture
        WRITE(numout,*) 'Buiding indecies for crop:'
        WRITE(numout,*) indall_crop
        WRITE(numout,*) indagec_crop
        WRITE(numout,*) indold_crop
        WRITE(numout,*) 'Buiding indecies for bioe1:'
        WRITE(numout,*) indall_bioe1
        WRITE(numout,*) indagec_bioe1
        WRITE(numout,*) indold_bioe1

        IF (gluc_use_harvest_biomass) THEN
          ! calculate the fuel fraction when input is the harvest biomass
          fuelfrac = harvest_biomass(:,3)
        ELSE
          fuelfrac = zero
        ENDIF

        IF (SingleAgeClass) THEN  
          CALL fharvest_SinAgeC (npts, dt_days, harvest_matrix,newvegfrac,   &
             fuelfrac, &
             glccZero,glccZero,glccZero,&
             def_fuel_1hr_remain, def_fuel_10hr_remain,        &
             def_fuel_100hr_remain, def_fuel_1000hr_remain,    &
             deforest_litter_remain, deforest_biomass_remain,  &
             convflux, cflux_prod10, cflux_prod100,                  &
             glccReal, IncreDeficit, glcc_pft, glcc_pftmtc,          &
             veget_cov_max, prod10, prod100, flux10, flux100,            &
             PFTpresent, senescence, moiavail_month, moiavail_week,  &
             gpp_week, ngd_minus5, resp_maint, resp_growth,          & 
             resp_hetero, npp_daily, when_growthinit, npp_longterm,  &
             ind, lm_lastyearmax, everywhere, age,                   &
             co2_to_bm, gpp_daily, co2_fire,                         &
             time_hum_min, gdd_midwinter, gdd_from_growthinit,       &
             gdd_m5_dormance, ncd_dormance,                          &
             lignin_struc, carbon, leaf_frac,                        &
             deepC_a, deepC_s, deepC_p,                              &
             leaf_age, bm_to_litter, biomass, litter,                &
             fuel_1hr, fuel_10hr, fuel_100hr, fuel_1000hr)

          CALL glcc_SinAgeC (npts, dt_days, glccZero,newvegfrac,   &
             glccSecondShift,glccPrimaryShift,glccNetLCC,&
             def_fuel_1hr_remain, def_fuel_10hr_remain,        &
             def_fuel_100hr_remain, def_fuel_1000hr_remain,    &
             deforest_litter_remain, deforest_biomass_remain,  &
             convflux, cflux_prod10, cflux_prod100,                  &
             glccReal, IncreDeficit, glcc_pft, glcc_pftmtc,          &
             veget_cov_max, prod10, prod100, flux10, flux100,            &
             PFTpresent, senescence, moiavail_month, moiavail_week,  &
             gpp_week, ngd_minus5, resp_maint, resp_growth,          & 
             resp_hetero, npp_daily, when_growthinit, npp_longterm,  &
             ind, lm_lastyearmax, everywhere, age,                   &
             co2_to_bm, gpp_daily, co2_fire,                         &
             time_hum_min, gdd_midwinter, gdd_from_growthinit,       &
             gdd_m5_dormance, ncd_dormance,                          &
             lignin_struc, carbon, leaf_frac,                        &
             deepC_a, deepC_s, deepC_p,                              &
             leaf_age, bm_to_litter, biomass, litter,                &
             fuel_1hr, fuel_10hr, fuel_100hr, fuel_1000hr)

        ! Multiple age class + wood harvest
        ELSE

          IF (allow_forestry_harvest) THEN

            IF (gluc_use_harvest_biomass) THEN
              CALL Get_harvest_matrix (npts,veget_cov_max,newvegfrac,   &
                            harvest_matrix,                             &
                            biomass, harvest_biomass, area_land_m2,     &
                            glcc_pft,glcc_pftmtc,                       &
                            Deficit_pf2yf_final, Deficit_sf2yf_final,   &
                            pf2yf_compen_sf2yf, sf2yf_compen_pf2yf)
            ENDIF

            CALL prepare_balance_check(outflux_sta,influx_sta,pool_sta,    &
                 veget_cov_max,                                            &
                 co2_to_bm,gpp_daily,npp_daily,                            &
                 biomass,litter,carbon,prod10,prod100,                     &
                 bm_to_litter,resp_maint,resp_growth,resp_hetero,          &
                 convflux,cflux_prod10,cflux_prod100,co2_fire)

            CALL fharvest_MulAgeC (npts, dt_days, harvest_matrix,newvegfrac,   &
               fuelfrac,&
               def_fuel_1hr_remain, def_fuel_10hr_remain,              &
               def_fuel_100hr_remain, def_fuel_1000hr_remain,          &
               deforest_litter_remain, deforest_biomass_remain,        &
               convflux,                   &
               glcc_pft, glcc_pftmtc,          &
               veget_cov_max, prod10, prod100,         &
               PFTpresent, senescence, moiavail_month, moiavail_week,  &
               gpp_week, ngd_minus5, resp_maint, resp_growth,          & 
               resp_hetero, npp_daily, when_growthinit, npp_longterm,  &
               ind, lm_lastyearmax, everywhere, age,                   &
               co2_to_bm, gpp_daily, co2_fire,                         &
               time_hum_min, gdd_midwinter, gdd_from_growthinit,       &
               gdd_m5_dormance, ncd_dormance,                          &
               lignin_struc, carbon, leaf_frac,                        &
               deepC_a, deepC_s, deepC_p,                              &
               leaf_age, bm_to_litter, biomass, litter,                &
               fuel_1hr, fuel_10hr, fuel_100hr, fuel_1000hr)

            CALL prepare_balance_check(outflux_end,influx_end,pool_end,    &
                 veget_cov_max,                                            &
                 co2_to_bm,gpp_daily,npp_daily,                            &
                 biomass,litter,carbon,prod10,prod100,                     &
                 bm_to_litter,resp_maint,resp_growth,resp_hetero,          &
                 convflux,cflux_prod10,cflux_prod100,co2_fire)

            CALL luc_balance_check(outflux_sta,influx_sta,pool_sta,        &
                 outflux_end,influx_end,pool_end,                          &
                 npts,lalo,'wood harvest')

          ENDIF

          CALL prepare_balance_check(outflux_sta,influx_sta,pool_sta,    &
               veget_cov_max,                                            &
               co2_to_bm,gpp_daily,npp_daily,                            &
               biomass,litter,carbon,prod10,prod100,                     &
               bm_to_litter,resp_maint,resp_growth,resp_hetero,          &
               convflux,cflux_prod10,cflux_prod100,co2_fire)

          CALL glcc_MulAgeC (npts, dt_days, newvegfrac,              &
             glccSecondShift,glccPrimaryShift,glccNetLCC,            &
             def_fuel_1hr_remain, def_fuel_10hr_remain,              &
             def_fuel_100hr_remain, def_fuel_1000hr_remain,          &
             deforest_litter_remain, deforest_biomass_remain,        &
             convflux,                  &
             glccReal, IncreDeficit, glcc_pft, glcc_pftmtc,          &
             veget_cov_max, prod10, prod100,         &
             PFTpresent, senescence, moiavail_month, moiavail_week,  &
             gpp_week, ngd_minus5, resp_maint, resp_growth,          & 
             resp_hetero, npp_daily, when_growthinit, npp_longterm,  &
             ind, lm_lastyearmax, everywhere, age,                   &
             co2_to_bm, gpp_daily, co2_fire,                         &
             time_hum_min, gdd_midwinter, gdd_from_growthinit,       &
             gdd_m5_dormance, ncd_dormance,                          &
             lignin_struc, carbon, leaf_frac,                        &
             deepC_a, deepC_s, deepC_p,                              &
             leaf_age, bm_to_litter, biomass, litter,                &
             fuel_1hr, fuel_10hr, fuel_100hr, fuel_1000hr)

          CALL prepare_balance_check(outflux_end,influx_end,pool_end,    &
               veget_cov_max,                                            &
               co2_to_bm,gpp_daily,npp_daily,                            &
               biomass,litter,carbon,prod10,prod100,                     &
               bm_to_litter,resp_maint,resp_growth,resp_hetero,          &
               convflux,cflux_prod10,cflux_prod100,co2_fire)

          CALL luc_balance_check(outflux_sta,influx_sta,pool_sta,        &
               outflux_end,influx_end,pool_end,                          &
               npts,lalo,'gross land cover change')
           
          IF (gluc_allow_trans_bioe) THEN

            CALL prepare_balance_check(outflux_sta,influx_sta,pool_sta,    &
                 veget_cov_max,                                            &
                 co2_to_bm,gpp_daily,npp_daily,                            &
                 biomass,litter,carbon,prod10,prod100,                     &
                 bm_to_litter,resp_maint,resp_growth,resp_hetero,          &
                 convflux,cflux_prod10,cflux_prod100,co2_fire)

            CALL glcc_bioe1 (npts, dt_days, newvegfrac,              &
             glccSecondShift*0,                                      &
             def_fuel_1hr_remain, def_fuel_10hr_remain,              &
             def_fuel_100hr_remain, def_fuel_1000hr_remain,          &
             deforest_litter_remain, deforest_biomass_remain,        &
             convflux, cflux_prod10, cflux_prod100,                  &
             glcc_pft, glcc_pftmtc,                                  &
             veget_cov_max, prod10, prod100, flux10, flux100,        &
             PFTpresent, senescence, moiavail_month, moiavail_week,  &
             gpp_week, ngd_minus5, resp_maint, resp_growth,          & 
             resp_hetero, npp_daily, when_growthinit, npp_longterm,  &
             ind, lm_lastyearmax, everywhere, age,                   &
             co2_to_bm, gpp_daily, co2_fire,                         &
             time_hum_min, gdd_midwinter, gdd_from_growthinit,       &
             gdd_m5_dormance, ncd_dormance,                          &
             lignin_struc, carbon, leaf_frac,                        &
             deepC_a, deepC_s, deepC_p,                              &
             leaf_age, bm_to_litter, biomass, litter,                &
             fuel_1hr, fuel_10hr, fuel_100hr, fuel_1000hr)

            CALL prepare_balance_check(outflux_end,influx_end,pool_end,    &
                 veget_cov_max,                                            &
                 co2_to_bm,gpp_daily,npp_daily,                            &
                 biomass,litter,carbon,prod10,prod100,                     &
                 bm_to_litter,resp_maint,resp_growth,resp_hetero,          &
                 convflux,cflux_prod10,cflux_prod100,co2_fire)

            CALL luc_balance_check(outflux_sta,influx_sta,pool_sta,        &
                 outflux_end,influx_end,pool_end,                          &
                 npts,lalo,'glcc with bioenergy PFT')

          ENDIF

        ENDIF !(SingleAgeClass)

        ! clear the constants specific for gross land use change
        CALL stomate_gluc_constants_init_clear

      ELSE ! .NOT. use_age_class
        IF (allow_deforest_fire) THEN
              ![chaoyue] veget_cov_max is used as the old veget_cov_max in lcchange_deffire
              ! veget_cov_max_new is used as the new veget_cov_max in lcchange_deffire
              CALL lcchange_deffire (npts, dt_days, veget_cov_max, veget_cov_max_new, &
                   biomass, ind, age, PFTpresent, senescence, when_growthinit, everywhere, &        
                   co2_to_bm, bm_to_litter, turnover_daily, bm_sapl, cn_ind, &
                   flux10,flux100, &
                   prod10,prod100,&
                   convflux,&
                   cflux_prod10,cflux_prod100, leaf_frac,&
                   npp_longterm, lm_lastyearmax, litter, litter_avail, litter_not_avail, &
                   carbon,&
                   deepC_a, deepC_s, deepC_p,&
                   fuel_1hr,fuel_10hr,fuel_100hr,fuel_1000hr,&
                   lcc,bafrac_deforest_accu,emideforest_litter_accu,emideforest_biomass_accu,&
                   deflitsup_total,defbiosup_total)
        ELSE
              ![chaoyue] veget_cov_max is used as veget_cov_max_old in lcchange_main
              ! veget_cov_max_new is used as the new veget_cov_max in lcchange_main
!!!qcj++ peatland
              IF (agri_peat) THEN
                 CALL lcchange_main_agripeat (npts, dt_days, veget_cov_max, veget_cov_max_new, &
                   biomass, ind, age, PFTpresent, senescence, when_growthinit, everywhere, &
                   co2_to_bm, bm_to_litter, turnover_daily, bm_sapl, cn_ind,&
                   flux10,flux100, &
                   prod10,prod100,convflux, &
                   cflux_prod10,cflux_prod100,leaf_frac,&
                   npp_longterm, lm_lastyearmax, litter, litter_avail, litter_not_avail, &
                   carbon, &
                   deepC_a, deepC_s, deepC_p,&
                   fuel_1hr,fuel_10hr,fuel_100hr,fuel_1000hr, &
                   veget_cov_max_adjusted,lalo,carbon_save,deepC_a_save,deepC_s_save,deepC_p_save,delta_fsave)
!,biomass_remove)
              ELSE
                 IF (.NOT. dyn_peat) THEN
              !WRITE (numout,*) 'QCJ check before lcchange_main'  
              CALL lcchange_main (npts, dt_days, veget_cov_max, veget_cov_max_new, &
                   biomass, ind, age, PFTpresent, senescence, when_growthinit, everywhere, &
                   co2_to_bm, bm_to_litter, turnover_daily, bm_sapl, cn_ind,&
                   flux10,flux100, &
                   prod10,prod100,convflux, &
                   cflux_prod10,cflux_prod100,leaf_frac,&
                   npp_longterm, lm_lastyearmax, litter, litter_avail, litter_not_avail, &
                   carbon, &
                   deepC_a, deepC_s, deepC_p,&
                   fuel_1hr,fuel_10hr,fuel_100hr,fuel_1000hr, &
!!!qcj++ peatland
                   veget_cov_max_new_peatdgvm,lalo )
               !WRITE (numout,*) 'QCJ check after lcchange_main'
                   IF (ok_dgvm_peat) THEN
                      veget_cov_max(:,:) = veget_cov_max_new_peatdgvm(:,:)
                   ENDIF
                 ENDIF
              ENDIF
        ENDIF
      ENDIF ! (use_age_class)

      IF (use_age_class) THEN
        !! Update the 10-year-turnover pool content following flux emission
        !!     (linear decay (10%) of the initial carbon input)
        DO  l = 0, 8
          m = 10 - l
          cflux_prod10(:,:) =  cflux_prod10(:,:) + flux10(:,m,:)
          prod10(:,m,:)     =  prod10(:,m-1,:)   - flux10(:,m-1,:)
          flux10(:,m,:)     =  flux10(:,m-1,:)
        ENDDO
        
        cflux_prod10(:,:) = cflux_prod10(:,:) + flux10(:,1,:) 
        flux10(:,1,:)     = 0.1 * prod10(:,0,:)
        prod10(:,1,:)     = prod10(:,0,:)
        
        !! Update the 100-year-turnover pool content following flux emission\n
        DO l = 0, 98
           m = 100 - l
           cflux_prod100(:,:)  =  cflux_prod100(:,:) + flux100(:,m,:)
           prod100(:,m,:)      =  prod100(:,m-1,:)   - flux100(:,m-1,:)
           flux100(:,m,:)      =  flux100(:,m-1,:)
        ENDDO
        
        cflux_prod100(:,:)  = cflux_prod100(:,:) + flux100(:,1,:) 
        flux100(:,1,:)      = 0.01 * prod100(:,0,:)
        prod100(:,1,:)      = prod100(:,0,:)
        prod10(:,0,:)        = zero
        prod100(:,0,:)       = zero 
      ENDIF ! (use_age_class)
      ! set the flag
      do_now_stomate_lcchange=.FALSE.

      ! Set the flag done_stomate_lcchange to be used in the end of sechiba_main to update the fractions.
      done_stomate_lcchange=.TRUE.

      CALL stomate_check_mass_values(npts,biomass(:,:,:,:), 'lpj: after land use change')

    ENDIF ! do_now_stomate_lcchange

    ! Update the status of age classes
    IF(EndOfYear) THEN

      ! start the mass balance check
      mass_before_2rd = zero
      mass_change_2rd = zero
      mass_after_2rd = zero
      mass_balance_2rd = zero
      pool_sta = zero
      influx_sta = zero
      outflux_sta = zero
      
      ind_biomass = 1
      ind_litter = 2
      ind_soil = 3
      pool_sta(:,ind_biomass) = SUM(SUM(biomass(:,:,:,icarbon),DIM=3) * veget_cov_max,DIM=2)
      pool_sta(:,ind_litter) = SUM(SUM(SUM(litter(:,:,:,:,icarbon),DIM=4),DIM=2) * veget_cov_max,DIM=2)
      pool_sta(:,ind_soil) = SUM(SUM(carbon(:,:,:),DIM=2) * veget_cov_max,DIM=2)
      mass_before_2rd(:,icarbon) = SUM(pool_sta(:,:),DIM=2)

      ind_co2tobm = 1
      ind_gpp = 2
      ind_npp = 3
      influx_sta(:,ind_co2tobm) = SUM(co2_to_bm*veget_cov_max,DIM=2)
      influx_sta(:,ind_gpp) = SUM(gpp_daily*veget_cov_max,DIM=2)
      influx_sta(:,ind_npp) = SUM(npp_daily*veget_cov_max,DIM=2)

      ind_bm2lit = 1
      ind_respm = 2
      ind_respg = 3
      ind_resph = 4
      outflux_sta(:,ind_bm2lit) = SUM(SUM(bm_to_litter(:,:,:,icarbon),DIM=3) * veget_cov_max,DIM=2)
      outflux_sta(:,ind_respm) = SUM(resp_maint*veget_cov_max,DIM=2)
      outflux_sta(:,ind_respg) = SUM(resp_growth*veget_cov_max,DIM=2)
      outflux_sta(:,ind_resph) = SUM(resp_hetero*veget_cov_max,DIM=2)


      IF (use_age_class) THEN
        CALL age_class_distr(npts, lalo, resolution, bound_spa, &
             biomass, veget_cov_max, ind, &
             lm_lastyearmax, leaf_frac, co2_to_bm, &
             fuel_1hr, fuel_10hr, fuel_100hr, fuel_1000hr, &
             everywhere, litter, carbon, lignin_struc, &
             deepC_a, deepC_s, deepC_p, &
             bm_to_litter, PFTpresent, when_growthinit,&
             senescence, npp_longterm, gpp_daily, leaf_age, age, &
             gdd_from_growthinit, gdd_midwinter, time_hum_min, hum_min_dormance, &
             gdd_m5_dormance, &
             ncd_dormance, moiavail_month, moiavail_week, ngd_minus5, &
             gpp_week, resp_maint, resp_growth, npp_daily, resp_hetero)
      ENDIF

      pool_end = zero
      influx_end = zero
      outflux_end = zero
      
      ind_biomass = 1
      ind_litter = 2
      ind_soil = 3
      pool_end(:,ind_biomass) = SUM(SUM(biomass(:,:,:,icarbon),DIM=3) * veget_cov_max,DIM=2)
      pool_end(:,ind_litter) = SUM(SUM(SUM(litter(:,:,:,:,icarbon),DIM=4),DIM=2) * veget_cov_max,DIM=2)
      pool_end(:,ind_soil) = SUM(SUM(carbon(:,:,:),DIM=2) * veget_cov_max,DIM=2)

      ind_co2tobm = 1
      ind_gpp = 2
      ind_npp = 3
      influx_end(:,ind_co2tobm) = SUM(co2_to_bm*veget_cov_max,DIM=2)
      influx_end(:,ind_gpp) = SUM(gpp_daily*veget_cov_max,DIM=2)
      influx_end(:,ind_npp) = SUM(npp_daily*veget_cov_max,DIM=2)


      ind_bm2lit = 1
      ind_respm = 2
      ind_respg = 3
      ind_resph = 4
      outflux_end(:,ind_bm2lit) = SUM(SUM(bm_to_litter(:,:,:,icarbon),DIM=3) * veget_cov_max,DIM=2)
      outflux_end(:,ind_respm) = SUM(resp_maint*veget_cov_max,DIM=2)
      outflux_end(:,ind_respg) = SUM(resp_growth*veget_cov_max,DIM=2)
      outflux_end(:,ind_resph) = SUM(resp_hetero*veget_cov_max,DIM=2)

      CALL stomate_check_mass_values(npts,biomass(:,:,:,:), 'lpj: after age_class_distr')

      ! mass_change_2rd is the mass increase
      mass_change_2rd(:,icarbon) = SUM(influx_end(:,:),DIM=2) - SUM(influx_sta(:,:),DIM=2) + &
                               + SUM(outflux_sta(:,:),DIM=2) - SUM(outflux_end(:,:),DIM=2)

      mass_after_2rd(:,icarbon) = SUM(pool_end(:,:),dim=2)
      mass_balance_2rd(:,icarbon) = mass_before_2rd(:,icarbon) - mass_after_2rd(:,icarbon) &
                               + mass_change_2rd(:,icarbon)
       
      DO ipts = 1,npts
        IF (ABS(mass_balance_2rd(ipts,icarbon)) .GE. 1e-3) THEN
          WRITE (numout,*) 'FATAL Error'
          WRITE (numout,*) 'Mass conservation failed after age_class_distr'
          WRITE (numout,*) ' limit: '       , 1e-3
          WRITE (numout,*) ' Mismatch :'    , mass_balance_2rd(ipts,icarbon) 
          WRITE (numout,*) ' Coordinates :' , lalo(ipts,:) 
          WRITE (numout,*) ' gridpoint: '   , ipts , ' of ngrids: ',npts
          STOP
        ENDIF
      ENDDO
      WRITE (numout,*) 'mass balance check successful after age_class_distr'

    ENDIF


    !! 16. Recalculate lai
    !!      This should be done whenever biomass is modified
    CALL setlai(biomass, sla_calc, slai)
    CALL setlai(biomass, sla_calc, lai)

    !! 17. Calculate vcmax 
    CALL vmax (npts, dt_days, &
         leaf_age, leaf_frac, &
         vcmax, &!)
         N_limfert, moiavail_month) !!!qcj++ peatland, following Druel et al. 2017 GMD paper,add humrel_month for dessication impact

    !MM deplacement pour initialisation correcte des grandeurs cumulees :
    cflux_prod_total(:) = convflux(:,iwplcc) + cflux_prod10(:,iwplcc) + cflux_prod100(:,iwplcc)
    prod10_total(:)=SUM(prod10(:,:,iwplcc),dim=2)
    prod100_total(:)=SUM(prod100(:,:,iwplcc),dim=2)

    cflux_prod_harvest_total = zero
    prod10_harvest_total = zero
    prod100_harvest_total = zero
    
  !! 18. Total heterotrophic respiration

    tot_soil_carb(:,:) = zero
    tot_litter_carb(:,:) = zero
    DO j=2,nvm

       tot_litter_carb(:,j) = tot_litter_carb(:,j) + (litter(:,istructural,j,iabove,icarbon) + &
            &          litter(:,imetabolic,j,iabove,icarbon) + &
            &          litter(:,istructural,j,ibelow,icarbon) + litter(:,imetabolic,j,ibelow,icarbon))

       tot_soil_carb(:,j) = tot_soil_carb(:,j) + (carbon(:,iactive,j) + &
            &          carbon(:,islow,j)+  carbon(:,ipassive,j))

    ENDDO
    tot_litter_soil_carb(:,:) = tot_litter_carb(:,:) + tot_soil_carb(:,:)

!!$     DO k = 1, nelements ! Loop over # elements
!!$        tot_live_biomass(:,:,k) = biomass(:,:,ileaf,k) + biomass(:,:,isapabove,k) + biomass(:,:,isapbelow,k) +&
!!$             &                    biomass(:,:,iheartabove,k) + biomass(:,:,iheartbelow,k) + &
!!$             &                    biomass(:,:,iroot,k) + biomass(:,:,ifruit,k) + biomass(:,:,icarbres,k)
!!$    END DO ! Loop over # elements

    tot_live_biomass(:,:,:) = biomass(:,:,ileaf,:) + biomass(:,:,isapabove,:) + biomass(:,:,isapbelow,:) +&
             &                    biomass(:,:,iheartabove,:) + biomass(:,:,iheartbelow,:) + &
             &                    biomass(:,:,iroot,:) + biomass(:,:,ifruit,:) + biomass(:,:,icarbres,:)


    tot_turnover(:,:,:) = turnover_daily(:,:,ileaf,:) + turnover_daily(:,:,isapabove,:) + &
         &         turnover_daily(:,:,isapbelow,:) + turnover_daily(:,:,iheartabove,:) + &
         &         turnover_daily(:,:,iheartbelow,:) + turnover_daily(:,:,iroot,:) + &
         &         turnover_daily(:,:,ifruit,:) + turnover_daily(:,:,icarbres,:)

    tot_bm_to_litter(:,:,:) = bm_to_litter(:,:,ileaf,:) + bm_to_litter(:,:,isapabove,:) +&
         &             bm_to_litter(:,:,isapbelow,:) + bm_to_litter(:,:,iheartbelow,:) +&
         &             bm_to_litter(:,:,iheartabove,:) + bm_to_litter(:,:,iroot,:) + &
         &             bm_to_litter(:,:,ifruit,:) + bm_to_litter(:,:,icarbres,:)

    carb_mass_variation(:)=-carb_mass_total(:)
    carb_mass_total(:)=SUM((tot_live_biomass(:,:,icarbon)+tot_litter_carb+tot_soil_carb)*veget_cov_max,dim=2) + &
         &                 (prod10_total + prod100_total) +  (prod10_harvest_total + prod100_harvest_total)
    carb_mass_variation(:)=carb_mass_total(:)+carb_mass_variation(:)
    
  !! 17. Write history

    CALL xios_orchidee_send_field("RESOLUTION_X",resolution(:,1))
    CALL xios_orchidee_send_field("RESOLUTION_Y",resolution(:,2))
    CALL xios_orchidee_send_field("CONTFRAC_STOMATE",contfrac(:))
    CALL xios_orchidee_send_field("T2M_MONTH",t2m_month)
    CALL xios_orchidee_send_field("T2M_WEEK",t2m_week)
    CALL xios_orchidee_send_field("HET_RESP",resp_hetero(:,:))
    CALL xios_orchidee_send_field("CO2_FIRE",co2_fire)
    CALL xios_orchidee_send_field("CO2_TAKEN",co2_to_bm)
    CALL xios_orchidee_send_field("LAI",lai)
    CALL xios_orchidee_send_field("VEGET_COV_MAX",veget_cov_max)
    CALL xios_orchidee_send_field("NPP_STOMATE",npp_daily)
    CALL xios_orchidee_send_field("GPP",gpp_daily)
    CALL xios_orchidee_send_field("IND",ind)
    CALL xios_orchidee_send_field("CN_IND",cn_ind)
    CALL xios_orchidee_send_field("WOODMASS_IND",woodmass_ind)
    CALL xios_orchidee_send_field("TOTAL_M",tot_live_biomass(:,:,icarbon))
    CALL xios_orchidee_send_field("MOISTRESS",moiavail_week)
    CALL xios_orchidee_send_field("LEAF_M",biomass(:,:,ileaf,icarbon))
    CALL xios_orchidee_send_field("SAP_M_AB",biomass(:,:,isapabove,icarbon))
    CALL xios_orchidee_send_field("SAP_M_BE",biomass(:,:,isapbelow,icarbon))
    CALL xios_orchidee_send_field("HEART_M_AB",biomass(:,:,iheartabove,icarbon))
    CALL xios_orchidee_send_field("HEART_M_BE",biomass(:,:,iheartbelow,icarbon))
    CALL xios_orchidee_send_field("ROOT_M",biomass(:,:,iroot,icarbon))
    CALL xios_orchidee_send_field("FRUIT_M",biomass(:,:,ifruit,icarbon))
    ! HERE WE ADD A OUTPUT VARIABLE, NAMED CROPYIELD 
    CALL xios_orchidee_send_field("CROPYIELD",biomass(:,:,ifruit,icarbon))
    DO i=1, npts
    DO j=2,nvm
       IF (is_peat(j)) THEN
    !      WRITE (numout,*) 'QC check output,pft',j
    !      WRITE (numout,*) 'QC check biomass,',biomass(i,j,ileaf,icarbon),biomass(i,j,isapabove,icarbon),biomass(i,j,isapbelow,icarbon),biomass(i,j,iheartabove,icarbon),biomass(i,j,iheartbelow,icarbon),biomass(i,j,iroot,icarbon),biomass(i,j,ifruit,icarbon),biomass(i,j,icarbres,icarbon)
    !      WRITE (numout,*) 'QC check litter,',litter(i,istructural,j,iabove,icarbon),litter(i,imetabolic,j,iabove,icarbon),litter(i,istructural,j,ibelow,icarbon),litter(i,imetabolic,j,ibelow,icarbon)
    !      WRITE (numout,*) 'QC check soc,',carbon(i,iactive,j),carbon(i,islow,j),carbon(i,ipassive,j)
       ENDIF
    ENDDO
    ENDDO

    CALL xios_orchidee_send_field("BIOMYIELD",tot_live_biomass(:,:,icarbon))
    ! END DEFINE
    CALL xios_orchidee_send_field("RESERVE_M",biomass(:,:,icarbres,icarbon))
!    CALL xios_orchidee_send_field("TOTAL_TURN",tot_turnover)
    CALL xios_orchidee_send_field("LEAF_TURN",turnover_daily(:,:,ileaf,icarbon))
    CALL xios_orchidee_send_field("MAINT_RESP",resp_maint)
    CALL xios_orchidee_send_field("GROWTH_RESP",resp_growth)
    CALL xios_orchidee_send_field("SAP_AB_TURN",turnover_daily(:,:,isapabove,icarbon))
    CALL xios_orchidee_send_field("ROOT_TURN",turnover_daily(:,:,iroot,icarbon))
    CALL xios_orchidee_send_field("FRUIT_TURN",turnover_daily(:,:,ifruit,icarbon))
    CALL xios_orchidee_send_field("TOTAL_BM_LITTER",tot_bm_to_litter(:,:,icarbon))
    CALL xios_orchidee_send_field("LEAF_BM_LITTER",bm_to_litter(:,:,ileaf,icarbon))
    CALL xios_orchidee_send_field("SAP_AB_BM_LITTER",bm_to_litter(:,:,isapabove,icarbon))
    CALL xios_orchidee_send_field("SAP_BE_BM_LITTER",bm_to_litter(:,:,isapbelow,icarbon))
    CALL xios_orchidee_send_field("HEART_AB_BM_LITTER",bm_to_litter(:,:,iheartabove,icarbon))
    CALL xios_orchidee_send_field("HEART_BE_BM_LITTER",bm_to_litter(:,:,iheartbelow,icarbon))
    CALL xios_orchidee_send_field("ROOT_BM_LITTER",bm_to_litter(:,:,iroot,icarbon))
    CALL xios_orchidee_send_field("FRUIT_BM_LITTER",bm_to_litter(:,:,ifruit,icarbon))
    CALL xios_orchidee_send_field("RESERVE_BM_LITTER",bm_to_litter(:,:,icarbres,icarbon))
    CALL xios_orchidee_send_field("LITTER_STR_AB",litter(:,istructural,:,iabove,icarbon))
    CALL xios_orchidee_send_field("LITTER_MET_AB",litter(:,imetabolic,:,iabove,icarbon))
    CALL xios_orchidee_send_field("LITTER_STR_BE",litter(:,istructural,:,ibelow,icarbon))
    CALL xios_orchidee_send_field("LITTER_MET_BE",litter(:,imetabolic,:,ibelow,icarbon))
    CALL xios_orchidee_send_field("DEADLEAF_COVER",deadleaf_cover)
    CALL xios_orchidee_send_field("TOTAL_SOIL_CARB",tot_litter_soil_carb)
    CALL xios_orchidee_send_field("CARBON_ACTIVE",carbon(:,iactive,:))
    CALL xios_orchidee_send_field("CARBON_SLOW",carbon(:,islow,:))
    CALL xios_orchidee_send_field("CARBON_PASSIVE",carbon(:,ipassive,:))
    CALL xios_orchidee_send_field("LITTERHUM",litterhum_daily)
    CALL xios_orchidee_send_field("TURNOVER_TIME",turnover_time)
!    CALL xios_orchidee_send_field("PROD10",prod10)
!    CALL xios_orchidee_send_field("FLUX10",flux10)
!    CALL xios_orchidee_send_field("PROD100",prod100)
!    CALL xios_orchidee_send_field("FLUX100",flux100)
    CALL xios_orchidee_send_field("CONVFLUX",convflux)
    CALL xios_orchidee_send_field("CFLUX_PROD10",cflux_prod10)
    CALL xios_orchidee_send_field("CFLUX_PROD100",cflux_prod100)
    CALL xios_orchidee_send_field("HARVEST_ABOVE",harvest_above)
!!!qcj++ peatland
    CALL xios_orchidee_send_field("HARVEST_BIO",harvest_bio)
    CALL xios_orchidee_send_field("VCMAX",vcmax)
    CALL xios_orchidee_send_field("AGE",age)
    CALL xios_orchidee_send_field("HEIGHT",height)
    CALL xios_orchidee_send_field("FIREINDEX",fireindex(:,:))
!gmjc
    CALL xios_orchidee_send_field("LITTER_STR_AVAIL",litter_avail(:,istructural,:))
    CALL xios_orchidee_send_field("LITTER_MET_AVAIL",litter_avail(:,imetabolic,:))
    CALL xios_orchidee_send_field("LITTER_STR_NAVAIL",litter_not_avail(:,istructural,:))
    CALL xios_orchidee_send_field("LITTER_MET_NAVAIL",litter_not_avail(:,imetabolic,:))
    CALL xios_orchidee_send_field("LITTER_STR_AVAILF",litter_avail_frac(:,istructural,:))
    CALL xios_orchidee_send_field("LITTER_MET_AVAILF",litter_avail_frac(:,imetabolic,:))
    IF (ANY(ok_LAIdev)) CALL xios_orchidee_send_field("N_LIMFERT",N_limfert)
!!!qcj++ peatland
    IF (ok_peat) THEN
       CALL xios_orchidee_send_field("TCARBON_ACRO", tcarbon_acro)
       CALL xios_orchidee_send_field("TCARBON_CATO", tcarbon_cato)
       CALL xios_orchidee_send_field("CARBON_ACRO",carbon_acro(:,:))
       CALL xios_orchidee_send_field("CARBON_CATO",carbon_cato(:,:))
       CALL xios_orchidee_send_field("HEIGHT_ACRO",height_acro)
       CALL xios_orchidee_send_field("HEIGHT_CATO",height_cato)
       CALL xios_orchidee_send_field("RESP_ACRO_oxic",resp_acro_oxic_d(:,:))
       CALL xios_orchidee_send_field("RESP_ACRO_anoxic",resp_acro_anoxic_d(:,:))
       CALL xios_orchidee_send_field("RESP_CATO",resp_cato_d(:,:))
       CALL xios_orchidee_send_field("LITTER_to_ACRO",litter_to_acro_d(:,:))
       CALL xios_orchidee_send_field("ACRO_to_CATO",acro_to_cato_d)
    ENDIF
    CALL xios_orchidee_send_field("deepC_a_save",deepC_a_save)
    CALL xios_orchidee_send_field("deepC_s_save",deepC_s_save)
    CALL xios_orchidee_send_field("deepC_p_save",deepC_p_save)
    CALL xios_orchidee_send_field("delta_fsave",delta_fsave)
    CALL xios_orchidee_send_field("carbon_a_save",carbon_save(:,iactive,:))
    CALL xios_orchidee_send_field("carbon_s_save",carbon_save(:,islow,:))
    CALL xios_orchidee_send_field("carbon_p_save",carbon_save(:,ipassive,:))
    CALL xios_orchidee_send_field('HUMREL_MONTH',moiavail_month)
    CALL xios_orchidee_send_field('NPP0_CUMUL',npp0_cumul)
    CALL xios_orchidee_send_field('fpeat_map',fpeat_map)
    CALL xios_orchidee_send_field("WTP_MONTH",wtp_month)
    CALL xios_orchidee_send_field("WTPMAX_MONTH",wtpmax_month)

 !   CALL xios_orchidee_send_field("biomass_remove_leaf",biomass_remove(:,:,ileaf,icarbon))
 !   CALL xios_orchidee_send_field("biomass_remove_sapabove",biomass_remove(:,:,isapabove,icarbon))
 !   CALL xios_orchidee_send_field("biomass_remove_sapbelow",biomass_remove(:,:,isapbelow,icarbon))
 !   CALL xios_orchidee_send_field("biomass_remove_heartabove",biomass_remove(:,:,iheartabove,icarbon))
 !   CALL xios_orchidee_send_field("biomass_remove_heartbelow",biomass_remove(:,:,iheartbelow,icarbon))
 !   CALL xios_orchidee_send_field("biomass_remove_root",biomass_remove(:,:,iroot,icarbon))
 !   CALL xios_orchidee_send_field("biomass_remove_fruit",biomass_remove(:,:,ifruit,icarbon))
 !   CALL xios_orchidee_send_field("biomass_remove_carbres",biomass_remove(:,:,icarbres,icarbon))
  
    !calculate grassland co2 fluxes
    DO j=2,nvm
      IF ((.NOT. is_tree(j)) .AND. natural(j)) THEN
        veget_cov_max_gm(:,j) = veget_cov_max(:,j)
      ENDIF
    END DO ! nvm
    veget_exist_gm(:) = SUM(veget_cov_max_gm,dim=2)
    WHERE (veget_exist_gm(:) .GT. 0.0)
      co2_gm(:) = SUM((gpp_daily-(resp_maint+resp_growth+resp_hetero)-co2_fire &
                    -harvest_gm-ranimal_gm-ch4_pft_gm+cinput_gm) &
                    *veget_cov_max_gm,dim=2)/veget_exist_gm
      ch4_gm(:) = SUM(ch4_pft_gm*veget_cov_max_gm,dim=2)/veget_exist_gm
      n2o_gm(:) = SUM(n2o_pft_gm*veget_cov_max_gm,dim=2)/veget_exist_gm
    ELSEWHERE
      co2_gm(:) = zero
      ch4_gm(:) = zero
      n2o_gm(:) = zero
    ENDWHERE
    CALL xios_orchidee_send_field("CO2_GM",co2_gm)
    CALL xios_orchidee_send_field("CH4_GM",ch4_gm)
    CALL xios_orchidee_send_field("N2O_GM",n2o_gm)
    CALL xios_orchidee_send_field("N2O_PFT_GM",n2o_pft_gm)
!end gmjc


    ! ipcc history
    ! Carbon stock transformed from gC/m2 into kgC/m2
    CALL xios_orchidee_send_field("cVeg",SUM(tot_live_biomass(:,:,icarbon)*veget_cov_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("cLitter",SUM(tot_litter_carb*veget_cov_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("cSoil",SUM(tot_soil_carb*veget_cov_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("cProduct",(prod10_total + prod100_total + prod10_harvest_total + prod100_harvest_total)/1e3)
    CALL xios_orchidee_send_field("cMassVariation",carb_mass_variation/1e3/one_day)

    CALL xios_orchidee_send_field("lai_ipcc",SUM(lai*veget_cov_max,dim=2)) ! m2/m2
    
    ! Carbon fluxes transformed from gC/m2/d into kgC/m2/s
    CALL xios_orchidee_send_field("gpp_ipcc",SUM(gpp_daily*veget_cov_max,dim=2)/1e3/one_day) 
    CALL xios_orchidee_send_field("ra",SUM((resp_maint+resp_growth)*veget_cov_max,dim=2)/1e3/one_day)
    CALL xios_orchidee_send_field("npp_ipcc",SUM(npp_daily*veget_cov_max,dim=2)/1e3/one_day)
    CALL xios_orchidee_send_field("rh",SUM(resp_hetero*veget_cov_max,dim=2)/1e3/one_day)
    CALL xios_orchidee_send_field("fFire",SUM(co2_fire*veget_cov_max,dim=2)/1e3/one_day)
!gmjc
    IF (enable_grazing) THEN
        vartmp(:)=(SUM((harvest_gm+ranimal_gm+ch4_pft_gm) & ! specific
         & *veget_cov_max,dim=2))/1e3/one_day
    ELSE
       vartmp(:)= 0.0
    ENDIF
    CALL xios_orchidee_send_field("fGrazing",vartmp(:))
    IF (enable_grazing) THEN
        vartmp(:)=(SUM((cinput_gm) & ! specific
         & *veget_cov_max,dim=2))/1e3/one_day
    ELSE
       vartmp(:)= 0.0
    ENDIF
    CALL xios_orchidee_send_field("fManureApl",vartmp(:))
!end gmjc 
    CALL xios_orchidee_send_field("fHarvest",harvest_above/1e3/one_day)
    CALL xios_orchidee_send_field("fLuc",cflux_prod_total/1e3/one_day)
    CALL xios_orchidee_send_field("fWoodharvest",cflux_prod_harvest_total/1e3/one_day)
!gmjc
    IF (enable_grazing) THEN
        vartmp(:)=(SUM((gpp_daily+co2_to_bm-(resp_maint+resp_growth+resp_hetero)-co2_fire &
                 -harvest_gm-ranimal_gm-ch4_pft_gm+cinput_gm) & ! specific            
         & *veget_cov_max,dim=2)-cflux_prod_total-cflux_prod_harvest_total-harvest_above)/1e3/one_day
    ELSE
       vartmp(:)=(SUM((gpp_daily+co2_to_bm-(resp_maint+resp_growth+resp_hetero)-co2_fire) &
         & *veget_cov_max,dim=2)-cflux_prod_total-cflux_prod_harvest_total-harvest_above)/1e3/one_day
    ENDIF
!end gmjc
    CALL xios_orchidee_send_field("nbp",vartmp(:))
    CALL xios_orchidee_send_field("fVegLitter",SUM((tot_bm_to_litter(:,:,icarbon) + tot_turnover(:,:,icarbon))*&
         veget_cov_max,dim=2)/1e3/one_day)
    CALL xios_orchidee_send_field("fLitterSoil",SUM(SUM(soilcarbon_input,dim=2)*veget_cov_max,dim=2)/1e3/one_day)

    ! Carbon stock transformed from gC/m2 into kgC/m2
    CALL xios_orchidee_send_field("cLeaf",SUM(biomass(:,:,ileaf,icarbon)*veget_cov_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("cWood",SUM((biomass(:,:,isapabove,icarbon)+biomass(:,:,iheartabove,icarbon))*&
          veget_cov_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("cRoot",SUM(( biomass(:,:,iroot,icarbon) + biomass(:,:,isapbelow,icarbon) + &
         biomass(:,:,iheartbelow,icarbon) )*veget_cov_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("cMisc",SUM(( biomass(:,:,icarbres,icarbon) + biomass(:,:,ifruit,icarbon))*&
         veget_cov_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("cLitterAbove",SUM((litter(:,istructural,:,iabove,icarbon)+&
         litter(:,imetabolic,:,iabove,icarbon))*veget_cov_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("cLitterBelow",SUM((litter(:,istructural,:,ibelow,icarbon)+&
         litter(:,imetabolic,:,ibelow,icarbon))*veget_cov_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("cSoilFast",SUM(carbon(:,iactive,:)*veget_cov_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("cSoilMedium",SUM(carbon(:,islow,:)*veget_cov_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("cSoilSlow",SUM(carbon(:,ipassive,:)*veget_cov_max,dim=2)/1e3)

    ! Vegetation fractions [0,100]
    CALL xios_orchidee_send_field("landCoverFrac",veget_cov_max*100)
    vartmp(:)=zero
    DO j = 2,nvm
       IF (is_deciduous(j)) THEN
          vartmp(:) = vartmp(:) + veget_cov_max(:,j)*100
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("treeFracPrimDec",vartmp)
    vartmp(:)=zero
    DO j = 2,nvm
       IF (is_evergreen(j)) THEN
          vartmp(:) = vartmp(:) + veget_cov_max(:,j)*100
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("treeFracPrimEver",vartmp)
    vartmp(:)=zero
    DO j = 2,nvm
       IF ( .NOT.(is_c4(j)) ) THEN
          vartmp(:) = vartmp(:) + veget_cov_max(:,j)*100
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("c3PftFrac",vartmp)
    vartmp(:)=zero
    DO j = 2,nvm
       IF ( is_c4(j) ) THEN
          vartmp(:) = vartmp(:) + veget_cov_max(:,j)*100
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("c4PftFrac",vartmp)

    ! Carbon fluxes transformed from gC/m2/d into kgC/m2/s
    CALL xios_orchidee_send_field("rGrowth",SUM(resp_growth*veget_cov_max,dim=2)/1e3/one_day)
    CALL xios_orchidee_send_field("rMaint",SUM(resp_maint*veget_cov_max,dim=2)/1e3/one_day)
    CALL xios_orchidee_send_field("nppLeaf",SUM(bm_alloc(:,:,ileaf,icarbon)*veget_cov_max,dim=2)/1e3/one_day)
    CALL xios_orchidee_send_field("nppWood",SUM(bm_alloc(:,:,isapabove,icarbon)*veget_cov_max,dim=2)/1e3/one_day)
    CALL xios_orchidee_send_field("nppRoot",SUM(( bm_alloc(:,:,isapbelow,icarbon) + bm_alloc(:,:,iroot,icarbon) ) * &
         veget_cov_max,dim=2)/1e3/one_day)


    CALL histwrite_p (hist_id_stomate, 'RESOLUTION_X', itime, &
         resolution(:,1), npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'RESOLUTION_Y', itime, &
         resolution(:,2), npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CONTFRAC', itime, &
         contfrac(:), npts, hori_index)

    CALL histwrite_p (hist_id_stomate, 'LITTER_STR_AB', itime, &
         litter(:,istructural,:,iabove,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_MET_AB', itime, &
         litter(:,imetabolic,:,iabove,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_STR_BE', itime, &
         litter(:,istructural,:,ibelow,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_MET_BE', itime, &
         litter(:,imetabolic,:,ibelow,icarbon), npts*nvm, horipft_index)
!spitfiretest
    !CALL histwrite (hist_id_stomate, 'fuel_1hr_met_b', itime, &
    !     fuel_1hr(:,:,imetabolic,icarbon), npts*nvm, horipft_index)
    !CALL histwrite (hist_id_stomate, 'fuel_1hr_str_b', itime, &
    !     fuel_1hr(:,:,istructural,icarbon), npts*nvm, horipft_index)
    !CALL histwrite (hist_id_stomate, 'fuel_10hr_met_b', itime, &
    !     fuel_10hr(:,:,imetabolic,icarbon), npts*nvm, horipft_index)
    !CALL histwrite (hist_id_stomate, 'fuel_10hr_str_b', itime, &
    !     fuel_10hr(:,:,istructural,icarbon), npts*nvm, horipft_index)
    !CALL histwrite (hist_id_stomate, 'fuel_100hr_met_b', itime, &
    !     fuel_100hr(:,:,imetabolic,icarbon), npts*nvm, horipft_index)
    !CALL histwrite (hist_id_stomate, 'fuel_100hr_str_b', itime, &
    !     fuel_100hr(:,:,istructural,icarbon), npts*nvm, horipft_index)
    !CALL histwrite (hist_id_stomate, 'fuel_1000hr_met_b', itime, &
    !     fuel_1000hr(:,:,imetabolic,icarbon), npts*nvm, horipft_index)
    !CALL histwrite (hist_id_stomate, 'fuel_1000hr_str_b', itime, &
    !     fuel_1000hr(:,:,istructural,icarbon), npts*nvm, horipft_index)
!endspittest

    CALL histwrite_p (hist_id_stomate, 'DEADLEAF_COVER', itime, &
         deadleaf_cover, npts, hori_index)

    CALL histwrite_p (hist_id_stomate, 'TOTAL_SOIL_CARB', itime, &
         tot_litter_soil_carb, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'CARBON_ACTIVE', itime, &
         carbon(:,iactive,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'CARBON_SLOW', itime, &
         carbon(:,islow,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'CARBON_PASSIVE', itime, &
         carbon(:,ipassive,:), npts*nvm, horipft_index)

    CALL histwrite_p (hist_id_stomate, 'CARBON_ACTIVE_SURF', itime, &
         carbon_surf(:,iactive,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'CARBON_SLOW_SURF', itime, &
         carbon_surf(:,islow,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'CARBON_PASSIVE_SURF', itime, &
         carbon_surf(:,ipassive,:), npts*nvm, horipft_index)


    CALL histwrite_p (hist_id_stomate, 'TSURF_YEAR', itime, &
                    tsurf_year, npts, hori_index)
!pss:-


    CALL histwrite_p (hist_id_stomate, 'T2M_MONTH', itime, &
         t2m_month, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'T2M_WEEK', itime, &
         t2m_week, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'TSEASON', itime, &
         Tseason, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'TMIN_SPRING_TIME', itime, &
         Tmin_spring_time, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'ONSET_DATE', itime, &
         onset_date(:,:), npts*nvm, horipft_index)

    CALL histwrite_p (hist_id_stomate, 'HET_RESP', itime, &
         resp_hetero(:,:), npts*nvm, horipft_index)
! gmjc
    CALL histwrite_p(hist_id_stomate ,'T2M_14'   ,itime, &
         t2m_14, npts, hori_index)
!    CALL histwrite (hist_id_stomate, 'LITTER_RESP', itime, &
!         resp_hetero_litter_d(:,:), npts*nvm, horipft_index)
!    CALL histwrite (hist_id_stomate, 'ACTIVE_RESP', itime, &
!         resp_hetero_soil_d(:,iactive,:), npts*nvm, horipft_index)
!    CALL histwrite (hist_id_stomate, 'SLOW_RESP', itime, &
!         resp_hetero_soil_d(:,islow,:), npts*nvm, horipft_index)
!    CALL histwrite (hist_id_stomate, 'PASSIVE_RESP', itime, &
!         resp_hetero_soil_d(:,ipassive,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_STR_AVAIL', itime, &
         litter_avail(:,istructural,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_MET_AVAIL', itime, &
         litter_avail(:,imetabolic,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_STR_NAVAIL', itime, &
         litter_not_avail(:,istructural,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_MET_NAVAIL', itime, &
         litter_not_avail(:,imetabolic,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_STR_AVAILF', itime, &
         litter_avail_frac(:,istructural,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_MET_AVAILF', itime, &
         litter_avail_frac(:,imetabolic,:), npts*nvm, horipft_index)

    ! Crop is enabled
    IF (ANY(ok_LAIdev)) THEN
      CALL histwrite_p (hist_id_stomate, 'N_LIMFERT', itime, &
         N_limfert, npts*nvm, horipft_index)
    ENDIF
! end gmjc
    CALL histwrite_p (hist_id_stomate, 'FIREINDEX', itime, &
         fireindex(:,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTERHUM', itime, &
         litterhum_daily, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CO2_FIRE', itime, &
         co2_fire, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'CO2_TAKEN', itime, &
         co2_to_bm, npts*nvm, horipft_index)
    ! land cover change
    CALL histwrite_p (hist_id_stomate, 'CONVFLUX_HAR', itime, &
         convflux(:,iwphar), npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CONVFLUX_LCC', itime, &
         convflux(:,iwplcc), npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CFLUX_PROD10_HAR', itime, &
         cflux_prod10(:,iwphar), npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CFLUX_PROD10_LCC', itime, &
         cflux_prod10(:,iwplcc), npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CFLUX_PROD100_HAR', itime, &
         cflux_prod100(:,iwphar), npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CFLUX_PROD100_LCC', itime, &
         cflux_prod100(:,iwplcc), npts, hori_index)

    CALL histwrite_p (hist_id_stomate, 'HARVEST_ABOVE', itime, &
         harvest_above, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'PROD10_HAR', itime, &
         prod10(:,:,iwphar), npts*11, horip11_index)
    CALL histwrite_p (hist_id_stomate, 'PROD10_LCC', itime, &
         prod10(:,:,iwplcc), npts*11, horip11_index)
    CALL histwrite_p (hist_id_stomate, 'PROD100_HAR', itime, &
         prod100(:,:,iwphar), npts*101, horip101_index)
    CALL histwrite_p (hist_id_stomate, 'PROD100_LCC', itime, &
         prod100(:,:,iwplcc), npts*101, horip101_index)
    CALL histwrite_p (hist_id_stomate, 'FLUX10_HAR', itime, &
         flux10(:,:,iwphar), npts*10, horip10_index)
    CALL histwrite_p (hist_id_stomate, 'FLUX10_LCC', itime, &
         flux10(:,:,iwplcc), npts*10, horip10_index)
    CALL histwrite_p (hist_id_stomate, 'FLUX100_HAR', itime, &
         flux100(:,:,iwphar), npts*100, horip100_index)
    CALL histwrite_p (hist_id_stomate, 'FLUX100_LCC', itime, &
         flux100(:,:,iwplcc), npts*100, horip100_index)
    CALL histwrite_p (hist_id_stomate, 'DefLitSurplus', itime, &
         deflitsup_total, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'DefBioSurplus', itime, &
         defbiosup_total, npts*nvm, horipft_index)

    IF (use_bound_spa) THEN
      CALL histwrite_p (hist_id_stomate, 'bound_spa', itime, &
         bound_spa, npts*nvm, horipft_index)
    ENDIF

    IF (do_now_stomate_lcchange) THEN
        CALL histwrite_p (hist_id_stomate, 'LCC', itime, &
             lcc, npts*nvm, horipft_index)
    ENDIF

    CALL histwrite_p (hist_id_stomate, 'LAI', itime, &
         lai, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'FPC_MAX', itime, &
         fpc_max, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'MAXFPC_LASTYEAR', itime, &
         maxfpc_lastyear, npts*nvm, horipft_index) 
    CALL histwrite_p (hist_id_stomate, 'VEGET_COV_MAX', itime, &
         veget_cov_max, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'NPP', itime, &
         npp_daily, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'GPP', itime, &
         gpp_daily, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'IND', itime, &
         ind, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'CN_IND', itime, &
         cn_ind, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'WOODMASS_IND', itime, &
         woodmass_ind, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'TOTAL_M', itime, &
         tot_live_biomass(:,:,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LEAF_M', itime, &
         biomass(:,:,ileaf,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SAP_M_AB', itime, &
         biomass(:,:,isapabove,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SAP_M_BE', itime, &
         biomass(:,:,isapbelow,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'HEART_M_AB', itime, &
         biomass(:,:,iheartabove,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'HEART_M_BE', itime, &
         biomass(:,:,iheartbelow,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'ROOT_M', itime, &
         biomass(:,:,iroot,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'FRUIT_M', itime, &
         biomass(:,:,ifruit,icarbon), npts*nvm, horipft_index)
!!!!! crop variables
    CALL histwrite_p (hist_id_stomate, 'CROPYIELD', itime, &
         biomass(:,:,ifruit,icarbon), npts*nvm, horipft_index)

    CALL histwrite_p (hist_id_stomate, 'BIOMYIELD', itime, &
         biomass(:,:,ileaf,icarbon)+biomass(:,:,isapabove,icarbon) &
        +biomass(:,:,ifruit,icarbon)+biomass(:,:,icarbres,icarbon), npts*nvm,horipft_index)

    CALL histwrite_p (hist_id_stomate, 'CROP_EXPORT', itime, &
         crop_export, npts*nvm, horipft_index)
!!!!! end crop variables, xuhui
    CALL histwrite_p (hist_id_stomate, 'RESERVE_M', itime, &
         biomass(:,:,icarbres,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'TOTAL_TURN', itime, &
         tot_turnover(:,:,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LEAF_TURN', itime, &
         turnover_daily(:,:,ileaf,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SAP_AB_TURN', itime, &
         turnover_daily(:,:,isapabove,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'ROOT_TURN', itime, &
         turnover_daily(:,:,iroot,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'FRUIT_TURN', itime, &
         turnover_daily(:,:,ifruit,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'TOTAL_BM_LITTER', itime, &
         tot_bm_to_litter(:,:,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LEAF_BM_LITTER', itime, &
         bm_to_litter(:,:,ileaf,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SAP_AB_BM_LITTER', itime, &
         bm_to_litter(:,:,isapabove,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SAP_BE_BM_LITTER', itime, &
         bm_to_litter(:,:,isapbelow,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'HEART_AB_BM_LITTER', itime, &
         bm_to_litter(:,:,iheartabove,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'HEART_BE_BM_LITTER', itime, &
         bm_to_litter(:,:,iheartbelow,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'ROOT_BM_LITTER', itime, &
         bm_to_litter(:,:,iroot,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'FRUIT_BM_LITTER', itime, &
         bm_to_litter(:,:,ifruit,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'RESERVE_BM_LITTER', itime, &
         bm_to_litter(:,:,icarbres,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'MAINT_RESP', itime, &
         resp_maint, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'GROWTH_RESP', itime, &
         resp_growth, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'AGE', itime, &
         age, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'HEIGHT', itime, &
         height, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'MOISTRESS', itime, &
         moiavail_week, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'VCMAX', itime, &
         vcmax, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'TURNOVER_TIME', itime, &
         turnover_time, npts*nvm, horipft_index)
!!DZADD
!    CALL histwrite_p (hist_id_stomate, 'LEAF_FRAC1', itime, leaf_frac(:,:,1), npts*nvm, horipft_index)
!    CALL histwrite_p (hist_id_stomate, 'LEAF_FRAC2', itime, leaf_frac(:,:,2), npts*nvm, horipft_index)
!    CALL histwrite_p (hist_id_stomate, 'LEAF_FRAC3', itime, leaf_frac(:,:,3), npts*nvm, horipft_index)
!    CALL histwrite_p (hist_id_stomate, 'LEAF_FRAC4', itime, leaf_frac(:,:,4), npts*nvm, horipft_index)
!!ENDDZADD

    IF ( hist_id_stomate_IPCC > 0 ) THEN
       vartmp(:)=SUM(tot_live_biomass(:,:,icarbon)*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cVeg", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(tot_litter_carb*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cLitter", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(tot_soil_carb*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cSoil", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=(prod10_total + prod100_total + prod10_harvest_total + prod100_harvest_total)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cProduct", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=carb_mass_variation/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "cMassVariation", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(lai*veget_cov_max,dim=2)
       CALL histwrite_p (hist_id_stomate_IPCC, "lai", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(gpp_daily*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "gpp", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((resp_maint+resp_growth)*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "ra", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(npp_daily*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "npp", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(resp_hetero*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "rh", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(co2_fire*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "fFire", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=harvest_above/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "fHarvest", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=cflux_prod_total/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "fLuc", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=cflux_prod_harvest_total/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "fWoodharvest", itime, &
        vartmp, npts, hori_index)
!gmjc
       IF (enable_grazing) THEN
           vartmp(:)=(SUM((gpp_daily+co2_to_bm-(resp_maint+resp_growth+resp_hetero)-co2_fire &
                    -harvest_gm-ranimal_gm-ch4_pft_gm+cinput_gm) & ! specific
            & *veget_cov_max,dim=2)-cflux_prod_total-cflux_prod_harvest_total-harvest_above)/1e3/one_day
       ELSE
          vartmp(:)=(SUM((gpp_daily+co2_to_bm-(resp_maint+resp_growth+resp_hetero)-co2_fire) &
            & *veget_cov_max,dim=2)-cflux_prod_total-cflux_prod_harvest_total-harvest_above)/1e3/one_day
       ENDIF
!end gmjc
!       vartmp(:)=(SUM((gpp_daily-(resp_maint+resp_growth+resp_hetero)-co2_fire) &
!            &        *veget_cov_max,dim=2)-cflux_prod_total-harvest_above)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "nbp", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((tot_bm_to_litter(:,:,icarbon) + tot_turnover(:,:,icarbon))*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "fVegLitter", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(SUM(soilcarbon_input,dim=2)*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "fLitterSoil", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(biomass(:,:,ileaf,icarbon)*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cLeaf", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((biomass(:,:,isapabove,icarbon)+biomass(:,:,iheartabove,icarbon))*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cWood", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(( biomass(:,:,iroot,icarbon) + biomass(:,:,isapbelow,icarbon) + biomass(:,:,iheartbelow,icarbon) ) &
            &        *veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cRoot", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(( biomass(:,:,icarbres,icarbon) + biomass(:,:,ifruit,icarbon))*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cMisc", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((litter(:,istructural,:,iabove,icarbon)+litter(:,imetabolic,:,iabove,icarbon))*&
            veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cLitterAbove", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((litter(:,istructural,:,ibelow,icarbon)+litter(:,imetabolic,:,ibelow,icarbon))*&
            veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cLitterBelow", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(carbon(:,iactive,:)*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cSoilFast", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(carbon(:,islow,:)*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cSoilMedium", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(carbon(:,ipassive,:)*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cSoilSlow", itime, &
            vartmp, npts, hori_index)
       DO j=1,nvm
          histvar(:,j)=veget_cov_max(:,j)*100
       ENDDO
       CALL histwrite_p (hist_id_stomate_IPCC, "landCoverFrac", itime, &
            histvar, npts*nvm, horipft_index)
       !-
       vartmp(:)=zero
       DO j = 2,nvm
          IF (is_deciduous(j)) THEN
             vartmp(:) = vartmp(:) + veget_cov_max(:,j)*100
          ENDIF
       ENDDO
       CALL histwrite_p (hist_id_stomate_IPCC, "treeFracPrimDec", itime, &
            vartmp, npts, hori_index)
       !-
       vartmp(:)=zero
       DO j = 2,nvm
          IF (is_evergreen(j)) THEN
             vartmp(:) = vartmp(:) + veget_cov_max(:,j)*100
          ENDIF
       ENDDO
       CALL histwrite_p (hist_id_stomate_IPCC, "treeFracPrimEver", itime, &
            vartmp, npts, hori_index)
       !-
       vartmp(:)=zero
       DO j = 2,nvm
          IF ( .NOT.(is_c4(j)) ) THEN
             vartmp(:) = vartmp(:) + veget_cov_max(:,j)*100
          ENDIF
       ENDDO
       CALL histwrite_p (hist_id_stomate_IPCC, "c3PftFrac", itime, &
            vartmp, npts, hori_index)
       !-
       vartmp(:)=zero
       DO j = 2,nvm
          IF ( is_c4(j) ) THEN
             vartmp(:) = vartmp(:) + veget_cov_max(:,j)*100
          ENDIF
       ENDDO
       CALL histwrite_p (hist_id_stomate_IPCC, "c4PftFrac", itime, &
            vartmp, npts, hori_index)
       !-
       vartmp(:)=SUM(resp_growth*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "rGrowth", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(resp_maint*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "rMaint", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(bm_alloc(:,:,ileaf,icarbon)*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "nppLeaf", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(bm_alloc(:,:,isapabove,icarbon)*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "nppWood", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(( bm_alloc(:,:,isapbelow,icarbon) + bm_alloc(:,:,iroot,icarbon) )*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "nppRoot", itime, &
            vartmp, npts, hori_index)

       CALL histwrite_p (hist_id_stomate_IPCC, 'RESOLUTION_X', itime, &
            resolution(:,1), npts, hori_index)
       CALL histwrite_p (hist_id_stomate_IPCC, 'RESOLUTION_Y', itime, &
            resolution(:,2), npts, hori_index)
       CALL histwrite_p (hist_id_stomate_IPCC, 'CONTFRAC', itime, &
            contfrac(:), npts, hori_index)

    ENDIF

    IF (printlev>=4) WRITE(numout,*) 'Leaving stomate_lpj'

  END SUBROUTINE StomateLpj


!! ================================================================================================================================
!! SUBROUTINE   : harvest
!!
!>\BRIEF        Harvest of croplands
!!
!! DESCRIPTION  : To take into account biomass harvest from crop (mainly to take 
!! into account for the reduced litter input and then decreased soil carbon. it is a 
!! constant (40\%) fraction of above ground biomass.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::harvest_above the harvested biomass
!!
!! REFERENCE(S) :
!! - Piao, S., P. Ciais, P. Friedlingstein, N. de Noblet-Ducoudre, P. Cadule, N. Viovy, and T. Wang. 2009. 
!!   Spatiotemporal patterns of terrestrial carbon cycle during the 20th century. Global Biogeochemical 
!!   Cycles 23:doi:10.1029/2008GB003339.
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE harvest(npts, dt_days, veget_cov_max, &
       bm_to_litter, turnover_daily, &
       harvest_above,harvest_bio)

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER, INTENT(in)                                    :: npts            !! Domain size (unitless) 
    REAL(r_std), INTENT(in)                                :: dt_days         !! Time step (days)                               
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)           :: veget_cov_max       !! new "maximal" coverage fraction of a PFT (LAI -> 
                                                                              !! infinity) on ground @tex $(m^2 m^{-2})$ @endtex 
    
   !! 0.2 Output variables
   
   !! 0.3 Modified variables

    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: bm_to_litter !! [DISPENSABLE] conversion of biomass to litter 
                                                                                     !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: turnover_daily   !! Turnover rates 
                                                                               !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts), INTENT(inout)            :: harvest_above    !! harvest above ground biomass for agriculture 
                                                                               !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
!!!qcj++ peatland
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)            :: harvest_bio
    !! 0.4 Local variables

    INTEGER(i_std)                                         :: i, j, k, l, m    !! indices                       
    REAL(r_std)                                            :: above_old        !! biomass of previous time step 
                                                                               !! @tex $(gC m^{-2})$ @endtex 
!_ ================================================================================================================================

  !! 1. Yearly initialisation

    above_old             = zero
    harvest_above         = zero
    harvest_bio           = zero 

    DO i = 1, npts
       DO j = 1,nvm
          IF ((.NOT. natural(j)) .AND. (.NOT. is_peat(j))) THEN
             above_old = turnover_daily(i,j,ileaf,icarbon) + turnover_daily(i,j,isapabove,icarbon) + &
                  &       turnover_daily(i,j,iheartabove,icarbon) + turnover_daily(i,j,ifruit,icarbon) + &
                  &       turnover_daily(i,j,icarbres,icarbon) + turnover_daily(i,j,isapbelow,icarbon) + &
                  &       turnover_daily(i,j,iheartbelow,icarbon) + turnover_daily(i,j,iroot,icarbon)

             turnover_daily(i,j,ileaf,icarbon) = turnover_daily(i,j,ileaf,icarbon)*frac_turnover_daily
             turnover_daily(i,j,isapabove,icarbon) = turnover_daily(i,j,isapabove,icarbon)*frac_turnover_daily
             turnover_daily(i,j,isapbelow,icarbon) = turnover_daily(i,j,isapbelow,icarbon)*frac_turnover_daily
             turnover_daily(i,j,iheartabove,icarbon) = turnover_daily(i,j,iheartabove,icarbon)*frac_turnover_daily
             turnover_daily(i,j,iheartbelow,icarbon) = turnover_daily(i,j,iheartbelow,icarbon)*frac_turnover_daily
             turnover_daily(i,j,iroot,icarbon) = turnover_daily(i,j,iroot,icarbon)*frac_turnover_daily
             turnover_daily(i,j,ifruit,icarbon) = turnover_daily(i,j,ifruit,icarbon)*frac_turnover_daily
             turnover_daily(i,j,icarbres,icarbon) = turnover_daily(i,j,icarbres,icarbon)*frac_turnover_daily
             harvest_above(i)  = harvest_above(i) + veget_cov_max(i,j) * above_old *(un - frac_turnover_daily)
             harvest_bio(i,j) = harvest_bio(i,j)+veget_cov_max(i,j) * above_old *(un- frac_turnover_daily) 
          ENDIF
       ENDDO
    ENDDO

!!$    harvest_above = harvest_above
  END SUBROUTINE harvest

!! ================================================================================================================================
!! SUBROUTINE   : lpj_cover_peat
!!
!>\BRIEF        Peatland area fraction is calculated by TOPMODEL
!               1. IF ok_dgvm, then non-peat vegetations' area is calculated by DGVM (bioclimatic conditions)
!               2. IF .NOT. ok_dgvm, then non-peat vegetations' area is read from vegetation maps 
!               We adjust non-peatland vegetations coverage according to peatland cover.
!                     a. Peatland initiates/expands ---> all non-peatland vegetations shrink, proportionally
!                     b. Peatland shrink/disappear ----> all non-peatland vegetations expand, proportionally
!               When peatland shrinks, change of area and C density will be saved (*_save). Then when pealtand expands, peatland will take
!               carbon from the *_save first, if fsave is not enough, extra carbon will be taken from natural PFTs
!!
!! \n
!_ ================================================================================================================================

  SUBROUTINE lpj_cover_peat(npts, lalo,cn_ind, ind, biomass, &
       veget_cov_max_new, veget_cov_max, veget_cov_max_old, &
       litter, litter_avail, litter_not_avail, carbon, &
       fuel_1hr, fuel_10hr, fuel_100hr, fuel_1000hr, &
       turnover_daily, bm_to_litter, &
       co2_to_bm, co2_fire, resp_hetero, resp_maint, resp_growth, gpp_daily, &
       deepC_a, deepC_s, deepC_p, &
       dt_days,age, PFTpresent, senescence, when_growthinit,&
       everywhere, leaf_frac, lm_lastyearmax, npp_longterm,&
       carbon_save,deepC_a_save,deepC_s_save,deepC_p_save,delta_fsave,liqwt_max_lastyear)

!! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                                  :: npts             !! Domain size (unitless)  
    REAL(r_std),DIMENSION(npts,2),INTENT(in)                   :: lalo
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)                :: cn_ind           !! Crown area 
                                                                                    !! @tex $(m^2)$ @endtex per individual
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)                :: ind              !! Number of individuals 
                                                                                    !! @tex $(m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)                :: veget_cov_max_new    !! saved fpeat 
                                                                                  !! @tex ($gC individual^{-1}$) @endtex
    REAL(r_std), INTENT(in)                                     :: dt_days !! Time step of vegetation dynamics for stomate
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)                :: veget_cov_max_old
    REAL(r_std), DIMENSION(npts), INTENT(in)                    :: liqwt_max_lastyear
    !! 0.2 Output variables

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(npts,nlitt,nvm,nlevs,nelements), INTENT(inout) :: litter    !! Metabolic and structural litter, above and 
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements), INTENT(inout)                 :: fuel_1hr
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements), INTENT(inout)                 :: fuel_10hr
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements), INTENT(inout)                 :: fuel_100hr
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements), INTENT(inout)                 :: fuel_1000hr
                                                                                       !! below ground @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nlitt,nvm), INTENT(inout):: litter_avail
    REAL(r_std), DIMENSION(npts,nlitt,nvm) , INTENT(inout):: litter_not_avail
    REAL(r_std), DIMENSION(npts,ncarb,nvm), INTENT(inout)             :: carbon        !! Carbon pool: active, slow, or passive 
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: biomass        !! Biomass @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)                  :: veget_cov_max      !! "Maximal" coverage fraction of a PFT (LAI->
                                                                                       !! infinity) on ground (unitless)
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: turnover_daily !! Turnover rates (gC m^{-2} day^{-1})
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: bm_to_litter   !! Conversion of biomass to litter 
                                                                                       !! @tex $(gC m^{-2} day^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: co2_to_bm             !! biomass up take for establishment           
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: co2_fire
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: resp_hetero
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: resp_maint
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: resp_growth
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: gpp_daily

    REAL(r_std), DIMENSION(npts,ndeep,nvm), INTENT(inout)         :: deepC_a           !! Permafrost soil carbon (g/m**3) active
    REAL(r_std), DIMENSION(npts,ndeep,nvm), INTENT(inout)         :: deepC_s           !! Permafrost soil carbon (g/m**3) slow
    REAL(r_std), DIMENSION(npts,ndeep,nvm), INTENT(inout)         :: deepC_p           !! Permafrost soil carbon (g/m**3) passive

    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: age !! mean age (years)
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)               :: senescence       !! plant senescent (only for deciduous trees) Set
                                                                                  !! to .FALSE. if PFT is introduced or killed
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)               :: PFTpresent       !! Is pft there (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: everywhere       !! is the PFT everywhere in the grid box or very 
                                                                                  !! localized (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: when_growthinit  !! how many days ago was the beginning of the 
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout) :: leaf_frac        !! fraction of leaves in leaf age class 
                                                                                  !! (unitless)
                                                                                  !! growing season (days)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: lm_lastyearmax   !! last year's maximum leaf mass for each PFT 
                                                                                  !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: npp_longterm     !! "long term" net primary productivity 
                                                                                  !! @tex ($gC m^{-2} year^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts,ncarb,nvm), INTENT(inout)     :: carbon_save
    REAL(r_std), DIMENSION(npts,ndeep), INTENT(inout)         :: deepC_a_save
    REAL(r_std), DIMENSION(npts,ndeep), INTENT(inout)         :: deepC_s_save
    REAL(r_std), DIMENSION(npts,ndeep), INTENT(inout)         :: deepC_p_save
    REAL(r_std), DIMENSION(npts), INTENT(inout)               :: delta_fsave

    !! 0.4 Local variables
     REAL(r_std), DIMENSION(npts,ncarb,nvm)                     :: carbon_save_tmp
    INTEGER(i_std)                                              :: i,j,k,m,l               !! Index (unitless)
    REAL(r_std), DIMENSION(npts,nlitt,nlevs,nelements)          :: dilu_lit              !! Litter dilution @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nparts,nelements)               :: dilu_bio              !! Biomass dilution 
    REAL(r_std), DIMENSION(npts,nparts,nelements)               :: dilu_turnover_daily
    REAL(r_std), DIMENSION(npts,nparts,nelements)               :: dilu_bm_to_litter
    REAL(r_std), DIMENSION(npts)                                :: dilu_co2flux_new
    REAL(r_std), DIMENSION(npts)                                :: dilu_gpp_daily
    REAL(r_std), DIMENSION(npts)                                :: dilu_resp_growth
    REAL(r_std), DIMENSION(npts)                                :: dilu_resp_maint
    REAL(r_std), DIMENSION(npts)                                :: dilu_resp_hetero
    REAL(r_std), DIMENSION(npts)                                :: dilu_co2_to_bm
    REAL(r_std), DIMENSION(npts)                                :: dilu_co2_fire
    REAL(r_std), DIMENSION(npts,nvm)                            :: co2flux_new
    REAL(r_std), DIMENSION(npts,nvm)                            :: co2flux_old
    REAL(r_std), DIMENSION(npts,ncarb,nvm)                       :: carbon_old
    REAL(r_std), DIMENSION(npts,ncarb,nvm)                       :: carbon_tmp

    REAL(r_std), DIMENSION(nvm)                                 :: delta_veg        !! Conversion factors (unitless)
    REAL(r_std)                                                 :: delta_veg_sum    !! Conversion factors (unitless)
    REAL(r_std)                                                 :: exp_nat_sum
    REAL(r_std)                                                 :: nat_sum
    REAL(r_std), DIMENSION(npts,nlitt,nelements)                :: dilu_f1hr        !! Litter dilution @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nlitt,nelements)                :: dilu_f10hr       !! Litter dilution @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nlitt,nelements)                :: dilu_f100hr      !! Litter dilution @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nlitt,nelements)                :: dilu_f1000hr     !! Litter dilution @tex $(gC m^{-2})$ @endtex
    REAL(r_std),DIMENSION(npts,nvm)                             :: delta_ind          !! change in number of individuals
    REAL(r_std), DIMENSION(npts)                                :: frac_nat
    REAL(r_std), DIMENSION(npts)                                :: sum_veget_natveg !! Conversion factors (unitless)

    REAL(r_std)                                                  :: sum_veg
    REAL(r_std)                                                  :: exp_nat_sum2
    REAL(r_std)                                                  :: sumvpeat_old, sumvpeat ! last an new sum of peatland vegetation
    REAL(r_std)                                                  :: rapport  ! (S-B) / (S-A)
    REAL(r_std)                                                  :: SUMveg
    REAL(r_std), DIMENSION(npts,nvm)                             :: veget_tmp
    REAL(r_std), DIMENSION(npts)                                 :: delta_fpeat
    REAL(r_std), DIMENSION(npts)                                 :: diff_fpeat
    REAL(r_std), DIMENSION(npts)                                 :: extra_fra
    REAL(r_std), DIMENSION(npts)                                 :: sum_nat
    REAL(r_std), DIMENSION(npts,ncarb,nvm)                       :: dilu_carbon
    REAL(r_std), DIMENSION(npts,ndeep,nvm)                       :: dilu_a
    REAL(r_std), DIMENSION(npts,ndeep,nvm)                       :: dilu_s
    REAL(r_std), DIMENSION(npts,ndeep,nvm)                       :: dilu_p
    REAL(r_std),DIMENSION(npts,nparts,nelements)                 :: bm_new
    REAL(r_std),DIMENSION(npts,nparts,nelements)                 :: bm_gain
    REAL(r_std), DIMENSION(npts,nvm)                             :: co2_take
    REAL(r_std), DIMENSION(npts,ncarb)                           :: carbon_obtain
    REAL(r_std), DIMENSION(npts)                                 :: biomass_before
    REAL(r_std), DIMENSION(npts)                                 :: biomass_after
    REAL(r_std), DIMENSION(npts)                                 :: litter_before
    REAL(r_std), DIMENSION(npts)                                 :: litter_after
    REAL(r_std), DIMENSION(npts)                                 :: soilc_before
    REAL(r_std), DIMENSION(npts)                                 :: soilc_after
    REAL(r_std), DIMENSION(npts)                                 :: delta_Cpool
    REAL(r_std), DIMENSION(npts,ncarb)                           :: excess_C
    REAL(r_std)                                                  :: excess_tmp
    REAL(r_std), DIMENSION(npts)                                 :: flux_before
    REAL(r_std), DIMENSION(npts)                                 :: flux_after
    REAL(r_std)                                                  :: delta_nonpeat_tmp
    ! peatland PFT number in nvm
    INTEGER(i_std)                                            :: c3crop, c4crop, peat_grass, peat_c3crop, peat_c4crop, peat_forest
!_  
!================================================================================================================================
   
   IF (printlev>=3) WRITE(numout,*) 'Entering lpj_cover_peat'

      ! WRITE(numout,*) 'chunjing Entering lpj_cover_peat'

       biomass_before(:)=zero
       biomass_after(:)=zero
       litter_before(:)=zero
       litter_after(:)=zero
       soilc_before(:)=zero
       soilc_after(:)=zero
       carbon_save_tmp(:,:,:)=zero  
       excess_C(:,:)=zero
 !    flux_before(:)= zero
   !    flux_after(:)=zero

       DO i=1, npts
          DO l=1,ncarb
             DO j=1,nvm 
                soilc_before(i)=soilc_before(i)+carbon(i,l,j)*veget_cov_max_old(i,j)
             ENDDO
          ENDDO
       ENDDO

       !! 1.1  Calculate initial values of vegetation cover
       frac_nat(:) = un
       sum_veget_natveg(:) = zero
       veget_cov_max(:,ibare_sechiba) = un
       co2flux_new = undef
       co2flux_old = undef

       carbon_old(:,:,:)=carbon(:,:,:)

       bm_new(:,:,:)=zero
       co2_take(:,:)=zero
       bm_gain(:,:,:)=zero
       carbon_obtain(:,:)=zero

       ! search for non-peatland crop
    DO j=1,nvm
       IF (.NOT. natural(j) .AND. is_croppeat(j) .AND. .NOT. is_c4(j)) THEN
         peat_c3crop = j  ! C3 peatland crop
       ENDIF
       IF (.NOT. natural(j) .AND. is_croppeat(j) .AND. is_c4(j)) THEN
         peat_c4crop = j  ! C4 peatland crop
       ENDIF
       IF (.NOT. natural(j) .AND. .NOT. is_croppeat(j) .AND. .NOT. is_c4(j)) THEN
         c3crop = j  ! C3 non-peatland crop
       ENDIF
       IF (.NOT. natural(j) .AND. .NOT. is_croppeat(j) .AND. is_c4(j)) THEN
         c4crop = j  ! C4 non-peatland crop
       ENDIF
       IF (is_peat(j) .AND. .NOT. is_croppeat(j) .AND. .NOT. is_tree(j)) THEN
         peat_grass = j  ! C3 peatland grass
       ENDIF
       IF (natural(j) .AND. is_peat(j) .AND. .NOT. is_croppeat(j) .AND. is_tree(j)) THEN
         peat_forest = j ! C3 peatland forest
       ENDIF
    ENDDO

    delta_fpeat(:) = zero

    IF (ok_dgvm) THEN
!!! 1. Fraction of non-peat vegetations are results of bioclimatic conditions
       DO j = 2,nvm ! loop over PFTs
          IF ( natural(j) .AND. .NOT. pasture(j) .AND. .NOT. is_peat(j)) THEN
        
             ! Summation of individual tree crown area to get total foliar projected coverage
             veget_cov_max(:,j) = ind(:,j) * cn_ind(:,j)
             sum_veget_natveg(:) = sum_veget_natveg(:) + veget_cov_max(:,j)
          ELSE

             !fraction occupied by agriculture needs to be substracted for the DGVM
             !this is used below to constrain veget for natural vegetation, see below
             frac_nat(:) = frac_nat(:) - veget_cov_max(:,j)

          ENDIF

       ENDDO ! loop over PFTs

       DO i = 1, npts ! loop over grid points

          ! Recalculation of vegetation projected coverage when ::frac_nat was below ::sum_veget_natveg
          ! It means that non-natural vegetation will recover ::veget_cov_max as natural vegetation
          IF (sum_veget_natveg(i) .GT. frac_nat(i) .AND. frac_nat(i) .GT. min_stomate) THEN
             DO j = 2,nvm ! loop over PFTs
                IF( natural(j) .AND. .NOT. pasture(j) .AND. .NOT. is_peat(j) ) THEN
                   veget_cov_max(i,j) =  veget_cov_max(i,j) * frac_nat(i) / sum_veget_natveg(i)
                ENDIF
             ENDDO ! loop over PFTs
          ENDIF
       ENDDO ! loop over grid points

       ! Renew veget_cov_max of bare soil as 0 to difference of veget_cov_max (ibare_sechiba) 
       ! to current veget_cov_max
       DO j = 2,nvm ! loop over PFTs
          veget_cov_max(:,ibare_sechiba) = veget_cov_max(:,ibare_sechiba) - veget_cov_max(:,j)
       ENDDO ! loop over PFTs
       veget_cov_max(:,ibare_sechiba) = MAX( veget_cov_max(:,ibare_sechiba), zero )

       veget_tmp(:,:)=veget_cov_max(:,:)

!!! 2. Fraction of peatland vegetation are calculated by TOPMODEL
!!! First establishment of pealtand vegetation, companied by increase in biomass
       DO i = 1, npts ! Loop over # pixels - domain size
         delta_veg(:) = veget_cov_max_new(i,:)-veget_tmp(i,:)
         DO j=1, nvm ! Loop over # PFTs
            IF (is_peat(j)) THEN
!!!change veget_cov_max to veget_cov_max_new
              IF (delta_veg(j) .LT. -min_stomate) THEN
  !!!peatland contracts, no limitation
                 IF (veget_cov_max_new(i,j) > min_stomate) THEN
                    veget_cov_max(i,j)=veget_cov_max_new(i,j)
                 ELSE
                    veget_cov_max(i,j)=zero
                 ENDIF
              ELSEIF(delta_veg(j) .GT. min_stomate) THEN
  !!!peatland expands
                 IF (veget_tmp(i,j) .LE. zero) THEN   !!!veget_tmp(i,j) .LE. zero
     !!!first initiation of peatland              
                    veget_cov_max(i,j)=veget_cov_max_new(i,j)
                 ELSE
     !!peatland expands from non-zero, the soil need to be wet enough to support growing of peat vegetations
                    IF (liqwt_max_lastyear(i) > 0.6) THEN 
                       veget_cov_max(i,j)=veget_cov_max_new(i,j)   
                    ENDIF
                 ENDIF
              ENDIF
              IF ( veget_cov_max(i,j)-veget_tmp(i,j) > min_stomate) THEN
            !! Initial setting of peatland vegetation new establishment
                IF (veget_tmp(i,j) .LE. zero)  THEN
                   IF (is_tree(j)) THEN
                       cn_ind(i,j) = cn_sapl(j)
                   ELSE
                       cn_ind(i,j)=un
                   ENDIF
                   ind(i,j)= veget_cov_max(i,j) / cn_ind(i,j)
                   PFTpresent(i,j) = .TRUE.
                   everywhere(i,j) = un
                   senescence(i,j) = .FALSE.
                   age(i,j) = zero
                   when_growthinit(i,j) = large_value
                   leaf_frac(i,j,1) = 1.0
                   npp_longterm(i,j) = npp_longterm_init
                   lm_lastyearmax(i,j) = bm_sapl(j,ileaf,icarbon) * ind(i,j)
                  ! bm_new(i,:,icarbon) = bm_sapl(j,:,icarbon) * ind(i,j) /veget_cov_max(i,j)
                  ! co2_take(i,j) = SUM(bm_new(i,:,icarbon))/ dt_days
                ENDIF
              ENDIF
            ENDIF
         ENDDO
       ENDDO

!!! 3. Adjust fractions of Non-peatland vegetations accordingly

       DO i = 1, npts ! Loop over # pixels - domain size
         sum_veg=SUM(veget_tmp(i,:))
         sumvpeat=zero
         sumvpeat_old=zero
         DO j = 1,nvm
           IF (is_peat(j)) THEN
             sumvpeat=sumvpeat+veget_cov_max(i,j)
             sumvpeat_old=sumvpeat_old+veget_tmp(i,j)
           ENDIF
         ENDDO

!!!!3.1 Peat vegetation area increase:
         IF (sumvpeat .GT. sumvpeat_old) THEN
! All non-pealtand PFTs (natural vegetation +baresoil) decreases, keep the the proportion of natural PFTs 
           rapport = ( sum_veg - sumvpeat ) / ( sum_veg - sumvpeat_old)
           DO j = 1, nvm
             IF ( natural(j) .AND. .NOT. is_peat(j)) THEN
               veget_cov_max(i,j) = veget_tmp(i,j) * rapport
             ENDIF
           ENDDO
         ELSE
! Peat vegetation decrease: natural vegetation fractions will not change, the decrease of peat is replaced by bare soil. 
! The DGVM will re-introduce natural PFT's.
           DO j = 1, nvm
             IF (j==1) THEN
                veget_cov_max(i,j)=veget_tmp(i,j)+sumvpeat_old-sumvpeat
             ENDIF
           ENDDO
         ENDIF
       ENDDO
    ELSE !(ok_dgvm=False)
!!!non-peat vegetations and crops' area are read from PFT maps, peatland fraction is calculated by TOPMODEL
!!!substract peatland fraction from natural vegetations (area of crops will not change and no crops growing in peatland)
!!! 1. Fraction of peatland vegetation are calculated by TOPMODEL
!!! First establishment of pealtand vegetation, companied by increase in biomass
       DO i = 1, npts ! Loop over # pixels - domain size
         delta_veg(:) = veget_cov_max_new(i,:)-veget_cov_max_old(i,:)
         DO j=1, nvm ! Loop over # PFTs
            IF (is_peat(j)) THEN
!!!change veget_cov_max to veget_cov_max_new
              IF (delta_veg(j) .LT. -min_stomate) THEN
  !!!peatland contracts, no limitation
                 IF (veget_cov_max_new(i,j) > min_stomate) THEN
                    veget_cov_max(i,j)=veget_cov_max_new(i,j)
                 ELSE
                    veget_cov_max(i,j)=zero
                 ENDIF
              ELSEIF(delta_veg(j) .GT. min_stomate) THEN
  !!!peatland expands
                 IF (veget_cov_max_old(i,j) .LE. zero) THEN   
     !!!first initiation of peatland              
                    veget_cov_max(i,j)=veget_cov_max_new(i,j)
                 ELSE
     !!peatland expands from non-zero, the soil need to be wet enough to support growing of peat vegetations
                    IF (liqwt_max_lastyear(i) > 0.6) THEN
                       veget_cov_max(i,j)=veget_cov_max_new(i,j)
                    ENDIF
                 ENDIF
              ENDIF
            ENDIF
         ENDDO
       ENDDO

!!! 2. adjust fractions of Non-peatland natural vegetations accordingly
!!!!!!! no crops growing in peatland
       DO j = 1,nvm ! loop over PFTs
          IF ( natural(j) .AND. .NOT. pasture(j) .AND. .NOT. is_peat(j)) THEN
             sum_veget_natveg(:) = sum_veget_natveg(:) + veget_cov_max_new(:,j)
          ENDIF
       ENDDO
       DO i=1, npts
          !!!area of crops (PFT12, PFT13) are from PFT maps and will not occupied by peatland
          veget_cov_max(i,c3crop)= veget_cov_max_new(i,c3crop)
          veget_cov_max(i,c4crop)= veget_cov_max_new(i,c4crop)
          IF (veget_cov_max(i,peat_grass) .LE. sum_veget_natveg(i)) THEN
             IF (sum_veget_natveg(i) .GT. min_stomate) THEN
               DO j=1,nvm
                  IF ( natural(j) .AND. .NOT. pasture(j) .AND. .NOT. is_peat(j)) THEN
                    veget_cov_max(i,j)=veget_cov_max_new(i,j)-veget_cov_max_new(i,j)/sum_veget_natveg(i)*veget_cov_max(i,peat_grass)
                  ENDIF
               ENDDO
             ELSE
               veget_cov_max(i,peat_grass)=zero
               DO j=1,nvm
                  IF (.NOT. is_peat(j)) THEN
                    veget_cov_max(i,j)=veget_cov_max_new(i,j)
                  ENDIF
               ENDDO
             ENDIF
          ELSE
             !! JC comment do we need to consider peat_forest? 
             veget_cov_max(i,peat_grass)=sum_veget_natveg(i)
             DO j=1,nvm
                IF ( natural(j) .AND. .NOT. pasture(j) .AND. .NOT. is_peat(j)) THEN
                   veget_cov_max(i,j)=zero
                ENDIF
             ENDDO 
             veget_cov_max(i,peat_c3crop)=veget_cov_max_new(i,peat_c3crop) !!!=0
             veget_cov_max(i,peat_c4crop)=veget_cov_max_new(i,peat_c4crop) !!!=0
          ENDIF
       ENDDO

       DO i=1, npts
          DO j = 1,nvm
             IF (is_peat(j)) THEN
                IF ( veget_cov_max(i,j)-veget_cov_max_old(i,j) > min_stomate) THEN
                   IF (veget_cov_max_old(i,j) .LE. zero)  THEN
                      IF (is_tree(j)) THEN
                        cn_ind(i,j) = cn_sapl(j)
                      ELSE
                        cn_ind(i,j)=un
                      ENDIF
                      ind(i,j)= veget_cov_max(i,j) / cn_ind(i,j)
                      PFTpresent(i,j) = .TRUE.
                      everywhere(i,j) = un
                      senescence(i,j) = .FALSE.
                      age(i,j) = zero
                      when_growthinit(i,j) = large_value
                      leaf_frac(i,j,1) = 1.0
                      npp_longterm(i,j) = npp_longterm_init
                      lm_lastyearmax(i,j) = bm_sapl(j,ileaf,icarbon) * ind(i,j)
                     ! bm_new(i,:,icarbon) = bm_sapl(j,:,icarbon) * ind(i,j)/veget_cov_max(i,j)
                     ! co2_take(i,j) = SUM(bm_new(i,:,icarbon))/ dt_days
                   ENDIF
                ENDIF    
             ENDIF
          ENDDO
       ENDDO
    ENDIF

!!! Correct vegetation fraction, normalize fractions of  veget_cov_max smaller than min_vegfrac 
!!! to avoid numerical error
    DO i = 1, npts ! loop over grid points
      DO j=1,nvm
        IF ( veget_cov_max(i,j) .LT. min_vegfrac ) THEN
           veget_cov_max(i,j) = zero
        ENDIF
      ENDDO 
      SUMveg =SUM(veget_cov_max(i,:))
      veget_cov_max(i,:) = veget_cov_max(i,:)/SUMveg
    ENDDO

!!!Adjust C fluxes and C pools
 
    DO i = 1, npts ! loop over grid points
      IF ( ABS( SUM(veget_cov_max(i,:)) - SUM(veget_cov_max_new(i,1:13))) > 1000.*min_stomate ) THEN
         WRITE(numout,*) 'qcj check lpj_cover_peat,veget',lalo(i,:)
         WRITE(numout,*) 'qcj check veget,sum_new',SUM(veget_cov_max(i,:))
      ENDIF

      ! Calculate the change in veget_cov_max between previous time step and current time step
       delta_veg(:) = veget_cov_max(i,:)-veget_cov_max_old(i,:)
       delta_veg_sum = SUM(delta_veg,MASK=delta_veg.LT.zero)
       DO j=1,nvm
         IF (is_peat(j)) THEN
           delta_fpeat(i)=delta_fpeat(i) + delta_veg(j)
         ENDIF
       ENDDO
       exp_nat_sum=zero
       DO j=1,nvm
          IF ((.NOT. is_peat(j)) .AND. (delta_veg(j) .GT. min_stomate)) THEN
             exp_nat_sum=exp_nat_sum+delta_veg(j)
          ENDIF
       ENDDO
 
       dilu_lit(i,:,:,:) = zero
       dilu_f1hr(i,:,:) = zero
       dilu_f10hr(i,:,:) = zero
       dilu_f100hr(i,:,:) = zero
       dilu_f1000hr(i,:,:) = zero

       dilu_turnover_daily(i,:,:)=zero
       dilu_bm_to_litter(i,:,:)=zero
       dilu_co2flux_new(i)=zero
       dilu_gpp_daily(i)=zero
       dilu_resp_growth(i)=zero
       dilu_resp_maint(i)=zero
       dilu_resp_hetero(i)=zero
       dilu_co2_to_bm(i)=zero
       dilu_co2_fire(i)=zero

       dilu_bio(i,:,:) = zero

       dilu_carbon(i,:,:)=zero
       dilu_a(i,:,:)=zero
       dilu_s(i,:,:)=zero
       dilu_p(i,:,:)=zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!Deal with Biomass C
!      IF (delta_fpeat(i) .GE. min_stomate) THEN
!         DO j=1, nvm
!            IF (is_peat(j)) THEN
!              bm_new(i,:,:)=delta_fpeat(i)* bm_sapl(j,:,:)
!              IF (veget_cov_max(i,j).GT. min_stomate) THEN
!                 co2_to_bm(i,j) = co2_to_bm(i,j) + (SUM(bm_new(i,:,icarbon))/(dt_days*veget_cov_max(i,j)))
!                 biomass(i,j,:,:)= (biomass(i,j,:,:)* veget_cov_max_old(i,j)+bm_new(i,:,:))/ veget_cov_max(i,j)
!              ENDIF
!            ELSE
!              IF(veget_cov_max(i,j).GT.min_stomate) THEN
!                biomass(i,j,:,:) = biomass(i,j,:,:) * veget_cov_max_old(i,j) /veget_cov_max(i,j)
!              ENDIF
!            ENDIF
!         ENDDO
!      ELSE
!         DO j=1, nvm
!            IF ( delta_veg(j) < -min_stomate ) THEN
!                dilu_bio(i,:,:)=dilu_bio(i,:,:)-delta_veg(j)*biomass(i,j,:,:)/delta_veg_sum
!            ENDIF
!         ENDDO
!         DO j=1, nvm
!            IF ( delta_veg(j) > min_stomate ) THEN
!               biomass(i,j,:,:)=(biomass(i,j,:,:)*veget_cov_max_old(i,j)+dilu_bio(i,:,:) * delta_veg(j)) / veget_cov_max(i,j)
!            ENDIF
!         ENDDO
!      ENDIF

      DO j=1, nvm
!! JC comment: Note there could be peat_grass or peat_forest, need to be changed later
!! Now there is risk of double counting if both peat_grass and peat_forest exists 
         IF (is_peat(j)) THEN !! JC comment: Note there could be peat_grass or peat_forest, need to be changed later
             IF (delta_fpeat(i) .GT. min_stomate) THEN
                bm_new(i,:,:)=delta_fpeat(i)* bm_sapl(j,:,:)
             ENDIF

             IF (veget_cov_max(i,j).GT.min_stomate) THEN
                co2_to_bm(i,j) = co2_to_bm(i,j) + (SUM(bm_new(i,:,icarbon))/ (dt_days*veget_cov_max(i,j)))
                biomass(i,j,:,:)= (biomass(i,j,:,:)* veget_cov_max_old(i,j)+ bm_new(i,:,:))/ veget_cov_max(i,j)
             ENDIF
         ELSE
            IF(veget_cov_max(i,j).GT.min_stomate) THEN
              biomass(i,j,:,:) = biomass(i,j,:,:) * veget_cov_max_old(i,j) / veget_cov_max(i,j)
            ENDIF
         ENDIF
      ENDDO



!!!!!!!!!!!!!!!!!!!!!!!!!!!!Deal with fluxes
      DO j=1, nvm
          co2flux_old(i,j)=resp_maint(i,j)+resp_growth(i,j)+resp_hetero(i,j)+co2_fire(i,j)-co2_to_bm(i,j)-gpp_daily(i,j)
          co2flux_new(i,j)=resp_maint(i,j)+resp_growth(i,j)+resp_hetero(i,j)+co2_fire(i,j)-co2_to_bm(i,j)-gpp_daily(i,j)
      ENDDO

!!!if vegetation coverage decreases, compute dilution of litter, biomass, turnover...
      DO j=1, nvm ! loop over PFTs
         IF ( delta_veg(j) < -min_stomate ) THEN
            dilu_lit(i,:,:,:) =  dilu_lit(i,:,:,:) + delta_veg(j) * litter(i,:,j,:,:) / delta_veg_sum
            dilu_f1hr(i,:,:) =  dilu_f1hr(i,:,:) + delta_veg(j) * fuel_1hr(i,j,:,:) / delta_veg_sum
            dilu_f10hr(i,:,:) =  dilu_f10hr(i,:,:) + delta_veg(j) * fuel_10hr(i,j,:,:) / delta_veg_sum
            dilu_f100hr(i,:,:) =  dilu_f100hr(i,:,:) + delta_veg(j) * fuel_100hr(i,j,:,:) / delta_veg_sum
            dilu_f1000hr(i,:,:) =  dilu_f1000hr(i,:,:) + delta_veg(j) * fuel_1000hr(i,j,:,:) / delta_veg_sum
            dilu_turnover_daily(i,:,:)=dilu_turnover_daily(i,:,:)+delta_veg(j)*turnover_daily(i,j,:,:)/delta_veg_sum
            dilu_bm_to_litter(i,:,:)=dilu_bm_to_litter(i,:,:)+delta_veg(j)*bm_to_litter(i,j,:,:)/delta_veg_sum
            dilu_co2flux_new(i)=dilu_co2flux_new(i)+delta_veg(j)*co2flux_old(i,j)/delta_veg_sum
            dilu_gpp_daily(i)=dilu_gpp_daily(i)+delta_veg(j)*gpp_daily(i,j)/delta_veg_sum
            dilu_resp_growth(i)=dilu_resp_growth(i)+delta_veg(j)*resp_growth(i,j)/delta_veg_sum
            dilu_resp_maint(i)=dilu_resp_maint(i)+delta_veg(j)*resp_maint(i,j)/delta_veg_sum
            dilu_resp_hetero(i)=dilu_resp_hetero(i)+delta_veg(j)*resp_hetero(i,j)/delta_veg_sum
            dilu_co2_to_bm(i)=dilu_co2_to_bm(i)+delta_veg(j)*co2_to_bm(i,j)/delta_veg_sum
            dilu_co2_fire(i)=dilu_co2_fire(i)+delta_veg(j)*co2_fire(i,j)/delta_veg_sum
         ENDIF
      ENDDO ! loop over PFTs

!!!PFTs that increases, recalculate the litter, biomass... with taking into accout the change in veget_cov_max
      DO j=1, nvm
         IF ( delta_veg(j) > min_stomate) THEN
           litter(i,:,j,:,:)=(litter(i,:,j,:,:) * veget_cov_max_old(i,j) + dilu_lit(i,:,:,:) * delta_veg(j)) / veget_cov_max(i,j)
           fuel_1hr(i,j,:,:)=(fuel_1hr(i,j,:,:) * veget_cov_max_old(i,j) + dilu_f1hr(i,:,:) * delta_veg(j)) / veget_cov_max(i,j)
           fuel_10hr(i,j,:,:)=(fuel_10hr(i,j,:,:) * veget_cov_max_old(i,j) + dilu_f10hr(i,:,:) * delta_veg(j)) / veget_cov_max(i,j)
           fuel_100hr(i,j,:,:)=(fuel_100hr(i,j,:,:) * veget_cov_max_old(i,j) + dilu_f100hr(i,:,:) * delta_veg(j)) / veget_cov_max(i,j)
           fuel_1000hr(i,j,:,:)=(fuel_1000hr(i,j,:,:) * veget_cov_max_old(i,j) + dilu_f1000hr(i,:,:) * delta_veg(j)) / veget_cov_max(i,j)
           IF (is_grassland_manag(j) .AND. is_grassland_grazed(j)) THEN
              litter_avail(i,:,j) = litter_avail(i,:,j) * veget_cov_max_old(i,j) / veget_cov_max(i,j)
              litter_not_avail(i,:,j) = litter(i,:,j,iabove,icarbon) - litter_avail(i,:,j)
           ENDIF
           turnover_daily(i,j,:,:)=(turnover_daily(i,j,:,:)*veget_cov_max_old(i,j)+dilu_turnover_daily(i,:,:)*delta_veg(j))/veget_cov_max(i,j)
           bm_to_litter(i,j,:,:)=(bm_to_litter(i,j,:,:)*veget_cov_max_old(i,j)+dilu_bm_to_litter(i,:,:)*delta_veg(j))/veget_cov_max(i,j)
           co2flux_new(i,j)=(co2flux_old(i,j)*veget_cov_max_old(i,j)+dilu_co2flux_new(i)*delta_veg(j))/veget_cov_max(i,j)
           gpp_daily(i,j)=(gpp_daily(i,j)*veget_cov_max_old(i,j)+dilu_gpp_daily(i)*delta_veg(j))/veget_cov_max(i,j)
           resp_growth(i,j)=(resp_growth(i,j)*veget_cov_max_old(i,j)+dilu_resp_growth(i)*delta_veg(j))/veget_cov_max(i,j)
           resp_maint(i,j)=(resp_maint(i,j)*veget_cov_max_old(i,j)+dilu_resp_maint(i)*delta_veg(j))/veget_cov_max(i,j)
           resp_hetero(i,j)=(resp_hetero(i,j)*veget_cov_max_old(i,j)+dilu_resp_hetero(i)*delta_veg(j))/veget_cov_max(i,j)
           co2_fire(i,j)=(co2_fire(i,j)*veget_cov_max_old(i,j)+dilu_co2_fire(i)*delta_veg(j))/veget_cov_max(i,j)
           co2_to_bm(i,j)=(co2_to_bm(i,j)*veget_cov_max_old(i,j)+dilu_co2_to_bm(i)*delta_veg(j))/veget_cov_max(i,j)
         ENDIF
      ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!Deal with Soil carbon pools    
      DO j=1, nvm
!!!peatland shrink, *_save record old_peat area and C (old_peat)
         IF ( delta_veg(j) .LE. zero) THEN
!! JC comment: Note there could be peat_grass or peat_forest, need to be changed later
!! Now there is risk of double counting if both peat_grass and peat_forest exists
            IF (is_peat(j)) THEN
             ! WRITE (numout,*) 'CQcheck ,peatland shrink'
             ! WRITE (numout,*) 'CQcheck ,old carbon_save',carbon_save(i,2,j),'fsave',delta_fsave(i)   
              carbon_save_tmp(i,:,j)=carbon_old(i,:,j)* (-delta_veg(j))
              carbon_save(i,:,j) = carbon_save(i,:,j) + carbon_save_tmp(i,:,j) 
              delta_fsave(i)=delta_fsave(i)-delta_veg(j)
              IF ( ok_pc ) THEN
                 deepC_a_save(i,:)=  deepC_a_save(i,:) + deepC_a(i,:,j)* (-delta_veg(j))
                 deepC_s_save(i,:)=  deepC_s_save(i,:) + deepC_s(i,:,j)* (-delta_veg(j))
                 deepC_p_save(i,:)=  deepC_p_save(i,:) + deepC_p(i,:,j)* (-delta_veg(j))
              ENDIF
             ! WRITE (numout,*) 'CQcheck ,new carbon_save',carbon_save(i,2,j),'fsave',delta_fsave(i)
            ELSE
              IF (delta_veg(j) < - min_stomate) THEN
!!!C from shrinked non-peat PFTs
                 dilu_carbon(i,:,j) =  dilu_carbon(i,:,j) - delta_veg(j) * carbon_old(i,:,j)
              ENDIF 
            ENDIF
         ENDIF
      ENDDO

      IF (delta_fpeat(i) .LE. zero) THEN
!!! peatland shrink
        DO j=1,nvm
          IF ((delta_veg(j) > min_stomate) .AND. (.NOT. is_peat(j)) .AND. (exp_nat_sum > min_stomate)) THEN
!!! expanding non-peat PFTs get C from shrinking non-peat PFTs and C of old_peat
               carbon(i,:,j)=(carbon_old(i,:,j) * veget_cov_max_old(i,j)+SUM(dilu_carbon(i,:,:),dim=2)*delta_veg(j)/exp_nat_sum+&
                             carbon_save_tmp(i,:,peat_grass)*delta_veg(j)/exp_nat_sum)/veget_cov_max(i,j)
 
               IF ( ok_pc ) THEN
                  IF (carbon_old(i,iactive,j).GT. min_stomate) THEN
                     deepC_a(i,:,j)=deepC_a(i,:,j)*carbon(i,iactive,j)/carbon_old(i,iactive,j)
                    ! WRITE (numout,*), 'CQ check deepC_a0,pft',j, deepC_a(i,:,j)
                  ENDIF
                  IF (carbon_old(i,islow,j).GT. min_stomate) THEN
                     deepC_s(i,:,j)=deepC_s(i,:,j)*carbon(i,islow,j)/carbon_old(i,islow,j)
                  ENDIF
                  IF (carbon_old(i,ipassive,j).GT. min_stomate) THEN
                     deepC_p(i,:,j)=deepC_p(i,:,j)*carbon(i,ipassive,j)/carbon_old(i,ipassive,j)
                  ENDIF
               ENDIF
          ENDIF    
        ENDDO
      ENDIF


!!! peatland expand, occupy old_peat first
      IF (delta_fpeat(i) .GT. zero) THEN
         IF (delta_fpeat(i) .LE. delta_fsave(i)) THEN
!!!if there is old_peat and old_peat is large enough for peatland expansion, peatland obtain C only from old_peat
            carbon_obtain(i,:)= MIN(un,delta_fpeat(i)/delta_fsave(i))*carbon_save(i,:,peat_grass)
            carbon_save(i,:,peat_grass)=MAX((un-delta_fpeat(i)/delta_fsave(i)),zero)*carbon_save(i,:,peat_grass)
            delta_fsave(i)=MAX(delta_fsave(i)-delta_fpeat(i),zero)
            IF ( ok_pc ) THEN
               deepC_a_save(i,:)=MAX((un-delta_fpeat(i)/delta_fsave(i)),zero)*deepC_a_save(i,:)
               deepC_s_save(i,:)=MAX((un-delta_fpeat(i)/delta_fsave(i)),zero)*deepC_s_save(i,:)
               deepC_p_save(i,:)=MAX((un-delta_fpeat(i)/delta_fsave(i)),zero)*deepC_p_save(i,:)
            ENDIF
         ELSE
!!!if there is no old_peat or old_peat is not large enough, get old_peat first
!then get some c from mineral soil
            diff_fpeat(i)=delta_fpeat(i)-delta_fsave(i)
            nat_sum=zero
            DO j=1,nvm
               IF (natural(j) .AND. (veget_cov_max_old(i,j) .GT. min_stomate) .AND. .NOT. is_peat(j) ) THEN
                  nat_sum=nat_sum+veget_cov_max_old(i,j)
               ENDIF
            ENDDO
            DO j=1,nvm
               IF (natural(j) .AND. (veget_cov_max_old(i,j) .GT. min_stomate) .AND. .NOT. is_peat(j)) THEN
                  carbon_obtain(i,:)=carbon_obtain(i,:)+diff_fpeat(i)*(veget_cov_max_old(i,j)/nat_sum)*carbon_old(i,:,j)
               ENDIF
            ENDDO
            carbon_obtain(i,:)= carbon_obtain(i,:)+carbon_save(i,:,peat_grass) 
            carbon_save(i,:,peat_grass)=zero
            delta_fsave(i)= zero
            IF (ok_pc) THEN
               deepC_a_save(i,:)=zero
               deepC_s_save(i,:)=zero
               deepC_p_save(i,:)=zero
            ENDIF
         ENDIF

!!!add *_obtain to peatland 
         DO j=1,nvm
            IF (is_peat(j)) THEN
                carbon(i,:,j)=(carbon_old(i,:,j)*veget_cov_max_old(i,j)+ carbon_obtain(i,:))/veget_cov_max(i,j)
                IF ( ok_pc ) THEN
                   IF (carbon_old(i,iactive,j).GT. min_stomate) THEN
                      deepC_a(i,:,j)=deepC_a(i,:,j)*carbon(i,iactive,j)/carbon_old(i,iactive,j)
                   ENDIF
                   IF (carbon_old(i,islow,j).GT. min_stomate) THEN
                      deepC_s(i,:,j)=deepC_s(i,:,j)*carbon(i,islow,j)/carbon_old(i,islow,j)
                   ENDIF
                   IF (carbon_old(i,ipassive,j).GT. min_stomate) THEN
                      deepC_p(i,:,j)=deepC_p(i,:,j)*carbon(i,ipassive,j)/carbon_old(i,ipassive,j)
                   ENDIF
                ENDIF
            ENDIF
         ENDDO

         carbon_tmp(i,:,:) = carbon(i,:,:)
!!!substract *_obtain from non-peat PFTs
         DO l=1, ncarb 
   !! if contracting non-peat PFTs can provide the *_obtain
            IF (SUM(dilu_carbon(i,l,:)) .GE. carbon_obtain(i,l)) THEN
        !! subtract *_obtain (has been given to peatland) from dilu_carbon,and then give the remaining dilu_carbon to expanding non-peat PFTs
               delta_nonpeat_tmp=zero
               DO j=1,nvm
                 IF (.NOT. is_peat(j) .AND. .NOT. is_croppeat(j)) THEN
                   delta_nonpeat_tmp=delta_nonpeat_tmp+delta_veg(j)
                 ENDIF 
               ENDDO
               IF (delta_nonpeat_tmp .GT. min_stomate) THEN
                  DO j=1,nvm
                    IF ((delta_veg(j) > min_stomate) .AND. (natural(j)) .AND. (exp_nat_sum > min_stomate) .AND. (.NOT. is_peat(j))) THEN
                        carbon(i,l,j)=(carbon_old(i,l,j) * veget_cov_max_old(i,j)+ &
                             (SUM(dilu_carbon(i,l,:))-carbon_obtain(i,l))*delta_veg(j)/exp_nat_sum)/veget_cov_max(i,j)
               !     WRITE (numout,*) 'CQ check, C substracted from nonpeat1',SUM(dilu_carbon(i,l,:))-carbon_obtain(i,l)
                      IF ( ok_pc ) THEN
                          IF ( (l .EQ. iactive) .AND. (carbon_old(i,iactive,j).GT. min_stomate)) THEN
                             deepC_a(i,:,j)=deepC_a(i,:,j)*carbon(i,iactive,j)/carbon_old(i,iactive,j)
                        !   WRITE (numout,*), 'CQ check deepC_a2,pft',j,deepC_a(i,:,j)
                          ENDIF
                          IF ( (l .EQ. islow) .AND. (carbon_old(i,islow,j).GT. min_stomate)) THEN
                             deepC_s(i,:,j)=deepC_s(i,:,j)*carbon(i,islow,j)/carbon_old(i,islow,j)
                          ENDIF
                          IF ( (l .EQ. ipassive) .AND. (carbon_old(i,ipassive,j).GT. min_stomate)) THEN
                             deepC_p(i,:,j)=deepC_p(i,:,j)*carbon(i,ipassive,j)/carbon_old(i,ipassive,j)
                          ENDIF
                      ENDIF
                    ENDIF
                  ENDDO
               ELSE  
        !! peatland is the only PFT that is expanding, add all dilu_carbon to peatland
                  DO j=1,nvm
                     IF (is_peat(j)) THEN
                        carbon(i,l,j)= carbon(i,l,j)+(SUM(dilu_carbon(i,l,:))-carbon_obtain(i,l))/veget_cov_max(i,j)
                        IF ( ok_pc ) THEN
                           IF ((carbon_tmp(i,iactive,j).GT. min_stomate).AND. (l .EQ. iactive)) THEN
                              deepC_a(i,:,j)=deepC_a(i,:,j)*carbon(i,iactive,j)/carbon_tmp(i,iactive,j)
                           ENDIF
                           IF ((carbon_tmp(i,islow,j).GT. min_stomate).AND. (l .EQ. islow)) THEN
                              deepC_s(i,:,j)=deepC_s(i,:,j)*carbon(i,islow,j)/carbon_tmp(i,islow,j)
                           ENDIF
                           IF ((carbon_tmp(i,ipassive,j).GT. min_stomate) .AND. (l .EQ. ipassive)) THEN
                              deepC_p(i,:,j)=deepC_p(i,:,j)*carbon(i,ipassive,j)/carbon_tmp(i,ipassive,j)
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF
            ELSE
   !! If contracting non-peat PFTs can't provide enough C to *_obtain, all dilu_carbon will be added to peat
    !!  expanding non-peat PFTs can't receive any C, adjust their C density
               DO j=1,nvm
                  IF ((delta_veg(j) > min_stomate) .AND. (natural(j)) .AND. (.NOT. is_peat(j))) THEN
                     carbon(i,l,j)=(carbon_old(i,l,j) * veget_cov_max_old(i,j))/veget_cov_max(i,j)
                 ENDIF
               ENDDO
     !! Then the difference between dilu_carbon and *_obtain need to be substracted from non-peat PFTs  
               excess_C(i,l)=SUM(dilu_carbon(i,l,:))-carbon_obtain(i,l)
               nat_sum=zero
               DO j=1,nvm
                 IF (natural(j) .AND. (veget_cov_max(i,j) .GT. min_stomate) .AND. (.NOT. is_peat(j))) THEN
                    nat_sum=nat_sum+veget_cov_max(i,j)*carbon(i,l,j)
                 ENDIF
               ENDDO
               DO j=1,nvm
                  IF ((natural(j) .AND. .NOT. is_peat(j) ).AND. (veget_cov_max(i,j) .GT. min_stomate) .AND. (excess_C(i,l) .LT. zero) .AND. (nat_sum .GT. zero)) THEN  
                     excess_tmp=MAX(-carbon(i,l,j) *veget_cov_max(i,j),excess_C(i,l)*veget_cov_max(i,j)*carbon(i,l,j)/nat_sum)
                     nat_sum=nat_sum-veget_cov_max(i,j)*carbon(i,l,j)                     
                     carbon(i,l,j)=(carbon(i,l,j)*veget_cov_max(i,j)+ excess_tmp)/veget_cov_max(i,j)
                     excess_C(i,l)=excess_C(i,l)-excess_tmp

                     IF (ok_pc) THEN
                        IF ( (l .EQ. iactive) .AND. (carbon_old(i,iactive,j).GT. min_stomate)) THEN
                           deepC_a(i,:,j)=deepC_a(i,:,j)*carbon(i,iactive,j)/carbon_old(i,iactive,j)
                         !  WRITE (numout,*), 'CQ check deepC_a3,pft',j, deepC_a(i,:,j)
                        ENDIF
                        IF ( (l .EQ. islow) .AND. (carbon_old(i,islow,j).GT. min_stomate)) THEN
                           deepC_s(i,:,j)=deepC_s(i,:,j)*carbon(i,islow,j)/carbon_old(i,islow,j)
                        ENDIF
                        IF ( (l .EQ. ipassive) .AND. (carbon_old(i,ipassive,j).GT. min_stomate)) THEN
                           deepC_p(i,:,j)=deepC_p(i,:,j)*carbon(i,ipassive,j)/carbon_old(i,ipassive,j)
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
               IF (excess_C(i,l) .LT. zero) THEN !!!!Sum of C from all mineral soils is smaller than old_peat C
                  DO j=1,nvm
                     IF (is_peat(j)) THEN
                        carbon(i,l,j)=(carbon(i,l,j)*veget_cov_max(i,j)+ excess_C(i,l))/veget_cov_max(i,j)
                        IF ( ok_pc ) THEN
                           IF ((carbon_tmp(i,iactive,j).GT. min_stomate) .AND. (l .EQ. iactive)) THEN
                              deepC_a(i,:,j)=deepC_a(i,:,j)*carbon(i,iactive,j)/carbon_tmp(i,iactive,j)
                           ENDIF
                           IF ((carbon_tmp(i,islow,j).GT. min_stomate) .AND. (l .EQ. islow)) THEN
                              deepC_s(i,:,j)=deepC_s(i,:,j)*carbon(i,islow,j)/carbon_tmp(i,islow,j)
                           ENDIF
                           IF ((carbon_tmp(i,ipassive,j).GT. min_stomate) .AND. (l .EQ. ipassive)) THEN
                              deepC_p(i,:,j)=deepC_p(i,:,j)*carbon(i,ipassive,j)/carbon_tmp(i,ipassive,j)
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF 
            ENDIF
         ENDDO !(ncarb)
      ENDIF !(delta_fpeat(i) > min_stomate)
    ENDDO !(npts)


       DO i=1, npts
          DO l=1,ncarb
             DO j=1,nvm
                soilc_after(i)=soilc_after(i)+carbon(i,l,j)*veget_cov_max(i,j)
             ENDDO
          ENDDO
       ENDDO


       DO i=1, npts
         ! IF (ABS(biomass_after(i)-biomass_before(i)) .GT. 1000.*min_stomate) THEN
         ! ENDIF
          IF ( ABS(soilc_after(i)-soilc_before(i)) .GT. 1000.*min_stomate) THEN
             WRITE(numout,*) ' qcj check lpj_cover_peat,C',lalo(i,:)
             WRITE(numout,*) ' qcj check C,after-before',soilc_after(i),soilc_before(i)
          ENDIF
       ENDDO

  END SUBROUTINE lpj_cover_peat
!_
!================================================================================================================================
END MODULE stomate_lpj
