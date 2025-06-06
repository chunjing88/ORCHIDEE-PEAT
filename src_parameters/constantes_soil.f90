! =================================================================================================================================
! MODULE 	: constantes_soil
!
! CONTACT       : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE       : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         "constantes_soil" module contains subroutine to initialize the parameters related to soil and hydrology.
!!
!!\n DESCRIPTION : "constantes_soil" module contains subroutine to initialize the parameters related to soil and hydrology.
!!                 This module alos USE constates_soil and can therfor be used to acces the subroutines and the constantes.
!!                 The constantes declarations can also be used seperatly with "USE constantes_soil_var".
!!
!! RECENT CHANGE(S): 
!!
!! REFERENCE(S)	:
!!
!! SVN          :
!! $HeadURL: $
!! $Date: $
!! $Revision: $
!! \n
!_ ================================================================================================================================

MODULE constantes_soil

  USE constantes_soil_var
  USE ioipsl_para 

  IMPLICIT NONE

CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : config_soil_parameters
!!
!>\BRIEF        This subroutine reads in the configuration file all the parameters related to soil and hydrology. 
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) :
!!
!! FLOWCHART    :
!! \n
!_ ================================================================================================================================

  SUBROUTINE config_soil_parameters()

    USE ioipsl

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.4 Local variables 

    INTEGER(i_std), PARAMETER      :: error_level = 3         !! Switch to 2 to turn fatal errors into warnings.(1-3, unitless)
    LOGICAL                        :: ok_freeze               !! Local variable used to set default values for all flags 
    !! controling the soil freezing scheme
    !_ ================================================================================================================================

    ! Following initializations are only done for option impose_param
    IF ( ok_sechiba .AND. impose_param ) THEN

       !Config Key   = DRY_SOIL_HEAT_CAPACITY
       !Config Desc  = Dry soil Heat capacity of soils
       !Config If    = OK_SECHIBA 
       !Config Def   = 1.80e+6
       !Config Help  = Values taken from : PIELKE,'MESOSCALE METEOROLOGICAL MODELING',P.384.
       !Config Units = [J.m^{-3}.K^{-1}] 
       CALL getin_p("DRY_SOIL_HEAT_CAPACITY",so_capa_dry)

       !! Check parameter value (correct range)
       IF ( so_capa_dry <= zero ) THEN
          CALL ipslerr_p(error_level, "config_soil_parameters.", &
               &     "Wrong parameter value for DRY_SOIL_HEAT_CAPACITY.", &
               &     "This parameter should be positive. ", &
               &     "Please, check parameter value in run.def. ")
       END IF


       !Config Key   = DRY_SOIL_HEAT_COND
       !Config Desc  = Dry soil Thermal Conductivity of soils
       !Config If    = OK_SECHIBA
       !Config Def   = 0.40 
       !Config Help  = Values taken from : PIELKE,'MESOSCALE METEOROLOGICAL MODELING',P.384.
       !Config Units = [W.m^{-2}.K^{-1}] 
       CALL getin_p("DRY_SOIL_HEAT_COND",so_cond_dry)

       !! Check parameter value (correct range)
       IF ( so_cond_dry <= zero ) THEN
          CALL ipslerr_p(error_level, "config_soil_parameters.", &
               &     "Wrong parameter value for DRY_SOIL_HEAT_COND.", &
               &     "This parameter should be positive. ", &
               &     "Please, check parameter value in run.def. ")
       END IF


       !Config Key   = WET_SOIL_HEAT_CAPACITY
       !Config Desc  = Wet soil Heat capacity of soils 
       !Config If    = OK_SECHIBA
       !Config Def   = 3.03e+6
       !Config Help  = 
       !Config Units = [J.m^{-3}.K^{-1}]
       CALL getin_p("WET_SOIL_HEAT_CAPACITY",so_capa_wet)

       !! Check parameter value (correct range)
       IF ( so_capa_wet <= zero ) THEN
          CALL ipslerr_p(error_level, "config_soil_parameters.", &
               &     "Wrong parameter value for WET_SOIL_HEAT_CAPACITY.", &
               &     "This parameter should be positive. ", &
               &     "Please, check parameter value in run.def. ")
       END IF


       !Config Key   = WET_SOIL_HEAT_COND
       !Config Desc  = Wet soil Thermal Conductivity of soils
       !Config If    = OK_SECHIBA 
       !Config Def   = 1.89 
       !Config Help  = 
       !Config Units = [W.m^{-2}.K^{-1}]
       CALL getin_p("WET_SOIL_HEAT_COND",so_cond_wet)

       !! Check parameter value (correct range)
       IF ( so_cond_wet <= zero ) THEN
          CALL ipslerr_p(error_level, "config_soil_parameters.", &
               &     "Wrong parameter value for WET_SOIL_HEAT_COND.", &
               &     "This parameter should be positive. ", &
               &     "Please, check parameter value in run.def. ")
       END IF


       !Config Key   = SNOW_HEAT_COND
       !Config Desc  = Thermal Conductivity of snow
       !Config If    = OK_SECHIBA  
       !Config Def   = 0.3
       !Config Help  = 
       !Config Units = [W.m^{-2}.K^{-1}]
       CALL getin_p("SNOW_HEAT_COND",sn_cond)

       !! Check
       IF ( sn_cond <= zero ) THEN
          CALL ipslerr_p(error_level, "config_soil_parameters.", &
               &     "Wrong parameter value for SNOW_HEAT_COND.", &
               &     "This parameter should be positive. ", &
               &     "Please, check parameter value in run.def. ")
       END IF


       !Config Key   = SNOW_DENSITY
       !Config Desc  = Snow density for the soil thermodynamics 
       !Config If    = OK_SECHIBA 
       !Config Def   = 330.0
       !Config Help  = 
       !Config Units = [-] 
       CALL getin_p("SNOW_DENSITY",sn_dens)

       !! Check parameter value (correct range)
       IF ( sn_dens <= zero ) THEN
          CALL ipslerr_p(error_level, "config_soil_parameters.", &
               &     "Wrong parameter value for SNOW_DENSITY.", &
               &     "This parameter should be positive. ", &
               &     "Please, check parameter value in run.def. ")
       END IF


       !! Calculation of snow capacity
       !! If sn_dens is redefined by the user, sn_capa needs to be reset
       sn_capa = 2100.0_r_std*sn_dens


       !Config Key   = NOBIO_WATER_CAPAC_VOLUMETRI
       !Config Desc  = 
       !Config If    = OK_SECHIBA and .NOT.(HYDROL_CWRR)
       !Config Def   = 150.
       !Config Help  = 
       !Config Units = [s/m^2]
       CALL getin_p('NOBIO_WATER_CAPAC_VOLUMETRI',mx_eau_nobio)

       !! Check parameter value (correct range)
       IF ( mx_eau_nobio <= zero ) THEN
          CALL ipslerr_p(error_level, "config_soil_parameters.", &
               &     "Wrong parameter value for NOBIO_WATER_CAPAC_VOLUMETRI.", &
               &     "This parameter should be positive. ", &
               &     "Please, check parameter value in run.def. ")
       END IF


       !Config Key   = SECHIBA_QSINT 
       !Config Desc  = Interception reservoir coefficient
       !Config If    = OK_SECHIBA 
       !Config Def   = 0.1
       !Config Help  = Transforms leaf area index into size of interception reservoir
       !Config         for slowproc_derivvar or stomate
       !Config Units = [m]
       CALL getin_p('SECHIBA_QSINT',qsintcst)

       !! Check parameter value (correct range)
       IF ( qsintcst <= zero ) THEN
          CALL ipslerr_p(error_level, "config_soil_parameters.", &
               &     "Wrong parameter value for SECHIBA_QSINT.", &
               &     "This parameter should be positive. ", &
               &     "Please, check parameter value in run.def. ")
       END IF

        !Config Key  = SOIL_LAYERS_DISCRE_METHOD
        !Config Desc = Select which soil layer discretization method use
        !Config If   = 
        !Config Def  = 0 (Thermix method)
        !Config Help = 0 = thermix, 1 = permafrost, any other value is not valid 
        !Config Units= [FLAG]
        SO_DISCRETIZATION_METHOD = 0
        CALL getin_p('SOIL_LAYERS_DISCRE_METHOD', SO_DISCRETIZATION_METHOD)
        !! Check parameter value (correct range)
        IF ( SO_DISCRETIZATION_METHOD < zero .OR. SO_DISCRETIZATION_METHOD > deux ) THEN
           CALL ipslerr_p(error_level, "config_soil_parameters.", &
                &     "Wrong parameter value for SOIL_LAYERS_DISCRE_METHOD.", &
                &     "Use 0 for thermix or 1 for permafrost method ", &
                &     "Please, check parameter value in run.def. ")
        END IF

       IF ( .NOT.(hydrol_cwrr) ) THEN

          !Config Key   = CHOISNEL_DIFF_MIN
          !Config Desc  = Diffusion constant for the slow regime
          !Config If    = OK_SECHIBA and .NOT.(HYDROL_CWRR)
          !Config Def   = 0.001
          !Config Help  = 
          !Config Units = [kg/m^2/dt]
          CALL getin_p('CHOISNEL_DIFF_MIN',min_drain)

          !! Check parameter value (correct range)
          IF ( min_drain <= zero ) THEN
             CALL ipslerr_p(error_level, "config_soil_parameters.", &
                  &     "Wrong parameter value for CHOISNEL_DIFF_MIN.", &
                  &     "This parameter should be positive. ", &
                  &     "Please, check parameter value in run.def. ")
          END IF


          !Config Key   = CHOISNEL_DIFF_MAX
          !Config Desc  = Diffusion constant for the fast regime
          !Config If    = OK_SECHIBA and .NOT.(HYDROL_CWRR)
          !Config Def   = 0.1
          !Config Help  = 
          !Config Units = [kg/m^2/dt]
          CALL getin_p('CHOISNEL_DIFF_MAX',max_drain)

          !! Check parameter value (correct range)
          IF (  ( max_drain <= zero ) .OR. ( max_drain <= min_drain ) ) THEN
             CALL ipslerr_p(error_level, "config_soil_parameters.", &
                  &     "Wrong parameter value for CHOISNEL_DIFF_MAX.", &
                  &     "This parameter should be positive or greater than CHOISNEL_DIFF_MIN.", &
                  &     "Please, check parameter value in run.def. ")
          END IF


          !Config Key   = CHOISNEL_DIFF_EXP
          !Config Desc  = The exponential in the diffusion law
          !Config If    = OK_SECHIBA and .NOT.(HYDROL_CWRR)
          !Config Def   = 1.5
          !Config Help  = 
          !Config Units = [-]
          CALL getin_p('CHOISNEL_DIFF_EXP',exp_drain)

          !! Check parameter value (correct range)
          IF ( exp_drain <= zero ) THEN
             CALL ipslerr_p(error_level, "config_soil_parameters.", &
                  &     "Wrong parameter value for CHOISNEL_DIFF_EXP.", &
                  &     "This parameter should be positive. ", &
                  &     "Please, check parameter value in run.def. ")
          END IF


          !Config Key   = CHOISNEL_RSOL_CSTE
          !Config Desc  = Constant in the computation of resistance for bare  soil evaporation 
          !Config If    = OK_SECHIBA and .NOT.(HYDROL_CWRR)
          !Config Def   = 33.E3
          !Config Help  = 
          !Config Units = [s/m^2]
          CALL getin_p('CHOISNEL_RSOL_CSTE',rsol_cste)

          !! Check parameter value (correct range)
          IF ( rsol_cste <= zero ) THEN
             CALL ipslerr_p(error_level, "config_soil_parameters.", &
                  &     "Wrong parameter value for CHOISNEL_RSOL_CSTE.", &
                  &     "This parameter should be positive. ", &
                  &     "Please, check parameter value in run.def. ")
          END IF


          !Config Key   = HCRIT_LITTER
          !Config Desc  = Scaling depth for litter humidity
          !Config If    = OK_SECHIBA and .NOT.(HYDROL_CWRR) 
          !Config Def   = 0.08 
          !Config Help  = 
          !Config Units = [m]
          CALL getin_p('HCRIT_LITTER',hcrit_litter)

          !! Check parameter value (correct range)
          IF ( hcrit_litter <= zero ) THEN
             CALL ipslerr_p(error_level, "config_soil_parameters.", &
                  &     "Wrong parameter value for HCRIT_LITTER.", &
                  &     "This parameter should be positive. ", &
                  &     "Please, check parameter value in run.def. ")
          END IF

       END IF

    END IF ! IF ( ok_sechiba .AND. impose_param ) THEN

!!!qcj++ peatland
    CALL getin_p("TAU_PEAT_GRASS",tau_peat_grass)
    CALL getin_p("TAU_PEAT_SHRUB",tau_peat_shrub)
    CALL getin_p("TAU_PEAT_MOSS",tau_peat_moss)

    CALL getin_p("Z_TAU",z_tau)
    CALL getin_p("TAU_agrPEAT",tau_agrpeat)
    !! Variables related to soil freezing in thermosoil module
    !
    !Config Key  = OK_FREEZE
    !Config Desc = Activate the complet soil freezing scheme
    !Config If   = OK_SECHIBA 
    !Config Def  = FALSE
    !Config Help = Activate soil freezing thermal effects. Activates soil freezing hydrological effects in CWRR scheme.
    !Config Units= [FLAG]

    ! ok_freeze is a flag that controls the default values for several flags controling 
    ! the different soil freezing processes
    ! Set ok_freeze=true for the complete soil freezing scheme
    ! ok_freeze is a local variable only used in this subroutine
    ok_freeze = .FALSE.
    CALL getin_p('OK_FREEZE',ok_freeze)


    !Config Key  = READ_REFTEMP
    !Config Desc = Initialize soil temperature using climatological temperature
    !Config If   = 
    !Config Def  = True/False depening on OK_FREEZE
    !Config Help = 
    !Config Units= [FLAG]

    IF (ok_freeze) THEN
       read_reftemp = .TRUE.
    ELSE
       read_reftemp = .FALSE.
    END IF
    CALL getin_p ('READ_REFTEMP',read_reftemp)

    !Config Key  = OK_FREEZE_THERMIX
    !Config Desc = Activate thermal part of the soil freezing scheme
    !Config If   = 
    !Config Def  = True if OK_FREEZE else false
    !Config Help = 
    !Config Units= [FLAG]
    IF (ok_freeze) THEN
       ok_freeze_thermix = .TRUE.
    ELSE
       ok_freeze_thermix = .FALSE.
    END IF
    CALL getin_p ('OK_FREEZE_THERMIX',ok_freeze_thermix)


    !Config Key  = OK_ECORR
    !Config Desc = Energy correction for freezing
    !Config If   = OK_FREEZE_THERMIX
    !Config Def  = True if OK_FREEZE else false
    !Config Help = Energy conservation : Correction to make sure that the same latent heat is 
    !Config        released and consumed during freezing and thawing
    !Config Units= [FLAG]
    IF (ok_freeze) THEN
       ok_Ecorr = .TRUE.
    ELSE
       ok_Ecorr = .FALSE.
    END IF
    CALL getin_p ('OK_ECORR',ok_Ecorr)
    IF (ok_Ecorr .AND. .NOT. ok_freeze_thermix) THEN
       CALL ipslerr_p(3,'thermosoil_init','OK_ECORR cannot be activated without OK_FREEZE_THERMIX', &
            'Adapt run parameters with OK_FREEZE_THERMIX=y','')
    END IF

    !Config Key = POROS
    !Config Desc = Soil porosity 
    !Config If = OK_SECHIBA
    !Config Def = 0.41
    !Config Help = From USDA classification, mean value
    !Config Units = [-] 
    poros=0.41
    CALL getin_p('POROS',poros)


    !Config Key = fr_dT
    !Config Desc = Freezing window    
    !Config If = OK_SECHIBA
    !Config Def = 2.0
    !Config Help = 
    !Config Units = [K] 
    fr_dT=2.0
    CALL getin_p('FR_DT',fr_dT)


    !Config Key  = QZ_USDA
    !Config Desc = quartz content
    !Config If   = 
    !Config Def  = 0.92, 0.82, 0.60, 0.25, 0.10, 0.40, 0.60, 0.10, 0.35, 0.52, 0.10, 0.25
    !Config Help =
    !Config Units= [-]
    CALL getin_p('QZ_USDA', QZ_usda)

    !Config Key   = SOILC_MAX
    !Config Desc  = soil carbon above which soil thermal properties equals to organic soil properties 
    !Config If    = 
    !Config Def   = 130000
    !Config Help  = 
    !Config Units = [gC/m3] 
    CALL getin_p("SOILC_MAX", soilc_max)

    !! Variables related to soil Freezing in hydrol module

    !Config Key  = SMCMAX_FAO 
    !Config Desc = Fao Porosity  
    !Config If   = 
    !Config Def  = 0.41_r_std, 0.43_r_std, 0.41_r_std
    !Config Help =
    !Config Units= [FLAG]
    CALL getin_p('SMCMAX_FAO', SMCMAX_fao)

    !Config Key  = OK_FREEZE_CWRR
    !Config Desc = CWRR freezing scheme by I. Gouttevin
    !Config If   = 
    !Config Def  = True if OK_FREEZE else false
    !Config Help =
    !Config Units= [FLAG]
    IF (ok_freeze) THEN
       ok_freeze_cwrr = .TRUE.
    ELSE
       ok_freeze_cwrr = .FALSE.
    END IF
    CALL getin_p('OK_FREEZE_CWRR',ok_freeze_cwrr)


    IF (ok_freeze_cwrr) THEN
       !Config Key  = OK_THERMODYNAMICAL_FREEZING
       !Config Desc = Calculate frozen fraction thermodynamically 
       !Config If   = HYDROL_CWRR .AND. OK_FREEZE_CWRR
       !Config Def  = True
       !Config Help = Calculate frozen fraction thermodynamically if true,
       !Config      = else calculate frozen fraction linearly 
       !Config Units= [FLAG]
       ok_thermodynamical_freezing = .TRUE.
       CALL getin_p('OK_THERMODYNAMICAL_FREEZING',ok_thermodynamical_freezing)
    END IF


    !! 1 Some initializations
    !
    !
    !Config Key   = CHECK_CWRR
    !Config Desc  = Check detailed CWRR water balance
    !Config Def   = n
    !Config If    = HYDROL_CWRR
    !Config Help  = This parameters allows the user to check
    !Config         the detailed water balance in each time step
    !Config         of CWRR and stop execution if not correct
    !Config Units = [FLAG]
    !
    check_cwrr = .FALSE.
    CALL getin_p('CHECK_CWRR', check_cwrr)

    !Config Key   = CHECK_CWRR2
    !Config Desc  = Caluculate diagnostics to check CWRR water balance
    !Config Def   = n
    !Config If    = HYDROL_CWRR2
    !Config Help  = The verifictaions are done in post-treatement
    !Config Units = [FLAG]
    !
    check_cwrr2 = .FALSE.
    CALL getin_p('CHECK_CWRR2', check_cwrr2)

  END SUBROUTINE config_soil_parameters


END MODULE constantes_soil
