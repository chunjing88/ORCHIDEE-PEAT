! ================================================================================================================================
!  MODULE       : diffuco
!
!  CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
!  LICENCE      : IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF   This module calculates the limiting coefficients, both aerodynamic
!! and hydrological, for the turbulent heat fluxes.
!!
!!\n DESCRIPTION: The aerodynamic resistance R_a is used to limit
!! the transport of fluxes from the surface layer of vegetation to the point in the atmosphere at which
!! interaction with the LMDZ atmospheric circulation model takes place. The aerodynamic resistance is
!! calculated either within the module r_aerod (if the surface drag coefficient is provided by the LMDZ, and 
!! if the flag 'ldq_cdrag_from_gcm' is set to TRUE) or r_aero (if the surface drag coefficient must be calculated).\n
!!
!! Within ORCHIDEE, evapotranspiration is a function of the Evaporation Potential, but is modulated by a
!! series of resistances (canopy and aerodynamic) of the surface layer, here represented by beta.\n
!!
!! DESCRIPTION	:
!! \latexonly 
!!     \input{diffuco_intro.tex}
!! \endlatexonly
!! \n
!!
!! This module calculates the beta for several different scenarios: 
!! - diffuco_snow calculates the beta coefficient for sublimation by snow, 
!! - diffuco_inter calculates the beta coefficient for interception loss by each type of vegetation, 
!! - diffuco_bare calculates the beta coefficient for bare soil, 
!! - diffuco_trans or diffuco_trans_co2 both calculate the beta coefficient for transpiration for each type
!!   of vegetation. Those routines differ by the formulation used to calculate the canopy resistance (Jarvis in 
!!   diffuco_trans, Farqhar in diffuco_trans_co2)
!! - chemistry_bvoc calculates the beta coefficient for emissions of biogenic compounds \n
!!
!! Finally, the module diffuco_comb computes the combined $\alpha$ and $\beta$ coefficients for use 
!! elsewhere in the module. \n

!! RECENT CHANGE(S): Nathalie le 28 mars 2006 - sur proposition de Fred Hourdin, ajout
!! d'un potentiometre pour regler la resistance de la vegetation (rveg is now in pft_parameters)
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_sechiba/diffuco.f90 $
!! $Date: 2018-04-11 10:59:44 +0200 (Wed, 11 Apr 2018) $
!! $Revision: 5188 $
!! \n
!_ ================================================================================================================================

MODULE diffuco

  ! modules used :
  USE constantes
  USE constantes_soil
  USE qsat_moisture
  USE sechiba_io_p
  USE ioipsl
  USE pft_parameters
  USE grid
  USE time, ONLY : one_day, dt_sechiba
  USE ioipsl_para 
  USE xios_orchidee
  USE chemistry, ONLY : chemistry_initialize, chemistry_bvoc, chemistry_clear
  IMPLICIT NONE

  ! public routines :
  PRIVATE
  PUBLIC :: diffuco_main, diffuco_initialize, diffuco_finalize, diffuco_clear

  INTERFACE Arrhenius_modified
     MODULE PROCEDURE Arrhenius_modified_0d, Arrhenius_modified_1d
  END INTERFACE

  !
  ! variables used inside diffuco module : declaration and initialisation
  !
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: wind                      !! Wind module (m s^{-1})
!$OMP THREADPRIVATE(wind)

CONTAINS


!!  =============================================================================================================================
!! SUBROUTINE:    diffuco_initialize
!!
!>\BRIEF	  Allocate module variables, read from restart file or initialize with default values
!!
!! DESCRIPTION:	  Allocate module variables, read from restart file or initialize with default values.
!!                Call chemistry_initialize for initialization of variables needed for the calculations of BVOCs.
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S): None
!! 
!! FLOWCHART: None
!! \n
!_ ==============================================================================================================================
  SUBROUTINE diffuco_initialize (kjit,    kjpindex, index,                  &
                                 rest_id, lalo,     neighbours, resolution, &
                                 rstruct, q_cdrag,  q_cdrag_pft)
    
    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                         :: kjit             !! Time step number (-) 
    INTEGER(i_std), INTENT(in)                         :: kjpindex         !! Domain size (-)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: index            !! Indeces of the points on the map (-)
    INTEGER(i_std),INTENT (in)                         :: rest_id          !! _Restart_ file identifier (-)
    REAL(r_std),DIMENSION (kjpindex,2),   INTENT (in)  :: lalo             !! Geographical coordinates
    INTEGER(i_std),DIMENSION (kjpindex,NbNeighb),INTENT (in):: neighbours  !! Vector of neighbours for each 
    REAL(r_std),DIMENSION (kjpindex,2), INTENT(in)     :: resolution       !! The size in km of each grid-box in X and Y
    
    !! 0.2 Output variables
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: rstruct          !! Structural resistance for the vegetation
    
    !! 0.3 Modified variables
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: q_cdrag          !! Surface drag coefficient  (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: q_cdrag_pft      !! pft specific cdrag
    
    !! 0.4 Local variables
    INTEGER                                            :: ilai
    INTEGER                                            :: jv
    INTEGER                                            :: ier
    CHARACTER(LEN=4)                                   :: laistring
    CHARACTER(LEN=80)                                  :: var_name        
    REAL(r_std),DIMENSION (kjpindex,nvm)              :: temppft
    !_ ================================================================================================================================
    
    !! 1. Define flag ldq_cdrag_from_gcm. This flag determines if the cdrag should be taken from the GCM or be calculated. 
    !!    The default value is true if the q_cdrag variables was already initialized. This is the case when coupled to the LMDZ.

    !Config Key   = CDRAG_FROM_GCM
    !Config Desc  = Keep cdrag coefficient from gcm.
    !Config If    = OK_SECHIBA
    !Config Def   = y
    !Config Help  = Set to .TRUE. if you want q_cdrag coming from GCM (if q_cdrag on initialization is non zero).
    !Config         Keep cdrag coefficient from gcm for latent and sensible heat fluxes.
    !Config Units = [FLAG]
    IF ( ABS(MAXVAL(q_cdrag)) .LE. EPSILON(q_cdrag)) THEN
       ldq_cdrag_from_gcm = .FALSE.
    ELSE
       ldq_cdrag_from_gcm = .TRUE.
    ENDIF
    CALL getin_p('CDRAG_from_GCM', ldq_cdrag_from_gcm)
    IF (printlev>=2) WRITE(numout,*) "ldq_cdrag_from_gcm = ",ldq_cdrag_from_gcm

    IF (printlev>=2) WRITE(numout,*) 'getting restart of q_cdrag_pft'
    var_name = 'cdrag_pft'
    q_cdrag_pft = zero
    CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., temppft, "gather", nbp_glo, index_g)
    IF (MINVAL(temppft) < MAXVAL(temppft) .OR. MAXVAL(temppft) .NE. val_exp) THEN
        q_cdrag_pft(:,:) = temppft(:,:)
    ENDIF
    IF (printlev>=2) WRITE(numout,*) 'done restart of q_cdrag_pft'

    !! 2. Allocate module variables
    ALLOCATE (wind(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'diffuco_initialize','Problem in allocate of variable wind','','')

    !! 3. Read variables from restart file
    IF (printlev>=3) WRITE (numout,*) 'Read DIFFUCO variables from restart file'

    CALL ioconf_setatt_p('UNITS', 's/m')
    CALL ioconf_setatt_p('LONG_NAME','Structural resistance')
    CALL restget_p (rest_id, 'rstruct', nbp_glo, nvm, 1, kjit, .TRUE., rstruct, "gather", nbp_glo, index_g)
    IF ( ALL(rstruct(:,:) == val_exp) ) THEN
       DO jv = 1, nvm
          rstruct(:,jv) = rstruct_const(jv)
       ENDDO
    ENDIF
    
    !! 4. Initialize chemistry module
    IF (printlev>=3) WRITE(numout,*) "ok_bvoc:",ok_bvoc
    IF ( ok_bvoc ) CALL chemistry_initialize(kjpindex, lalo, neighbours, resolution)
    
  END SUBROUTINE diffuco_initialize



!! ================================================================================================================================
!! SUBROUTINE    : diffuco_main
!!
!>\BRIEF	 The root subroutine for the module, which calls all other required
!! subroutines.
!! 
!! DESCRIPTION   : 

!! This is the main subroutine for the module. 
!! First it calculates the surface drag coefficient (via a call to diffuco_aero), using available parameters to determine
!! the stability of air in the surface layer by calculating the Richardson Nubmber. If a value for the 
!! surface drag coefficient is passed down from the atmospheric model and and if the flag 'ldq_cdrag_from_gcm' 
!! is set to TRUE, then the subroutine diffuco_aerod is called instead. This calculates the aerodynamic coefficient. \n
!!
!! Following this, an estimation of the saturated humidity at the surface is made (via a call
!! to qsatcalc in the module qsat_moisture). Following this the beta coefficients for sublimation (via 
!! diffuco_snow), interception (diffuco_inter), bare soil (diffuco_bare), and transpiration (via 
!! diffuco_trans_co2 if co2 is considered, diffuco_trans otherwise) are calculated in sequence. Finally 
!! the alpha and beta coefficients are combined (diffuco_comb). \n
!!
!! The surface drag coefficient is calculated for use within the module enerbil. It is required to to
!! calculate the aerodynamic coefficient for every flux. \n
!!
!! The various beta coefficients are used within the module enerbil for modifying the rate of evaporation, 
!! as appropriate for the surface. As explained in Chapter 2 of Guimberteau (2010), that module (enerbil) 
!! calculates the rate of evaporation essentially according to the expression $E = /beta E_{pot}$, where
!! E is the total evaporation and $E_{pot}$ the evaporation potential. If $\beta = 1$, there would be
!! essentially no resistance to evaporation, whereas should $\beta = 0$, there would be no evaporation and
!! the surface layer would be subject to some very stong hydrological stress. \n
!!
!! The following processes are calculated:
!! - call diffuco_aero for aerodynamic transfer coeficient
!! - call diffuco_snow for partial beta coefficient: sublimation
!! - call diffuco_inter for partial beta coefficient: interception for each type of vegetation
!! - call diffuco_bare for partial beta coefficient: bare soil
!! - call diffuco_trans for partial beta coefficient: transpiration for each type of vegetation, using Jarvis formula
!! - call diffuco_trans_co2 for partial beta coefficient: transpiration for each type of vegetation, using Farqhar's formula
!! - call diffuco_comb for alpha and beta coefficient
!! - call chemistry_bvoc for alpha and beta coefficients for biogenic emissions
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): humrel, q_cdrag, vbeta, vbeta1, vbeta4,
!! vbeta2, vbeta3, rveget, cimean   
!!
!! REFERENCE(S) :				        
!! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
!! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
!! Circulation Model. Journal of Climate, 6, pp.248-273.
!! - de Rosnay, P, 1999. Représentation des interactions sol-plante-atmosphère dans le modèle de circulation générale
!! du LMD, 1999. PhD Thesis, Université Paris 6, available (25/01/12): 
!! http://www.ecmwf.int/staff/patricia_de_rosnay/publications.html#8
!! - Ducharne, A, 1997. Le cycle de l'eau: modélisation de l'hydrologie continentale, étude de ses interactions avec 
!! le climat, PhD Thesis, Université Paris 6
!! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
!! sur le cycle de l'eau, PhD Thesis, available (25/01/12):
!! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
!! - Lathière, J, 2005. Evolution des émissions de composés organiques et azotés par la biosphère continentale dans le 
!! modèle LMDz-INCA-ORCHIDEE, Université Paris 6
!!
!! FLOWCHART	:
!! \latexonly 
!!     \includegraphics[scale=0.5]{diffuco_main_flowchart.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE diffuco_main (kjit, kjpindex, index, indexveg, indexlai, indexlai0, u, v, &
     & zlev, z0m, z0h, roughheight, roughheight_pft, temp_sol, temp_sol_pft,  temp_air, temp_growth, rau, q_cdrag, q_cdrag_pft, &
     & qsurf, qair, q2m, t2m, pb, &
     & rsol, evap_bare_lim, evapot, evapot_corr, snow, flood_frac, flood_res, frac_nobio, snow_nobio, totfrac_nobio, &
!     & rsol, evap_bare_lim, evap_bare_lim_pft, evapot, evapot_corr, snow, flood_frac, flood_res, frac_nobio, snow_nobio, totfrac_nobio, &
     & swnet, swdown, coszang, ccanopy, humrel, veget, veget_max, lai, qsintveg, qsintmax, assim_param, &
     & vbeta , vbeta_pft, vbeta1, vbeta2, vbeta3, vbeta3pot, vbeta4, vbeta4_pft, vbeta5, gsmean, rveget, rstruct, cimean, gpp, &
     & lalo, neighbours, resolution, ptnlev1, precip_rain, frac_age, tot_bare_soil, frac_snow_veg, frac_snow_nobio, &
     & hist_id, hist2_id,wtp_year)

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                         :: kjit             !! Time step number (-) 
    INTEGER(i_std), INTENT(in)                         :: kjpindex         !! Domain size (-)
    INTEGER(i_std),INTENT (in)                         :: hist_id          !! _History_ file identifier (-)
    INTEGER(i_std),INTENT (in)                         :: hist2_id         !! _History_ file 2 identifier (-)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)     :: index          !! Indeces of the points on the map (-)
    INTEGER(i_std),DIMENSION (kjpindex*(nlai+1)), INTENT (in) :: indexlai  !! Indeces of the points on the 3D map
    INTEGER(i_std),DIMENSION (kjpindex*(nlai)), INTENT (in) :: indexlai0  !! Indeces of the points on the 3D map
    INTEGER(i_std),DIMENSION (kjpindex*nvm), INTENT (in) :: indexveg       !! Indeces of the points on the 3D map (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: u                !! Eastward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: v                !! Northward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: zlev             !! Height of first layer (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: z0m              !! Surface roughness Length for momentum (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: z0h              !! Surface roughness Length for heat (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: roughheight      !! Effective height for roughness (m)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: roughheight_pft  !! Effective height for roughness (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: temp_sol         !! Skin temperature (K)
    REAL(r_std),DIMENSION (kjpindex, nvm), INTENT (in) :: temp_sol_pft         !! Skin temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: temp_air         !! Lowest level temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: temp_growth      !! Growth temperature (°C) - Is equal to t2m_month
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: rau              !! Air Density (kg m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: qsurf            !! Near surface air specific humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: qair             !! Lowest level air specific humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: q2m              !! 2m air specific humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: t2m              !! 2m air temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: snow             !! Snow mass (kg)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: flood_frac       !! Fraction of floodplains
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: flood_res        !! Reservoir in floodplains (estimation to avoid over-evaporation)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: pb               !! Surface level pressure (hPa)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: rsol             !! Bare soil evaporation resistance (s m^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: evap_bare_lim    !! Limit to the bare soil evaporation when the 
!    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: evap_bare_lim_pft    !! Limit to the bare soil evaporation when the 
                                                                           !! 11-layer hydrology is used (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: evapot           !! Soil Potential Evaporation (mm day^{-1}) 
                                                                           !! NdN This variable does not seem to be used at 
                                                                           !! all in diffuco
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: evapot_corr      !! Soil Potential Evaporation
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in) :: frac_nobio     !! Fraction of ice,lakes,cities,... (-)
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in) :: snow_nobio     !! Snow on ice,lakes,cities,... (kg m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: totfrac_nobio    !! Total fraction of ice+lakes+cities+... (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: swnet            !! Net surface short-wave flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: swdown           !! Down-welling surface short-wave flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: coszang          !! Cosine of the solar zenith angle (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: ccanopy          !! CO2 concentration inside the canopy (ppm)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget            !! Fraction of vegetation type (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget_max        !! Max. fraction of vegetation type (LAI->infty)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: lai              !! Leaf area index (m^2 m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: qsintveg         !! Water on vegetation due to interception (kg m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: qsintmax         !! Maximum water on vegetation for interception 
                                                                           !! (kg m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm,npco2), INTENT (in) :: assim_param !! min+max+opt temps, vcmax, vjmax
                                                                           !! for photosynthesis (K ??)
    REAL(r_std),DIMENSION (kjpindex,2),   INTENT (in)  :: lalo               !! Geographical coordinates
    INTEGER(i_std),DIMENSION (kjpindex,NbNeighb),INTENT (in):: neighbours    !! Vector of neighbours for each 
                                                                             !! grid point (1=N, 2=E, 3=S, 4=W)
    REAL(r_std),DIMENSION (kjpindex,2), INTENT(in)     :: resolution         !! The size in km of each grid-box in X and Y
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: ptnlev1            !! 1st level of soil temperature (Kelvin)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: precip_rain        !! Rain precipitation expressed in mm/tstep
    REAL(r_std),DIMENSION (kjpindex,nvm,nleafages), INTENT (in)  :: frac_age !! Age efficiency for isoprene emissions (from STOMATE)
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: tot_bare_soil      !! Total evaporating bare soil fraction 
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: frac_snow_veg      !! Snow cover fraction on vegeted area
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in):: frac_snow_nobio    !! Snow cover fraction on non-vegeted area

!!qcj++ peatland
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: wtp_year
    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: vbeta            !! Total beta coefficient (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: vbeta_pft            !! Total beta coefficient (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: vbeta1           !! Beta for sublimation (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: vbeta4           !! Beta for bare soil evaporation (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)  :: vbeta4_pft       !! Beta for bare soil evaporation (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: vbeta5           !! Beta for floodplains
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: gsmean           !! Mean stomatal conductance to CO2 (mol m-2 s-1) 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: vbeta2           !! Beta for interception loss (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: vbeta3           !! Beta for transpiration (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: vbeta3pot        !! Beta for potential transpiration
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: rveget           !! Stomatal resistance for the whole canopy (s m^{-1})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: rstruct          !! Structural resistance for the vegetation
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: cimean           !! Mean leaf Ci (ppm)  
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(out)  :: gpp              !! Assimilation ((gC m^{-2} s^{-1}), total area)  

    !! 0.3 Modified variables
 
    REAL(r_std),DIMENSION (kjpindex, nvm), INTENT (inout) :: humrel        !! Soil moisture stress (within range 0 to 1)
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)      :: q_cdrag       !! Surface drag coefficient  (-)
    REAL(r_std),DIMENSION (kjpindex, nvm), INTENT (inout) :: q_cdrag_pft   !!

    !! 0.4 Local variables

    REAL(r_std),DIMENSION (kjpindex,nvm)     :: vbeta23            !! Beta for fraction of wetted foliage that will
                                                                   !! transpire once intercepted water has evaporated (-)
    REAL(r_std),DIMENSION (kjpindex)         :: raero              !! Aerodynamic resistance (s m^{-1})
    INTEGER(i_std)                           :: ilaia, jv
    CHARACTER(LEN=4)                         :: laistring
    CHARACTER(LEN=80)                        :: var_name           !! To store variables names for I/O
    REAL(r_std),DIMENSION(kjpindex)          :: qsatt              !! Surface saturated humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex,nvm)     :: cim                !! Intercellular CO2 over nlai 

!_ ================================================================================================================================
    wind(:) = SQRT (u(:)*u(:) + v(:)*v(:))
  
  !! 1. Calculate the different coefficients

    IF (.NOT.ldq_cdrag_from_gcm) THEN
        ! Case 1a)
       CALL diffuco_aero (kjpindex, kjit, u, v, zlev, z0h, z0m, roughheight, roughheight_pft, temp_sol, temp_sol_pft, temp_air, &
                          qsurf, qair, snow, q_cdrag, q_cdrag_pft)
    ELSE
        !!! in coupled case, we are not able to distinguish cdrag for each PFT
        DO jv = 1,nvm
            q_cdrag_pft(:,jv) = q_cdrag
        ENDDO
    ENDIF

    ! Case 1b)
    CALL diffuco_raerod (kjpindex, u, v, q_cdrag, raero)

  !! 2. Make an estimation of the saturated humidity at the surface

    CALL qsatcalc (kjpindex, temp_sol, pb, qsatt)

  !! 3. Calculate the beta coefficient for sublimation
  
    CALL diffuco_snow (kjpindex, qair, qsatt, rau, u, v, q_cdrag, &
         snow, frac_nobio, totfrac_nobio, snow_nobio, frac_snow_veg, frac_snow_nobio, &
         vbeta1)


    CALL diffuco_flood (kjpindex, qair, qsatt, rau, u, v, q_cdrag, evapot, evapot_corr, &
         & flood_frac, flood_res, vbeta5)

  !! 4. Calculate the beta coefficient for interception

    CALL diffuco_inter (kjpindex, qair, qsatt, rau, u, v, q_cdrag, q_cdrag_pft, humrel, veget, &
       & qsintveg, qsintmax, rstruct, vbeta2, vbeta23) 


  !! 5. Calculate the beta coefficient for transpiration

    IF ( ok_co2 ) THEN

      ! case 5a)
      CALL diffuco_trans_co2 (kjpindex, lalo, swdown, pb, qsurf, q2m, t2m, temp_growth, rau, u, v, q_cdrag, q_cdrag_pft, humrel, &
                              assim_param, ccanopy, &
                              veget, veget_max, lai, qsintveg, qsintmax, vbeta3, vbeta3pot, &
                              rveget, rstruct, cimean, gsmean, gpp, vbeta23, hist_id, indexveg, indexlai, indexlai0, &
                              index, kjit, cim, wtp_year) !!!qcj++ peatland         
    ELSE

      ! case 5b) 
      CALL diffuco_trans (kjpindex, swnet, temp_air, pb, qair, rau, u, v, q_cdrag, humrel, &
           & veget, veget_max, lai, qsintveg, qsintmax, vbeta3, vbeta3pot, rveget, rstruct, cimean, &
           & gsmean, vbeta23)

    ENDIF


    !
    !biogenic emissions
    !
    IF ( ok_bvoc ) THEN
       CALL chemistry_bvoc (kjpindex, swdown, coszang, temp_air, &
            temp_sol, ptnlev1, precip_rain, humrel, veget_max, &
            lai, frac_age, lalo, ccanopy, cim, wind, snow, &
            veget, hist_id, hist2_id, kjit, index, &
            indexlai, indexveg)
    ENDIF
    !
    ! combination of coefficient : alpha and beta coefficient
    ! beta coefficient for bare soil
    !

!    CALL diffuco_bare (kjpindex, u, v, q_cdrag, rsol, evap_bare_lim, evap_bare_lim_pft, humrel, &
    CALL diffuco_bare (kjpindex, u, v, q_cdrag, rsol, evap_bare_lim, humrel, &
         veget, veget_max, tot_bare_soil, vbeta2, vbeta3, vbeta4, vbeta4_pft)

  !! 6. Combine the alpha and beta coefficients

    ! Ajout qsintmax dans les arguments de la routine.... Nathalie / le 13-03-2006
    CALL diffuco_comb (kjpindex, humrel, rau, u, v, q_cdrag, pb, qair, temp_sol, temp_air, snow, &
       & veget, veget_max, lai, tot_bare_soil, vbeta1, vbeta2, vbeta3, vbeta4, vbeta4_pft, vbeta, vbeta_pft, qsintmax)

    CALL xios_orchidee_send_field("q_cdrag",q_cdrag)
    CALL xios_orchidee_send_field("cdrag_pft",q_cdrag_pft)
    CALL xios_orchidee_send_field("raero",raero)
    CALL xios_orchidee_send_field("wind",wind)
    CALL xios_orchidee_send_field("qsatt",qsatt)
    CALL xios_orchidee_send_field("coszang",coszang)
    IF ( ok_co2 ) CALL xios_orchidee_send_field('cim', cim)

    IF ( .NOT. almaoutput ) THEN
       CALL histwrite_p(hist_id, 'raero', kjit, raero, kjpindex, index)
       CALL histwrite_p(hist_id, 'cdrag', kjit, q_cdrag, kjpindex, index)
       CALL histwrite_p(hist_id, 'cdrag_pft', kjit, q_cdrag_pft, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'Wind', kjit, wind, kjpindex, index)
       CALL histwrite_p(hist_id, 'qsatt', kjit, qsatt, kjpindex, index)
       IF ( ok_co2 ) CALL histwrite_p(hist_id, 'cim', kjit, cim, kjpindex*nvm, indexveg)

       IF ( hist2_id > 0 ) THEN
          CALL histwrite_p(hist2_id, 'raero', kjit, raero, kjpindex, index)
          CALL histwrite_p(hist2_id, 'cdrag', kjit, q_cdrag, kjpindex, index)
!          CALL histwrite_p(hist2_id, 'cdrag_pft', kjit, q_cdrag_pft, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'Wind', kjit, wind, kjpindex, index)
          CALL histwrite_p(hist2_id, 'qsatt', kjit, qsatt, kjpindex, index)
       ENDIF
    ELSE
       IF ( ok_co2 ) CALL histwrite_p(hist_id, 'cim', kjit, cim, kjpindex*nvm, indexveg)
    ENDIF

    IF (printlev>=3) WRITE (numout,*) ' diffuco_main done '

  END SUBROUTINE diffuco_main

!!  =============================================================================================================================
!! SUBROUTINE: diffuco_finalize
!!
!>\BRIEF          Write to restart file
!!
!! DESCRIPTION:   This subroutine writes the module variables and variables calculated in diffuco
!!                to restart file
!!
!! RECENT CHANGE(S): None
!! REFERENCE(S): None
!! FLOWCHART: None
!! \n
!_ ==============================================================================================================================
  SUBROUTINE diffuco_finalize (kjit, kjpindex, rest_id, rstruct, q_cdrag_pft )

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                         :: kjit             !! Time step number (-) 
    INTEGER(i_std), INTENT(in)                         :: kjpindex         !! Domain size (-)
    INTEGER(i_std),INTENT (in)                         :: rest_id          !! _Restart_ file identifier (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: rstruct          !! Structural resistance for the vegetation
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: q_cdrag_pft      !! drag coefficient for the vegetation

    !! 0.4 Local variables
    INTEGER                                            :: ilai
    CHARACTER(LEN=4)                                   :: laistring
    CHARACTER(LEN=80)                                  :: var_name        

!_ ================================================================================================================================
    
  !! 1. Prepare the restart file for the next simulation
    IF (printlev>=3) WRITE (numout,*) 'Complete restart file with DIFFUCO variables '

    CALL restput_p (rest_id, 'rstruct', nbp_glo, nvm, 1, kjit, rstruct, 'scatter',  nbp_glo, index_g)

    CALL restput_p (rest_id, 'cdrag_pft', nbp_glo, nvm, 1, kjit, q_cdrag_pft, 'scatter', nbp_glo, index_g)
    
  END SUBROUTINE diffuco_finalize


!! ================================================================================================================================
!! SUBROUTINE		 			: diffuco_clear
!!
!>\BRIEF					Housekeeping module to deallocate the variables
!! rstruct and raero
!!
!! DESCRIPTION				        : Housekeeping module to deallocate the variables
!! rstruct and raero
!!
!! RECENT CHANGE(S)                             : None
!!
!! MAIN OUTPUT VARIABLE(S)	                : None
!!
!! REFERENCE(S)				        : None
!!
!! FLOWCHART                                    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE diffuco_clear()

    ! Deallocate and reset variables in chemistry module
    CALL chemistry_clear

  END SUBROUTINE diffuco_clear


!! ================================================================================================================================
!! SUBROUTINE	: diffuco_aero
!!
!>\BRIEF	This module first calculates the surface drag 
!! coefficient, for cases in which the surface drag coefficient is NOT provided by the coupled 
!! atmospheric model LMDZ or when the flag ldq_cdrag_from_gcm is set to FALSE 
!!
!! DESCRIPTION	: Computes the surface drag coefficient, for cases 
!! in which it is NOT provided by the coupled atmospheric model LMDZ. The module first uses the 
!! meteorolgical input to calculate the Richardson Number, which is an indicator of atmospheric 
!! stability in the surface layer. The formulation used to find this surface drag coefficient is 
!! dependent on the stability determined. \n
!!
!! Designation of wind speed
!! \latexonly 
!!     \input{diffucoaero1.tex}
!! \endlatexonly
!!
!! Calculation of geopotential. This is the definition of Geopotential height (e.g. Jacobson 
!! eqn.4.47, 2005). (required for calculation of the Richardson Number)
!! \latexonly 
!!     \input{diffucoaero2.tex}
!! \endlatexonly
!! 
!! \latexonly 
!!     \input{diffucoaero3.tex}
!! \endlatexonly
!!
!! Calculation of the virtual air temperature at the surface (required for calculation
!! of the Richardson Number)
!! \latexonly 
!!     \input{diffucoaero4.tex}
!! \endlatexonly
!!
!! Calculation of the virtual surface temperature (required for calculation of th
!! Richardson Number)
!! \latexonly 
!!     \input{diffucoaero5.tex}
!! \endlatexonly
!!
!! Calculation of the squared wind shear (required for calculation of the Richardson
!! Number)
!! \latexonly 
!!     \input{diffucoaero6.tex}
!! \endlatexonly
!! 
!! Calculation of the Richardson Number. The Richardson Number is defined as the ratio 
!! of potential to kinetic energy, or, in the context of atmospheric science, of the
!! generation of energy by wind shear against consumption
!! by static stability and is an indicator of flow stability (i.e. for when laminar flow 
!! becomes turbulent and vise versa). It is approximated using the expression below:
!! \latexonly 
!!     \input{diffucoaero7.tex}
!! \endlatexonly
!!
!! The Richardson Number hence calculated is subject to a minimum value:
!! \latexonly 
!!     \input{diffucoaero8.tex}
!! \endlatexonly
!! 
!! Computing the drag coefficient. We add the add the height of the vegetation to the 
!! level height to take into account that the level 'seen' by the vegetation is actually 
!! the top of the vegetation. Then we we can subtract the displacement height.
!! \latexonly 
!!     \input{diffucoaero9.tex}
!! \endlatexonly
!! 
!! For the stable case (i.e $R_i$ $\geq$ 0)
!! \latexonly 
!!     \input{diffucoaero10.tex}
!! \endlatexonly
!!
!! \latexonly 
!!     \input{diffucoaero11.tex}
!! \endlatexonly
!!          
!! For the unstable case (i.e. $R_i$ < 0)
!! \latexonly 
!!     \input{diffucoaero12.tex}
!! \endlatexonly
!!
!! \latexonly 
!!     \input{diffucoaero13.tex}
!! \endlatexonly
!!               
!! If the Drag Coefficient becomes too small than the surface may uncouple from the atmosphere.
!! To prevent this, a minimum limit to the drag coefficient is defined as:
!!
!! \latexonly 
!!     \input{diffucoaero14.tex}
!! \endlatexonly
!! 
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): q_cdrag
!!
!! REFERENCE(S)	: 
!! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
!! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
!! Circulation Model. Journal of Climate, 6, pp.248-273
!! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
!! sur le cycle de l'eau, PhD Thesis, available from:
!! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
!! - Jacobson M.Z., Fundamentals of Atmospheric Modeling (2nd Edition), published Cambridge 
!! University Press, ISBN 0-521-54865-9
!!
!! FLOWCHART	:
!! \latexonly 
!!     \includegraphics[scale=0.5]{diffuco_aero_flowchart.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE diffuco_aero (kjpindex, kjit, u, v, zlev, z0h, z0m, roughheight, roughheight_pft, temp_sol, temp_sol_pft, temp_air, &
                           qsurf, qair, snow, q_cdrag, q_cdrag_pft)

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                          :: kjpindex, kjit   !! Domain size
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: u                !! Eastward Lowest level wind speed (m s^{-1}) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: v                !! Northward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: zlev             !! Height of first atmospheric layer (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: z0h               !! Surface roughness Length for heat (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: z0m               !! Surface roughness Length for momentum (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: roughheight      !! Effective roughness height (m)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: roughheight_pft  
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: temp_sol         !! Ground temperature (K)
    REAL(r_std),DIMENSION (kjpindex, nvm), INTENT (in)  :: temp_sol_pft     !! Ground temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: temp_air         !! Lowest level temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: qsurf            !! near surface specific air humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: qair             !! Lowest level specific air humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: snow             !! Snow mass (kg)

    !! 0.2 Output variables
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: q_cdrag          !! Surface drag coefficient  (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: q_cdrag_pft     

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                      :: ji, jv
    REAL(r_std)                                         :: speed, zg, zdphi, ztvd, ztvs, zdu2
    REAL(r_std)                                         :: zri, cd_neut, zscf, cd_tmp
    REAL(r_std),DIMENSION(nvm)                          :: cd_neut_pft, cd_tmp_pft, ztvs_pft, zri_pft, zscf_pft
    REAL(r_std)                                         :: snowfact
!_ ================================================================================================================================
    snowfact = 1
    zri_pft = zero
    cd_neut_pft = zero
    zscf_pft = zero

  !! 1. Initialisation

    ! test if we have to work with q_cdrag or to calcul it
    DO ji=1,kjpindex
       
       !! 1a).1 Designation of wind speed
       !! \latexonly 
       !!     \input{diffucoaero1.tex}
       !! \endlatexonly
       speed = wind(ji)
    
       !! 1a).2 Calculation of geopotentiel
       !! This is the definition of Geopotential height (e.g. Jacobson eqn.4.47, 2005). (required
       !! for calculation of the Richardson Number)
       !! \latexonly 
       !!     \input{diffucoaero2.tex}
       !! \endlatexonly
       zg = zlev(ji) * cte_grav
      
       !! \latexonly 
       !!     \input{diffucoaero3.tex}
       !! \endlatexonly
       zdphi = zg/cp_air
       
       !! 1a).3 Calculation of the virtual air temperature at the surface 
       !! required for calculation of the Richardson Number
       !! \latexonly 
       !!     \input{diffucoaero4.tex}
       !! \endlatexonly
       ztvd = (temp_air(ji) + zdphi / (un + rvtmp2 * qair(ji))) * (un + retv * qair(ji)) 
       
       !! 1a).4 Calculation of the virtual surface temperature 
       !! required for calculation of the Richardson Number
       !! \latexonly 
       !!     \input{diffucoaero5.tex}
       !! \endlatexonly
       ztvs = temp_sol(ji) * (un + retv * qsurf(ji))
       DO jv = 1,nvm
           IF (ok_LAIdev(jv)) THEN
               ztvs_pft(jv) = temp_sol_pft(ji,jv) * (un + retv * qsurf(ji))
           ENDIF
       ENDDO
     
       !! 1a).5 Calculation of the squared wind shear 
       !! required for calculation of the Richardson Number
       !! \latexonly 
       !!     \input{diffucoaero6.tex}
       !! \endlatexonly
       zdu2 = MAX(cepdu2,speed**2)
       
       !! 1a).6 Calculation of the Richardson Number
       !!  The Richardson Number is defined as the ratio of potential to kinetic energy, or, in the 
       !!  context of atmospheric science, of the generation of energy by wind shear against consumption
       !!  by static stability and is an indicator of flow stability (i.e. for when laminar flow 
       !!  becomes turbulent and vise versa).\n
       !!  It is approximated using the expression below:
       !!  \latexonly 
       !!     \input{diffucoaero7.tex}
       !! \endlatexonly
       zri = zg * (ztvd - ztvs) / (zdu2 * ztvd)
       DO jv = 2,nvm
           IF (ok_LAIdev(jv)) THEN
               zri_pft(jv) = zg * (ztvd - ztvs_pft(jv)) / (zdu2 * ztvd)
               zri_pft(jv) = MAX(MIN(zri_pft(jv),5.),-5.)
           ENDIF
       ENDDO
      
       !! The Richardson Number hence calculated is subject to a minimum value:
       !! \latexonly 
       !!     \input{diffucoaero8.tex}
       !! \endlatexonly       
       zri = MAX(MIN(zri,5.),-5.)
       
       !! 1a).7 Computing the drag coefficient
       !!  We add the add the height of the vegetation to the level height to take into account
       !!  that the level 'seen' by the vegetation is actually the top of the vegetation. Then we 
       !!  we can subtract the displacement height.
       !! \latexonly 
       !!     \input{diffucoaero9.tex}
       !! \endlatexonly

       !! 7.0 Snow smoothering
       !! Snow induces low levels of turbulence.
       !! Sensible heat fluxes can therefore be reduced of ~1/3. Pomeroy et al., 1998
       cd_neut = ct_karman ** 2. / ( LOG( (zlev(ji) + roughheight(ji)) / z0m(ji) ) * LOG( (zlev(ji) + roughheight(ji)) / z0h(ji) ) )
       DO jv = 2,nvm
           cd_neut_pft(jv) = ct_karman ** 2. / ( LOG( (zlev(ji) + roughheight_pft(ji, jv)) / z0m(ji) ) * LOG( (zlev(ji) + roughheight_pft(ji, jv)) / z0h(ji) ) )
       ENDDO
       
       !! 1a).7.1 - for the stable case (i.e $R_i$ $\geq$ 0)
       IF (zri .GE. zero) THEN
          
          !! \latexonly 
          !!     \input{diffucoaero10.tex}
          !! \endlatexonly
          zscf = SQRT(un + cd * ABS(zri))
         
          !! \latexonly 
          !!     \input{diffucoaero11.tex}
          !! \endlatexonly          
          cd_tmp=cd_neut/(un + trois * cb * zri * zscf)
!          DO jv = 2,nvm
!              cd_tmp_pft(jv) = cd_neut_pft(jv)/(un + trois * cb * zri * zscf)
!          ENDDO
       ELSE
          
          !! 1a).7.2 - for the unstable case (i.e. $R_i$ < 0)
          !! \latexonly 
          !!     \input{diffucoaero12.tex}
          !! \endlatexonly
          zscf = un / (un + trois * cb * cc * cd_neut * SQRT(ABS(zri) * &
               & ((zlev(ji) + roughheight(ji)) / z0m(ji))))

          !! \latexonly 
          !!     \input{diffucoaero13.tex}
          !! \endlatexonly               
          cd_tmp=cd_neut * (un - trois * cb * zri * zscf)
!          DO jv = 2,nvm
!              zscftemp = un / (un + trois * cb * cc * cd_neut_pft(jv) * SQRT(ABS(zri) * &
!                       & ((zlev(ji) + roughheight_pft(ji,jv)) / z0(ji))))
!              cd_tmp_pft(jv) = cd_neut_pft(jv) * (un - trois * cb * zri * zscftemp)
!          ENDDO

       ENDIF

       DO jv = 2,nvm 
           IF (ok_LAIdev(jv)) THEN
               IF (zri_pft(jv) .GE. zero) THEN ! stable case
                   zscf_pft(jv) = SQRT(un + cd*ABS(zri_pft(jv)))
                   cd_tmp_pft(jv) = cd_neut_pft(jv)/(un + trois * cb * zri_pft(jv) * zscf_pft(jv))
               ELSE ! unstable case
                   zscf_pft(jv) = un / (un + trois * cb * cc * cd_neut_pft(jv) * SQRT(ABS(zri_pft(jv)) * &
                                & ((zlev(ji) + roughheight_pft(ji,jv)) / z0m(ji)/snowfact)))
                   cd_tmp_pft(jv) = cd_neut_pft(jv) * (un - trois * cb * zri_pft(jv) * zscf_pft(jv))
               ENDIF
           ENDIF
       ENDDO
       
       !! If the Drag Coefficient becomes too small than the surface may uncouple from the atmosphere.
       !! To prevent this, a minimum limit to the drag coefficient is defined as:
       
       !! \latexonly 
       !!     \input{diffucoaero14.tex}
       !! \endlatexonly
       !!
       q_cdrag(ji) = MAX(cd_tmp, 1.e-4/MAX(speed,min_wind))
       DO jv = 2,nvm
           IF (ok_LAIdev(jv)) THEN
               q_cdrag_pft(ji,jv) = MAX(cd_tmp_pft(jv), 1.e-4/MAX(speed,min_wind))
           ELSE
               q_cdrag_pft(ji,jv) = q_cdrag(ji)
           ENDIF
       ENDDO

       ! In some situations it might be useful to give an upper limit on the cdrag as well. 
       ! The line here should then be uncommented.
      !q_cdrag(ji) = MIN(q_cdrag(ji), 0.5/MAX(speed,min_wind))

    END DO

    IF (printlev>=3) WRITE (numout,*) ' not ldqcdrag_from_gcm : diffuco_aero done '

  END SUBROUTINE diffuco_aero


!! ================================================================================================================================
!! SUBROUTINE    : diffuco_snow
!!
!>\BRIEF         This subroutine computes the beta coefficient for snow sublimation.
!!
!! DESCRIPTION   : This routine computes beta coefficient for snow sublimation, which
!! integrates the snow on both vegetation and other surface types (e.g. ice, lakes,
!! cities etc.) \n
!!
!! A critical depth of snow (snowcri) is defined to calculate the fraction of each grid-cell
!! that is covered with snow (snow/snowcri) while the remaining part is snow-free.
!! We also carry out a first calculation of sublimation (subtest) to lower down the beta
!! coefficient if necessary (if subtest > snow). This is a predictor-corrector test. 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::vbeta1 
!!
!! REFERENCE(S) :
!! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
!! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
!! Circulation Model. Journal of Climate, 6, pp. 248-273
!! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
!! sur le cycle de l'eau, PhD Thesis, available from:
!! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  
SUBROUTINE diffuco_snow (kjpindex, qair, qsatt, rau, u, v,q_cdrag, &
       & snow, frac_nobio, totfrac_nobio, snow_nobio, frac_snow_veg, frac_snow_nobio, &
       vbeta1)

  !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables
 
    INTEGER(i_std), INTENT(in)                           :: kjpindex       !! Domain size (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: qair           !! Lowest level specific air humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: qsatt          !! Surface saturated humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: rau            !! Air density (kg m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: u              !! Eastward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: v              !! Northward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: q_cdrag        !! Surface drag coefficient  (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: snow           !! Snow mass (kg m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in) :: frac_nobio     !! Fraction of ice, lakes, cities etc. (-)
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in) :: snow_nobio     !! Snow on ice, lakes, cities etc. (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: totfrac_nobio  !! Total fraction of ice, lakes, cities etc. (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)         :: frac_snow_veg  !! Snow cover fraction on vegeted area
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in)  :: frac_snow_nobio!! Snow cover fraction on non-vegeted area
    
    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)       :: vbeta1         !! Beta for sublimation (dimensionless ratio) 
    
    !! 0.3 Modified variables

    !! 0.4 Local variables

    REAL(r_std)                                          :: subtest        !! Sublimation for test (kg m^{-2})
    REAL(r_std)                                          :: zrapp          !! Modified factor (ratio)
    REAL(r_std)                                          :: speed          !! Wind speed (m s^{-1})
    REAL(r_std)                                          :: vbeta1_add     !! Beta for sublimation (ratio)
    INTEGER(i_std)                                       :: ji, jv         !! Indices (-)
!_ ================================================================================================================================

  !! 1. Calculate beta coefficient for snow sublimation on the vegetation\n

    DO ji=1,kjpindex  ! Loop over # pixels - domain size

       ! Fraction of mesh that can sublimate snow
       vbeta1(ji) = (un - totfrac_nobio(ji)) * frac_snow_veg(ji)

       ! Limitation of sublimation in case of snow amounts smaller than the atmospheric demand. 
       speed = MAX(min_wind, wind(ji))

       subtest = dt_sechiba * vbeta1(ji) * speed * q_cdrag(ji) * rau(ji) * &
               & ( qsatt(ji) - qair(ji) )

       IF ( subtest .GT. min_sechiba ) THEN
          zrapp = snow(ji) / subtest
          IF ( zrapp .LT. un ) THEN
             vbeta1(ji) = vbeta1(ji) * zrapp
          ENDIF
       ENDIF

    END DO ! Loop over # pixels - domain size

  !! 2. Add the beta coefficients calculated from other surfaces types (snow on ice,lakes, cities...)

    DO jv = 1, nnobio ! Loop over # other surface types
!!$      !
!!$      IF ( jv .EQ. iice ) THEN
!!$        !
!!$        !  Land ice is of course a particular case
!!$        !
!!$        DO ji=1,kjpindex
!!$          vbeta1(ji) = vbeta1(ji) + frac_nobio(ji,jv)
!!$        ENDDO
!!$        !
!!$      ELSE
        !
        DO ji=1,kjpindex ! Loop over # pixels - domain size

           vbeta1_add = frac_nobio(ji,jv) * frac_snow_nobio(ji, jv)

           ! Limitation of sublimation in case of snow amounts smaller than
           ! the atmospheric demand. 
           speed = MAX(min_wind, wind(ji))
            
            !!     Limitation of sublimation by the snow accumulated on the ground 
            !!     A first approximation is obtained with the old values of
            !!     qair and qsol_sat: function of temp-sol and pb. (see call of qsatcalc)
           subtest = dt_sechiba * vbeta1_add * speed * q_cdrag(ji) * rau(ji) * &
                & ( qsatt(ji) - qair(ji) )

           IF ( subtest .GT. min_sechiba ) THEN
              zrapp = snow_nobio(ji,jv) / subtest
              IF ( zrapp .LT. un ) THEN
                 vbeta1_add = vbeta1_add * zrapp
              ENDIF
           ENDIF

           vbeta1(ji) = vbeta1(ji) + vbeta1_add

        ENDDO ! Loop over # pixels - domain size

!!$      ENDIF
      
    ENDDO ! Loop over # other surface types

    IF (printlev>=3) WRITE (numout,*) ' diffuco_snow done '

  END SUBROUTINE diffuco_snow


!! ================================================================================================================================
!! SUBROUTINE		 			: diffuco_flood 
!!
!>\BRIEF				       	This routine computes partial beta coefficient : floodplains
!!
!! DESCRIPTION				        : 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S)	                : vbeta5
!!
!! REFERENCE(S)				        : None
!!
!! FLOWCHART                                    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE diffuco_flood (kjpindex, qair, qsatt, rau, u, v, q_cdrag, evapot, evapot_corr, &
       & flood_frac, flood_res, vbeta5)

    ! interface description
    ! input scalar 
    INTEGER(i_std), INTENT(in)                               :: kjpindex   !! Domain size
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: qair       !! Lowest level specific humidity
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: qsatt      !! Surface saturated humidity
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: rau        !! Density
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: u          !! Lowest level wind speed 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: v          !! Lowest level wind speed
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: q_cdrag    !! Surface drag coefficient  (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: flood_res  !! water mass in flood reservoir
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: flood_frac !! fraction of floodplains
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: evapot     !! Potential evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: evapot_corr!! Potential evaporation2
    ! output fields
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: vbeta5     !! Beta for floodplains

    ! local declaration
    REAL(r_std)                                              :: subtest, zrapp, speed
    INTEGER(i_std)                                           :: ji, jv

!_ ================================================================================================================================
    !
    ! beta coefficient for sublimation for floodplains
    !
    DO ji=1,kjpindex
       !
       IF (evapot(ji) .GT. min_sechiba) THEN
          vbeta5(ji) = flood_frac(ji) *evapot_corr(ji)/evapot(ji)
       ELSE
          vbeta5(ji) = flood_frac(ji)
       ENDIF
       !
       ! -- Limitation of evaporation in case of water amounts smaller than
       !    the atmospheric demand. 
       
       !
       speed = MAX(min_wind, SQRT (u(ji)*u(ji) + v(ji)*v(ji)))
       !
       subtest = dt_sechiba * vbeta5(ji) * speed * q_cdrag(ji) * rau(ji) * &
               & ( qsatt(ji) - qair(ji) )
       !  
       IF ( subtest .GT. min_sechiba ) THEN
          zrapp = flood_res(ji) / subtest
          IF ( zrapp .LT. un ) THEN
             vbeta5(ji) = vbeta5(ji) * zrapp
          ENDIF
       ENDIF
       !
    END DO

    IF (printlev>=3) WRITE (numout,*) ' diffuco_flood done '

  END SUBROUTINE diffuco_flood


!! ================================================================================================================================
!! SUBROUTINE    : diffuco_inter
!!
!>\BRIEF	 This routine computes the partial beta coefficient
!! for the interception for each type of vegetation
!!
!! DESCRIPTION   : We first calculate the dry and wet parts of each PFT (wet part = qsintveg/qsintmax).
!! The former is submitted to transpiration only (vbeta3 coefficient, calculated in 
!! diffuco_trans or diffuco_trans_co2), while the latter is first submitted to interception loss 
!! (vbeta2 coefficient) and then to transpiration once all the intercepted water has been evaporated 
!! (vbeta23 coefficient). Interception loss is also submitted to a predictor-corrector test, 
!! as for snow sublimation. \n
!!
!! \latexonly 
!!     \input{diffucointer1.tex}
!! \endlatexonly
!! Calculate the wet fraction of vegetation as  the ration between the intercepted water and the maximum water 
!! on the vegetation. This ratio defines the wet portion of vegetation that will be submitted to interception loss.
!!
!! \latexonly 
!!     \input{diffucointer2.tex}
!! \endlatexonly
!!
!! Calculation of $\beta_3$, the canopy transpiration resistance
!! \latexonly 
!!     \input{diffucointer3.tex}
!! \endlatexonly            
!! 
!! We here determine the limitation of interception loss by the water stored on the leaf. 
!! A first approximation of interception loss is obtained using the old values of
!! qair and qsol_sat, which are functions of temp-sol and pb. (see call of 'qsatcalc')
!! \latexonly 
!!     \input{diffucointer4.tex}
!! \endlatexonly
!!
!! \latexonly
!!     \input{diffucointer5.tex}
!! \endlatexonly
!!
!! \latexonly 
!!     \input{diffucointer6.tex}
!! \endlatexonly
!!
!! Once the whole water stored on foliage has evaporated, transpiration can take place on the fraction
!! 'zqsvegrap'.
!! \latexonly 
!!     \input{diffucointer7.tex}
!! \endlatexonly
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::vbeta2, ::vbeta23
!!
!! REFERENCE(S) :
!! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
!! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
!! Circulation Model. Journal of Climate, 6, pp. 248-273
!! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
!! sur le cycle de l'eau, PhD Thesis, available from:
!! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
!! - Perrier, A, 1975. Etude physique de l'évaporation dans les conditions naturelles. Annales 
!! Agronomiques, 26(1-18): pp. 105-123, pp. 229-243
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE diffuco_inter (kjpindex, qair, qsatt, rau, u, v, q_cdrag, q_cdrag_pft, humrel, veget, &
     & qsintveg, qsintmax, rstruct, vbeta2, vbeta23)
   
  !! 0 Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                           :: kjpindex   !! Domain size (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: qair       !! Lowest level specific air humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: qsatt      !! Surface saturated humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: rau        !! Air Density (kg m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: u          !! Eastward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: v          !! Northward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: q_cdrag    !! Surface drag coefficient  (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: q_cdrag_pft   !!Product of Surface drag coefficient and wind
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: humrel     !! Soil moisture stress (within range 0 to 1)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: veget      !! vegetation fraction for each type (fraction)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: qsintveg   !! Water on vegetation due to interception (kg m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: qsintmax   !! Maximum water on vegetation (kg m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: rstruct    !! architectural resistance (s m^{-1})
    
    !! 0.2 Output variables
    
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)   :: vbeta2     !! Beta for interception loss (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)   :: vbeta23    !! Beta for fraction of wetted foliage that will 
                                                                       !! transpire (-)

    !! 0.4 Local variables

    INTEGER(i_std)                                       :: ji, jv                               !! (-), (-)
    REAL(r_std)                                          :: zqsvegrap, ziltest, zrapp, speed     !!
!_ ================================================================================================================================

  !! 1. Initialize

    vbeta2(:,:) = zero
    vbeta23(:,:) = zero
   
  !! 2. The beta coefficient for interception by vegetation. 
    
    DO jv = 2,nvm

      DO ji=1,kjpindex

         IF (veget(ji,jv) .GT. min_sechiba .AND. qsintveg(ji,jv) .GT. zero ) THEN

            zqsvegrap = zero
            IF (qsintmax(ji,jv) .GT. min_sechiba ) THEN

            !! \latexonly 
            !!     \input{diffucointer1.tex}
            !! \endlatexonly
            !!
            !! We calculate the wet fraction of vegetation as  the ration between the intercepted water and the maximum water 
            !! on the vegetation. This ratio defines the wet portion of vegetation that will be submitted to interception loss.
            !!
                zqsvegrap = MAX(zero, qsintveg(ji,jv) / qsintmax(ji,jv))
            END IF

            !! \latexonly 
            !!     \input{diffucointer2.tex}
            !! \endlatexonly
            speed = MAX(min_wind, wind(ji))

            !! Calculation of $\beta_3$, the canopy transpiration resistance
            !! \latexonly 
            !!     \input{diffucointer3.tex}
            !! \endlatexonly
            IF (.NOT. ok_LAIdev(jv)) THEN
                vbeta2(ji,jv) = veget(ji,jv) * zqsvegrap * (un / (un + speed * q_cdrag(ji) * rstruct(ji,jv)))
            ELSE
                vbeta2(ji,jv) = veget(ji,jv) * zqsvegrap * (un / (un + speed * q_cdrag_pft(ji,jv) * rstruct(ji,jv)))
                ! for crops, assuming separate field with separate drag
                ! coefficient for each crop
                ! this will affect the estimates of transpiration, xuhui
            ENDIF
            
            !! We here determine the limitation of interception loss by the water stored on the leaf. 
            !! A first approximation of interception loss is obtained using the old values of
            !! qair and qsol_sat, which are functions of temp-sol and pb. (see call of 'qsatcalc')
            !! \latexonly 
            !!     \input{diffucointer4.tex}
            !! \endlatexonly
            ziltest = dt_sechiba * vbeta2(ji,jv) * speed * q_cdrag(ji) * rau(ji) * &
               & ( qsatt(ji) - qair(ji) )

            IF ( ziltest .GT. min_sechiba ) THEN

                !! \latexonly 
                !!     \input{diffucointer5.tex}
                !! \endlatexonly
                zrapp = qsintveg(ji,jv) / ziltest
                IF ( zrapp .LT. un ) THEN
                   
                    !! \latexonly 
                    !!     \input{diffucointer6.tex}
                    !! \endlatexonly
                    !!
		    !! Once the whole water stored on foliage has evaporated, transpiration can take place on the fraction
                    !! 'zqsvegrap'.
                   IF ( humrel(ji,jv) >= min_sechiba ) THEN
                      vbeta23(ji,jv) = MAX(vbeta2(ji,jv) - vbeta2(ji,jv) * zrapp, zero)
                   ELSE
                      ! We don't want transpiration when the soil cannot deliver it
                      vbeta23(ji,jv) = zero
                   ENDIF
                    
                    !! \latexonly 
                    !!     \input{diffucointer7.tex}
                    !! \endlatexonly
                    vbeta2(ji,jv) = vbeta2(ji,jv) * zrapp
                ENDIF
            ENDIF
        END IF
!        ! Autre formulation possible pour l'evaporation permettant une transpiration sur tout le feuillage
!        !commenter si formulation Nathalie sinon Tristan
!        speed = MAX(min_wind, wind(ji))
!        
!        vbeta23(ji,jv) = MAX(zero, veget(ji,jv) * (un / (un + speed * q_cdrag(ji) * rstruct(ji,jv))) - vbeta2(ji,jv))

      END DO

    END DO

    IF (printlev>=3) WRITE (numout,*) ' diffuco_inter done '

  END SUBROUTINE diffuco_inter


!! ==============================================================================================================================
!! SUBROUTINE      : diffuco_bare
!!
!>\BRIEF	   This routine computes the partial beta coefficient corresponding to
!! bare soil
!!
!! DESCRIPTION	   : Bare soil evaporation is either limited by a soil resistance 
!! (rsol) that is calculated in hydrolc.f90, when Choisnel hydrology is used or 
!! submitted to a maximum possible flow (evap_bare_lim) if the 11-layer hydrology is used.\n
!! 
!! Calculation of wind speed
!! \latexonly 
!!     \input{diffucobare1.tex}
!! \endlatexonly
!!             
!! The calculation of $\beta_4$
!! \latexonly 
!!     \input{diffucobare2.tex}
!! \endlatexonly
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::vbeta4
!!
!! REFERENCE(S)	 :
!! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
!! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
!! Circulation Model. Journal of Climate, 6, pp.248-273
!! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
!! sur le cycle de l'eau, PhD Thesis, available from:
!! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

!  SUBROUTINE diffuco_bare (kjpindex, u, v, q_cdrag, rsol, evap_bare_lim, evap_bare_lim_pft, humrel, &
  SUBROUTINE diffuco_bare (kjpindex, u, v, q_cdrag, rsol, evap_bare_lim, humrel, &
       & veget, veget_max, tot_bare_soil, vbeta2, vbeta3, vbeta4, vbeta4_pft)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                         :: kjpindex       !! Domain size (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: u              !! Eastward Lowest level wind speed (m s^{-1}) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: v              !! Northward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: q_cdrag        !!  Surface drag coefficient  (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: rsol           !! resistance for bare soil evaporation  (s m^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: evap_bare_lim  !! limiting factor for bare soil evaporation when the 
!    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: evap_bare_lim_pft  !! limiting factor for bare soil evaporation when the 
                                                                         !! 11-layer hydrology is used (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: humrel         !! Soil moisture stress (within range 0 to 1)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget          !! Type of vegetation fraction (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget_max      !! Type of vegetation max fraction
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: vbeta2         !! Beta for Interception 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: vbeta3         !! Beta for Transpiration 
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)      :: tot_bare_soil  !! Total evaporating bare soil fraction 

    !! 0.2 Output variables
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: vbeta4         !! Beta for bare soil evaporation (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: vbeta4_pft     !! Beta for bare soil evaporation (-)
    
    !! 0.3 Modified variables
        
    !! 0.4 Local variables
    REAL(r_std)                                    :: humveg_prod

    INTEGER(i_std)                                     :: ji, jv
    REAL(r_std)                                        :: speed          !! Surface wind speed (m s^{-1})
!_ ================================================================================================================================

    vbeta4 = zero
    vbeta4_pft = zero

  !! 1. Calculation of the soil resistance and the beta (beta_4) for bare soil

    IF ( .NOT. hydrol_cwrr ) THEN
       DO ji = 1, kjpindex
          
          vbeta4(ji) = zero
          !     
          ! 1.   Soil resistance and beta for bare soil
          !      note: tot_bare_soil contains the fraction of bare soil
          !            see slowproc module
          !
          speed = MAX(min_wind, wind(ji))
          !
          humveg_prod = tot_bare_soil(ji) * humrel(ji,1)
          !
          DO jv = 2, nvm
             humveg_prod = humveg_prod + veget(ji,jv) * humrel(ji,jv)
          ENDDO
            
             !! \latexonly 
             !!     \input{diffucobare1.tex}
             !! \endlatexonly
          IF (tot_bare_soil(ji) .GE. min_sechiba) THEN
             
             ! Correction Nathalie de Noblet - le 27 Mars 2006
             ! Selon recommandation de Frederic Hourdin: supprimer humrel dans formulation vbeta4
             !vbeta4(ji) = tot_bare_soil(ji) *humrel(ji,1)* (un / (un + speed * q_cdrag(ji) * rsol(ji)))
             ! Nathalie - le 28 mars 2006 - vbeta4 n'etait pas calcule en fonction de
             ! rsol mais de rsol_cste * hdry! Dans ce cas inutile de calculer rsol(ji)!!
             vbeta4(ji) = tot_bare_soil(ji) * (un / (un + speed * q_cdrag(ji) * rsol(ji)))
             
          ENDIF
          !Commenter la ligne ci-dessous si calcul Nathalie sinon Tristan
!          vbeta4(ji) = MIN(humveg_prod * (un / (un + speed * q_cdrag(ji) * rsol(ji))), &
!               & un - SUM(vbeta2(ji,:)+vbeta3(ji,:)))
          
       END DO
    ELSE
       DO ji = 1, kjpindex

          ! The limitation by 1-beta2-beta3 is due to the fact that evaporation under vegetation is possible
          !! \latexonly 
          !!     \input{diffucobare3.tex}
          !! \endlatexonly
          vbeta4(ji) = MIN(evap_bare_lim(ji), un - SUM(vbeta2(ji,:)+vbeta3(ji,:)))
          DO jv = 1,nvm
              IF ( veget_max(ji,jv) .GT. 0 ) THEN
!                  vbeta4_pft(ji,jv) = MIN(evap_bare_lim_pft(ji,jv), veget_max(ji,jv) - (vbeta2(ji,jv)+vbeta3(ji,jv)))
                  vbeta4_pft(ji,jv) = MIN(evap_bare_lim(ji), veget_max(ji,jv) - (vbeta2(ji,jv)+vbeta3(ji,jv)))
              ELSE
                  vbeta4_pft(ji,jv) = 0
              ENDIF
          ENDDO
       END DO
    ENDIF
    
    IF (printlev>=3) WRITE (numout,*) ' diffuco_bare done '
    
  END SUBROUTINE diffuco_bare


!! ================================================================================================================================
!! SUBROUTINE	: diffuco_trans 
!!
!>\BRIEF        This routine computes the partial beta coefficient 
!! corresponding to transpiration for each vegetation type.
!!
!! DESCRIPTION  : Beta coefficients for transpiration are calculated 
!! here using Jarvis formulation for stomatal resistance and
!! the structural resistance to represent the vertical gradient of 
!! transpiration within the canopy. \n
!!
!! The Jarvis formulation as used here is derived by Lohanner et al. (1980) from Jarvis (1976). This formulation is
!! semi-empirical: \n
!!
!! \latexonly 
!!     \input{diffucotrans4.tex}
!! \endlatexonly
!! \n
!!            
!! where in this expression LAI is the single sided Leaf Area Index, R_{new}^{SW} the net shortwave radiation,
!! R_{SO} the half light saturation factor, \delta c the water vapour concentration deficit, and a, k_0 and \lambda
!! are all parameters that are derived from extensive measurement of surface layer vegetation. \n
!!
!! Structural resistance (or architectural resistance) is a function of vegetation type and is assigned based on the
!! particular Plant Functional Type (PFT) in question. The range of values for the structural resistance are listed
!! in the module 'pft_parameters', and are described in de Noblet-Ducoudré et al (1993). \n
!!
!! vbetaco2 is here set to zero as this way to compute canopy resistances is only used 
!! without STOMATE, and there is therefore no photosynthesis. \n
!!
!! Moisture concentration at the leaf level.
!! \latexonly 
!!     \input{diffucotrans1.tex}
!! \endlatexonly
!!   
!! Calulation of the beta coefficient for vegetation transpiration, beta_3.
!! \latexonly 
!!     \input{diffucotrans2.tex}
!! \endlatexonly
!! \latexonly             
!!     \input{diffucotrans3.tex}
!! \endlatexonly
!!
!! \latexonly 
!!     \input{diffucotrans4.tex}
!! \endlatexonly            
!!
!! This is the formulation for beta_3.
!! \latexonly 
!!     \input{diffucotrans5.tex}
!! \endlatexonly
!!
!! \latexonly 
!!     \input{diffucotrans6.tex}
!! \endlatexonly
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::vbeta3, ::rveget, ::cimean and ::vbetaco2
!!
!! REFERENCE(S) :
!! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
!! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
!! Circulation Model. Journal of Climate, 6, pp.248-273
!! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
!! sur le cycle de l'eau, PhD Thesis, available from:
!! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
!! - Jarvis, PG, 1976. The interpretation of the variations in leaf water potential and stomatal
!! conductance found in canopies in the fields. Philosophical Transactions of the Royal Society of
!! London, Series B, 273, pp. 593-610
!! - Lohammer T, Larsson S, Linder S & Falk SO, 1980. Simulation models of gaseous exchange in Scotch
!! pine. Structure and function of Northern Coniferous Forest, Ecological Bulletin, 32, pp. 505-523
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE diffuco_trans (kjpindex, swnet, temp_air, pb, qair, rau, u, v, q_cdrag, humrel, &
                            veget, veget_max, lai, qsintveg, qsintmax, vbeta3, vbeta3pot, rveget, rstruct, &
                            cimean, vbetaco2, vbeta23)  

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    
    INTEGER(i_std), INTENT(in)                         :: kjpindex   !! Domain size (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: swnet      !! Short wave net flux at surface (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: temp_air   !! Air temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: pb         !! Lowest level pressure (hPa)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: qair       !! Lowest level specific air humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: rau        !! Air Density (kg m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: u          !! Eastward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: v          !! Northward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: q_cdrag    !! Surface drag coefficient  (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: humrel     !! Soil moisture stress (within range 0 to 1)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget      !! Type of vegetation (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget_max  !! Max. vegetation fraction (LAI->infty) (fraction)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: lai        !! Leaf area index (m^2 m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: qsintveg   !! Water on vegetation due to interception (kg m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: qsintmax   !! Maximum water on vegetation (kg m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: rstruct    !! Structural resistance (s m^{-1})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: vbeta23    !! Beta for wetted foliage fraction that will transpire (-)
    
    !! 0.2 Output variables
    
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: vbeta3     !! Beta for Transpiration (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: vbeta3pot  !! Beta for Potential Transpiration
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: rveget     !! Stomatal resistance of the whole canopy (s m^{-1})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: cimean     !! STOMATE: mean intercellular ci (see enerbil) 
                                                                     !! (\mumol m^{-2} s^{-1}) 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: vbetaco2   !! STOMATE: Beta for CO2 (-)

    !! 0.3 Modified variables
  
    !! 0.4 Local variables

    INTEGER(i_std)                                     :: ji, jv
    REAL(r_std)                                        :: speed
    REAL(r_std), DIMENSION(kjpindex)                   :: zdefconc, zqsvegrap
    REAL(r_std), DIMENSION(kjpindex)                   :: qsatt
    REAL(r_std),DIMENSION (kjpindex,nvm)               :: rveget_min           !! Minimal Surface resistance of vegetation
!_ ================================================================================================================================


  !! 1.  Moisture concentration at the leaf level.
    
    CALL qsatcalc (kjpindex, temp_air, pb, qsatt)
    
    !! \latexonly 
    !!     \input{diffucotrans1.tex}
    !! \endlatexonly
    zdefconc(:) = rau(:) * MAX( qsatt(:) - qair(:), zero )

   
  !! 2. Calulation of the beta coefficient for vegetation transpiration, beta_3.

    rveget(:,:) = undef_sechiba
    rveget_min(:,:) = undef_sechiba
    vbeta3(:,:) = zero
    vbeta3pot(:,:) = zero

    DO jv = 2,nvm

      zqsvegrap(:) = zero

      DO ji = 1, kjpindex

         !! \latexonly 
         !!     \input{diffucotrans2.tex}
         !! \endlatexonly
         speed = MAX(min_wind, wind(ji))

         IF (qsintmax(ji,jv) .GT. min_sechiba) THEN
        
            !! \latexonly 
            !!     \input{diffucotrans3.tex}
            !! \endlatexonly
            zqsvegrap(ji) = MAX(zero, qsintveg(ji,jv) / qsintmax(ji,jv))
         ENDIF
         
         IF ( ( veget(ji,jv)*lai(ji,jv) .GT. min_sechiba ) .AND. &
              ( kzero(jv) .GT. min_sechiba ) .AND. &
              ( swnet(ji) .GT. min_sechiba ) ) THEN

            !! \latexonly 
            !!     \input{diffucotrans4.tex}
            !! \endlatexonly            
            rveget(ji,jv) = (( swnet(ji) + rayt_cste ) / swnet(ji) ) &
                 * ((defc_plus + (defc_mult * zdefconc(ji) )) / kzero(jv)) * (un / lai(ji,jv))

            rveget_min(ji,jv) = (defc_plus / kzero(jv)) * (un / lai(ji,jv))

            ! Corrections Nathalie - le 28 mars 2006 - sur conseils Fred Hourdin
            ! Introduction d'un potentiometre (rveg_pft) pour regler la somme rveg+rstruct
            ! vbeta3(ji,jv) = veget(ji,jv) * (un - zqsvegrap(ji)) * humrel(ji,jv) * &
            !     (un / (un + speed * q_cdrag(ji) * (rveget(ji,jv) + rstruct(ji,jv))))
            
            !! This is the formulation for $beta_3$.
            !! \latexonly 
            !!     \input{diffucotrans5.tex}
            !! \endlatexonly
            vbeta3(ji,jv) = veget(ji,jv) * (un - zqsvegrap(ji)) * humrel(ji,jv) * &
                 (un / (un + speed * q_cdrag(ji) * (rveg_pft(jv)*(rveget(ji,jv) + rstruct(ji,jv)))))
            
            ! Fin ajout Nathalie
            ! Ajout Nathalie - Juin 2006

            !! \latexonly 
            !!     \input{diffucotrans6.tex}
            !! \endlatexonly
            vbeta3(ji,jv) = vbeta3(ji,jv) + MIN( vbeta23(ji,jv), &
                 veget(ji,jv) * zqsvegrap(ji) * humrel(ji,jv) * &
                 (un / (un + speed * q_cdrag(ji) * (rveg_pft(jv)*(rveget(ji,jv) + rstruct(ji,jv))))))
            ! Fin ajout Nathalie
            ! Autre possibilite permettant la transpiration sur toute la canopee
            !Commenter si formulation Nathalie sinon Tristan
!            vbeta3(ji,jv) = MAX(zero, MIN(vbeta23(ji,jv), &
!                 & veget_max(ji,jv) * humrel(ji,jv) / &
!                 & (un + speed * q_cdrag(ji) * (rveg_pft(jv)*(rveget(ji,jv) + rstruct(ji,jv))))))

           ! vbeta3pot for computation of potential transpiration (needed for irrigation)
            vbeta3pot(ji,jv) = &
                 &  MAX(zero, veget_max(ji,jv) / &
                 & (un + speed * q_cdrag(ji) * (rveg_pft(jv)*(rveget_min(ji,jv) + rstruct(ji,jv)))))
         ENDIF

      ENDDO

    ENDDO

    ! STOMATE
    cimean(:,:) = zero
    vbetaco2(:,:) = zero

    IF (printlev>=3) WRITE (numout,*) ' diffuco_trans done '

  END SUBROUTINE diffuco_trans


!! ==============================================================================================================================
!! SUBROUTINE   : diffuco_trans_co2
!!
!>\BRIEF        This subroutine computes carbon assimilation and stomatal 
!! conductance, following respectively Farqhuar et al. (1980) and Ball et al. (1987).
!!
!! DESCRIPTION  :\n
!! *** General:\n 
!! The equations are different depending on the photosynthesis mode (C3 versus C4).\n 
!! Assimilation and conductance are computed over 20 levels of LAI and then 
!! integrated at the canopy level.\n 
!! This routine also computes partial beta coefficient: transpiration for each 
!! type of vegetation.\n
!! There is a main loop on the PFTs, then inner loops on the points where 
!! assimilation has to be calculated.\n
!! This subroutine is called by diffuco_main only if photosynthesis is activated
!! for sechiba (flag STOMATE_OK_CO2=TRUE), otherwise diffuco_trans is called.\n
!! This subroutine is called at each sechiba time step by sechiba_main.\n
!! *** Details:
!! - Integration at the canopy level\n
!! \latexonly
!! \input{diffuco_trans_co2_1.1.tex}
!! \endlatexonly
!! - Light''s extinction \n
!! The available light follows a simple Beer extinction law. 
!! The extinction coefficients (ext_coef) are PFT-dependant constants and are defined in constant_co2.f90.\n
!! \latexonly
!! \input{diffuco_trans_co2_1.2.tex}
!! \endlatexonly
!! - Estimation of relative humidity of air (for calculation of the stomatal conductance)\n
!! \latexonly
!! \input{diffuco_trans_co2_1.3.tex}
!! \endlatexonly
!! - Calculation of the water limitation factor\n
!! \latexonly
!! \input{diffuco_trans_co2_2.1.tex}
!! \endlatexonly
!! - Calculation of temperature dependent parameters for C4 plants\n
!! \latexonly
!! \input{diffuco_trans_co2_2.2.tex}
!! \endlatexonly
!! - Calculation of temperature dependent parameters for C3 plants\n
!! \latexonly
!! \input{diffuco_trans_co2_2.3.tex}
!! \endlatexonly
!! - Vmax scaling\n 
!! Vmax is scaled into the canopy due to reduction of nitrogen 
!! (Johnson and Thornley,1984).\n
!! \latexonly
!! \input{diffuco_trans_co2_2.4.1.tex}
!! \endlatexonly
!! - Assimilation for C4 plants (Collatz et al., 1992)\n
!! \latexonly
!! \input{diffuco_trans_co2_2.4.2.tex}
!! \endlatexonly         
!! - Assimilation for C3 plants (Farqhuar et al., 1980)\n
!! \latexonly
!! \input{diffuco_trans_co2_2.4.3.tex}
!! \endlatexonly
!! - Estimation of the stomatal conductance (Ball et al., 1987)\n
!! \latexonly
!! \input{diffuco_trans_co2_2.4.4.tex}
!! \endlatexonly
!!
!! RECENT CHANGE(S): N. de Noblet          2006/06
!!                - addition of q2m and t2m as input parameters for the 
!!                calculation of Rveget
!!                - introduction of vbeta23
!!
!! MAIN OUTPUT VARIABLE(S): beta coefficients, resistances, CO2 intercellular 
!! concentration
!!
!! REFERENCE(S) :
!! - Ball, J., T. Woodrow, and J. Berry (1987), A model predicting stomatal 
!! conductance and its contribution to the control of photosynthesis under 
!! different environmental conditions, Prog. Photosynthesis, 4, 221– 224.
!! - Collatz, G., M. Ribas-Carbo, and J. Berry (1992), Coupled photosynthesis 
!! stomatal conductance model for leaves of C4 plants, Aust. J. Plant Physiol.,
!! 19, 519–538.
!! - Farquhar, G., S. von Caemmener, and J. Berry (1980), A biochemical model of 
!! photosynthesis CO2 fixation in leaves of C3 species, Planta, 149, 78–90.
!! - Johnson, I. R., and J. Thornley (1984), A model of instantaneous and daily
!! canopy photosynthesis, J Theor. Biol., 107, 531 545
!! - McMurtrie, R.E., Rook, D.A. and Kelliher, F.M., 1990. Modelling the yield of Pinus radiata on a
!! site limited by water and nitrogen. For. Ecol. Manage., 30: 381-413
!! - Bounoua, L., Hall, F. G., Sellers, P. J., Kumar, A., Collatz, G. J., Tucker, C. J., and Imhoff, M. L. (2010), Quantifying the 
!! negative feedback of vegetation to greenhouse warming: A modeling approach, Geophysical Research Letters, 37, Artn L23701, 
!! Doi 10.1029/2010gl045338
!! - Bounoua, L., Collatz, G. J., Sellers, P. J., Randall, D. A., Dazlich, D. A., Los, S. O., Berry, J. A., Fung, I., 
!! Tucker, C. J., Field, C. B., and Jensen, T. G. (1999), Interactions between vegetation and climate: Radiative and physiological 
!! effects of doubled atmospheric co2, Journal of Climate, 12, 309-324, Doi 10.1175/1520-0442(1999)012<0309:Ibvacr>2.0.Co;2
!! - Sellers, P. J., Bounoua, L., Collatz, G. J., Randall, D. A., Dazlich, D. A., Los, S. O., Berry, J. A., Fung, I., 
!! Tucker, C. J., Field, C. B., and Jensen, T. G. (1996), Comparison of radiative and physiological effects of doubled atmospheric
!! co2 on climate, Science, 271, 1402-1406, DOI 10.1126/science.271.5254.1402
!! - Lewis, J. D., Ward, J. K., and Tissue, D. T. (2010), Phosphorus supply drives nonlinear responses of cottonwood 
!! (populus deltoides) to increases in co2 concentration from glacial to future concentrations, New Phytologist, 187, 438-448, 
!! DOI 10.1111/j.1469-8137.2010.03307.x
!! - Kattge, J., Knorr, W., Raddatz, T., and Wirth, C. (2009), Quantifying photosynthetic capacity and its relationship to leaf 
!! nitrogen content for global-scale terrestrial biosphere models, Global Change Biology, 15, 976-991, 
!! DOI 10.1111/j.1365-2486.2008.01744.x
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

SUBROUTINE diffuco_trans_co2 (kjpindex, lalo, swdown, pb, qsurf, q2m, t2m, temp_growth, rau, u, v, q_cdrag, q_cdrag_pft,  humrel, &
                                assim_param, Ca, &
                                veget, veget_max, lai, qsintveg, qsintmax, vbeta3, vbeta3pot, rveget, rstruct, &
                                cimean, gsmean, gpp, vbeta23, hist_id, indexveg, indexlai, indexlai0, & 
                                index, kjit, cim, wtp_year) !!!qcj++ peatland 

    !
    !! 0. Variable and parameter declaration
    !

    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                               :: kjpindex         !! Domain size (unitless)
    REAL(r_std), DIMENSION(kjpindex,2), INTENT(in)           :: lalo             !! Geographical coordinates for pixels (degrees)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: swdown           !! Downwelling short wave flux 
                                                                                 !! @tex ($W m^{-2}$) @endtex 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: pb               !! Lowest level pressure (hPa)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: qsurf             !! Near surface specific humidity 
                                                                                 !! @tex ($kg kg^{-1}$) @endtex
! N. de Noblet - 2006/06 - addition of q2m and t2m
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: q2m              !! 2m specific humidity 
                                                                                 !! @tex ($kg kg^{-1}$) @endtex
! In off-line mode q2m and qair are the same.
! In off-line mode t2m and temp_air are the same.
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: t2m              !! 2m air temperature (K)
! N. de Noblet - 2006/06 - addition of q2m and t2m
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: temp_growth      !! Growth temperature (°C) - Is equal to t2m_month
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: rau              !! air density @tex ($kg m^{-3}$) @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: u                !! Lowest level wind speed 
                                                                                 !! @tex ($m s^{-1}$) @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: v                !! Lowest level wind speed 
                                                                                 !! @tex ($m s^{-1}$) @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: q_cdrag          !! Surface drag coefficient  (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: q_cdrag_pft      !!
    REAL(r_std),DIMENSION (kjpindex,nvm,npco2), INTENT (in)  :: assim_param      !! min+max+opt temps (K), vcmax, vjmax for 
                                                                                 !! photosynthesis 
                                                                                 !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: Ca               !! CO2 concentration inside the canopy
                                                                                 !! @tex ($\mu mol mol^{-1}$) @endtex
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: humrel           !! Soil moisture stress (0-1,unitless) 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: veget            !! Coverage fraction of vegetation for each PFT 
                                                                                 !! depending on LAI (0-1, unitless) 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: veget_max        !! Maximum vegetation fraction of each PFT inside 
                                                                                 !! the grid box (0-1, unitless) 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: lai              !! Leaf area index @tex ($m^2 m^{-2}$) @endtex
                                                                                 !! @tex ($m s^{-1}$) @endtex
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: qsintveg         !! Water on vegetation due to interception 
                                                                                 !! @tex ($kg m^{-2}$) @endte
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: qsintmax         !! Maximum water on vegetation
                                                                                 !! @tex ($kg m^{-2}$) @endtex
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: vbeta23          !! Beta for fraction of wetted foliage that will 
                                                                                 !! transpire (unitless) 
    INTEGER(i_std),INTENT (in)                               :: hist_id          !! _History_ file identifier (-)    
    INTEGER(i_std),DIMENSION (kjpindex*nvm), INTENT (in) :: indexveg       !! Indeces of the points on the 3D map (-)   
    INTEGER(i_std),DIMENSION (kjpindex*(nlai+1)), INTENT (in) :: indexlai  !! Indeces of the points on the 3D map
    INTEGER(i_std),DIMENSION (kjpindex*(nlai)), INTENT (in) :: indexlai0  !! Indeces of the points on the 3D map
INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: index            !! Indeces of the points on the map (-)
    INTEGER(i_std), INTENT(in)                         :: kjit             !! Time step number (-)        

!!!qcj++ peatland
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: wtp_year    
  
    !
    !! 0.2 Output variables
    !
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: vbeta3           !! Beta for Transpiration (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: vbeta3pot        !! Beta for Potential Transpiration
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: rveget           !! stomatal resistance of vegetation 
                                                                                 !! @tex ($s m^{-1}$) @endtex
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: rstruct          !! structural resistance @tex ($s m^{-1}$) @endtex
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: cimean           !! mean intercellular CO2 concentration 
                                                                                 !! @tex ($\mu mol mol^{-1}$) @endtex
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: gsmean           !! mean stomatal conductance to CO2 (umol m-2 s-1)
    REAL(r_Std),DIMENSION (kjpindex,nvm), INTENT (out)       :: gpp
 !! Assimilation ((gC m^{-2} s^{-1}), total area)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: cim              !! Intercellular CO2 over nlai 
    !
    !! 0.3 Modified variables
    !
 
    !
    !! 0.4 Local variables
    !
    REAL(r_std),DIMENSION (kjpindex,nvm)  :: vcmax                               !! maximum rate of carboxylation 
                                                                                 !! @tex ($\mu mol CO2 m^{-2} s^{-1}$) @endtex
!!!qcj++ peatland
    REAL(r_std),DIMENSION (kjpindex,nvm)  :: vcmax_tmp
    REAL(r_std),DIMENSION (kjpindex,nvm)  :: vcmax_wt

    INTEGER(i_std)                        :: ji, jv, jl, limit_photo             !! indices (unitless)
    REAL(r_std), DIMENSION(kjpindex,nlai+1) :: info_limitphoto
    REAL(r_std), DIMENSION(kjpindex,nvm,nlai)  :: leaf_ci                        !! intercellular CO2 concentration (ppm)
    REAL(r_std), DIMENSION(kjpindex)      :: leaf_ci_lowest                      !! intercellular CO2 concentration at the lowest 
                                                                                 !! LAI level
                                                                                 !! @tex ($\mu mol mol^{-1}$) @endtex
    INTEGER(i_std), DIMENSION(kjpindex)   :: ilai                                !! counter for loops on LAI levels (unitless)
    REAL(r_std), DIMENSION(kjpindex)      :: zqsvegrap                           !! relative water quantity in the water 
                                                                                 !! interception reservoir (0-1,unitless) 
    REAL(r_std)                           :: speed                               !! wind speed @tex ($m s^{-1}$) @endtex
    ! Assimilation
    LOGICAL, DIMENSION(kjpindex)          :: assimilate                          !! where assimilation is to be calculated 
                                                                                 !! (unitless) 
    LOGICAL, DIMENSION(kjpindex)          :: calculate                           !! where assimilation is to be calculated for 
                                                                                 !! in the PFTs loop (unitless) 
    INTEGER(i_std)                        :: nic,inic,icinic                     !! counter/indices (unitless)
    INTEGER(i_std), DIMENSION(kjpindex)   :: index_calc                          !! index (unitless)
    INTEGER(i_std)                        :: nia,inia,nina,inina,iainia          !! counter/indices (unitless)
    INTEGER(i_std), DIMENSION(kjpindex)   :: index_assi,index_non_assi           !! indices (unitless)
    REAL(r_std), DIMENSION(kjpindex, nlai+1)      :: vc2                                 !! rate of carboxylation (at a specific LAI level) 
                                                                                 !! @tex ($\mu mol CO2 m^{-2} s^{-1}$) @endtex 
    REAL(r_std), DIMENSION(kjpindex, nlai+1)      :: vj2                                 !! rate of Rubisco regeneration (at a specific LAI 
                                                                                 !! level) @tex ($\mu mol e- m^{-2} s^{-1}$) @endtex 
    REAL(r_std), DIMENSION(kjpindex, nlai+1)      :: assimi                              !! assimilation (at a specific LAI level) 
                                                                                 !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
                                                                                 !! (temporary variables)
    REAL(r_std), DIMENSION(kjpindex)             :: gstop                        !! stomatal conductance to H2O at topmost level 
                                                                                 !! @tex ($m s^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex,nlai+1)      :: gs                                  !! stomatal conductance to CO2 
                                                                                 !! @tex ($\mol m^{-2} s^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex,nlai)      :: templeafci  


    REAL(r_std), DIMENSION(kjpindex)      :: gamma_star                          !! CO2 compensation point (ppm)
                                                                                 !! @tex ($\mu mol mol^{-1}$) @endtex

    REAL(r_std), DIMENSION(kjpindex)      :: air_relhum                          !! air relative humidity at 2m 
                                                                                 !! @tex ($kg kg^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex)      :: VPD                                 !! Vapor Pressure Deficit (kPa)
    REAL(r_std), DIMENSION(kjpindex)      :: water_lim                           !! water limitation factor (0-1,unitless)

    REAL(r_std), DIMENSION(kjpindex)      :: gstot                               !! total stomatal conductance to H2O
                                                                                 !! Final unit is
                                                                                 !! @tex ($m s^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex)      :: assimtot                            !! total assimilation 
                                                                                 !! @tex ($\mu mol CO2 m^{-2} s^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex)      :: Rdtot                               !! Total Day respiration (respiratory CO2 release other than by photorespiration) (mumol CO2 m−2 s−1) 
    REAL(r_std), DIMENSION(kjpindex)      :: leaf_gs_top                         !! leaf stomatal conductance to H2O at topmost level 
                                                                                 !! @tex ($\mol H2O m^{-2} s^{-1}$) @endtex
    REAL(r_std), DIMENSION(nlai+1)        :: laitab                              !! tabulated LAI steps @tex ($m^2 m^{-2}$) @endtex
    REAL(r_std), DIMENSION(kjpindex)      :: qsatt                               !! surface saturated humidity at 2m (??) 
                                                                                 !! @tex ($g g^{-1}$) @endtex
    REAL(r_std), DIMENSION(nvm,nlai)      :: light                               !! fraction of light that gets through upper LAI    
                                                                                 !! levels (0-1,unitless)
    REAL(r_std), DIMENSION(kjpindex)      :: T_Vcmax                             !! Temperature dependance of Vcmax (unitless)
    REAL(r_std), DIMENSION(kjpindex)      :: S_Vcmax_acclim_temp                 !! Entropy term for Vcmax 
                                                                                 !! accounting for acclimation to temperature (J K-1 mol-1)
    REAL(r_std), DIMENSION(kjpindex)      :: T_Jmax                              !! Temperature dependance of Jmax
    REAL(r_std), DIMENSION(kjpindex)      :: S_Jmax_acclim_temp                  !! Entropy term for Jmax 
                                                                                 !! accounting for acclimation toxs temperature (J K-1 mol-1)
    REAL(r_std), DIMENSION(kjpindex)      :: T_gm                                !! Temperature dependance of gmw
    REAL(r_std), DIMENSION(kjpindex)      :: T_Rd                                !! Temperature dependance of Rd (unitless)
    REAL(r_std), DIMENSION(kjpindex)      :: T_Kmc                               !! Temperature dependance of KmC (unitless)
    REAL(r_std), DIMENSION(kjpindex)      :: T_KmO                               !! Temperature dependance of KmO (unitless)
    REAL(r_std), DIMENSION(kjpindex)      :: T_Sco                               !! Temperature dependance of Sco
    REAL(r_std), DIMENSION(kjpindex)      :: T_gamma_star                        !! Temperature dependance of gamma_star (unitless)    
    REAL(r_std), DIMENSION(kjpindex)      :: vc                                  !! Maximum rate of Rubisco activity-limited carboxylation (mumol CO2 m−2 s−1)
    REAL(r_std), DIMENSION(kjpindex)      :: vj                                  !! Maximum rate of e- transport under saturated light (mumol CO2 m−2 s−1)
    REAL(r_std), DIMENSION(kjpindex)      :: gm                                  !! Mesophyll diffusion conductance (molCO2 m−2 s−1 bar−1)
    REAL(r_std), DIMENSION(kjpindex)      :: g0var
    REAL(r_std), DIMENSION(kjpindex,nlai+1)      :: Rd                                  !! Day respiration (respiratory CO2 release other than by photorespiration) (mumol CO2 m−2 s−1)
    REAL(r_std), DIMENSION(kjpindex)      :: Kmc                                 !! Michaelis–Menten constant of Rubisco for CO2 (mubar)
    REAL(r_std), DIMENSION(kjpindex)      :: KmO                                 !! Michaelis–Menten constant of Rubisco for O2 (mubar)
    REAL(r_std), DIMENSION(kjpindex)      :: Sco                                 !! Relative CO2 /O2 specificity factor for Rubisco (bar bar-1)
    REAL(r_std), DIMENSION(kjpindex)      :: gb_co2                              !! Boundary-layer conductance (molCO2 m−2 s−1 bar−1)
    REAL(r_std), DIMENSION(kjpindex)      :: gb_h2o                              !! Boundary-layer conductance (molH2O m−2 s−1 bar−1)
    REAL(r_std), DIMENSION(kjpindex)      :: fvpd                                !! Factor for describing the effect of leaf-to-air vapour difference on gs (-)
    REAL(r_std), DIMENSION(kjpindex)      :: low_gamma_star                      !! Half of the reciprocal of Sc/o (bar bar-1)
    REAL(r_std)                           :: N_Vcmax                             !! Nitrogen level dependance of Vcmacx and Jmax 
    REAL(r_std)                           :: fcyc                                !! Fraction of electrons at PSI that follow cyclic transport around PSI (-)
    REAL(r_std)                           :: z                                   !! A lumped parameter (see Yin et al. 2009) ( mol mol-1)                          
    REAL(r_std)                           :: Rm                                  !! Day respiration in the mesophyll (umol CO2 m−2 s−1)
    REAL(r_std)                           :: Cs_star                             !! Cs -based CO2 compensation point in the absence of Rd (ubar)
    REAL(r_std), DIMENSION(kjpindex)      :: Iabs                                !! Photon flux density absorbed by leaf photosynthetic pigments (umol photon m−2 s−1)
    REAL(r_std), DIMENSION(kjpindex)      :: Jmax                                !! Maximum value of J under saturated light (umol e− m−2 s−1)
    REAL(r_std), DIMENSION(kjpindex, nlai+1)      :: JJ                                  !! Rate of e− transport (umol e− m−2 s−1)
    REAL(r_std)                           :: J2                                  !! Rate of all e− transport through PSII (umol e− m−2 s−1)
    REAL(r_std)                           :: VpJ2                                !! e− transport-limited PEP carboxylation rate (umol CO2 m−2 s−1)
    REAL(r_std)                           :: A_1, A_3                            !! Lowest First and third roots of the analytical solution for a general cubic equation (see Appendix A of Yin et al. 2009) (umol CO2 m−2 s−1)
    REAL(r_std)                           :: A_1_tmp, A_3_tmp                            !! Temporary First and third roots of the analytical solution for a general cubic equation (see Appendix A of Yin et al. 2009) (umol CO2 m−2 s−1)
    REAL(r_std)                           :: Obs                                 !! Bundle-sheath oxygen partial pressure (ubar)
    REAL(r_std), DIMENSION(kjpindex, nlai+1)      :: Cc                                  !! Chloroplast CO2 partial pressure (ubar)
    REAL(r_std)                           :: ci_star                             !! Ci -based CO2 compensation point in the absence of Rd (ubar)        
    REAL(r_std)                           :: a,b,c,d,m,f,j,g,h,i,l,p,q,r         !! Variables used for solving the cubic equation (see Yin et al. (2009))
    REAL(r_std)                           :: QQ,UU,PSI,x1,x2,x3                      !! Variables used for solving the cubic equation (see Yin et al. (2009))
                                        
    REAL(r_std)                           :: cresist                             !! coefficient for resistances (??)
    CHARACTER(LEN=150)                    :: printstr, printstr2                 !! For temporary uses
    REAL(r_std), DIMENSION(kjpindex)      :: laisum                              !! when calculating cim over nlai

! @defgroup Photosynthesis Photosynthesis
! @{   
    ! 1. Preliminary calculations\n
!_ ================================================================================================================================

    cim(:,:)=zero
    leaf_ci = zero

    !
    ! 1.1 Calculate LAI steps\n
    ! The integration at the canopy level is done over nlai fixed levels.
    !! \latexonly
    !! \input{diffuco_trans_co2_1.1.tex}
    !! \endlatexonly
! @}
! @codeinc
    DO jl = 1, nlai+1
      laitab(jl) = laimax*(EXP(lai_level_depth*REAL(jl-1,r_std))-1.)/(EXP(lai_level_depth*REAL(nlai,r_std))-un)
    ENDDO
! @endcodeinc

! @addtogroup Photosynthesis
! @{   
    !
    ! 1.2 Calculate light fraction for each LAI step\n
    ! The available light follows a simple Beer extinction law. 
    ! The extinction coefficients (ext_coef) are PFT-dependant constants and are defined in constant_co2.f90.
    !! \latexonly
    !! \input{diffuco_trans_co2_1.2.tex}
    !! \endlatexonly
! @}
! @codeinc
    DO jl = 1, nlai
      DO jv = 1, nvm
        light(jv,jl) = exp( -ext_coeff(jv)*laitab(jl) )
      ENDDO
    ENDDO 
! @endcodeinc
    !
    ! Photosynthesis parameters
    !

    ! Choice of downregulation. Note that downregulation_co2_new excludes
    IF (downregulation_co2_new) THEN
       ! Option used for CMIP6 version from 6.1.11
       ! A minimum value is used for the CO2 concentration(Ca)
       ! For low CO2 concentration values(under downregulation_co2_minimum) the
       ! parametrization allows a behavior in the same way as for the
       ! preindustral period.
       DO jv= 1, nvm
           vcmax(:,jv) = assim_param(:,jv,ivcmax)*(un-downregulation_co2_coeff_new(jv) * &
                         (MAX(Ca(:),downregulation_co2_minimum)-downregulation_co2_baselevel) / &
                         (MAX(Ca(:),downregulation_co2_minimum)+20.))
       ENDDO
    ELSE IF (downregulation_co2) THEN
       ! Option used for CMIP6 version 6.1.0 up to 6.1.10
       DO jv= 1, nvm
          vcmax(:,jv) = assim_param(:,jv,ivcmax)*(un-downregulation_co2_coeff(jv) * &
                          log(Ca(:)/downregulation_co2_baselevel))
       ENDDO
    ELSE
       vcmax(:,:) = assim_param(:,:,ivcmax)
    ENDIF

!!!qcj++ peatland

   IF (wtp_constraint) THEN
      vcmax_tmp=vcmax(:,:)
      DO jv= 1, nvm
!mosses and tropical peatland trees not constrainted by wtp
         IF ( is_peat(jv) .AND. (.NOT. is_mosspeat(jv)) .AND. (.NOT. (leaf_tab(jv)==1 .AND. pheno_type(jv)==1) ) .AND. (.NOT. (leaf_tab(jv)==1 .AND. pheno_type(jv)==3) ) ) THEN
            vcmax_wt(:,jv)=wtp_year(:,jv)*factor_wtp
            WHERE (vcmax_wt(:,jv) .GT. 0.5*factor_wtp) ! If WT>0.5, no constraint on Vcmax 
                vcmax_wt(:,jv) = 0.5*factor_wtp
            ENDWHERE   
            vcmax(:,jv) = vcmax_tmp(:,jv) + vcmax_wt(:,jv)
         ENDIF
      ENDDO
   ENDIF

!    DO jv = 1, nvm
!       vcmax(:,:) = Vcmax25(jv)
!    ENDDO

! @addtogroup Photosynthesis
! @{   
    !
    ! 1.3 Estimate relative humidity of air (for calculation of the stomatal conductance).\n
    !! \latexonly
    !! \input{diffuco_trans_co2_1.3.tex}
    !! \endlatexonly
! @}
    !
! N. de Noblet - 2006/06 - We use q2m/t2m instead of qair.
!    CALL qsatcalc (kjpindex, temp_air, pb, qsatt)
!    air_relhum(:) = &
!      ( qair(:) * pb(:) / (0.622+qair(:)*0.378) ) / &
!      ( qsatt(:)*pb(:) / (0.622+qsatt(:)*0.378 ) )
! @codeinc
    CALL qsatcalc (kjpindex, t2m, pb, qsatt)
    air_relhum(:) = &
      ( qsurf(:) * pb(:) / (Tetens_1+qsurf(:)* Tetens_2) ) / &
      ( qsatt(:)*pb(:) / (Tetens_1+qsatt(:)*Tetens_2 ) )


    VPD(:) = ( qsatt(:)*pb(:) / (Tetens_1+qsatt(:)*Tetens_2 ) ) &
         - ( qsurf(:) * pb(:) / (Tetens_1+qsurf(:)* Tetens_2) )
    ! VPD is needed in kPa
    VPD(:) = VPD(:)/10.

! @endcodeinc
! N. de Noblet - 2006/06 
    !
    !
    ! 2. beta coefficient for vegetation transpiration
    !
    rstruct(:,1) = rstruct_const(1)
    rveget(:,:) = undef_sechiba
    !
    vbeta3(:,:) = zero
    vbeta3pot(:,:) = zero
    gsmean(:,:) = zero
    gpp(:,:) = zero
    !
    cimean(:,1) = Ca(:)
    !
    ! @addtogroup Photosynthesis
    ! @{   
    ! 2. Loop over vegetation types\n
    ! @} 
    !
    DO jv = 2,nvm
       gamma_star(:)=zero
       Kmo(:)=zero
       Kmc(:)=zero
       gm(:)=zero
       g0var(:) =zero

       Cc(:,:)=zero
       Vc2(:,:)=zero
       JJ(:,:)=zero
       info_limitphoto(:,:)=zero
       gs(:,:)=zero
       templeafci(:,:)=zero
       assimi(:,:)=zero
       Rd(:,:)=zero

      !
      ! @addtogroup Photosynthesis
      ! @{   
      !
      ! 2.1 Initializations\n
      !! \latexonly
      !! \input{diffuco_trans_co2_2.1.tex}
      !! \endlatexonly
      ! @}      
      !
      ! beta coefficient for vegetation transpiration
      !
      rstruct(:,jv) = rstruct_const(jv)
      cimean(:,jv) = Ca(:)
      !
      !! mask that contains points where there is photosynthesis
      !! For the sake of vectorisation [DISPENSABLE], computations are done only for convenient points.
      !! nia is the number of points where the assimilation is calculated and nina the number of points where photosynthesis is not
      !! calculated (based on criteria on minimum or maximum values on LAI, vegetation fraction, shortwave incoming radiation, 
      !! temperature and relative humidity).
      !! For the points where assimilation is not calculated, variables are initialized to specific values. 
      !! The assimilate(kjpindex) array contains the logical value (TRUE/FALSE) relative to this photosynthesis calculation.

      nia=0
      nina=0
      !
      DO ji=1,kjpindex
         !
!!         IF ( ( lai(ji,jv) .GT. 0.01 ) .AND. & !! original
!!         IF ( ( lai(ji,jv) .GT. 0.001 ) .AND. & !! I think this one is better
!in precision, but it will affect the general performance
           IF ( ( (ok_LAIdev(jv) .AND. (lai(ji,jv) .GT. 0.001) ) .OR. ( (.NOT. ok_LAIdev(jv)) .AND. (lai(ji,jv) .GT. 0.01) ) ) .AND. &
              ( veget_max(ji,jv) .GT. min_sechiba ) ) THEN

#ifdef STRICT_CHECK 
            IF (vcmax(ji, jv) <= 0) THEN
               WRITE(printstr, *)  'coordinates lat=', lalo(ji, 1),' long=', lalo(ji, 2), 'PFT=', jv
               WRITE(printstr2, *) 'lai=', lai(ji,jv), '  vcmax=', vcmax(ji,jv)
               CALL ipslerr_p(3, 'diffuco_trans_co2', 'vcmax must be bigger than 0 to be consistent with lai', & 
                              TRIM(printstr), TRIM(printstr2))
            ENDIF
#endif

            IF ( ( veget(ji,jv) .GT. min_sechiba ) .AND. &
                 ( swdown(ji) .GT. min_sechiba )   .AND. &
                 ( humrel(ji,jv) .GT. min_sechiba) .AND. &
                 ( temp_growth(ji) .GT. tphoto_min(jv) ) .AND. &
                 ( temp_growth(ji) .LT. tphoto_max(jv) ) ) THEN
               !
               assimilate(ji) = .TRUE.
               nia=nia+1
               index_assi(nia)=ji
               !
            ELSE
               !
               assimilate(ji) = .FALSE.
               nina=nina+1
               index_non_assi(nina)=ji
               !
            ENDIF
         ELSE
            !
            assimilate(ji) = .FALSE.
            nina=nina+1
            index_non_assi(nina)=ji
            !
         ENDIF
         !

      ENDDO
      !

      gstot(:) = zero
      gstop(:) = zero
      assimtot(:) = zero
      Rdtot(:)=zero
      leaf_gs_top(:) = zero
      !
      zqsvegrap(:) = zero
      WHERE (qsintmax(:,jv) .GT. min_sechiba)
      !! relative water quantity in the water interception reservoir
          zqsvegrap(:) = MAX(zero, qsintveg(:,jv) / qsintmax(:,jv))
      ENDWHERE
      !
      !! Calculates the water limitation factor.
      water_lim(:) = humrel(:,jv)

      ! give a default value of ci for all pixel that do not assimilate
      DO jl=1,nlai
         DO inina=1,nina
            leaf_ci(index_non_assi(inina),jv,jl) = Ca(index_non_assi(inina)) 
         ENDDO
      ENDDO
      !
      ilai(:) = 1
      !
      ! Here is the calculation of assimilation and stomatal conductance
      ! based on the work of Farquahr, von Caemmerer and Berry (FvCB model) 
      ! as described in Yin et al. 2009
      ! Yin et al. developed a extended version of the FvCB model for C4 plants
      ! and proposed an analytical solution for both photosynthesis pathways (C3 and C4)
      ! Photosynthetic parameters used are those reported in Yin et al. 
      ! Except For Vcmax25, relationships between Vcmax25 and Jmax25 for which we use 
      ! Medlyn et al. (2002) and Kattge & Knorr (2007)
      ! Because these 2 references do not consider mesophyll conductance, we neglect this term
      ! in the formulations developed by Yin et al. 
      ! Consequently, gm (the mesophyll conductance) tends to the infinite
      ! This is of importance because as stated by Kattge & Knorr and Medlyn et al.,
      ! values of Vcmax and Jmax derived with different model parametrizations are not 
      ! directly comparable and the published values of Vcmax and Jmax had to be standardized
      ! to one consistent formulation and parametrization

      ! See eq. 6 of Yin et al. (2009)
      ! Parametrization of Medlyn et al. (2002) - from Bernacchi et al. (2001)
      T_KmC(:)        = Arrhenius(kjpindex,t2m,298.,E_KmC(jv))
      T_KmO(:)        = Arrhenius(kjpindex,t2m,298.,E_KmO(jv))
      T_Sco(:)        = Arrhenius(kjpindex,t2m,298.,E_Sco(jv))
      T_gamma_star(:) = Arrhenius(kjpindex,t2m,298.,E_gamma_star(jv))


      ! Parametrization of Yin et al. (2009) - from Bernacchi et al. (2001)
      T_Rd(:)         = Arrhenius(kjpindex,t2m,298.,E_Rd(jv))


      ! For C3 plants, we assume that the Entropy term for Vcmax and Jmax 
      ! acclimates to temperature as shown by Kattge & Knorr (2007) - Eq. 9 and 10
      ! and that Jmax and Vcmax respond to temperature following a modified Arrhenius function
      ! (with a decrease of these parameters for high temperature) as in Medlyn et al. (2002) 
      ! and Kattge & Knorr (2007).
      ! In Yin et al. (2009), temperature dependance to Vcmax is based only on a Arrhenius function
      ! Concerning this apparent unconsistency, have a look to the section 'Limitation of 
      ! Photosynthesis by gm' of Bernacchi (2002) that may provide an explanation
      
      ! Growth temperature tested by Kattge & Knorr range from 11 to 35°C
      ! So, we limit the relationship between these lower and upper limits
      S_Jmax_acclim_temp(:) = aSJ(jv) + bSJ(jv) * MAX(11., MIN(temp_growth(:),35.))      
      T_Jmax(:)  = Arrhenius_modified(kjpindex,t2m,298.,E_Jmax(jv),D_Jmax(jv),S_Jmax_acclim_temp)

      S_Vcmax_acclim_temp(:) = aSV(jv) + bSV(jv) * MAX(11., MIN(temp_growth(:),35.))   
      T_Vcmax(:) = Arrhenius_modified(kjpindex,t2m,298.,E_Vcmax(jv),D_Vcmax(jv),S_Vcmax_acclim_temp)
       

      
      vc(:) = vcmax(:,jv) * T_Vcmax(:)
      ! As shown by Kattge & Knorr (2007), we make use
      ! of Jmax25/Vcmax25 ratio (rJV) that acclimates to temperature for C3 plants
      ! rJV is written as a function of the growth temperature
      ! rJV = arJV + brJV * T_month 
      ! See eq. 10 of Kattge & Knorr (2007)
      ! and Table 3 for Values of arJV anf brJV 
      ! Growth temperature is monthly temperature (expressed in °C) - See first paragraph of
      ! section Methods/Data of Kattge & Knorr
      vj(:) = ( arJV(jv) + brJV(jv) *  MAX(11., MIN(temp_growth(:),35.)) ) * vcmax(:,jv) * T_Jmax(:)

      T_gm(:)    = Arrhenius_modified(kjpindex,t2m,298.,E_gm(jv),D_gm(jv),S_gm(jv))
      gm(:) = gm25(jv) * T_gm(:) * MAX(1-stress_gm(jv), water_lim(:))

      g0var(:) = g0(jv)* MAX(1-stress_gs(jv), water_lim(:))
      ! @endcodeinc
      !
      KmC(:)=KmC25(jv)*T_KmC(:)
      KmO(:)=KmO25(jv)*T_KmO(:)
      Sco(:)=Sco25(jv)*T_sco(:)
      gamma_star(:) = gamma_star25(jv)*T_gamma_star(:)



      ! low_gamma_star is defined by Yin et al. (2009)
      ! as the half of the reciprocal of Sco - See Table 2
      low_gamma_star(:) = 0.5 / Sco(:)

      ! VPD expressed in kPa
      ! Note : MIN(1.-min_sechiba,MAX(min_sechiba,(a1(jv) - b1(jv) * VPD(:)))) is always between 0-1 not including 0 and 1
      fvpd(:) = 1. / ( 1. / MIN(1.-min_sechiba,MAX(min_sechiba,(a1(jv) - b1(jv) * VPD(:)))) - 1. ) &
                * MAX(1-stress_gs(jv), water_lim(:))

      ! leaf boundary layer conductance 
      ! conversion from a conductance in (m s-1) to (mol H2O m-2 s-1)
      ! from Pearcy et al. (1991, see below)
      gb_h2o(:) = gb_ref * 44.6 * (tp_00/t2m(:)) * (pb(:)/pb_std) 

      ! conversion from (mol H2O m-2 s-1) to (mol CO2 m-2 s-1)
      gb_co2(:) = gb_h2o(:) / ratio_H2O_to_CO2

      !
      ! @addtogroup Photosynthesis
      ! @{   
      !
      ! 2.4 Loop over LAI discretized levels to estimate assimilation and conductance\n
      ! @}           
      !
      !! The calculate(kjpindex) array is of type logical to indicate wether we have to sum over this LAI fixed level (the LAI of
      !! the point for the PFT is lower or equal to the LAI level value). The number of such points is incremented in nic and the 
      !! corresponding point is indexed in the index_calc array.
      JJ(:,:)=zero
      vc2(:,:)=zero
      vj2(:,:)=zero
      Cc(:,:)=zero
      gs(:,:)=zero
      assimi(:,:)=zero
      Rd(:,:)=zero

      DO jl = 1, nlai
         !
         nic=0
         calculate(:) = .FALSE.
         !
         IF (nia .GT. 0) then
            DO inia=1,nia
               calculate(index_assi(inia)) = (laitab(jl) .LE. lai(index_assi(inia),jv) )
               IF ( calculate(index_assi(inia)) ) THEN
                  nic=nic+1
                  index_calc(nic)=index_assi(inia)
               ENDIF
            ENDDO
         ENDIF
         !
         ! @addtogroup Photosynthesis
         ! @{   
         !
         ! 2.4.1 Vmax is scaled into the canopy due to reduction of nitrogen 
         !! (Johnson and Thornley,1984).\n
         !! \latexonly
         !! \input{diffuco_trans_co2_2.4.1.tex}
         !! \endlatexonly
         ! @}           
         !
         N_Vcmax = ( un - .7_r_std * ( un - light(jv,jl) ) )
         !

         vc2(:,jl) = vc(:) * N_Vcmax * MAX(1-stress_vcmax(jv), water_lim(:))
         vj2(:,jl) = vj(:) * N_Vcmax * MAX(1-stress_vcmax(jv), water_lim(:))

         ! see Comment in legend of Fig. 6 of Yin et al. (2009)
         ! Rd25 is assumed to equal 0.01 Vcmax25 
         Rd(:,jl) = vcmax(:,jv) * N_Vcmax * 0.01 * T_Rd(:)  * MAX(1-stress_vcmax(jv), water_lim(:))

         Iabs(:)=swdown(:)*W_to_mol*RG_to_PAR*ext_coeff(jv)*light(jv,jl)
         
         ! eq. 4 of Yin et al (2009)
         Jmax(:)=vj2(:,jl)
         JJ(:,jl) = ( alpha_LL(jv) * Iabs(:) + Jmax(:) - sqrt((alpha_LL(jv) * Iabs(:) + Jmax(:) )**2. &
              - 4 * theta(jv) * Jmax(:) * alpha_LL(jv) * Iabs(:)) ) &
              / ( 2 * theta(jv))

         !
         IF ( is_c4(jv) )  THEN
            !
            ! @addtogroup Photosynthesis
            ! @{   
            !
            ! 2.4.2 Assimilation for C4 plants (Collatz et al., 1992)\n
            !! \latexonly
            !! \input{diffuco_trans_co2_2.4.2.tex}
            !! \endlatexonly
            ! @}           
            !
            !
            !
            IF (nic .GT. 0) THEN
               DO inic=1,nic

                  ! Analytical resolution of the Assimilation based Yin et al. (2009)
                  icinic=index_calc(inic)

                  ! Eq. 28 of Yin et al. (2009)
                  fcyc= 1. - ( 4.*(1.-fpsir(jv))*(1.+fQ(jv)) + 3.*h_protons(jv)*fpseudo(jv) ) / &
                       ( 3.*h_protons(jv) - 4.*(1.-fpsir(jv)))
                                    
                  ! See paragraph after eq. (20b) of Yin et al.
                  Rm=Rd(icinic,jl)/2.
                                
                  ! We assume that cs_star equals ci_star (see Comment in legend of Fig. 6 of Yin et al. (2009)
                  ! Equation 26 of Yin et al. (2009)
                  Cs_star = (gbs(jv) * low_gamma_star(icinic) * Oi - &
                       ( 1. + low_gamma_star(icinic) * alpha(jv) / 0.047) * Rd(icinic,jl) + Rm ) &
                       / ( gbs(jv) + kp(jv) ) 

                  ! eq. 11 of Yin et al (2009)
                  J2 = JJ(icinic,jl) / ( 1. - fpseudo(jv) / ( 1. - fcyc ) )

                  ! Equation right after eq. (20d) of Yin et al. (2009)
                  z = ( 2. + fQ(jv) - fcyc ) / ( h_protons(jv) * (1. - fcyc ))

                  VpJ2 = fpsir(jv) * J2 * z / 2.

                  A_3=9999.

                  ! See eq. right after eq. 18 of Yin et al. (2009)
                  DO limit_photo=1,2
                     ! Is Vc limiting the Assimilation
                     IF ( limit_photo .EQ. 1 ) THEN
                        a = 1. + kp(jv) / gbs(jv)
                        b = 0.
                        x1 = Vc2(icinic,jl)
                        x2 = KmC(icinic)/KmO(icinic)
                        x3 = KmC(icinic)
                        ! Is J limiting the Assimilation
                     ELSE
                        a = 1.
                        b = VpJ2
                        x1 = (1.- fpsir(jv)) * J2 * z / 3.
                        x2 = 7. * low_gamma_star(icinic) / 3.
                        x3 = 0.
                     ENDIF

                     m=fvpd(icinic)-g0var(icinic)/gb_co2(icinic)
                     d=g0var(icinic)*(Ca(icinic)-Cs_star) + fvpd(icinic)*Rd(icinic,jl)
                     f=(b-Rm-low_gamma_star(icinic)*Oi*gbs(jv))*x1*d + a*gbs(jv)*x1*Ca(icinic)*d
                     j=(b-Rm+gbs(jv)*x3 + x2*gbs(jv)*Oi)*m + (alpha(jv)*x2/0.047-1.)*d &
                          + a*gbs(jv)*(Ca(icinic)*m - d/gb_co2(icinic) - (Ca(icinic) - Cs_star ))
 
                     g=(b-Rm-low_gamma_star(icinic)*Oi*gbs(jv))*x1*m - (alpha(jv)*low_gamma_star(icinic)/0.047+1.)*x1*d &
                          + a*gbs(jv)*x1*(Ca(icinic)*m - d/gb_co2(icinic) - (Ca(icinic)-Cs_star ))
 
                     h=-((alpha(jv)*low_gamma_star(icinic)/0.047+1.)*x1*m + (a*gbs(jv)*x1*(m-1.))/gb_co2(icinic) )
                     i= ( b-Rm + gbs(jv)*x3 + x2*gbs(jv)*Oi )*d + a*gbs(jv)*Ca(icinic)*d
                     l= ( alpha(jv)*x2/0.047 - 1.)*m - (a*gbs(jv)*(m-1.))/gb_co2(icinic)
 
                     p = (j-(h-l*Rd(icinic,jl))) / l
                     q = (i+j*Rd(icinic,jl)-g) / l
                     r = -(f-i*Rd(icinic,jl)) / l 
 
                     ! See Yin et al. (2009) and  Baldocchi (1994)
                     QQ = ( (p**2._r_std) - 3._r_std * q) / 9._r_std
                     UU = ( 2._r_std* (p**3._r_std) - 9._r_std *p*q + 27._r_std *r) /54._r_std

                     IF ( (QQ .GE. 0._r_std) .AND. (ABS(UU/(QQ**1.5_r_std) ) .LE. 1._r_std) ) THEN
                        PSI = ACOS(UU/(QQ**1.5_r_std))
                        A_3_tmp = -2._r_std * SQRT(QQ) * COS(( PSI + 4._r_std * PI)/3._r_std ) - p / 3._r_std
                        IF (( A_3_tmp .LT. A_3 )) THEN
                           A_3 = A_3_tmp
                           info_limitphoto(icinic,jl)=2.
                        ELSE
                        ! In case, J is not limiting the assimilation
                        ! we have to re-initialise a, b, x1, x2 and x3 values
                        ! in agreement with a Vc-limited assimilation 
                           a = 1. + kp(jv) / gbs(jv)
                           b = 0.
                           x1 = Vc2(icinic,jl)
                           x2 = KmC(icinic)/KmO(icinic)
                           x3 = KmC(icinic)
                           info_limitphoto(icinic,jl)=1.
                        ENDIF
                     ENDIF

                     IF ( ( A_3 .EQ. 9999. ) .OR. ( A_3 .LT. (-Rd(icinic,jl)) ) ) THEN
                        IF ( printlev>=4 ) THEN
                           WRITE(numout,*) 'We have a problem in diffuco_trans_co2 for A_3'
                           WRITE(numout,*) 'no real positive solution found for pft:',jv
                           WRITE(numout,*) 't2m:',t2m(icinic)
                           WRITE(numout,*) 'vpd:',vpd(icinic)
                        END IF
                        A_3 = -Rd(icinic,jl)
                     ENDIF
                     assimi(icinic,jl) = A_3

                     IF ( ABS( assimi(icinic,jl) + Rd(icinic,jl) ) .LT. min_sechiba ) THEN
                        gs(icinic,jl) = g0var(icinic)
                        !leaf_ci keeps its initial value (Ca).
                     ELSE
                        ! Eq. 24 of Yin et al. (2009) 
                        Obs = ( alpha(jv) * assimi(icinic,jl) ) / ( 0.047 * gbs(jv) ) + Oi
                        ! Eq. 23 of Yin et al. (2009)
                        Cc(icinic,jl) = ( ( assimi(icinic,jl) + Rd(icinic,jl) ) * ( x2 * Obs + x3 ) + low_gamma_star(icinic) &
                             * Obs * x1 ) &
                             / MAX(min_sechiba, x1 - ( assimi(icinic,jl) + Rd(icinic,jl) ))
                        ! Eq. 22 of Yin et al. (2009)
                        leaf_ci(icinic,jv,jl) = ( Cc(icinic,jl) - ( b - assimi(icinic,jl) - Rm ) / gbs(jv) ) / a
                        ! Eq. 25 of Yin et al. (2009)
                        ! It should be Cs instead of Ca but it seems that 
                        ! other equations in Appendix C make use of Ca
                        gs(icinic,jl) = g0var(icinic) + ( assimi(icinic,jl) + Rd(icinic,jl) ) / &
                             ( Ca(icinic) - Cs_star ) * fvpd(icinic)             
                     ENDIF
                  ENDDO                 ! lim_photo 
               ENDDO
            ENDIF
         ELSE
            !
            ! @addtogroup Photosynthesis
            ! @{   
            !
            ! 2.4.3 Assimilation for C3 plants (Farqhuar et al., 1980)\n
            !! \latexonly
            !! \input{diffuco_trans_co2_2.4.3.tex}
            !! \endlatexonly
            ! @}           
            !
            !
            IF (nic .GT. 0) THEN
               DO inic=1,nic
                  icinic=index_calc(inic)
			
                  A_1=9999.

                  ! See eq. right after eq. 18 of Yin et al. (2009)
                  DO limit_photo=1,2
                     ! Is Vc limiting the Assimilation
                     IF ( limit_photo .EQ. 1 ) THEN
                        x1 = vc2(icinic,jl)
                        ! It should be O not Oi (comment from Vuichard)
                        x2 = KmC(icinic) * ( 1. + 2*gamma_star(icinic)*Sco(icinic) / KmO(icinic) )
                        ! Is J limiting the Assimilation
                     ELSE
                        x1 = JJ(icinic,jl)/4.
                        x2 = 2. * gamma_star(icinic)
                     ENDIF


                     ! See Appendix B of Yin et al. (2009)
                     a = g0var(icinic) * ( x2 + gamma_star(icinic) ) + &
                          ( g0var(icinic) / gm(icinic) + fvpd(icinic) ) * ( x1 - Rd(icinic,jl) )
                     b = Ca(icinic) * ( x1 - Rd(icinic,jl) ) - gamma_star(icinic) * x1 - Rd(icinic,jl) * x2
                     c = Ca(icinic) + x2 + ( 1./gm(icinic) + 1./gb_co2(icinic) ) * ( x1 - Rd(icinic,jl) ) 
                     d = x2 + gamma_star(icinic) + ( x1 - Rd(icinic,jl) ) / gm(icinic)
                     m = 1./gm(icinic) + ( g0var(icinic)/gm(icinic) + fvpd(icinic) ) * ( 1./gm(icinic) + 1./gb_co2(icinic) )  
   
                     p = - ( d + (x1 - Rd(icinic,jl) ) / gm(icinic) + a * ( 1./gm(icinic) + 1./gb_co2(icinic) ) + &
                          ( g0var(icinic)/gm(icinic) + fvpd(icinic) ) * c ) / m
   
                     q = ( d * ( x1 - Rd(icinic,jl) ) + a*c + ( g0var(icinic)/gm(icinic) + fvpd(icinic) ) * b ) / m
                     r = - a * b / m
   
                     ! See Yin et al. (2009) 
                     QQ = ( (p**2._r_std) - 3._r_std * q) / 9._r_std
                     UU = ( 2._r_std* (p**3._r_std) - 9._r_std *p*q + 27._r_std *r) /54._r_std
               
                     IF ( (QQ .GE. 0._r_std) .AND. (ABS(UU/(QQ**1.5_r_std) ) .LE. 1._r_std) ) THEN
                        PSI = ACOS(UU/(QQ**1.5_r_std))
                        A_1_tmp = -2._r_std * SQRT(QQ) * COS( PSI / 3._r_std ) - p / 3._r_std
                        IF (( A_1_tmp .LT. A_1 )) THEN
                           A_1 = A_1_tmp
                           info_limitphoto(icinic,jl)=2.
                        ELSE
                        ! In case, J is not limiting the assimilation
                        ! we have to re-initialise x1 and x2 values
                        ! in agreement with a Vc-limited assimilation 
                           x1 = vc2(icinic,jl)
                           ! It should be O not Oi (comment from Vuichard)
                           x2 = KmC(icinic) * ( 1. + 2*gamma_star(icinic)*Sco(icinic) / KmO(icinic) )                           
                           info_limitphoto(icinic,jl)=1.
                        ENDIF
                     ENDIF
                  ENDDO
                  IF ( (A_1 .EQ. 9999.) .OR. ( A_1 .LT. (-Rd(icinic,jl)) ) ) THEN
                     IF ( printlev>=4 ) THEN
                        WRITE(numout,*) 'We have a problem in diffuco_trans_co2 for A_1'
                        WRITE(numout,*) 'no real positive solution found for pft:',jv
                        WRITE(numout,*) 't2m:',t2m(icinic)
                        WRITE(numout,*) 'vpd:',vpd(icinic)
                     END IF
                     A_1 = -Rd(icinic,jl)
                  ENDIF
                  assimi(icinic,jl) = A_1

                  IF ( ABS( assimi(icinic,jl) + Rd(icinic,jl) ) .LT. min_sechiba ) THEN
                     gs(icinic,jl) = g0var(icinic)
                  ELSE
                     ! Eq. 18 of Yin et al. (2009)
                     Cc(icinic,jl) = ( gamma_star(icinic) * x1 + ( assimi(icinic,jl) + Rd(icinic,jl) ) * x2 )  &
                          / MAX( min_sechiba, x1 - ( assimi(icinic,jl) + Rd(icinic,jl) ) )
                     ! Eq. 17 of Yin et al. (2009)
                     leaf_ci(icinic,jv,jl) = Cc(icinic,jl) + assimi(icinic,jl) / gm(icinic) 
                     ! See eq. right after eq. 15 of Yin et al. (2009)
                     ci_star = gamma_star(icinic) - Rd(icinic,jl) / gm(icinic)
                     ! 
                     ! Eq. 15 of Yin et al. (2009)
                     gs(icinic,jl) = g0var(icinic) + ( assimi(icinic,jl) + Rd(icinic,jl) ) / ( leaf_ci(icinic,jv,jl) &
                          - ci_star ) * fvpd(icinic)
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
         !
         IF (nic .GT. 0) THEN
            !
            DO inic=1,nic
               !
               ! @addtogroup Photosynthesis
               ! @{   
               !
               !! 2.4.4 Estimatation of the stomatal conductance (Ball et al., 1987).\n
               !! \latexonly
               !! \input{diffuco_trans_co2_2.4.4.tex}
               !! \endlatexonly
               ! @}           
               !
               icinic=index_calc(inic)
               !
               ! keep stomatal conductance of topmost level
               !
               IF ( jl .EQ. 1 ) THEN
                  leaf_gs_top(icinic) = gs(icinic,jl)
                  !
               ENDIF
               !
               ! @addtogroup Photosynthesis
               ! @{   
               !
               !! 2.4.5 Integration at the canopy level\n
               !! \latexonly
               !! \input{diffuco_trans_co2_2.4.5.tex}
               !! \endlatexonly
               ! @}           
               ! total assimilation and conductance
               assimtot(icinic) = assimtot(icinic) + &
                    assimi(icinic,jl) * (laitab(jl+1)-laitab(jl))
               Rdtot(icinic) = Rdtot(icinic) + &
                    Rd(icinic,jl) * (laitab(jl+1)-laitab(jl))
               gstot(icinic) = gstot(icinic) + &
                    gs(icinic,jl) * (laitab(jl+1)-laitab(jl))
               !
               ilai(icinic) = jl
               !
            ENDDO
            !
         ENDIF
      ENDDO  ! loop over LAI steps
 
      IF(jv==testpft) THEN
         templeafci(:,:)=leaf_ci(:,testpft,:)
         CALL histwrite_p(hist_id, 'Cc', kjit, Cc, kjpindex*(nlai+1), indexlai)
         CALL histwrite_p(hist_id, 'Vc', kjit, Vc2, kjpindex*(nlai+1), indexlai)
         CALL histwrite_p(hist_id, 'Vj', kjit, JJ, kjpindex*(nlai+1), indexlai)
         CALL histwrite_p(hist_id, 'limitphoto', kjit, info_limitphoto, kjpindex*(nlai+1), indexlai)
         CALL histwrite_p(hist_id, 'gammastar', kjit, gamma_star, kjpindex,index)
         CALL histwrite_p(hist_id, 'Kmo', kjit, Kmo, kjpindex,index)
         CALL histwrite_p(hist_id, 'Kmc', kjit, Kmc, kjpindex,index)
         CALL histwrite_p(hist_id, 'gm', kjit, gm, kjpindex, index)
         CALL histwrite_p(hist_id, 'gs', kjit, gs, kjpindex*(nlai+1), indexlai)
         CALL histwrite_p(hist_id, 'leafci', kjit, templeafci, kjpindex*(nlai), indexlai0)
         CALL histwrite_p(hist_id, 'assimi', kjit, assimi, kjpindex*(nlai+1), indexlai)
         CALL histwrite_p(hist_id, 'Rd', kjit, Rd, kjpindex*(nlai+1), indexlai)
      ENDIF
      !! Calculated intercellular CO2 over nlai needed for the chemistry module
      cim(:,jv)=0.
      laisum(:)=0
      DO jl=1,nlai
         WHERE (laitab(jl) .LE. lai(:,jv) )
            cim(:,jv)= cim(:,jv)+leaf_ci(:,jv,jl)*(laitab(jl+1)-laitab(jl))
            laisum(:)=laisum(:)+ (laitab(jl+1)-laitab(jl))
         ENDWHERE
      ENDDO
      WHERE (laisum(:)>0)
         cim(:,jv)= cim(:,jv)/laisum(:)
      ENDWHERE


      !
      !! 2.5 Calculate resistances
      !
      IF (nia .GT. 0) THEN
         !
         DO inia=1,nia
            !
            iainia=index_assi(inia)

            !! Mean stomatal conductance for CO2 (mol m-2 s-1)
            gsmean(iainia,jv) = gstot(iainia)
            !
            ! cimean is the "mean ci" calculated in such a way that assimilation 
            ! calculated in enerbil is equivalent to assimtot
            !
            IF ( ABS(gsmean(iainia,jv)-g0var(iainia)*laisum(iainia)) .GT. min_sechiba) THEN
               cimean(iainia,jv) = (fvpd(iainia)*(assimtot(iainia)+Rdtot(iainia))) /&
                 (gsmean(iainia,jv)-g0var(iainia)*laisum(iainia)) + gamma_star(iainia) 
            ELSE
               cimean(iainia,jv) = gamma_star(iainia) 
            ENDIF
                 
            ! conversion from umol m-2 (PFT) s-1 to gC m-2 (mesh area) tstep-1
            gpp(iainia,jv) = assimtot(iainia)*12e-6*veget_max(iainia,jv)*dt_sechiba
            
            !
            ! conversion from mol/m^2/s to m/s
            !
            ! As in Pearcy, Schulze and Zimmermann
            ! Measurement of transpiration and leaf conductance
            ! Chapter 8 of Plant Physiological Ecology
            ! Field methods and instrumentation, 1991
            ! Editors:
            !
            !    Robert W. Pearcy,
            !    James R. Ehleringer,
            !    Harold A. Mooney,
            !    Philip W. Rundel
            !
            ! ISBN: 978-0-412-40730-7 (Print) 978-94-010-9013-1 (Online)

            gstot(iainia) =  mol_to_m_1 *(t2m(iainia)/tp_00)*&
                 (pb_std/pb(iainia))*gstot(iainia)*ratio_H2O_to_CO2
            gstop(iainia) =  mol_to_m_1 * (t2m(iainia)/tp_00)*&
                 (pb_std/pb(iainia))*leaf_gs_top(iainia)*ratio_H2O_to_CO2*&
                 laitab(ilai(iainia)+1)
            !
            rveget(iainia,jv) = un/gstop(iainia)

            !
            !
            ! rstruct is the difference between rtot (=1./gstot) and rveget
            !
            ! Correction Nathalie - le 27 Mars 2006 - Interdire a rstruct d'etre negatif
            !rstruct(iainia,jv) = un/gstot(iainia) - &
            !     rveget(iainia,jv)
            rstruct(iainia,jv) = MAX( un/gstot(iainia) - &
                 rveget(iainia,jv), min_sechiba)
            !
            !
            !! wind is a global variable of the diffuco module.
            speed = MAX(min_wind, wind(iainia))
            !
            ! beta for transpiration
            !
            ! Corrections Nathalie - 28 March 2006 - on advices of Fred Hourdin
            !! Introduction of a potentiometer rveg_pft to settle the rveg+rstruct sum problem in the coupled mode.
            !! rveg_pft=1 in the offline mode. rveg_pft is a global variable declared in the diffuco module.
            !vbeta3(iainia,jv) = veget_max(iainia,jv) * &
            !  (un - zqsvegrap(iainia)) * &
            !  (un / (un + speed * q_cdrag(iainia) * (rveget(iainia,jv) + &
            !   rstruct(iainia,jv))))
            !! Global resistance of the canopy to evaporation
            IF (.NOT. ok_LAIdev(jv)) THEN
                cresist=(un / (un + speed * q_cdrag(iainia) * &
                     veget(iainia,jv)/veget_max(iainia,jv) * &
                     (rveg_pft(jv)*(rveget(iainia,jv) + rstruct(iainia,jv)))))
            ELSE  ! croplands, xuhui
                ! using pft specific drag coefficient 
                cresist=(un / (un + speed * q_cdrag_pft(iainia, jv) * &
                     veget(iainia,jv)/veget_max(iainia,jv) * &
                     (rveg_pft(jv)*(rveget(iainia,jv) + rstruct(iainia,jv)))))
            ENDIF

            IF ( humrel(iainia,jv) >= min_sechiba ) THEN
               vbeta3(iainia,jv) = veget(iainia,jv) * &
                    (un - zqsvegrap(iainia)) * cresist + &
                    MIN( vbeta23(iainia,jv), veget(iainia,jv) * &
                    zqsvegrap(iainia) * cresist )
            ELSE
               ! Because of a minimum conductance g0, vbeta3 cannot be zero even if humrel=0
               ! in the above equation.
               ! Here, we force transpiration to be zero when the soil cannot deliver it
               vbeta3(iainia,jv) = zero
            END IF

            ! vbeta3pot for computation of potential transpiration (needed for irrigation)
            vbeta3pot(iainia,jv) = MAX(zero, veget(iainia,jv) * cresist)
            !
            !
         ENDDO
         !
      ENDIF
      !
   END DO         ! loop over vegetation types
   !
      
   IF (printlev>=3) WRITE (numout,*) ' diffuco_trans_co2 done '


END SUBROUTINE diffuco_trans_co2


!! ================================================================================================================================
!! SUBROUTINE	   : diffuco_comb
!!
!>\BRIEF           This routine combines the previous partial beta 
!! coefficients and calculates the total alpha and complete beta coefficients.
!!
!! DESCRIPTION	   : Those integrated coefficients are used to calculate (in enerbil.f90) the total evapotranspiration 
!! from the grid-cell. \n
!!
!! In the case that air is more humid than surface, dew deposition can occur (negative latent heat flux). 
!! In this instance, for temperature above zero, all of the beta coefficients are set to 0, except for 
!! interception (vbeta2) and bare soil (vbeta4 with zero soil resistance). The amount of water that is 
!! intercepted by leaves is calculated based on the value of LAI of the surface. In the case of freezing 
!! temperatures, water is added to the snow reservoir, and so vbeta4 and vbeta2 are set to 0, and the 
!! total vbeta is set to 1.\n
!!
!! \latexonly 
!!     \input{diffucocomb1.tex}
!! \endlatexonly
!!
!! The beta and alpha coefficients are initially set to 1.
!! \latexonly 
!!     \input{diffucocomb2.tex}
!! \endlatexonly
!!
!! If snow is lower than the critical value:
!! \latexonly 
!!     \input{diffucocomb3.tex}
!! \endlatexonly
!! If in the presence of dew:
!! \latexonly 
!!     \input{diffucocomb4.tex}
!! \endlatexonly
!!
!! Determine where the water goes (soil, vegetation, or snow)
!! when air moisture exceeds saturation.
!! \latexonly 
!!     \input{diffucocomb5.tex}
!! \endlatexonly
!!
!! If it is not freezing dew is put into the interception reservoir and onto the bare soil. If it is freezing, 
!! water is put into the snow reservoir. 
!! Now modify vbetas where necessary: for soil and snow
!! \latexonly 
!!     \input{diffucocomb6.tex}
!! \endlatexonly
!!
!! and for vegetation
!! \latexonly 
!!     \input{diffucocomb7.tex}
!! \endlatexonly
!!
!! Then compute part of dew that can be intercepted by leafs.
!!
!! There will be no transpiration when air moisture is too high, under any circumstance
!! \latexonly 
!!     \input{diffucocomb8.tex}
!! \endlatexonly
!!
!! There will also be no interception loss on bare soil, under any circumstance.
!! \latexonly 
!!     \input{diffucocomb9.tex}
!! \endlatexonly
!!
!! The flowchart details the 'decision tree' which underlies the module. 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): vbeta1, vbeta4, humrel, vbeta2, vbeta3, vbeta
!!
!! REFERENCE(S) :
!! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
!! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
!! Circulation Model. Journal of Climate, 6, pp.248-273
!! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
!! sur le cycle de l'eau, PhD Thesis, available from:
!! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
!!
!! FLOWCHART    :
!! \latexonly 
!!     \includegraphics[scale=0.25]{diffuco_comb_flowchart.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE diffuco_comb (kjpindex, humrel, rau, u, v, q_cdrag, pb, qair, temp_sol, temp_air, &
       & snow, veget, veget_max, lai, tot_bare_soil, vbeta1, vbeta2, vbeta3 , vbeta4, vbeta4_pft, vbeta, vbeta_pft, qsintmax)    
    
    ! Ajout qsintmax dans les arguments de la routine Nathalie / le 13-03-2006

  !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables
    
    INTEGER(i_std), INTENT(in)                           :: kjpindex   !! Domain size (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: rau        !! Air Density (kg m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: u          !! Eastward Lowest level wind speed (m s^{-1}) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: v          !! Nortward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: q_cdrag    !! Surface drag coefficient  (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: pb         !! Lowest level pressure (hPa)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: qair       !! Lowest level specific air humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: temp_sol   !! Skin temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: temp_air   !! Lower air temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: snow       !! Snow mass (kg)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: veget      !! Fraction of vegetation type (fraction)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: veget_max  !! Maximum fraction of vegetation type (fraction)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: lai        !! Leaf area index (m^2 m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: qsintmax   !! Maximum water on vegetation (kg m^{-2})
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)        :: tot_bare_soil!! Total evaporating bare soil fraction 

    !! 0.2 Output variables
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)       :: vbeta      !! Total beta coefficient (-)
    REAL(r_std),DIMENSION (kjpindex, nvm), INTENT (out)  :: vbeta_pft    !! Total beta coefficient (-)

    !! 0.3 Modified variables 
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)     :: vbeta1     !! Beta for sublimation (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)     :: vbeta4     !! Beta for Bare soil evaporation (-) 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)  :: vbeta4_pft !! Beta for Bare soil evaporation (-) 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout) :: humrel     !! Soil moisture stress (within range 0 to 1)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout) :: vbeta2     !! Beta for interception loss (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout) :: vbeta3     !! Beta for Transpiration (-)
    
    !! 0.4 Local variables
    
    INTEGER(i_std)                                       :: ji, jv
    REAL(r_std)                                          :: zevtest, zsoil_moist, zrapp
    REAL(r_std), DIMENSION(kjpindex)                     :: qsatt
    LOGICAL, DIMENSION(kjpindex)                         :: toveg, tosnow
    REAL(r_std)                                          :: coeff_dew_veg
!_ ================================================================================================================================

    vbeta_pft = zero
    
    !! 1 If we are in presence of dew
     
    CALL qsatcalc (kjpindex, temp_sol, pb, qsatt)

    
    !! 1.1 Determine where the water goes 
    !! Determine where the water goes (soil, vegetation, or snow)
    !! when air moisture exceeds saturation.
    !! \latexonly 
    !!     \input{diffucocomb5.tex}
    !! \endlatexonly
    toveg(:) = .FALSE.
    tosnow(:) = .FALSE.
    DO ji = 1, kjpindex
      IF ( qsatt(ji) .LT. qair(ji) ) THEN
          IF (temp_air(ji) .GT. tp_00) THEN
              !! If it is not freezing dew is put into the 
              !! interception reservoir and onto the bare soil.
              toveg(ji) = .TRUE.
          ELSE
              !! If it is freezing water is put into the 
              !! snow reservoir.
              tosnow(ji) = .TRUE.
          ENDIF
      ENDIF
    END DO

    !! 1.2 Now modify vbetas where necessary.
    
    !! 1.2.1 Soil and snow
    !! \latexonly 
    !!     \input{diffucocomb6.tex}
    !! \endlatexonly
    DO ji = 1, kjpindex
      IF ( toveg(ji) ) THEN
        vbeta1(ji) = zero
        vbeta4(ji) = tot_bare_soil(ji)
        ! Correction Nathalie - le 13-03-2006: le vbeta ne sera calcule qu'une fois tous les vbeta2 redefinis
        !vbeta(ji) = vegetsum(ji)
        vbeta(ji) = vbeta4(ji)
        DO jv = 1,nvm
            IF (tot_bare_soil(ji)>0) THEN
                vbeta_pft(ji,jv) = vbeta4_pft(ji,jv)
                !vbeta_pft(ji,jv) = vbeta4(ji)*(veget_max(ji,jv) - veget(ji,jv))/tot_bare_soil(ji)
            ELSE 
                vbeta_pft(ji,jv) = 0   
            ENDIF
        ENDDO
      ENDIF
      IF ( tosnow(ji) ) THEN
        vbeta1(ji) = un
        vbeta4(ji) = zero
        vbeta(ji) = un
        DO jv = 1,nvm
            vbeta_pft(ji,jv) = un * veget_max(ji,jv)
        ENDDO
      ENDIF
    ENDDO

    !! 1.2.2 Vegetation and interception loss
    !! \latexonly 
    !!     \input{diffucocomb7.tex}
    !! \endlatexonly
    DO jv = 1, nvm
      
      DO ji = 1, kjpindex
               
        IF ( toveg(ji) ) THEN
           IF (qsintmax(ji,jv) .GT. min_sechiba) THEN
              
              ! Compute part of dew that can be intercepted by leafs.
              IF ( lai(ji,jv) .GT. min_sechiba) THEN
                IF (lai(ji,jv) .GT. 1.5) THEN
                   coeff_dew_veg= &
                         &   dew_veg_poly_coeff(6)*lai(ji,jv)**5 &
                         & - dew_veg_poly_coeff(5)*lai(ji,jv)**4 &
                         & + dew_veg_poly_coeff(4)*lai(ji,jv)**3 &
                         & - dew_veg_poly_coeff(3)*lai(ji,jv)**2 &
                         & + dew_veg_poly_coeff(2)*lai(ji,jv) &
                         & + dew_veg_poly_coeff(1)
                 ELSE
                    coeff_dew_veg=un
                 ENDIF
              ELSE
                 coeff_dew_veg=zero
              ENDIF
              IF (jv .EQ. 1) THEN
                 ! This line may not work with CWRR when frac_bare is distributed among three soiltiles
                 ! Fortunately, qsintmax(ji,1)=0 (LAI=0 in PFT1) so we never pass here
                 vbeta2(ji,jv) = coeff_dew_veg*tot_bare_soil(ji)
              ELSE
                 vbeta2(ji,jv) = coeff_dew_veg*veget(ji,jv)
              ENDIF
           ELSE
              vbeta2(ji,jv) = zero ! if qsintmax=0, vbeta2=0
           ENDIF
           vbeta_pft(ji,jv) = vbeta_pft(ji,jv) + vbeta2(ji,jv)
        ENDIF
        IF ( tosnow(ji) ) vbeta2(ji,jv) = zero
        
      ENDDO
      
    ENDDO

    !! 1.2.3 Vegetation and transpiration  
    !! There will be no transpiration when air moisture is too high, under any circumstance
    !! \latexonly 
    !!     \input{diffucocomb8.tex}
    !! \endlatexonly
    DO jv = 1, nvm
      DO ji = 1, kjpindex
        IF ( qsatt(ji) .LT. qair(ji) ) THEN
          vbeta3(ji,jv) = zero
          humrel(ji,jv) = zero
        ENDIF
      ENDDO
    ENDDO
    
    
    !! 1.2.4 Overrules 1.2.2
    !! There will also be no interception loss on bare soil, under any circumstance.
    !! \latexonly 
    !!     \input{diffucocomb9.tex}
    !! \endlatexonly
    DO ji = 1, kjpindex
       IF ( qsatt(ji) .LT. qair(ji) ) THEN
          vbeta2(ji,1) = zero
       ENDIF
    ENDDO

    !! 2  Now calculate vbeta in all cases (the equality needs to hold for enerbil to be consistent)

    DO ji = 1, kjpindex
          vbeta(ji) = vbeta4(ji) + SUM(vbeta2(ji,:)) + SUM(vbeta3(ji,:))

          IF (vbeta(ji) .LT. min_sechiba) THEN
             vbeta(ji) = zero
             vbeta4(ji) = zero
             vbeta2(ji,:)= zero
             vbeta3(ji,:)= zero
          END IF
    ENDDO 

    IF (printlev>=3) WRITE (numout,*) ' diffuco_comb done '

  END SUBROUTINE diffuco_comb


!! ================================================================================================================================
!! SUBROUTINE	: diffuco_raerod
!!
!>\BRIEF	Computes the aerodynamic resistance, for cases in which the
!! surface drag coefficient is provided by the coupled atmospheric model LMDZ and  when the flag
!! 'ldq_cdrag_from_gcm' is set to TRUE
!!
!! DESCRIPTION	: Simply computes the aerodynamic resistance, for cases in which the
!! surface drag coefficient is provided by the coupled atmospheric model LMDZ. If the surface drag coefficient
!! is not provided by the LMDZ or signalled by the flag 'ldq_cdrag_from_gcm' set to FALSE, then the subroutine
!! diffuco_aero is called instead of this one.
!!
!! Calculation of the aerodynamic resistance, for diganostic purposes. First calculate wind speed:
!! \latexonly 
!!     \input{diffucoaerod1.tex}
!! \endlatexonly       
!!
!! next calculate ::raero
!! \latexonly 
!!     \input{diffucoaerod2.tex}
!! \endlatexonly
!! 
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::raero
!!
!! REFERENCE(S)	:
!! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
!! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
!! Circulation Model. Journal of Climate, 6, pp.248-273
!! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influence de l'irrigation
!! sur le cycle de l'eau, PhD Thesis, available from:
!! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
!!
!! FLOWCHART    :  None
!! \n
!_ ================================================================================================================================

  SUBROUTINE diffuco_raerod (kjpindex, u, v, q_cdrag, raero)
    
    IMPLICIT NONE
    
  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                     :: kjpindex     !! Domain size (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)  :: u            !! Eastward Lowest level wind velocity (m s^{-1}) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)  :: v            !! Northward Lowest level wind velocity (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)  :: q_cdrag      !! Surface drag coefficient  (-)
    
    !! 0.2 Output variables 
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (out) :: raero        !! Aerodynamic resistance (s m^{-1})
     
    !! 0.3 Modified variables

    !! 0.4 Local variables
    
    INTEGER(i_std)                                 :: ji           !! (-)
    REAL(r_std)                                    :: speed        !! (m s^{-1})
!_ ================================================================================================================================
   
  !! 1. Simple calculation of the aerodynamic resistance, for diganostic purposes.

    DO ji=1,kjpindex

       !! \latexonly 
       !!     \input{diffucoaerod1.tex}
       !! \endlatexonly       
       speed = MAX(min_wind, wind(ji))

       !! \latexonly 
       !!     \input{diffucoaerod2.tex}
       !! \endlatexonly
       raero(ji) = un / (q_cdrag(ji)*speed)
       
    ENDDO
  
  END SUBROUTINE diffuco_raerod


  FUNCTION Arrhenius (kjpindex,temp,ref_temp,energy_act) RESULT ( val_arrhenius )
    !! 0.1 Input variables

    INTEGER(i_std),INTENT(in)                     :: kjpindex          !! Domain size (-)
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)    :: temp              !! Temperature (K)
    REAL(r_std), INTENT(in)                       :: ref_temp          !! Temperature of reference (K)
    REAL(r_std),INTENT(in)                        :: energy_act        !! Activation Energy (J mol-1)
    
    !! 0.2 Result

    REAL(r_std), DIMENSION(kjpindex)              :: val_arrhenius     !! Temperature dependance based on
                                                                       !! a Arrhenius function (-)
    
    val_arrhenius(:)=EXP(((temp(:)-ref_temp)*energy_act)/(ref_temp*RR*(temp(:))))
  END FUNCTION Arrhenius

  FUNCTION Arrhenius_modified_1d (kjpindex,temp,ref_temp,energy_act,energy_deact,entropy) RESULT ( val_arrhenius )
    !! 0.1 Input variables

    INTEGER(i_std),INTENT(in)                     :: kjpindex          !! Domain size (-)
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)    :: temp              !! Temperature (K)
    REAL(r_std), INTENT(in)                       :: ref_temp          !! Temperature of reference (K)
    REAL(r_std),INTENT(in)                        :: energy_act        !! Activation Energy (J mol-1)
    REAL(r_std),INTENT(in)                        :: energy_deact      !! Deactivation Energy (J mol-1)
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)    :: entropy           !! Entropy term (J K-1 mol-1)
        
    !! 0.2 Result

    REAL(r_std), DIMENSION(kjpindex)              :: val_arrhenius     !! Temperature dependance based on
                                                                       !! a Arrhenius function (-)
    
    val_arrhenius(:)=EXP(((temp(:)-ref_temp)*energy_act)/(ref_temp*RR*(temp(:))))  &
         * (1. + EXP( (ref_temp * entropy(:) - energy_deact) / (ref_temp * RR ))) &
         / (1. + EXP( (temp(:) * entropy(:) - energy_deact) / ( RR*temp(:))))
         
  END FUNCTION Arrhenius_modified_1d

  FUNCTION Arrhenius_modified_0d (kjpindex,temp,ref_temp,energy_act,energy_deact,entropy) RESULT ( val_arrhenius )
    !! 0.1 Input variables

    INTEGER(i_std),INTENT(in)                     :: kjpindex          !! Domain size (-)
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)    :: temp              !! Temperature (K)
    REAL(r_std), INTENT(in)                       :: ref_temp          !! Temperature of reference (K)
    REAL(r_std),INTENT(in)                        :: energy_act        !! Activation Energy (J mol-1)
    REAL(r_std),INTENT(in)                        :: energy_deact      !! Deactivation Energy (J mol-1)
    REAL(r_std),INTENT(in)                        :: entropy           !! Entropy term (J K-1 mol-1)
        
    !! 0.2 Result

    REAL(r_std), DIMENSION(kjpindex)              :: val_arrhenius     !! Temperature dependance based on
                                                                       !! a Arrhenius function (-)
    
    val_arrhenius(:)=EXP(((temp(:)-ref_temp)*energy_act)/(ref_temp*RR*(temp(:))))  &
         * (1. + EXP( (ref_temp * entropy - energy_deact) / (ref_temp * RR ))) &
         / (1. + EXP( (temp(:) * entropy - energy_deact) / ( RR*temp(:))))
         
  END FUNCTION Arrhenius_modified_0d


END MODULE diffuco
