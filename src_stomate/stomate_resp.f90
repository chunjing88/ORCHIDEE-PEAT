! =================================================================================================================================
! MODULE           : stomate_resp
!
! CONTACT          : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE          : IPSL (2006)
!                  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF           Calculates maintenance respiration for different plant components
!!
!!\n DESCRIPTION   : None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	   :
!!- McCree KJ. An equation for the respiration of white clover plants grown under controlled conditions. 
!! In: Setlik I, editor. Prediction and measurement of photosynthetic productivity. Wageningen, The Netherlands: Pudoc; 1970. p. 221-229.
!! - Krinner G, Viovy N, de Noblet-Ducoudre N, Ogee J, Polcher J, Friedlingstein P,
!! Ciais P, Sitch S, Prentice I C (2005) A dynamic global vegetation model for studies
!! of the coupled atmosphere-biosphere system. Global Biogeochemical Cycles, 19, GB1015,
!! doi: 10.1029/2003GB002199.\n
!! Ruimy A., Dedieu G., Saugier B. (1996), TURC: A diagnostic model
!! of continental gross primary productivity and net primary productivity,
!! Global Biogeochemical Cycles, 10, 269-285.\n

!! SVN :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_stomate/stomate_resp.f90 $
!! $Date: 2018-02-07 12:28:38 +0100 (Wed, 07 Feb 2018) $
!! $Revision: 4956 $
!! \n
!_ ================================================================================================================================
 
MODULE stomate_resp

  ! modules used:
  USE stomate_data
  USE pft_parameters
  USE constantes  
  USE constantes_soil 

  IMPLICIT NONE

  ! private & public routines
  PRIVATE
  PUBLIC maint_respiration,maint_respiration_clear

  LOGICAL, SAVE                                              :: firstcall_resp = .TRUE.                 !! first call
!$OMP THREADPRIVATE(firstcall_resp)

CONTAINS


!! ================================================================================================================================
!! SUBROUTINE 	: maint_respiration_clear
!!
!>\BRIEF        : Set the flag ::firstcall_resp to .TRUE. and as such activate section 
!!                1.1 of the subroutine maint_respiration (see below).
!_ ================================================================================================================================

  SUBROUTINE maint_respiration_clear
    firstcall_resp=.TRUE.
  END SUBROUTINE maint_respiration_clear


!! ================================================================================================================================
!! SUBROUTINE 	: maint_respiration
!!
!>\BRIEF         Calculate PFT maintenance respiration of each living plant part by 
!! multiplying the biomass of plant part by maintenance respiration coefficient which
!! depends on long term mean annual temperature. PFT maintenance respiration is carbon flux 
!! with the units @tex $(gC.m^{-2}dt_sechiba^{-1})$ @endtex, and the convention is from plants to the 
!! atmosphere.
!!
!! DESCRIPTION : The maintenance respiration of each plant part for each PFT is the biomass of the plant 
!! part multiplied by maintenance respiration coefficient. The biomass allocation to different 
!! plant parts is done in routine stomate_alloc.f90. The maintenance respiration coefficient is 
!! calculated in this routine.\n
!!
!! The maintenance respiration coefficient is the fraction of biomass that is lost during 
!! each time step, which increases linearly with temperature (2-meter air temperature for aboveground plant
!! tissues; root-zone temperature for below-ground tissues). Air temperature is an input forcing variable. 
!! Root-zone temperature is a convolution of root and soil temperature profiles and also calculated 
!! in this routine.\n
!!
!! The calculation of maintenance respiration coefficient (fraction of biomass respired) depends linearly
!! on temperature:
!! - the relevant temperature for different plant parts (air temperature or root-zone temperature)\n
!! - intercept: prescribed maintenance respiration coefficients at 0 Degree Celsius for 
!!   different plant parts for each PFT in routine stomate_constants.f90\n
!! - slope: calculated with a quadratic polynomial with the multi-annual mean air temperature 
!! (the constants are in routine stomate_constants.f90) as follows\n 
!!    \latexonly
!!      \input{resp3.tex} 
!!    \endlatexonly
!!   Where, maint_resp_slope1, maint_resp_slope2, maint_resp_slope3 are constant in stomate_constants.f90.
!!   Then coeff_maint is calculated as follows:\n
!!    \latexonly
!!      \input{resp4.tex} 
!!    \endlatexonly  
!! If the calculation result is negative, maintenance respiration coefficient will take the value 0.
!! Therefore the maintenance respiration will also be 0.\n
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): PFT maintenance respiration of different plant parts (::resp_maint_part_radia)
!!
!! REFERENCE(S)	:
!! McCree KJ. An equation for the respiration of white clover plants grown under controlled conditions. In: 
!! Setlik I, editor. Prediction and measurement of photosynthetic productivity. Wageningen, 
!! The Netherlands: Pudoc; 1970. p. 221-229.
!! Krinner G, Viovy N, de Noblet-Ducoudre N, Ogee J, Polcher J, Friedlingstein P,
!! Ciais P, Sitch S, Prentice I C (2005) A dynamic global vegetation model for studies
!! of the coupled atmosphere-biosphere system. Global Biogeochemical Cycles, 19, GB1015,
!! doi: 10.1029/2003GB002199.\n
!! Ruimy A., Dedieu G., Saugier B. (1996), TURC: A diagnostic model
!! of continental gross primary productivity and net primary productivity,
!! Global Biogeochemical Cycles, 10, 269-285.\n
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE maint_respiration ( npts,lai, t2m,t2m_longterm,stempdiag,height,&
       rprof,biomass,resp_maint_part_radia, &
!gmjc
       sla_calc, humrel_month) !!!qcj++ peatland, following Druel et al. 2017 GMD paper, add humrel_month for dessication impact
!end gmjc
!! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                         :: npts      !! Domain size - number of grid cells (unitless)
    REAL(r_std), DIMENSION(npts), INTENT(in)           :: t2m       !! 2 meter air temperature - forcing variable (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)           :: t2m_longterm !! Long term annual mean 2 meter reference air temperatures 
                                                                       !! calculated in stomate_season.f90 (K)
    REAL(r_std), DIMENSION(npts,nslm), INTENT (in)     :: stempdiag !! Soil temperature of each soil layer (K)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)       :: height    !! height of vegetation (m)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)       :: rprof     !! PFT root depth as calculated in stomate.f90 from parameter 
                                                                    !! humcste which is root profile for different PFTs 
                                                                    !! in slowproc.f90 (m)
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements),INTENT(in) :: biomass   !! PFT total biomass calculated in stomate_alloc.f90 
                                                                    !! @tex $(gC.m^{-2})$ @endtex

    !! 0.2 Output variables

!!!qcj++ peatland
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)         :: humrel_month

    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)        :: lai                   !! PFT leaf area index @tex $(m^2 m^{-2})$ @endtex

    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(out) :: resp_maint_part_radia !! PFT maintenance respiration of different plant
                                                                                  !! parts @tex $(gC.m^{-2}dt_sechiba^{-1} )$ @endtex
!gmjc
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: sla_calc
!end gmjc
    !! 0.3 Modified variables
 
    !! 0.4 Local variables

    REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)    :: z_soil       !! Variable to store depth of the different soil layers (m)
!$OMP THREADPRIVATE(z_soil)
    REAL(r_std), DIMENSION(npts,nvm)        :: t_root               !! PFT root temperature (convolution of root and soil 
                                                                    !! temperature profiles) (K)
    REAL(r_std), DIMENSION(npts,nvm,nparts) :: coeff_maint          !! PFT maintenance respiration coefficients of different 
                                                                    !! plant compartments at 0 deg C 
                                                                    !! @tex $(g.g^{-1}dt_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts)            :: rpc                  !! Scaling factor for integrating vertical soil 
                                                                    !! profiles (unitless)
    REAL(r_std), DIMENSION(npts,nparts)     :: t_maint_radia        !! Temperature which is pertinent for maintenance respiration, 
                                                                    !! which is air/root temperature for above/below-ground 
                                                                    !! compartments (K) 
    REAL(r_std), DIMENSION(npts)            :: tl                   !! Long term reference temperature in degrees Celcius 
                                                                    !! (= t2m_longterm - 273.15) (C)
    REAL(r_std), DIMENSION(npts)            :: slope                !! slope of the temperature dependence of maintenance 
                                                                    !! respiration coefficient (1/K)
    INTEGER(i_std)                          :: i,j,k,l,m            !! Indeces (unitless)
    INTEGER(i_std)                          :: ier                  !! Error handling 

!_ ================================================================================================================================
    
    
    IF (printlev>=3) WRITE(numout,*) 'Entering respiration'
    
 !! 1. Initializations
    
    IF ( firstcall_resp ) THEN

       !! 1.1. Soil levels (first call only)
       !       Set the depth of the different soil layers (number of layers: nslm) 
       !       previously calculated as variable diaglev in routines sechiba.f90 and slowproc.f90
       ALLOCATE(z_soil(0:nslm), stat=ier)
       IF ( ier /= 0 ) CALL ipslerr_p(3,'maint_respiration','Pb in allocate of z_soil','','')
       z_soil(0) = zero
       z_soil(1:nslm) = diaglev(1:nslm)

       firstcall_resp = .FALSE.
    ENDIF

    
    
    !! 1.2. Calculate root temperature
    !       Calculate root temperature as the convolution of root and soil temperature profiles
    DO j = 2,nvm ! Loop over # PFTs

       !! 1.2.1 Calculate rpc
       !  - rpc is an integration constant to make the integral over the root profile is equal 'one', 
       !    calculated as follows:\n
       !  \latexonly
       !    \input{resp1.tex} 
       !  \endlatexonly
       rpc(:) = un / ( un - EXP( -z_soil(nslm) / rprof(:,j) ) )

       !! 1.2.2 Calculate root temperature
       !        - Integrate root profile temperature (K) over soil layers (number of layers = nslm)
       !          with rpc and soil temperature (K) of each soil layer as follows:\n
       !        \latexonly
       !          \input{resp2.tex} 
       !        \endlatexonly
       !        Where, stempdiag is diagnostic temperature profile of soil (K)\n
       t_root(:,j) = zero

       DO l = 1, nslm ! Loop over # soil layers

          t_root(:,j) = &
               t_root(:,j) + stempdiag(:,l) * rpc(:) * &
               ( EXP( -z_soil(l-1)/rprof(:,j) ) - EXP( -z_soil(l)/rprof(:,j) ) )

       ENDDO ! Loop over # soil layers

    ENDDO ! Loop over # PFTs

 !! 2. Define maintenance respiration coefficients

    DO j = 2,nvm ! Loop over # PFTs

       !! 2.1 Temperature for maintenanace respiration
       !      Temperature which is used to calculate maintenance respiration for different plant compartments
       !      (above- and belowground)\n
       !      - for aboveground parts, we use 2-meter air temperature, t2m\n
       !      - for belowground parts, we use root temperature calculated in section 1.2 of this subroutine\n
       
       ! 2.1.1 Aboveground biomass
       t_maint_radia(:,ileaf) = t2m(:)
       t_maint_radia(:,isapabove) = t2m(:)
       t_maint_radia(:,ifruit) = t2m(:)

       ! 2.1.2 Belowground biomass
       t_maint_radia(:,isapbelow) = t_root(:,j)
       t_maint_radia(:,iroot) = t_root(:,j)

       !! 2.1.3 Heartwood biomass
       !        Heartwood does does not respire (coeff_maint_zero is set to zero)

       t_maint_radia(:,iheartbelow) = t_root(:,j)
       t_maint_radia(:,iheartabove) = t2m(:)

       !! 2.1.4 Reserve biomass
       !        Use aboveground temperature for trees and belowground temeperature for grasses
       IF ( is_tree(j) ) THEN
          t_maint_radia(:,icarbres) = t2m(:)
       ELSE
          t_maint_radia(:,icarbres) = t_root(:,j)
       ENDIF

       
       !! 2.2 Calculate maintenance respiration coefficients (coeff_maint)
       !      Maintenance respiration is a fraction of biomass defined by the coefficient 
       !      coeff_maint [Mc Cree, 1969]. Coeff_maint is defined through a linear relationship of temperature [Ruimy et al, 1996]
       !      which slope is the coefficient 'slope' and which intercept is 'coeff_maint_zero'.
       !     - Coeff_maint_zero is defined in stomate_data to cm_zero_plantpartname
       !     - Slope is calculated here through a second-degree polynomial [Krinner et al, 2005] 
       !    equation that makes it dependent on the long term temperature (to represent adaptation
       !    of the ecosystem to long term temperature).
       !         \latexonly
       !           \input{resp3.tex} 
       !         \endlatexonly
       !        Where, maint_resp_slope1, maint_resp_slope2, maint_resp_slope3 are constant in stomate_constants.f90.
       !        Then coeff_maint is calculated as follows:\n
       !         \latexonly
       !           \input{resp4.tex} 
       !         \endlatexonly
       ! If the calculation result is negative, coeff_maint will take the value 0.\n	
       tl(:) = t2m_longterm(:) - ZeroCelsius
       slope(:) = maint_resp_slope(j,1) + tl(:) * maint_resp_slope(j,2) + &
            tl(:)*tl(:) * maint_resp_slope(j,3)

       DO k = 1, nparts ! Loop over # plant parts

          coeff_maint(:,j,k) = &
               MAX( (coeff_maint_zero(j,k)*dt_sechiba/one_day) * &
               ( un + slope(:) * (t_maint_radia(:,k)-ZeroCelsius) ), zero )

       ENDDO ! Loop over # plant parts

    ENDDO ! Loop over # PFTs
    
 !! 3. Calculate maintenance respiration

    ! The maintenance respiration @tex $(gC.m^{-2}dt_sechiba^{-1})$ @endtex of each plant compartment for each PFT is 
    ! the biomass @tex $(gC.m^{-2})$ @endtex of the plant part multiplied by maintenance respiration 
    ! coefficient @tex $(g.g^{-1}dt_sechiba^{-1})$ @endtex, except that the maintenance respiration of leaves is 
    ! corrected by leaf area index (LAI) as follows:\n
    ! \latexonly     
    !   \input{resp5.tex} 
    ! \endlatexonly

    ! ibare_sechiba = 1, which means the there is only bare soil but not any PFT, consequently no LAI and
    !  maintenance respiration
    lai(:,ibare_sechiba) = zero
    resp_maint_part_radia(:,ibare_sechiba,:) = zero
    
    DO j = 2,nvm ! Loop over # PFTs
       
       ! 3.1 Maintenance respiration of the different plant parts
!JCMODIF
       IF ( .NOT. ok_LAIdev(j) ) THEN
           lai(:,j) = biomass(:,j,ileaf,icarbon) * sla_calc(:,j)
    !       lai(:,j) = biomass(:,j,ileaf) * sla(j)
       ENDIF

!ENDJCMODIF
       DO k = 1, nparts ! Loop over # plant parts

          IF ( k .EQ. ileaf ) THEN

             ! Leaves: respiration depends on leaf mass AND LAI.
!!$                WHERE ( (biomass(:,j,ileaf) > min_stomate) .AND. (lai(:,j) > 0.0) .AND. (lai(:,j) < val_exp) )
!!$                resp_maint_part_radia(:,j,k) = coeff_maint(:,j,k) * biomass(:,j,k) * &
!!$                        ( .3*lai(:,j) + 1.4*(1.-exp(-.5*lai(:,j))) ) / lai(:,j)
!!$             ELSEWHERE
!!$                resp_maint_part_radia(:,j,k) = 0.0
!!$             ENDWHERE
             DO i = 1, npts ! Loop over # pixels
                IF ( (biomass(i,j,ileaf,icarbon) > min_stomate) .AND. (lai(i,j) > min_stomate) ) THEN

!$                         IF (lai(i,j) < 100._r_std) THEN
!$                            resp_maint_part_radia(i,j,k) = coeff_maint(i,j,k) * biomass(i,j,k,icarbon) * &
!$                                 ( .3*lai(i,j) + 1.4*(1.-exp(-.5*lai(i,j))) ) / lai(i,j)
!$                         ELSE
!$                            resp_maint_part_radia(i,j,k) = coeff_maint(i,j,k) * biomass(i,j,k,icarbon) * &
!$                                 ( .3*lai(i,j) + 1.4 ) / lai(i,j)
!$                         ENDIF

                   ! Maintenance respiration is calculated as a fraction of biomass as defined by coeff_maint and 
                   ! is adjusted for the nitrogen effect through a third factor depending on LAI. The hypothesis 
                   ! here is that the vcmax (i.e. the nitrogen distribution) in the canopy decreases exponentially 
                   ! with LAI following the Beer-Lambert law with an asymptote defining the minimum of the function
                   ! at 30% of the LAI. The 1.4 parameter is an integration constant.
                   ! This method is also used in diffuco_trans_co2 2.4.1 for scaling vmax based on nitrogen reduction 
                   ! in the canopy.
                   resp_maint_part_radia(i,j,k) = coeff_maint(i,j,k) * biomass(i,j,k,icarbon) * &
                        ( maint_resp_min_vmax*lai(i,j) + maint_resp_coeff*(un - exp(-ext_coeff(j)*lai(i,j))) ) / lai(i,j)

!!!qcj++ peatland, following Druel et al. 2017 GMD paper, Eq.6
                   IF ( is_mosspeat(j) .AND. (humrel_month(i,j) .LT.humrel_mmin) ) THEN
                      resp_maint_part_radia(i,j,k)=resp_maint_part_radia(i,j,k) * &
                            ( vcmax_offset + ((1.0 - vcmax_offset) / humrel_mmin ) * humrel_month(i,j))
                   ENDIF

                   IF (resp_maint_part_radia(i,j,k)<0) THEN
                       WRITE(numout,*) "xuhui, resp_maint<0:"
                       WRITE(numout,*) 'k ',k
                       WRITE(numout,*) 'coeff_maint ',coeff_maint(i,j,k)
                       WRITE(numout,*) 'lai(i,j) ',lai(i,j)
                       WRITE(numout,*) 'biomass(i,j,k,icarbon)',biomass(i,j,k,icarbon)
                   ENDIF
                ELSE
                   resp_maint_part_radia(i,j,k) = zero
                ENDIF
             ENDDO ! Loop over # pixels
          ELSE

             resp_maint_part_radia(:,j,k) = coeff_maint(:,j,k) * biomass(:,j,k,icarbon)

          ENDIF

       ENDDO ! Loop over # plant parts

       ! 3.2 Total maintenance respiration of the plant
       !     VPP killer:
       !     resp_maint(:,j) = SUM( resp_maint_part(:,:), DIM=2 )

    ENDDO ! Loop over # PFTs

!    WRITE(numout,*) 'lai after stomate_resp: ',lai(1,12:14)

  END SUBROUTINE maint_respiration

END MODULE stomate_resp
