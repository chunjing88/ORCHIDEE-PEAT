! =================================================================================================================================
! MODULE       : stomate_turnover.f90
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        This module manages the end of the growing season and calculates herbivory and turnover of leaves, fruits, fine roots.
!! 
!!\n DESCRIPTION: This subroutine calculates leaf senescence due to climatic conditions or as a 
!! function of leaf age and new LAI, and subsequent turnover of the different plant biomass compartments (sections 1 to 6), 
!! herbivory (section 7), fruit turnover for trees (section 8) and sapwood conversion (section 9). 
!!
!! RECENT CHANGE(S): None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_stomate/stomate_turnover.f90 $
!! $Date: 2018-05-30 15:40:04 +0200 (Wed, 30 May 2018) $
!! $Revision: 5268 $
!! \n
!_ ================================================================================================================================

MODULE stomate_turnover

  ! modules used:
  USE xios_orchidee
  USE ioipsl_para
  USE stomate_data
  USE constantes
  USE pft_parameters

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC turn, turn_clear

  LOGICAL, SAVE                          :: firstcall_turnover = .TRUE.           !! first call (true/false)
!$OMP THREADPRIVATE(firstcall_turnover)

CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : turn_clear
!!
!>\BRIEF        Set flag ::firstcall_turnover to .TRUE., and therefore activate section 1 
!!              of subroutine turn which writes a message to the output.
!!                
!_ ================================================================================================================================

  SUBROUTINE turn_clear
    firstcall_turnover=.TRUE.
  END SUBROUTINE turn_clear


!! ================================================================================================================================
!! SUBROUTINE    : turn
!!
!>\BRIEF         Calculate turnover of leaves, roots, fruits and sapwood due to aging or climatic 
!!               induced senescence. Calculate herbivory.
!!
!! DESCRIPTION : This subroutine determines the turnover of leaves and fine roots (and stems for grasses)
!! and simulates following processes:
!! 1. Mean leaf age is calculated from leaf ages of separate leaf age classes. Should actually 
!!    be recalculated at the end of the routine, but it does not change too fast. The mean leaf 
!!    age is calculated using the following equation:
!!    \latexonly
!!    \input{turnover_lma_update_eqn1.tex}
!!    \endlatexonly
!!    \n 
!! 2. Meteorological senescence: the detection of the end of the growing season and shedding 
!!    of leaves, fruits and fine roots due to unfavourable meteorological conditions.
!!    The model distinguishes three different types of "climatic" leaf senescence, that do not 
!!    change the age structure: sensitivity to cold temperatures, to lack of water, or both. 
!!    If meteorological conditions are fulfilled, a flag ::senescence is set to TRUE. Note
!!    that evergreen species do not experience climatic senescence.
!!    Climatic senescence is triggered by sensitivity to cold temperatures where the critical 
!!    temperature for senescence is calculated using the following equation:
!!    \latexonly
!!    \input{turnover_temp_crit_eqn2.tex}
!!    \endlatexonly
!!    \n
!!    Climatic senescence is triggered by sensitivity to lack of water availability where the 
!!    moisture availability critical level is calculated using the following equation:
!!    \latexonly
!!    \input{turnover_moist_crit_eqn3.tex}
!!    \endlatexonly
!!    \n
!!    Climatic senescence is triggered by sensitivity to temperature or to lack of water where
!!    critical temperature and moisture availability are calculated as above.\n
!!    Trees in climatic senescence lose their fine roots at the same rate as they lose their leaves. 
!!    The rate of biomass loss of both fine roots and leaves is presribed through the equation:
!!    \latexonly
!!    \input{turnover_clim_senes_biomass_eqn4.tex}
!!    \endlatexonly
!!    \n
!!    with ::leaffall(j) a PFT-dependent time constant which is given in 
!!    ::stomate_constants. In grasses, leaf senescence is extended to the whole plant 
!!    (all carbon pools) except to its carbohydrate reserve.    
!! 3. Senescence due to aging: the loss of leaves, fruits and  biomass due to aging
!!    At a certain age, leaves fall off, even if the climate would allow a green plant
!!    all year round. Even if the meteorological conditions are favorable for leaf maintenance,
!!    plants, and in particular, evergreen trees, have to renew their leaves simply because the 
!!    old leaves become inefficient. Roots, fruits (and stems for grasses) follow leaves. 
!!    The ??senescence?? rate varies with leaf age. Note that plant is not declared senescent 
!!    in this case (wchich is important for allocation: if the plant loses leaves because of 
!!    their age, it can renew them). The leaf turnover rate due to aging of leaves is calculated
!!    using the following equation:
!!    \latexonly
!!    \input{turnover_age_senes_biomass_eqn5.tex}
!!    \endlatexonly
!!    \n
!!    Drop all leaves if there is a very low leaf mass during senescence. After this, the biomass 
!!    of different carbon pools both for trees and grasses is set to zero and the mean leaf age 
!!    is reset to zero. Finally, the leaf fraction and leaf age of the different leaf age classes 
!!    is set to zero. For deciduous trees: next to leaves, also fruits and fine roots are dropped.
!!    For grasses: all aboveground carbon pools, except the carbohydrate reserves are affected:
!! 4. Update the leaf biomass, leaf age class fraction and the LAI
!!    Older leaves will fall more frequently than younger leaves and therefore the leaf age 
!!    distribution needs to be recalculated after turnover. The fraction of biomass in each 
!!    leaf class is updated using the following equation:
!!    \latexonly
!!    \input{turnover_update_LeafAgeDistribution_eqn6.tex}
!!    \endlatexonly
!!    \n
!! 5. Simulate herbivory activity and update leaf and fruits biomass. Herbivore activity 
!!    affects the biomass of leaves and fruits as well as stalks (only for grasses).
!!    However, herbivores do not modify leaf age structure.
!! 6. Calculates fruit turnover for trees. Trees simply lose their fruits with a time 
!!    constant ::tau_fruit(j), that is set to 90 days for all PFTs in ::stomate_constants 
!! 7. Convert sapwood to heartwood for trees and update heart and sapwood above and 
!!    belowground biomass. Sapwood biomass is converted into heartwood biomass 
!!    with a time constant tau ::tau_sap(j) of 1 year. Note that this biomass conversion 
!!    is not added to "turnover" as the biomass is not lost. For the updated heartwood, 
!!    the sum of new heartwood above and new heartwood below after converting sapwood to 
!!    heartwood, is saved as ::hw_new(:). Creation of new heartwood decreases the age of 
!!    the plant ??carbon?? with a factor that is determined by: old heartwood ::hw_old(:) 
!!    divided by the new heartwood ::hw_new(:)
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLES: ::Biomass of leaves, fruits, fine roots and sapwood above (latter for grasses only), 
!!                        ::Update LAI, ::Update leaf age distribution with new leaf age class fraction 
!!
!! REFERENCE(S) : 
!! - Krinner, G., N. Viovy, N. de Noblet-Ducoudre, J. Ogee, J. Polcher, P. 
!! Friedlingstein, P. Ciais, S. Sitch and I.C. Prentice (2005), A dynamic global
!! vegetation model for studies of the coupled atmosphere-biosphere system, Global
!! Biogeochemical Cycles, 19, doi:10.1029/2003GB002199. 
!! - McNaughton, S. J., M. Oesterheld, D. A. Frank and K. J. Williams (1989), 
!! Ecosystem-level patterns of primary productivity and herbivory in terrestrial habitats, 
!! Nature, 341, 142-144, 1989. 
!! - Sitch, S., C. Huntingford, N. Gedney, P. E. Levy, M. Lomas, S. L. Piao, , Betts, R., Ciais, P., Cox, P., 
!! Friedlingstein, P., Jones, C. D., Prentice, I. C. and F. I. Woodward : Evaluation of the terrestrial carbon  
!! cycle, future plant geography and climate-carbon cycle feedbacks using 5 dynamic global vegetation 
!! models (dgvms), Global Change Biology, 14(9), 2015–2039, 2008. 
!!
!! FLOWCHART    : 
!! \latexonly
!! \includegraphics[scale=0.5]{turnover_flowchart_1.png}
!! \includegraphics[scale=0.5]{turnover_flowchart_2.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE turn (npts, dt, PFTpresent, &
       herbivores, &
       maxmoiavail_lastyear, minmoiavail_lastyear, &
       moiavail_week, moiavail_month, t2m_longterm, t2m_month, t2m_week, veget_cov_max, &
       gdd_from_growthinit, leaf_age, leaf_frac, age, lai, biomass, &
       turnover, senescence,turnover_time, &
!!! crops
       nrec, c_export, &
!!! end crops, xuhui
!gmjc
       sla_calc ,&
!end gmjc
       npp0_cumul) !!!qcj++ peatland
    !! 0. Variable and parameter declaration

    !! 0.1 Input variables 

    INTEGER(i_std), INTENT(in)                                 :: npts                 !! Domain size - number of grid cells 
                                                                                       !! (unitless) 
    REAL(r_std), INTENT(in)                                    :: dt                   !! time step (dt_days)
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                   :: PFTpresent           !! PFT exists (true/false)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: herbivores           !! time constant of probability of a leaf to 
                                                                                       !! be eaten by a herbivore (days) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: maxmoiavail_lastyear !! last year's maximum moisture availability 
                                                                                       !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: minmoiavail_lastyear !! last year's minimum moisture availability 
                                                                                       !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: moiavail_week        !! "weekly" moisture availability 
                                                                                       !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: moiavail_month       !! "monthly" moisture availability 
                                                                                       !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: t2m_longterm         !! "long term" 2 meter reference 
                                                                                       !! temperatures (K) 
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: t2m_month            !! "monthly" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: t2m_week             !! "weekly" 2 meter temperatures (K)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: veget_cov_max        !! "maximal" coverage fraction of a PFT (LAI 
                                                                                       !! -> infinity) on ground (unitless) 
    INTEGER(i_std), DIMENSION(npts, nvm), INTENT(in)            :: nrec
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: gdd_from_growthinit  !! gdd senescence for crop

    !! 0.2 Output variables

    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(out) :: turnover         !! Turnover @tex ($gC m^{-2}$) @endtex
    LOGICAL, DIMENSION(npts,nvm), INTENT(out)                  :: senescence           !! is the plant senescent? (true/false) 
                                                                                       !! (interesting only for deciduous trees: 
                                                                                       !! carbohydrate reserve) 
    REAL(r_std), DIMENSION(npts,nvm),INTENT(out)               :: c_export             !! c export (fruit & straws) from croplands
    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout)  :: leaf_age             !! age of the leaves (days)
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout)  :: leaf_frac            !! fraction of leaves in leaf age class 
                                                                                       !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: age                  !! age (years)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: lai                  !! leaf area index @tex ($m^2 m^{-2}$) 
                                                                                       !! @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: biomass        !! biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: turnover_time        !! turnover_time of grasses (days)
!gmjc
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: sla_calc
!end gmjc
!!!qcj++ peatland
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: npp0_cumul
    !! 0.4 Local  variables

    REAL(r_std), DIMENSION(npts,nvm)                           :: leaf_meanage         !! mean age of the leaves (days)
    REAL(r_std), DIMENSION(npts)                               :: dturnover            !! Intermediate variable for turnover ??
                                                                                       !! @tex ($gC m^{-2}$) @endtex 
    REAL(r_std), DIMENSION(npts)                               :: moiavail_crit        !! critical moisture availability, function 
                                                                                       !! of last year's moisture availability 
                                                                                       !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts)                               :: tl                   !! long term annual mean temperature, (C)
    REAL(r_std), DIMENSION(npts)                               :: t_crit               !! critical senescence temperature, function 
                                                                                       !! of long term annual temperature (K) 
    LOGICAL, DIMENSION(npts)                                   :: shed_rest            !! shed the remaining leaves? (true/false)
    REAL(r_std), DIMENSION(npts)                               :: sapconv              !! Sapwood conversion @tex ($gC m^{-2}$) 
                                                                                       !! @endtex 
    REAL(r_std), DIMENSION(npts)                               :: hw_old               !! old heartwood mass @tex ($gC m^{-2}$) 
                                                                                       !! @endtex 
    REAL(r_std), DIMENSION(npts)                               :: hw_new               !! new heartwood mass @tex ($gC m^{-2}$) 
                                                                                       !! @endtex 
    REAL(r_std), DIMENSION(npts)                               :: lm_old               !! old leaf mass @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nleafages)                     :: delta_lm             !! leaf mass change for each age class @tex 
                                                                                       !! ($gC m^{-2}$) @endtex 
    REAL(r_std), DIMENSION(npts)                               :: turnover_rate        !! turnover rate (unitless) 
    REAL(r_std), DIMENSION(npts,nvm)                           :: leaf_age_crit        !! critical leaf age (days)
    REAL(r_std), DIMENSION(npts,nvm)                           :: new_turnover_time    !! instantaneous turnover time (days)
    INTEGER(i_std)                                             :: j,m,k                !! Index (unitless)
    REAL(r_std), DIMENSION(npts,nvm)                           :: histvar              !! controls the history output
 
!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering turnover'

    !! 1. first call - output messages

    IF ( firstcall_turnover ) THEN

       IF (printlev >=2 ) THEN
          WRITE(numout,*) 'turnover:'
          WRITE(numout,*) ' > minimum mean leaf age for senescence (days) (::min_leaf_age_for_senescence) : '&
               ,min_leaf_age_for_senescence
       END IF
       firstcall_turnover = .FALSE.


    ENDIF

    !! 2. Initializations 

    !! 2.1 set output to zero
    turnover(:,:,:,:) = zero
    c_export(:,:) = zero
    new_turnover_time(:,:) = zero
    senescence(:,:) = .FALSE.

    !! 2.2 Recalculate mean leaf age
    !      Mean leaf age is recalculated from leaf ages of separate leaf age classes. Should actually be recalculated at the 
    !      end of this routine, but it does not change too fast.
    !      The mean leaf age is calculated using the following equation:
    !      \latexonly
    !      \input{turnover_lma_update_eqn1.tex}
    !      \endlatexonly
    !      \n
    leaf_meanage(:,:) = zero

    DO m = 1, nleafages
       leaf_meanage(:,:) = leaf_meanage(:,:) + leaf_age(:,:,m) * leaf_frac(:,:,m)
    ENDDO

    !! 3. Climatic senescence

    !     Three different types of "climatic" leaf senescence,
    !     that do not change the age structure. 
    DO j = 2,nvm ! Loop over # PFTs

       !! 3.1 Determine if there is climatic senescence. 
       !      The climatic senescence can be of three types:
       !      sensitivity to cold temperatures, to lack of water, or both. If meteorological conditions are 
       !      fulfilled, a flag senescence is set to TRUE.
       !      Evergreen species do not experience climatic senescence.

       SELECT CASE ( senescence_type(j) )


       CASE ('crop' )!for crop senescence is based on a GDD criterium as in crop models
          WHERE ( ( biomass(:,j,ileaf,icarbon) .GT. zero ) .AND. &
               ( leaf_meanage(:,j) .GT. min_leaf_age_for_senescence(j) ) .AND. &
               ( gdd_from_growthinit(:,j) .GT.  gdd_senescence(j)))
            ! note that orchidee-crop does not utilize this any longer, xuhui
             senescence(:,j) = .TRUE.
          ENDWHERE

       CASE ( 'cold' )

          !! 3.1.1 Summergreen species: Climatic senescence is triggered by sensitivity to cold temperatures
          !        Climatic senescence is triggered by sensitivity to cold temperatures as follows: 
          !        If biomass is large enough (i.e. when it is greater than zero), 
          !        AND (i.e. when leaf mean age is above a certain PFT-dependent treshold ::min_leaf_age_for_senescence,
          !        which is given in ::stomate_constants),      
          !        AND the monthly temperature is low enough (i.e. when monthly temperature ::t2m_month(:) is below a critical 
          !        temperature ::t_crit(:), which is calculated in this module),
          !        AND the temperature tendency is negative (i.e. when weekly temperatures ::t2m_week(:) are lower than monthly 
          !        temperatures ::t2m_month(:))
          !        If these conditions are met, senescence is set to TRUE.
          !
          !        The critical temperature for senescence is calculated using the following equation:
          !        \latexonly
          !        \input{turnover_temp_crit_eqn2.tex}
          !        \endlatexonly
          !        \n
          !
          ! Critical temperature for senescence may depend on long term annual mean temperature
          tl(:) = t2m_longterm(:) - ZeroCelsius
          t_crit(:) = ZeroCelsius + senescence_temp(j,1) + &
               tl(:) * senescence_temp(j,2) + &
               tl(:)*tl(:) * senescence_temp(j,3)

          WHERE ( ( biomass(:,j,ileaf,icarbon) .GT. zero ) .AND. &
               ( leaf_meanage(:,j) .GT. min_leaf_age_for_senescence(j) ) .AND. &
               ( t2m_month(:) .LT. t_crit(:) ) .AND. ( t2m_week(:) .LT. t2m_month(:) ) )


             senescence(:,j) = .TRUE.

          ENDWHERE

       CASE ( 'dry' )

          !! 3.1.2 Raingreen species: Climatic senescence is triggered by sensitivity to lack of water availability 
          !        Climatic senescence is triggered by sensitivity to lack of water availability as follows:  
          !        If biomass is large enough (i.e. when it is greater than zero), 
          !        AND (i.e. when leaf mean age is above a certain PFT-dependent treshold ::min_leaf_age_for_senescence,
          !        which is given in ::stomate_constants),      
          !        AND the moisture availability drops below a critical level (i.e. when weekly moisture availability 
          !        ::moiavail_week(:,j) is below a critical moisture availability ::moiavail_crit(:),
          !        which is calculated in this module), 
          !        If these conditions are met, senescence is set to TRUE.
          !
          !        The moisture availability critical level is calculated using the following equation:
          !        \latexonly
          !        \input{turnover_moist_crit_eqn3.tex}
          !        \endlatexonly
          !        \n
          moiavail_crit(:) = &
               MIN( MAX( minmoiavail_lastyear(:,j) + hum_frac(j) * &
               ( maxmoiavail_lastyear(:,j) - minmoiavail_lastyear(:,j) ), &
               senescence_hum(j) ), &
               nosenescence_hum(j) )

          WHERE ( ( biomass(:,j,ileaf,icarbon) .GT. zero ) .AND. &
               ( leaf_meanage(:,j) .GT. min_leaf_age_for_senescence(j) ) .AND. &
               ( moiavail_week(:,j) .LT. moiavail_crit(:) ) )

             senescence(:,j) = .TRUE.

          ENDWHERE

       CASE ( 'mixed' )

          !! 3.1.3 Mixed criterion: Climatic senescence is triggered by sensitivity to temperature or to lack of water  
          !        Climatic senescence is triggered by sensitivity to temperature or to lack of water availability as follows:
          !        If biomass is large enough (i.e. when it is greater than zero), 
          !        AND (i.e. when leaf mean age is above a certain PFT-dependent treshold ::min_leaf_age_for_senescence,
          !        which is given in ::stomate_constants),      
          !        AND the moisture availability drops below a critical level (i.e. when weekly moisture availability 
          !        ::moiavail_week(:,j) is below a critical moisture availability ::moiavail_crit(:), calculated in this module), 
          !        OR 
          !        the monthly temperature is low enough (i.e. when monthly temperature ::t2m_month(:) is below a critical 
          !        temperature ::t_crit(:), calculated in this module),
          !        AND the temperature tendency is negative (i.e. when weekly temperatures ::t2m_week(:) are lower than 
          !        monthly temperatures ::t2m_month(:)).
          !        If these conditions are met, senescence is set to TRUE.
          moiavail_crit(:) = &
               MIN( MAX( minmoiavail_lastyear(:,j) + hum_frac(j) * &
               (maxmoiavail_lastyear(:,j) - minmoiavail_lastyear(:,j) ), &
               senescence_hum(j) ), &
               nosenescence_hum(j) )

          tl(:) = t2m_longterm(:) - ZeroCelsius
          t_crit(:) = ZeroCelsius + senescence_temp(j,1) + &
               tl(:) * senescence_temp(j,2) + &
               tl(:)*tl(:) * senescence_temp(j,3)

          IF ( is_tree(j) ) THEN
             ! critical temperature for senescence may depend on long term annual mean temperature
             WHERE ( ( biomass(:,j,ileaf,icarbon) .GT. zero ) .AND. &
                  ( leaf_meanage(:,j) .GT. min_leaf_age_for_senescence(j) ) .AND. &
                  ( ( moiavail_week(:,j) .LT. moiavail_crit(:) ) .OR. &
                  ( ( t2m_month(:) .LT. t_crit(:) ) .AND. ( t2m_week(:) .LT. t2m_month(:) ) ) ) )
                senescence(:,j) = .TRUE.
             ENDWHERE
          ELSE

            turnover_time(:,j) = max_turnover_time(j)

            WHERE ( ( biomass(:,j,ileaf,icarbon) .GT. zero ) .AND. &
                 ( leaf_meanage(:,j) .GT. min_leaf_age_for_senescence(j) ) .AND. &
                 ( ( moiavail_week(:,j) .LT. moiavail_crit(:) )))
                turnover_time(:,j) = max_turnover_time(j) * &
                     (1.-   (1.- (moiavail_week(:,j)/  moiavail_crit(:)))**2)		
            ENDWHERE
            WHERE ( turnover_time(:,j) .LT. min_turnover_time(j) )		 
	       turnover_time(:,j) = min_turnover_time(j)			 
            ENDWHERE 									 

            WHERE ((( biomass(:,j,ileaf,icarbon) .GT. zero ) .AND. &
                ( leaf_meanage(:,j) .GT. min_leaf_age_for_senescence(j) ) .AND. &
                ((t2m_month(:) .LT. t_crit(:)) .AND. (lai(:,j) .GT. lai_max(j)/4.) .OR. &
                (t2m_month(:) .LT. ZeroCelsius)) .AND. ( t2m_week(:) .LT. t2m_month(:) )))
               turnover_time(:,j)= leaffall(j)
            ENDWHERE
!gmjc
          IF (is_grassland_manag(j)) THEN
            !WHERE (lai(:,j) .LT. 0.25)
            !   turnover_time(:,j)= max_turnover_time(j)
            !ENDWHERE
            WHERE (lai(:,j) .GT. 2.55)
               turnover_time(:,j)= MAX(45.0,(85.0-lai(:,j)*10.0))
            ENDWHERE
          ENDIF
!end gmjc
         ENDIF
       !! Evergreen species do not experience climatic senescence
       CASE ( 'none' )

          
       !! In case no climatic senescence type is recognized.
       CASE default

          WRITE(numout,*) '  turnover: don''t know how to treat this PFT.'
          WRITE(numout,*) '  number (::j) : ',j
          WRITE(numout,*) '  senescence type (::senescence_type(j)) : ',senescence_type(j)

          CALL ipslerr_p(3,"turn","Dont know how to treat this PFT.","","")

       END SELECT

       !! 3.2 Drop leaves and roots, plus stems and fruits for grasses

       ! for the new routine of crop.
       ! xuhui: whatever the senescence_type is chosen, when initiated ORCHIDEE-crop
       ! the following turnover process will be used 
 
       IF (ok_LAIdev(j)) THEN           ! using the new routine of crops
          
          ! for those pixels where harvesting, senescence occurs
          WHERE ((nrec(:, j) .gt. 0))  ! where harvested 
             senescence(:, j) = .TRUE. 
             turnover_time(:, j) = dt      !how many days for litters
          ELSEWHERE
             turnover_time(:, j) = max_turnover_time(j)
          ENDWHERE
       ENDIF

       IF ( is_tree(j) ) THEN

          !! 3.2.1 Trees in climatic senescence lose their fine roots at the same rate as they lose their leaves. 
          !        The rate of biomass loss of both fine roots and leaves is presribed through the equation:
          !        \latexonly
          !        \input{turnover_clim_senes_biomass_eqn4.tex}
          !        \endlatexonly
          !        \n
          !         with ::leaffall(j) a PFT-dependent time constant which is given in ::stomate_constants),
          WHERE ( senescence(:,j) )

             turnover(:,j,ileaf,icarbon) = biomass(:,j,ileaf,icarbon) * dt / leaffall(j)
             turnover(:,j,iroot,icarbon) = biomass(:,j,iroot,icarbon) * dt / leaffall(j)

          ENDWHERE

       ELSEIF (ok_LAIdev(j)) THEN ! for STICS

          WHERE (turnover_time(:,j) .LT. max_turnover_time(j)) 
             turnover(:,j,ileaf,icarbon) = biomass(:,j,ileaf,icarbon) * dt / turnover_time(:,j)
             turnover(:,j,isapabove,icarbon) = biomass(:,j,isapabove,icarbon) * dt / turnover_time(:,j)
             turnover(:,j,iroot,icarbon) = biomass(:,j,iroot,icarbon) * dt / turnover_time(:,j) 
             turnover(:,j,ifruit,icarbon) = biomass(:,j,ifruit,icarbon) * dt / turnover_time(:,j)
          ELSEWHERE
             turnover(:,j,ileaf,icarbon)= zero
             turnover(:,j,isapabove,icarbon) = zero
             turnover(:,j,iroot,icarbon) = zero
             turnover(:,j,ifruit,icarbon) = zero
          ENDWHERE

          IF ( ANY(biomass(1,j,:,icarbon)<0) ) THEN
              WRITE(numout,*) 'biomass low0 in turnover '
              WRITE(numout,*) 'biomass(1,j,:,icarbon): ',biomass(1,j,:,icarbon) 
              WRITE(numout,*) 'turnover(1,j,:,icarbon): ',turnover(1,j,:,icarbon)
          ENDIF

       ELSE  !grass or super-grass

          !! 3.2.2 In grasses, leaf senescence is extended to the whole plant 
          !        In grasses, leaf senescence is extended to the whole plant (all carbon pools) except to its
          !        carbohydrate reserve.      

          IF (senescence_type(j) .EQ. 'crop') THEN
             ! 3.2.2.1 crops with 'crop' phenological model
             WHERE ( senescence(:,j) )
                turnover(:,j,ileaf,icarbon) = biomass(:,j,ileaf,icarbon) * dt / leaffall(j)
                turnover(:,j,iroot,icarbon) = biomass(:,j,iroot,icarbon) * dt / leaffall(j)
                turnover(:,j,isapabove,icarbon) = biomass(:,j,isapabove,icarbon) * dt / leaffall(j)
                turnover(:,j,ifruit,icarbon) = biomass(:,j,ifruit,icarbon) * dt /leaffall(j)
             ENDWHERE
          ELSE
          ! 3.2.2.2 grass or crops based on 'mixed' phenological model
             WHERE (turnover_time(:,j) .LT. max_turnover_time(j)) 
                turnover(:,j,ileaf,icarbon) = biomass(:,j,ileaf,icarbon) * dt / turnover_time(:,j)
                turnover(:,j,isapabove,icarbon) = biomass(:,j,isapabove,icarbon) * dt / turnover_time(:,j)
                turnover(:,j,iroot,icarbon) = biomass(:,j,iroot,icarbon) * dt / turnover_time(:,j) 
                turnover(:,j,ifruit,icarbon) = biomass(:,j,ifruit,icarbon) * dt / turnover_time(:,j)
             ENDWHERE
          ENDIF
       ENDIF      ! tree/grass
       biomass(:,j,ileaf,icarbon) = biomass(:,j,ileaf,icarbon) - turnover(:,j,ileaf,icarbon)
       biomass(:,j,isapabove,icarbon) = biomass(:,j,isapabove,icarbon) - turnover(:,j,isapabove,icarbon)
       biomass(:,j,iroot,icarbon) = biomass(:,j,iroot,icarbon) - turnover(:,j,iroot,icarbon)
       biomass(:,j,ifruit,icarbon) = biomass(:,j,ifruit,icarbon) - turnover(:,j,ifruit,icarbon)

       IF (ok_LAIdev(j)) THEN
           ! part of the "turnovered" carbon ( 1-prc_residual of sapabove & leaf and all fruit) are exported from the ecosystem
           c_export(:,j) =  turnover(:,j,ifruit,icarbon) + (1 - prc_residual) * &
                            (turnover(:,j,ileaf,icarbon)+turnover(:,j,isapabove,icarbon))
           turnover(:,j,ifruit,icarbon) = 0.
           turnover(:,j,ileaf,icarbon) = prc_residual * turnover(:,j,ileaf,icarbon)
           turnover(:,j,isapabove,icarbon) = prc_residual * turnover(:,j,isapabove,icarbon)
       ENDIF
    ENDDO        ! loop over PFTs

    !! 4. Leaf fall
    !     At a certain age, leaves fall off, even if the climate would allow a green plant
    !     all year round. Even if the meteorological conditions are favorable for leaf maintenance,
    !     plants, and in particular, evergreen trees, have to renew their leaves simply because the 
    !     old leaves become inefficient.   
    !     Roots, fruits (and stems) follow leaves. The decay rate varies with leaf age.
    !     Note that plant is not declared senescent in this case (wchich is important for allocation:
    !     if the plant loses leaves because of their age, it can renew them).
    !
    !     The leaf turnover rate due to aging of leaves is calculated using the following equation:
    !     \latexonly
    !     \input{turnover_age_senes_biomass_eqn5.tex}
    !     \endlatexonly
    !     \n
    DO j = 2,nvm ! Loop over # PFTs

       !! save old leaf mass
       lm_old(:) = biomass(:,j,ileaf,icarbon)

       !! initialize leaf mass change in age class
       delta_lm(:,:) = zero

       IF ( is_tree(j) .OR. (.NOT. natural(j)) ) THEN

          !! 4.1 Trees: leaves, roots, fruits roots and fruits follow leaves.

          !! 4.1.1 Critical age: prescribed for trees
          leaf_age_crit(:,j) = leafagecrit(j)

       ELSE

          !! 4.2 Grasses: leaves, roots, fruits, sap follow leaves.

          !! 4.2.1 Critical leaf age depends on long-term temperature
          !        Critical leaf age depends on long-term temperature
          !        generally, lower turnover in cooler climates.
          leaf_age_crit(:,j) = &
               MIN( leafagecrit(j) * leaf_age_crit_coeff(1) , &
               MAX( leafagecrit(j) * leaf_age_crit_coeff(2) , &
               leafagecrit(j) - leaf_age_crit_coeff(3) * &
               ( t2m_longterm(:)-ZeroCelsius - leaf_age_crit_tref ) ) )

       END IF
       

       IF (.NOT. ok_LAIdev(j)) THEN   
           ! 4.2.2 Loop over leaf age classes
           DO m = 1, nleafages
    
              turnover_rate(:) = zero
    
              WHERE ( leaf_age(:,j,m) .GT. leaf_age_crit(:,j)/2. )
    
                 turnover_rate(:) =  &
                      MIN( 0.99_r_std, dt / ( leaf_age_crit(:,j) * &
                      ( leaf_age_crit(:,j) / leaf_age(:,j,m) )**quatre ) )
    
              ENDWHERE
              
              dturnover(:) = biomass(:,j,ileaf,icarbon) * leaf_frac(:,j,m) * turnover_rate(:)
              turnover(:,j,ileaf,icarbon) = turnover(:,j,ileaf,icarbon) + dturnover(:)
              biomass(:,j,ileaf,icarbon) = biomass(:,j,ileaf,icarbon) - dturnover(:)
    
              ! save leaf mass change
              delta_lm(:,m) = - dturnover(:)
              
              dturnover(:) = biomass(:,j,iroot,icarbon) * leaf_frac(:,j,m) * turnover_rate(:)
              turnover(:,j,iroot,icarbon) = turnover(:,j,iroot,icarbon) + dturnover(:)
              biomass(:,j,iroot,icarbon) = biomass(:,j,iroot,icarbon) - dturnover(:)
              
              dturnover(:) = biomass(:,j,ifruit,icarbon) * leaf_frac(:,j,m) * turnover_rate(:)
              turnover(:,j,ifruit,icarbon) = turnover(:,j,ifruit,icarbon) + dturnover(:)
              biomass(:,j,ifruit,icarbon) = biomass(:,j,ifruit,icarbon) - dturnover(:)
              
              IF (.NOT. is_tree(j)) THEN
                 dturnover(:) = biomass(:,j,isapabove,icarbon) * leaf_frac(:,j,m) * turnover_rate(:)
                 turnover(:,j,isapabove,icarbon) = turnover(:,j,isapabove,icarbon) + dturnover(:)
                 biomass(:,j,isapabove,icarbon) = biomass(:,j,isapabove,icarbon) - dturnover(:)
              ENDIF
              
           ENDDO
       ENDIF

       !! 4.3 Recalculate the fraction of leaf biomass in each leaf age class.
       !      Older leaves will fall more fast than younger leaves and therefore 
       !      the leaf age distribution needs to be recalculated after turnover. 
       !      The fraction of biomass in each leaf class is updated using the following equation:
       !      \latexonly
       !      \input{turnover_update_LeafAgeDistribution_eqn6.tex}
       !      \endlatexonly
       !      \n
       !
       !      new fraction = new leaf mass of that fraction / new total leaf mass
       !                   = (old fraction*old total leaf mass ::lm_old(:) + biomass change of that fraction ::delta_lm(:,m)  ) /
       !                     new total leaf mass ::biomass(:,j,ileaf
       DO m = 1, nleafages
          
          WHERE ( biomass(:,j,ileaf,icarbon) .GT. min_sechiba )
             leaf_frac(:,j,m) = ( leaf_frac(:,j,m)*lm_old(:) + delta_lm(:,m) ) / biomass(:,j,ileaf,icarbon)
          ELSEWHERE
             leaf_frac(:,j,m) = zero
          ENDWHERE

       ENDDO

!!!qcj++ peatland, following Druel et al. 2017 GMD paper, Eq.4
       IF ( is_mosspeat(j) ) THEN
          turnover_rate(:) = zero
          WHERE ( lai(:,j) .gt. llaimax(j) )
             turnover_rate(:) = EXP(llai_coef * (lai(:,j) - llaimax(j)))-1.0
          ENDWHERE
          dturnover(:) = biomass(:,j,ileaf,icarbon) * turnover_rate(:)
          turnover(:,j,ileaf,icarbon) = turnover(:,j,ileaf,icarbon) + dturnover(:)
          biomass(:,j,ileaf,icarbon) = biomass(:,j,ileaf,icarbon) - dturnover(:)

          dturnover(:) = biomass(:,j,ifruit,icarbon) * turnover_rate(:)
          turnover(:,j,ifruit,icarbon) = turnover(:,j,ifruit,icarbon) + dturnover(:)
          biomass(:,j,ifruit,icarbon) = biomass(:,j,ifruit,icarbon) - dturnover(:)
       ENDIF

       IF ( ANY(biomass(1,j,:,icarbon)<0) ) THEN
           WRITE(numout,*) 'biomass low0 in turnover after leaf fall '
           WRITE(numout,*) 'biomass(1,j,:,icarbon): ',biomass(1,j,:,icarbon) 
           WRITE(numout,*) 'turnover(1,j,:,icarbon): ',turnover(1,j,:,icarbon)
       ENDIF
    ENDDO         ! loop over PFTs


    !! 5. New (provisional) LAI 
    !     ::lai(:,j) is determined from the leaf biomass ::biomass(:,j,ileaf,icarbon) and the 
    !     specific leaf surface :: sla(j) (m^2 gC^{-1})
    !     The leaf area index is updated using the following equation:
    !     \latexonly
    !     \input{turnover_update_LAI_eqn7.tex}
    !     \endlatexonly
    !     \n

    !    lai(:,ibare_sechiba) = zero
    !    DO j = 2, nvm ! Loop over # PFTs
    !        lai(:,j) = biomass(:,j,ileaf,icarbon) * sla(j)
    !    ENDDO

    !! 6. Definitely drop all leaves if there is a very low leaf mass during senescence.

    !     Both for deciduous trees and grasses same conditions are checked:
    !     If biomass is large enough (i.e. when it is greater than zero), 
    !     AND when senescence is set to true
    !     AND the leaf biomass drops below a critical minimum biomass level (i.e. when it is lower than half
    !     the minimum initial LAI ::lai_initmin(j) divided by the specific leaf area ::sla(j),
    !     ::lai_initmin(j) is set to 0.3 in stomate_data.f90 and sla is a constant that is set to 0.015366 m2/gC), 
    !     If these conditions are met, the flag ::shed_rest(:) is set to TRUE.
    !
    !     After this, the biomass of different carbon pools both for trees and grasses is set to zero
    !     and the mean leaf age is reset to zero.
    !     Finally, the leaf fraction and leaf age of the different leaf age classes is set to zero.
    DO j = 2,nvm ! Loop over # PFTs

       shed_rest(:) = .FALSE.

       !! 6.1 For deciduous trees: next to leaves, also fruits and fine roots are dropped 
       !      For deciduous trees: next to leaves, also fruits and fine roots are dropped: fruit ::biomass(:,j,ifruit) 
       !      and fine root ::biomass(:,j,iroot) carbon pools are set to zero.
       IF ( is_tree(j) .AND. ( senescence_type(j) .NE. 'none' ) .AND. (.NOT. ok_LAIdev(j)) ) THEN

          ! check whether we shed the remaining leaves
          WHERE ( ( biomass(:,j,ileaf,icarbon) .GT. zero ) .AND. senescence(:,j) .AND. &
!JCMODIF
!               ( biomass(:,j,ileaf) .LT. (lai_initmin(j) / 2.)/sla(j) )             )
               ( biomass(:,j,ileaf,icarbon) .LT. (lai_initmin(j) / 2.)/sla_calc(:,j) ) )
!ENDJCMODIF
             shed_rest(:) = .TRUE.

             turnover(:,j,ileaf,icarbon)  = turnover(:,j,ileaf,icarbon) + biomass(:,j,ileaf,icarbon)
             turnover(:,j,iroot,icarbon)  = turnover(:,j,iroot,icarbon) + biomass(:,j,iroot,icarbon)
             turnover(:,j,ifruit,icarbon) = turnover(:,j,ifruit,icarbon) + biomass(:,j,ifruit,icarbon)

             biomass(:,j,ileaf,icarbon)  = zero
             biomass(:,j,iroot,icarbon)  = zero
             biomass(:,j,ifruit,icarbon) = zero

             ! reset leaf age and lai
             leaf_meanage(:,j) = zero
             lai(:,j) = zero
          ENDWHERE

       ENDIF

       !! 6.2 For grasses: all aboveground carbon pools, except the carbohydrate reserves are affected: 
       !      For grasses: all aboveground carbon pools, except the carbohydrate reserves are affected: 
       !      fruit ::biomass(:,j,ifruit,icarbon), fine root ::biomass(:,j,iroot,icarbon) and sapwood above 
       !      ::biomass(:,j,isapabove,icarbon) carbon pools are set to zero. 
       IF (( .NOT. is_tree(j)) .AND. (.NOT. ok_LAIdev(j)) ) THEN

          ! Shed the remaining leaves if LAI very low.
          WHERE ( ( biomass(:,j,ileaf,icarbon) .GT. zero ) .AND. senescence(:,j) .AND. &
!JCMODIF
!               (  biomass(:,j,ileaf) .LT. (lai_initmin(j) / 2.)/sla(j) ))
               (  biomass(:,j,ileaf,icarbon) .LT. (lai_initmin(j) / 2.)/sla_calc(:,j) ))
!ENDJCMODIF
             shed_rest(:) = .TRUE.

             turnover(:,j,ileaf,icarbon) = turnover(:,j,ileaf,icarbon) + biomass(:,j,ileaf,icarbon)
             turnover(:,j,isapabove,icarbon) = turnover(:,j,isapabove,icarbon) + biomass(:,j,isapabove,icarbon)
             turnover(:,j,iroot,icarbon) = turnover(:,j,iroot,icarbon) + biomass(:,j,iroot,icarbon)
             turnover(:,j,ifruit,icarbon) = turnover(:,j,ifruit,icarbon) + biomass(:,j,ifruit,icarbon)

             biomass(:,j,ileaf,icarbon) = zero
             biomass(:,j,isapabove,icarbon) = zero
             biomass(:,j,iroot,icarbon) = zero
             biomass(:,j,ifruit,icarbon) = zero

             ! reset leaf age and lai 
             leaf_meanage(:,j) = zero
             lai(:,j) = zero
          ENDWHERE

       ENDIF
       IF (printlev>=4) THEN
          IF ( ANY(biomass(1,j,:,icarbon)<0) ) THEN
              WRITE(numout,*) 'biomass low0 in turnover after leaf drop '
              WRITE(numout,*) 'biomass(1,j,:,icarbon): ',biomass(1,j,:,icarbon)
              WRITE(numout,*) 'turnover(1,j,:,icarbon):',turnover(1,j,:,icarbon)
          ENDIF
       ENDIF
       !! 6.3 Reset the leaf age structure: the leaf fraction and leaf age of the different leaf age classes is set to zero.
      
       DO m = 1, nleafages

          WHERE ( shed_rest(:) )

             leaf_age(:,j,m) = zero
             leaf_frac(:,j,m) = zero

          ENDWHERE

       ENDDO

!!!qcj++ peatland
       !! 6.4 For moss, following Druel et al., 2017 GMD paper, Eq.3
       IF ( is_mosspeat(j) ) THEN
          WHERE ( (npp0_cumul(:,j) .GT. npp0_c(1)) .AND. (npp0_cumul(:,j) .LE. npp0_c(2)) )
             dturnover(:) = biomass(:,j,ileaf,icarbon) * npp0_c(4) * &
                            & (npp0_cumul(:,j)-npp0_c(1))/ (npp0_c(2)-npp0_c(1))
             turnover(:,j,ileaf,icarbon) = turnover(:,j,ileaf,icarbon) + dturnover(:)
             biomass(:,j,ileaf,icarbon) = biomass(:,j,ileaf,icarbon) - dturnover(:)

             dturnover(:) = biomass(:,j,ifruit,icarbon) * npp0_c(4) * &
                            & (npp0_cumul(:,j)-npp0_c(1))/ (npp0_c(2)-npp0_c(1))
             turnover(:,j,ifruit,icarbon) = turnover(:,j,ifruit,icarbon) + dturnover(:)
             biomass(:,j,ifruit,icarbon) = biomass(:,j,ifruit,icarbon) - dturnover(:)
          ELSEWHERE ( (npp0_cumul(:,j) .GT. npp0_c(2)) .AND. (npp0_cumul(:,j) .LT. npp0_c(3)) )
             dturnover(:) = biomass(:,j,ileaf,icarbon) * npp0_c(4) * &
                            & (npp0_c(3)-npp0_cumul(:,j))/ (npp0_c(3)-npp0_c(2))
             turnover(:,j,ileaf,icarbon) = turnover(:,j,ileaf,icarbon) + dturnover(:)
             biomass(:,j,ileaf,icarbon) = biomass(:,j,ileaf,icarbon) - dturnover(:)

             dturnover(:) = biomass(:,j,ifruit,icarbon) * npp0_c(4) * &
                            & (npp0_c(3)-npp0_cumul(:,j))/ (npp0_c(3)-npp0_c(2))
             turnover(:,j,ifruit,icarbon) = turnover(:,j,ifruit,icarbon) + dturnover(:)
             biomass(:,j,ifruit,icarbon) = biomass(:,j,ifruit,icarbon) - dturnover(:)
          ENDWHERE
       ENDIF

    ENDDO          ! loop over PFTs
    
    !! 7. Herbivore activity: elephants, cows, gazelles but no lions.
 
    !     Herbivore activity affects the biomass of leaves and fruits as well 
    !     as stalks (only for grasses). Herbivore activity does not modify leaf 
    !     age structure. Herbivores ::herbivores(:,j) is the time constant of 
    !     probability of a leaf to be eaten by a herbivore, and is calculated in 
    !     ::stomate_season. following Mc Naughton et al. [1989].

    IF ( ok_herbivores ) THEN

       ! If the herbivore activity is allowed (if ::ok_herbivores is true, which is set in run.def), 
       ! remove the amount of biomass consumed by herbivory from the leaf biomass ::biomass(:,j,ileaf,icarbon) and 
       ! the fruit biomass ::biomass(:,j,ifruit,icarbon).
       ! The daily amount consumed equals the biomass multiplied by 1 day divided by the time constant ::herbivores(:,j).
       DO j = 2,nvm ! Loop over # PFTs

          IF ( is_tree(j) ) THEN

             !! For trees: only the leaves and fruit carbon pools are affected

             WHERE (biomass(:,j,ileaf,icarbon) .GT. zero)
                ! added by shilong
                WHERE (herbivores(:,j).GT. min_sechiba)
                   dturnover(:) = biomass(:,j,ileaf,icarbon) * dt / herbivores(:,j)
                   turnover(:,j,ileaf,icarbon) = turnover(:,j,ileaf,icarbon) + dturnover(:)
                   biomass(:,j,ileaf,icarbon) = biomass(:,j,ileaf,icarbon) - dturnover(:)

                   dturnover(:) = biomass(:,j,ifruit,icarbon) * dt / herbivores(:,j)
                   turnover(:,j,ifruit,icarbon) = turnover(:,j,ifruit,icarbon) + dturnover(:)
                   biomass(:,j,ifruit,icarbon) = biomass(:,j,ifruit,icarbon) - dturnover(:)
                ENDWHERE
             ENDWHERE

          ELSE

             ! For grasses: all aboveground carbon pools are affected: leaves, fruits and sapwood above
             WHERE ( biomass(:,j,ileaf,icarbon) .GT. zero )
                ! added by shilong
                WHERE (herbivores(:,j) .GT. min_sechiba)
                   dturnover(:) = biomass(:,j,ileaf,icarbon) * dt / herbivores(:,j)
                   turnover(:,j,ileaf,icarbon) = turnover(:,j,ileaf,icarbon) + dturnover(:)
                   biomass(:,j,ileaf,icarbon) = biomass(:,j,ileaf,icarbon) - dturnover(:)

                   dturnover(:) = biomass(:,j,isapabove,icarbon) * dt / herbivores(:,j)
                   turnover(:,j,isapabove,icarbon) = turnover(:,j,isapabove,icarbon) + dturnover(:)
                   biomass(:,j,isapabove,icarbon) = biomass(:,j,isapabove,icarbon) - dturnover(:)

                   dturnover(:) = biomass(:,j,ifruit,icarbon) * dt / herbivores(:,j)
                   turnover(:,j,ifruit,icarbon) = turnover(:,j,ifruit,icarbon) + dturnover(:)
                   biomass(:,j,ifruit,icarbon) = biomass(:,j,ifruit,icarbon) - dturnover(:)
                ENDWHERE

             ENDWHERE

          ENDIF  ! tree/grass?

       ENDDO    ! loop over PFTs

    ENDIF ! end herbivores

    !! 8. Fruit turnover for trees

    !     Fruit turnover for trees: trees simply lose their fruits with a time constant ::tau_fruit(j), 
    !     that is set to 90 days for all PFTs in ::stomate_constants

    DO k = 1,nelements 
       DO j = 2,nvm ! Loop over # PFTs
          IF ( is_tree(j) ) THEN

             dturnover(:) = biomass(:,j,ifruit,k) * dt / tau_fruit(j)
             turnover(:,j,ifruit,k) = turnover(:,j,ifruit,k) + dturnover(:)
             biomass(:,j,ifruit,k) = biomass(:,j,ifruit,k) - dturnover(:)
             
          ENDIF
       ENDDO       ! loop over PFTs
    END DO

    !! 9 Conversion of sapwood to heartwood both for aboveground and belowground sapwood and heartwood.

    !   Following LPJ (Sitch et al., 2003), sapwood biomass is converted into heartwood biomass 
    !   with a time constant tau ::tau_sap(j) of 1 year.
    !   Note that this biomass conversion is not added to "turnover" as the biomass is not lost!
    DO j = 2,nvm ! Loop over # PFTs

       IF ( is_tree(j) ) THEN

          !! For the recalculation of age in 9.2 (in case the vegetation is not dynamic ie. ::ok_dgvm is FALSE), 
          !! the heartwood above and below is stored in ::hw_old(:).
!!!qcj++ peatland
          IF ( (.NOT. ok_dgvm) .AND. (.NOT. ok_dgvm_peat) ) THEN 
             hw_old(:) = biomass(:,j,iheartabove,icarbon) + biomass(:,j,iheartbelow,icarbon)
          ELSEIF ( (.NOT. ok_dgvm) .AND. ok_dgvm_peat .AND. (.NOT. is_peat(j)) )  THEN
             hw_old(:) = biomass(:,j,iheartabove,icarbon) + biomass(:,j,iheartbelow,icarbon)  
          ENDIF

          !! 9.1 Calculate the rate of sapwood to heartwood conversion 
          !      Calculate the rate of sapwood to heartwood conversion with the time constant ::tau_sap(j) 
          !      and update aboveground and belowground sapwood ::biomass(:,j,isapabove) and ::biomass(:,j,isapbelow)
          !      and heartwood ::biomass(:,j,iheartabove) and ::biomass(:,j,iheartbelow).

          DO k = 1,nelements

             ! Above the ground
             sapconv(:) = biomass(:,j,isapabove,k) * dt / tau_sap(j)
             biomass(:,j,isapabove,k) = biomass(:,j,isapabove,k) - sapconv(:)
             biomass(:,j,iheartabove,k) =  biomass(:,j,iheartabove,k) + sapconv(:)
             
             ! Below the ground
             sapconv(:) = biomass(:,j,isapbelow,k) * dt / tau_sap(j)
             biomass(:,j,isapbelow,k) = biomass(:,j,isapbelow,k) - sapconv(:)
             biomass(:,j,iheartbelow,k) =  biomass(:,j,iheartbelow,k) + sapconv(:)

          END DO

          !! 9.2 If the vegetation is not dynamic, the age of the plant is decreased. 
          !      The updated heartwood, the sum of new heartwood above and new heartwood below after 
          !      converting sapwood to heartwood, is saved as ::hw_new(:) .
          !      Creation of new heartwood decreases the age of the plant with a factor that is determined by: 
          !      old heartwood ::hw_old(:) divided by the new heartwood ::hw_new(:)
!!!qcj++ peatland
          IF ( (.NOT. ok_dgvm) .AND. (.NOT. ok_dgvm_peat) ) THEN
             hw_new(:) = biomass(:,j,iheartabove,icarbon) + biomass(:,j,iheartbelow,icarbon)
             WHERE ( hw_new(:) .GT. min_sechiba )
                age(:,j) = age(:,j) * hw_old(:)/hw_new(:)
             ENDWHERE
          ELSEIF ( (.NOT. ok_dgvm) .AND. ok_dgvm_peat .AND. (.NOT. is_peat(j)) ) THEN
             hw_new(:) = biomass(:,j,iheartabove,icarbon) + biomass(:,j,iheartbelow,icarbon)
             WHERE ( hw_new(:) .GT. min_sechiba )
                age(:,j) = age(:,j) * hw_old(:)/hw_new(:)
             ENDWHERE
          ENDIF

       ENDIF

    ENDDO       ! loop over PFTs


    CALL xios_orchidee_send_field("HERBIVORES",herbivores)
    CALL xios_orchidee_send_field("LEAF_AGE",leaf_meanage)
    

    ! Write mean leaf age and time constant of probability of a leaf to be eaten by a herbivore 
    ! to the stomate output file.
    CALL histwrite_p (hist_id_stomate, 'LEAF_AGE', itime, &
         leaf_meanage, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'HERBIVORES', itime, &
         herbivores, npts*nvm, horipft_index)
    WHERE(senescence)
       histvar=un
    ELSEWHERE
       histvar=zero
    ENDWHERE    
    CALL histwrite_p (hist_id_stomate, 'SENESCENCE', itime, histvar, npts*nvm, horipft_index)
    IF (printlev>=4) WRITE(numout,*) 'Leaving turnover'

  END SUBROUTINE turn

END MODULE stomate_turnover
