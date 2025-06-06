! =================================================================================================================================
! MODULE       : stomate_prescribe
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         Initialize and update density, crown area.
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	:
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_stomate/stomate_prescribe.f90 $
!! $Date: 2017-10-26 15:35:04 +0200 (Thu, 26 Oct 2017) $
!! $Revision: 4719 $
!! \n
!_ ================================================================================================================================

MODULE stomate_prescribe

  ! modules used:

  USE ioipsl_para
  USE stomate_data
  USE pft_parameters
  USE constantes

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC prescribe,prescribe_clear

    ! first call
    LOGICAL, SAVE                                              :: firstcall_prescribe = .TRUE.
!$OMP THREADPRIVATE(firstcall_prescribe)

CONTAINS

! =================================================================================================================================
!! SUBROUTINE   : prescribe_clear
!!
!>\BRIEF        : Set the firstcall_prescribe flag back to .TRUE. to prepare for the next simulation.
!_=================================================================================================================================

  SUBROUTINE prescribe_clear
    firstcall_prescribe=.TRUE.
  END SUBROUTINE prescribe_clear

!! ================================================================================================================================
!! SUBROUTINE   : prescribe
!!
!>\BRIEF         Works only with static vegetation and agricultural PFT. Initialize biomass,
!!               density, crown area in the first call and update them in the following.
!!
!! DESCRIPTION (functional, design, flags): \n
!! This module works only with static vegetation and agricultural PFT.
!! In the first call, initialize density of individuals, biomass, crown area,
!! and leaf age distribution to some reasonable value. In the following calls,
!! these variables are updated.
!!
!! To fulfill these purposes, pipe model are used:
!! \latexonly 
!!     \input{prescribe1.tex}
!!     \input{prescribe2.tex}
!! \endlatexonly
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLES(S): ::ind, ::cn_ind, ::leaf_frac
!!
!! REFERENCES   :
!! - Krinner, G., N. Viovy, et al. (2005). "A dynamic global vegetation model 
!!   for studies of the coupled atmosphere-biosphere system." Global 
!!   Biogeochemical Cycles 19: GB1015, doi:1010.1029/2003GB002199.
!! - Sitch, S., B. Smith, et al. (2003), Evaluation of ecosystem dynamics,
!!   plant geography and terrestrial carbon cycling in the LPJ dynamic 
!!   global vegetation model, Global Change Biology, 9, 161-185.
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

 SUBROUTINE prescribe (npts, &
                        veget_cov_max, dt, PFTpresent, everywhere, when_growthinit, &
                        biomass, leaf_frac, ind, cn_ind, co2_to_bm)

!! 0. Parameters and variables declaration

   !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                                :: npts            !! Domain size (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: veget_cov_max   !! "maximal" coverage fraction of a PFT (LAI -> infinity) on ground (unitless;0-1)
    REAL(r_std), INTENT(in)                                   :: dt              !! time step (dt_days)
    !! 0.2 Output variables 

    !! 0.3 Modified variables

    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)               :: PFTpresent      !! PFT present (0 or 1)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: everywhere      !! is the PFT everywhere in the grid box or very localized (after its introduction) (?)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: when_growthinit !! how many days ago was the beginning of the growing season (days)
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: biomass   !! biomass (gC/(m^2 of ground))
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout) :: leaf_frac       !! fraction of leaves in leaf age class (unitless;0-1)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: ind             !! density of individuals (1/(m^2 of ground))
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: cn_ind          !! crown area per individual (m^2)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: co2_to_bm       !! co2 taken up by carbohydrate 
                                                                                 !! reserve at the beginning of the
                                                                                 !! growing season @tex ($gC m^{-2} 
                                                                                 !! of total ground/day$) @endtex 

    !! 0.4 Local variables

    REAL(r_std), DIMENSION(npts)                              :: dia             !! stem diameter (m)
    REAL(r_std), DIMENSION(npts)                              :: woodmass        !! woodmass (gC/(m^2 of ground))
    REAL(r_std), DIMENSION(npts)                              :: woodmass_ind    !! woodmass of an individual (gC)
    INTEGER(i_std)                                            :: i,j             !! index (unitless)

!_ ================================================================================================================================

 IF (.NOT. ok_dgvm_peat) THEN !!!qcj++ peatland

    DO j = 2,nvm

      ! only when the DGVM is not activated or agricultural PFT.

      IF ( ( .NOT. ok_dgvm .AND. lpj_gap_const_mort ) .OR. ( .NOT. natural(j) .OR. pasture(j)) ) THEN

        !
        !! 1.Update crown area
        !

        cn_ind(:,j) = zero

        IF ( is_tree(j) ) THEN

          !
          !! 1.1 treat for trees
          !

          dia(:) = zero

          DO i = 1, npts ! loop over grid points

            IF ( veget_cov_max(i,j) .GT. zero ) THEN

              !! 1.1.1 calculate wood mass on an area basis, which include sapwood and heartwood aboveground and belowground.

              woodmass(i) = (biomass(i,j,isapabove,icarbon) + biomass(i,j,isapbelow,icarbon) + &
                   biomass(i,j,iheartabove,icarbon) + biomass(i,j,iheartbelow,icarbon)) * veget_cov_max(i,j)  

              IF ( woodmass(i) .GT. min_stomate ) THEN

                !! 1.1.2 calculate critical individual density
!?? the logic for 1.1.3 and 1.1.2 is strange, it should be the case that first to calculate critical woodmass per individual,
!?? then calculate critical density.


                ! how to derive the following equation:
                ! first, TreeHeight=pipe_tune2 * Diameter^{pipe_tune3}
                ! we assume the tree is an ideal cylinder, so it volume is: Volume = pi*(Dia/2)^2*Height = pi/4 * Dia * pipe_tune2*Dia^{pipe_tune3} 
                !                                                                  = pi/4 * pipe_tune2 * Dia^{2+pipe_tune3}
                ! last, the woodmass_per_individual = pipe_density * Volume = pipe_density*pi/4.*pipe_tune2 * Dia^{2+pipe_tune3}             
                ind(i,j) = woodmass(i) / &
                           ( pipe_density*pi/4.*pipe_tune2 * maxdia(j)**(2.+pipe_tune3) )

                !! 1.1.3 individual biomass corresponding to this critical density of individuals

                woodmass_ind(i) = woodmass(i) / ind(i,j)

                !! 1.1.4 calculate stem diameter per individual tree

                dia(i) = ( woodmass_ind(i) / ( pipe_density * pi/4. * pipe_tune2 ) ) ** &
                         ( un / ( 2. + pipe_tune3 ) )

                !! 1.1.5 calculate provisional tree crown area for per individual tree

                ! equation: CrownArea=pipe_tune1 * Diameter^{1.6}
                cn_ind(i,j) = pipe_tune1 * MIN( maxdia(j), dia(i) ) ** pipe_tune_exp_coeff

                !! 1.1.6 If total tree crown area for this tree PFT exceeds its veget_cov_max, tree density is recalculated.

                IF ( cn_ind(i,j) * ind(i,j) .GT. 1.002* veget_cov_max(i,j) ) THEN

                  ind(i,j) = veget_cov_max(i,j) / cn_ind(i,j)

                ELSE

                   ind(i,j) = ( veget_cov_max(i,j) / &
                        &     ( pipe_tune1 * (woodmass(i)/(pipe_density*pi/4.*pipe_tune2)) &
                        &     **(pipe_tune_exp_coeff/(2.+pipe_tune3)) ) ) &
                        &     ** (1./(1.-(pipe_tune_exp_coeff/(2.+pipe_tune3))))
                   

                  woodmass_ind(i) = woodmass(i) / ind(i,j)

                  dia(i) = ( woodmass_ind(i) / ( pipe_density * pi/4. * pipe_tune2 ) ) ** &
                           ( un / ( 2. + pipe_tune3 ) )

                  ! final crown area
                  cn_ind(i,j) = pipe_tune1 * MIN( maxdia(j), dia(i) ) ** pipe_tune_exp_coeff

                ENDIF

              ELSE !woodmas=0  => impose some value

                dia(:) = maxdia(j)

                cn_ind(i,j) = pipe_tune1 * MIN( maxdia(j), dia(i) ) ** pipe_tune_exp_coeff

              ENDIF ! IF ( woodmass(i) .GT. min_stomate )

            ENDIF    ! veget_cov_max .GT. 0.

          ENDDO      ! loop over grid points

        ELSEIF ( is_mosspeat(j) ) THEN  !!!qcj++ peatland,following Druel et al. 2017 GMD paper

          WHERE ( veget_cov_max(:,j) .GT. zero )
             cn_ind(:,j) = un
          ENDWHERE

        ELSE !grass

          !
          !! 1.2 grasses: crown area always set to 1m**2
          !

          WHERE ( veget_cov_max(:,j) .GT. zero )
            cn_ind(:,j) = un
          ENDWHERE

        ENDIF   !IF ( is_tree(j) )

        !
        !! 2 density of individuals
        !
        
        WHERE ( veget_cov_max(:,j) .GT. zero )

          ind(:,j) = veget_cov_max(:,j) / cn_ind(:,j)  

        ELSEWHERE

          ind(:,j) = zero

        ENDWHERE

      ENDIF     ! IF ( ( .NOT. ok_dgvm .AND. lpj_gap_const_mort ) .OR. ( .NOT. natural(j) ) )

    ENDDO       ! loop over PFTs

    !
    !!? it's better to move the code for first call at the beginning of the module.
    !! 2 If it's the first call for this module, 
    !

    IF (( firstcall_prescribe ) .AND. (TRIM(stom_restname_in) == 'NONE')) THEN

       IF (printlev >= 2) THEN
          WRITE(numout,*) 'prescribe:'
          ! impose some biomass if zero and PFT prescribed
          WRITE(numout,*) '   > Imposing initial biomass for prescribed trees, '// &
               'initial reserve mass for prescribed grasses.'
          WRITE(numout,*) '   > Declaring prescribed PFTs present.'
       ENDIF

      DO j = 2,nvm ! loop over PFTs
        DO i = 1, npts ! loop over grid points

          ! is vegetation static or PFT agricultural?
          ! Static vegetation or agricultural PFT
          IF ( ( .NOT. ok_dgvm ) .OR. &
               ( ( .NOT. natural(j) .OR. pasture(j)) .AND. ( veget_cov_max(i,j) .GT. min_stomate ) ) ) THEN

            !
            !! 2.1 if tree biomass is extremely small, prescribe the biomass by assuming they have sapling biomass, which is a constant in the model.
            !!     then set all the leaf age as 1.
            !
            ! if tree PFT and biomass too small, prescribe the biomass to a value.
            IF ( is_tree(j) .AND. &
                 ( veget_cov_max(i,j) .GT. min_stomate ) .AND. &
                 ( SUM( biomass(i,j,:,icarbon) ) .LE. min_stomate ) ) THEN
               !!? here the code is redundant, as "veget_cov_max(i,j) .GT. min_stomate" is already met in the above if condition.
               IF (veget_cov_max(i,j) .GT. min_stomate) THEN
                  biomass(i,j,:,:) = (bm_sapl_rescale * bm_sapl(j,:,:) * ind(i,j)) / veget_cov_max(i,j)
               ELSE
                  biomass(i,j,:,:) = zero
               ENDIF

              ! set leaf age classes
              leaf_frac(i,j,:) = zero
              leaf_frac(i,j,1) = un

              ! set time since last beginning of growing season
              when_growthinit(i,j) = large_value

              ! seasonal trees: no leaves at beginning

              IF ( pheno_model(j) .NE. 'none' ) THEN

                biomass(i,j,ileaf,icarbon) = zero
                leaf_frac(i,j,1) = zero

              ENDIF

              co2_to_bm(i,j) = co2_to_bm(i,j) + ( SUM(biomass(i,j,:,icarbon))  / dt )
            ENDIF

            !
            !! 2.2 for grasses, set only the carbon reserve pool to "sapling" carbon reserve pool.
            !!     and set all leaf age to 1.

            IF ( ( .NOT. is_tree(j) ) .AND. &
                 ( veget_cov_max(i,j) .GT. min_stomate ) .AND. &
                 ( SUM( biomass(i,j,:,icarbon) ) .LE. min_stomate ) ) THEN

              biomass(i,j,icarbres,:) = bm_sapl(j,icarbres,:) * ind(i,j) / veget_cov_max(i,j)

              IF ( is_mosspeat(j) ) THEN   !!!qcj++ peatland, following Druel et al. 2017 GMD paper
                 biomass(i,j,ileaf,:) = bm_sapl(j,ileaf,:) * ind(i,j) /veget_cov_max(i,j)
              ENDIF


              ! set leaf age classes
              leaf_frac(i,j,:) = zero
              leaf_frac(i,j,1) = un

              ! set time since last beginning of growing season
              when_growthinit(i,j) = large_value
        
              co2_to_bm(i,j) = co2_to_bm(i,j) + ( biomass(i,j,icarbres,icarbon)  / dt )
              IF ( is_mosspeat(j) ) THEN   !!!qcj++ peatland, following Druel et al. 2017 GMD paper
                  co2_to_bm(i,j) = co2_to_bm(i,j) + ( biomass(i,j,ileaf,icarbon)  / dt )
              ENDIF
            ENDIF

            !
            !! 2.3 declare all PFTs with positive veget_cov_max as present everywhere in that grid box
            !

            IF ( veget_cov_max(i,j) .GT. min_stomate ) THEN
              PFTpresent(i,j) = .TRUE.
              everywhere(i,j) = un
            ENDIF

          ENDIF   ! not ok_dgvm  or agricultural

        ENDDO ! loop over grid points
      ENDDO ! loop over PFTs

    ENDIF !IF (( firstcall_prescribe ) .AND. (TRIM(stom_restname_in) == 'NONE')
  
 ELSE !IF ok_dgvm_peat

   IF (ok_dgvm) THEN
       CALL ipslerr_p(3,'prescribe','Coherence error','ok_dgvm_peat and ok_dgvm can not be true at the same time.','')
   ENDIF
   
   DO j = 2,nvm
      IF ( .NOT. is_peat(j) ) THEN

        !
        !! 1.Update crown area
        !

        cn_ind(:,j) = zero
        IF ( is_tree(j) ) THEN
          !! 1.1 treat for trees
          dia(:) = zero
          DO i = 1, npts ! loop over grid points
            IF ( veget_cov_max(i,j) .GT. zero ) THEN
              !! 1.1.1 calculate wood mass on an area basis, which include sapwood and heartwood aboveground and belowground.
              woodmass(i) = (biomass(i,j,isapabove,icarbon) + biomass(i,j,isapbelow,icarbon) + &
                   biomass(i,j,iheartabove,icarbon) + biomass(i,j,iheartbelow,icarbon)) * veget_cov_max(i,j)

              IF ( woodmass(i) .GT. min_stomate ) THEN
                !! 1.1.2 calculate critical individual density
                ! how to derive the following equation:
                ! first, TreeHeight=pipe_tune2 * Diameter^{pipe_tune3}
                ! we assume the tree is an ideal cylinder, so it volume is: Volume = pi*(Dia/2)^2*Height = pi/4 * Dia * pipe_tune2*Dia^{pipe_tune3}
                !                                                                  = pi/4 * pipe_tune2 * Dia^{2+pipe_tune3}
                ! last, the woodmass_per_individual = pipe_density * Volume = pipe_density*pi/4.*pipe_tune2 * Dia^{2+pipe_tune3}
                ind(i,j) = woodmass(i) / &
                           ( pipe_density*pi/4.*pipe_tune2 * maxdia(j)**(2.+pipe_tune3) )

                !! 1.1.3 individual biomass corresponding to this critical density of individuals

                woodmass_ind(i) = woodmass(i) / ind(i,j)
                !! 1.1.4 calculate stem diameter per individual tree

                dia(i) = ( woodmass_ind(i) / ( pipe_density * pi/4. * pipe_tune2 ) ) ** &
                         ( un / ( 2. + pipe_tune3 ) )

                !! 1.1.5 calculate provisional tree crown area for per individual tree

                ! equation: CrownArea=pipe_tune1 * Diameter^{1.6}
                cn_ind(i,j) = pipe_tune1 * MIN( maxdia(j), dia(i) ) ** pipe_tune_exp_coeff
                !! 1.1.6 If total tree crown area for this tree PFT exceeds its veget_cov_max, tree density is recalculated.
                IF ( cn_ind(i,j) * ind(i,j) .GT. 1.002* veget_cov_max(i,j) ) THEN

                  ind(i,j) = veget_cov_max(i,j) / cn_ind(i,j)

                ELSE

                   ind(i,j) = ( veget_cov_max(i,j) / &
                        &     ( pipe_tune1 * (woodmass(i)/(pipe_density*pi/4.*pipe_tune2)) &
                        &     **(pipe_tune_exp_coeff/(2.+pipe_tune3)) ) ) &
                        &     ** (1./(1.-(pipe_tune_exp_coeff/(2.+pipe_tune3))))


                  woodmass_ind(i) = woodmass(i) / ind(i,j)

                  dia(i) = ( woodmass_ind(i) / ( pipe_density * pi/4. * pipe_tune2 ) ) ** &
                           ( un / ( 2. + pipe_tune3 ) )

                  ! final crown area
                  cn_ind(i,j) = pipe_tune1 * MIN( maxdia(j), dia(i) ) ** pipe_tune_exp_coeff

                ENDIF

              ELSE !woodmas=0  => impose some value

                dia(:) = maxdia(j)

                cn_ind(i,j) = pipe_tune1 * MIN( maxdia(j), dia(i) ) ** pipe_tune_exp_coeff

              ENDIF ! IF ( woodmass(i) .GT. min_stomate )
            ENDIF    ! veget_cov_max .GT. 0.

          ENDDO      ! loop over grid points
        ELSE !grass

          !
          !! 1.2 grasses: crown area always set to 1m**2
          !

          WHERE ( veget_cov_max(:,j) .GT. zero )
            cn_ind(:,j) = un
          ENDWHERE

        ENDIF   !IF ( is_tree(j) )

        WHERE ( veget_cov_max(:,j) .GT. zero )

          ind(:,j) = veget_cov_max(:,j) / cn_ind(:,j)

        ELSEWHERE

          ind(:,j) = zero

        ENDWHERE

      ENDIF     ! IF (.NOT. is_peat(j)) 

   ENDDO       ! loop over PFTs

    !! 2 If it's the first call for this module,
    !

    IF (( firstcall_prescribe ) .AND. (TRIM(stom_restname_in) == 'NONE')) THEN

       IF (printlev >= 2) THEN
          WRITE(numout,*) 'prescribe:'
          ! impose some biomass if zero and PFT prescribed
          WRITE(numout,*) '   > Imposing initial biomass for prescribed trees, '// &
               'initial reserve mass for prescribed grasses.'
          WRITE(numout,*) '   > Declaring prescribed PFTs present.'
       ENDIF

      DO j = 2,nvm ! loop over PFTs
        DO i = 1, npts ! loop over grid points

         IF ( .NOT. is_peat(j) ) THEN

            !! 2.1 if tree biomass is extremely small, prescribe the biomass by assuming they have sapling biomass, which is a constant in the model.
            !!     then set all the leaf age as 1.
            !
            ! if tree PFT and biomass too small, prescribe the biomass to a value.
            IF ( is_tree(j) .AND. &
                 ( veget_cov_max(i,j) .GT. min_stomate ) .AND. &
                 ( SUM( biomass(i,j,:,icarbon) ) .LE. min_stomate ) ) THEN
               !!? here the code is redundant, as "veget_cov_max(i,j) .GT. min_stomate" is already met in the above if condition.
               IF (veget_cov_max(i,j) .GT. min_stomate) THEN
                  biomass(i,j,:,:) = (bm_sapl_rescale * bm_sapl(j,:,:) * ind(i,j)) / veget_cov_max(i,j)
               ELSE
                  biomass(i,j,:,:) = zero
               ENDIF

              ! set leaf age classes
              leaf_frac(i,j,:) = zero
              leaf_frac(i,j,1) = un

              ! set time since last beginning of growing season
              when_growthinit(i,j) = large_value

              ! seasonal trees: no leaves at beginning

              IF ( pheno_model(j) .NE. 'none' ) THEN

                biomass(i,j,ileaf,icarbon) = zero
                leaf_frac(i,j,1) = zero

              ENDIF

              co2_to_bm(i,j) = co2_to_bm(i,j) + ( SUM(biomass(i,j,:,icarbon))  / dt )
            ENDIF

            !! 2.2 for grasses, set only the carbon reserve pool to "sapling" carbon reserve pool.
            !!     and set all leaf age to 1.

            IF ( ( .NOT. is_tree(j) ) .AND. &
                 ( veget_cov_max(i,j) .GT. min_stomate ) .AND. &
                 ( SUM( biomass(i,j,:,icarbon) ) .LE. min_stomate ) ) THEN

              biomass(i,j,icarbres,:) = bm_sapl(j,icarbres,:) * ind(i,j) / veget_cov_max(i,j)

              ! set leaf age classes
              leaf_frac(i,j,:) = zero
              leaf_frac(i,j,1) = un

              ! set time since last beginning of growing season
              when_growthinit(i,j) = large_value

              co2_to_bm(i,j) = co2_to_bm(i,j) + ( biomass(i,j,icarbres,icarbon)  / dt )
            ENDIF

            !
            !! 2.3 declare all PFTs with positive veget_cov_max as present everywhere in that grid box
            !

            IF ( veget_cov_max(i,j) .GT. min_stomate ) THEN
              PFTpresent(i,j) = .TRUE.
              everywhere(i,j) = un
            ENDIF

         ENDIF   ! not is_peat

        ENDDO ! loop over grid points
      ENDDO ! loop over PFTs

    ENDIF !IF (( firstcall_prescribe )

 ENDIF  ! ok_dgvm_peat

   firstcall_prescribe = .FALSE.

  END SUBROUTINE prescribe

END MODULE stomate_prescribe
