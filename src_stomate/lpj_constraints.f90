! =================================================================================================================================
! MODULE       : lpj_constraints
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
!              This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Groups the subroutines that: (1) initialize all variables in 
!! lpj_constraints and (2) check the temperature threshold to decide for each PFT 
!! if it can adapt to and regenerate under prevailing climate conditions**1
!! 
!!\n RECENT CHANGE(S) : None
!!
!! REFERENCE(S) :
!! - Sitch, S., B. Smith, et al. (2003), Evaluation of ecosystem dynamics,
!!   plant geography and terrestrial carbon cycling in the LPJ dynamic 
!!   global vegetation model, Global Change Biology, 9, 161-185.\n
!! - Smith, B., I. C. Prentice, et al. (2001), Representation of vegetation
!!   dynamics in the modelling of terrestrial ecosystems: comparing two
!!   contrasting approaches within European climate space,
!!   Global Ecology and Biogeography, 10, 621-637.\n
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_stomate/lpj_constraints.f90 $
!! $Date: 2017-10-26 15:35:04 +0200 (Thu, 26 Oct 2017) $
!! $Revision: 4719 $
!! \n
!_ ================================================================================================================================

MODULE lpj_constraints

  ! modules used:
  USE xios_orchidee
  USE ioipsl_para
  USE stomate_data
  USE constantes
  USE pft_parameters

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC constraints,constraints_clear

  REAL(r_std), PARAMETER          :: grow_limit = 7+ZeroCelsius          !! Growing-season temperature limit to tree extension (K)
  LOGICAL, SAVE                   :: firstcall_constraints = .TRUE.      !! first call
!$OMP THREADPRIVATE(firstcall_constraints)

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE   : constraints_clear
!!
!>\BRIEF        Set the flag ::firstcall_constraints to .TRUE. and as such activate section 
!!              1.1 of the subroutine constraints (see subroutine constraints).
!!
!_ ================================================================================================================================

  SUBROUTINE constraints_clear
    firstcall_constraints = .TRUE. 
  END SUBROUTINE constraints_clear


!! ================================================================================================================================
!! SUBROUTINE   : constraints
!!
!>\BRIEF        Determine whether each PFT can adapt to and regenerate under the prevailing climate
!! conditions. Climate conditions are characterised by different threshold values for air temperature.
!!
!! DESCRIPTION : PFTs are adapted to the climate conditions if the daily air temperature does not drop below 
!! the treshold values ::tcm_crit. Some PFT's do not have a ::tcm_crit treshold. Seasonal trees die if leafage 
!! does not show a clear seasonality. (i.e. if the start of the growing season is never detected)
!! 
!! If the monthly temperature is below ::tcm_crit i.e. the critical temperature of the coldest month, the 
!! PFT will be able to regenerate. If minimum temperatures do not drop below ::tcm_crit, its regenerative 
!! capacity decreases with time. Hence, plants that need vernalization die after a few years if they don't
!! vernalize (even if they would not loose their leaves).
!!
!! The treshold values ::t_min_crit, the critical temperature of the coldest month is defined in 
!! stomate_constants.f90'. However, ::regenerate_min, the critical temperature to support regeneration is 
!! calculated in this routine from parameters of which none depenent on the PFT. The value for 
!! ::large_value is defined as 1.E33 in stomata_constraints
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : ::adapted (0-1, unitless) and ::regenerate (0-1, unitless)
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART : 
!! \latexonly 
!!     \includegraphics[scale=0.3]{lpj_constraints_flowchart.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================
  
    SUBROUTINE constraints (npts, dt, &
       t2m_month, t2m_min_daily, when_growthinit, &
       adapted, regenerate, Tseason, &
!qcj++ peatland
       wtpmax_month,Kwt)
    
    !! 0. Variable and parameter declaration
    
    !! 0.1 Input variable

    INTEGER(i_std), INTENT(in)                      :: npts            !! Domain size (unitless)
    REAL(r_std), INTENT(in)                         :: dt              !! Time step   (days)
    REAL(r_std), DIMENSION(npts), INTENT(in)        :: t2m_month       !! "Monthly" 2-meter temperature by defualt 
                                                                       !! monthly spans 20 days a (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)        :: Tseason         !! "seasonal" 2-meter temperature (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)        :: t2m_min_daily   !! Daily minimum 2-meter temperature (K)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)    :: when_growthinit !! Days since beginning of growing season (days)
!qcj++ peatland
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)    :: wtpmax_month
    REAL(r_std), DIMENSION(npts,nvm),INTENT(out)    :: Kwt
    !! 0.2 Output variables

    !! 0.3 Modified variables
  
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout) :: adapted         !! plant adaptation to climate: progressively decrease if temperature is too low
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout) :: regenerate      !! plant regeneration: if vernalization criteria is not reached (i.e not allow
                                                                       !! plant to correctly seed and reproduce) then regeneration progressively decrease

    !! 0.4 Local variables
    REAL(r_std)                                     :: tau_adapt       !! Memory length for adaption (days)
    REAL(r_std)                                     :: tau_regenerate  !! Memory length for regeneration (days)
    REAL(r_std)                                     :: regenerate_min  !! Critical value of "regenerate" below which plant 
                                                                       !! dies (unitless)
    INTEGER(i_std)                                  :: i,j               !! Index
!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering constraints' ! diagnostic level in implimentation

  !! 1. Initializations

    tau_adapt = one_year
    tau_regenerate = one_year

    !! 1.1 Print parameter settings
    IF ( firstcall_constraints ) THEN
       IF (printlev >= 2) THEN
          WRITE(numout,*) 'constraints:'
       
          WRITE(numout,*) '   > Memory length for adaption (d): ',tau_adapt
          WRITE(numout,*) '   > Memory length for regeneration (d): ',tau_regenerate
          WRITE(numout,*) '   > Longest sustainable time without vernalization (y):', too_long
          WRITE(numout,*) '   > For trees, longest sustainable time without growth init (y):', too_long
       END IF
       firstcall_constraints = .FALSE.
       
    ENDIF

    !! 1.2 Calculate critical value for "regenerate"
    !      Critical value for "regenerate", below this value, the last vernalization
    !      happened too far in the past, The PFT is can not regenerate under the 
    !      prevailing climate conditions.
    regenerate_min = exp ( - too_long * one_year / tau_regenerate )

  !! 2. Calculate ::adapted and ::regenerate
    Kwt(:,:)=un
    DO j = 2,nvm ! Loop over # PFTs

       !! 2.1 PFT mask for natural or agriculture vegetations
       IF ( natural(j) .OR. agriculture ) THEN

          !! 2.1.1 Climate criteria
          !! There is no critical temperature for the PFT 
          !  Frost restistant PFT do not have a critical temperature for growth
          IF ( tmin_crit(j) .EQ. undef ) THEN
             adapted(:,j) = un 
          !!qcj++ peatland: test run shows that tropical PFTs occur in the arctic region
          !!the minimum temperature criterion described in Krinner et al. 2005
          !![WHERE (t2m_min_daily(:) .LT. tmin_crit(j)) adapted(:,j)=0] was removed 
          !!I reintroduce the minimum temperature criterion for tropical trees
          !!while the mortality rate of boreal and temperate trees is increased when t2m_min_daily<tmin_crit in lpj_gap!
          ELSEIF ( is_peat(j) .AND. is_tropical(j) ) THEN
             WHERE ( t2m_min_daily(:) .LT. tmin_crit(j) )
                adapted(:,j) = zero
             ENDWHERE

          ENDIF

          !! 2.1.2 Seasonal trees die if leafage does not show a clear seasonality.
          !        Seasonal trees die if leafage does not show a clear seasonality
          !        (i.e. if the start of the growing season is never detected).
          IF ( is_tree(j) .AND. ( pheno_model(j) .NE. 'none' ) ) THEN

             WHERE ( when_growthinit(:,j) .GT. too_long*one_year .AND. when_growthinit(:,j).LT. large_value)
                adapted(:,j) = zero
             ENDWHERE

          ENDIF

          ! "seasonal" temperature Tseason must exceed gro_limit degree (grow_limit) for trees to be declared adapted.
          IF ( is_tree(j) ) THEN
             WHERE ( Tseason(:) .LT. grow_limit )
                adapted(:,j)=zero
             ENDWHERE
          ENDIF

          !qcj++ peatland, peatland boreal needleleaf(leaf_tab=2) evergreen trees, 
          !and broadleaf deciduous shrubs (leaf_tab=1 and pheno_type=2)
          !constrained by water table position[can not establish if too wet]
          IF ( is_tree(j) .AND. is_peat(j) .AND. leaf_tab(j)==2 ) THEN
            ! WRITE (numout,*) 'QCJ check wtpmax_month,',wtpmax_month(:,j)  
            Kwt(:,j)=(un-zero)/(wtpwet_crit(j)-zero)*(wtpmax_month(:,j)-wtpwet_crit(j))  
            Kwt(:,j)=MAX( zero, MIN( un, Kwt(:,j) ) )
          ENDIF

          IF ( is_tree(j) .AND. is_peat(j) .AND. (leaf_tab(j)==1 .AND. pheno_type(j)==2) ) THEN
            Kwt(:,j)=(un-zero)/(wtpwet_crit(j)-zero)*(wtpmax_month(:,j)-wtpwet_crit(j))
            Kwt(:,j)=MAX( zero, MIN( un, Kwt(:,j) ) )
          ENDIF
  

          !qcj++ peatland, peatland moss and grass constrained by water table position
          ![adataption decreases if too dry]

          IF ( natural(j) .AND. is_peat(j) .AND. (.NOT. is_tree(j)) ) THEN
             Kwt(:,j) = un + (zero-un)/(wtpdry_crit(j)-wtp_crit(j))*(wtpmax_month(:,j)-wtp_crit(j)) 
             Kwt(:,j)=MAX( zero, MIN( un, Kwt(:,j) ) )
          ENDIF

          ! if Tseason was less than grow_limit from previous year then the adataption decrease
          adapted(:,j) = un - ( un - adapted(:,j) ) * (tau_adapt- dt)/tau_adapt 

          !! 2.1.3 Test if PFT is regenerative

          !! 2.1.3.1 Check PFT vernalization. 
          !          If sufficiently cold, PFT will be able to regenerate for some time.
          !          Several PFTs (ex: evergreen) don't need vernalization 
          IF ( tcm_crit(j) .EQ. undef ) THEN

             regenerate(:,j) = un

          !! 2.1.3.2 PFT needs vernalization
          ELSE
             ! if the vernalization is reach, the plant flowering will be correct and then the production of seed
             ! will be correct, making the plant able to reproduce and regenerate (i.e fitness)
             WHERE ( t2m_month(:) .LE. tcm_crit(j) )
                regenerate(:,j) = un
             ENDWHERE

             ! each time vernalization is not reached, the fitness decrease
             !  hence, with time ::regenerate approaches 0
             regenerate(:,j) = regenerate(:,j) * (tau_regenerate-dt)/tau_regenerate

          ENDIF

          !! 2.1.4 
          !        Plants that need vernalization die after a few years if they don't
          !        vernalize (even if they would not loose their leaves).
          WHERE ( regenerate(:,j) .LE. regenerate_min )
             adapted(:,j) = zero
          ENDWHERE

       !! 2.1 PFT except natural and agriculture vegetation
       !      Should be developed if needed 
       ELSE

          adapted(:,j) = zero

          regenerate(:,j) = zero

       ENDIF ! PFT of natural or agriculture

    ENDDO ! Loop over # PFTs

  !! 3. Write history files
    CALL xios_orchidee_send_field("ADAPTATION",adapted)
    CALL xios_orchidee_send_field("REGENERATION",regenerate)

    CALL histwrite_p (hist_id_stomate, 'ADAPTATION', itime, &
         adapted, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'REGENERATION', itime, &
         regenerate, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'T2M_MIN_DAILY', itime, &
         t2m_min_daily, npts, horipft_index)

    CALL xios_orchidee_send_field("T2M_MIN_DAILY",t2m_min_daily)
!!!qcj++ peatland
    CALL xios_orchidee_send_field("KWT_ADAPTATION",Kwt)
    IF (printlev>=4) WRITE(numout,*) 'Leaving constraints'

  END SUBROUTINE constraints

END MODULE lpj_constraints
