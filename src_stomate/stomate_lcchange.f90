! =================================================================================================================================
! MODULE       : stomate_lcchange
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Impact of land cover change on carbon stocks
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): Including permafrost carbon
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_stomate/stomate_lcchange.f90 $
!! $Date: 2017-10-26 15:35:04 +0200 (Thu, 26 Oct 2017) $
!! $Revision: 4719 $
!! \n
!_ ================================================================================================================================


MODULE stomate_lcchange

  ! modules used:
  
  USE ioipsl_para
  USE stomate_data
  USE pft_parameters
  USE constantes
  USE constantes_soil_var
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC lcchange_main
  PUBLIC lcchange_deffire
  PUBLIC lcchange_main_agripeat
  PUBLIC agripeat_adjust_fractions
CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : lcchange_main
!!
!>\BRIEF        Impact of land cover change on carbon stocks
!!
!! DESCRIPTION  : This subroutine is always activate if VEGET_UPDATE>0Y in the configuration file, which means that the 
!! vegetation map is updated regulary. lcchange_main is called from stomateLpj the first time step after the vegetation 
!! map has been changed. 
!! The impact of land cover change on carbon stocks is computed in this subroutine. The land cover change is written
!! by the difference of current and previous "maximal" coverage fraction of a PFT. 
!! On the basis of this difference, the amount of 'new establishment'/'biomass export',
!! and increase/decrease of each component, are estimated.\n
!!
!! Main structure of lpj_establish.f90 is:
!! 1. Initialization
!! 2. Calculation of changes in carbon stocks and biomass by land cover change
!! 3. Update 10 year- and 100 year-turnover pool contents
!! 4. History
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : ::prod10, ::prod100, ::flux10, ::flux100,
!!   :: cflux_prod10 and :: cflux_prod100 
!!
!! REFERENCES   : None
!!
!! FLOWCHART    : 
!! \latexonly 
!!     \includegraphics[scale=0.5]{lcchange.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE lcchange_main ( npts, dt_days, veget_cov_max_old, veget_cov_max_new, &
       biomass, ind, age, PFTpresent, senescence, when_growthinit, everywhere, &        
       co2_to_bm, bm_to_litter, turnover_daily, bm_sapl, cn_ind,flux10,flux100, &
       prod10, prod100, &
       convflux, &
       cflux_prod10, cflux_prod100, leaf_frac,&
       npp_longterm, lm_lastyearmax, litter, litter_avail, litter_not_avail, &
       carbon, &
       deepC_a, deepC_s, deepC_p, &
       fuel_1hr, fuel_10hr, fuel_100hr, fuel_1000hr, &
!!!qcj++ peatland
       veget_cov_max_new_peatdgvm,lalo)

    
    IMPLICIT NONE
    
  !! 0. Variable and parameter declaration 
    
    !! 0.1 Input variables
    
    INTEGER, INTENT(in)                                       :: npts             !! Domain size - number of pixels (unitless)
    REAL(r_std), INTENT(in)                                   :: dt_days          !! Time step of vegetation dynamics for stomate
                                                                                  !! (days)
    REAL(r_std), DIMENSION(nvm, nparts,nelements), INTENT(in) :: bm_sapl          !! biomass of sapling 
                                                                                  !! @tex ($gC individual^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: veget_cov_max_old!! Current "maximal" coverage fraction of a PFT (LAI
                                                                                  !! -> infinity) on ground
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: veget_cov_max_new!! New "maximal" coverage fraction of a PFT (LAI ->
                                                                                  !! infinity) on ground (unitless) 
 
!!!qcj++ peatland
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)             :: veget_cov_max_new_peatdgvm
    REAL(r_std),DIMENSION(npts,2),INTENT(in)                    :: lalo
    REAL(r_std), DIMENSION(npts)                              :: sumnat
    REAL(r_std), DIMENSION(npts)                              :: sumcrop
    REAL(r_std), DIMENSION(npts)                              :: sumcroppeat
    REAL(r_std), DIMENSION(npts)                              :: sumnat_adjust
    REAL(r_std), DIMENSION(npts)                              :: sumgridPFT
    REAL(r_std), DIMENSION(npts)                              :: sum_peat_old
    REAL(r_std), DIMENSION(npts)                              :: sum_peat_new
    REAL(r_std), DIMENSION(npts)                              :: sum_peatBare
    REAL(r_std), DIMENSION(npts)                              :: sum_peatBare_adjust
    REAL(r_std), DIMENSION(npts)                              :: delta_peat
    REAL(r_std), DIMENSION(npts)                              :: sumnatpeat

    !! 0.2 Output variables

    REAL(r_std), DIMENSION(:,:), INTENT(out)                  :: convflux         !! release during first year following land cover
                                                                                  !! change  (npts, nwp)
    REAL(r_std), DIMENSION(:,:), INTENT(out)                  :: cflux_prod10     !! total annual release from the 10 year-turnover
                                                                                  !! pool @tex ($gC m^{-2}$) @endtex
                                                                                  !! change  (npts, nwp)
    REAL(r_std), DIMENSION(:,:), INTENT(out)                  :: cflux_prod100    !! total annual release from the 100 year-
                                                                                  !! turnover pool @tex ($gC m^{-2}$) @endtex
                                                                                  !! change  (npts, nwp)
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout):: turnover_daily   !! Turnover rates 
 

    !! 0.3 Modified variables   
    
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout):: biomass    !! biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: ind              !! Number of individuals @tex ($m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: age              !! mean age (years)
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)               :: senescence       !! plant senescent (only for deciduous trees) Set
                                                                                  !! to .FALSE. if PFT is introduced or killed
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)               :: PFTpresent       !! Is pft there (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: everywhere       !! is the PFT everywhere in the grid box or very 
                                                                                  !! localized (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: when_growthinit  !! how many days ago was the beginning of the 
                                                                                  !! growing season (days)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: co2_to_bm        !! biomass uptaken 
                                                                                  !! @tex ($gC m^{-2} day^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: bm_to_litter !! conversion of biomass to litter 
                                                                                  !! @tex ($gC m^{-2} day^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: cn_ind           !! crown area of individuals 
                                                                                  !! @tex ($m^{2}$) @endtex
    REAL(r_std), DIMENSION(npts,0:10,nwp), INTENT(inout)              :: prod10           !! products remaining in the 10 year-turnover
                                                                                  !! pool after the annual release for each 
                                                                                  !! compartment (10 + 1 : input from year of land
                                                                                  !! cover change)  (npts,0:10,nwp)
    REAL(r_std), DIMENSION(npts,0:100,nwp), INTENT(inout)              :: prod100          !! products remaining in the 100 year-turnover
                                                                                  !! pool after the annual release for each 
                                                                                  !! compartment (100 + 1 : input from year of land
                                                                                  !! cover change) (npts,0:100,nwp)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)              :: flux10           !! annual release from the 10/100 year-turnover 
                                                                                  !! pool compartments (npts,10,nwp)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)              :: flux100          !! annual release from the 100/100 year-turnover
                                                                                  !! pool compartments (npts,100,nwp)
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout) :: leaf_frac        !! fraction of leaves in leaf age class 
                                                                                  !! (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: lm_lastyearmax   !! last year's maximum leaf mass for each PFT 
                                                                                  !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: npp_longterm     !! "long term" net primary productivity 
                                                                                  !! @tex ($gC m^{-2} year^{-1}$) @endtex
    REAL(r_std),DIMENSION(npts,nlitt,nvm,nlevs,nelements), INTENT(inout):: litter !! metabolic and structural litter, above and 
                                                                                  !! below ground @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nlitt,nvm), INTENT(inout):: litter_avail
    REAL(r_std), DIMENSION(npts,nlitt,nvm) , INTENT(inout):: litter_not_avail
    REAL(r_std),DIMENSION(npts,ncarb,nvm), INTENT(inout)      :: carbon           !! carbon pool: active, slow, or passive 
    REAL(r_std), DIMENSION(npts,ndeep,nvm), INTENT(inout)         :: deepC_a      !! Permafrost soil carbon (g/m**3) active
    REAL(r_std), DIMENSION(npts,ndeep,nvm), INTENT(inout)         :: deepC_s      !! Permafrost soil carbon (g/m**3) slow
    REAL(r_std), DIMENSION(npts,ndeep,nvm), INTENT(inout)         :: deepC_p      !! Permafrost soil carbon (g/m**3) passive
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements),INTENT(inout)        :: fuel_1hr
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements),INTENT(inout)        :: fuel_10hr
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements),INTENT(inout)        :: fuel_100hr
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements),INTENT(inout)        :: fuel_1000hr
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements)       :: fuel_all_type
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements,4)     :: fuel_type_frac

    !! 0.4 Local variables

    INTEGER(i_std)                                            :: i, j, k, l, m    !! indices (unitless)
    REAL(r_std),DIMENSION(npts,nelements)                     :: bm_new           !! biomass increase @tex ($gC m^{-2}$) @endtex
    REAL(r_std),DIMENSION(npts,nparts,nelements)              :: biomass_loss     !! biomass loss @tex ($gC m^{-2}$) @endtex
    REAL(r_std)                                               :: above            !! aboveground biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std),DIMENSION(npts,nlitt,nlevs,nelements)         :: dilu_lit         !! Litter dilution @tex ($gC m^{-2}$) @endtex
    REAL(r_std),DIMENSION(npts,ncarb)                         :: dilu_soil_carbon !! Soil Carbondilution @tex ($gC m^{-2}$) @endtex
    REAL(r_std),DIMENSION(npts,ndeep,ncarb)                   :: dilu_soil_carbon_vertres !!vertically-resolved Soil Carbondilution (gC/m²)

    REAL(r_std),DIMENSION(nvm)                                :: delta_veg        !! changes in "maximal" coverage fraction of PFT 
    REAL(r_std)                                                 :: delta_veg_sum    !! sum of delta_veg
    REAL(r_std),DIMENSION(npts,nvm)                             :: delta_ind        !! change in number of individuals  
    REAL(r_std), DIMENSION(npts,3)                              :: C_start
    REAL(r_std), DIMENSION(npts,3)                              :: C_end
    REAL(r_std), DIMENSION(npts,3)                              :: C_check
    REAL(r_std), DIMENSION(npts,ndeep,3)                        :: deepC_start
    REAL(r_std), DIMENSION(npts,ndeep,3)                        :: deepC_end
    REAL(r_std), DIMENSION(npts,ndeep,3)                        :: deepC_check
    REAL(r_std),DIMENSION(npts,4)                               :: Fin_start
    REAL(r_std),DIMENSION(npts,4)                               :: Fin_end
    REAL(r_std),DIMENSION(npts,4)                               :: Fout_start
    REAL(r_std),DIMENSION(npts,4)                               :: Fout_end
    REAL(r_std),DIMENSION(npts,4)                               :: Fin_check
    REAL(r_std),DIMENSION(npts,4)                               :: Fout_check


!_ ================================================================================================================================

  IF (.NOT. ok_dgvm_peat) THEN !!!qcj++ peatland

    IF (printlev>=3) WRITE(numout,*) 'Entering lcchange_main'
    
  !! 1. initialization
    
    prod10(:,0,:)         = zero
    prod100(:,0,:)        = zero   
    above               = zero
    convflux(:,:)         = zero
    cflux_prod10(:,:)     = zero
    cflux_prod100(:,:)    = zero
    delta_ind(:,:)      = zero
    delta_veg(:)        = zero
    dilu_soil_carbon_vertres(:,:,:) = zero

  !! 3. calculation of changes in carbon stocks and biomass by land cover change\n
    
    DO i = 1, npts ! Loop over # pixels - domain size
       
       !! 3.1 initialization of carbon stocks\n
       delta_veg(:) = veget_cov_max_new(i,:)-veget_cov_max_old(i,:)
       delta_veg_sum = SUM(delta_veg,MASK=delta_veg.LT.0.)
    
       dilu_lit(i,:,:,:) = zero
       dilu_soil_carbon(i,:) = zero
       biomass_loss(i,:,:) = zero
       
       !! 3.2 if vegetation coverage decreases, compute dilution of litter, soil carbon, and biomass.\n
       DO j=2, nvm
          IF ( delta_veg(j) < -min_stomate ) THEN 
             dilu_lit(i,:,:,:) = dilu_lit(i,:,:,:) + delta_veg(j)*litter(i,:,j,:,:) / delta_veg_sum
             biomass_loss(i,:,:) = biomass_loss(i,:,:) + biomass(i,j,:,:)*delta_veg(j) / delta_veg_sum
             IF ( ok_pc ) THEN
                    dilu_soil_carbon_vertres(i,:,iactive)=dilu_soil_carbon_vertres(i,:,iactive) + &
                         delta_veg(j) * deepC_a(i,:,j) / delta_veg_sum
                    dilu_soil_carbon_vertres(i,:,islow)=dilu_soil_carbon_vertres(i,:,islow) + &
                         delta_veg(j) * deepC_s(i,:,j) / delta_veg_sum
                    dilu_soil_carbon_vertres(i,:,ipassive)=dilu_soil_carbon_vertres(i,:,ipassive) + &
                         delta_veg(j) * deepC_p(i,:,j) / delta_veg_sum
             ELSE
                    dilu_soil_carbon(i,:) =  dilu_soil_carbon(i,:) + delta_veg(j) * carbon(i,:,j) / delta_veg_sum
             ENDIF
          ENDIF
       ENDDO
       
       !! 3.3 
       DO j=2, nvm ! Loop over # PFTs

          !! 3.3.1 The case that vegetation coverage of PFTj increases
          IF ( delta_veg(j) > min_stomate) THEN

             !! 3.3.1.1 Initial setting of new establishment
             IF (veget_cov_max_old(i,j) .LT. min_stomate) THEN 
                IF (is_tree(j)) THEN

                   ! cn_sapl(j)=0.5; stomate_data.f90
                   cn_ind(i,j) = cn_sapl(j) 
                ELSE
                   cn_ind(i,j) = un
                ENDIF
                ind(i,j)= delta_veg(j) / cn_ind(i,j)
                PFTpresent(i,j) = .TRUE.
                everywhere(i,j) = 1.
                senescence(i,j) = .FALSE.
                age(i,j) = zero

                ! large_value = 1.E33_r_std
                when_growthinit(i,j) = large_value 
                leaf_frac(i,j,1) = 1.0
                npp_longterm(i,j) = npp_longterm_init
                lm_lastyearmax(i,j) = bm_sapl(j,ileaf,icarbon) * ind(i,j)
             ENDIF
             IF ( cn_ind(i,j) > min_stomate ) THEN
                delta_ind(i,j) = delta_veg(j) / cn_ind(i,j) 
             ENDIF
             
             !! 3.3.1.2 Update of biomass in each each carbon stock component 
             !!         Update of biomass in each each carbon stock component (leaf, sapabove, sapbelow,
             !>         heartabove, heartbelow, root, fruit, and carbres)\n
             DO k = 1, nparts ! loop over # carbon stock components, nparts = 8; stomate_constant.f90 
                DO l = 1,nelements ! loop over # elements

                   bm_new(i,l) = delta_ind(i,j) * bm_sapl(j,k,l) 
                   IF (veget_cov_max_old(i,j) .GT. min_stomate) THEN

                      ! in the case that bm_new is overestimated compared with biomass?
                      IF ((bm_new(i,l)/delta_veg(j)) > biomass(i,j,k,l)) THEN
                         bm_new(i,l) = biomass(i,j,k,l)*delta_veg(j)
                      ENDIF
                   ENDIF
                   biomass(i,j,k,l) = ( biomass(i,j,k,l) * veget_cov_max_old(i,j) + bm_new(i,l) ) / veget_cov_max_new(i,j)
                   co2_to_bm(i,j) = co2_to_bm(i,j) + (bm_new(i,icarbon)* dt_days) / (one_year * veget_cov_max_new(i,j))
                END DO ! loop over # elements
             ENDDO ! loop over # carbon stock components

             !! 3.3.1.3 Calculation of dilution in litter, soil carbon, and  input of litter
             !!        In this 'IF statement', dilu_* is zero. Formulas for litter and soil carbon
             !!         could be shortend?? Are the following formulas correct?

             ! Litter
             litter(i,:,j,:,:)=(litter(i,:,j,:,:) * veget_cov_max_old(i,j) + &
                  dilu_lit(i,:,:,:) * delta_veg(j)) / veget_cov_max_new(i,j)
                !gmjc available and not available litter for grazing
                ! only not available litter increase/decrease, available litter will not
                ! change, due to tree litter can not be eaten
               IF (is_grassland_manag(j) .AND. is_grassland_grazed(j)) THEN
                 litter_avail(i,:,j) = litter_avail(i,:,j) * veget_cov_max_old(i,j) / veget_cov_max_new(i,j)
                 litter_not_avail(i,:,j) = litter(i,:,j,iabove,icarbon) - litter_avail(i,:,j)
               ENDIF
                !end gmjc            
             IF ( ok_pc ) THEN 
                deepC_a(i,:,j)=(deepC_a(i,:,j) * veget_cov_max_old(i,j) + &
                     dilu_soil_carbon_vertres(i,:,iactive) * delta_veg(j)) / veget_cov_max_new(i,j)
                deepC_s(i,:,j)=(deepC_s(i,:,j) * veget_cov_max_old(i,j) + &
                     dilu_soil_carbon_vertres(i,:,islow) * delta_veg(j)) / veget_cov_max_new(i,j)
                deepC_p(i,:,j)=(deepC_p(i,:,j) * veget_cov_max_old(i,j) + &
                     dilu_soil_carbon_vertres(i,:,ipassive) * delta_veg(j)) / veget_cov_max_new(i,j)
             ELSE
                ! Soil carbon
                carbon(i,:,j)=(carbon(i,:,j) * veget_cov_max_old(i,j) + dilu_soil_carbon(i,:) * delta_veg(j)) / veget_cov_max_new(i,j)
             ENDIF

             DO l = 1,nelements

                ! Litter input
                bm_to_litter(i,j,isapbelow,l) = (bm_to_litter(i,j,isapbelow,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,isapbelow,l)*delta_veg(j)) / veget_cov_max_new(i,j)
                bm_to_litter(i,j,iheartbelow,l) = (bm_to_litter(i,j,iheartbelow,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,iheartbelow,l) *delta_veg(j)) / veget_cov_max_new(i,j)
                bm_to_litter(i,j,iroot,l) = (bm_to_litter(i,j,iroot,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,iroot,l)*delta_veg(j)) / veget_cov_max_new(i,j)
                bm_to_litter(i,j,ifruit,l) = (bm_to_litter(i,j,ifruit,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,ifruit,l)*delta_veg(j)) / veget_cov_max_new(i,j)
                bm_to_litter(i,j,icarbres,l) = (bm_to_litter(i,j,icarbres,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,icarbres,l)   *delta_veg(j)) / veget_cov_max_new(i,j)
                bm_to_litter(i,j,ileaf,l) = (bm_to_litter(i,j,ileaf,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,ileaf,l)*delta_veg(j)) / veget_cov_max_new(i,j)

             END DO

             age(i,j)=age(i,j)*veget_cov_max_old(i,j)/veget_cov_max_new(i,j)
             
          !! 3.3.2 The case that vegetation coverage of PFTj is no change or decreases
          ELSE 
 
             !! 3.3.2.1 Biomass export
             ! coeff_lcchange_*:  Coeff of biomass export for the year, decade, and century
             above = biomass(i,j,isapabove,icarbon) + biomass(i,j,iheartabove,icarbon)
             convflux(i,iwplcc)  = convflux(i,iwplcc)  - ( coeff_lcchange_1(j) * above * delta_veg(j) ) 
             prod10(i,0,iwplcc)  = prod10(i,0,iwplcc)  - ( coeff_lcchange_10(j) * above * delta_veg(j) )
             prod100(i,0,iwplcc) = prod100(i,0,iwplcc) - ( coeff_lcchange_100(j) * above * delta_veg(j) )

             !! If the vegetation is to small, it has been set to 0.
             IF ( veget_cov_max_new(i,j) .LT. min_stomate ) THEN

                ind(i,j) = zero
                biomass(i,j,:,:) = zero
                PFTpresent(i,j) = .FALSE.
                senescence(i,j) = .FALSE.
                age(i,j) = zero
                when_growthinit(i,j) = undef
                everywhere(i,j) = zero
                carbon(i,:,j) = zero
                litter(i,:,j,:,:) = zero
                bm_to_litter(i,j,:,:) = zero

                IF (ok_pc) THEN
                   deepC_a(i,:,j) = zero
                   deepC_s(i,:,j) = zero
                   deepC_p(i,:,j) = zero
                ENDIF

             ENDIF

          ENDIF ! End if PFT's coverage reduction
          
       ENDDO ! Loop over # PFTs
       
       !! 2.4 update 10 year-turnover pool content following flux emission
       !!     (linear decay (10%) of the initial carbon input)
       DO  l = 0, 8
          m = 10 - l
          cflux_prod10(i,iwplcc) =  cflux_prod10(i,iwplcc) + flux10(i,m,iwplcc)
          prod10(i,m,iwplcc)     =  prod10(i,m-1,iwplcc)   - flux10(i,m-1,iwplcc)
          flux10(i,m,iwplcc)     =  flux10(i,m-1,iwplcc)
          
          IF (prod10(i,m,iwplcc) .LT. 1.0) prod10(i,m,iwplcc) = zero
       ENDDO
       
       cflux_prod10(i,iwplcc) = cflux_prod10(i,iwplcc) + flux10(i,1,iwplcc) 
       flux10(i,1,iwplcc)     = 0.1 * prod10(i,0,iwplcc)
       prod10(i,1,iwplcc)     = prod10(i,0,iwplcc)
       
       !! 2.5 update 100 year-turnover pool content following flux emission\n
       DO   l = 0, 98
          m = 100 - l
          cflux_prod100(i,iwplcc)  =  cflux_prod100(i,iwplcc) + flux100(i,m,iwplcc)
          prod100(i,m,iwplcc)      =  prod100(i,m-1,iwplcc)   - flux100(i,m-1,iwplcc)
          flux100(i,m,iwplcc)      =  flux100(i,m-1,iwplcc)
          
          IF (prod100(i,m,iwplcc).LT.1.0) prod100(i,m,iwplcc) = zero
       ENDDO
       
       cflux_prod100(i,iwplcc)  = cflux_prod100(i,iwplcc) + flux100(i,1,iwplcc) 
       flux100(i,1,iwplcc)      = 0.01 * prod100(i,0,iwplcc)
       prod100(i,1,iwplcc)      = prod100(i,0,iwplcc)
       prod10(i,0,iwplcc)        = zero
       prod100(i,0,iwplcc)       = zero 
       


    ENDDO ! Loop over # pixels - domain size
    
    !! We redistribute the updated litter into four fuel classes, so that
    !! the balance between aboveground litter and fuel is mainted. The subtraction
    !! of fuel burned by land cover change fires from the fuel pool is made here. 
    fuel_all_type(:,:,:,:) = fuel_1hr(:,:,:,:) + fuel_10hr(:,:,:,:) + &
                               fuel_100hr(:,:,:,:) + fuel_1000hr(:,:,:,:)
    fuel_type_frac(:,:,:,:,:) = 0.25
    WHERE(fuel_all_type(:,:,:,:) > min_stomate)
      fuel_type_frac(:,:,:,:,1) = fuel_1hr(:,:,:,:)/fuel_all_type(:,:,:,:)
      fuel_type_frac(:,:,:,:,2) = fuel_10hr(:,:,:,:)/fuel_all_type(:,:,:,:)
      fuel_type_frac(:,:,:,:,3) = fuel_100hr(:,:,:,:)/fuel_all_type(:,:,:,:)
      fuel_type_frac(:,:,:,:,4) = fuel_1000hr(:,:,:,:)/fuel_all_type(:,:,:,:)
    ENDWHERE
    DO j=1,nvm
      fuel_1hr(:,j,:,:) = litter(:,:,j,iabove,:) * fuel_type_frac(:,j,:,:,1) 
      fuel_10hr(:,j,:,:) = litter(:,:,j,iabove,:) * fuel_type_frac(:,j,:,:,2) 
      fuel_100hr(:,j,:,:) = litter(:,:,j,iabove,:) * fuel_type_frac(:,j,:,:,3) 
      fuel_1000hr(:,j,:,:) = litter(:,:,j,iabove,:) * fuel_type_frac(:,j,:,:,4) 
    END DO 

  !! 3. history
    convflux        = convflux/one_year*dt_days
    cflux_prod10    = cflux_prod10/one_year*dt_days
    cflux_prod100   = cflux_prod100/one_year*dt_days
    
    IF (printlev>=4) WRITE(numout,*) 'Leaving lcchange_main'

  ELSE  !!!ok_dgvm_peat

    IF (printlev>=3) WRITE(numout,*) 'Entering lcchange_main, ok_dgvm_peat'

  !! 1. initialization

    prod10(:,0,:)         = zero
    prod100(:,0,:)        = zero
    above               = zero
    convflux(:,:)         = zero
    cflux_prod10(:,:)     = zero
    cflux_prod100(:,:)    = zero
    delta_ind(:,:)      = zero
    dilu_soil_carbon_vertres(:,:,:) = zero

    C_start(:,:)=zero
    C_end(:,:)=zero

    deepC_start(:,:,:)=zero
    deepC_end(:,:,:)=zero
    deepC_check(:,:,:)=zero

    Fin_start(:,:)=zero
    Fin_end(:,:)=zero
    Fout_start(:,:)=zero
    Fout_end(:,:)=zero

    C_start(:,1)=zero !SUM(SUM(biomass(:,:,:,icarbon),DIM=3) * veget_cov_max_old,DIM=2)
    C_start(:,2)=SUM(SUM(SUM(litter(:,:,:,:,icarbon),DIM=4),DIM=2) *veget_cov_max_old,DIM=2)
    C_start(:,3)=zero!SUM(SUM(carbon(:,:,:),DIM=2) * veget_cov_max_old,DIM=2)

    Fin_start(:,1)=zero !SUM(SUM(bm_to_litter(:,:,:,icarbon),DIM=3) * veget_cov_max_old,DIM=2)
    Fin_start(:,2)=zero 
    Fin_start(:,3)=zero !SUM(co2_to_bm * veget_cov_max_old,DIM=2)
    Fin_start(:,4)=zero 

    DO l=1,ndeep
      deepC_start(:,l,1)=SUM(deepC_a(:,l,:)*veget_cov_max_old(:,:),DIM=2)
      deepC_start(:,l,2)=SUM(deepC_s(:,l,:)*veget_cov_max_old(:,:),DIM=2)
      deepC_start(:,l,3)=SUM(deepC_p(:,l,:)*veget_cov_max_old(:,:),DIM=2)
    ENDDO
!!! non-peatland
  !! 3. calculation of changes in carbon stocks and biomass by land cover change\n
    DO i = 1, npts ! Loop over # pixels - domain size
      ! IF (lalo(i,1)==27 .AND. lalo(i,2)==23) THEN
      ! WRITE (numout,*) 'QCJ check lalo',lalo(i,1),lalo(i,2) 
      ! WRITE (numout,*) 'QCJ check veget_cov_max_old sum',SUM(veget_cov_max_old(i,:)) 
      ! WRITE (numout,*) 'QCJ check veget_cov_max_new sum, PFTmap',SUM(veget_cov_max_new(i,:))
      ! ENDIF

       veget_cov_max_new_peatdgvm(i,:)=zero

       sumnat(i)=zero
       sumcrop(i)=zero
       sumcroppeat(i)=zero
       sumnat_adjust(i)=zero
       sumnatpeat(i)=zero
       delta_veg(:) = zero
 ! peatland
       DO j=2, nvm
          IF ( natural(j) .AND. .NOT. pasture(j) .AND. ( is_peat(j)) ) THEN
            veget_cov_max_new_peatdgvm(i,j)=veget_cov_max_old(i,j)
            sumnatpeat(i)=sumnatpeat(i)+veget_cov_max_old(i,j)
            !IF (lalo(i,1)==27 .AND. lalo(i,2)==23) THEN
            !WRITE (numout,*) 'QCJ check natural peat,',j,veget_cov_max_old(i,j),veget_cov_max_new(i,j)
            !ENDIF
          ENDIF
       ENDDO

       DO j=2, nvm
         IF ( natural(j) .AND. (.NOT. pasture(j)) .AND. (.NOT. is_peat(j)) ) THEN
           sumnat(i)=sumnat(i)+veget_cov_max_new(i,j)
           !IF (lalo(i,1)==27 .AND. lalo(i,2)==23) THEN
           !WRITE (numout,*) 'QCJ check natural upland,',j,veget_cov_max_old(i,j),veget_cov_max_new(i,j)
           !ENDIF   
         ENDIF

         IF ( ( .NOT. natural(j) .OR. pasture(j) ) .AND. (.NOT. is_croppeat(j)) ) THEN
           veget_cov_max_new_peatdgvm(i,j)=veget_cov_max_new(i,j)
           !IF (lalo(i,1)==27 .AND. lalo(i,2)==23) THEN
           !WRITE (numout,*) 'QCJ check crop,',j,veget_cov_max_old(i,j),veget_cov_max_new(i,j)
           !ENDIF    
           sumcrop(i)=sumcrop(i)+veget_cov_max_new(i,j)
         ENDIF 
         IF ( is_croppeat(j) ) THEN
           veget_cov_max_new_peatdgvm(i,j)=veget_cov_max_new(i,j) 
           sumcroppeat(i)=sumcroppeat(i)+veget_cov_max_new(i,j)
         ENDIF  
       ENDDO
       
       veget_cov_max_new_peatdgvm(i,ibare_sechiba)=veget_cov_max_old(i,ibare_sechiba)
       !IF (lalo(i,1)==27 .AND. lalo(i,2)==23) THEN
       !WRITE (numout,*) 'QCJ check baresoil,old,',veget_cov_max_old(i,ibare_sechiba)
       !ENDIF 
       sumnat_adjust(i)=un-sumnatpeat(i)-sumcrop(i)-sumcroppeat(i)-veget_cov_max_new_peatdgvm(i,ibare_sechiba)
       !IF (lalo(i,1)==27 .AND. lalo(i,2)==23) THEN
       ! WRITE (numout,*) 'QCJ check sumnat',sumnat(i)
       ! WRITE (numout,*) 'QCJ check sumnat_adjust',sumnat_adjust(i)
       !ENDIF    
       IF (sumnat_adjust(i) .LT. zero) THEN
          veget_cov_max_new_peatdgvm(i,ibare_sechiba)=veget_cov_max_new_peatdgvm(i,ibare_sechiba)+sumnat_adjust(i)
          sumnat_adjust(i)=zero  
        !  IF (lalo(i,1)==27 .AND. lalo(i,2)==23) THEN 
        !  WRITE (numout,*) 'QCJ check baresoil,new1',veget_cov_max_new_peatdgvm(i,ibare_sechiba)
        !  ENDIF
       ENDIF
        
       IF (sumnat(i) .NE. sumnat_adjust(i)) THEN
         IF (sumnat(i) .GT. zero) THEN
           DO j=2, nvm 
              IF ( natural(j) .AND. (.NOT. pasture(j)) .AND. (.NOT. is_peat(j)) ) THEN
                veget_cov_max_new_peatdgvm(i,j)=veget_cov_max_new(i,j)/sumnat(i)*sumnat_adjust(i)
              ENDIF 
           ENDDO
         ELSE 
           DO j=2, nvm
              IF ( natural(j) .AND. (.NOT. pasture(j)) .AND. (.NOT. is_peat(j)) ) THEN
                veget_cov_max_new_peatdgvm(i,j)=veget_cov_max_new(i,j)
              ENDIF
           ENDDO
           veget_cov_max_new_peatdgvm(i,ibare_sechiba)=veget_cov_max_new_peatdgvm(i,ibare_sechiba)+ &
                  (sumnat_adjust(i)-sumnat(i))   
         ENDIF
       ELSE   
         DO j=2, nvm
            IF ( natural(j) .AND. (.NOT. pasture(j)) .AND. (.NOT. is_peat(j)) ) THEN
               veget_cov_max_new_peatdgvm(i,j)=veget_cov_max_new(i,j)
            ENDIF
         ENDDO     
       ENDIF  

       !IF (lalo(i,1)==27 .AND. lalo(i,2)==23) THEN
       !   WRITE (numout,*) 'QCJ check,veget_cov_max_new_peatdgvm1,',veget_cov_max_new_peatdgvm(i,:),'sum',SUM(veget_cov_max_new_peatdgvm(i,:))
       !ENDIF

       IF ( ANY(veget_cov_max_new_peatdgvm(i,:) .LT. zero) ) THEN
          !WRITE(numout,*) 'lcchange_main veget_cov_max_new_peatdgvm<0,','i=',i,'veget_cov_max_new_peatdgvm=',veget_cov_max_new_peatdgvm(i,:)
          WRITE(numout,*) 'QCJ check veget_cov_max_new little than 0 in lcchange, at:',lalo(i,:)
          WRITE(numout,*) 'QCJ check veget_cov_max_new little than 0 ',veget_cov_max_new(i,:)
          !CALL ipslerr_p(3,'lcchange_main','ok_dgvm_peat= T','veget_cov_max_new_peatdgvm<0','')
       ENDIF

       delta_veg(:)=veget_cov_max_new_peatdgvm(i,:)-veget_cov_max_old(i,:)

       delta_veg_sum=zero
       DO j=1, nvm
         IF ( (.NOT. is_peat(j)) .AND. (.NOT. is_croppeat(j)) .AND. (delta_veg(j) .LT. -min_stomate) ) THEN
           delta_veg_sum = delta_veg_sum+delta_veg(j)
         ENDIF  
       ENDDO 

       dilu_lit(i,:,:,:) = zero
       dilu_soil_carbon(i,:) = zero
       biomass_loss(i,:,:) = zero

       !! 3.2 if vegetation coverage decreases, compute dilution of litter, soil carbon, and biomass.\n
       DO j=1, nvm
          IF ( (.NOT. is_peat(j)) .AND. (.NOT. is_croppeat(j)) .AND. (delta_veg(j) < -min_stomate) ) THEN
             dilu_lit(i,:,:,:) = dilu_lit(i,:,:,:) + delta_veg(j)*litter(i,:,j,:,:) / delta_veg_sum
             biomass_loss(i,:,:) = biomass_loss(i,:,:) + biomass(i,j,:,:)*delta_veg(j) / delta_veg_sum
             IF ( ok_pc ) THEN
                    dilu_soil_carbon_vertres(i,:,iactive)=dilu_soil_carbon_vertres(i,:,iactive) + &
                         delta_veg(j) * deepC_a(i,:,j) / delta_veg_sum
                    dilu_soil_carbon_vertres(i,:,islow)=dilu_soil_carbon_vertres(i,:,islow) + &
                         delta_veg(j) * deepC_s(i,:,j) / delta_veg_sum
                    dilu_soil_carbon_vertres(i,:,ipassive)=dilu_soil_carbon_vertres(i,:,ipassive) + &
                         delta_veg(j) * deepC_p(i,:,j) / delta_veg_sum
             ELSE
                    dilu_soil_carbon(i,:) =  dilu_soil_carbon(i,:) + delta_veg(j) * carbon(i,:,j) / delta_veg_sum
             ENDIF
          ENDIF
       ENDDO

       !! 3.3
       DO j=1, nvm ! Loop over # PFTs

         IF (.NOT. is_peat (j) .AND. (.NOT. is_croppeat(j)) ) THEN 
          !! 3.3.1 The case that vegetation coverage of PFTj increases
           IF ( delta_veg(j) > min_stomate ) THEN

             !! 3.3.1.1 Initial setting of new establishment
             IF (veget_cov_max_old(i,j) .LT. min_stomate) THEN
                IF (is_tree(j)) THEN

                   ! cn_sapl(j)=0.5; stomate_data.f90
                   cn_ind(i,j) = cn_sapl(j)
                ELSE
                   cn_ind(i,j) = un
                ENDIF
                ind(i,j)= delta_veg(j) / cn_ind(i,j)
                PFTpresent(i,j) = .TRUE.
                everywhere(i,j) = 1.
                senescence(i,j) = .FALSE.
                age(i,j) = zero

                ! large_value = 1.E33_r_std
                when_growthinit(i,j) = large_value
                leaf_frac(i,j,1) = 1.0
                npp_longterm(i,j) = npp_longterm_init
                lm_lastyearmax(i,j) = bm_sapl(j,ileaf,icarbon) * ind(i,j)
             ENDIF

             IF ( cn_ind(i,j) > min_stomate ) THEN
                delta_ind(i,j) = delta_veg(j) / cn_ind(i,j)
             ENDIF

             !! 3.3.1.2 Update of biomass in each each carbon stock component
             !!         Update of biomass in each each carbon stock component (leaf, sapabove, sapbelow,
             !>         heartabove, heartbelow, root, fruit, and carbres)\n

             DO k = 1, nparts ! loop over # carbon stock components, nparts = 8; stomate_constant.f90
                DO l = 1,nelements ! loop over # elements

                   bm_new(i,l) = delta_ind(i,j) * bm_sapl(j,k,l)
                   IF (veget_cov_max_old(i,j) .GT. min_stomate) THEN

                      ! in the case that bm_new is overestimated compared with biomass?
                      IF ((bm_new(i,l)/delta_veg(j)) > biomass(i,j,k,l)) THEN
                         bm_new(i,l) = biomass(i,j,k,l)*delta_veg(j)
                      ENDIF
                   ENDIF
                   biomass(i,j,k,l) = ( biomass(i,j,k,l) * veget_cov_max_old(i,j) + bm_new(i,l) ) / veget_cov_max_new_peatdgvm(i,j)
                   co2_to_bm(i,j) = co2_to_bm(i,j) + (bm_new(i,icarbon)* dt_days) / (one_year * veget_cov_max_new_peatdgvm(i,j))
                END DO ! loop over # elements
             ENDDO ! loop over # carbon stock components

             !! 3.3.1.3 Calculation of dilution in litter, soil carbon, and  input of litter
             !!        In this 'IF statement', dilu_* is zero. Formulas for litter and soil carbon
             !!         could be shortend?? Are the following formulas correct?
             ! Litter
             litter(i,:,j,:,:)=(litter(i,:,j,:,:) * veget_cov_max_old(i,j) + &
                  dilu_lit(i,:,:,:) * delta_veg(j)) / veget_cov_max_new_peatdgvm(i,j)
                !gmjc available and not available litter for grazing
                ! only not available litter increase/decrease, available litter will not
                ! change, due to tree litter can not be eaten
               IF (is_grassland_manag(j) .AND. is_grassland_grazed(j)) THEN
                 litter_avail(i,:,j) = litter_avail(i,:,j) * veget_cov_max_old(i,j) / veget_cov_max_new_peatdgvm(i,j)
                 litter_not_avail(i,:,j) = litter(i,:,j,iabove,icarbon) - litter_avail(i,:,j)
               ENDIF
                !end gmjc
             IF ( ok_pc ) THEN
                deepC_a(i,:,j)=(deepC_a(i,:,j) * veget_cov_max_old(i,j) + &
                     dilu_soil_carbon_vertres(i,:,iactive) * delta_veg(j)) / veget_cov_max_new_peatdgvm(i,j)
                deepC_s(i,:,j)=(deepC_s(i,:,j) * veget_cov_max_old(i,j) + &
                     dilu_soil_carbon_vertres(i,:,islow) * delta_veg(j)) / veget_cov_max_new_peatdgvm(i,j)
                deepC_p(i,:,j)=(deepC_p(i,:,j) * veget_cov_max_old(i,j) + &
                     dilu_soil_carbon_vertres(i,:,ipassive) * delta_veg(j)) / veget_cov_max_new_peatdgvm(i,j)
             ELSE
                ! Soil carbon
                carbon(i,:,j)=(carbon(i,:,j) * veget_cov_max_old(i,j) + dilu_soil_carbon(i,:) * delta_veg(j)) / veget_cov_max_new_peatdgvm(i,j)
             ENDIF

             DO l = 1,nelements

                ! Litter input
                bm_to_litter(i,j,isapbelow,l) = (bm_to_litter(i,j,isapbelow,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,isapbelow,l)*delta_veg(j)) / veget_cov_max_new_peatdgvm(i,j)
                bm_to_litter(i,j,iheartbelow,l) = (bm_to_litter(i,j,iheartbelow,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,iheartbelow,l) *delta_veg(j)) / veget_cov_max_new_peatdgvm(i,j)
                bm_to_litter(i,j,iroot,l) = (bm_to_litter(i,j,iroot,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,iroot,l)*delta_veg(j)) / veget_cov_max_new_peatdgvm(i,j)
                bm_to_litter(i,j,ifruit,l) = (bm_to_litter(i,j,ifruit,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,ifruit,l)*delta_veg(j)) / veget_cov_max_new_peatdgvm(i,j)
                bm_to_litter(i,j,icarbres,l) = (bm_to_litter(i,j,icarbres,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,icarbres,l)   *delta_veg(j)) / veget_cov_max_new_peatdgvm(i,j)
                bm_to_litter(i,j,ileaf,l) = (bm_to_litter(i,j,ileaf,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,ileaf,l)*delta_veg(j)) / veget_cov_max_new_peatdgvm(i,j)

             END DO

             age(i,j)=age(i,j)*veget_cov_max_old(i,j)/veget_cov_max_new_peatdgvm(i,j)

          !! 3.3.2 The case that vegetation coverage of PFTj is no change or decreases
           ELSE
             !! 3.3.2.1 Biomass export
             ! coeff_lcchange_*:  Coeff of biomass export for the year, decade, and century
             above = biomass(i,j,isapabove,icarbon) + biomass(i,j,iheartabove,icarbon)
             convflux(i,iwplcc)  = convflux(i,iwplcc)  - ( coeff_lcchange_1(j) * above * delta_veg(j) )
             prod10(i,0,iwplcc)  = prod10(i,0,iwplcc)  - ( coeff_lcchange_10(j) * above * delta_veg(j) )
             prod100(i,0,iwplcc) = prod100(i,0,iwplcc) - ( coeff_lcchange_100(j) * above * delta_veg(j) )

             !! If the vegetation is to small, it has been set to 0.
             IF ( veget_cov_max_new_peatdgvm(i,j) .LT. min_stomate ) THEN

                ind(i,j) = zero
                biomass(i,j,:,:) = zero
                PFTpresent(i,j) = .FALSE.
                senescence(i,j) = .FALSE.
                age(i,j) = zero
                when_growthinit(i,j) = undef
                everywhere(i,j) = zero
                carbon(i,:,j) = zero
                litter(i,:,j,:,:) = zero
                bm_to_litter(i,j,:,:) = zero

                IF (ok_pc) THEN
                   deepC_a(i,:,j) = zero
                   deepC_s(i,:,j) = zero
                   deepC_p(i,:,j) = zero
                ENDIF

             ENDIF

          ENDIF ! PFT's coverage reduction

         ENDIF ! .NOT. is_peat  .AND. (.NOT. is_croppeat(j))

       ENDDO  !Loop over # PFTs

!!! peatland (natural peat+agricultural peat)
    IF (agri_peat) THEN     

       sum_peat_old(i) = zero
       sum_peat_new(i) = zero
       sum_peatBare(i)=zero
       sum_peatBare_adjust(i)=zero
       DO j=1,nvm
          IF (is_peat(j) .OR. is_croppeat(j) .OR. j==ibare_sechiba) THEN
             sum_peat_old(i) = sum_peat_old(i) + veget_cov_max_old(i,j)
             sum_peat_new(i) = sum_peat_new(i) + veget_cov_max_new_peatdgvm(i,j)
          ENDIF
       ENDDO

       DO j=1,nvm
          IF (is_peat(j) .OR. j==ibare_sechiba) THEN
             sum_peatBare(i)=sum_peatBare(i)+veget_cov_max_new_peatdgvm(i,j)
          ENDIF
       ENDDO
       delta_peat(i) = sum_peat_new(i)-sum_peat_old(i)
       sum_peatBare_adjust(i)=sum_peatBare(i)-delta_peat(i)

       DO j=1, nvm
          IF ( ( is_peat(j) .OR. j==ibare_sechiba) .AND. (sum_peatBare(i) .GT. zero) ) THEN
             veget_cov_max_new_peatdgvm(i,j)=veget_cov_max_new_peatdgvm(i,j)*sum_peatBare_adjust(i)/sum_peatBare(i)
           
          ENDIF
       ENDDO

       delta_veg(:)=zero 
       DO j=1, nvm
         IF ( is_peat(j) .OR. is_croppeat(j) .OR. j==ibare_sechiba) THEN
           delta_veg(j) = veget_cov_max_new_peatdgvm(i,j)-veget_cov_max_old(i,j)
         ENDIF
       ENDDO

       delta_veg_sum=zero
       DO j=1, nvm
         IF ( ( is_peat(j) .OR. is_croppeat(j) .OR. j==ibare_sechiba) .AND. (delta_veg(j) .LT. -min_stomate) ) THEN
           delta_veg_sum = delta_veg_sum+delta_veg(j)
         ENDIF
       ENDDO

       dilu_lit(i,:,:,:) = zero
       dilu_soil_carbon(i,:) = zero
       biomass_loss(i,:,:) = zero

       !! if vegetation coverage decreases, compute dilution of litter, soil carbon, and biomass.\n
       DO j=1, nvm
          IF ( ( is_peat(j) .OR. is_croppeat(j) .OR. j==ibare_sechiba ) .AND. (delta_veg(j) < -min_stomate) ) THEN
             dilu_lit(i,:,:,:) = dilu_lit(i,:,:,:) + delta_veg(j)*litter(i,:,j,:,:) / delta_veg_sum
             biomass_loss(i,:,:) = biomass_loss(i,:,:) + biomass(i,j,:,:)*delta_veg(j) / delta_veg_sum
             IF ( ok_pc ) THEN
                    dilu_soil_carbon_vertres(i,:,iactive)=dilu_soil_carbon_vertres(i,:,iactive) + &
                         delta_veg(j) * deepC_a(i,:,j) / delta_veg_sum
                    dilu_soil_carbon_vertres(i,:,islow)=dilu_soil_carbon_vertres(i,:,islow) + &
                         delta_veg(j) * deepC_s(i,:,j) / delta_veg_sum
                    dilu_soil_carbon_vertres(i,:,ipassive)=dilu_soil_carbon_vertres(i,:,ipassive) + &
                         delta_veg(j) * deepC_p(i,:,j) / delta_veg_sum
             ELSE
                    dilu_soil_carbon(i,:) =  dilu_soil_carbon(i,:) + delta_veg(j) * carbon(i,:,j) / delta_veg_sum
             ENDIF
          ENDIF
       ENDDO
       DO j=1, nvm ! Loop over # PFTs

         IF ( is_peat(j) .OR. is_croppeat(j) .OR. j==ibare_sechiba ) THEN
          !! 3.3.1 The case that vegetation coverage of PFTj increases
          IF ( delta_veg(j) > min_stomate ) THEN

             !! 3.3.1.1 Initial setting of new establishment
             IF ( (veget_cov_max_old(i,j) .LT. min_stomate) .AND. (j/=ibare_sechiba) )  THEN
                IF (is_tree(j)) THEN
                   cn_ind(i,j) = cn_sapl(j)
                ELSE
                   cn_ind(i,j) = un
                ENDIF
                ind(i,j)= delta_veg(j) / cn_ind(i,j)
                PFTpresent(i,j) = .TRUE.
                everywhere(i,j) = 1.
                senescence(i,j) = .FALSE.
                age(i,j) = zero

               when_growthinit(i,j) = large_value
               leaf_frac(i,j,1) = 1.0
               npp_longterm(i,j) = npp_longterm_init
               lm_lastyearmax(i,j) = bm_sapl(j,ileaf,icarbon) * ind(i,j)
             ENDIF

             IF ( cn_ind(i,j) > min_stomate .AND. (j/=ibare_sechiba) ) THEN
                delta_ind(i,j) = delta_veg(j) / cn_ind(i,j)
             ENDIF


             !! 3.3.1.2 Update of biomass in each each carbon stock component
             !!         Update of biomass in each each carbon stock component (leaf, sapabove, sapbelow,
             !>         heartabove, heartbelow, root, fruit, and carbres)\n

             DO k = 1, nparts ! loop over # carbon stock components, nparts = 8; stomate_constant.f90
                DO l = 1,nelements ! loop over # elements

                   bm_new(i,l) = delta_ind(i,j) * bm_sapl(j,k,l)
                   IF (veget_cov_max_old(i,j) .GT. min_stomate) THEN

                      ! in the case that bm_new is overestimated compared with biomass?
                      IF ((bm_new(i,l)/delta_veg(j)) > biomass(i,j,k,l)) THEN
                         bm_new(i,l) = biomass(i,j,k,l)*delta_veg(j)
                      ENDIF
                   ENDIF
                   biomass(i,j,k,l) = ( biomass(i,j,k,l) * veget_cov_max_old(i,j) + bm_new(i,l) ) / veget_cov_max_new_peatdgvm(i,j)
                   co2_to_bm(i,j) = co2_to_bm(i,j) + (bm_new(i,icarbon)* dt_days) / (one_year * veget_cov_max_new_peatdgvm(i,j))
                END DO ! loop over # elements
             ENDDO ! loop over # carbon stock components

             !! 3.3.1.3 Calculation of dilution in litter, soil carbon, and  input of litter
             !!        In this 'IF statement', dilu_* is zero. Formulas for litter and soil carbon
             !!         could be shortend?? Are the following formulas correct?
             ! Litter
             litter(i,:,j,:,:)=(litter(i,:,j,:,:) * veget_cov_max_old(i,j) + &
                  dilu_lit(i,:,:,:) * delta_veg(j)) / veget_cov_max_new_peatdgvm(i,j)
                !gmjc available and not available litter for grazing
                ! only not available litter increase/decrease, available litter will not
                ! change, due to tree litter can not be eaten
               IF (is_grassland_manag(j) .AND. is_grassland_grazed(j)) THEN
                 litter_avail(i,:,j) = litter_avail(i,:,j) * veget_cov_max_old(i,j) / veget_cov_max_new_peatdgvm(i,j)
                 litter_not_avail(i,:,j) = litter(i,:,j,iabove,icarbon) - litter_avail(i,:,j)
               ENDIF
                !end gmjc

             IF ( ok_pc ) THEN
                deepC_a(i,:,j)=(deepC_a(i,:,j) * veget_cov_max_old(i,j) + &
                     dilu_soil_carbon_vertres(i,:,iactive) * delta_veg(j)) / veget_cov_max_new_peatdgvm(i,j)
                deepC_s(i,:,j)=(deepC_s(i,:,j) * veget_cov_max_old(i,j) + &
                     dilu_soil_carbon_vertres(i,:,islow) * delta_veg(j)) / veget_cov_max_new_peatdgvm(i,j)
                deepC_p(i,:,j)=(deepC_p(i,:,j) * veget_cov_max_old(i,j) + &
                     dilu_soil_carbon_vertres(i,:,ipassive) * delta_veg(j)) / veget_cov_max_new_peatdgvm(i,j)
             ELSE
                ! Soil carbon
                carbon(i,:,j)=(carbon(i,:,j) * veget_cov_max_old(i,j) + dilu_soil_carbon(i,:) * delta_veg(j)) / veget_cov_max_new_peatdgvm(i,j)
             ENDIF

             DO l = 1,nelements

                ! Litter input
                bm_to_litter(i,j,isapbelow,l) = (bm_to_litter(i,j,isapbelow,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,isapbelow,l)*delta_veg(j)) / veget_cov_max_new_peatdgvm(i,j)
                bm_to_litter(i,j,iheartbelow,l) = (bm_to_litter(i,j,iheartbelow,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,iheartbelow,l) *delta_veg(j)) / veget_cov_max_new_peatdgvm(i,j)
                bm_to_litter(i,j,iroot,l) = (bm_to_litter(i,j,iroot,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,iroot,l)*delta_veg(j)) / veget_cov_max_new_peatdgvm(i,j)
                bm_to_litter(i,j,ifruit,l) = (bm_to_litter(i,j,ifruit,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,ifruit,l)*delta_veg(j)) / veget_cov_max_new_peatdgvm(i,j)
                bm_to_litter(i,j,icarbres,l) = (bm_to_litter(i,j,icarbres,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,icarbres,l)   *delta_veg(j)) / veget_cov_max_new_peatdgvm(i,j)
                bm_to_litter(i,j,ileaf,l) = (bm_to_litter(i,j,ileaf,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,ileaf,l)*delta_veg(j)) / veget_cov_max_new_peatdgvm(i,j)

             END DO

             age(i,j)=age(i,j)*veget_cov_max_old(i,j)/veget_cov_max_new_peatdgvm(i,j)

          !! 3.3.2 The case that vegetation coverage of PFTj is no change or decreases
          ELSE
             !! 3.3.2.1 Biomass export
             ! coeff_lcchange_*:  Coeff of biomass export for the year, decade, and century
             above = biomass(i,j,isapabove,icarbon) + biomass(i,j,iheartabove,icarbon)
             convflux(i,iwplcc)  = convflux(i,iwplcc)  - ( coeff_lcchange_1(j) * above * delta_veg(j) )
             prod10(i,0,iwplcc)  = prod10(i,0,iwplcc)  - ( coeff_lcchange_10(j) * above * delta_veg(j) )
             prod100(i,0,iwplcc) = prod100(i,0,iwplcc) - ( coeff_lcchange_100(j) * above * delta_veg(j) )

             !! If the vegetation is to small, it has been set to 0.
             IF ( veget_cov_max_new_peatdgvm(i,j) .LT. min_stomate ) THEN

                ind(i,j) = zero
                biomass(i,j,:,:) = zero
                PFTpresent(i,j) = .FALSE.
                senescence(i,j) = .FALSE.
                age(i,j) = zero
                when_growthinit(i,j) = undef
                everywhere(i,j) = zero
                carbon(i,:,j) = zero
                litter(i,:,j,:,:) = zero
                bm_to_litter(i,j,:,:) = zero

                IF (ok_pc) THEN
                   deepC_a(i,:,j) = zero
                   deepC_s(i,:,j) = zero
                   deepC_p(i,:,j) = zero
                ENDIF

             ENDIF

          ENDIF ! PFT's coverage reduction

         ENDIF ! is_peat  .OR. is_croppeat(j)

       ENDDO  !Loop over # PFTs

    ENDIF
 
       !! 2.4 update 10 year-turnover pool content following flux emission
       !!     (linear decay (10%) of the initial carbon input)

       DO  l = 0, 8
          m = 10 - l
          cflux_prod10(i,iwplcc) =  cflux_prod10(i,iwplcc) + flux10(i,m,iwplcc)
          prod10(i,m,iwplcc)     =  prod10(i,m-1,iwplcc)   - flux10(i,m-1,iwplcc)
          flux10(i,m,iwplcc)     =  flux10(i,m-1,iwplcc)

          IF (prod10(i,m,iwplcc) .LT. 1.0) prod10(i,m,iwplcc) = zero
       ENDDO

       cflux_prod10(i,iwplcc) = cflux_prod10(i,iwplcc) + flux10(i,1,iwplcc)
       flux10(i,1,iwplcc)     = 0.1 * prod10(i,0,iwplcc)
       prod10(i,1,iwplcc)     = prod10(i,0,iwplcc)

       !! 2.5 update 100 year-turnover pool content following flux emission\n
       DO   l = 0, 98
          m = 100 - l
          cflux_prod100(i,iwplcc)  =  cflux_prod100(i,iwplcc) + flux100(i,m,iwplcc)
          prod100(i,m,iwplcc)      =  prod100(i,m-1,iwplcc)   - flux100(i,m-1,iwplcc)
          flux100(i,m,iwplcc)      =  flux100(i,m-1,iwplcc)

          IF (prod100(i,m,iwplcc).LT.1.0) prod100(i,m,iwplcc) = zero
       ENDDO

       cflux_prod100(i,iwplcc)  = cflux_prod100(i,iwplcc) + flux100(i,1,iwplcc)
       flux100(i,1,iwplcc)      = 0.01 * prod100(i,0,iwplcc)
       prod100(i,1,iwplcc)      = prod100(i,0,iwplcc)
       prod10(i,0,iwplcc)        = zero
       prod100(i,0,iwplcc)       = zero

    ENDDO   ! Loop over # pixels


    sumgridPFT=SUM(veget_cov_max_new_peatdgvm,DIM=2)
    DO i= 1, npts
       !WRITE (numout,*) 'QCJ check2,veget_cov_max_new_peatdgvm1,',veget_cov_max_new_peatdgvm(i,:),'sum',SUM(veget_cov_max_new_peatdgvm(i,:))
       IF ( (sumgridPFT(i) .LT. un-min_stomate) .OR. (sumgridPFT(i) .GT. un+min_stomate) ) THEN
          WRITE (numout,*) 'QCJ check sumgridPFT not equal to 1 in lcchange,at:',lalo(i,:)           
          WRITE (numout,*) 'QCJ check sumgridPFT=',sumgridPFT(i),'veget_cov_max=',veget_cov_max_new_peatdgvm(i,:)
          !CALL ipslerr_p(3,'stomate_lcchange','ok_dgvm_peat=T','sum of veget_cov_max_new_peatdgvm not equal to 1','') 
       ENDIF
    ENDDO    
   
    !! We redistribute the updated litter into four fuel classes, so that
    !! the balance between aboveground litter and fuel is mainted. The subtraction
    !! of fuel burned by land cover change fires from the fuel pool is made here.
    fuel_all_type(:,:,:,:) = fuel_1hr(:,:,:,:) + fuel_10hr(:,:,:,:) + &
                               fuel_100hr(:,:,:,:) + fuel_1000hr(:,:,:,:)
    fuel_type_frac(:,:,:,:,:) = 0.25
    WHERE(fuel_all_type(:,:,:,:) > min_stomate)
      fuel_type_frac(:,:,:,:,1) = fuel_1hr(:,:,:,:)/fuel_all_type(:,:,:,:)
      fuel_type_frac(:,:,:,:,2) = fuel_10hr(:,:,:,:)/fuel_all_type(:,:,:,:)
      fuel_type_frac(:,:,:,:,3) = fuel_100hr(:,:,:,:)/fuel_all_type(:,:,:,:)
      fuel_type_frac(:,:,:,:,4) = fuel_1000hr(:,:,:,:)/fuel_all_type(:,:,:,:)
    ENDWHERE

    DO j=1,nvm
      fuel_1hr(:,j,:,:) = litter(:,:,j,iabove,:) * fuel_type_frac(:,j,:,:,1)
      fuel_10hr(:,j,:,:) = litter(:,:,j,iabove,:) * fuel_type_frac(:,j,:,:,2)
      fuel_100hr(:,j,:,:) = litter(:,:,j,iabove,:) * fuel_type_frac(:,j,:,:,3)
      fuel_1000hr(:,j,:,:) = litter(:,:,j,iabove,:) * fuel_type_frac(:,j,:,:,4)
    END DO

  !! 3. history
    convflux        = convflux/one_year*dt_days
    cflux_prod10    = cflux_prod10/one_year*dt_days
    cflux_prod100   = cflux_prod100/one_year*dt_days

    C_end(:,1)= zero !SUM(SUM(biomass(:,:,:,icarbon),DIM=3) * veget_cov_max_new_peatdgvm,DIM=2)##Biomass,end should equal to Biomass,start+co2_to_bm
    C_end(:,2)=SUM(SUM(SUM(litter(:,:,:,:,icarbon),DIM=4),DIM=2)*veget_cov_max_new_peatdgvm,DIM=2)
    C_end(:,3)=zero !SUM(SUM(carbon(:,:,:),DIM=2) * veget_cov_max_new_peatdgvm,DIM=2) #carbon has not been adjusted when ok_pc=T

    DO l=1,ndeep
      deepC_end(:,l,1)=SUM(deepC_a(:,l,:)*veget_cov_max_new_peatdgvm(:,:),DIM=2)
      deepC_end(:,l,2)=SUM(deepC_s(:,l,:)*veget_cov_max_new_peatdgvm(:,:),DIM=2)
      deepC_end(:,l,3)=SUM(deepC_p(:,l,:)*veget_cov_max_new_peatdgvm(:,:),DIM=2)
    ENDDO

    Fin_end(:,1)=zero !SUM(SUM(bm_to_litter(:,:,:,icarbon),DIM=3) * veget_cov_max_new_peatdgvm,DIM=2) ##bm_to_litter,end should equal to bm_ti_litter,start+biomass_loss
    Fin_end(:,2)=zero
    Fin_end(:,3)=zero !SUM(co2_to_bm * veget_cov_max_new_peatdgvm,DIM=2) ##co2_to_bm,end should equal to co2_to_bm,start+bm_new
    Fin_end(:,4)=zero

    DO i = 1, npts ! loop over grid points
        C_check(i,:)=C_end(i,:)-C_start(i,:)
        !WRITE (numout,*) 'QCJ check lcchange_main, veget_cov_max_old,',veget_cov_max_old(i,:)
        !WRITE (numout,*) 'QCJ check lcchange_main, veget_cov_max_new_peatdgvm,',veget_cov_max_new_peatdgvm(i,:) 
        !WRITE (numout,*) 'QCJ check lcchange_main, delta vege_max,',veget_cov_max_new_peatdgvm(i,:)-veget_cov_max_old(i,:)

        !IF ( ANY(C_check(i,:) .GT. min_stomate*100.) ) THEN
           ! WRITE (numout,*) 'QCJ NOTE, C pool mass not conserved after adjust VEGET_COV_MAX,at grid',i,'Coordinates:',lalo(i,:)
           ! WRITE (numout,*) 'QCJ lcchange_main, biomass, litter, SOC',C_check(i,:)
            !WRITE (numout,*) 'QCJ delta vege_max,',veget_cov_max_new_peatdgvm(i,:)-veget_cov_max_old(i,:)
        !ENDIF

        deepC_check(i,:,:)=deepC_end(i,:,:)-deepC_start(i,:,:)
        !IF ( ANY(deepC_check(i,:,:) .GT. min_stomate*100.) ) THEN
          ! WRITE (numout,*) 'QCJ NOTE, deepC pool mass not conserved after adjust VEGET_COV_MAX,at grid',i,'Coordinates:',lalo(i,:)
          ! WRITE (numout,*) 'QCJ lcchange_main, active first 10 layers',deepC_check(i,:10,1)
          ! WRITE (numout,*) 'QCJ lcchange_main, slow first 10 layers',deepC_check(i,:10,2)
          ! WRITE (numout,*) 'QCJ lcchange_main, passive first 10 layers',deepC_check(i,:10,3) 
        !ENDIF

        Fin_check(i,:)=Fin_end(i,:)-Fin_start(i,:)
        !IF ( ANY(Fin_check(i,:) .GT. min_stomate*100.) ) THEN
         !  WRITE (numout,*) 'QCJ NOTE, Fin mass not conserved after adjust VEGET_COV_MAX, at grid',i
         !  WRITE (numout,*) 'QCJ lcchange_main, bm_to_litter, turnover_daily, co2_to_bm, gpp_daily',Fin_check(i,:)
        !ENDIF
    ENDDO

    IF (printlev>=4) WRITE(numout,*) 'Leaving lcchange_main, ok_dgvm_peat'

  ENDIF
   
  END SUBROUTINE lcchange_main
  

  !! The lcchange modelling including consideration of deforestation fires
  SUBROUTINE lcchange_deffire ( npts, dt_days, veget_cov_max, veget_cov_max_new,&
       biomass, ind, age, PFTpresent, senescence, when_growthinit, everywhere, &        
       co2_to_bm, bm_to_litter, turnover_daily, bm_sapl, cn_ind,flux10,flux100, &
       prod10,prod100,&
       convflux,&
       cflux_prod10,cflux_prod100, leaf_frac,&
       npp_longterm, lm_lastyearmax, litter, litter_avail, litter_not_avail, &
       carbon,&
       deepC_a, deepC_s, deepC_p,&
       fuel_1hr,fuel_10hr,fuel_100hr,fuel_1000hr,&
       lcc,bafrac_deforest_accu,emideforest_litter_accu,emideforest_biomass_accu,&
       deflitsup_total,defbiosup_total)
    
    IMPLICIT NONE
    
    !! 0. Variable and parameter declaration 
    
    !! 0.1 Input variables
    
    INTEGER, INTENT(in)                                       :: npts             !! Domain size - number of pixels (unitless)
    REAL(r_std), INTENT(in)                                   :: dt_days          !! Time step of vegetation dynamics for stomate
                                                                                  !! (days)
    REAL(r_std), DIMENSION(nvm, nparts,nelements), INTENT(in) :: bm_sapl          !! biomass of sapling 
                                                                                  !! @tex ($gC individual^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)                    :: bafrac_deforest_accu !!cumulative deforestation fire burned fraction, unitless
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements), INTENT(in)    :: emideforest_litter_accu !!cumulative deforestation fire carbon emissions from litter
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(in)   :: emideforest_biomass_accu !!cumulative deforestation fire carbon emissions from tree biomass
    REAL(r_std), DIMENSION(npts,nvm),INTENT(in)                     :: lcc !! land cover change happened at this day

    !! 0.2 Output variables

    REAL(r_std), DIMENSION(:,:), INTENT(out)                 :: convflux       !! release during first year following land cover
                                                                               !! change (npts,nwp)
    REAL(r_std), DIMENSION(:,:), INTENT(out)                 :: cflux_prod10   !! total annual release from the 10 year-turnover
                                                                               !! pool @tex ($gC m^{-2}$) @endtex
                                                                               !! (npts,nwp)
    REAL(r_std), DIMENSION(:,:), INTENT(out)                 :: cflux_prod100  !! total annual release from the 100 year-
                                                                               !! turnover pool @tex ($gC m^{-2}$) @endtex
                                                                               !! (npts,nwp)
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout):: turnover_daily   !! Turnover rates 
                                                                                      !! @tex ($gC m^{-2} day^{-1}$) @endtex

    !! 0.3 Modified variables   
    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: veget_cov_max        !! "maximal" coverage fraction of a PFT (LAI ->
                                                                                  !! infinity) on ground (unitless) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: veget_cov_max_new    !! new "maximal" coverage fraction of a PFT (LAI
                                                                                  !! -> infinity) on ground
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout):: biomass    !! biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: ind              !! Number of individuals @tex ($m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: age              !! mean age (years)
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)               :: senescence       !! plant senescent (only for deciduous trees) Set
                                                                                  !! to .FALSE. if PFT is introduced or killed
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)               :: PFTpresent       !! Is pft there (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: everywhere       !! is the PFT everywhere in the grid box or very 
                                                                                  !! localized (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: when_growthinit  !! how many days ago was the beginning of the 
                                                                                  !! growing season (days)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: co2_to_bm        !! biomass uptaken 
                                                                                  !! @tex ($gC m^{-2} day^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: bm_to_litter !! conversion of biomass to litter 
                                                                                  !! @tex ($gC m^{-2} day^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: cn_ind           !! crown area of individuals 
                                                                                  !! @tex ($m^{2}$) @endtex
    REAL(r_std), DIMENSION(npts,0:10,nwp), INTENT(inout)          :: prod10               !! products remaining in the 10 year-turnover
                                                                                  !! pool after the annual release for each 
                                                                                  !! compartment (10 + 1 : input from year of land
                                                                                  !! cover change) (npts,0:10,nwp)
    REAL(r_std), DIMENSION(npts,0:100,nwp), INTENT(inout)         :: prod100               !! products remaining in the 100 year-turnover
                                                                                  !! pool after the annual release for each 
                                                                                  !! compartment (100 + 1 : input from year of land
                                                                                  !! cover change)  (npts,0:100,nwp)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)            :: flux10             !! annual release from the 10/100 year-turnover 
                                                                                  !! pool compartments (npts,0:10,nwp)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)           :: flux100             !! annual release from the 10/100 year-turnover
                                                                                  !! pool compartments (npts,0:100,nwp)
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout) :: leaf_frac        !! fraction of leaves in leaf age class 
                                                                                  !! (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: lm_lastyearmax   !! last year's maximum leaf mass for each PFT 
                                                                                  !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: npp_longterm     !! "long term" net primary productivity 
                                                                                  !! @tex ($gC m^{-2} year^{-1}$) @endtex
    REAL(r_std),DIMENSION(npts,nlitt,nvm,nlevs,nelements), INTENT(inout):: litter !! metabolic and structural litter, above and 
                                                                                  !! below ground @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nlitt,nvm), INTENT(inout):: litter_avail
    REAL(r_std), DIMENSION(npts,nlitt,nvm) , INTENT(inout):: litter_not_avail
    REAL(r_std),DIMENSION(npts,ncarb,nvm), INTENT(inout)      :: carbon           !! carbon pool: active, slow, or passive 
                                                                                  !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,ndeep,nvm), INTENT(inout)         :: deepC_a      !! Permafrost soil carbon (g/m**3) active
    REAL(r_std), DIMENSION(npts,ndeep,nvm), INTENT(inout)         :: deepC_s      !! Permafrost soil carbon (g/m**3) slow
    REAL(r_std), DIMENSION(npts,ndeep,nvm), INTENT(inout)         :: deepC_p      !! Permafrost soil carbon (g/m**3) passive

    REAL(r_std),DIMENSION(npts,nvm), INTENT(inout)                :: deflitsup_total
    REAL(r_std),DIMENSION(npts,nvm), INTENT(inout)                :: defbiosup_total
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements),INTENT(inout)        :: fuel_1hr
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements),INTENT(inout)        :: fuel_10hr
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements),INTENT(inout)        :: fuel_100hr
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements),INTENT(inout)        :: fuel_1000hr

    !! 0.4 Local variables

    INTEGER(i_std)                                            :: i, j, k, l, m, ilit, ipart    !! indices (unitless)
    REAL(r_std),DIMENSION(npts,nelements)                     :: bm_new           !! biomass increase @tex ($gC m^{-2}$) @endtex
    REAL(r_std),DIMENSION(npts,nparts,nelements)              :: biomass_loss     !! biomass loss @tex ($gC m^{-2}$) @endtex
    REAL(r_std)                                               :: above            !! aboveground biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std),DIMENSION(npts,nlitt,nlevs,nelements)         :: dilu_lit         !! Litter dilution @tex ($gC m^{-2}$) @endtex
    REAL(r_std),DIMENSION(npts,ncarb)                         :: dilu_soil_carbon !! Soil Carbondilution @tex ($gC m^{-2}$) @endtex
    REAL(r_std),DIMENSION(npts,ndeep,ncarb)                   :: dilu_soil_carbon_vertres !!vertically-resolved Soil Carbondilution (gC/m²)

    REAL(r_std),DIMENSION(nvm)                                :: delta_veg        !! changes in "maximal" coverage fraction of PFT 
    REAL(r_std)                                               :: delta_veg_sum    !! sum of delta_veg
    REAL(r_std),DIMENSION(npts,nvm)                           :: delta_ind        !! change in number of individuals  
    REAL(r_std),DIMENSION(npts,nvm,nlitt)                     :: deforest_litter_surplus !! Surplus in ground litter for deforested land after 
                                                                                         !! accounting for fire emissions
    REAL(r_std),DIMENSION(npts,nvm,nparts)                    :: deforest_biomass_surplus !!Surplus in live biomass for deforested forest
                                                                                          !!after accounting for fire emissions 
    REAL(r_std),DIMENSION(npts,nvm,nlitt)                     :: deforest_litter_deficit
    REAL(r_std),DIMENSION(npts,nvm,nparts)                    :: deforest_biomass_deficit
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements)       :: fuel_all_type
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements,4)     :: fuel_type_frac

    REAL(r_std),DIMENSION(npts,nvm)                           :: pool_start_pft        !! change in number of individuals  
    REAL(r_std),DIMENSION(npts)                           :: pool_start        !! change in number of individuals  
    REAL(r_std),DIMENSION(npts,nvm)                           :: pool_end_pft        !! change in number of individuals  
    REAL(r_std),DIMENSION(npts)                           :: pool_end        !! change in number of individuals  
    REAL(r_std),DIMENSION(npts)                           :: outflux        !! change in number of individuals  
!_ ================================================================================================================================



    pool_start_pft(:,:) = SUM(biomass(:,:,:,icarbon),DIM=3) &
                    + SUM(SUM(litter(:,:,:,:,icarbon),DIM=2),DIM=3) &
                    + SUM(carbon(:,:,:),DIM=2) &
                    + SUM(bm_to_litter(:,:,:,icarbon),DIM=3) &
                    + SUM(turnover_daily(:,:,:,icarbon),DIM=3)
                    
    pool_start(:) = SUM(pool_start_pft(:,:)*veget_cov_max(:,:),DIM=2) &
                    + SUM(prod10(:,:,iwplcc),DIM=2) + SUM(prod100(:,:,iwplcc),DIM=2)
                    

    deforest_biomass_surplus(:,:,:) = zero
    deforest_litter_surplus(:,:,:) = zero
    deforest_biomass_deficit(:,:,:) = zero
    deforest_litter_deficit(:,:,:) = zero

    IF (printlev>=3) WRITE(numout,*) 'Entering lcchange_main'
    
  !! 1. initialization
    
    prod10(:,0,:)       = zero
    prod100(:,0,:)      = zero   
    above               = zero
    convflux(:,:)       = zero
    cflux_prod10(:,:)   = zero
    cflux_prod100(:,:)  = zero
    delta_ind(:,:)      = zero
    delta_veg(:)        = zero
    dilu_soil_carbon_vertres(:,:,:) =zero
  !! 2. calculation of changes in carbon stocks and biomass by land cover change\n
    
    DO i = 1, npts ! Loop over # pixels - domain size
       
       !! 2.1 initialization of carbon stocks\n
       delta_veg(:) = veget_cov_max_new(i,:)-veget_cov_max(i,:)
       delta_veg_sum = SUM(delta_veg,MASK=delta_veg.LT.0.) !note `delta_veg_sum` is a negative number
       
       dilu_lit(i,:,:,:) = zero
       dilu_soil_carbon(i,:) = zero
       biomass_loss(i,:,:) = zero
       
       !! 2.2 Compute dilution pool of litter, soil carbon, and biomass for 
       !! decreasing PFTs.
       DO j=2, nvm
          IF ( delta_veg(j) < -min_stomate ) THEN 

             ! We make distinction between tree and grass because tree cover reduction might be due to fires. 
             ! The litter that is burned in fire should be excluded from diluting litter pool.
             IF (is_tree(j)) THEN
                deforest_litter_surplus(i,j,:) = -1*delta_veg(j)*litter(i,:,j,iabove,icarbon) - emideforest_litter_accu(i,j,:,icarbon)

                ! Here we compensate the litter burned by deforestation fire if it's higher than the litter available for
                ! burning. It follows the same logic as biomass which is described below.
                DO ilit = 1,nlitt
                  IF (deforest_litter_surplus(i,j,ilit) < zero) THEN
                     IF (veget_cov_max_new(i,j) < min_stomate) THEN
                        !WRITE (numout,*) 'Cumulative deforestation fire emission exceeds litter for point',i,',PFT ',j, &
                        !                 'However the new veget_cov_max is zero, there is not remaining litter to be diluted'
                        !STOP
                        deforest_litter_deficit(i,j,ilit) = deforest_litter_surplus(i,j,ilit) 

                     ELSE IF (litter(i,ilit,j,iabove,icarbon)*veget_cov_max_new(i,j) < -deforest_litter_surplus(i,j,ilit)) THEN
                        !WRITE (numout,*) 'Cumulative deforestation fire emission exceeds litter for point',i,',PFT ',j, &
                        !                 'However the remaing litter is not engough for diluting'
                        !STOP 
                        deforest_litter_deficit(i,j,ilit) = deforest_litter_surplus(i,j,ilit) 
                     ELSE
                        litter(i,ilit,j,iabove,icarbon) = ( litter(i,ilit,j,iabove,icarbon)*veget_cov_max_new(i,j) &
                               + deforest_litter_surplus(i,j,ilit) )/veget_cov_max_new(i,j)
                     END IF
                  ELSE
                     dilu_lit(i,ilit,iabove,icarbon) = dilu_lit(i,ilit,iabove,icarbon) -1 * deforest_litter_surplus(i,j,ilit)
                  END IF
                END DO
                dilu_lit(i,:,ibelow,:) = dilu_lit(i,:,ibelow,:) + delta_veg(j)*litter(i,:,j,ibelow,:) 
             ELSE
                dilu_lit(i,:,:,:) = dilu_lit(i,:,:,:) + delta_veg(j)*litter(i,:,j,:,:) 
             END IF

             IF (is_tree(j)) THEN
                deforest_biomass_surplus(i,j,:) = -1*delta_veg(j)*biomass(i,j,:,icarbon) - emideforest_biomass_accu(i,j,:,icarbon)
                ! Here we check if the biomass burned by deforestation fires is higher than the amount
                ! that could be deforested, if yes, the extra burned biomass is compensated by the biomass
                ! that is not deforested. Here we assume that if this happens for one deforested tree PFT,
                ! it happens for all deforested tree PFTs, so that we don't assume this extra burned biomass
                ! could be compenstated by other tree PFTs.
                DO ipart = 1,nparts
                  IF (deforest_biomass_surplus(i,j,ipart) < zero) THEN
                     IF (veget_cov_max_new(i,j) < min_stomate) THEN
                        !WRITE (numout,*) 'Cumulative deforestation fire emission exceeds biomass for point',i,',PFT ',j, &
                        !                 'However the new veget_cov_max is zero, there is not remaining biomass to be diluted'
                        !STOP 
                        deforest_biomass_deficit(i,j,ipart) = deforest_biomass_surplus(i,j,ipart)

                     ELSE IF (biomass(i,j,ipart,icarbon)*veget_cov_max_new(i,j) < -deforest_biomass_surplus(i,j,ipart)) THEN
                        !WRITE (numout,*) 'Cumulative deforestation fire emission exceeds biomass for point',i,',PFT ',j, &
                        !                 'However the remaing biomass is not engough for diluting'
                        !STOP 
                        deforest_biomass_deficit(i,j,ipart) = deforest_biomass_surplus(i,j,ipart)
                     ELSE
                        biomass(i,j,ipart,icarbon) = ( biomass(i,j,ipart,icarbon)*veget_cov_max_new(i,j) &
                               + deforest_biomass_surplus(i,j,ipart) )/veget_cov_max_new(i,j)
                     END IF
                  ELSE
                     biomass_loss(i,ipart,icarbon) = biomass_loss(i,ipart,icarbon) -1 * deforest_biomass_surplus(i,j,ipart)
                  END IF
                END DO
             ELSE 
                biomass_loss(i,:,:) = biomass_loss(i,:,:) + biomass(i,j,:,:)*delta_veg(j) 
             END IF 

             !IF (ANY( deforest_biomass_surplus(i,j,:) .LT. 0.0 ) .OR. ANY( deforest_litter_surplus(i,j,:) .LT. 0.0 ) ) THEN
             !   STOP 'Negative biomass or litter surplus'
             !ENDIF 

             IF ( ok_pc ) THEN
                    dilu_soil_carbon_vertres(i,:,iactive)=dilu_soil_carbon_vertres(i,:,iactive) + &
                         delta_veg(j) * deepC_a(i,:,j) / delta_veg_sum
                    dilu_soil_carbon_vertres(i,:,islow)=dilu_soil_carbon_vertres(i,:,islow) + &
                         delta_veg(j) * deepC_s(i,:,j) / delta_veg_sum
                    dilu_soil_carbon_vertres(i,:,ipassive)=dilu_soil_carbon_vertres(i,:,ipassive) + &
                         delta_veg(j) * deepC_p(i,:,j) / delta_veg_sum
             ELSE
                    dilu_soil_carbon(i,:) =  dilu_soil_carbon(i,:) + delta_veg(j) * carbon(i,:,j) / delta_veg_sum
             ENDIF
          ENDIF
       ENDDO !nbpts


       ! Note here `biomass_loss` and `dilu_lit` will change their sign from negative to positive
       IF ( delta_veg_sum < -min_stomate ) THEN
         biomass_loss(i,:,:) = biomass_loss(i,:,:) / delta_veg_sum
         dilu_lit(i,:,:,:) = dilu_lit(i,:,:,:) / delta_veg_sum
       END IF

       
       !! 2.3 Dilut the litter, soil carbon from decreasing PFTs to increasing ones.
       !! Establish new biomass for increasing PFTs.
       DO j=2, nvm ! Loop over # PFTs

          !! 2.3.1 The case that vegetation coverage of PFTj increases
          IF ( delta_veg(j) > min_stomate) THEN

             !! 2.3.1.1 The PFTj increased from zero to non-zeor, we have to 
             !! initialize it by setting new establishment
             IF (veget_cov_max(i,j) .LT. min_stomate) THEN 
                IF (is_tree(j)) THEN
                   cn_ind(i,j) = cn_sapl(j) ! cn_sapl(j)=0.5; stomate_data.f90
                ELSE
                   cn_ind(i,j) = un
                ENDIF

                ind(i,j)= delta_veg(j) / cn_ind(i,j)
                PFTpresent(i,j) = .TRUE.
                everywhere(i,j) = 1.
                senescence(i,j) = .FALSE.
                age(i,j) = zero
                when_growthinit(i,j) = large_value ! large_value = 1.E33_r_std
                leaf_frac(i,j,1) = 1.0
                npp_longterm(i,j) = npp_longterm_init
                lm_lastyearmax(i,j) = bm_sapl(j,ileaf,icarbon) * ind(i,j)
             ENDIF

             
             ! Calculate individual density increase because of coverage increase
             IF ( cn_ind(i,j) > min_stomate ) THEN
                delta_ind(i,j) = delta_veg(j) / cn_ind(i,j) 
             ENDIF
             !! 2.3.1.2 The increase in `ind` should be companied by increase in 
             !! biomass, we do this by assuming increased `ind` are saplings.
             DO k = 1, nparts ! loop over # carbon stock components, nparts = 8; stomate_constant.f90 
                DO l = 1,nelements ! loop over # elements
                   bm_new(i,l) = delta_ind(i,j) * bm_sapl(j,k,l) 
                   IF (veget_cov_max(i,j) .GT. min_stomate) THEN
                      ! Adjust bm_new equal to existing biomass if it's 
                      ! larger than the latter 
                      IF ((bm_new(i,l)/delta_veg(j)) > biomass(i,j,k,l)) THEN
                         bm_new(i,l) = biomass(i,j,k,l)*delta_veg(j)
                      ENDIF
                   ENDIF
                   biomass(i,j,k,l) = ( biomass(i,j,k,l) * veget_cov_max(i,j) + bm_new(i,l) ) / veget_cov_max_new(i,j)
                   co2_to_bm(i,j) = co2_to_bm(i,j) + (bm_new(i,icarbon)* dt_days) / (one_year * veget_cov_max_new(i,j))
                END DO ! loop over # elements
             ENDDO ! loop over # carbon stock components

             !! 2.3.1.3 Tow tasks are done here:
             !! A. We transfer the litter and soil carbon from the 
             !! reduced PFTs to the increases PFTs.

             ! Litter
             litter(i,:,j,:,:)=(litter(i,:,j,:,:) * veget_cov_max(i,j) + &
                  dilu_lit(i,:,:,:) * delta_veg(j)) / veget_cov_max_new(i,j)

               !!######################This part needs to be discussed with JinFeng ############
                !gmjc available and not available litter for grazing
                ! only not available litter increase/decrease, available litter will not
                ! change, due to tree litter can not be eaten
               IF (is_grassland_manag(j) .AND. is_grassland_grazed(j)) THEN
                 litter_avail(i,:,j) = litter_avail(i,:,j) * veget_cov_max(i,j) / veget_cov_max_new(i,j)
                 litter_not_avail(i,:,j) = litter(i,:,j,iabove,icarbon) - litter_avail(i,:,j)
               ENDIF
                !end gmjc            
               !!###############################################################################

             ! Soil carbon
             IF ( ok_pc ) THEN 
                deepC_a(i,:,j)=(deepC_a(i,:,j) * veget_cov_max(i,j) + &
                     dilu_soil_carbon_vertres(i,:,iactive) * delta_veg(j)) / veget_cov_max_new(i,j)
                deepC_s(i,:,j)=(deepC_s(i,:,j) * veget_cov_max(i,j) + &
                     dilu_soil_carbon_vertres(i,:,islow) * delta_veg(j)) / veget_cov_max_new(i,j)
                deepC_p(i,:,j)=(deepC_p(i,:,j) * veget_cov_max(i,j) + &
                     dilu_soil_carbon_vertres(i,:,ipassive) * delta_veg(j)) / veget_cov_max_new(i,j)
             ELSE
                carbon(i,:,j)=(carbon(i,:,j) * veget_cov_max(i,j) + dilu_soil_carbon(i,:) * delta_veg(j)) / veget_cov_max_new(i,j)
             ENDIF

             !! B. For the biomass pool of reducing PFTs, we cannot transfer them directly to the 
             !! increasing PFTs, because the latter ones are treated with new sapling estalishement
             !! in section 2.3.1.2. So we assume the non-harvestable biomass of reducing PFTs will
             !! go to litter pool via `bm_to_litter`, and these are further directly transferred to
             !! the increasing PFTs.
             !! 
             !! The non-harvestable parts are: isapbelow,iheartbelow,iroot,icarbres,ileaf,ifruit
             !! Note that the icarbres,ileaf,ifruit could be burned in deforestation fires, the 
             !! emissions from these parts are already subtracted from `biomass_loss`, as done 
             !! in section 2.2. The harvestable biomass parts go to harvest pool and this will done
             !! in the section for the reducing PFTs.
             DO l = 1,nelements

                bm_to_litter(i,j,isapbelow,l) = bm_to_litter(i,j,isapbelow,l) + &
                                                & biomass_loss(i,isapbelow,l)*delta_veg(j) / veget_cov_max_new(i,j)
                bm_to_litter(i,j,iheartbelow,l) = bm_to_litter(i,j,iheartbelow,l) + biomass_loss(i,iheartbelow,l) *delta_veg(j) &
                                                  &  / veget_cov_max_new(i,j)
                bm_to_litter(i,j,iroot,l) = bm_to_litter(i,j,iroot,l) + biomass_loss(i,iroot,l)*delta_veg(j) / veget_cov_max_new(i,j)
                bm_to_litter(i,j,ifruit,l) = bm_to_litter(i,j,ifruit,l) + biomass_loss(i,ifruit,l)*delta_veg(j) / veget_cov_max_new(i,j)
                bm_to_litter(i,j,icarbres,l) = bm_to_litter(i,j,icarbres,l) + &
                                               & biomass_loss(i,icarbres,l)   *delta_veg(j) / veget_cov_max_new(i,j)
                bm_to_litter(i,j,ileaf,l) = bm_to_litter(i,j,ileaf,l) + biomass_loss(i,ileaf,l)*delta_veg(j) / veget_cov_max_new(i,j)
             END DO

             age(i,j)=age(i,j)*veget_cov_max(i,j)/veget_cov_max_new(i,j)
             
          !! 2.3.2 The case that vegetation coverage of PFTj has no change or decreases.
          ELSE 
             
             !! 2.3.2.1 Complete disappearing of PFTj, i.e., changes from non-zero
             !! to zero.
             IF ( veget_cov_max_new(i,j) .LT. min_stomate ) THEN 
                veget_cov_max_new(i,j)= zero
                ind(i,j) = zero
                biomass(i,j,:,:) = zero
                PFTpresent(i,j) = .FALSE.
                senescence(i,j) = .FALSE.
                age(i,j) = zero
                when_growthinit(i,j) = undef
                everywhere(i,j) = zero
                carbon(i,:,j) = zero
                litter(i,:,j,:,:) = zero
                litter_avail(i,:,j) = zero
                litter_not_avail(i,:,j) = zero
                bm_to_litter(i,j,:,:) = zero
                turnover_daily(i,j,:,:) = zero
                deepC_a(i,:,j) = zero
                deepC_s(i,:,j) = zero
                deepC_p(i,:,j) = zero
             ENDIF
          ENDIF ! The end the two cases: PFT-coverage reduction versus 
                ! non-change-or-increase
       ENDDO ! 2.3 Loop over # PFTs

       !! 2.4 Biomass harvest and turnover of different harvest pools

       !!?? Here we have some problem regarding grassland/cropland area dereasing,
       !!?? Because their sapwood/heartwood aboveground are also treated as 
       !!?? wood products.

       !! 2.4.1 We have already deforestation fire fluxes from sapwood/hearwood aboveground,
       !! now we just assume the remaining unburned parts are harvested, as 10-year and
       !! 100-year product pool.

    print *,'delta_veg_sum',delta_veg_sum
    print *,'prod10_in_lcc_before_assign',prod10(:,:,iwplcc)
    print *,'biomass_loss',biomass_loss(:,:,:)
       ! Note before we divide biomass_loss by `delta_veg_sum` to convert it based on PFT area,
       ! Now we multiply it again by `delta_veg_sum` to convert it back based on grid cell area.
       ! Also note `delta_veg_sum` is negative, so we should multiply again by (-1)
       above = (biomass_loss(i,isapabove,icarbon) + biomass_loss(i,iheartabove,icarbon))*delta_veg_sum*(-1)
       convflux(i,iwplcc)  = SUM(emideforest_biomass_accu(i,:,isapabove,icarbon)+emideforest_biomass_accu(i,:,iheartabove,icarbon))
       prod10(i,0,iwplcc)  = 0.4* above
       prod100(i,0,iwplcc) = 0.6 * above 
       print *,'above_in_lcc_before_assign',above

       !! 2.4.2 update 10 year-turnover pool content following flux emission
       !!     (linear decay (10%) of the initial carbon input)
       DO  l = 0, 8
          m = 10 - l
          cflux_prod10(i,iwplcc) =  cflux_prod10(i,iwplcc) + flux10(i,m,iwplcc)
          prod10(i,m,iwplcc)     =  prod10(i,m-1,iwplcc)   - flux10(i,m-1,iwplcc)
          flux10(i,m,iwplcc)     =  flux10(i,m-1,iwplcc)
          IF (prod10(i,m,iwplcc) .LT. 1.0) prod10(i,m,iwplcc) = zero

       ENDDO
       
       cflux_prod10(i,iwplcc) = cflux_prod10(i,iwplcc) + flux10(i,1,iwplcc) 
       flux10(i,1,iwplcc)     = 0.1 * prod10(i,0,iwplcc)
       prod10(i,1,iwplcc)     = prod10(i,0,iwplcc)
       
       !! 3.5 update 100 year-turnover pool content following flux emission\n
       DO   l = 0, 98
          m = 100 - l
          cflux_prod100(i,iwplcc)  =  cflux_prod100(i,iwplcc) + flux100(i,m,iwplcc)
          prod100(i,m,iwplcc)      =  prod100(i,m-1,iwplcc)   - flux100(i,m-1,iwplcc)
          flux100(i,m,iwplcc)      =  flux100(i,m-1,iwplcc)
          
          IF (prod100(i,m,iwplcc).LT.1.0) prod100(i,m,iwplcc) = zero
       ENDDO
       
       cflux_prod100(i,iwplcc)  = cflux_prod100(i,iwplcc) + flux100(i,1,iwplcc) 
       flux100(i,1,iwplcc)      = 0.01 * prod100(i,0,iwplcc)
       prod100(i,1,iwplcc)      = prod100(i,0,iwplcc)
       prod10(i,0,iwplcc)        = zero
       prod100(i,0,iwplcc)       = zero 

    ENDDO ! Loop over # pixels - domain size
    print *,'prod10_in_lcc_after_assign',prod10(:,:,iwplcc)

    !!Jinfeng's grassland management module might should also be put here.

    !! We redistribute the updated litter into four fuel classes, so that
    !! the balance between aboveground litter and fuel is mainted. The subtraction
    !! of fuel burned by land cover change fires from the fuel pool is made here. 
    fuel_all_type(:,:,:,:) = fuel_1hr(:,:,:,:) + fuel_10hr(:,:,:,:) + &
                               fuel_100hr(:,:,:,:) + fuel_1000hr(:,:,:,:)
    fuel_type_frac(:,:,:,:,:) = 0.25
    WHERE(fuel_all_type(:,:,:,:) > min_stomate)
      fuel_type_frac(:,:,:,:,1) = fuel_1hr(:,:,:,:)/fuel_all_type(:,:,:,:)
      fuel_type_frac(:,:,:,:,2) = fuel_10hr(:,:,:,:)/fuel_all_type(:,:,:,:)
      fuel_type_frac(:,:,:,:,3) = fuel_100hr(:,:,:,:)/fuel_all_type(:,:,:,:)
      fuel_type_frac(:,:,:,:,4) = fuel_1000hr(:,:,:,:)/fuel_all_type(:,:,:,:)
    ENDWHERE
    DO j=1,nvm
      fuel_1hr(:,j,:,:) = litter(:,:,j,iabove,:) * fuel_type_frac(:,j,:,:,1) 
      fuel_10hr(:,j,:,:) = litter(:,:,j,iabove,:) * fuel_type_frac(:,j,:,:,2) 
      fuel_100hr(:,j,:,:) = litter(:,:,j,iabove,:) * fuel_type_frac(:,j,:,:,3) 
      fuel_1000hr(:,j,:,:) = litter(:,:,j,iabove,:) * fuel_type_frac(:,j,:,:,4) 
    END DO 
    
  !! 4. history
    
    veget_cov_max(:,:) = veget_cov_max_new(:,:)
    convflux(:,iwplcc)       = convflux(:,iwplcc)/one_year*dt_days
    cflux_prod10(:,iwplcc)    = cflux_prod10(:,iwplcc)/one_year*dt_days
    cflux_prod100(:,iwplcc)   = cflux_prod100(:,iwplcc)/one_year*dt_days

    pool_end_pft(:,:) = SUM(biomass(:,:,:,icarbon),DIM=3) &
                    + SUM(SUM(litter(:,:,:,:,icarbon),DIM=2),DIM=3) &
                    + SUM(carbon(:,:,:),DIM=2) &
                    + SUM(bm_to_litter(:,:,:,icarbon),DIM=3) &
                    + SUM(turnover_daily(:,:,:,icarbon),DIM=3)
                    
    pool_end(:) = SUM(pool_end_pft(:,:)*veget_cov_max(:,:),DIM=2) &
                    + SUM(prod10(:,:,iwplcc),DIM=2) + SUM(prod100(:,:,iwplcc),DIM=2)


    outflux(:) = SUM(SUM(emideforest_biomass_accu(:,:,:,icarbon),DIM=3),DIM=2) &
                 + SUM(SUM(emideforest_litter_accu(:,:,:,icarbon),DIM=3),DIM=2) &
                 + SUM(flux10(:,:,iwplcc),DIM=2) + SUM(flux100(:,:,iwplcc),DIM=2) &
                 - SUM(co2_to_bm(:,:)*veget_cov_max(:,:),DIM=2)

    print *,"pool_start: ",pool_start(:)
    print *,"pool_end: ",pool_end(:)
    print *,"outflux: ",outflux(:)
    print *,"pool_change: ",pool_start(:)-pool_end(:)
    print *,'prod10_end_lcc',prod10(:,:,:)

    deflitsup_total(:,:) = SUM(deforest_litter_surplus(:,:,:),dim=3)
    defbiosup_total(:,:) = SUM(deforest_biomass_surplus(:,:,:),dim=3)

    CALL histwrite (hist_id_stomate, 'dilu_lit_met', itime, &
                    dilu_lit(:,imetabolic,iabove,icarbon), npts, hori_index)
    CALL histwrite (hist_id_stomate, 'dilu_lit_str', itime, &
                    dilu_lit(:,istructural,iabove,icarbon), npts, hori_index)


    CALL histwrite (hist_id_stomate, 'SurpBioLEAF', itime, &
         deforest_biomass_surplus(:,:,ileaf), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'SurpBioRESERVE', itime, &
         deforest_biomass_surplus(:,:,icarbres), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'SurpBioFRUIT', itime, &
         deforest_biomass_surplus(:,:,ifruit), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'SurpBioSapABOVE', itime, &
         deforest_biomass_surplus(:,:,isapabove), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'SurpBioHeartABOVE', itime, &
         deforest_biomass_surplus(:,:,iheartabove), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'SurpBioSapBELOW', itime, &
         deforest_biomass_surplus(:,:,isapbelow), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'SurpBioHeartBELOW', itime, &
         deforest_biomass_surplus(:,:,iheartbelow), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'SurpBioROOT', itime, &
         deforest_biomass_surplus(:,:,iroot), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'SurpLitMET', itime, &
         deforest_litter_surplus(:,:,imetabolic), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'SurpLitSTR', itime, &
         deforest_litter_surplus(:,:,istructural), npts*nvm, horipft_index)

    CALL histwrite (hist_id_stomate, 'DefiBioLEAF', itime, &
         deforest_biomass_deficit(:,:,ileaf), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'DefiBioRESERVE', itime, &
         deforest_biomass_deficit(:,:,icarbres), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'DefiBioFRUIT', itime, &
         deforest_biomass_deficit(:,:,ifruit), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'DefiBioSapABOVE', itime, &
         deforest_biomass_deficit(:,:,isapabove), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'DefiBioHeartABOVE', itime, &
         deforest_biomass_deficit(:,:,iheartabove), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'DefiBioSapBELOW', itime, &
         deforest_biomass_deficit(:,:,isapbelow), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'DefiBioHeartBELOW', itime, &
         deforest_biomass_deficit(:,:,iheartbelow), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'DefiBioROOT', itime, &
         deforest_biomass_deficit(:,:,iroot), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'DefiLitMET', itime, &
         deforest_litter_deficit(:,:,imetabolic), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'DefiLitSTR', itime, &
         deforest_litter_deficit(:,:,istructural), npts*nvm, horipft_index)


    IF (printlev>=4) WRITE(numout,*) 'Leaving lcchange_main'
    
  END SUBROUTINE lcchange_deffire


  !SUBROUTINE lcc_neighbour_shift(ipts,neighbours,veget_cov_max,lcc,veget_cov_max_new)  
  !  INTEGER(i_std), DIMENSION(npts,8), INTENT(in)             :: neighbours      !! indices of the 8 neighbours of each grid point 
  !                                                                               !! (unitless);  
  !                                                                               !! [1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW]  


  !END SUBROUTINE lcc_neighbour_shift

    !print *,'end_biomass',SUM(SUM(biomass(:,:,:,icarbon),DIM=3)*veget_cov_max(:,:),DIM=2)
    !print *,'end_litter',SUM(SUM(SUM(litter(:,:,:,:,icarbon),DIM=2),DIM=3)*veget_cov_max(:,:),DIM=2)
    !print *, 'end_soil',SUM(SUM(carbon(:,:,:),DIM=2)*veget_cov_max(:,:),DIM=2)
    !print *,'end_bm2lit',sum(SUM(bm_to_litter(:,:,:,icarbon),DIM=3)*veget_cov_max(:,:),dim=2)
    !print *,'end_turnover',sum(SUM(turnover_daily(:,:,:,icarbon),DIM=3)*veget_cov_max(:,:),dim=2)
    !print *,'end_prod10', SUM(prod10(:,:),DIM=2)
    !print *,'end_prod100',SUM(prod100(:,:),DIM=2)

!    !!block to check
!    pool_end_pft(:,:) = SUM(biomass(:,:,:,icarbon),DIM=3) &
!                    + SUM(SUM(litter(:,:,:,:,icarbon),DIM=2),DIM=3) &
!                    + SUM(carbon(:,:,:),DIM=2) &
!                    + SUM(bm_to_litter(:,:,:,icarbon),DIM=3) &
!                    + SUM(turnover_daily(:,:,:,icarbon),DIM=3)
!                    
!    pool_end(:) = SUM(pool_end_pft(:,:)*veget_cov_max(:,:),DIM=2) &
!                    + SUM(prod10(:,:),DIM=2) + SUM(prod100(:,:),DIM=2)
!
!    outflux(:) = SUM(SUM(emideforest_biomass_accu(:,:,:,icarbon),DIM=3),DIM=2) &
!                 + SUM(SUM(emideforest_litter_accu(:,:,:,icarbon),DIM=3),DIM=2) &
!                 + SUM(flux10(:,:),DIM=2) + SUM(flux100,DIM=2) &
!                 - SUM(co2_to_bm(:,:)*veget_cov_max(:,:),DIM=2)
!
!    print *,"pool_start: ",pool_start(:)
!    print *,"pool_end: ",pool_end(:)
!    print *,"outflux: ",outflux(:)
!    print *,"pool_change: ",pool_start(:)-pool_end(:)
!    !!end block to check


  

!! ================================================================================================================================
!! SUBROUTINE   : lcchange_main_agripeat
!!
!>\BRIEF        Impact of land cover change on carbon stocks, when peatland modules are activated
!!
!! DESCRIPTION  : This subroutine is activate if DGVM=T during spinups, and then during historical simulations VEGET_UPDATE>0Y and AGRI_PEAT=True
!! lcchange_main_agripeat is called from stomateLpj the first time step after the vegetation map has been changed. 
!! After the vegetation map has been read, natural peatland area occupy each PFT in proportion to grid cell area of the PFT.
!! Crops occupied by peatland become new PFTs (PFT15, PFT16) (CALL SUBROUTINE agripeat_adjust_fractions).
!! The impact of land cover change on carbon stocks is computed in this subroutine. 
!! On the basis of this difference, the amount of 'new establishment'/'biomass export',and increase/decrease of each component, are estimated.
!!
!_ ================================================================================================================================
  SUBROUTINE lcchange_main_agripeat ( npts, dt_days, veget_cov_max_old, veget_cov_max_new, &
       biomass, ind, age, PFTpresent, senescence, when_growthinit, everywhere, &
       co2_to_bm, bm_to_litter, turnover_daily, bm_sapl, cn_ind,flux10,flux100, &
       prod10, prod100, &
       convflux, &
       cflux_prod10, cflux_prod100, leaf_frac,&
       npp_longterm, lm_lastyearmax, litter, litter_avail, litter_not_avail, &
       carbon, &
       deepC_a, deepC_s, deepC_p, &
       fuel_1hr, fuel_10hr, fuel_100hr, fuel_1000hr, &
       veget_cov_max_adjusted,lalo,carbon_save,deepC_a_save,deepC_s_save,deepC_p_save,delta_fsave)
!,biomass_remove)

    IMPLICIT NONE

  !! 0. Variable and parameter declaration 

    !! 0.1 Input variables

    INTEGER, INTENT(in)                                       :: npts             !! Domain size - number of pixels (unitless)
    REAL(r_std), INTENT(in)                                   :: dt_days          !! Time step of vegetation dynamics for stomate
                                                                                  !! (days)
    REAL(r_std), DIMENSION(nvm, nparts,nelements), INTENT(in) :: bm_sapl          !! biomass of sapling 
                                                                                  !! @tex ($gC individual^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: veget_cov_max_old!! Current "maximal" coverage fraction of a PFT (LAI
                                                                                  !! -> infinity) on ground
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: veget_cov_max_new!! New "maximal" coverage fraction of a PFT (LAI ->
                                                                                  !! infinity) on ground (unitless) 

    REAL(r_std),DIMENSION(npts,2),INTENT(in)                   :: lalo
 
    !! 0.2 Output variables

    REAL(r_std), DIMENSION(:,:), INTENT(out)                  :: convflux         !! release during first year following land cover
                                                                                  !! change  (npts, nwp)
    REAL(r_std), DIMENSION(:,:), INTENT(out)                  :: cflux_prod10     !! total annual release from the 10 year-turnover
                                                                                  !! pool @tex ($gC m^{-2}$) @endtex
                                                                                  !! change  (npts, nwp)
    REAL(r_std), DIMENSION(:,:), INTENT(out)                  :: cflux_prod100    !! total annual release from the 100 year-
                                                                                  !! turnover pool @tex ($gC m^{-2}$) @endtex
                                                                                  !! change  (npts, nwp)
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout):: turnover_daily   !! Turnover rates 
                                                                                      !! @tex ($gC m^{-2} day^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm),INTENT(out)              :: veget_cov_max_adjusted    !! fraction of vegetation after adjustment


    !! 0.3 Modified variables   

    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout):: biomass    !! biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: ind              !! Number of individuals @tex ($m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: age              !! mean age (years)
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)               :: senescence       !! plant senescent (only for deciduous trees) Set
                                                                                  !! to .FALSE. if PFT is introduced or killed
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)               :: PFTpresent       !! Is pft there (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: everywhere       !! is the PFT everywhere in the grid box or very 
                                                                                  !! localized (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: when_growthinit  !! how many days ago was the beginning of the 
                                                                                  !! growing season (days)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: co2_to_bm        !! biomass uptaken 
                                                                                  !! @tex ($gC m^{-2} day^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: bm_to_litter !! conversion of biomass to litter 
                                                                                  !! @tex ($gC m^{-2} day^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: cn_ind           !! crown area of individuals 
                                                                                  !! @tex ($m^{2}$) @endtex
    REAL(r_std), DIMENSION(npts,0:10,nwp), INTENT(inout)              :: prod10           !! products remaining in the 10 year-turnover
                                                                                  !! pool after the annual release for each 
                                                                                  !! compartment (10 + 1 : input from year of land
                                                                                  !! cover change)  (npts,0:10,nwp)
    REAL(r_std), DIMENSION(npts,0:100,nwp), INTENT(inout)              :: prod100          !! products remaining in the 100 year-turnover
                                                                                  !! pool after the annual release for each 
                                                                                  !! compartment (100 + 1 : input from year of land
                                                                                  !! cover change) (npts,0:100,nwp)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)              :: flux10           !! annual release from the 10/100 year-turnover 
                                                                                  !! pool compartments (npts,10,nwp)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)              :: flux100          !! annual release from the 100/100 year-turnover
                                                                                  !! pool compartments (npts,100,nwp)
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout) :: leaf_frac        !! fraction of leaves in leaf age class 
                                                                                  !! (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: lm_lastyearmax   !! last year's maximum leaf mass for each PFT 
                                                                                  !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: npp_longterm     !! "long term" net primary productivity 
                                                                                  !! @tex ($gC m^{-2} year^{-1}$) @endtex
    REAL(r_std),DIMENSION(npts,nlitt,nvm,nlevs,nelements), INTENT(inout):: litter !! metabolic and structural litter, above and 
                                                                                  !! below ground @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nlitt,nvm), INTENT(inout):: litter_avail
    REAL(r_std), DIMENSION(npts,nlitt,nvm) , INTENT(inout):: litter_not_avail
    REAL(r_std),DIMENSION(npts,ncarb,nvm), INTENT(inout)      :: carbon           !! carbon pool: active, slow, or passive 
                                                                                  !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,ndeep,nvm), INTENT(inout)         :: deepC_a      !! Permafrost soil carbon (g/m**3) active
    REAL(r_std), DIMENSION(npts,ndeep,nvm), INTENT(inout)         :: deepC_s      !! Permafrost soil carbon (g/m**3) slow
    REAL(r_std), DIMENSION(npts,ndeep,nvm), INTENT(inout)         :: deepC_p      !! Permafrost soil carbon (g/m**3) passive
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements),INTENT(inout)        :: fuel_1hr
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements),INTENT(inout)        :: fuel_10hr
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements),INTENT(inout)        :: fuel_100hr
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements),INTENT(inout)        :: fuel_1000hr

    REAL(r_std), DIMENSION(npts,ncarb,nvm), INTENT(inout)     :: carbon_save
    REAL(r_std), DIMENSION(npts,ndeep), INTENT(inout)         :: deepC_a_save
    REAL(r_std), DIMENSION(npts,ndeep), INTENT(inout)         :: deepC_s_save
    REAL(r_std), DIMENSION(npts,ndeep), INTENT(inout)         :: deepC_p_save
    REAL(r_std), DIMENSION(npts), INTENT(inout)               :: delta_fsave

   ! REAL(r_std),DIMENSION(npts,nvm,nparts,nelements), INTENT(out) :: biomass_remove

    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements)       :: fuel_all_type
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements,4)     :: fuel_type_frac

    !! 0.4 Local variables

    INTEGER(i_std)                                            :: i, j, k, l, m, n    !! indices (unitless)
    REAL(r_std),DIMENSION(npts,nelements)                     :: bm_new           !! biomass increase @tex ($gC m^{-2}$) @endtex
    REAL(r_std),DIMENSION(npts,nparts,nelements)              :: biomass_loss     !! biomass loss @tex ($gC m^{-2}$) @endtex
    REAL(r_std),DIMENSION(npts,nparts,nelements)              :: biomass_loss_peat
    REAL(r_std)                                               :: above            !! aboveground biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std),DIMENSION(npts,nlitt,nlevs,nelements)         :: dilu_lit         !! Litter dilution @tex ($gC m^{-2}$) @endtex
    REAL(r_std),DIMENSION(npts,nlitt,nlevs,nelements)         :: dilu_lit_peat
    REAL(r_std),DIMENSION(npts,ncarb)                         :: dilu_soil_carbon !! Soil Carbondilution @tex ($gC m^{-2}$) @endtex
    REAL(r_std),DIMENSION(npts,ndeep,ncarb)                   :: dilu_soil_carbon_vertres !!vertically-resolved Soil Carbondilution (gC/m²)
    REAL(r_std),DIMENSION(npts)                               :: delta_nat_peat
    REAL(r_std),DIMENSION(nvm)                                :: delta_veg        !! changes in "maximal" coverage fraction of PFT 
    REAL(r_std)                                               :: delta_veg_sum    !! sum of delta_veg, vegetations
    REAL(r_std),DIMENSION(npts,nvm)                           :: delta_ind        !! change in number of individuals  

    REAL(r_std), DIMENSION(npts)                              :: soilc_before
    REAL(r_std), DIMENSION(npts)                              :: soilc_after
    REAL(r_std), DIMENSION(npts,ndeep)                        :: delta_a      !! Permafrost soil carbon (g/m**3) active
    REAL(r_std), DIMENSION(npts,ndeep)                        :: delta_s      !! Permafrost soil carbon (g/m**3) slow
    REAL(r_std), DIMENSION(npts,ndeep)                        :: delta_p      !!  
    REAL(r_std),DIMENSION(npts,ncarb)                         :: delta_carbon
    ! peatland PFT number in nvm
    INTEGER(i_std)                                            :: c3crop, c4crop, peat_grass, peat_c3crop, peat_c4crop, peat_forest 
!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering lcchange_main_agripeat'

    soilc_before(:)=zero
    soilc_after(:)=zero
   
    DO i=1, npts
       DO j=1,nvm
           DO l=1,ncarb
              soilc_before(i)=soilc_before(i)+carbon(i,l,j)*veget_cov_max_old(i,j)
           ENDDO
       ENDDO 
    ENDDO

    
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

    ! adjust vegetation fractions with different assumptions on the peatland
    ! occupation by new natural and croplands
    CALL agripeat_adjust_fractions (npts, lalo,veget_cov_max_new, veget_cov_max_old, veget_cov_max_adjusted)

  !! 1. initialization
    
    prod10(:,0,:)         = zero
    prod100(:,0,:)        = zero   
    above               = zero
    convflux(:,:)         = zero
    cflux_prod10(:,:)     = zero
    cflux_prod100(:,:)    = zero
    delta_ind(:,:)      = zero
    delta_veg(:)        = zero
    dilu_soil_carbon_vertres(:,:,:) = zero

    delta_nat_peat(:)  = zero
    delta_a(:,:) =zero
    delta_s(:,:) =zero
    delta_p(:,:) =zero
    delta_carbon(:,:) =zero

  !! 3. calculation of changes in carbon stocks and biomass by land cover change\n

    DO i = 1, npts ! Loop over # pixels - domain size

       !! 3.1 initialization of carbon stocks\n
       delta_veg(:) = veget_cov_max_adjusted(i,:)-veget_cov_max_old(i,:)
       ! will be calculated later
       delta_veg_sum = zero
       !delta_veg_sum = SUM(delta_veg,MASK=delta_veg.LT.0.)
 
       dilu_lit(i,:,:,:) = zero
       dilu_soil_carbon(i,:) = zero
       biomass_loss(i,:,:) = zero

       ! changes in natural peatland
       DO j=2,nvm
         IF (is_peat(j)) THEN
           delta_nat_peat(i) = delta_nat_peat(i) + delta_veg(j)
         ENDIF
       ENDDO
       dilu_lit_peat(i,:,:,:) = zero
       biomass_loss_peat(i,:,:)=zero
       ! biomass_remove(i,:,:,:)=zero

       !! With peatland, the dilution need to be do for peatland and non-peatland seperately
       ! !! 3.2 if vegetation coverage decreases, compute dilution of litter, soil carbon, and biomass.\n
       ! DO j=2, nvm
       !    IF ( delta_veg(j) < -min_stomate ) THEN
       !       dilu_lit(i,:,:,:) = dilu_lit(i,:,:,:) + delta_veg(j)*litter(i,:,j,:,:) / delta_veg_sum
       !       biomass_loss(i,:,:) = biomass_loss(i,:,:) + biomass(i,j,:,:)*delta_veg(j) / delta_veg_sum
       !       IF ( ok_pc ) THEN
       !              dilu_soil_carbon_vertres(i,:,iactive)=dilu_soil_carbon_vertres(i,:,iactive) + &
       !                   delta_veg(j) * deepC_a(i,:,j) / delta_veg_sum
       !              dilu_soil_carbon_vertres(i,:,islow)=dilu_soil_carbon_vertres(i,:,islow) + &
       !                   delta_veg(j) * deepC_s(i,:,j) / delta_veg_sum
       !              dilu_soil_carbon_vertres(i,:,ipassive)=dilu_soil_carbon_vertres(i,:,ipassive) + &
       !                   delta_veg(j) * deepC_p(i,:,j) / delta_veg_sum
       !       ELSE
       !              dilu_soil_carbon(i,:) =  dilu_soil_carbon(i,:) + delta_veg(j) * carbon(i,:,j) / delta_veg_sum
       !       ENDIF
       !    ENDIF
       ! ENDDO

       !! 3.3
       !!expanding vegetations (peat, and non-peat vegetations), biomass
       DO j=2,nvm

          !! 3.3.1 The case that vegetation coverage of PFTj increases
          IF ( delta_veg(j) > min_stomate) THEN

             !! 3.3.1.1 Initial setting of new establishment
             IF (veget_cov_max_old(i,j) .LT. min_stomate) THEN 
                IF (is_tree(j)) THEN

                   ! cn_sapl(j)=0.5; stomate_data.f90
                   cn_ind(i,j) = cn_sapl(j) 
                ELSE
                   cn_ind(i,j) = un
                ENDIF
                ind(i,j)= delta_veg(j) / cn_ind(i,j)
                PFTpresent(i,j) = .TRUE.
                everywhere(i,j) = 1.
                senescence(i,j) = .FALSE.
                age(i,j) = zero

                ! large_value = 1.E33_r_std
                when_growthinit(i,j) = large_value 
                leaf_frac(i,j,1) = 1.0
                npp_longterm(i,j) = npp_longterm_init
                lm_lastyearmax(i,j) = bm_sapl(j,ileaf,icarbon) * ind(i,j)
             ENDIF
             IF ( cn_ind(i,j) > min_stomate ) THEN
                delta_ind(i,j) = delta_veg(j) / cn_ind(i,j) 
             ENDIF
             
             !! 3.3.1.2 Update of biomass in each each carbon stock component 
             !!         Update of biomass in each each carbon stock component (leaf, sapabove, sapbelow,
             !>         heartabove, heartbelow, root, fruit, and carbres)\n
             DO k = 1, nparts ! loop over # carbon stock components, nparts = 8; stomate_constant.f90 
                DO l = 1,nelements ! loop over # elements

                   bm_new(i,l) = delta_ind(i,j) * bm_sapl(j,k,l)
                   IF (veget_cov_max_old(i,j) .GT. min_stomate) THEN
                     ! in the case that bm_new is overestimated compared with biomass?
                     IF ((bm_new(i,l)/delta_veg(j)) > biomass(i,j,k,l)) THEN
                        bm_new(i,l) = biomass(i,j,k,l)*delta_veg(j)
                     ENDIF
                  ENDIF
                  biomass(i,j,k,l) = ( biomass(i,j,k,l) * veget_cov_max_old(i,j) + bm_new(i,l) ) / veget_cov_max_adjusted(i,j)
                  co2_to_bm(i,j) = co2_to_bm(i,j) + (bm_new(i,icarbon)*dt_days) / (one_year * veget_cov_max_adjusted(i,j))
                END DO ! loop over # elements
             ENDDO ! loop over # carbon stock components

             ! With peatland, dilution need to be calculated separately 
             !! 3.3.1.3 Calculation of dilution in litter, soil carbon, and input of litter
             !!        In this 'IF statement', dilu_* is zero. Formulas for litter and soil carbon
             !!         could be shortend?? Are the following formulas correct?


         !ELSE ! delta_veg <= 0

         ENDIF ! End if PFT's coverage reduction
       ENDDO !  Loop over # PFTs

       !! when natural peatland decrease, the only possibility is to be drained
       !for croplands
       !! 2.2 Crops on peatland increase, natural peatland decrease.
       IF (delta_nat_peat(i) < -min_stomate) THEN
          !! crops on peatland obtain soilC from shrinking natural peat
          DO j=1,nvm
             IF (is_peat(j)) THEN ! natural northern peatland (grass) or tropical peatland (forest)
                !dead biomass due to agricultural use is removed from the field
                ! biomass_remove(i,j,:,:) =  biomass(i,j,:,:)*delta_veg(j)*0.75
                dilu_lit_peat(i,:,:,:) = dilu_lit_peat(i,:,:,:) - delta_veg(j)*litter(i,:,j,:,:) 
                biomass_loss_peat(i,:,:) = biomass_loss_peat(i,:,:) - biomass(i,j,:,:)*delta_veg(j)
                delta_a(i,:)=delta_a(i,:)-deepC_a(i,:,j)*delta_veg(j)
                delta_s(i,:)=delta_s(i,:)-deepC_s(i,:,j)*delta_veg(j)
                delta_p(i,:)=delta_p(i,:)-deepC_p(i,:,j)*delta_veg(j)
                delta_carbon(i,:)= delta_carbon(i,:)-carbon(i,:,j) * delta_veg(j)
                ! Biomass export
                above = biomass(i,j,isapabove,icarbon) + biomass(i,j,iheartabove,icarbon)
                convflux(i,iwplcc)  = convflux(i,iwplcc)  - ( coeff_lcchange_1(j) * above *delta_veg(j) )
                prod10(i,0,iwplcc)  = prod10(i,0,iwplcc)  - ( coeff_lcchange_10(j) * above *delta_veg(j) )
                prod100(i,0,iwplcc) = prod100(i,0,iwplcc) - ( coeff_lcchange_100(j) * above *delta_veg(j) )
             ENDIF 
          ENDDO
          IF ((delta_veg(peat_c3crop) > min_stomate) .AND. (delta_veg(peat_c4crop) > min_stomate)) THEN        
              !! both PFT15 and PFT16 increase, give C from natural peat to them proportionally 
              DO n=2,nvm 
                IF (n .EQ. peat_c3crop .OR. n .EQ. peat_c4crop) THEN   
                  litter(i,:,n,:,:)=(litter(i,:,n,:,:) * veget_cov_max_old(i,n) - &
                        & dilu_lit_peat(i,:,:,:) * (delta_veg(n)/delta_nat_peat(i))) / veget_cov_max_adjusted(i,n)
                  DO l = 1,nelements
                      bm_to_litter(i,n,isapbelow,l) = (bm_to_litter(i,n,isapbelow,l)*veget_cov_max_old(i,n) - &
                           & biomass_loss_peat(i,isapbelow,l)*delta_veg(n)/delta_nat_peat(i)) / veget_cov_max_adjusted(i,n)
                      bm_to_litter(i,n,iheartbelow,l) =(bm_to_litter(i,n,iheartbelow,l)*veget_cov_max_old(i,n)-& 
                            & biomass_loss_peat(i,iheartbelow,l) *delta_veg(n)/delta_nat_peat(i))/ veget_cov_max_adjusted(i,n)
                      bm_to_litter(i,n,iroot,l) =(bm_to_litter(i,n,iroot,l)*veget_cov_max_old(i,n) - &
                            & biomass_loss_peat(i,iroot,l)*delta_veg(n)/delta_nat_peat(i))/ veget_cov_max_adjusted(i,n)
                      bm_to_litter(i,n,ifruit,l) = (bm_to_litter(i,n,ifruit,l)*veget_cov_max_old(i,n) - & 
                            & biomass_loss_peat(i,ifruit,l)*delta_veg(n)/delta_nat_peat(i))/ veget_cov_max_adjusted(i,n)
                      bm_to_litter(i,n,icarbres,l) = (bm_to_litter(i,n,icarbres,l)*veget_cov_max_old(i,n) - &
                            &  biomass_loss_peat(i,icarbres,l) *delta_veg(n)/delta_nat_peat(i)) / veget_cov_max_adjusted(i,n)
                      bm_to_litter(i,n,ileaf,l) =(bm_to_litter(i,n,ileaf,l)*veget_cov_max_old(i,n) - & 
                            & biomass_loss_peat(i,ileaf,l)* delta_veg(n)/delta_nat_peat(i)) / veget_cov_max_adjusted(i,n)
                  ENDDO
                  IF (ok_pc) THEN
                     deepC_a(i,:,n)=(deepC_a(i,:,n) * veget_cov_max_old(i,n) - &
                             delta_a(i,:)*delta_veg(n)/delta_nat_peat(i)) / veget_cov_max_adjusted(i,n)
                     deepC_s(i,:,n)=(deepC_s(i,:,n) * veget_cov_max_old(i,n) - &
                             delta_s(i,:)*delta_veg(n)/delta_nat_peat(i)) / veget_cov_max_adjusted(i,n)
                     deepC_p(i,:,n)=(deepC_p(i,:,n) * veget_cov_max_old(i,n) - &
                             delta_p(i,:)*delta_veg(n)/delta_nat_peat(i)) /veget_cov_max_adjusted(i,n)     
                  ENDIF 
                     carbon(i,:,n)=(carbon(i,:,n) * veget_cov_max_old(i,n) - delta_carbon(i,:) * delta_veg(n)/delta_nat_peat(i))/veget_cov_max_adjusted(i,n)
                ENDIF
              ENDDO      
          ELSE 
              !! one of PFT15 and PFT16 increase, the other one decrease
              !! the increasing one obtain C from the decreasing one and from shrinking nature peat
              DO n=2,nvm
                IF (n .EQ. peat_c3crop .OR. n .EQ. peat_c4crop) THEN
                  IF (delta_veg(n) .LE. 0.) THEN
                    ! biomass_remove(i,n,:,:)= biomass(i,n,:,:)*delta_veg(n)*0.75
                     biomass_loss_peat(i,:,:) = biomass_loss_peat(i,:,:) - biomass(i,n,:,:)*delta_veg(n)
                     dilu_lit_peat(i,:,:,:) = dilu_lit_peat(i,:,:,:) - delta_veg(n)*litter(i,:,n,:,:)
                     delta_a(i,:)=delta_a(i,:)-deepC_a(i,:,n)*delta_veg(n)
                     delta_s(i,:)=delta_s(i,:)-deepC_s(i,:,n)*delta_veg(n)
                     delta_p(i,:)=delta_p(i,:)-deepC_p(i,:,n)*delta_veg(n)
                     delta_carbon(i,:)= delta_carbon(i,:)-carbon(i,:,n) * delta_veg(n)
                     above = biomass(i,n,isapabove,icarbon) + biomass(i,n,iheartabove,icarbon)
                     convflux(i,iwplcc)  = convflux(i,iwplcc)  - ( coeff_lcchange_1(n) * above *delta_veg(n) )
                     prod10(i,0,iwplcc)  = prod10(i,0,iwplcc)  - ( coeff_lcchange_10(n) * above *delta_veg(n) )
                     prod100(i,0,iwplcc) = prod100(i,0,iwplcc) - ( coeff_lcchange_100(n) * above*delta_veg(n) )
                  ENDIF
                ENDIF
              ENDDO

              DO n=2,nvm
                IF (n .EQ. peat_c3crop .OR. n .EQ. peat_c4crop) THEN
                  IF (delta_veg(n) .GT. 0.) THEN   
                    litter(i,:,n,:,:)=(litter(i,:,n,:,:) * veget_cov_max_old(i,n) + &
                       dilu_lit_peat(i,:,:,:)) / veget_cov_max_adjusted(i,n)
                    DO l = 1,nelements
                       bm_to_litter(i,n,isapbelow,l) = (bm_to_litter(i,n,isapbelow,l)*veget_cov_max_old(i,n) + &
                          & biomass_loss_peat(i,isapbelow,l))/ veget_cov_max_adjusted(i,n)
                       bm_to_litter(i,n,iheartbelow,l) = (bm_to_litter(i,n,iheartbelow,l)*veget_cov_max_old(i,n) +&
                          & biomass_loss_peat(i,iheartbelow,l))/ veget_cov_max_adjusted(i,n)
                       bm_to_litter(i,n,iroot,l) = (bm_to_litter(i,n,iroot,l)*veget_cov_max_old(i,n) + &
                          & biomass_loss_peat(i,iroot,l))/ veget_cov_max_adjusted(i,n)
                       bm_to_litter(i,n,ifruit,l) = (bm_to_litter(i,n,ifruit,l)*veget_cov_max_old(i,n) + & 
                          & biomass_loss_peat(i,ifruit,l))/ veget_cov_max_adjusted(i,n)
                       bm_to_litter(i,n,icarbres,l) = (bm_to_litter(i,n,icarbres,l)*veget_cov_max_old(i,n) + &
                          &  biomass_loss_peat(i,icarbres,l))/ veget_cov_max_adjusted(i,n)
                       bm_to_litter(i,n,ileaf,l) = (bm_to_litter(i,n,ileaf,l)*veget_cov_max_old(i,n) + &
                          & biomass_loss_peat(i,ileaf,l))/ veget_cov_max_adjusted(i,n)
                    END DO
                    IF (ok_pc) THEN
                       deepC_a(i,:,n)=(deepC_a(i,:,n) * veget_cov_max_old(i,n) + delta_a(i,:)) / veget_cov_max_adjusted(i,n)
                       deepC_s(i,:,n)=(deepC_s(i,:,n) * veget_cov_max_old(i,n) + delta_s(i,:)) / veget_cov_max_adjusted(i,n)
                       deepC_p(i,:,n)=(deepC_p(i,:,n) * veget_cov_max_old(i,n) + delta_p(i,:)) /veget_cov_max_adjusted(i,n)
                    ENDIF
                       carbon(i,:,n)=(carbon(i,:,n) * veget_cov_max_old(i,n) + delta_carbon(i,:)) / veget_cov_max_adjusted(i,n)  
                  ENDIF
                ENDIF 
              ENDDO
          ENDIF ! Crops on peatland increase, natural peatland decrease.

          !!non-peat vegetations, dilute within non-peat vegetations
          DO j=1,nvm
            IF (.NOT. is_peat(j) .AND. .NOT. is_croppeat(j)) THEN
              IF ( delta_veg(j) .LT. zero ) THEN
                 delta_veg_sum = delta_veg_sum+delta_veg(j)
              ENDIF
            ENDIF
          ENDDO
          DO j=1,nvm
            IF (.NOT. is_peat(j) .AND. .NOT. is_croppeat(j)) THEN
              IF ( delta_veg(j) < -min_stomate ) THEN
                 dilu_lit(i,:,:,:) = dilu_lit(i,:,:,:) + delta_veg(j)*litter(i,:,j,:,:) / delta_veg_sum
                 biomass_loss(i,:,:) = biomass_loss(i,:,:) + biomass(i,j,:,:)*delta_veg(j) / delta_veg_sum
                 IF ( ok_pc ) THEN
                      dilu_soil_carbon_vertres(i,:,iactive)=dilu_soil_carbon_vertres(i,:,iactive) + &
                           delta_veg(j) * deepC_a(i,:,j) / delta_veg_sum
                      dilu_soil_carbon_vertres(i,:,islow)=dilu_soil_carbon_vertres(i,:,islow) + &
                           delta_veg(j) * deepC_s(i,:,j) / delta_veg_sum
                      dilu_soil_carbon_vertres(i,:,ipassive)=dilu_soil_carbon_vertres(i,:,ipassive) + &
                           delta_veg(j) * deepC_p(i,:,j) / delta_veg_sum
                 ENDIF
                 dilu_soil_carbon(i,:) =  dilu_soil_carbon(i,:) + delta_veg(j) * carbon(i,:,j) / delta_veg_sum
              ENDIF
            ENDIF
          ENDDO
          DO j=1,nvm
            IF (.NOT. is_peat(j) .AND. .NOT. is_croppeat(j)) THEN
              IF ( delta_veg(j) > min_stomate) THEN
                 litter(i,:,j,:,:)=(litter(i,:,j,:,:) * veget_cov_max_old(i,j) + &
                        dilu_lit(i,:,:,:) * delta_veg(j)) / veget_cov_max_adjusted(i,j)
                 IF (is_grassland_manag(j) .AND. is_grassland_grazed(j)) THEN
                   litter_avail(i,:,j) = litter_avail(i,:,j) * veget_cov_max_old(i,j) / veget_cov_max_adjusted(i,j)
                   litter_not_avail(i,:,j) = litter(i,:,j,iabove,icarbon) - litter_avail(i,:,j)
                 ENDIF
                IF ( ok_pc ) THEN
                   deepC_a(i,:,j)=(deepC_a(i,:,j) * veget_cov_max_old(i,j) + &
                      dilu_soil_carbon_vertres(i,:,iactive) * delta_veg(j)) /veget_cov_max_adjusted(i,j)
                   deepC_s(i,:,j)=(deepC_s(i,:,j) * veget_cov_max_old(i,j) + &
                      dilu_soil_carbon_vertres(i,:,islow) * delta_veg(j)) /veget_cov_max_adjusted(i,j)
                   deepC_p(i,:,j)=(deepC_p(i,:,j) * veget_cov_max_old(i,j) + &
                      dilu_soil_carbon_vertres(i,:,ipassive) * delta_veg(j)) /veget_cov_max_adjusted(i,j)
                ENDIF
                carbon(i,:,j)=(carbon(i,:,j) * veget_cov_max_old(i,j) +dilu_soil_carbon(i,:) * delta_veg(j)) / veget_cov_max_adjusted(i,j)

                DO l = 1,nelements
                   bm_to_litter(i,j,isapbelow,l) = (bm_to_litter(i,j,isapbelow,l)*veget_cov_max_old(i,j) + &
                      &  biomass_loss(i,isapbelow,l)*delta_veg(j)) / veget_cov_max_adjusted(i,j)
                   bm_to_litter(i,j,iheartbelow,l) = (bm_to_litter(i,j,iheartbelow,l)*veget_cov_max_old(i,j) + &
                      &  biomass_loss(i,iheartbelow,l) *delta_veg(j))/ veget_cov_max_adjusted(i,j) 
                   bm_to_litter(i,j,iroot,l) = (bm_to_litter(i,j,iroot,l)*veget_cov_max_old(i,j) + &
                      &  biomass_loss(i,iroot,l)*delta_veg(j)) / veget_cov_max_adjusted(i,j)
                   bm_to_litter(i,j,ifruit,l) = (bm_to_litter(i,j,ifruit,l)*veget_cov_max_old(i,j) + &
                      &  biomass_loss(i,ifruit,l)*delta_veg(j)) / veget_cov_max_adjusted(i,j)
                   bm_to_litter(i,j,icarbres,l) = (bm_to_litter(i,j,icarbres,l)*veget_cov_max_old(i,j) + &
                      &  biomass_loss(i,icarbres,l) *delta_veg(j)) / veget_cov_max_adjusted(i,j)
                   bm_to_litter(i,j,ileaf,l) = (bm_to_litter(i,j,ileaf,l)*veget_cov_max_old(i,j) + &
                      &  biomass_loss(i,ileaf,l)*delta_veg(j))/ veget_cov_max_adjusted(i,j)
                ENDDO

                age(i,j)=age(i,j)*veget_cov_max_old(i,j)/veget_cov_max_adjusted(i,j)
              ELSE
                above = biomass(i,j,isapabove,icarbon) + biomass(i,j,iheartabove,icarbon)
                convflux(i,iwplcc)  = convflux(i,iwplcc)  - ( coeff_lcchange_1(j) * above * delta_veg(j) )
                prod10(i,0,iwplcc)  = prod10(i,0,iwplcc)  - ( coeff_lcchange_10(j) * above *delta_veg(j) )
                prod100(i,0,iwplcc) = prod100(i,0,iwplcc) - ( coeff_lcchange_100(j) * above *delta_veg(j) )
              ENDIF ! End if PFT's coverage reduction
            ENDIF
          ENDDO
       ELSE
         !!2.3 Crops on peatland decrease,abandoned agricultural peat become
         !natural vegetations, natural peatland will not change
         delta_veg_sum = SUM(delta_veg,MASK=delta_veg.LT. 0.)
         !! if vegetation coverage decreases, compute dilution of litter, soil carbon, and biomass.\n
         DO j=1, nvm
           IF ( delta_veg(j) < -min_stomate ) THEN
              dilu_lit(i,:,:,:) = dilu_lit(i,:,:,:) + delta_veg(j)*litter(i,:,j,:,:) / delta_veg_sum
              biomass_loss(i,:,:) = biomass_loss(i,:,:) + biomass(i,j,:,:)*delta_veg(j) / delta_veg_sum
              IF ( ok_pc ) THEN
                 dilu_soil_carbon_vertres(i,:,iactive)=dilu_soil_carbon_vertres(i,:,iactive) + &
                     delta_veg(j) * deepC_a(i,:,j) / delta_veg_sum
                 dilu_soil_carbon_vertres(i,:,islow)=dilu_soil_carbon_vertres(i,:,islow) + &
                     delta_veg(j) * deepC_s(i,:,j) / delta_veg_sum
                 dilu_soil_carbon_vertres(i,:,ipassive)=dilu_soil_carbon_vertres(i,:,ipassive) + &
                     delta_veg(j) * deepC_p(i,:,j) / delta_veg_sum
              ENDIF
              dilu_soil_carbon(i,:) =  dilu_soil_carbon(i,:) + delta_veg(j) * carbon(i,:,j) / delta_veg_sum
           ENDIF
         ENDDO
         DO j=1, nvm ! Loop over # PFTs
          !! The case that vegetation coverage of PFTj increases
          IF ( delta_veg(j) > min_stomate) THEN
             litter(i,:,j,:,:)=(litter(i,:,j,:,:) * veget_cov_max_old(i,j) + &
                  dilu_lit(i,:,:,:) * delta_veg(j)) / veget_cov_max_adjusted(i,j)
               IF (is_grassland_manag(j) .AND. is_grassland_grazed(j)) THEN
                 litter_avail(i,:,j) = litter_avail(i,:,j) * veget_cov_max_old(i,j) / veget_cov_max_adjusted(i,j)
                 litter_not_avail(i,:,j) = litter(i,:,j,iabove,icarbon) - litter_avail(i,:,j)
               ENDIF
             IF ( ok_pc ) THEN
                deepC_a(i,:,j)=(deepC_a(i,:,j) * veget_cov_max_old(i,j) + &
                     dilu_soil_carbon_vertres(i,:,iactive) * delta_veg(j)) / veget_cov_max_adjusted(i,j)
                deepC_s(i,:,j)=(deepC_s(i,:,j) * veget_cov_max_old(i,j) + &
                     dilu_soil_carbon_vertres(i,:,islow) * delta_veg(j)) / veget_cov_max_adjusted(i,j)
                deepC_p(i,:,j)=(deepC_p(i,:,j) * veget_cov_max_old(i,j) + &
                     dilu_soil_carbon_vertres(i,:,ipassive) * delta_veg(j)) / veget_cov_max_adjusted(i,j)
             ENDIF
             carbon(i,:,j)=(carbon(i,:,j) * veget_cov_max_old(i,j) + dilu_soil_carbon(i,:) * delta_veg(j)) / veget_cov_max_adjusted(i,j)

             DO l = 1,nelements
                bm_to_litter(i,j,isapbelow,l) = (bm_to_litter(i,j,isapbelow,l)* veget_cov_max_old(i,j) + &
                     &                         biomass_loss(i,isapbelow,l)*delta_veg(j))/ veget_cov_max_adjusted(i,j)
                bm_to_litter(i,j,iheartbelow,l) = (bm_to_litter(i,j,iheartbelow,l)* veget_cov_max_old(i,j) +&
                     &                          biomass_loss(i,iheartbelow,l) *delta_veg(j))/ veget_cov_max_adjusted(i,j)
                bm_to_litter(i,j,iroot,l) = (bm_to_litter(i,j,iroot,l)* veget_cov_max_old(i,j) + &
                     &                          biomass_loss(i,iroot,l)*delta_veg(j))/ veget_cov_max_adjusted(i,j)
                bm_to_litter(i,j,ifruit,l) =(bm_to_litter(i,j,ifruit,l)* veget_cov_max_old(i,j) + &
                     &                          biomass_loss(i,ifruit,l)*delta_veg(j)) / veget_cov_max_adjusted(i,j)
                bm_to_litter(i,j,icarbres,l) = (bm_to_litter(i,j,icarbres,l)* veget_cov_max_old(i,j) + &
                     &                         biomass_loss(i,icarbres,l) *delta_veg(j)) / veget_cov_max_adjusted(i,j)
                bm_to_litter(i,j,ileaf,l) = (bm_to_litter(i,j,ileaf,l)* veget_cov_max_old(i,j) + & 
                     &                         biomass_loss(i,ileaf,l)*delta_veg(j))/ veget_cov_max_adjusted(i,j)
             END DO
             age(i,j)=age(i,j)*veget_cov_max_old(i,j)/veget_cov_max_adjusted(i,j)
          ELSE
             above = biomass(i,j,isapabove,icarbon) + biomass(i,j,iheartabove,icarbon)
             convflux(i,iwplcc)  = convflux(i,iwplcc)  - ( coeff_lcchange_1(j) * above * delta_veg(j) )
             prod10(i,0,iwplcc)  = prod10(i,0,iwplcc)  - ( coeff_lcchange_10(j) * above * delta_veg(j) )
             prod100(i,0,iwplcc) = prod100(i,0,iwplcc) - ( coeff_lcchange_100(j) * above * delta_veg(j) )
          ENDIF ! End if PFT's coverage reduction
         ENDDO ! Loop over # PFTs
       ENDIF
 
       DO j=2,nvm ! Loop over PFTs
           !! If the vegetation is to small, it has been set to 0.
          IF (veget_cov_max_adjusted(i,j) .LT. min_stomate ) THEN

            ind(i,j) = zero
            biomass(i,j,:,:) = zero
            PFTpresent(i,j) = .FALSE.
            senescence(i,j) = .FALSE.
            age(i,j) = zero
            when_growthinit(i,j) = undef
            everywhere(i,j) = zero
            carbon(i,:,j) = zero
            litter(i,:,j,:,:) = zero
            bm_to_litter(i,j,:,:) = zero
            turnover_daily(i,j,:,:) = zero

            IF (ok_pc) THEN
               deepC_a(i,:,j)=zero
               deepC_s(i,:,j)=zero
               deepC_p(i,:,j)=zero
            ENDIF

          ENDIF
       ENDDO

       !! 2.4 update 10 year-turnover pool content following flux emission
       !!     (linear decay (10%) of the initial carbon input)
       DO  l = 0, 8
          m = 10 - l
          cflux_prod10(i,iwplcc) =  cflux_prod10(i,iwplcc) + flux10(i,m,iwplcc)
          prod10(i,m,iwplcc)     =  prod10(i,m-1,iwplcc)   - flux10(i,m-1,iwplcc)
          flux10(i,m,iwplcc)     =  flux10(i,m-1,iwplcc)
          
          IF (prod10(i,m,iwplcc) .LT. 1.0) prod10(i,m,iwplcc) = zero
       ENDDO
       
       cflux_prod10(i,iwplcc) = cflux_prod10(i,iwplcc) + flux10(i,1,iwplcc) 
       flux10(i,1,iwplcc)     = 0.1 * prod10(i,0,iwplcc)
       prod10(i,1,iwplcc)     = prod10(i,0,iwplcc)
       
       !! 2.5 update 100 year-turnover pool content following flux emission\n
       DO   l = 0, 98
          m = 100 - l
          cflux_prod100(i,iwplcc)  =  cflux_prod100(i,iwplcc) + flux100(i,m,iwplcc)
          prod100(i,m,iwplcc)      =  prod100(i,m-1,iwplcc)   - flux100(i,m-1,iwplcc)
          flux100(i,m,iwplcc)      =  flux100(i,m-1,iwplcc)
          
          IF (prod100(i,m,iwplcc).LT.1.0) prod100(i,m,iwplcc) = zero
       ENDDO
       
       cflux_prod100(i,iwplcc)  = cflux_prod100(i,iwplcc) + flux100(i,1,iwplcc) 
       flux100(i,1,iwplcc)      = 0.01 * prod100(i,0,iwplcc)
       prod100(i,1,iwplcc)      = prod100(i,0,iwplcc)
       prod10(i,0,iwplcc)        = zero
       prod100(i,0,iwplcc)       = zero 
       


    ENDDO ! Loop over # pixels - domain size
    
    !! We redistribute the updated litter into four fuel classes, so that
    !! the balance between aboveground litter and fuel is mainted. The subtraction
    !! of fuel burned by land cover change fires from the fuel pool is made here. 
    fuel_all_type(:,:,:,:) = fuel_1hr(:,:,:,:) + fuel_10hr(:,:,:,:) + &
                               fuel_100hr(:,:,:,:) + fuel_1000hr(:,:,:,:)
    fuel_type_frac(:,:,:,:,:) = 0.25
    WHERE(fuel_all_type(:,:,:,:) > min_stomate)
      fuel_type_frac(:,:,:,:,1) = fuel_1hr(:,:,:,:)/fuel_all_type(:,:,:,:)
      fuel_type_frac(:,:,:,:,2) = fuel_10hr(:,:,:,:)/fuel_all_type(:,:,:,:)
      fuel_type_frac(:,:,:,:,3) = fuel_100hr(:,:,:,:)/fuel_all_type(:,:,:,:)
      fuel_type_frac(:,:,:,:,4) = fuel_1000hr(:,:,:,:)/fuel_all_type(:,:,:,:)
    ENDWHERE
    DO j=1,nvm
      fuel_1hr(:,j,:,:) = litter(:,:,j,iabove,:) * fuel_type_frac(:,j,:,:,1) 
      fuel_10hr(:,j,:,:) = litter(:,:,j,iabove,:) * fuel_type_frac(:,j,:,:,2) 
      fuel_100hr(:,j,:,:) = litter(:,:,j,iabove,:) * fuel_type_frac(:,j,:,:,3) 
      fuel_1000hr(:,j,:,:) = litter(:,:,j,iabove,:) * fuel_type_frac(:,j,:,:,4) 
    END DO 

  !! 3. history
    convflux        = convflux/one_year*dt_days
    cflux_prod10    = cflux_prod10/one_year*dt_days
    cflux_prod100   = cflux_prod100/one_year*dt_days

    DO i=1, npts
       DO j=1,nvm
           DO l=1,ncarb
              soilc_after(i)=soilc_after(i)+carbon(i,l,j)*veget_cov_max_adjusted(i,j)
           ENDDO
       ENDDO
       IF ( ABS(soilc_after(i)-soilc_before(i)) .GT. 1000.*min_stomate) THEN
             WRITE(numout,*) ' qcj check lcchange_main_agripeat,C',lalo(i,:)
             WRITE (numout,*) 'qcj check lcchange_main_agripeat,before-after', soilc_before(i),soilc_after(i)
       ENDIF
    ENDDO


    IF (printlev>=4) WRITE(numout,*) 'Leaving lcchange_main_agripeat'



  END SUBROUTINE lcchange_main_agripeat

!!
!================================================================================================================================
!! SUBROUTINE   : agripeat_adjust_fractions
!! DESCRIPTION: Adjust fractions of vegetations according to peatland area
!_
!================================================================================================================================
  SUBROUTINE agripeat_adjust_fractions (npts, lalo,veget_cov_max_new, veget_cov_max_old, veget_cov_max_adjusted)

    IMPLICIT NONE
  !! 0. Variable and parameter declaration 

    !! Input variables
    INTEGER, INTENT(in)                                       :: npts !! Domain size - number of pixels (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: veget_cov_max_old
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: veget_cov_max_new
    REAL(r_std),DIMENSION(npts,2),INTENT(in)                   :: lalo

    !! Output variables
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)             :: veget_cov_max_adjusted  

   !!  Local variables
    REAL(r_std), DIMENSION(npts,nvm)                          :: veget_cov_max_tmp
    REAL(r_std), DIMENSION(npts,nvm)                          :: veget_cov_max_tmp2
!!!natural vegetations (not including natural peatland)
    REAL(r_std), DIMENSION(npts)                              :: sum_veget_natold
    REAL(r_std), DIMENSION(npts)                              :: sum_veget_natnew
    REAL(r_std), DIMENSION(npts)                              :: sum_veget_nattmp
    REAL(r_std), DIMENSION(npts)                              :: sum_veget_nattmp_adjust
!!! crops (peat+non-peat)
    REAL(r_std), DIMENSION(npts)                              :: sum_crops_old
    REAL(r_std), DIMENSION(npts)                              :: sum_crops_new
!!!crops on peatland
    REAL(r_std), DIMENSION(npts)                              :: sum_peat_crops_old
    REAL(r_std), DIMENSION(npts)                              :: sum_peat_crops_new
!!!natural peatland+ crops peatland
    REAL(r_std), DIMENSION(npts)                              :: sum_peat_old
    INTEGER(i_std)                                            :: i, j, k, l

    ! peatland PFT number in nvm
    INTEGER(i_std)                                            :: c3crop, c4crop, peat_grass, peat_c3crop, peat_c4crop, peat_forest
!_
!================================================================================================================================

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

    veget_cov_max_adjusted(:,:) = veget_cov_max_new(:,:)

    sum_veget_natold(:) = zero
    sum_veget_natnew(:) = zero
    sum_veget_nattmp(:) = zero
    
    sum_veget_nattmp_adjust(:) = zero
      
    sum_peat_crops_old(:)=zero
    sum_peat_crops_new(:)=zero
 
    sum_crops_old(:)=zero
    sum_crops_new(:)=zero

    sum_peat_old(:)=zero
    DO i = 1, npts
       DO j=1,nvm
          ! Natural Peatland + other Natural land 
          IF ( natural(j) .AND. .NOT. pasture(j) .AND. .NOT. is_peat(j) .AND. .NOT. is_croppeat(j) ) THEN
             sum_veget_natold(i) = sum_veget_natold(i) + veget_cov_max_old(i,j)
          ENDIF
          ! Crop Peatland
          IF (is_croppeat(j)) THEN
             sum_peat_crops_old(i) = sum_peat_crops_old(i) + veget_cov_max_old(i,j)
          ENDIF
       !ENDDO
       !DO j=1,nvm
          IF (is_peat(j) .OR. is_croppeat(j)) THEN
            sum_peat_old(i) = sum_peat_old(i) + veget_cov_max_old(i,j) !+sum_peat_crops_old(i)
          ENDIF
       ENDDO
    ENDDO

!!!!!!!!!!!!!!!!!!!! Hard-coded, may need to be improved. qcj++

  IF (agri_peat_prop) THEN  
!!!Overlay vegetation maps with peatland, assume that peatland occupy each non-peat PFT in proportion to grid cell area of the PFT
    DO i= 1, npts
       DO j=1,nvm !! PFT1-13, non-peatland, natural vegetations
          ! rescale all natural and crop non-peatland PFTs
          IF (.NOT. is_peat(j) .AND. .NOT. is_croppeat(j)) THEN
            veget_cov_max_tmp(i,j) = veget_cov_max_new(i,j) * (un - sum_peat_old(i))
          ENDIF
       ENDDO 
  !! peatland overlay with crops creating two new PFTs
       ! proportionally replaced JC comment: I do not get this
       veget_cov_max_tmp(i,peat_c3crop) = sum_peat_old(i)*veget_cov_max_new(i,c3crop)
       veget_cov_max_tmp(i,peat_c4crop) = sum_peat_old(i)*veget_cov_max_new(i,c4crop)
  !! substract PFT15 and PFT16 from PFT12 and PFT13, respectively
  !! subtract PFT15 and PFT16 from natural peat (PFT14)
       veget_cov_max_tmp(i,peat_grass) = sum_peat_old(i)-veget_cov_max_tmp(i,peat_c3crop)-veget_cov_max_tmp(i,peat_c4crop)
    ENDDO

    DO i = 1, npts
!! IF first occur of agriculture on peatland
      IF ( (sum_peat_old(i) .GT. zero) .AND. (sum_peat_crops_old(i) .EQ. zero)) THEN
         veget_cov_max_adjusted(i,:) = veget_cov_max_tmp(i,:)
      ENDIF

!! Change vegetation map 
      IF (sum_peat_crops_old(i) .GT. zero) THEN
   !! compare veget_cov_max_tmp with veget_cov_max_old
        DO j=1,nvm
           IF ( natural(j) .AND. .NOT. pasture(j)) THEN
              sum_veget_natnew(i) = sum_veget_natnew(i) + veget_cov_max_new(i,j)
           ENDIF
        ENDDO
        IF ((sum_veget_natnew(i) .LE. sum_veget_natold(i))) THEN
        ! sum of natural vegetations decrease, sum of PFT12 and PFT13 in the vegetation map increases
        ! natural peatland decrease. crops on peatland increase
            veget_cov_max_adjusted(i,:) = veget_cov_max_tmp(i,:)
        ELSE
        ! natural vegetations increase
        ! natural peatland should not increase 
           veget_cov_max_adjusted(i,peat_grass) = veget_cov_max_old(i,peat_grass)
          !! adjust non-peatland natural vegetations accordingly
           sum_veget_nattmp(i) = sum_veget_natnew(i)-veget_cov_max_old(i,peat_grass)
           DO j=1,nvm
              IF ( natural(j) .AND. .NOT. pasture(j)) THEN
                veget_cov_max_adjusted(i,j)=veget_cov_max_new(i,j)*sum_veget_nattmp(i)/sum_veget_natnew(i)
              ENDIF
           ENDDO  
          !! crops decrease                   
           veget_cov_max_adjusted(i,peat_c3crop) = veget_cov_max_tmp(i,peat_c3crop)
           veget_cov_max_adjusted(i,peat_c4crop) = veget_cov_max_tmp(i,peat_c4crop)
           veget_cov_max_adjusted(i,c3crop) = veget_cov_max_tmp(i,c3crop)
           veget_cov_max_adjusted(i,c4crop) = veget_cov_max_tmp(i,c4crop)
        ENDIF
      ENDIF
    ENDDO !(npts)
  ENDIF

  IF (agri_peat_MINcrop) THEN
!!! crops will not be planted on peatland, as long as there is enough non-peat mineral soil exists in the grid cell 
    DO i = 1, npts
        DO j=1,nvm
           IF ( natural(j) .AND. .NOT. pasture(j) .AND. .NOT. is_peat(j) .AND. .NOT. is_croppeat(j)) THEN
              sum_veget_natnew(i) = sum_veget_natnew(i) + veget_cov_max_new(i,j)
           ENDIF
        ENDDO
        sum_crops_new(i)=veget_cov_max_new(i,c3crop)+veget_cov_max_new(i,c4crop)

        IF (sum_veget_natnew(i) .GE. veget_cov_max_old(i,peat_grass)) THEN
        ! there is enough natural areas for peatland 
           veget_cov_max_adjusted(i,peat_grass) = veget_cov_max_old(i,peat_grass)
        ! substarct PFT14 from natural non-peat vegetations
           sum_veget_nattmp(i) = sum_veget_natnew(i)- veget_cov_max_adjusted(i,peat_grass)
           DO j=1,nvm
              IF ( natural(j) .AND. .NOT. pasture(j) .AND. &
                 .NOT. is_peat(j) .AND. .NOT. is_croppeat(j) .AND. &
                 (sum_veget_natnew(i) .GT. min_stomate) ) THEN
                 veget_cov_max_adjusted(i,j)= veget_cov_max_new(i,j)*sum_veget_nattmp(i)/sum_veget_natnew(i)
              ENDIF
           ENDDO

           sum_peat_crops_new(i) = MIN(sum_crops_new(i),sum_peat_old(i)-veget_cov_max_adjusted(i,peat_grass))  
           IF (sum_crops_new(i) .GT. min_stomate) THEN
              veget_cov_max_adjusted(i,peat_c3crop) = MAX (zero,sum_peat_crops_new(i)*veget_cov_max_new(i,c3crop)/sum_crops_new(i))
              veget_cov_max_adjusted(i,peat_c4crop) = MAX (zero,sum_peat_crops_new(i)*veget_cov_max_new(i,c4crop)/sum_crops_new(i))
           ENDIF
           veget_cov_max_adjusted(i,c3crop) = veget_cov_max_new(i,c3crop)-veget_cov_max_adjusted(i,peat_c3crop) 
           veget_cov_max_adjusted(i,c4crop) = veget_cov_max_new(i,c4crop)-veget_cov_max_adjusted(i,peat_c4crop)       
        ELSE 
        ! sum of all natural vegetaions < PFT14,thus all natural vegetations are PFT14 and part of peatland are occupied by crops  
          veget_cov_max_adjusted(i,peat_grass) = sum_veget_natnew(i)
        ! substarct PFT14 from natural non-peat vegetations
          DO j=1,nvm
             IF ( natural(j) .AND. .NOT. pasture(j)) THEN
                 veget_cov_max_adjusted(i,j)= zero
             ENDIF
          ENDDO

          sum_peat_crops_new(i) = sum_peat_old(i)-veget_cov_max_adjusted(i,peat_grass)
          IF (sum_crops_new(i) .GT. min_stomate) THEN
             veget_cov_max_adjusted(i,peat_c3crop) = sum_peat_crops_new(i)*veget_cov_max_new(i,c3crop)/sum_crops_new(i)
             veget_cov_max_adjusted(i,peat_c4crop) = sum_peat_crops_new(i)*veget_cov_max_new(i,c4crop)/sum_crops_new(i)
          ENDIF
          veget_cov_max_adjusted(i,c3crop) = veget_cov_max_new(i,c3crop)-veget_cov_max_adjusted(i,peat_c3crop)
          veget_cov_max_adjusted(i,c4crop) = veget_cov_max_new(i,c4crop)-veget_cov_max_adjusted(i,peat_c4crop)
        ENDIF
    ENDDO
  ENDIF

  IF (agri_peat_MAXcrop) THEN
!!! For a grid cell that has crops, crops will be planted on peatland first, then to be planted on mineral soils
    DO i = 1, npts
       DO j=1,nvm
          IF ( natural(j) .AND. .NOT. pasture(j) ) THEN
            sum_veget_natnew(i) = sum_veget_natnew(i) + veget_cov_max_new(i,j)
          ENDIF
       ENDDO
       sum_crops_new(i)=veget_cov_max_new(i,c3crop)+veget_cov_max_new(i,c4crop)

       IF (sum_crops_new(i) .GE. sum_peat_old(i)) THEN
       ! all peatland are disturbed
          IF (sum_crops_new(i) .GT. min_stomate) THEN
            veget_cov_max_adjusted(i,peat_c3crop)=sum_peat_old(i)*veget_cov_max_new(i,c3crop)/sum_crops_new(i)          
            veget_cov_max_adjusted(i,peat_c4crop)=sum_peat_old(i)*veget_cov_max_new(i,c4crop)/sum_crops_new(i)
          ENDIF
          veget_cov_max_adjusted(i,peat_grass)=zero
          veget_cov_max_adjusted(i,c3crop)=veget_cov_max_new(i,c3crop)-veget_cov_max_adjusted(i,peat_c3crop)
          veget_cov_max_adjusted(i,c4crop)=veget_cov_max_new(i,c4crop)-veget_cov_max_adjusted(i,peat_c4crop)
       ! non-peat PFTs, unchanged from the map
          DO j=1,nvm
            IF ( natural(j) .AND. .NOT. pasture(j)) THEN
               veget_cov_max_adjusted(i,j)= veget_cov_max_new(i,j)
            ENDIF
          ENDDO
       ELSE
       ! all crops are planted on peatland, and there is still some natural peatland
          veget_cov_max_adjusted(i,peat_c3crop)= veget_cov_max_new(i,c3crop)
          veget_cov_max_adjusted(i,peat_c4crop)= veget_cov_max_new(i,c4crop)
          veget_cov_max_adjusted(i,c3crop)= zero
          veget_cov_max_adjusted(i,c4crop)= zero
          ! natural peatland should not increase
          veget_cov_max_adjusted(i,peat_grass)= MIN(sum_peat_old(i)-veget_cov_max_adjusted(i,peat_c3crop)-veget_cov_max_adjusted(i,peat_c4crop), &
                                    veget_cov_max_old(i,peat_grass))
        ! adjust natural non-peat PFTs according to PFT14
          sum_veget_nattmp(i) = sum_veget_natnew(i)- veget_cov_max_adjusted(i,peat_grass)
          DO j=1,nvm
            IF ( natural(j) .AND. .NOT. pasture(j)) THEN
               veget_cov_max_adjusted(i,j)= veget_cov_max_new(i,j)*sum_veget_nattmp(i)/sum_veget_natnew(i)  
            ENDIF
          ENDDO
       ENDIF
    ENDDO
  ENDIF

  DO i = 1, npts
    IF (ABS(SUM(veget_cov_max_adjusted(i,:))-SUM(veget_cov_max_new(i,:))) .GT. min_stomate) THEN
       WRITE (numout,*) 'qcj check agripeat_adjust_fractions,veget', lalo(i,:)
       WRITE (numout,*) 'qcj check agripeat_adjust_fractions,after-before',SUM(veget_cov_max_adjusted(i,:)),SUM(veget_cov_max_new(i,:))
    ENDIF
  ENDDO


  END SUBROUTINE agripeat_adjust_fractions
!_
!================================================================================================================================

END MODULE stomate_lcchange
