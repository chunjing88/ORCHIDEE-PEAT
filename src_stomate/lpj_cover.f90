! =================================================================================================================================
! MODULE       : lpj_cover
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
!                This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        Recalculate vegetation cover and LAI
!!
!!\n DESCRIPTION : None
!!
!! RECENT CHANGE(S) : Including permafrost carbon
!!
!! REFERENCE(S) : 
!!        Sitch, S., B. Smith, et al. (2003), Evaluation of ecosystem dynamics,
!!        plant geography and terrestrial carbon cycling in the LPJ dynamic 
!!        global vegetation model, Global Change Biology, 9, 161-185.\n
!!        Smith, B., I. C. Prentice, et al. (2001), Representation of vegetation
!!        dynamics in the modelling of terrestrial ecosystems: comparing two
!!        contrasting approaches within European climate space,
!!        Global Ecology and Biogeography, 10, 621-637.\n
!!
!! SVN :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_stomate/lpj_cover.f90 $
!! $Date: 2017-10-26 15:35:04 +0200 (Thu, 26 Oct 2017) $
!! $Revision: 4719 $
!! \n
!_ ================================================================================================================================

MODULE lpj_cover

  ! modules used:
  USE ioipsl_para
  USE stomate_data
  USE pft_parameters
  USE constantes_soil_var
  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC cover

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE     : lpj_cover
!!
!>\BRIEF          Recalculate vegetation cover and LAI
!!
!!\n DESCRIPTION : Veget_cov_max is first renewed here according to newly calculated foliage biomass in this calculation step 
!! Then, litter, soil carbon, and biomass are also recalcuted with taking into account the changes in Veget_cov_max (i.e. delta_veg)
!! Grid-scale fpc (foliage projected coverage) is calculated to obtain the shadede ground area by leaf's light capture
!! Finally, grid-scale fpc is adjusted not to exceed 1.0
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : ::lai (leaf area index, @tex $(m^2 m^{-2})$ @endtex), 
!! :: veget (fractional vegetation cover, unitless)
!!
!! REFERENCE(S)   : None
!! 
!! FLOWCHART :
!! \latexonly 
!!     \includegraphics[scale=0.5]{lpj_cover_flowchart.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE cover (npts, cn_ind, ind, biomass, &
       veget_cov_max, veget_cov_max_old, lai, & 
       litter, litter_avail, litter_not_avail, carbon, &
       fuel_1hr, fuel_10hr, fuel_100hr, fuel_1000hr, &
       turnover_daily, bm_to_litter, &
       co2_to_bm, co2_fire, resp_hetero, resp_maint, resp_growth, gpp_daily, &
       deepC_a, deepC_s, deepC_p, &
       fpeat_map,veget_cov_max_new,lalo,wettile_dgvm,date) !!!qcj++ peatland

!! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                                  :: npts             !! Domain size (unitless)  
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)                :: cn_ind           !! Crown area 
                                                                                    !! @tex $(m^2)$ @endtex per individual
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)                :: ind              !! Number of individuals 
                                                                                    !! @tex $(m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)                :: veget_cov_max_old!! "Maximal" coverage fraction of a PFT (LAI-> 
                                                                                    !! infinity) on ground at beginning of time 
!!!qcj++ peatland
    REAL(r_std), DIMENSION(npts,nstm), INTENT(in)               :: fpeat_map
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)                :: veget_cov_max_new
    REAL(r_std),DIMENSION(npts,2),INTENT(in)                    :: lalo
    LOGICAL,DIMENSION (nstm), INTENT (in)                       :: wettile_dgvm
    INTEGER(i_std),INTENT(in)                                   :: date
    REAL(r_std), DIMENSION(npts)                                :: sum_files
    REAL(r_std), DIMENSION(npts,nstm)                           :: fpeat_map_adjust
    REAL(r_std), DIMENSION(npts)                                :: sum_fpeat_map
    REAL(r_std), DIMENSION(npts)                                :: sum_fpeat_map_adjust
    REAL(r_std), DIMENSION(npts,nvm)                            :: veget_peat
    REAL(r_std), DIMENSION(npts,nvm)                            :: veget_peat_old
    REAL(r_std), DIMENSION(npts,3)                              :: C_start
    REAL(r_std), DIMENSION(npts,3)                              :: C_end 
    REAL(r_std), DIMENSION(npts,3)                              :: C_check
    REAL(r_std),DIMENSION(npts,4)                               :: Fin_start
    REAL(r_std),DIMENSION(npts,4)                               :: Fin_end
    REAL(r_std),DIMENSION(npts,4)                               :: Fout_start
    REAL(r_std),DIMENSION(npts,4)                               :: Fout_end
    REAL(r_std),DIMENSION(npts,4)                               :: Fin_check
    REAL(r_std),DIMENSION(npts,4)                               :: Fout_check
    REAL(r_std), DIMENSION(npts)                                :: sum_veget
    REAL(r_std), DIMENSION(npts,nstm)                           :: sum_peat
    REAL(r_std), DIMENSION(npts)                                :: sum_peat_old
    REAL(r_std), DIMENSION(npts)                                :: sum_peat_new
    REAL(r_std), DIMENSION(npts,nstm)                           :: sum_peat_tiles
    REAL(r_std), DIMENSION(npts)                                :: sum_crop
    REAL(r_std), DIMENSION(npts)                                :: sum_peatBare
    REAL(r_std), DIMENSION(npts)                                :: sum_peatBare_old
    REAL(r_std), DIMENSION(npts)                                :: sum_peatBare_adjust
    REAL(r_std), DIMENSION(npts)                                :: delta_peatBare
    REAL(r_std)                                                 :: delta_veg_shrink
    REAL(r_std)                                                 :: delta_veg_expand

    !! 0.2 Output variables

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)             :: lai                 !! Leaf area index OF AN INDIVIDUAL PLANT 
                                                                                       !! @tex $(m^2 m^{-2})$ @endtex
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
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)                  :: veget_cov_max  !! "Maximal" coverage fraction of a PFT (LAI->
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

    !! 0.4 Local variables

    INTEGER(i_std)                                              :: i,j,k,m,jst           !! Index (unitless)
    REAL(r_std), DIMENSION(npts,nlitt,nlevs,nelements)          :: dilu_lit              !! Litter dilution @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,ncarb)                          :: dilu_soil_carbon      !! Soil Carbondilution 
                                                                                         !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nparts,nelements)               :: dilu_bio              !! Biomass dilution 
    REAL(r_std), DIMENSION(npts)                                :: dilu_TCarbon
    REAL(r_std), DIMENSION(npts,nparts,nelements)               :: dilu_turnover_daily
    REAL(r_std), DIMENSION(npts,nparts,nelements)               :: dilu_bm_to_litter
    REAL(r_std), DIMENSION(npts)                                :: dilu_co2flux_new
    REAL(r_std), DIMENSION(npts)                                :: dilu_gpp_daily
    REAL(r_std), DIMENSION(npts)                                :: dilu_resp_growth
    REAL(r_std), DIMENSION(npts)                                :: dilu_resp_maint
    REAL(r_std), DIMENSION(npts)                                :: dilu_resp_hetero
    REAL(r_std), DIMENSION(npts)                                :: dilu_co2_to_bm
    REAL(r_std), DIMENSION(npts)                                :: dilu_co2_fire
    REAL(r_std), DIMENSION(npts,nvm)                            :: TCarbon
    REAL(r_std), DIMENSION(npts,nvm)                            :: co2flux_new
    REAL(r_std), DIMENSION(npts,nvm)                            :: co2flux_old
    REAL(r_std), DIMENSION(npts,ncarb,nvm)                       :: carbon_old
    REAL(r_std),DIMENSION(npts,ndeep,ncarb)                     :: dilu_soil_carbon_vertres !!vertically-resolved Soil Carbondilution (gC/m²)

    REAL(r_std), DIMENSION(nvm)                                 :: delta_veg        !! Conversion factors (unitless)
    REAL(r_std), DIMENSION(nvm)                                 :: reduct           !! Conversion factors (unitless)
    REAL(r_std)                                                 :: delta_veg_sum    !! Conversion factors (unitless)
    REAL(r_std)                                                 :: diff             !! Conversion factors (unitless)
    REAL(r_std)                                                 :: sr               !! Conversion factors (unitless)
    REAL(r_std), DIMENSION(npts)                                :: frac_nat         !! Conversion factors (unitless)
    REAL(r_std), DIMENSION(npts)                                :: sum_vegettree    !! Conversion factors (unitless)
    REAL(r_std), DIMENSION(npts)                                :: sum_vegetgrass   !! Conversion factors (unitless) 
    REAL(r_std), DIMENSION(npts)                                :: sum_veget_natveg !! Conversion factors (unitless)
    REAL(r_std), DIMENSION(npts)                                :: vartmp           !! Temporary variable used to add history
    REAL(r_std), DIMENSION(npts,nlitt,nelements)                :: dilu_f1hr        !! Litter dilution @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nlitt,nelements)                :: dilu_f10hr       !! Litter dilution @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nlitt,nelements)                :: dilu_f100hr      !! Litter dilution @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nlitt,nelements)                :: dilu_f1000hr     !! Litter dilution @tex $(gC m^{-2})$ @endtex
!_ ================================================================================================================================

 !! 1. If the vegetation is dynamic, calculate new maximum vegetation cover for natural plants

  IF (.NOT. ok_dgvm_peat) THEN !!!qcj++ peatland
  
    IF ( ok_dgvm ) THEN

       !! 1.1  Calculate initial values of vegetation cover
       frac_nat(:) = un
       sum_veget_natveg(:) = zero
       veget_cov_max(:,ibare_sechiba) = un
       co2flux_new = undef
       co2flux_old = undef
       TCarbon = undef

       carbon_old(:,:,:)=carbon(:,:,:)

       DO j = 2,nvm ! loop over PFTs

          IF ( natural(j) .AND. .NOT. pasture(j)) THEN
	     
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
                IF( natural(j) .AND. .NOT. pasture(j)) THEN
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

       !! 1.2 Calculate carbon fluxes between PFTs to maintain mass balance
       !! Assure carbon closure when veget_cov_max changes(delta_veg): if veget_cov_max of some PFTs decrease, we use "dilu" to 
       !! record the corresponding lost in carbon (biomass, litter, soil carbon, gpp, respiration etc.) for 
       !! these PFTs, and re-allocate "dilu" to those PFTs with increasing veget_cov_max.
       DO i = 1, npts ! loop over grid points
          
          ! Calculate the change in veget_cov_max between previous time step and current time step
          delta_veg(:) = veget_cov_max(i,:)-veget_cov_max_old(i,:)
          delta_veg_sum = SUM(delta_veg,MASK=delta_veg.LT.zero)

          dilu_lit(i,:,:,:) = zero
          dilu_f1hr(i,:,:) = zero
          dilu_f10hr(i,:,:) = zero
          dilu_f100hr(i,:,:) = zero
          dilu_f1000hr(i,:,:) = zero
          dilu_soil_carbon(i,:) = zero
          dilu_soil_carbon_vertres(i,:,:) = zero

          dilu_bio(i,:,:) = zero
          dilu_TCarbon(i)=zero

          dilu_turnover_daily(i,:,:)=zero
          dilu_bm_to_litter(i,:,:)=zero
          dilu_co2flux_new(i)=zero
          dilu_gpp_daily(i)=zero
          dilu_resp_growth(i)=zero
          dilu_resp_maint(i)=zero
          dilu_resp_hetero(i)=zero
          dilu_co2_to_bm(i)=zero
          dilu_co2_fire(i)=zero

          ! Calculate TCarbon: total carbon including biomass, litter and soil carbon, as well as "today's" turnover and 
          ! bm_to_litter due to mortality, because today's turnover and bm_to_litter are not yet added into "litter" until tomorrow. 
          DO j=1, nvm
                TCarbon(i,j)=SUM(biomass(i,j,:,icarbon))+SUM(carbon(i,:,j))+SUM(litter(i,:,j,:,icarbon))+SUM(turnover_daily(i,j,:,icarbon))+SUM(bm_to_litter(i,j,:,icarbon))
                co2flux_old(i,j)=resp_maint(i,j)+resp_growth(i,j)+resp_hetero(i,j)+co2_fire(i,j)-co2_to_bm(i,j)-gpp_daily(i,j)
                co2flux_new(i,j)=resp_maint(i,j)+resp_growth(i,j)+resp_hetero(i,j)+co2_fire(i,j)-co2_to_bm(i,j)-gpp_daily(i,j)
          ENDDO

          DO j=1, nvm ! loop over PFTs
             IF ( delta_veg(j) < -min_stomate ) THEN
                dilu_lit(i,:,:,:) =  dilu_lit(i,:,:,:) + delta_veg(j) * litter(i,:,j,:,:) / delta_veg_sum
                dilu_f1hr(i,:,:) =  dilu_f1hr(i,:,:) + delta_veg(j) * fuel_1hr(i,j,:,:) / delta_veg_sum
                dilu_f10hr(i,:,:) =  dilu_f10hr(i,:,:) + delta_veg(j) * fuel_10hr(i,j,:,:) / delta_veg_sum
                dilu_f100hr(i,:,:) =  dilu_f100hr(i,:,:) + delta_veg(j) * fuel_100hr(i,j,:,:) / delta_veg_sum
                dilu_f1000hr(i,:,:) =  dilu_f1000hr(i,:,:) + delta_veg(j) * fuel_1000hr(i,j,:,:) / delta_veg_sum
                dilu_soil_carbon(i,:) =  dilu_soil_carbon(i,:) + delta_veg(j) * carbon(i,:,j) / delta_veg_sum
                dilu_TCarbon(i)= dilu_TCarbon(i) + delta_veg(j) * TCarbon(i,j) / delta_veg_sum
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

          DO j=1, nvm ! loop over PFTs
             IF ( delta_veg(j) > min_stomate) THEN

                ! Dilution of reservoirs
                ! Recalculate the litter and soil carbon with taking into accout the change in 
                ! veget_cov_max (delta_veg)
                ! Litter
                litter(i,:,j,:,:)=(litter(i,:,j,:,:) * veget_cov_max_old(i,j) + dilu_lit(i,:,:,:) * delta_veg(j)) &
                                  / veget_cov_max(i,j)
                fuel_1hr(i,j,:,:)=(fuel_1hr(i,j,:,:) * veget_cov_max_old(i,j) + dilu_f1hr(i,:,:) * delta_veg(j)) / veget_cov_max(i,j)
                fuel_10hr(i,j,:,:)=(fuel_10hr(i,j,:,:) * veget_cov_max_old(i,j) + dilu_f10hr(i,:,:) * delta_veg(j)) / veget_cov_max(i,j)
                fuel_100hr(i,j,:,:)=(fuel_100hr(i,j,:,:) * veget_cov_max_old(i,j) + dilu_f100hr(i,:,:) * delta_veg(j)) / veget_cov_max(i,j)
                fuel_1000hr(i,j,:,:)=(fuel_1000hr(i,j,:,:) * veget_cov_max_old(i,j) + dilu_f1000hr(i,:,:) * delta_veg(j)) / veget_cov_max(i,j)
                !JCADD available and not available litter for grazing
                ! only not available litter change, available litter will not
                ! change, because tree litter can not be eaten
               IF (is_grassland_manag(j) .AND. is_grassland_grazed(j)) THEN
                 litter_avail(i,:,j) = litter_avail(i,:,j) * veget_cov_max_old(i,j) / veget_cov_max(i,j)
                 litter_not_avail(i,:,j) = litter(i,:,j,iabove,icarbon) - litter_avail(i,:,j)
               ENDIF
                !ENDJCADD    
                !IF ( ok_pc ) THEN
                !   deepC_a(i,:,j)=(deepC_a(i,:,j) * veget_cov_max_old(i,j) + &
                !        dilu_soil_carbon_vertres(i,:,iactive) * delta_veg(j)) / veget_cov_max(i,j)
                !   deepC_s(i,:,j)=(deepC_s(i,:,j) * veget_cov_max_old(i,j) + &
                !        dilu_soil_carbon_vertres(i,:,islow) * delta_veg(j)) / veget_cov_max(i,j)
                !   deepC_p(i,:,j)=(deepC_p(i,:,j) * veget_cov_max_old(i,j) + &
                !        dilu_soil_carbon_vertres(i,:,ipassive) * delta_veg(j)) / veget_cov_max(i,j)
                !ENDIF
                ! Soil carbon
                carbon(i,:,j)=(carbon(i,:,j) * veget_cov_max_old(i,j) + dilu_soil_carbon(i,:) * delta_veg(j)) / veget_cov_max(i,j)
                IF ( ok_pc ) THEN
                   IF (carbon_old(i,iactive,j) .GT. min_stomate) THEN
                      deepC_a(i,:,j)=deepC_a(i,:,j)*carbon(i,iactive,j)/carbon_old(i,iactive,j)
                   ENDIF
                   IF (carbon_old(i,islow,j) .GT. min_stomate) THEN
                      deepC_s(i,:,j)=deepC_s(i,:,j)*carbon(i,islow,j)/carbon_old(i,islow,j)
                   ENDIF
                   IF (carbon_old(i,ipassive,j) .GT. min_stomate) THEN
                      deepC_p(i,:,j)=deepC_p(i,:,j)*carbon(i,ipassive,j)/carbon_old(i,ipassive,j)
                   ENDIF
                ENDIF

                !biomass(i,j,:,:)=(biomass(i,j,:,:) * veget_cov_max_old(i,j) + dilu_bio(i,:,:) * delta_veg(j)) / veget_cov_max(i,j)
                TCarbon(i,j)=(TCarbon(i,j) * veget_cov_max_old(i,j) + dilu_TCarbon(i) * delta_veg(j)) / veget_cov_max(i,j)

                turnover_daily(i,j,:,:)=(turnover_daily(i,j,:,:)*veget_cov_max_old(i,j)+dilu_turnover_daily(i,:,:)*delta_veg(j))/veget_cov_max(i,j)
                bm_to_litter(i,j,:,:)=(bm_to_litter(i,j,:,:)*veget_cov_max_old(i,j)+dilu_bm_to_litter(i,:,:)*delta_veg(j))/veget_cov_max(i,j)
                co2flux_new(i,j)=(co2flux_old(i,j)*veget_cov_max_old(i,j)+dilu_co2flux_new(i)*delta_veg(j))/veget_cov_max(i,j)
                gpp_daily(i,j)=(gpp_daily(i,j)*veget_cov_max_old(i,j)+dilu_gpp_daily(i)*delta_veg(j))/veget_cov_max(i,j)
                resp_growth(i,j)=(resp_growth(i,j)*veget_cov_max_old(i,j)+dilu_resp_growth(i)*delta_veg(j))/veget_cov_max(i,j)
                resp_maint(i,j)=(resp_maint(i,j)*veget_cov_max_old(i,j)+dilu_resp_maint(i)*delta_veg(j))/veget_cov_max(i,j)
                resp_hetero(i,j)=(resp_hetero(i,j)*veget_cov_max_old(i,j)+dilu_resp_hetero(i)*delta_veg(j))/veget_cov_max(i,j)
                co2_to_bm(i,j)=(co2_to_bm(i,j)*veget_cov_max_old(i,j)+dilu_co2_to_bm(i)*delta_veg(j))/veget_cov_max(i,j)
                co2_fire(i,j)=(co2_fire(i,j)*veget_cov_max_old(i,j)+dilu_co2_fire(i)*delta_veg(j))/veget_cov_max(i,j)

             ENDIF

             IF(veget_cov_max(i,j).GT.min_stomate) THEN

                ! Correct biomass densities to conserve mass
                ! since it's defined on veget_cov_max
                biomass(i,j,:,:) = biomass(i,j,:,:) * veget_cov_max_old(i,j) / veget_cov_max(i,j)

             ENDIF

          ENDDO ! loop over PFTs
      ENDDO ! loop over grid points

      vartmp(:)=SUM(co2flux_new*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tCO2FLUX", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(co2flux_old*veget_cov_max_old,dim=2)
      CALL histwrite_p (hist_id_stomate, "tCO2FLUX_OLD", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(TCarbon*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tCARBON", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(gpp_daily*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tGPP", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(resp_growth*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tRESP_GROWTH", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(resp_maint*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tRESP_MAINT", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(resp_hetero*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tRESP_HETERO", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(co2_to_bm*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tCO2_TAKEN", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(co2_fire*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tCO2_FIRE", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(SUM(biomass(:,:,:,icarbon),dim=3)*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tBIOMASS", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(SUM(SUM(litter(:,:,:,:,icarbon),dim=4),dim=2)*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tLITTER", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(SUM(fuel_1hr(:,:,:,icarbon),dim=3)*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tFUEL1HR", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(SUM(fuel_10hr(:,:,:,icarbon),dim=3)*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tFUEL10HR", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(SUM(fuel_100hr(:,:,:,icarbon),dim=3)*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tFUEL100HR", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(SUM(fuel_1000hr(:,:,:,icarbon),dim=3)*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tFUEL1000HR", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(SUM(carbon,dim=2)*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tSOILC", itime, vartmp, npts, hori_index)

      IF ( ok_pc ) THEN
        vartmp(:)=SUM(SUM(deepC_a,dim=2)*veget_cov_max,dim=2)
        CALL histwrite_p (hist_id_stomate, "tDEEPCa", itime, vartmp, npts, hori_index)
        vartmp(:)=SUM(SUM(deepC_s,dim=2)*veget_cov_max,dim=2)
        CALL histwrite_p (hist_id_stomate, "tDEEPCs", itime, vartmp, npts, hori_index)
        vartmp(:)=SUM(SUM(deepC_p,dim=2)*veget_cov_max,dim=2)
        CALL histwrite_p (hist_id_stomate, "tDEEPCp", itime, vartmp, npts, hori_index)
      ENDIF

    ENDIF

  ELSE !IF ok_dgvm_peat

       C_start(:,:)=zero
       C_end(:,:)=zero

       Fin_start(:,:)=zero
       Fin_end(:,:)=zero
       Fout_start(:,:)=zero
       Fout_end(:,:)=zero

       C_start(:,1)=SUM(SUM(biomass(:,:,:,icarbon),DIM=3) * veget_cov_max_old,DIM=2)
       C_start(:,2)=SUM(SUM(SUM(litter(:,:,:,:,icarbon),DIM=4),DIM=2) *veget_cov_max_old,DIM=2)
       C_start(:,3)=SUM(SUM(carbon(:,:,:),DIM=2) * veget_cov_max_old,DIM=2)

       Fin_start(:,1)=SUM(SUM(bm_to_litter(:,:,:,icarbon),DIM=3) * veget_cov_max_old,DIM=2)
       Fin_start(:,2)=SUM(SUM(turnover_daily(:,:,:,icarbon),DIM=3) * veget_cov_max_old,DIM=2)
       Fin_start(:,3)=SUM(co2_to_bm * veget_cov_max_old,DIM=2)
       Fin_start(:,4)=SUM(gpp_daily* veget_cov_max_old,DIM=2)

       Fout_start(:,1)=SUM(co2_fire* veget_cov_max_old,DIM=2)
       Fout_start(:,2)=SUM(resp_hetero* veget_cov_max_old,DIM=2)
       Fout_start(:,3)=SUM(resp_maint* veget_cov_max_old,DIM=2)
       Fout_start(:,4)=SUM(resp_growth* veget_cov_max_old,DIM=2)

       !! 1.1  Calculate initial values of vegetation cover
       !sum_veget_natveg(:) = zero
       co2flux_new(:,:) = zero
       co2flux_old(:,:) = zero
       TCarbon(:,:) = zero
       carbon_old(:,:,:)=carbon(:,:,:)

       sum_peat(:,:)= zero
       sum_crop(:)= zero
       sum_veget(:)=zero
       sum_peatBare(:)=zero
       sum_peatBare_adjust(:)=zero
       sum_peat_old(:)=zero
       sum_peat_new(:)=zero 
       sum_peat_tiles(:,:)=zero

       veget_cov_max(:,:)=veget_cov_max_old(:,:)

       sum_fpeat_map(:)=zero
       fpeat_map_adjust(:,:)=zero

       DO jst=1,nstm  
          IF (wettile_dgvm(jst)) THEN
           sum_fpeat_map(:)=sum_fpeat_map(:)+fpeat_map(:,jst)
          ENDIF
       ENDDO
    
       DO j=2,nvm
         IF ( (.NOT. natural(j)) .OR. pasture(j) ) THEN
           sum_crop(:)=sum_crop(:)+veget_cov_max_old(:,j)
         ENDIF
       ENDDO

!!!natural wetlands, baresoil, croplands are read from input files
!!!make sure the two files compatible with each other 
!!!(natural wetlands+croplands+baresoil <=1.0)
       DO j = 2,nvm
         IF ( natural(j) .AND. .NOT. pasture(j) .AND. is_peat(j) ) THEN
           sum_peat_old(:)=sum_peat_old(:)+veget_cov_max_old(:,j)
         ENDIF
       ENDDO

       IF ( date == 1 ) THEN
          sum_files(:)=sum_fpeat_map(:)+sum_crop(:)+veget_cov_max_old(:,ibare_sechiba)      
       ELSE
          sum_files(:)=sum_peat_old(:)+sum_crop(:)+veget_cov_max_old(:,ibare_sechiba)
       ENDIF      

       DO i =1,npts

          !WRITE (numout,*) 'QCJ check veget_cov_max_old',veget_cov_max_old(i,:),'sum_fpeat_map',sum_fpeat_map(i)  
          IF (date==1) THEN  
            IF (sum_files(i) .GT. un) THEN
              sum_fpeat_map_adjust(i)=un-sum_crop(i)-veget_cov_max_old(i,ibare_sechiba)
             !WRITE (numout,*) 'QCJ check1,sum_fpeat_map_adjust',sum_fpeat_map_adjust(i),'sum_crop',sum_crop(i),'veget_cov_max_old bare',veget_cov_max_old(i,ibare_sechiba)
            ELSE
              sum_fpeat_map_adjust(i)=sum_fpeat_map(i)
             !WRITE (numout,*) 'QCJ check2,sum_fpeat_map_adjust',sum_fpeat_map_adjust(i)
            ENDIF
          ELSE
              sum_fpeat_map_adjust(i)=sum_peat_old(i)+veget_cov_max_old(i,ibare_sechiba)    
              !WRITE (numout,*) 'QCJ check3,sum_fpeat_map_adjust',sum_fpeat_map_adjust(i)
          ENDIF
            
          DO jst=1,nstm
            IF ( wettile_dgvm(jst) .AND. (sum_fpeat_map(i) .GT. zero) ) THEN  
              fpeat_map_adjust(i,jst)=fpeat_map(i,jst)/sum_fpeat_map(i)*sum_fpeat_map_adjust(i)
            ENDIF
          ENDDO 
       ENDDO
      
       IF ( date == 1 ) THEN
         sum_peatBare(:)=sum_fpeat_map_adjust(:)+veget_cov_max_old(:,ibare_sechiba)
       ELSE 
         sum_peatBare(:)=sum_peat_old(:)+veget_cov_max_old(:,ibare_sechiba)
       ENDIF  
  
       DO jst=1,nstm
         IF ( wettile_dgvm(jst) ) THEN
           DO j = 2,nvm ! loop over PFTs
             IF ( natural(j) .AND. .NOT. pasture(j) .AND. is_peat(j) .AND. pref_soil_veg(j)==jst ) THEN
             ! Summation of individual tree crown area to get total foliar projected coverage
               veget_cov_max(:,j) = ind(:,j) * cn_ind(:,j)
               sum_peat_tiles(:,jst)=sum_peat_tiles(:,jst)+veget_cov_max(:,j)
             ENDIF
           ENDDO
         ENDIF
       ENDDO      
         
       DO i = 1, npts
         DO jst=1,nstm
           IF ( wettile_dgvm(jst) .AND. (sum_peat_tiles(i,jst) .GT. fpeat_map_adjust(i,jst)) ) THEN
             !WRITE (numout,*) 'QCJ check sum_peat_tiles,',sum_peat_tiles(i,jst),'fpeat_map_adjust',fpeat_map_adjust(i,jst)
             DO j = 2,nvm ! loop over PFTs
               IF ( natural(j) .AND. .NOT. pasture(j) .AND. is_peat(j) .AND. pref_soil_veg(j)==jst ) THEN
                 veget_cov_max(i,j)=veget_cov_max(i,j)/sum_peat_tiles(i,jst)*fpeat_map_adjust(i,jst)
               ENDIF 
             ENDDO
           ENDIF    
         ENDDO     
       ENDDO             

       DO j = 2,nvm
         IF ( natural(j) .AND. .NOT. pasture(j) .AND. is_peat(j) ) THEN
           sum_peat_new(:)=sum_peat_new(:)+veget_cov_max(:,j)
         ENDIF
       ENDDO
         
       DO i = 1, npts
         veget_cov_max(i,ibare_sechiba)=MAX( zero, sum_peatBare(i)-sum_peat_new(i))
       ENDDO  

       sum_veget(:)=sum_peat_new(:)+sum_crop(:)+veget_cov_max(:,ibare_sechiba)

       DO i=1,npts
         IF (sum_veget(i) .GT. un+min_stomate) THEN
           WRITE (numout,*) 'QCJ check sum veget greater than 1.0, in lpj_cover, at:',lalo(i,:)
           WRITE (numout,*) 'QCJ check sum_veget',sum_veget(i),'sum_peat_new,sum_crop,veget_cov_max',sum_peat_new(i),sum_crop(i),veget_cov_max(i,ibare_sechiba)
           WRITE (numout,*) 'QCJ check veget in lpg_cover,veget_cov_max',veget_cov_max(i,:)
           WRITE (numout,*) 'QCJ check veget in lpg_cover,veget_cov_max_old',veget_cov_max_old(i,:)
           CALL ipslerr_p(3,'lpj_cover','ok_dgvm_peat=T','sum veget_cov_max>1.0','')
         ENDIF
       ENDDO

       !! 1.2 Calculate carbon fluxes between PFTs to maintain mass balance
       !! Assure carbon closure when veget_cov_max changes(delta_veg): if veget_cov_max of some PFTs decrease, we use "dilu" to
       !! record the corresponding lost in carbon (biomass, litter, soil carbon, gpp, respiration etc.) for
       !! these PFTs, and re-allocate "dilu" to those PFTs with increasing veget_cov_max.

     DO jst=1,nstm
      IF (wettile_dgvm(jst)) THEN

        DO i = 1, npts ! loop over grid points

          ! Calculate the change in veget_cov_max between previous time step and current time step
          delta_veg(:)=veget_cov_max(i,:)-veget_cov_max_old(i,:)

          delta_veg_shrink=zero
          delta_veg_expand=zero
  
          DO j = 1,nvm ! loop over PFTs
             IF ( natural(j) .AND. .NOT. pasture(j) .AND. is_peat(j) .AND. pref_soil_veg(j)==jst ) THEN
                IF (delta_veg(j) .GT. min_stomate) THEN  
                  delta_veg_expand = delta_veg_expand + delta_veg(j)
                ELSEIF (delta_veg(j) .LT. -min_stomate) THEN
                  delta_veg_shrink = delta_veg_shrink + delta_veg(j)
                ENDIF
             ENDIF
          ENDDO
 
          dilu_lit(i,:,:,:) = zero
          dilu_f1hr(i,:,:) = zero
          dilu_f10hr(i,:,:) = zero
          dilu_f100hr(i,:,:) = zero
          dilu_f1000hr(i,:,:) = zero
          dilu_soil_carbon(i,:) = zero
          dilu_soil_carbon_vertres(i,:,:) = zero

          dilu_bio(i,:,:) = zero
          dilu_TCarbon(i)=zero

          dilu_turnover_daily(i,:,:)=zero
          dilu_bm_to_litter(i,:,:)=zero
          dilu_co2flux_new(i)=zero
          dilu_gpp_daily(i)=zero
          dilu_resp_growth(i)=zero
          dilu_resp_maint(i)=zero
          dilu_resp_hetero(i)=zero
          dilu_co2_to_bm(i)=zero
          dilu_co2_fire(i)=zero

          ! Calculate TCarbon: total carbon including biomass, litter and soil carbon, as well as "today's" turnover and
          ! bm_to_litter due to mortality, because today's turnover and bm_to_litter are not yet added into "litter" until tomorrow.
          DO j=1, nvm
             IF ( natural(j) .AND. .NOT. pasture(j) .AND. is_peat(j) .AND. pref_soil_veg(j)==jst ) THEN
                TCarbon(i,j)=SUM(biomass(i,j,:,icarbon))+SUM(carbon(i,:,j))+SUM(litter(i,:,j,:,icarbon))+SUM(turnover_daily(i,j,:,icarbon))+SUM(bm_to_litter(i,j,:,icarbon))
                co2flux_old(i,j)=resp_maint(i,j)+resp_growth(i,j)+resp_hetero(i,j)+co2_fire(i,j)-co2_to_bm(i,j)-gpp_daily(i,j)
                co2flux_new(i,j)=resp_maint(i,j)+resp_growth(i,j)+resp_hetero(i,j)+co2_fire(i,j)-co2_to_bm(i,j)-gpp_daily(i,j)
             ENDIF
          ENDDO

          IF ( delta_veg_shrink .LT. -min_stomate ) THEN

            DO j=1, nvm ! loop over PFTs
             IF ( (delta_veg(j) < -min_stomate) .AND. is_peat(j) .AND. pref_soil_veg(j)==jst .AND. natural(j) ) THEN !-min_stomate
                dilu_lit(i,:,:,:) =  dilu_lit(i,:,:,:) - delta_veg(j) * litter(i,:,j,:,:) 
                dilu_f1hr(i,:,:) =  dilu_f1hr(i,:,:) - delta_veg(j) * fuel_1hr(i,j,:,:) 
                dilu_f10hr(i,:,:) =  dilu_f10hr(i,:,:) - delta_veg(j) * fuel_10hr(i,j,:,:)  
                dilu_f100hr(i,:,:) =  dilu_f100hr(i,:,:) - delta_veg(j) * fuel_100hr(i,j,:,:) 
                dilu_f1000hr(i,:,:) =  dilu_f1000hr(i,:,:) - delta_veg(j) * fuel_1000hr(i,j,:,:) 
                dilu_soil_carbon(i,:) =  dilu_soil_carbon(i,:) - delta_veg(j) * carbon(i,:,j) 
                dilu_TCarbon(i)= dilu_TCarbon(i) - delta_veg(j) * TCarbon(i,j) 
                dilu_turnover_daily(i,:,:)=dilu_turnover_daily(i,:,:) - delta_veg(j)*turnover_daily(i,j,:,:)
                dilu_bm_to_litter(i,:,:)=dilu_bm_to_litter(i,:,:) - delta_veg(j)*bm_to_litter(i,j,:,:)
                dilu_co2flux_new(i)=dilu_co2flux_new(i) - delta_veg(j)*co2flux_old(i,j)
                dilu_gpp_daily(i)=dilu_gpp_daily(i) - delta_veg(j)*gpp_daily(i,j)
                dilu_resp_growth(i)=dilu_resp_growth(i) - delta_veg(j)*resp_growth(i,j)
                dilu_resp_maint(i)=dilu_resp_maint(i) - delta_veg(j)*resp_maint(i,j)
                dilu_resp_hetero(i)=dilu_resp_hetero(i) - delta_veg(j)*resp_hetero(i,j)
                dilu_co2_to_bm(i)=dilu_co2_to_bm(i) - delta_veg(j)*co2_to_bm(i,j)
                dilu_co2_fire(i)=dilu_co2_fire(i) - delta_veg(j)*co2_fire(i,j)

                IF ( ok_pc ) THEN
                   dilu_soil_carbon_vertres(i,:,iactive)=dilu_soil_carbon_vertres(i,:,iactive) - delta_veg(j)*deepC_a(i,:,j)
                   dilu_soil_carbon_vertres(i,:,islow)=dilu_soil_carbon_vertres(i,:,islow) - delta_veg(j)*deepC_s(i,:,j)
                   dilu_soil_carbon_vertres(i,:,ipassive)=dilu_soil_carbon_vertres(i,:,ipassive) - delta_veg(j)*deepC_p(i,:,j)
                ENDIF
             ENDIF
            ENDDO ! loop over PFTs
          ENDIF ! delta_veg_shrink<-min_stomate
      
          IF ( delta_veg_expand .GT. min_stomate) THEN  

           IF (delta_veg_shrink .LT. -min_stomate) THEN !Expanded PFT get C from shrinked PFT   

             IF ( delta_veg(ibare_sechiba) .LT. -min_stomate ) THEN
                delta_veg_shrink=delta_veg_shrink+delta_veg(ibare_sechiba)
                dilu_lit(i,:,:,:) =  dilu_lit(i,:,:,:) - delta_veg(ibare_sechiba) * litter(i,:,ibare_sechiba,:,:)
                dilu_f1hr(i,:,:) =  dilu_f1hr(i,:,:) - delta_veg(ibare_sechiba) * fuel_1hr(i,ibare_sechiba,:,:)
                dilu_f10hr(i,:,:) =  dilu_f10hr(i,:,:) - delta_veg(ibare_sechiba) * fuel_10hr(i,ibare_sechiba,:,:)
                dilu_f100hr(i,:,:) =  dilu_f100hr(i,:,:) - delta_veg(ibare_sechiba) * fuel_100hr(i,ibare_sechiba,:,:)
                dilu_f1000hr(i,:,:) =  dilu_f1000hr(i,:,:) - delta_veg(ibare_sechiba) * fuel_1000hr(i,ibare_sechiba,:,:)
                dilu_soil_carbon(i,:) =  dilu_soil_carbon(i,:) - delta_veg(ibare_sechiba) * carbon(i,:,ibare_sechiba)
                dilu_TCarbon(i)= dilu_TCarbon(i) - delta_veg(ibare_sechiba) * TCarbon(i,ibare_sechiba)
                dilu_turnover_daily(i,:,:)=dilu_turnover_daily(i,:,:) - delta_veg(ibare_sechiba)*turnover_daily(i,ibare_sechiba,:,:)
                dilu_bm_to_litter(i,:,:)=dilu_bm_to_litter(i,:,:) - delta_veg(ibare_sechiba)*bm_to_litter(i,ibare_sechiba,:,:)
                dilu_co2flux_new(i)=dilu_co2flux_new(i) - delta_veg(ibare_sechiba)*co2flux_old(i,ibare_sechiba)
                dilu_gpp_daily(i)=dilu_gpp_daily(i) - delta_veg(ibare_sechiba)*gpp_daily(i,ibare_sechiba)
                dilu_resp_growth(i)=dilu_resp_growth(i) - delta_veg(ibare_sechiba)*resp_growth(i,ibare_sechiba)
                dilu_resp_maint(i)=dilu_resp_maint(i) - delta_veg(ibare_sechiba)*resp_maint(i,ibare_sechiba)
                dilu_resp_hetero(i)=dilu_resp_hetero(i) - delta_veg(ibare_sechiba)*resp_hetero(i,ibare_sechiba)
                dilu_co2_to_bm(i)=dilu_co2_to_bm(i) - delta_veg(ibare_sechiba)*co2_to_bm(i,ibare_sechiba)
                dilu_co2_fire(i)=dilu_co2_fire(i) - delta_veg(ibare_sechiba)*co2_fire(i,ibare_sechiba)

                IF ( ok_pc ) THEN
                   dilu_soil_carbon_vertres(i,:,iactive)=dilu_soil_carbon_vertres(i,:,iactive) - delta_veg(ibare_sechiba)*deepC_a(i,:,ibare_sechiba)            
                   dilu_soil_carbon_vertres(i,:,islow)=dilu_soil_carbon_vertres(i,:,islow) - delta_veg(ibare_sechiba)*deepC_s(i,:,ibare_sechiba)
                   dilu_soil_carbon_vertres(i,:,ipassive)=dilu_soil_carbon_vertres(i,:,ipassive) - delta_veg(ibare_sechiba)*deepC_p(i,:,ibare_sechiba)
                ENDIF
             ELSEIF ( delta_veg(ibare_sechiba) .GT. min_stomate ) THEN
                litter(i,:,ibare_sechiba,:,:)=litter(i,:,ibare_sechiba,:,:) * veget_cov_max_old(i,ibare_sechiba)/veget_cov_max(i,ibare_sechiba)
                fuel_1hr(i,ibare_sechiba,:,:)=fuel_1hr(i,ibare_sechiba,:,:) * veget_cov_max_old(i,ibare_sechiba)/veget_cov_max(i,ibare_sechiba)
                fuel_10hr(i,ibare_sechiba,:,:)=fuel_10hr(i,ibare_sechiba,:,:) * veget_cov_max_old(i,ibare_sechiba)/veget_cov_max(i,ibare_sechiba)
                fuel_100hr(i,ibare_sechiba,:,:)=fuel_100hr(i,ibare_sechiba,:,:) * veget_cov_max_old(i,ibare_sechiba)/veget_cov_max(i,ibare_sechiba)
                fuel_1000hr(i,ibare_sechiba,:,:)=fuel_1000hr(i,ibare_sechiba,:,:) * veget_cov_max_old(i,ibare_sechiba)/veget_cov_max(i,ibare_sechiba)
                carbon(i,:,ibare_sechiba)=carbon(i,:,ibare_sechiba) * veget_cov_max_old(i,ibare_sechiba)/veget_cov_max(i,ibare_sechiba)
                IF (ok_pc) THEN
                   deepC_a(i,:,ibare_sechiba)=deepC_a(i,:,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba)/veget_cov_max(i,ibare_sechiba)
                   deepC_s(i,:,ibare_sechiba)=deepC_s(i,:,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba)/veget_cov_max(i,ibare_sechiba)
                   deepC_p(i,:,ibare_sechiba)=deepC_p(i,:,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba)/veget_cov_max(i,ibare_sechiba)
                ENDIF
                TCarbon(i,ibare_sechiba)=TCarbon(i,ibare_sechiba) * veget_cov_max_old(i,ibare_sechiba)/veget_cov_max(i,ibare_sechiba)
                turnover_daily(i,ibare_sechiba,:,:)=turnover_daily(i,ibare_sechiba,:,:)*veget_cov_max_old(i,ibare_sechiba)/veget_cov_max(i,ibare_sechiba)
                bm_to_litter(i,ibare_sechiba,:,:)=bm_to_litter(i,ibare_sechiba,:,:)*veget_cov_max_old(i,ibare_sechiba)/veget_cov_max(i,ibare_sechiba)
                co2flux_new(i,ibare_sechiba)=co2flux_old(i,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba)/veget_cov_max(i,ibare_sechiba)
                gpp_daily(i,ibare_sechiba)=gpp_daily(i,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba)/veget_cov_max(i,ibare_sechiba)
                resp_growth(i,ibare_sechiba)=resp_growth(i,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba)/veget_cov_max(i,ibare_sechiba)
                resp_maint(i,ibare_sechiba)=resp_maint(i,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba)/veget_cov_max(i,ibare_sechiba)
                resp_hetero(i,ibare_sechiba)=resp_hetero(i,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba)/veget_cov_max(i,ibare_sechiba)
                co2_to_bm(i,ibare_sechiba)=co2_to_bm(i,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba)/veget_cov_max(i,ibare_sechiba)
                co2_fire(i,ibare_sechiba)=co2_fire(i,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba)/veget_cov_max(i,ibare_sechiba)              
             ENDIF    

             DO j=1, nvm ! loop over PFTs
              IF ( (delta_veg(j) > min_stomate) .AND. is_peat(j) .AND. pref_soil_veg(j)==jst .AND. natural(j) ) THEN !min_stomate

                ! Dilution of reservoirs
                ! Recalculate the litter and soil carbon with taking into accout the change in
                ! veget_cov_max (delta_veg)
                ! Litter
                litter(i,:,j,:,:)=(litter(i,:,j,:,:) * veget_cov_max_old(i,j) + dilu_lit(i,:,:,:) * delta_veg(j)/delta_veg_expand) &
                                  / veget_cov_max(i,j)
                fuel_1hr(i,j,:,:)=(fuel_1hr(i,j,:,:) * veget_cov_max_old(i,j) + dilu_f1hr(i,:,:) * delta_veg(j)/delta_veg_expand) & 
                                  / veget_cov_max(i,j)
                fuel_10hr(i,j,:,:)=(fuel_10hr(i,j,:,:) * veget_cov_max_old(i,j) + dilu_f10hr(i,:,:) * delta_veg(j)/delta_veg_expand) &
                                  / veget_cov_max(i,j)
                fuel_100hr(i,j,:,:)=(fuel_100hr(i,j,:,:) * veget_cov_max_old(i,j) + dilu_f100hr(i,:,:) * delta_veg(j)/delta_veg_expand) &
                                  / veget_cov_max(i,j)
                fuel_1000hr(i,j,:,:)=(fuel_1000hr(i,j,:,:) * veget_cov_max_old(i,j) + dilu_f1000hr(i,:,:) * delta_veg(j)/delta_veg_expand) &
                                  / veget_cov_max(i,j)
                !JCADD available and not available litter for grazing
                ! only not available litter change, available litter will not
                ! change, because tree litter can not be eaten

               carbon(i,:,j)=(carbon(i,:,j) * veget_cov_max_old(i,j) + dilu_soil_carbon(i,:) * delta_veg(j)/delta_veg_expand) &
                              / veget_cov_max(i,j)
               IF ( ok_pc ) THEN
                 !IF (carbon_old(i,iactive,j) .GT. min_stomate) THEN !min_stomate
                 deepC_a(i,:,j)=(deepC_a(i,:,j)*veget_cov_max_old(i,j)+dilu_soil_carbon_vertres(i,:,iactive)*delta_veg(j) &
                                    /delta_veg_expand)/veget_cov_max(i,j)
                 !ENDIF
                 !IF (carbon_old(i,islow,j) .GT. min_stomate) THEN !min_stomate
                 deepC_s(i,:,j)=(deepC_s(i,:,j)*veget_cov_max_old(i,j)+dilu_soil_carbon_vertres(i,:,islow)*delta_veg(j) &
                                   /delta_veg_expand)/veget_cov_max(i,j)
                 !ENDIF
                 !IF (carbon_old(i,ipassive,j) .GT. min_stomate) THEN !min_stomate
                 deepC_p(i,:,j)=(deepC_p(i,:,j)*veget_cov_max_old(i,j)+dilu_soil_carbon_vertres(i,:,ipassive)*delta_veg(j) &
                                  /delta_veg_expand)/veget_cov_max(i,j)
                 !ENDIF
               ENDIF

               TCarbon(i,j)=(TCarbon(i,j) * veget_cov_max_old(i,j) + dilu_TCarbon(i) * delta_veg(j)/delta_veg_expand) &
                            / veget_cov_max(i,j)

               turnover_daily(i,j,:,:)=(turnover_daily(i,j,:,:)*veget_cov_max_old(i,j)+dilu_turnover_daily(i,:,:)*delta_veg(j)&
                           /delta_veg_expand)/veget_cov_max(i,j)


               bm_to_litter(i,j,:,:)=(bm_to_litter(i,j,:,:)*veget_cov_max_old(i,j)+dilu_bm_to_litter(i,:,:)*delta_veg(j) &
                          /delta_veg_expand)/veget_cov_max(i,j)
               co2flux_new(i,j)=(co2flux_old(i,j)*veget_cov_max_old(i,j)+dilu_co2flux_new(i)*delta_veg(j) &
                         /delta_veg_expand)/veget_cov_max(i,j)
               gpp_daily(i,j)=(gpp_daily(i,j)*veget_cov_max_old(i,j)+dilu_gpp_daily(i)*delta_veg(j)/delta_veg_expand) &
                         /veget_cov_max(i,j)
               resp_growth(i,j)=(resp_growth(i,j)*veget_cov_max_old(i,j)+dilu_resp_growth(i)*delta_veg(j) & 
                        /delta_veg_expand)/veget_cov_max(i,j)
               resp_maint(i,j)=(resp_maint(i,j)*veget_cov_max_old(i,j)+dilu_resp_maint(i)*delta_veg(j)/delta_veg_expand) &
                        /veget_cov_max(i,j)
               resp_hetero(i,j)=(resp_hetero(i,j)*veget_cov_max_old(i,j)+dilu_resp_hetero(i)*delta_veg(j)/delta_veg_expand) &
                        /veget_cov_max(i,j)
  
               co2_to_bm(i,j)=(co2_to_bm(i,j)*veget_cov_max_old(i,j)+dilu_co2_to_bm(i)*delta_veg(j)/delta_veg_expand) &
                        /veget_cov_max(i,j)
               co2_fire(i,j)=(co2_fire(i,j)*veget_cov_max_old(i,j)+dilu_co2_fire(i)*delta_veg(j)/delta_veg_expand) &
                        /veget_cov_max(i,j)

              ENDIF  !delta_veg(j) > min_stomate
             ENDDO !j=1, nvm 

           ELSE   !delta_veg_shrink .GT. -min_stomate, no peat PFT shrinks
            IF ( veget_cov_max_old(i,ibare_sechiba) .GT. min_stomate) THEN ! expanded PFT get C from baresoil
              IF ( (veget_cov_max_old(i,ibare_sechiba) .GT. delta_veg_expand ) .AND. (veget_cov_max(i,ibare_sechiba) .GT. min_stomate) )THEN
                DO j=1, nvm
                  IF ( natural(j) .AND. .NOT. pasture(j) .AND. is_peat(j) .AND. pref_soil_veg(j)==jst .AND. (delta_veg(j) > min_stomate) ) THEN
                    litter(i,:,j,:,:)=(litter(i,:,j,:,:) * veget_cov_max_old(i,j) + litter(i,:,ibare_sechiba,:,:) & 
                                       *delta_veg(j) )/veget_cov_max(i,j)
                    fuel_1hr(i,j,:,:)=(fuel_1hr(i,j,:,:) * veget_cov_max_old(i,j) + fuel_1hr(i,ibare_sechiba,:,:) &
                                       *delta_veg(j) )/veget_cov_max(i,j)
                    fuel_10hr(i,j,:,:)=(fuel_10hr(i,j,:,:) * veget_cov_max_old(i,j) + fuel_10hr(i,ibare_sechiba,:,:) &
                                       *delta_veg(j) )/veget_cov_max(i,j)
                    fuel_100hr(i,j,:,:)=(fuel_100hr(i,j,:,:) * veget_cov_max_old(i,j) + fuel_100hr(i,ibare_sechiba,:,:) &
                                       *delta_veg(j) )/veget_cov_max(i,j)
                    fuel_1000hr(i,j,:,:)=(fuel_1000hr(i,j,:,:) * veget_cov_max_old(i,j) + fuel_1000hr(i,ibare_sechiba,:,:) &
                                       *delta_veg(j) )/veget_cov_max(i,j)
                    carbon(i,:,j)=(carbon(i,:,j) * veget_cov_max_old(i,j) + carbon(i,:,ibare_sechiba) &
                                       *delta_veg(j) )/veget_cov_max(i,j)
                    IF (ok_pc) THEN 
                      deepC_a(i,:,j)=(deepC_a(i,:,j)*veget_cov_max_old(i,j)+ deepC_a(i,:,ibare_sechiba) &
                                       *delta_veg(j) )/veget_cov_max(i,j)   
                      deepC_s(i,:,j)=(deepC_s(i,:,j)*veget_cov_max_old(i,j)+ deepC_s(i,:,ibare_sechiba) &
                                       *delta_veg(j) )/veget_cov_max(i,j)
                      deepC_p(i,:,j)=(deepC_p(i,:,j)*veget_cov_max_old(i,j)+ deepC_p(i,:,ibare_sechiba) &
                                       *delta_veg(j) )/veget_cov_max(i,j)
                    ENDIF
                    TCarbon(i,j)=(TCarbon(i,j) * veget_cov_max_old(i,j) + TCarbon(i,ibare_sechiba) &
                                       *delta_veg(j) )/veget_cov_max(i,j)
                    turnover_daily(i,j,:,:)=(turnover_daily(i,j,:,:)*veget_cov_max_old(i,j)+ &
                              turnover_daily(i,ibare_sechiba,:,:) *delta_veg(j) )/veget_cov_max(i,j) 
                    bm_to_litter(i,j,:,:)=(bm_to_litter(i,j,:,:)*veget_cov_max_old(i,j) + &
                              bm_to_litter(i,ibare_sechiba,:,:) *delta_veg(j) )/veget_cov_max(i,j)
                    co2flux_new(i,j)=(co2flux_old(i,j)*veget_cov_max_old(i,j)+ &
                              co2flux_old(i,ibare_sechiba) *delta_veg(j) )/veget_cov_max(i,j)
                    gpp_daily(i,j)=(gpp_daily(i,j)*veget_cov_max_old(i,j)+ &
                              gpp_daily(i,ibare_sechiba) *delta_veg(j) )/veget_cov_max(i,j)
                    resp_growth(i,j)=(resp_growth(i,j)*veget_cov_max_old(i,j)+ &
                              resp_growth(i,ibare_sechiba) *delta_veg(j) )/veget_cov_max(i,j)
                    resp_maint(i,j)=(resp_maint(i,j)*veget_cov_max_old(i,j)+ &
                              resp_maint(i,ibare_sechiba) *delta_veg(j) )/veget_cov_max(i,j)
                    resp_hetero(i,j)=(resp_hetero(i,j)*veget_cov_max_old(i,j)+ &
                              resp_hetero(i,ibare_sechiba) *delta_veg(j) )/veget_cov_max(i,j)
                    co2_to_bm(i,j)=(co2_to_bm(i,j)*veget_cov_max_old(i,j)+ &
                              co2_to_bm(i,ibare_sechiba) *delta_veg(j) )/veget_cov_max(i,j)
                    co2_fire(i,j)=(co2_fire(i,j)*veget_cov_max_old(i,j)+ &
                             co2_fire(i,ibare_sechiba) *delta_veg(j) )/veget_cov_max(i,j)          
                  ENDIF 
                ENDDO
                litter(i,:,ibare_sechiba,:,:) =(litter(i,:,ibare_sechiba,:,:)* veget_cov_max_old(i,ibare_sechiba)- &
                                       litter(i,:,ibare_sechiba,:,:)*delta_veg_expand)/veget_cov_max(i,ibare_sechiba)
                fuel_1hr(i,ibare_sechiba,:,:)=(fuel_1hr(i,ibare_sechiba,:,:)* veget_cov_max_old(i,ibare_sechiba)- &
                                       fuel_1hr(i,ibare_sechiba,:,:)*delta_veg_expand)/veget_cov_max(i,ibare_sechiba)
                fuel_10hr(i,ibare_sechiba,:,:)=(fuel_10hr(i,ibare_sechiba,:,:)* veget_cov_max_old(i,ibare_sechiba)- &
                                       fuel_10hr(i,ibare_sechiba,:,:)*delta_veg_expand)/veget_cov_max(i,ibare_sechiba)
                fuel_100hr(i,ibare_sechiba,:,:)=(fuel_100hr(i,ibare_sechiba,:,:)* veget_cov_max_old(i,ibare_sechiba)- &
                                       fuel_100hr(i,ibare_sechiba,:,:)*delta_veg_expand)/veget_cov_max(i,ibare_sechiba)
                fuel_1000hr(i,ibare_sechiba,:,:)=(fuel_1000hr(i,ibare_sechiba,:,:)* veget_cov_max_old(i,ibare_sechiba)- &
                                       fuel_1000hr(i,ibare_sechiba,:,:)*delta_veg_expand)/veget_cov_max(i,ibare_sechiba)

                carbon(i,:,ibare_sechiba)=(carbon(i,:,ibare_sechiba)* veget_cov_max_old(i,ibare_sechiba)- &
                                          carbon(i,:,ibare_sechiba)*delta_veg_expand)/veget_cov_max(i,ibare_sechiba)
                IF (ok_pc) THEN
                  deepC_a(i,:,ibare_sechiba)=(deepC_a(i,:,ibare_sechiba)* veget_cov_max_old(i,ibare_sechiba)- &
                                          deepC_a(i,:,ibare_sechiba)*delta_veg_expand)/veget_cov_max(i,ibare_sechiba)
   
                  deepC_s(i,:,ibare_sechiba)=(deepC_s(i,:,ibare_sechiba)* veget_cov_max_old(i,ibare_sechiba)- &
                                          deepC_s(i,:,ibare_sechiba)*delta_veg_expand)/veget_cov_max(i,ibare_sechiba)
                  deepC_p(i,:,ibare_sechiba)=(deepC_p(i,:,ibare_sechiba)* veget_cov_max_old(i,ibare_sechiba)- &
                                          deepC_p(i,:,ibare_sechiba)*delta_veg_expand)/veget_cov_max(i,ibare_sechiba)
                ENDIF
                TCarbon(i,ibare_sechiba)=(TCarbon(i,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba)- &
                                         TCarbon(i,ibare_sechiba)*delta_veg_expand)/veget_cov_max(i,ibare_sechiba)
                turnover_daily(i,ibare_sechiba,:,:)=(turnover_daily(i,ibare_sechiba,:,:)*veget_cov_max_old(i,ibare_sechiba)- & 
                                         turnover_daily(i,ibare_sechiba,:,:)*delta_veg_expand)/veget_cov_max(i,ibare_sechiba)
                bm_to_litter(i,ibare_sechiba,:,:)=(bm_to_litter(i,ibare_sechiba,:,:)*veget_cov_max_old(i,ibare_sechiba)- &
                                         bm_to_litter(i,ibare_sechiba,:,:)*delta_veg_expand)/veget_cov_max(i,ibare_sechiba)
                co2flux_new(i,ibare_sechiba)=(co2flux_new(i,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba)- &
                                         co2flux_new(i,ibare_sechiba)*delta_veg_expand)/veget_cov_max(i,ibare_sechiba)
                gpp_daily(i,ibare_sechiba)=(gpp_daily(i,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba)- &
                                         gpp_daily(i,ibare_sechiba)*delta_veg_expand)/veget_cov_max(i,ibare_sechiba)
                resp_growth(i,ibare_sechiba)=(resp_growth(i,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba)- &
                                         resp_growth(i,ibare_sechiba)*delta_veg_expand)/veget_cov_max(i,ibare_sechiba)
                resp_maint(i,ibare_sechiba)=(resp_maint(i,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba)- &
                                         resp_maint(i,ibare_sechiba)*delta_veg_expand)/veget_cov_max(i,ibare_sechiba)
                resp_hetero(i,ibare_sechiba)=(resp_hetero(i,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba) - &
                                         resp_hetero(i,ibare_sechiba)*delta_veg_expand)/veget_cov_max(i,ibare_sechiba)
                co2_to_bm(i,ibare_sechiba)=(co2_to_bm(i,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba) - &
                                         co2_to_bm(i,ibare_sechiba)*delta_veg_expand)/veget_cov_max(i,ibare_sechiba)
                co2_fire(i,ibare_sechiba)=(co2_fire(i,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba) - &
                                         co2_fire(i,ibare_sechiba)*delta_veg_expand)/veget_cov_max(i,ibare_sechiba)
              ELSE ! expanded PFT get C from baresoil, baresoil C depleted

                DO j=1, nvm
                  IF ( natural(j) .AND. .NOT. pasture(j) .AND. is_peat(j) .AND. pref_soil_veg(j)==jst .AND. (delta_veg(j) > min_stomate) ) THEN
                    litter(i,:,j,:,:)=(litter(i,:,j,:,:) * veget_cov_max_old(i,j) + litter(i,:,ibare_sechiba,:,:) &
                                       *veget_cov_max_old(i,ibare_sechiba)/delta_veg_expand*delta_veg(j) )/veget_cov_max(i,j)
                    fuel_1hr(i,j,:,:)=(fuel_1hr(i,j,:,:) * veget_cov_max_old(i,j) + fuel_1hr(i,ibare_sechiba,:,:) &
                                       *veget_cov_max_old(i,ibare_sechiba)/delta_veg_expand*delta_veg(j) )/veget_cov_max(i,j)
                    fuel_10hr(i,j,:,:)=(fuel_10hr(i,j,:,:) * veget_cov_max_old(i,j) + fuel_10hr(i,ibare_sechiba,:,:) &
                                       *veget_cov_max_old(i,ibare_sechiba)/delta_veg_expand*delta_veg(j) )/veget_cov_max(i,j)
                    fuel_100hr(i,j,:,:)=(fuel_100hr(i,j,:,:) * veget_cov_max_old(i,j) + fuel_100hr(i,ibare_sechiba,:,:) &
                                       *veget_cov_max_old(i,ibare_sechiba)/delta_veg_expand*delta_veg(j) )/veget_cov_max(i,j)
                    fuel_1000hr(i,j,:,:)=(fuel_1000hr(i,j,:,:) * veget_cov_max_old(i,j) + fuel_1000hr(i,ibare_sechiba,:,:) &
                                       *veget_cov_max_old(i,ibare_sechiba)/delta_veg_expand*delta_veg(j) )/veget_cov_max(i,j)
                    carbon(i,:,j)=(carbon(i,:,j) * veget_cov_max_old(i,j) + carbon(i,:,ibare_sechiba) &
                                       *veget_cov_max_old(i,ibare_sechiba)/delta_veg_expand*delta_veg(j) )/veget_cov_max(i,j)
                    IF (ok_pc) THEN
                      deepC_a(i,:,j)=(deepC_a(i,:,j)*veget_cov_max_old(i,j)+ deepC_a(i,:,ibare_sechiba) &
                                       *veget_cov_max_old(i,ibare_sechiba)/delta_veg_expand*delta_veg(j) )/veget_cov_max(i,j)
                      deepC_s(i,:,j)=(deepC_s(i,:,j)*veget_cov_max_old(i,j)+ deepC_s(i,:,ibare_sechiba) &
                                       *veget_cov_max_old(i,ibare_sechiba)/delta_veg_expand*delta_veg(j) )/veget_cov_max(i,j)
                      deepC_p(i,:,j)=(deepC_p(i,:,j)*veget_cov_max_old(i,j)+ deepC_p(i,:,ibare_sechiba) &
                                       *veget_cov_max_old(i,ibare_sechiba)/delta_veg_expand*delta_veg(j) )/veget_cov_max(i,j)
                    ENDIF
                    TCarbon(i,j)=(TCarbon(i,j) * veget_cov_max_old(i,j) + TCarbon(i,ibare_sechiba) &
                                       *veget_cov_max_old(i,ibare_sechiba)/delta_veg_expand*delta_veg(j) )/veget_cov_max(i,j)
                    turnover_daily(i,j,:,:)=(turnover_daily(i,j,:,:)*veget_cov_max_old(i,j) + turnover_daily(i,ibare_sechiba,:,:)&
                                       *veget_cov_max_old(i,ibare_sechiba)/delta_veg_expand*delta_veg(j) )/veget_cov_max(i,j)
                    bm_to_litter(i,j,:,:)=(bm_to_litter(i,j,:,:)*veget_cov_max_old(i,j) + bm_to_litter(i,ibare_sechiba,:,:) &
                                       *veget_cov_max_old(i,ibare_sechiba)/delta_veg_expand*delta_veg(j) )/veget_cov_max(i,j)
                    co2flux_new(i,j)=(co2flux_old(i,j)*veget_cov_max_old(i,j)+ co2flux_old(i,ibare_sechiba) &
                                       *veget_cov_max_old(i,ibare_sechiba)/delta_veg_expand*delta_veg(j) )/veget_cov_max(i,j)
                    gpp_daily(i,j)=(gpp_daily(i,j)*veget_cov_max_old(i,j)+ gpp_daily(i,ibare_sechiba) &
                                       *veget_cov_max_old(i,ibare_sechiba)/delta_veg_expand*delta_veg(j) )/veget_cov_max(i,j)
                    resp_growth(i,j)=(resp_growth(i,j)*veget_cov_max_old(i,j)+ resp_growth(i,ibare_sechiba) &
                                       *veget_cov_max_old(i,ibare_sechiba)/delta_veg_expand*delta_veg(j) )/veget_cov_max(i,j)
                    resp_maint(i,j)=(resp_maint(i,j)*veget_cov_max_old(i,j)+ resp_maint(i,ibare_sechiba) &
                                       *veget_cov_max_old(i,ibare_sechiba)/delta_veg_expand*delta_veg(j) )/veget_cov_max(i,j)
                    resp_hetero(i,j)=(resp_hetero(i,j)*veget_cov_max_old(i,j)+ resp_hetero(i,ibare_sechiba) &
                                       *veget_cov_max_old(i,ibare_sechiba)/delta_veg_expand*delta_veg(j) )/veget_cov_max(i,j)
                    co2_to_bm(i,j)=(co2_to_bm(i,j)*veget_cov_max_old(i,j)+ co2_to_bm(i,ibare_sechiba) &
                                       *veget_cov_max_old(i,ibare_sechiba)/delta_veg_expand*delta_veg(j) )/veget_cov_max(i,j)
                    co2_fire(i,j)=(co2_fire(i,j)*veget_cov_max_old(i,j)+ co2_fire(i,ibare_sechiba) &
                                       *veget_cov_max_old(i,ibare_sechiba)/delta_veg_expand*delta_veg(j) )/veget_cov_max(i,j)
                  ENDIF
                ENDDO
                litter(i,:,ibare_sechiba,:,:) = zero
                fuel_1hr(i,ibare_sechiba,:,:) = zero
                fuel_10hr(i,ibare_sechiba,:,:) = zero
                fuel_100hr(i,ibare_sechiba,:,:) = zero
                fuel_1000hr(i,ibare_sechiba,:,:) = zero
                carbon(i,:,ibare_sechiba) = zero
                IF (ok_pc) THEN
                  deepC_a(i,:,ibare_sechiba) = zero
                  deepC_s(i,:,ibare_sechiba) = zero
                  deepC_p(i,:,ibare_sechiba) = zero
                ENDIF
                TCarbon(i,ibare_sechiba) = zero
                turnover_daily(i,ibare_sechiba,:,:) = zero
                bm_to_litter(i,ibare_sechiba,:,:) = zero
                co2flux_new(i,ibare_sechiba) = zero
                gpp_daily(i,ibare_sechiba) = zero
                resp_growth(i,ibare_sechiba) = zero
                resp_maint(i,ibare_sechiba) = zero
                resp_hetero(i,ibare_sechiba) = zero
                co2_to_bm(i,ibare_sechiba) = zero
                co2_fire(i,ibare_sechiba) = zero
              ENDIF   
            ELSE  !expanded PFT do not get C

             DO j=1, nvm
               IF ( natural(j) .AND. .NOT. pasture(j) .AND. is_peat(j) .AND. pref_soil_veg(j)==jst .AND. (delta_veg(j) > min_stomate) ) THEN
                  litter(i,:,j,:,:)=litter(i,:,j,:,:) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                  fuel_1hr(i,j,:,:)=fuel_1hr(i,j,:,:) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                  fuel_10hr(i,j,:,:)=fuel_10hr(i,j,:,:) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                  fuel_100hr(i,j,:,:)=fuel_100hr(i,j,:,:) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                  fuel_1000hr(i,j,:,:)=fuel_1000hr(i,j,:,:) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                  carbon(i,:,j)=carbon(i,:,j) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                  IF (ok_pc) THEN  
                    deepC_a(i,:,j)=deepC_a(i,:,j) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                    deepC_s(i,:,j)=deepC_s(i,:,j) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                    deepC_p(i,:,j)=deepC_p(i,:,j) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                  ENDIF
                  TCarbon(i,j)=TCarbon(i,j) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                  turnover_daily(i,j,:,:)=turnover_daily(i,j,:,:) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                  bm_to_litter(i,j,:,:)=bm_to_litter(i,j,:,:) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                  co2flux_new(i,j)=co2flux_new(i,j) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                  gpp_daily(i,j)=gpp_daily(i,j) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                  resp_growth(i,j)=resp_growth(i,j) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                  resp_maint(i,j)=resp_maint(i,j) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                  resp_hetero(i,j)=resp_hetero(i,j) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                  co2_to_bm(i,j)=co2_to_bm(i,j) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                  co2_fire(i,j)=co2_fire(i,j) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
               ENDIF
             ENDDO

            ENDIF   

           ENDIF 
  
          ELSE !delta_veg_expand .LT. min_stomate, no peat PFT expands
           IF ( veget_cov_max(i,ibare_sechiba) .GT. min_stomate) THEN !shrinked PFT give C to baresoil
               litter(i,:,ibare_sechiba,:,:)=(litter(i,:,ibare_sechiba,:,:) * veget_cov_max_old(i,ibare_sechiba) + dilu_lit(i,:,:,:) ) &
                                  / veget_cov_max(i,ibare_sechiba)
               fuel_1hr(i,ibare_sechiba,:,:)=(fuel_1hr(i,ibare_sechiba,:,:) * veget_cov_max_old(i,ibare_sechiba) + dilu_f1hr(i,:,:) ) &
                                  / veget_cov_max(i,ibare_sechiba)
               fuel_10hr(i,ibare_sechiba,:,:)=(fuel_10hr(i,ibare_sechiba,:,:) * veget_cov_max_old(i,ibare_sechiba) + dilu_f10hr(i,:,:) ) &
                                  / veget_cov_max(i,ibare_sechiba)
               fuel_100hr(i,ibare_sechiba,:,:)=(fuel_100hr(i,ibare_sechiba,:,:) * veget_cov_max_old(i,ibare_sechiba) + dilu_f100hr(i,:,:) ) &
                                  / veget_cov_max(i,ibare_sechiba)
               fuel_1000hr(i,ibare_sechiba,:,:)=(fuel_1000hr(i,ibare_sechiba,:,:) * veget_cov_max_old(i,ibare_sechiba) + dilu_f1000hr(i,:,:) ) &
                                  / veget_cov_max(i,ibare_sechiba)

               carbon(i,:,ibare_sechiba)=(carbon(i,:,ibare_sechiba) * veget_cov_max_old(i,ibare_sechiba) + dilu_soil_carbon(i,:) ) &
                              / veget_cov_max(i,ibare_sechiba)
               IF ( ok_pc ) THEN
                 !IF (carbon_old(i,iactive,ibare_sechiba) .GT. min_stomate) THEN !min_stomate
                 deepC_a(i,:,ibare_sechiba)=(deepC_a(i,:,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba)+dilu_soil_carbon_vertres(i,:,iactive) ) &
                                    /veget_cov_max(i,ibare_sechiba)
              !   ENDIF
              !   IF (carbon_old(i,islow,ibare_sechiba) .GT. min_stomate) THEN !min_stomate
                 deepC_s(i,:,ibare_sechiba)=(deepC_s(i,:,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba)+dilu_soil_carbon_vertres(i,:,islow) ) &
                                   /veget_cov_max(i,ibare_sechiba)
              !   ENDIF
              !   IF (carbon_old(i,ipassive,ibare_sechiba) .GT. min_stomate) THEN !min_stomate
                 deepC_p(i,:,ibare_sechiba)=(deepC_p(i,:,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba)+dilu_soil_carbon_vertres(i,:,ipassive) ) &
                                  /veget_cov_max(i,ibare_sechiba)
              !   ENDIF
               ENDIF

               TCarbon(i,ibare_sechiba)=(TCarbon(i,ibare_sechiba) * veget_cov_max_old(i,ibare_sechiba) + dilu_TCarbon(i) ) &
                            /veget_cov_max(i,ibare_sechiba)

               turnover_daily(i,ibare_sechiba,:,:)=(turnover_daily(i,ibare_sechiba,:,:)*veget_cov_max_old(i,ibare_sechiba)+dilu_turnover_daily(i,:,:) )&
                           /veget_cov_max(i,ibare_sechiba)


               bm_to_litter(i,ibare_sechiba,:,:)=(bm_to_litter(i,ibare_sechiba,:,:)*veget_cov_max_old(i,ibare_sechiba)+dilu_bm_to_litter(i,:,:) ) &
                          /veget_cov_max(i,ibare_sechiba)
               co2flux_new(i,ibare_sechiba)=(co2flux_old(i,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba)+dilu_co2flux_new(i) ) &
                         /veget_cov_max(i,ibare_sechiba)
               gpp_daily(i,ibare_sechiba)=(gpp_daily(i,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba)+dilu_gpp_daily(i) ) &
                         /veget_cov_max(i,ibare_sechiba)
               resp_growth(i,ibare_sechiba)=(resp_growth(i,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba)+dilu_resp_growth(i) ) &
                         /veget_cov_max(i,ibare_sechiba)
               resp_maint(i,ibare_sechiba)=(resp_maint(i,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba)+dilu_resp_maint(i) ) &
                        /veget_cov_max(i,ibare_sechiba)
               resp_hetero(i,ibare_sechiba)=(resp_hetero(i,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba)+dilu_resp_hetero(i) ) &
                        /veget_cov_max(i,ibare_sechiba)

               co2_to_bm(i,ibare_sechiba)=(co2_to_bm(i,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba)+dilu_co2_to_bm(i) ) &
                        /veget_cov_max(i,ibare_sechiba)
               co2_fire(i,ibare_sechiba)=(co2_fire(i,ibare_sechiba)*veget_cov_max_old(i,ibare_sechiba)+dilu_co2_fire(i) ) &
                        /veget_cov_max(i,ibare_sechiba)
  
           ELSE 
             DO j=1, nvm
               IF ( natural(j) .AND. .NOT. pasture(j) .AND. is_peat(j) .AND. pref_soil_veg(j)==jst .AND. (delta_veg(j) < -min_stomate) ) THEN
                 IF (veget_cov_max(i,j) .GT. min_stomate) THEN !shrinked PFT keep their own C
                    litter(i,:,j,:,:)=litter(i,:,j,:,:) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                    fuel_1hr(i,j,:,:)=fuel_1hr(i,j,:,:) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                    fuel_10hr(i,j,:,:)=fuel_10hr(i,j,:,:) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                    fuel_100hr(i,j,:,:)=fuel_100hr(i,j,:,:) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                    fuel_1000hr(i,j,:,:)=fuel_1000hr(i,j,:,:) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                    carbon(i,:,j)=carbon(i,:,j) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                    deepC_a(i,:,j)=deepC_a(i,:,j) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                    deepC_s(i,:,j)=deepC_s(i,:,j) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                    deepC_p(i,:,j)=deepC_p(i,:,j) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                    TCarbon(i,j)=TCarbon(i,j) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                    turnover_daily(i,j,:,:)=turnover_daily(i,j,:,:) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                    bm_to_litter(i,j,:,:)=bm_to_litter(i,j,:,:) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                    co2flux_new(i,j)=co2flux_new(i,j) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                    gpp_daily(i,j)=gpp_daily(i,j) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                    resp_growth(i,j)=resp_growth(i,j) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                    resp_maint(i,j)=resp_maint(i,j) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                    resp_hetero(i,j)=resp_hetero(i,j) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                    co2_to_bm(i,j)=co2_to_bm(i,j) * veget_cov_max_old(i,j)/veget_cov_max(i,j)
                    co2_fire(i,j)=co2_fire(i,j) * veget_cov_max_old(i,j)/veget_cov_max(i,j)

                 ELSE
                      WRITE (numout,*) 'QCJ check Peat PFTs shrinked to 0, in lpj_cover, at: ',lalo(i,:)
                      CALL ipslerr_p(3,'lpj_cover','ok_dgvm_peat=T','Peat PFTs shrinked to zero,but no peat PFT expanded, no baresoil, shrinked C have no where to go','')    
                 ENDIF       
 
               ENDIF   
             ENDDO
           ENDIF 

           DO j=1, nvm
             IF ((delta_veg(j) > min_stomate) .AND. is_grassland_manag(j) .AND. is_grassland_grazed(j)) THEN
                litter_avail(i,:,j) = litter_avail(i,:,j) * veget_cov_max_old(i,j) / veget_cov_max(i,j)
                litter_not_avail(i,:,j) = litter(i,:,j,iabove,icarbon) - litter_avail(i,:,j)
             ENDIF
           ENDDO

          ENDIF !delta_veg_expand
    
          DO j=1, nvm 
            IF ( (veget_cov_max(i,j).GT. min_stomate) .AND. is_peat(j) .AND. natural(j) .AND. pref_soil_veg(j)==jst ) THEN !min_stomate

                ! Correct biomass densities to conserve mass
                ! since it's defined on veget_cov_max
                biomass(i,j,:,:) = biomass(i,j,:,:) * veget_cov_max_old(i,j) / veget_cov_max(i,j)
            ENDIF
          ENDDO ! loop over PFTs

       ENDDO ! loop over grid points
  
      ENDIF ! wettile_dgvm

     ENDDO ! jst=1,nstm    
     
     !DO i=1,npts
      ! WRITE (numout,*) 'QCJ check lpj_cover,veget_cov_max,',veget_cov_max(i,:)
     !ENDDO


      C_end(:,1)=SUM(SUM(biomass(:,:,:,icarbon),DIM=3) * veget_cov_max,DIM=2)
      C_end(:,2)=SUM(SUM(SUM(litter(:,:,:,:,icarbon),DIM=4),DIM=2)*veget_cov_max,DIM=2)
      C_end(:,3)=SUM(SUM(carbon(:,:,:),DIM=2) * veget_cov_max,DIM=2)

      Fin_end(:,1)=SUM(SUM(bm_to_litter(:,:,:,icarbon),DIM=3) * veget_cov_max,DIM=2)
      Fin_end(:,2)=SUM(SUM(turnover_daily(:,:,:,icarbon),DIM=3) * veget_cov_max,DIM=2)
      Fin_end(:,3)=SUM(co2_to_bm * veget_cov_max,DIM=2)
      Fin_end(:,4)=SUM(gpp_daily* veget_cov_max,DIM=2)

      Fout_end(:,1)=SUM(co2_fire* veget_cov_max,DIM=2)
      Fout_end(:,2)=SUM(resp_hetero* veget_cov_max,DIM=2)
      Fout_end(:,3)=SUM(resp_maint* veget_cov_max,DIM=2)
      Fout_end(:,4)=SUM(resp_growth* veget_cov_max,DIM=2)

      DO i = 1, npts ! loop over grid points
         C_check(i,:)=C_end(i,:)-C_start(i,:)
         !WRITE (numout,*) 'QCJ check lpj_cover, veget_cov_max old:',veget_cov_max_old(i,:)
         !WRITE (numout,*) 'QCJ check lpj_cover, veget_cov_max new:',veget_cov_max(i,:)
         !WRITE (numout,*) 'QCJ check lpj_cover, delta veget_cov_max:',veget_cov_max(i,:)-veget_cov_max_old(i,:)
  
         !IF ( ANY(C_check(i,:) .GT. min_stomate*100.) ) THEN
            !WRITE (numout,*) 'QCJ NOTE, C pool mass not conserved after adjust VEGET_COV_MAX,at grid',i
            !WRITE (numout,*) 'QCJ lpj_cover, biomass, litter, SOC',C_check(i,:)
         !ENDIF 

         Fin_check(i,:)=Fin_end(i,:)-Fin_start(i,:)
         !IF ( ANY(Fin_check(i,:) .GT. min_stomate*100.) ) THEN
            !WRITE (numout,*) 'QCJ NOTE, Fin mass not conserved after adjust VEGET_COV_MAX, at grid',i,'Coordinates:',lalo(i,:)
            !WRITE (numout,*) 'QCJ lpj_cover, bm_to_litter, turnover_daily, co2_to_bm, gpp_daily',Fin_check(i,:)
         !ENDIF

         !Fout_check(i,:)=Fout_end(i,:)-Fout_start(i,:)
         !IF ( ANY(Fout_check(i,:) .GT. min_stomate*100.) ) THEN
            !WRITE (numout,*) 'QCJ NOTE, Fout mass not conserved after adjust VEGET_COV_MAX, at grid',i
            !WRITE (numout,*) 'QCJ lpj_cover, co2_fire, resp_hetero, resp_maint,resp_growth ',Fout_check(i,:)
         !ENDIF
      ENDDO 

      vartmp(:)=SUM(co2flux_new*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tCO2FLUX", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(co2flux_old*veget_cov_max_old,dim=2)
      CALL histwrite_p (hist_id_stomate, "tCO2FLUX_OLD", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(TCarbon*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tCARBON", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(gpp_daily*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tGPP", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(resp_growth*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tRESP_GROWTH", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(resp_maint*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tRESP_MAINT", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(resp_hetero*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tRESP_HETERO", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(co2_to_bm*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tCO2_TAKEN", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(co2_fire*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tCO2_FIRE", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(SUM(biomass(:,:,:,icarbon),dim=3)*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tBIOMASS", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(SUM(SUM(litter(:,:,:,:,icarbon),dim=4),dim=2)*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tLITTER", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(SUM(fuel_1hr(:,:,:,icarbon),dim=3)*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tFUEL1HR", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(SUM(fuel_10hr(:,:,:,icarbon),dim=3)*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tFUEL10HR", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(SUM(fuel_100hr(:,:,:,icarbon),dim=3)*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tFUEL100HR", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(SUM(fuel_1000hr(:,:,:,icarbon),dim=3)*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tFUEL1000HR", itime, vartmp, npts, hori_index)
      vartmp(:)=SUM(SUM(carbon,dim=2)*veget_cov_max,dim=2)
      CALL histwrite_p (hist_id_stomate, "tSOILC", itime, vartmp, npts, hori_index)

      IF ( ok_pc ) THEN
        vartmp(:)=SUM(SUM(deepC_a,dim=2)*veget_cov_max,dim=2)
        CALL histwrite_p (hist_id_stomate, "tDEEPCa", itime, vartmp, npts, hori_index)
        vartmp(:)=SUM(SUM(deepC_s,dim=2)*veget_cov_max,dim=2)
        CALL histwrite_p (hist_id_stomate, "tDEEPCs", itime, vartmp, npts, hori_index)
        vartmp(:)=SUM(SUM(deepC_p,dim=2)*veget_cov_max,dim=2)
        CALL histwrite_p (hist_id_stomate, "tDEEPCp", itime, vartmp, npts, hori_index)
      ENDIF

  ENDIF  !IF ok_dgvm_peat

  END SUBROUTINE cover

END MODULE lpj_cover
