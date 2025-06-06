! =================================================================================================================================
! MODULE       : pft_parameters
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2011)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        This module initializes all the pft parameters in function of the
!!              number of vegetation types and of the values chosen by the user.
!!
!!\n DESCRIPTION:  This module allocates and initializes the pft parameters in function of the number of pfts
!!                 and the values of the parameters. \n
!!                 The number of PFTs is read in control.f90 (subroutine control_initialize). \n
!!                 Then we can initialize the parameters. \n
!!                 This module is the result of the merge of constantes_co2, constantes_veg, stomate_constants.\n
!!
!! RECENT CHANGE(S): Josefine Ghattas 2013 : The declaration part has been extracted and moved to module pft_parameters_var
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: $
!! $Date: 2018-05-30 15:40:04 +0200 (Wed, 30 May 2018) $
!! $Revision: 5268 $
!! \n
!_ ================================================================================================================================

MODULE pft_parameters

  USE pft_parameters_var
  USE vertical_soil_var
  USE constantes_mtc
  USE constantes_soil_var !!! only for nstm & pref_soil_veg
  USE ioipsl
  USE ioipsl_para 
  USE defprec

  IMPLICIT NONE

CONTAINS
  

!! ================================================================================================================================
!! SUBROUTINE   : pft_parameters_main
!!
!>\BRIEF          This subroutine initializes all the pft parameters in function of the
!! number of vegetation types chosen by the user.
!!
!! DESCRIPTION  : This subroutine is called after the reading of the number of PFTS and the options 
!!                activated by the user in the configuration files. \n
!!                The allocation is done just before reading the correspondence table  between PFTs and MTCs
!!                defined by the user in the configuration file.\n
!!                With the correspondence table, the subroutine can initialize the pft parameters in function
!!                of the flags activated (ok_sechiba, ok_stomate, ok_co2, routing, new_hydrol...) in order to
!!                optimize the memory allocation. \n
!!                If the number of PFTs and pft_to_mtc are not found, the standard configuration will be used
!!                (13 PFTs, PFT = MTC). \n 
!!                Some restrictions : the pft 1 can only be the bare soil and it is unique. \n
!!                Algorithm : Build new PFT from 13 generic-PFT or meta-classes.
!!                1. Read the number of PFTs in "run.def". If nothing is found, it is assumed that the user intend to use 
!!                   the standard of PFTs (13).
!!                2. Read the index vector in "run.def". The index vector associates one PFT to one meta-classe (or generic PFT).
!!                   When the association is done, the PFT defined by the user inherited the default values from the meta classe.
!!                   If nothing is found, it is assumed to use the standard index vector (PFT = MTC).
!!                3. Check consistency
!!                4. Memory allocation and initialization.
!!                5. The parameters are read in the configuration file in config_initialize (control module).
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE pft_parameters_main()

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.4 Local variables  

    INTEGER(i_std) :: j                             !! Index (unitless)

    !_ ================================================================================================================================ 

    !
    ! PFT global
    !

    IF(l_first_pft_parameters) THEN

       !! 1. First time step
       IF(printlev>=3) THEN
          WRITE(numout,*) 'l_first_pft_parameters :we read the parameters from the def files'
       ENDIF

       !! 2. Memory allocation for the pfts-parameters
       CALL pft_parameters_alloc()

       !! 3. Correspondance table 

       !! 3.1 Initialisation of the correspondance table
       !! Initialisation of the correspondance table
       DO j = 1, nvm
          pft_to_mtc(j) = j
       ENDDO ! j=1, nvm

       !! 3.2 Reading of the conrrespondance table in the .def file
       !
       !Config Key   = PFT_TO_MTC
       !Config Desc  = correspondance array linking a PFT to MTC
       !Config if    = OK_SECHIBA or OK_STOMATE
       !Config Def   = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13
       !Config Help  =
       !Config Units = [-]
       CALL getin_p('PFT_TO_MTC',pft_to_mtc)

       !! 3.3 If the user want to use the standard configuration, he needn't to fill the correspondance array
       !!     If the configuration is wrong, send a error message to the user.
       IF(nvm /= nvmc ) THEN
          !
          IF(pft_to_mtc(1) == undef_int) THEN
             STOP ' The array PFT_TO_MTC is empty : we stop'
          ENDIF !(pft_to_mtc(1) == undef_int)
          !
       ENDIF !(nvm /= nvmc )

       !! 3.4 Some error messages

       !! 3.4.1 What happened if pft_to_mtc(j) > nvmc or pft_to_mtc(j) <=0 (if the mtc doesn't exist)?
       DO j = 1, nvm ! Loop over # PFTs  
          !
          IF( (pft_to_mtc(j) > nvmc) .OR. (pft_to_mtc(j) <= 0) ) THEN
             WRITE(numout,*) 'the metaclass chosen does not exist'
             STOP 'we stop reading pft_to_mtc'
          ENDIF !( (pft_to_mtc(j) > nvmc) .OR. (pft_to_mtc(j) <= 0) )
          !
       ENDDO  ! Loop over # PFTs  


       !! 3.4.2 Check if pft_to_mtc(1) = 1 
       IF(pft_to_mtc(1) /= 1) THEN
          !
          WRITE(numout,*) 'the first pft has to be the bare soil'
          STOP 'we stop reading next values of pft_to_mtc'
          !
       ELSE
          !
          DO j = 2,nvm ! Loop over # PFTs different from bare soil
             !
             IF(pft_to_mtc(j) == 1) THEN
                WRITE(numout,*) 'only pft_to_mtc(1) has to be the bare soil'
                STOP 'we stop reading pft_to_mtc'
             ENDIF ! (pft_to_mtc(j) == 1)
             !
          ENDDO ! Loop over # PFTs different from bare soil
          !
       ENDIF !(pft_to_mtc(1) /= 1)


       !! 4.Initialisation of the pfts-parameters
       CALL pft_parameters_init()

       !! 5. Useful data

       !! 5.1 Read the name of the PFTs given by the user
       !
       !Config Key   = PFT_NAME
       !Config Desc  = Name of a PFT
       !Config if    = OK_SECHIBA or OK_STOMATE
       !Config Def   = bare ground, tropical broad-leaved evergreen, tropical broad-leaved raingreen, 
       !Config         temperate needleleaf evergreen, temperate broad-leaved evergreen temperate broad-leaved summergreen,
       !Config         boreal needleleaf evergreen, boreal broad-leaved summergreen, boreal needleleaf summergreen,
       !Config         C3 grass, C4 grass, C3 agriculture, C4 agriculture    
       !Config Help  = the user can name the new PFTs he/she introducing for new species
       !Config Units = [-]
       CALL getin_p('PFT_NAME',pft_name)

       !! 5.2 A useful message to the user: correspondance between the number of the pft
       !! and the name of the associated mtc 
       IF (printlev >=1 ) THEN
          WRITE(numout,*) ''
          DO j = 2,nvm ! Loop over # PFTs
             WRITE(numout,*) 'The PFT',j, 'called ', trim(PFT_name(j)),' corresponds to the MTC : ',trim(MTC_name(pft_to_mtc(j)))
          END DO
          WRITE(numout,*) ''
       END IF


       !! 6. End message
       IF (printlev>=3) WRITE(numout,*) 'pft_parameters_done'

       !! 8. Reset flag
       l_first_pft_parameters = .FALSE.

    ELSE 

       RETURN

    ENDIF !(l_first_pft_parameters)

  END SUBROUTINE pft_parameters_main


!! ================================================================================================================================
!! SUBROUTINE   : pft_parameters_init 
!!
!>\BRIEF          This subroutine initializes all the pft parameters by the default values
!! of the corresponding metaclasse. 
!!
!! DESCRIPTION  : This subroutine is called after the reading of the number of PFTS and the correspondence
!!                table defined by the user in the configuration files. \n
!!                With the correspondence table, the subroutine can search the default values for the parameter
!!                even if the PFTs are classified in a random order (except bare soil). \n
!!                With the correspondence table, the subroutine can initialize the pft parameters in function
!!                of the flags activated (ok_sechiba, ok_stomate, ok_co2, routing, new_hydrol...).\n
!!
!! RECENT CHANGE(S): Didier Solyga : Simplified PFT loops : use vector notation. 
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE pft_parameters_init()

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.1 Input variables
    INTEGER(i_std)                :: jv            !! Index (unitless) 
    !_ ================================================================================================================================ 

   ! 1.1 For parameters used anytime
   PFT_name(:) = MTC_name(pft_to_mtc(:))
   !
   ! Vegetation structure 
   !
   veget_ori_fixed_test_1(:) = veget_ori_fixed_mtc(pft_to_mtc(:))
   llaimax(:) = llaimax_mtc(pft_to_mtc(:))
   llaimin(:) = llaimin_mtc(pft_to_mtc(:))
   height_presc(:) = height_presc_mtc(pft_to_mtc(:))
   z0_over_height(:) = z0_over_height_mtc(pft_to_mtc(:))
   ratio_z0m_z0h(:) = ratio_z0m_z0h_mtc(pft_to_mtc(:))
   type_of_lai(:) = type_of_lai_mtc(pft_to_mtc(:))
   natural(:) = natural_mtc(pft_to_mtc(:))
! dgvmjc
    pasture(:) = pasture_mtc(pft_to_mtc(:))
! end dgvmjc
   !
   ! Water - sechiba
   !
   IF (zmaxh == 2.0) THEN
      WRITE(numout,*)'Initialize humcst using reference values for 2m soil depth'
      humcste(:) = humcste_ref2m(pft_to_mtc(:))  ! values for 2m soil depth
   ELSE IF (zmaxh == 4.0) THEN
      WRITE(numout,*)'Initialize humcst using reference values for 4m soil depth'
      humcste(:) = humcste_ref4m(pft_to_mtc(:))  ! values for 4m soil depth 
   ELSE
      WRITE(numout,*)'Note that humcste is initialized with values for 2m soil depth bur zmaxh=', zmaxh
      humcste(:) = humcste_ref2m(pft_to_mtc(:))  ! values for 2m soil depth
   END IF

   irrig_threshold(:) = irrig_threshold_mtc(pft_to_mtc(:))
   irrig_fulfill(:) = irrig_fulfill_mtc(pft_to_mtc(:))
   !
   ! Soil - vegetation
   !
   pref_soil_veg(:) = pref_soil_veg_mtc(pft_to_mtc(:))
   !
   !
   ! Vegetation - age classes
   !
   agec_group(:) = agec_group_mtc(pft_to_mtc(:))
   !
   ! Photosynthesis
   !
!!!qcj++ peatland
   is_peat(:)=is_peat_mtc(pft_to_mtc(:))
   is_croppeat(:)=is_croppeat_mtc(pft_to_mtc(:))
   is_shrubpeat(:)=is_shrubpeat_mtc(pft_to_mtc(:))
   is_mosspeat(:)=is_mosspeat_mtc(pft_to_mtc(:))
   is_mineralwet(:)=is_mineralwet_mtc(pft_to_mtc(:))

   is_c4(:) = is_c4_mtc(pft_to_mtc(:))
   vcmax_fix(:) = vcmax_fix_mtc(pft_to_mtc(:))
   downregulation_co2_coeff(:) = downregulation_co2_coeff_mtc(pft_to_mtc(:))
   downregulation_co2_coeff_new(:) = downregulation_co2_coeff_new_mtc(pft_to_mtc(:))
   E_KmC(:)      = E_KmC_mtc(pft_to_mtc(:))
   E_KmO(:)      = E_KmO_mtc(pft_to_mtc(:))
   E_Sco(:)      = E_Sco_mtc(pft_to_mtc(:))
   E_gamma_star(:) = E_gamma_star_mtc(pft_to_mtc(:))
   E_Vcmax(:)    = E_Vcmax_mtc(pft_to_mtc(:))
   E_Jmax(:)     = E_Jmax_mtc(pft_to_mtc(:))
   aSV(:)        = aSV_mtc(pft_to_mtc(:))
   bSV(:)        = bSV_mtc(pft_to_mtc(:))
   tphoto_min(:) = tphoto_min_mtc(pft_to_mtc(:))
   tphoto_max(:) = tphoto_max_mtc(pft_to_mtc(:))
   aSJ(:)        = aSJ_mtc(pft_to_mtc(:))
   bSJ(:)        = bSJ_mtc(pft_to_mtc(:))
   D_Vcmax(:)     = D_Vcmax_mtc(pft_to_mtc(:))
   D_Jmax(:)     = D_Jmax_mtc(pft_to_mtc(:))
   E_gm(:)       = E_gm_mtc(pft_to_mtc(:))  
   S_gm(:)       = S_gm_mtc(pft_to_mtc(:))  
   D_gm(:)       = D_gm_mtc(pft_to_mtc(:))
   E_Rd(:)       = E_Rd_mtc(pft_to_mtc(:))
   Vcmax25(:)    = Vcmax25_mtc(pft_to_mtc(:))
   arJV(:)       = arJV_mtc(pft_to_mtc(:))
   brJV(:)       = brJV_mtc(pft_to_mtc(:))
   KmC25(:)      = KmC25_mtc(pft_to_mtc(:))
   KmO25(:)      = KmO25_mtc(pft_to_mtc(:))
   Sco25(:)      = Sco25_mtc(pft_to_mtc(:)) 
   gm25(:)       = gm25_mtc(pft_to_mtc(:)) 
   gamma_star25(:)  = gamma_star25_mtc(pft_to_mtc(:))
   a1(:)         = a1_mtc(pft_to_mtc(:))
   b1(:)         = b1_mtc(pft_to_mtc(:))
   g0(:)         = g0_mtc(pft_to_mtc(:))
   h_protons(:)  = h_protons_mtc(pft_to_mtc(:))
   fpsir(:)      = fpsir_mtc(pft_to_mtc(:))
   fQ(:)         = fQ_mtc(pft_to_mtc(:))     
   fpseudo(:)    = fpseudo_mtc(pft_to_mtc(:))    
   kp(:)         = kp_mtc(pft_to_mtc(:))
   alpha(:)      = alpha_mtc(pft_to_mtc(:))
   gbs(:)        = gbs_mtc(pft_to_mtc(:))
   theta(:)      = theta_mtc(pft_to_mtc(:))        
   alpha_LL(:)   = alpha_LL_mtc(pft_to_mtc(:))
   stress_vcmax(:) = stress_vcmax_mtc(pft_to_mtc(:)) 
   stress_gs(:)    = stress_gs_mtc(pft_to_mtc(:)) 
   stress_gm(:)    = stress_gm_mtc(pft_to_mtc(:)) 
   ext_coeff(:) = ext_coeff_mtc(pft_to_mtc(:))
   ext_coeff_vegetfrac(:) = ext_coeff_vegetfrac_mtc(pft_to_mtc(:))
   !
   !! Define labels from physiologic characteristics 
   !
   leaf_tab(:) = leaf_tab_mtc(pft_to_mtc(:)) 
   pheno_model(:) = pheno_model_mtc(pft_to_mtc(:))   
   !
   is_tree(:) = .FALSE.
   DO jv = 1,nvm
      IF ( leaf_tab(jv) <= 2 ) is_tree(jv) = .TRUE.
   END DO
      !
   is_deciduous(:) = .FALSE.
   DO jv = 1,nvm
      IF ( is_tree(jv) .AND. (pheno_model(jv) /= "none") ) is_deciduous(jv) = .TRUE.
   END DO
   !
   is_evergreen(:) = .FALSE.
   DO jv = 1,nvm
      IF ( is_tree(jv) .AND. (pheno_model(jv) == "none") ) is_evergreen(jv) = .TRUE.
   END DO
   !
   is_needleleaf(:) = .FALSE.
   DO jv = 1,nvm
      IF ( leaf_tab(jv) == 2 ) is_needleleaf(jv) = .TRUE.
   END DO


    !
    ! 1. Correspondance between the PFTs values and thes MTCs values 
    !

   IF (ok_sechiba) THEN
      !
      ! Vegetation structure - sechiba
      !
      rveg_pft(:) = rveg_mtc(pft_to_mtc(:))
      !
      ! Evapotranspiration -  sechiba
      !
      rstruct_const(:) = rstruct_const_mtc(pft_to_mtc(:))
      kzero(:) = kzero_mtc(pft_to_mtc(:))
      !
      ! Water - sechiba
      !
      wmax_veg(:) = wmax_veg_mtc(pft_to_mtc(:))
      IF ( hydrol_cwrr .AND. OFF_LINE_MODE ) THEN
         throughfall_by_pft(:) = 0.
      ELSE
         throughfall_by_pft(:) = throughfall_by_mtc(pft_to_mtc(:))
      ENDIF
      !
      ! Albedo - sechiba
      !
      snowa_aged_vis(:) = snowa_aged_vis_mtc(pft_to_mtc(:))
      snowa_aged_nir(:) = snowa_aged_nir_mtc(pft_to_mtc(:))
      snowa_dec_vis(:) = snowa_dec_vis_mtc(pft_to_mtc(:)) 
      snowa_dec_nir(:) = snowa_dec_nir_mtc(pft_to_mtc(:)) 
      alb_leaf_vis(:) = alb_leaf_vis_mtc(pft_to_mtc(:))  
      alb_leaf_nir(:) = alb_leaf_nir_mtc(pft_to_mtc(:))
      !-

      !chaoyue+
      ! Permafrost - sechiba
      permafrost_veg_exists(:)= permafrost_veg_exists_mtc(pft_to_mtc(:))
      !chaoyue-
      
   ENDIF !(ok_sechiba)

    ! 1.1 For parameters used anytime

    PFT_name(:) = MTC_name(pft_to_mtc(:))
    !
    ! Vegetation structure 
    !
    veget_ori_fixed_test_1(:) = veget_ori_fixed_mtc(pft_to_mtc(:))
    llaimax(:) = llaimax_mtc(pft_to_mtc(:))
    llaimin(:) = llaimin_mtc(pft_to_mtc(:))
    height_presc(:) = height_presc_mtc(pft_to_mtc(:))
    z0_over_height(:) = z0_over_height_mtc(pft_to_mtc(:))
    ratio_z0m_z0h(:) = ratio_z0m_z0h_mtc(pft_to_mtc(:))
    type_of_lai(:) = type_of_lai_mtc(pft_to_mtc(:))
    natural(:) = natural_mtc(pft_to_mtc(:))
    !
    ! Water - sechiba
    !
    IF (zmaxh == 2.0) THEN
       IF (printlev>=2) WRITE(numout,*)'Initialize humcst using reference values for 2m soil depth'
       humcste(:) = humcste_ref2m(pft_to_mtc(:))  ! values for 2m soil depth
    ELSE IF (zmaxh == 4.0) THEN
       IF (printlev>=2) WRITE(numout,*)'Initialize humcst using reference values for 4m soil depth'
       humcste(:) = humcste_ref4m(pft_to_mtc(:))  ! values for 4m soil depth 
    ELSE
       IF (printlev>=2) WRITE(numout,*)'Note that humcste is initialized with values for 2m soil depth bur zmaxh=', zmaxh
       humcste(:) = humcste_ref2m(pft_to_mtc(:))  ! values for 2m soil depth
    END IF
    !
    ! Soil - vegetation
    !
    pref_soil_veg(:) = pref_soil_veg_mtc(pft_to_mtc(:))
    !
    ! Photosynthesis
    !
!!!qcj++ peatland
    is_peat(:)=is_peat_mtc(pft_to_mtc(:))
    is_croppeat(:)=is_croppeat_mtc(pft_to_mtc(:))
    is_shrubpeat(:)=is_shrubpeat_mtc(pft_to_mtc(:))
    is_mosspeat(:)=is_mosspeat_mtc(pft_to_mtc(:))
    is_mineralwet(:)=is_mineralwet_mtc(pft_to_mtc(:))

    is_c4(:) = is_c4_mtc(pft_to_mtc(:))
    vcmax_fix(:) = vcmax_fix_mtc(pft_to_mtc(:))
    downregulation_co2_coeff(:) = downregulation_co2_coeff_mtc(pft_to_mtc(:))
   downregulation_co2_coeff_new(:) = downregulation_co2_coeff_new_mtc(pft_to_mtc(:))
    E_KmC(:)      = E_KmC_mtc(pft_to_mtc(:))
    E_KmO(:)      = E_KmO_mtc(pft_to_mtc(:))
    E_Sco(:)      = E_Sco_mtc(pft_to_mtc(:))
    E_gamma_star(:) = E_gamma_star_mtc(pft_to_mtc(:))
    E_Vcmax(:)    = E_Vcmax_mtc(pft_to_mtc(:))
    E_Jmax(:)     = E_Jmax_mtc(pft_to_mtc(:))
    aSV(:)        = aSV_mtc(pft_to_mtc(:))
    bSV(:)        = bSV_mtc(pft_to_mtc(:))
    tphoto_min(:) = tphoto_min_mtc(pft_to_mtc(:))
    tphoto_max(:) = tphoto_max_mtc(pft_to_mtc(:))
    aSJ(:)        = aSJ_mtc(pft_to_mtc(:))
    bSJ(:)        = bSJ_mtc(pft_to_mtc(:))
    D_Vcmax(:)     = D_Vcmax_mtc(pft_to_mtc(:))
    D_Jmax(:)     = D_Jmax_mtc(pft_to_mtc(:))
    E_gm(:)       = E_gm_mtc(pft_to_mtc(:)) 
    S_gm(:)       = S_gm_mtc(pft_to_mtc(:)) 
    D_gm(:)       = D_gm_mtc(pft_to_mtc(:)) 
    E_Rd(:)       = E_Rd_mtc(pft_to_mtc(:))
    Vcmax25(:)    = Vcmax25_mtc(pft_to_mtc(:))
    arJV(:)       = arJV_mtc(pft_to_mtc(:))
    brJV(:)       = brJV_mtc(pft_to_mtc(:))
    KmC25(:)      = KmC25_mtc(pft_to_mtc(:))
    KmO25(:)      = KmO25_mtc(pft_to_mtc(:))
    Sco25(:)      = Sco25_mtc(pft_to_mtc(:))
    gm25(:)       = gm25_mtc(pft_to_mtc(:)) 
    gamma_star25(:)  = gamma_star25_mtc(pft_to_mtc(:))
    a1(:)         = a1_mtc(pft_to_mtc(:))
    b1(:)         = b1_mtc(pft_to_mtc(:))
    g0(:)         = g0_mtc(pft_to_mtc(:))
    h_protons(:)  = h_protons_mtc(pft_to_mtc(:))
    fpsir(:)      = fpsir_mtc(pft_to_mtc(:))
    fQ(:)         = fQ_mtc(pft_to_mtc(:))     
    fpseudo(:)    = fpseudo_mtc(pft_to_mtc(:))    
    kp(:)         = kp_mtc(pft_to_mtc(:))
    alpha(:)      = alpha_mtc(pft_to_mtc(:))
    gbs(:)        = gbs_mtc(pft_to_mtc(:))
    theta(:)      = theta_mtc(pft_to_mtc(:))        
    alpha_LL(:)   = alpha_LL_mtc(pft_to_mtc(:))
    stress_vcmax(:) = stress_vcmax_mtc(pft_to_mtc(:))
    stress_gs(:)    = stress_gs_mtc(pft_to_mtc(:))
    stress_gm(:)    = stress_gm_mtc(pft_to_mtc(:))
    ext_coeff(:) = ext_coeff_mtc(pft_to_mtc(:))
    ext_coeff_vegetfrac(:) = ext_coeff_vegetfrac_mtc(pft_to_mtc(:))
    !
    !! Define labels from physiologic characteristics 
    !
    leaf_tab(:) = leaf_tab_mtc(pft_to_mtc(:)) 
    pheno_model(:) = pheno_model_mtc(pft_to_mtc(:))   
    !
    is_tree(:) = .FALSE.
    DO jv = 1,nvm
       IF ( leaf_tab(jv) <= 2 ) is_tree(jv) = .TRUE.
    END DO
    !
    is_deciduous(:) = .FALSE.
    DO jv = 1,nvm
       IF ( is_tree(jv) .AND. (pheno_model(jv) /= "none") ) is_deciduous(jv) = .TRUE.
    END DO
    !
    is_evergreen(:) = .FALSE.
    DO jv = 1,nvm
       IF ( is_tree(jv) .AND. (pheno_model(jv) == "none") ) is_evergreen(jv) = .TRUE.
    END DO
    !
    is_needleleaf(:) = .FALSE.
    DO jv = 1,nvm
       IF ( leaf_tab(jv) == 2 ) is_needleleaf(jv) = .TRUE.
    END DO


    ! 1.2 For sechiba parameters

   IF (ok_stomate) THEN
      bm_sapl(:,:,:) = zero
      maxdia(:) = undef
      migrate(:) = undef
      cn_sapl(:) = undef
      leaf_timecst(:) = undef
      lai_initmin(:) = undef
      !
      ! Vegetation structure - stomate
      !
      sla(:) = sla_mtc(pft_to_mtc(:))
      availability_fact(:) = availability_fact_mtc(pft_to_mtc(:))
      !
      ! Allocation - stomate
      !
      R0(:) = R0_mtc(pft_to_mtc(:)) 
      S0(:) = S0_mtc(pft_to_mtc(:)) 
      !
      !pss+:Wetland CH4 methane
      !
      rdepth_v(:) = rdepth_v_mtc(pft_to_mtc(:))
      sdepth_v(:) = sdepth_v_mtc(pft_to_mtc(:))
      tveg_v(:) = tveg_v_mtc(pft_to_mtc(:))
      !pss-

      !
      ! Respiration - stomate
      !
      frac_growthresp(:) = frac_growthresp_mtc(pft_to_mtc(:))  
      maint_resp_slope_c(:) = maint_resp_slope_c_mtc(pft_to_mtc(:))               
      maint_resp_slope_b(:) = maint_resp_slope_b_mtc(pft_to_mtc(:))
      maint_resp_slope_a(:) = maint_resp_slope_a_mtc(pft_to_mtc(:))
      cm_zero_leaf(:) = cm_zero_leaf_mtc(pft_to_mtc(:))
      cm_zero_sapabove(:) = cm_zero_sapabove_mtc(pft_to_mtc(:))
      cm_zero_sapbelow(:) = cm_zero_sapbelow_mtc(pft_to_mtc(:)) 
      cm_zero_heartabove(:) = cm_zero_heartabove_mtc(pft_to_mtc(:)) 
      cm_zero_heartbelow(:) = cm_zero_heartbelow_mtc(pft_to_mtc(:))
      cm_zero_root(:) = cm_zero_root_mtc(pft_to_mtc(:))
      cm_zero_fruit(:) = cm_zero_fruit_mtc(pft_to_mtc(:))
      cm_zero_carbres(:) = cm_zero_carbres_mtc(pft_to_mtc(:))
      !
      ! Fire - stomate
      !
      flam(:) = flam_mtc(pft_to_mtc(:))
      resist(:) = resist_mtc(pft_to_mtc(:))
      !spitfire
      dens_fuel(:) = dens_fuel_mtc(pft_to_mtc(:))
      f_sh(:) = f_sh_mtc(pft_to_mtc(:))
      crown_length(:) = crown_length_mtc(pft_to_mtc(:))
      BTpar1(:) = BTpar1_mtc(pft_to_mtc(:))
      BTpar2(:) = BTpar2_mtc(pft_to_mtc(:))
      r_ck(:) = r_ck_mtc(pft_to_mtc(:))
      p_ck(:) = p_ck_mtc(pft_to_mtc(:))
      ef_CO2(:) = ef_CO2_mtc(pft_to_mtc(:))
      ef_CO(:) = ef_CO_mtc(pft_to_mtc(:))
      ef_CH4(:) = ef_CH4_mtc(pft_to_mtc(:))
      ef_VOC(:) = ef_VOC_mtc(pft_to_mtc(:))
      ef_TPM(:) = ef_TPM_mtc(pft_to_mtc(:))
      ef_NOx(:) = ef_NOx_mtc(pft_to_mtc(:))
      me(:) = me_mtc(pft_to_mtc(:))
      fire_max_cf_100hr(:) = fire_max_cf_100hr_mtc(pft_to_mtc(:))
      fire_max_cf_1000hr(:) = fire_max_cf_1000hr_mtc(pft_to_mtc(:))
      !endspit
      !
      ! grassland management 
      !
      !gmjc
      is_grassland_manag(:) = is_grassland_manag_mtc(pft_to_mtc(:))
      is_grassland_cut(:) = is_grassland_cut_mtc(pft_to_mtc(:))
      is_grassland_grazed(:) = is_grassland_grazed_mtc(pft_to_mtc(:))
      management_intensity(:) = management_intensity_mtc(pft_to_mtc(:))
      management_start(:) = management_start_mtc(pft_to_mtc(:))
      deposition_start(:) = deposition_start_mtc(pft_to_mtc(:))
      nb_year_management(:) = nb_year_management_mtc(pft_to_mtc(:))
      sla_min(:) = sla_min_mtc(pft_to_mtc(:))
      sla_max(:) = sla_max_mtc(pft_to_mtc(:))
      !end gmjc

      !gluc plus bioenergy
      is_bioe1(:) = is_bioe1_mtc(pft_to_mtc(:))

      !
      ! Flux - LUC
      !
      coeff_lcchange_1(:) = coeff_lcchange_1_mtc(pft_to_mtc(:))
      coeff_lcchange_10(:) = coeff_lcchange_10_mtc(pft_to_mtc(:))
      coeff_lcchange_100(:) = coeff_lcchange_100_mtc(pft_to_mtc(:))

      coeff_indwood_1(:) = coeff_indwood_1_mtc(pft_to_mtc(:))
      coeff_indwood_10(:) = coeff_indwood_10_mtc(pft_to_mtc(:))
      coeff_indwood_100(:) = coeff_indwood_100_mtc(pft_to_mtc(:))
      !
      ! Phenology
      !
      !
      ! 1. Stomate
      !
      lai_max_to_happy(:) = lai_max_to_happy_mtc(pft_to_mtc(:))  
      lai_max(:) = lai_max_mtc(pft_to_mtc(:))
      pheno_type(:) = pheno_type_mtc(pft_to_mtc(:))
      !
      ! 2. Leaf Onset
      !
      pheno_gdd_crit_c(:) = pheno_gdd_crit_c_mtc(pft_to_mtc(:))
      pheno_gdd_crit_b(:) = pheno_gdd_crit_b_mtc(pft_to_mtc(:))         
      pheno_gdd_crit_a(:) = pheno_gdd_crit_a_mtc(pft_to_mtc(:))
      pheno_moigdd_t_crit(:) = pheno_moigdd_t_crit_mtc(pft_to_mtc(:))
      ngd_crit(:) =  ngd_crit_mtc(pft_to_mtc(:))
      ncdgdd_temp(:) = ncdgdd_temp_mtc(pft_to_mtc(:)) 
      hum_frac(:) = hum_frac_mtc(pft_to_mtc(:))
      hum_min_time(:) = hum_min_time_mtc(pft_to_mtc(:))
      tau_sap(:) = tau_sap_mtc(pft_to_mtc(:))
      tau_leafinit(:) = tau_leafinit_mtc(pft_to_mtc(:))  
      tau_fruit(:) = tau_fruit_mtc(pft_to_mtc(:))
      ecureuil(:) = ecureuil_mtc(pft_to_mtc(:))
      alloc_min(:) = alloc_min_mtc(pft_to_mtc(:))
      alloc_max(:) = alloc_max_mtc(pft_to_mtc(:))
      demi_alloc(:) = demi_alloc_mtc(pft_to_mtc(:))
      leaflife_tab(:) = leaflife_mtc(pft_to_mtc(:))
      !
      ! 3. Senescence
      !
      leaffall(:) = leaffall_mtc(pft_to_mtc(:))
      leafagecrit(:) = leafagecrit_mtc(pft_to_mtc(:))
      senescence_type(:) = senescence_type_mtc(pft_to_mtc(:)) 
      senescence_hum(:) = senescence_hum_mtc(pft_to_mtc(:)) 
      nosenescence_hum(:) = nosenescence_hum_mtc(pft_to_mtc(:)) 
      max_turnover_time(:) = max_turnover_time_mtc(pft_to_mtc(:))
      min_turnover_time(:) = min_turnover_time_mtc(pft_to_mtc(:))
      min_leaf_age_for_senescence(:) = min_leaf_age_for_senescence_mtc(pft_to_mtc(:))
      senescence_temp_c(:) = senescence_temp_c_mtc(pft_to_mtc(:))
      senescence_temp_b(:) = senescence_temp_b_mtc(pft_to_mtc(:))
      senescence_temp_a(:) = senescence_temp_a_mtc(pft_to_mtc(:))
      gdd_senescence(:) = gdd_senescence_mtc(pft_to_mtc(:))
      !
      ! DGVM
      !
      residence_time(:) = residence_time_mtc(pft_to_mtc(:))
      tmin_crit(:) = tmin_crit_mtc(pft_to_mtc(:))
      tcm_crit(:) = tcm_crit_mtc(pft_to_mtc(:))
      !qcj++ peatland
      wtpwet_crit(:) = wtpwet_crit_mtc(pft_to_mtc(:))
      wtpdry_crit(:) = wtpdry_crit_mtc(pft_to_mtc(:))
      wtp_crit(:) = wtp_crit_mtc(pft_to_mtc(:))
      wt_mortality(:) = wt_mortality_mtc(pft_to_mtc(:))

      !-

      !!!!! crop parameters

      ! STICS:: main LAIdev
      ok_LAIdev(:) = ok_LAIdev_mtc(pft_to_mtc(:))
      ! STICS:: 
      SP_codeplante(:) = SP_codeplante_mtc(pft_to_mtc(:))
      SP_stade0(:) = SP_stade0_mtc(pft_to_mtc(:))
      SP_iplt0(:) = SP_iplt0_mtc(pft_to_mtc(:))
      SP_nbox(:) = SP_nbox_mtc(pft_to_mtc(:))
      SP_iwater(:) = SP_iwater_mtc(pft_to_mtc(:))
      SP_codesimul(:) = SP_codesimul_mtc(pft_to_mtc(:))
      SP_codelaitr(:) = SP_codelaitr_mtc(pft_to_mtc(:))
      SP_slamax(:) = SP_slamax_mtc(pft_to_mtc(:))
      SP_slamin(:) = SP_slamin_mtc(pft_to_mtc(:))
      SP_codeperenne(:) = SP_codeperenne_mtc(pft_to_mtc(:))

      SP_codcueille(:) = SP_codcueille_mtc(pft_to_mtc(:))
      SP_codegdh(:) = SP_codegdh_mtc(pft_to_mtc(:))
      SP_codetemp(:) = SP_codetemp_mtc(pft_to_mtc(:))
      SP_coderetflo(:) = SP_coderetflo_mtc(pft_to_mtc(:))
      SP_codeinnact(:) = SP_codeinnact_mtc(pft_to_mtc(:))
      SP_codeh2oact(:) = SP_codeh2oact_mtc(pft_to_mtc(:))
      SP_stressdev(:) = SP_stressdev_mtc(pft_to_mtc(:))
      SP_innlai(:) = SP_innlai_mtc(pft_to_mtc(:))
      SP_innsenes(:) = SP_innsenes_mtc(pft_to_mtc(:))
      SP_codebfroid(:) = SP_codebfroid_mtc(pft_to_mtc(:))

      SP_codephot(:) = SP_codephot_mtc(pft_to_mtc(:))
      SP_codedormance(:) = SP_codedormance_mtc(pft_to_mtc(:))
      SP_codefauche(:) = SP_codefauche_mtc(pft_to_mtc(:))
      SP_codetempfauche(:) = SP_codetempfauche_mtc(pft_to_mtc(:))
      SP_codlainet(:) = SP_codlainet_mtc(pft_to_mtc(:))
      SP_codeindetermin(:) = SP_codeindetermin_mtc(pft_to_mtc(:))
      SP_codeinitprec(:) = SP_codeinitprec_mtc(pft_to_mtc(:))
      SP_culturean(:) = SP_culturean_mtc(pft_to_mtc(:))

      SP_jvc(:) = SP_jvc_mtc(pft_to_mtc(:))
      SP_tfroid(:) = SP_tfroid_mtc(pft_to_mtc(:))
      SP_ampfroid(:) = SP_ampfroid_mtc(pft_to_mtc(:))
      SP_jvcmini(:) = SP_jvcmini_mtc(pft_to_mtc(:))
      SP_tgmin(:) = SP_tgmin_mtc(pft_to_mtc(:))
      SP_stpltger(:) = SP_stpltger_mtc(pft_to_mtc(:))
      SP_profsem(:) = SP_profsem_mtc(pft_to_mtc(:))
      SP_propjgermin(:) = SP_propjgermin_mtc(pft_to_mtc(:))

      SP_tdmax(:) = SP_tdmax_mtc(pft_to_mtc(:))
      SP_nbjgerlim(:) = SP_nbjgerlim_mtc(pft_to_mtc(:))
      SP_densitesem(:) = SP_densitesem_mtc(pft_to_mtc(:))
      SP_vigueurbat(:) = SP_vigueurbat_mtc(pft_to_mtc(:))
      SP_codepluiepoquet(:) = SP_codepluiepoquet_mtc(pft_to_mtc(:))
      SP_codehypo(:) = SP_codehypo_mtc(pft_to_mtc(:))
      SP_elmax(:) = SP_elmax_mtc(pft_to_mtc(:))
      SP_belong(:) = SP_belong_mtc(pft_to_mtc(:))

      SP_celong(:) = SP_celong_mtc(pft_to_mtc(:))
      SP_nlevlim1(:) = SP_nlevlim1_mtc(pft_to_mtc(:))
      SP_nlevlim2(:) = SP_nlevlim2_mtc(pft_to_mtc(:))
      SP_codrecolte(:) = SP_codrecolte_mtc(pft_to_mtc(:))
      SP_variete(:) = SP_variete_mtc(pft_to_mtc(:))
      SP_codegermin(:) = SP_codegermin_mtc(pft_to_mtc(:))

      S_codeulaivernal(:) = S_codeulaivernal_mtc(pft_to_mtc(:))
      SP_swfacmin(:) = SP_swfacmin_mtc(pft_to_mtc(:))
      SP_neffmax(:) = SP_neffmax_mtc(pft_to_mtc(:))
      SP_nsatrat(:) = SP_nsatrat_mtc(pft_to_mtc(:))


      ! STICS:: LAI CALCULATION
      SP_laiplantule(:) = SP_laiplantule_mtc(pft_to_mtc(:))
      SP_vlaimax(:) = SP_vlaimax_mtc(pft_to_mtc(:))
      SP_stlevamf(:) = SP_stlevamf_mtc(pft_to_mtc(:))
      SP_stdrpmat(:) = SP_stdrpmat_mtc(pft_to_mtc(:))
      SP_stamflax(:) = SP_stamflax_mtc(pft_to_mtc(:))
      SP_udlaimax(:) = SP_udlaimax_mtc(pft_to_mtc(:))
      SP_laicomp(:) = SP_laicomp_mtc(pft_to_mtc(:))
      SP_adens(:) = SP_adens_mtc(pft_to_mtc(:))
      SP_bdens(:) = SP_bdens_mtc(pft_to_mtc(:))

      SP_tcxstop(:) = SP_tcxstop_mtc(pft_to_mtc(:))
      SP_tcmax(:) = SP_tcmax_mtc(pft_to_mtc(:))
      SP_tcmin(:) = SP_tcmin_mtc(pft_to_mtc(:))
      SP_dlaimax(:) = SP_dlaimax_mtc(pft_to_mtc(:))
      SP_dlaimin(:) = SP_dlaimin_mtc(pft_to_mtc(:))
      SP_pentlaimax(:) = SP_pentlaimax_mtc(pft_to_mtc(:))
      SP_tigefeuil(:) = SP_tigefeuil_mtc(pft_to_mtc(:))
     
      SP_stlaxsen(:) = SP_stlaxsen_mtc(pft_to_mtc(:))
      SP_stsenlan(:) = SP_stsenlan_mtc(pft_to_mtc(:))
      SP_stlevdrp(:) = SP_stlevdrp_mtc(pft_to_mtc(:))
      SP_stflodrp(:) = SP_stflodrp_mtc(pft_to_mtc(:))
      SP_stdrpdes(:) = SP_stdrpdes_mtc(pft_to_mtc(:))
     
      SP_phyllotherme(:) = SP_phyllotherme_mtc(pft_to_mtc(:))
      SP_lai0(:) = SP_lai0_mtc(pft_to_mtc(:))
      SP_tustressmin(:) = SP_tustressmin_mtc(pft_to_mtc(:))

 
      ! STICS:: LAI SENESCENCE
      SP_nbfgellev(:) = SP_nbfgellev_mtc(pft_to_mtc(:))
      SP_ratiodurvieI(:) = SP_ratiodurvieI_mtc(pft_to_mtc(:))
      SP_durvieF(:) = SP_durvieF_mtc(pft_to_mtc(:))
      SP_ratiosen(:) = SP_ratiosen_mtc(pft_to_mtc(:))
      SP_tdmin(:) = SP_tdmin_mtc(pft_to_mtc(:))
      
      ! STICS:: F_humerac
 
      SP_sensrsec(:) = SP_sensrsec_mtc(pft_to_mtc(:))

      ! STICS:: gel

      SP_codgellev(:) = SP_codgellev_mtc(pft_to_mtc(:))
      SP_tletale(:) = SP_tletale_mtc(pft_to_mtc(:))
      SP_tdebgel(:) = SP_tdebgel_mtc(pft_to_mtc(:))
      SP_tgellev10(:) = SP_tgellev10_mtc(pft_to_mtc(:))
      SP_tgellev90(:) = SP_tgellev90_mtc(pft_to_mtc(:))

      SP_tgeljuv10(:) = SP_tgeljuv10_mtc(pft_to_mtc(:))
      SP_tgeljuv90(:) = SP_tgeljuv90_mtc(pft_to_mtc(:))
      SP_tgelveg10(:) = SP_tgelveg10_mtc(pft_to_mtc(:))
      SP_tgelveg90(:) = SP_tgelveg90_mtc(pft_to_mtc(:))


      ! STICS:: PHOTOPERIOD

      SP_sensiphot(:) = SP_sensiphot_mtc(pft_to_mtc(:))
      SP_phosat(:) = SP_phosat_mtc(pft_to_mtc(:))
      SP_phobase(:) = SP_phobase_mtc(pft_to_mtc(:))
 
      ! STICS:: CARBON ALLOCATION
      
      SP_stoprac(:) = SP_stoprac_mtc(pft_to_mtc(:))
      SP_zracplantule(:) = SP_zracplantule_mtc(pft_to_mtc(:))
      SP_codtrophrac(:) = SP_codtrophrac_mtc(pft_to_mtc(:))
      SP_repracpermax(:) = SP_repracpermax_mtc(pft_to_mtc(:))
      SP_repracpermin(:) = SP_repracpermin_mtc(pft_to_mtc(:))
      SP_krepracperm(:) = SP_krepracperm_mtc(pft_to_mtc(:))
      SP_repracseumax(:) = SP_repracseumax_mtc(pft_to_mtc(:))
      SP_repracseumin(:) = SP_repracseumin_mtc(pft_to_mtc(:))
      SP_krepracseu(:) = SP_krepracseu_mtc(pft_to_mtc(:))
      SP_codetemprac(:) = SP_codetemprac_mtc(pft_to_mtc(:))
      SP_codedyntalle(:) = SP_codedyntalle_mtc(pft_to_mtc(:))
      SP_nbjgrain(:) = SP_nbjgrain_mtc(pft_to_mtc(:))
      SP_maxgs(:) = SP_maxgs_mtc(pft_to_mtc(:))
      SP_codgelflo(:) = SP_codgelflo_mtc(pft_to_mtc(:))
      SP_tgelflo10(:) = SP_tgelflo10_mtc(pft_to_mtc(:))
      SP_tgelflo90(:) = SP_tgelflo90_mtc(pft_to_mtc(:))
      SP_cgrain(:) = SP_cgrain_mtc(pft_to_mtc(:))
      SP_cgrainv0(:) = SP_cgrainv0_mtc(pft_to_mtc(:))
      SP_nbgrmax(:) = SP_nbgrmax_mtc(pft_to_mtc(:))
      SP_nbgrmin(:) = SP_nbgrmin_mtc(pft_to_mtc(:))
      SP_codazofruit(:) = SP_codazofruit_mtc(pft_to_mtc(:))
      SP_codeir(:) = SP_codeir_mtc(pft_to_mtc(:))
      SP_vitircarb(:) = SP_vitircarb_mtc(pft_to_mtc(:))
      SP_irmax(:) = SP_irmax_mtc(pft_to_mtc(:))
      SP_vitircarbT(:) = SP_vitircarbT_mtc(pft_to_mtc(:))
      SP_codetremp(:) = SP_codetremp_mtc(pft_to_mtc(:))
      SP_tminremp(:) = SP_tminremp_mtc(pft_to_mtc(:))
      SP_tmaxremp(:) = SP_tmaxremp_mtc(pft_to_mtc(:))
      SP_pgrainmaxi(:) = SP_pgrainmaxi_mtc(pft_to_mtc(:))
      
      !! SPECIFIC FOR DYNAMIC INN STRATEGY
      
      SP_DY_INN(:) = SP_DY_INN_mtc(pft_to_mtc(:))
      SP_avenfert(:) = SP_avenfert_mtc(pft_to_mtc(:))
      
      !!!!! end crop parameters
   ENDIF !(ok_stomate)

    ! 1.3 For BVOC parameters

    IF (ok_bvoc) THEN
       !
       ! Biogenic Volatile Organic Compounds
       !
       em_factor_isoprene(:) = em_factor_isoprene_mtc(pft_to_mtc(:))
       em_factor_monoterpene(:) = em_factor_monoterpene_mtc(pft_to_mtc(:))
       LDF_mono = LDF_mono_mtc 
       LDF_sesq = LDF_sesq_mtc 
       LDF_meth = LDF_meth_mtc 
       LDF_acet = LDF_acet_mtc 

       em_factor_apinene(:) = em_factor_apinene_mtc(pft_to_mtc(:))
       em_factor_bpinene(:) = em_factor_bpinene_mtc(pft_to_mtc(:))
       em_factor_limonene(:) = em_factor_limonene_mtc(pft_to_mtc(:))
       em_factor_myrcene(:) = em_factor_myrcene_mtc(pft_to_mtc(:))
       em_factor_sabinene(:) = em_factor_sabinene_mtc(pft_to_mtc(:))
       em_factor_camphene(:) = em_factor_camphene_mtc(pft_to_mtc(:))
       em_factor_3carene(:) = em_factor_3carene_mtc(pft_to_mtc(:))
       em_factor_tbocimene(:) = em_factor_tbocimene_mtc(pft_to_mtc(:))
       em_factor_othermonot(:) = em_factor_othermonot_mtc(pft_to_mtc(:))
       em_factor_sesquiterp(:) = em_factor_sesquiterp_mtc(pft_to_mtc(:))

       beta_mono = beta_mono_mtc
       beta_sesq = beta_sesq_mtc
       beta_meth = beta_meth_mtc
       beta_acet = beta_acet_mtc
       beta_oxyVOC = beta_oxyVOC_mtc

       em_factor_ORVOC(:) = em_factor_ORVOC_mtc(pft_to_mtc(:)) 
       em_factor_OVOC(:) = em_factor_OVOC_mtc(pft_to_mtc(:))
       em_factor_MBO(:) = em_factor_MBO_mtc(pft_to_mtc(:))
       em_factor_methanol(:) = em_factor_methanol_mtc(pft_to_mtc(:))
       em_factor_acetone(:) = em_factor_acetone_mtc(pft_to_mtc(:)) 
       em_factor_acetal(:) = em_factor_acetal_mtc(pft_to_mtc(:))
       em_factor_formal(:) = em_factor_formal_mtc(pft_to_mtc(:))
       em_factor_acetic(:) = em_factor_acetic_mtc(pft_to_mtc(:))
       em_factor_formic(:) = em_factor_formic_mtc(pft_to_mtc(:))
       em_factor_no_wet(:) = em_factor_no_wet_mtc(pft_to_mtc(:))
       em_factor_no_dry(:) = em_factor_no_dry_mtc(pft_to_mtc(:))
       Larch(:) = Larch_mtc(pft_to_mtc(:)) 
       !-
    ENDIF !(ok_bvoc)

    ! 1.4 For stomate parameters

    IF (ok_stomate) THEN
       !
       ! Vegetation structure - stomate
       !
       sla(:) = sla_mtc(pft_to_mtc(:))
       availability_fact(:) = availability_fact_mtc(pft_to_mtc(:))
       !
       ! Allocation - stomate
       !
       R0(:) = R0_mtc(pft_to_mtc(:)) 
       S0(:) = S0_mtc(pft_to_mtc(:)) 
       !
       ! Respiration - stomate
       !
       frac_growthresp(:) = frac_growthresp_mtc(pft_to_mtc(:))  
       maint_resp_slope_c(:) = maint_resp_slope_c_mtc(pft_to_mtc(:))               
       maint_resp_slope_b(:) = maint_resp_slope_b_mtc(pft_to_mtc(:))
       maint_resp_slope_a(:) = maint_resp_slope_a_mtc(pft_to_mtc(:))
       cm_zero_leaf(:) = cm_zero_leaf_mtc(pft_to_mtc(:))
       cm_zero_sapabove(:) = cm_zero_sapabove_mtc(pft_to_mtc(:))
       cm_zero_sapbelow(:) = cm_zero_sapbelow_mtc(pft_to_mtc(:)) 
       cm_zero_heartabove(:) = cm_zero_heartabove_mtc(pft_to_mtc(:)) 
       cm_zero_heartbelow(:) = cm_zero_heartbelow_mtc(pft_to_mtc(:))
       cm_zero_root(:) = cm_zero_root_mtc(pft_to_mtc(:))
       cm_zero_fruit(:) = cm_zero_fruit_mtc(pft_to_mtc(:))
       cm_zero_carbres(:) = cm_zero_carbres_mtc(pft_to_mtc(:))
       !
       ! Fire - stomate
       !
       flam(:) = flam_mtc(pft_to_mtc(:))
       resist(:) = resist_mtc(pft_to_mtc(:))
       !
       ! Flux - LUC
       !
       coeff_lcchange_1(:) = coeff_lcchange_1_mtc(pft_to_mtc(:))
       coeff_lcchange_10(:) = coeff_lcchange_10_mtc(pft_to_mtc(:))
       coeff_lcchange_100(:) = coeff_lcchange_100_mtc(pft_to_mtc(:))

       coeff_indwood_1(:) = coeff_indwood_1_mtc(pft_to_mtc(:))
       coeff_indwood_10(:) = coeff_indwood_10_mtc(pft_to_mtc(:))
       coeff_indwood_100(:) = coeff_indwood_100_mtc(pft_to_mtc(:))
       !
       ! Phenology
       !
       !
       ! 1. Stomate
       !
       lai_max_to_happy(:) = lai_max_to_happy_mtc(pft_to_mtc(:))  
       lai_max(:) = lai_max_mtc(pft_to_mtc(:))
       pheno_type(:) = pheno_type_mtc(pft_to_mtc(:))
       !
       ! 2. Leaf Onset
       !
       pheno_gdd_crit_c(:) = pheno_gdd_crit_c_mtc(pft_to_mtc(:))
       pheno_gdd_crit_b(:) = pheno_gdd_crit_b_mtc(pft_to_mtc(:))         
       pheno_gdd_crit_a(:) = pheno_gdd_crit_a_mtc(pft_to_mtc(:))
       pheno_moigdd_t_crit(:) = pheno_moigdd_t_crit_mtc(pft_to_mtc(:))
       ngd_crit(:) =  ngd_crit_mtc(pft_to_mtc(:))
       ncdgdd_temp(:) = ncdgdd_temp_mtc(pft_to_mtc(:)) 
       hum_frac(:) = hum_frac_mtc(pft_to_mtc(:))
       hum_min_time(:) = hum_min_time_mtc(pft_to_mtc(:))
       tau_sap(:) = tau_sap_mtc(pft_to_mtc(:))
       tau_leafinit(:) = tau_leafinit_mtc(pft_to_mtc(:))  
       tau_fruit(:) = tau_fruit_mtc(pft_to_mtc(:))
       ecureuil(:) = ecureuil_mtc(pft_to_mtc(:))
       alloc_min(:) = alloc_min_mtc(pft_to_mtc(:))
       alloc_max(:) = alloc_max_mtc(pft_to_mtc(:))
       demi_alloc(:) = demi_alloc_mtc(pft_to_mtc(:))
       leaflife_tab(:) = leaflife_mtc(pft_to_mtc(:))
       !
       ! 3. Senescence
       !
       leaffall(:) = leaffall_mtc(pft_to_mtc(:))
       leafagecrit(:) = leafagecrit_mtc(pft_to_mtc(:))
       senescence_type(:) = senescence_type_mtc(pft_to_mtc(:)) 
       senescence_hum(:) = senescence_hum_mtc(pft_to_mtc(:)) 
       nosenescence_hum(:) = nosenescence_hum_mtc(pft_to_mtc(:)) 
       max_turnover_time(:) = max_turnover_time_mtc(pft_to_mtc(:))
       min_turnover_time(:) = min_turnover_time_mtc(pft_to_mtc(:))
       min_leaf_age_for_senescence(:) = min_leaf_age_for_senescence_mtc(pft_to_mtc(:))
       senescence_temp_c(:) = senescence_temp_c_mtc(pft_to_mtc(:))
       senescence_temp_b(:) = senescence_temp_b_mtc(pft_to_mtc(:))
       senescence_temp_a(:) = senescence_temp_a_mtc(pft_to_mtc(:))
       gdd_senescence(:) = gdd_senescence_mtc(pft_to_mtc(:))
       !
       ! DGVM
       !
       residence_time(:) = residence_time_mtc(pft_to_mtc(:))
       tmin_crit(:) = tmin_crit_mtc(pft_to_mtc(:))
       tcm_crit(:) = tcm_crit_mtc(pft_to_mtc(:))
      !qcj++ peatland
       wtpwet_crit(:) = wtpwet_crit_mtc(pft_to_mtc(:))
       wtpdry_crit(:) = wtpdry_crit_mtc(pft_to_mtc(:))
       wtp_crit(:) = wtp_crit_mtc(pft_to_mtc(:))
       wt_mortality(:) = wt_mortality_mtc(pft_to_mtc(:))
       !-
    ENDIF !(ok_stomate)

  END SUBROUTINE pft_parameters_init


!! ================================================================================================================================
!! SUBROUTINE   : pft_parameters_alloc
!!
!>\BRIEF         This subroutine allocates memory needed for the PFT parameters 
!! in function  of the flags activated.  
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE pft_parameters_alloc()

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.1 Input variables 

    !! 0.4 Local variables

    LOGICAL :: l_error                             !! Diagnostic boolean for error allocation (true/false) 
    INTEGER :: ier                                !! Return value for memory allocation (0-N, unitless)
    CHARACTER(LEN=100) :: str !! Temporary variable

    !_ ================================================================================================================================

    !
    ! 1. Parameters used anytime
    !

    l_error = .FALSE.

    ALLOCATE(pft_to_mtc(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for pft_to_mtc. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(PFT_name(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for PFT_name. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(height_presc(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for height_presc. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(z0_over_height(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for z0_over_height. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(ratio_z0m_z0h(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for ratio_z0m_z0h. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF
! dgvmjc
    ALLOCATE(pasture(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for pasture. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF
! end dgvmjc
!!! crop irrig
   ALLOCATE(irrig_threshold(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for irrig_threshold. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(irrig_fulfill(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for irrig_fulfill. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF
!!! end crop irrig, xuhui
    ALLOCATE(natural(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for natural. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(is_tree(nvm),stat=ier) 
       l_error = l_error .OR. (ier /= 0) 
       IF (l_error) THEN 
          WRITE(numout,*) ' Memory allocation error for is_tree. We stop. We need nvm words = ',nvm 
       STOP 'pft_parameters_alloc' 
    END IF 

!!!qcj++ peatland
    ALLOCATE(is_peat(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for is_peat. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(is_croppeat(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for is_croppeat. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(is_shrubpeat(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for is_shrubpeat. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(is_mosspeat(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for is_mosspeat. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(is_mineralwet(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for is_mineralwet. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(is_c4(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for is_c4. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(vcmax_fix(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for vcmax_fix. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(humcste(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for humcste. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(downregulation_co2_coeff(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for downregulation_co2_coeff. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(downregulation_co2_coeff_new(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for downregulation_co2_coeff_new. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(E_KmC(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for E_KmC. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(E_KmO(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for E_KmO. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(E_Sco(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for E_Sco. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(E_gamma_star(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for E_gamma_star. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(E_vcmax(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for E_Vcmax. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(E_Jmax(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for E_Jmax. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(aSV(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for aSV. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(bSV(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for bSV. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(tphoto_min(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for tphoto_min. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(tphoto_max(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for tphoto_max. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(aSJ(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for aSJ. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(bSJ(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for bSJ. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(D_Vcmax(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for D_Vcmax. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(D_Jmax(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for D_Jmax. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(E_gm(nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for E_gm. We stop. We need nvm words = ',nvm 
       STOP 'pft_parameters_alloc' 
    END IF
    
    ALLOCATE(S_gm(nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for S_gm. We stop. We need nvm words = ',nvm 
       STOP 'pft_parameters_alloc' 
    END IF
    
    ALLOCATE(D_gm(nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for D_gm. We stop. We need nvm words = ',nvm 
       STOP 'pft_parameters_alloc' 
    END IF
    
    ALLOCATE(E_Rd(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for E_Rd. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(Vcmax25(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for Vcmax25. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(arJV(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for arJV. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(brJV(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for brJV. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(KmC25(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for KmC25. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(KmO25(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for KmO25. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(Sco25(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for Sco25. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF
    
    ALLOCATE(gm25(nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for gm25. We stop. We need nvm words = ',nvm 
       STOP 'pft_parameters_alloc' 
    END IF

    ALLOCATE(gamma_star25(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for gamma_star25. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(a1(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for a1. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(b1(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for b1. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(g0(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for g0. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(h_protons(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for h_protons. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(fpsir(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for fpsir. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(fQ(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for fQ. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(fpseudo(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for fpseudo. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(kp(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for kp. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(alpha(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for alpha. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(gbs(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for gbs. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(theta(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for theta. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(alpha_LL(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for alpha_LL. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(stress_vcmax(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for stress_vcmax. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF
    
    ALLOCATE(stress_gs(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for stress_gs. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF
    
    ALLOCATE(stress_gm(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for stress_gm. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(ext_coeff(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for ext_coeff. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(ext_coeff_vegetfrac(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for ext_coeff_vegetfrac. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(veget_ori_fixed_test_1(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for veget_ori_fixed_test_1. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(llaimax(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for llaimax. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(llaimin(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for llaimin. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(type_of_lai(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for type_of_lai. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

   ALLOCATE(agec_group(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for agec_group. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(age_class_bound(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for age_class_bound. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF
   age_class_bound(:)=0.

   ALLOCATE(start_index(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for start_index. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(nagec_pft(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for nagec_pft. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(leaf_tab(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for leaf_tab. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

    ALLOCATE(pref_soil_veg(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for pref_soil_veg. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(pheno_model(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for pheno_model. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(is_deciduous(nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for is_deciduous. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(is_evergreen(nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for is_evergreen. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(is_needleleaf(nvm),stat=ier)  
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for is_needleleaf. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(is_tropical(nvm),stat=ier)   
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for is_tropical. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF


    !
    ! 2. Parameters used if ok_sechiba only
    !
    IF ( ok_sechiba ) THEN

       l_error = .FALSE.

       ALLOCATE(rstruct_const(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for rstruct_const. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(kzero(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for kzero. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(rveg_pft(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for rveg_pft. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(wmax_veg(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for wmax_veg. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(throughfall_by_pft(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for throughfall_by_pft. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(snowa_dec_vis(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for snowa_dec_nir. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(snowa_dec_nir(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for snowa_dec_nir. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(snowa_aged_vis(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for snowa_aged_vis. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(snowa_aged_nir(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for snowa_aged_nir. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(alb_leaf_vis(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for alb_leaf_vis. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(alb_leaf_nir(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for alb_leaf_nir. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

      !chaoyue+
      ALLOCATE(permafrost_veg_exists(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for permafrost_veg_exists. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF
      !chaoyue-

      IF( ok_bvoc ) THEN

          l_error = .FALSE.

          ALLOCATE(em_factor_isoprene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_isoprene. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_monoterpene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_monoterpene. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_apinene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_apinene. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_bpinene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_bpinene. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_limonene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_limonene. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_myrcene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_myrcene. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_sabinene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_sabinene. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_camphene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_camphene. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_3carene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_3carene. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_tbocimene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_tbocimene. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_othermonot(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_othermonot. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_sesquiterp(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_sesquiterp. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF


          ALLOCATE(em_factor_ORVOC(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_ORVOC. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_OVOC(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0)       
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_OVOC. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_MBO(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_MBO. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_methanol(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_methanol. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_acetone(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_acetone. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_acetal(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_acetal. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_formal(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_formal. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_acetic(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0)       
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_acetic. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_formic(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_formic. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_no_wet(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0)
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_no_wet. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_no_dry(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0)       
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_no_dry. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(Larch(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for Larch. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

       ENDIF ! (ok_bvoc) 

       !!!!! crop parameters

       ALLOCATE(ok_LAIdev(nvm),stat=ier)
       IF (ier /= 0) THEN
         WRITE(str, *) 'Error code=', ier
         CALL ipslerr_p(3, 'config_pft_parameters', 'Memory error location for variable', ' ok_LAIdev', str)
       ENDIF
       !!!! end crop parameters

    ENDIF !(ok_sechiba)

    !
    ! 3. Parameters used if ok_stomate only
    !
    IF ( ok_stomate ) THEN

       l_error = .FALSE.

       ALLOCATE(sla(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for sla. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(availability_fact(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for availability_fact. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(R0(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for R0. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(S0(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for S0. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(L0(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for L0. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(pheno_gdd_crit_c(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for pheno_gdd_crit_c. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(pheno_gdd_crit_b(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for pheno_gdd_crit_b. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(pheno_gdd_crit_a(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for pheno_gdd_crit_a. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(pheno_gdd_crit(nvm,3),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for pheno_gdd_crit. We stop. We need nvm words = ',nvm*3
          STOP 'pft_parameters_alloc'
       END IF
       pheno_gdd_crit(:,:) = zero

       ALLOCATE(pheno_moigdd_t_crit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for pheno_moigdd_t_crit. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(ngd_crit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for ngd_crit. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(ncdgdd_temp(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for ncdgdd_temp. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(hum_frac(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for hum_frac. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(hum_min_time(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for hum_min_time. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(tau_sap(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for tau_sap. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(tau_leafinit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for tau_leafinit. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(tau_fruit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for tau_fruit. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(ecureuil(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for ecureuil. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(alloc_min(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for alloc_min. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(alloc_max(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for alloc_max. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(demi_alloc(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for . We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(frac_growthresp(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for frac_growthresp. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(maint_resp_slope(nvm,3),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for maint_resp_slope. We stop. We need nvm*3 words = ',nvm*3
          STOP 'pft_parameters_alloc'
       END IF
       maint_resp_slope(:,:) = zero

       ALLOCATE(maint_resp_slope_c(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for maint_resp_slope_c. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(maint_resp_slope_b(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for maint_resp_slope_b. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(maint_resp_slope_a(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for maint_resp_slope_a. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(coeff_maint_zero(nvm,nparts),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for coeff_maint_zero. We stop. We need nvm*nparts words = ',nvm*nparts
          STOP 'pft_parameters_alloc'
       END IF
       coeff_maint_zero(:,:) = zero

       ALLOCATE(cm_zero_leaf(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for cm_zero_leaf. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(cm_zero_sapabove(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for cm_zero_sapabove. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(cm_zero_sapbelow(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for cm_zero_sapbelow. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(cm_zero_heartabove(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for cm_zero_heartabove. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF
       
       ALLOCATE(cm_zero_heartbelow(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for cm_zero_heartbelow. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(cm_zero_root(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for cm_zero_root. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(cm_zero_fruit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for cm_zero_fruit. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(cm_zero_carbres(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for cm_zero_carbres. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

      !spitfire
      ALLOCATE(dens_fuel(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for dens_fuel. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(f_sh(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for f_sh. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(crown_length(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for crown_length. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(BTpar1(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for BTpar1. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(BTpar2(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for BTpar2. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(r_ck(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for r_ck. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(p_ck(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for p_ck. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(ef_CO2(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for ef_CO2. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(ef_CO(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for ef_CO. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(ef_CH4(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for ef_CH4. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(ef_VOC(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for ef_VOC. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(ef_TPM(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for ef_TPM. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(ef_NOx(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for ef_NOx. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(me(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for me. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(fire_max_cf_100hr(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for fire_max_cf_100hr. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(fire_max_cf_1000hr(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for fire_max_cf_1000hr. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF
      !endspit

      ! grassland management
!gmjc
      ALLOCATE(is_grassland_manag(nvm),stat=ier)
      l_error = l_error .OR. (ier .NE. 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for is_grassland_manag. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF
      ALLOCATE(is_grassland_cut(nvm),stat=ier)
      l_error = l_error .OR. (ier .NE. 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for is_grassland_cut. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF
      ALLOCATE(is_grassland_grazed(nvm),stat=ier)
      l_error = l_error .OR. (ier .NE. 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for is_grassland_grazed. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF
      ALLOCATE(management_intensity(nvm),stat=ier)
      l_error = l_error .OR. (ier .NE. 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for management_intensity. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF
      ALLOCATE(management_start(nvm),stat=ier)
      l_error = l_error .OR. (ier .NE. 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for management_start. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF
      ALLOCATE(deposition_start(nvm),stat=ier)
      l_error = l_error .OR. (ier .NE. 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for deposition_start. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF
      ALLOCATE(nb_year_management(nvm),stat=ier)
      l_error = l_error .OR. (ier .NE. 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for nb_year_management. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF
      ALLOCATE(sla_max(nvm),stat=ier)
      l_error = l_error .OR. (ier .NE. 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for sla_max. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF
      ALLOCATE(sla_min(nvm),stat=ier)
      l_error = l_error .OR. (ier .NE. 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for sla_min. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF
!end gmjc

      ALLOCATE(is_bioe1(nvm),stat=ier)
      l_error = l_error .OR. (ier .NE. 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for is_bioe1. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

       ALLOCATE(flam(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for . We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF
       ALLOCATE(resist(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for resist. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(coeff_lcchange_1(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for coeff_lcchange_1. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(coeff_lcchange_10(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for coeff_lcchange_10. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(coeff_lcchange_100(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for coeff_lcchange_100. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(coeff_indwood_1(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for coeff_indwood_1. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(coeff_indwood_10(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for coeff_indwood_10. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(coeff_indwood_100(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for coeff_indwood_100. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(lai_max_to_happy(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for lai_max_to_happy. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(lai_max(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for lai_max. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(pheno_type(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for pheno_type. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(leaffall(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for leaffall. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(leafagecrit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for leafagecrit. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(senescence_type(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for . We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(senescence_hum(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for senescence_hum. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(nosenescence_hum(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for nosenescence_hum. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(max_turnover_time(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for max_turnover_time. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(min_turnover_time(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for min_turnover_time. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(min_leaf_age_for_senescence(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for min_leaf_age_for_senescence. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(senescence_temp_c(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for senescence_temp_c. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(senescence_temp_b(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for senescence_temp_b. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(senescence_temp_a(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for senescence_temp_a. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(senescence_temp(nvm,3),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for senescence_temp. We stop. We need nvm*3 words = ',nvm*3
          STOP 'pft_parameters_alloc'
       END IF
       senescence_temp(:,:) = zero

       ALLOCATE(gdd_senescence(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for gdd_senescence. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(residence_time(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for residence_time. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(tmin_crit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for tmin_crit. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(tcm_crit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for tcm_crit. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

!qcj++ peatland
       ALLOCATE(wtpwet_crit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for wtpwet_crit. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(wtpdry_crit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for wtpdry_crit. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(wtp_crit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for wtp_crit. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(wt_mortality(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for wt_mortality. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(lai_initmin(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for . We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(bm_sapl(nvm,nparts,nelements),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for bm_sapl. We stop. We need nvm*nparts*nelements words = ',& 
               &  nvm*nparts*nelements
          STOP 'pft_parameters_alloc'
       END IF

!pss+
      ALLOCATE(rdepth_v(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for rdepth_v. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(sdepth_v(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for sdepth_v. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(tveg_v(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for tveg_v. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF
!pss- 

!!!!! crop parameters

      ALLOCATE(SP_codeplante(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)   
      ALLOCATE(SP_stade0(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)   
      ALLOCATE(SP_iplt0(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0) 
      IF (cyc_rot_max .GT. 1) THEN
         ALLOCATE(SP_iplt1(nvm), stat=ier)
         l_error = l_error .OR. (ier /= 0)
      ENDIF
      IF (cyc_rot_max .GT. 2) THEN
         ALLOCATE(SP_iplt2(nvm), stat=ier)
         l_error = l_error .OR. (ier /= 0)
      ENDIF
      ALLOCATE(SP_nbox(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0) 
      ALLOCATE(SP_iwater(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)

      ALLOCATE(SP_codesimul(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0) 
      ALLOCATE(SP_codelaitr(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)

      ALLOCATE(SP_slamax(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0) 
      ALLOCATE(SP_slamin(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)

      ALLOCATE(SP_codeperenne(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0) 
      ALLOCATE(SP_codcueille(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)

      ALLOCATE(SP_codegdh(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0) 
      ALLOCATE(SP_codetemp(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)

      ALLOCATE(SP_coderetflo(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0) 
      ALLOCATE(SP_codeinnact(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)

      ALLOCATE(SP_codeh2oact(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0) 
      ALLOCATE(SP_stressdev(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)

      ALLOCATE(SP_innlai(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_innsenes(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_codebfroid(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)

      ALLOCATE(SP_codephot(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_codedormance(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_codefauche(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_codetempfauche(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_codlainet(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_codeindetermin(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_codeinitprec(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_culturean(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_jvc(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_tfroid(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_ampfroid(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_jvcmini(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_tgmin(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_stpltger(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_profsem(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_propjgermin(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_tdmax(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_nbjgerlim(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_densitesem(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_vigueurbat(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_codepluiepoquet(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_codehypo(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_elmax(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_belong(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_celong(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_nlevlim1(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_nlevlim2(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_codrecolte(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_variete(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_codegermin(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)

      ALLOCATE(S_codeulaivernal(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_swfacmin(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_neffmax(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_nsatrat(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)


      ALLOCATE(SP_laiplantule(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)

      ALLOCATE(SP_vlaimax(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_stlevamf(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_stdrpmat(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_stamflax(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_udlaimax(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_laicomp(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_adens(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_bdens(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_tcxstop(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_tcmax(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_tcmin(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_dlaimax(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_dlaimin(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_pentlaimax(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_tigefeuil(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      
      ALLOCATE(SP_stlaxsen(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_stsenlan(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_stlevdrp(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_stflodrp(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_stdrpdes(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_phyllotherme(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)

      ALLOCATE(SP_lai0(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_tustressmin(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)



      ! STICS:: LAI SENESCENCE
      ALLOCATE(SP_nbfgellev(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_ratiodurvieI(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_durvieF(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_ratiosen(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_tdmin(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
    
      ! STICS:: F_humerac

      ALLOCATE(SP_sensrsec(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ! STICS:: GEL

      ALLOCATE(SP_codgellev(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_codgeljuv(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_codgelveg(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_tletale(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_tdebgel(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_tgellev10(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_tgellev90(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)


      ALLOCATE(SP_tgeljuv10(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_tgeljuv90(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_tgelveg10(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_tgelveg90(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)






      ! STICS:: Photoperiod

      ALLOCATE(SP_sensiphot(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_phosat(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_phobase(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)

      ! STICS:: CARBON ALLOCATION

      ALLOCATE(SP_stoprac(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_zracplantule(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_codtrophrac(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_repracpermax(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_repracpermin(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_krepracperm(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_repracseumax(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_repracseumin(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_krepracseu(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_codetemprac(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_codedyntalle(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_nbjgrain(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_maxgs(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_codgelflo(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_tgelflo10(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_tgelflo90(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_cgrain(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_cgrainv0(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_nbgrmax(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_nbgrmin(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_codazofruit(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_codeir(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_vitircarb(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_irmax(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_vitircarbT(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_codetremp(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_tminremp(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_tmaxremp(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      ALLOCATE(SP_pgrainmaxi(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)

      !! for dynamic nitrogen process

      ALLOCATE(SP_DY_INN(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)

      ALLOCATE(SP_avenfert(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)

      IF (l_error) THEN
          STOP 'pft_alloc : error in memory allocation of crop pft parameters'
      ENDIF
!!!!! end crop parameters

       ALLOCATE(migrate(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for migrate. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(maxdia(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for maxdia. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(cn_sapl(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for cn_sapl. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(leaf_timecst(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for leaf_timecst. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(leaflife_tab(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for leaflife_tab. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

    ENDIF ! (ok_stomate)

  END SUBROUTINE pft_parameters_alloc

!! ================================================================================================================================
!! SUBROUTINE   : config_pft_parameters 
!!
!>\BRIEF          This subroutine will read the imposed values for the global pft
!! parameters (sechiba + stomate). It is not called if IMPOSE_PARAM is set to NO.
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE config_pft_parameters

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.4 Local variable

    INTEGER(i_std) :: jv, ivm                   !! Index (untiless)

    !_ ================================================================================================================================ 


    !
    ! Vegetation structure
    !

     !Config  Key  = PERMAFROST_VEG_EXISTS
     !Config  Desc = Is the vegetation type a permafrost vegetation
     !adaptation ?
     !Config  if  = OK_SECHIBA
     !Config  Def  = y, y, y, y, y, y, y, y, y, y, y, y, y
     !Config  Help =
     !Config  Units = NONE
     CALL getin_p('PERMAFROST_VEG_EXISTS', permafrost_veg_exists)
     
      !Config Key   = SECHIBA_LAI
      !Config Desc  = laimax for maximum lai(see also type of lai interpolation)
      !Config if    = OK_SECHIBA or IMPOSE_VEG
      !Config Def   = 0., 8., 8., 4., 4.5, 4.5, 4., 4.5, 4., 2., 2., 2., 2.
      !Config Help  = Maximum values of lai used for interpolation of the lai map
      !Config Units = [m^2/m^2]
      CALL getin_p('SECHIBA_LAI',llaimax)

    !! Redefine the values for is_tree, is_deciduous, is_needleleaf, is_evergreen if values have been modified
    !! in run.def

    is_tree(:) = .FALSE.
    DO jv = 1,nvm
       IF ( leaf_tab(jv) <= 2 ) is_tree(jv) = .TRUE.
    END DO
    !
    is_deciduous(:) = .FALSE.
    DO jv = 1,nvm
       IF ( is_tree(jv) .AND. (pheno_model(jv) /= "none") ) is_deciduous(jv) = .TRUE.
    END DO
    !
    is_evergreen(:) = .FALSE.
    DO jv = 1,nvm
       IF ( is_tree(jv) .AND. (pheno_model(jv) == "none") ) is_evergreen(jv) = .TRUE.
    END DO
    !
    is_needleleaf(:) = .FALSE.
    DO jv = 1,nvm
       IF ( leaf_tab(jv) == 2 ) is_needleleaf(jv) = .TRUE.
    END DO


    !Config Key   = SECHIBA_LAI
    !Config Desc  = laimax for maximum lai(see also type of lai interpolation)
    !Config if    = OK_SECHIBA or IMPOSE_VEG
    !Config Def   = 0., 8., 8., 4., 4.5, 4.5, 4., 4.5, 4., 2., 2., 2., 2.
    !Config Help  = Maximum values of lai used for interpolation of the lai map
    !Config Units = [m^2/m^2]
    CALL getin_p('SECHIBA_LAI',llaimax)

    !Config Key   = LLAIMIN
    !Config Desc  = laimin for minimum lai(see also type of lai interpolation)
    !Config if    = OK_SECHIBA or IMPOSE_VEG
    !Config Def   = 0., 8., 0., 4., 4.5, 0., 4., 0., 0., 0., 0., 0., 0.
    !Config Help  = Minimum values of lai used for interpolation of the lai map
    !Config Units = [m^2/m^2]
    CALL getin_p('LLAIMIN',llaimin)

    !Config Key   = SLOWPROC_HEIGHT
    !Config Desc  = prescribed height of vegetation 
    !Config if    = OK_SECHIBA
    !Config Def   = 0., 30., 30., 20., 20., 20., 15., 15., 15., .5, .6, 1., 1.
    !Config Help  =
    !Config Units = [m] 
    CALL getin_p('SLOWPROC_HEIGHT',height_presc)

      !Config Key   = NATURAL
      !Config Desc  = natural? 
      !Config if    = OK_SECHIBA, OK_STOMATE
      !Config Def   = y, y, y, y, y, y, y, y, y, y, y, n, n, n 
      !Config Help  =
      !Config Units = [BOOLEAN]
      CALL getin_p('NATURAL',natural)
! dgvmjc
    !Config Key   = PASTURE
    !Config Desc  = pasture?
    !Config if    = OK_SECHIBA, OK_STOMATE
    !Config Def   = y, y, y, y, y, y, y, y, y, y, y, n, n
    !Config Help  =
    !Config Units = [BOOLEAN]
    CALL getin_p('PASTURE',pasture)
! end dgvmjc

    !
    !Config Key   = RATIO_Z0M_Z0H
    !Config Desc  = Ratio between z0m and z0h
    !Config Def   = 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 
    !Config if    = OK_SECHIBA
    !Config Help  = 
    !Config Units = [-]
    CALL getin_p('RATIO_Z0M_Z0H',ratio_z0m_z0h)

      !Config Key   = IS_C4
      !Config Desc  = flag for C4 vegetation types
      !Config if    = OK_SECHIBA or OK_STOMATE
      !Config Def   = n, n, n, n, n, n, n, n, n, n, n, y, n, y, n
      !Config Help  =
      !Config Units = [BOOLEAN]
      CALL getin_p('IS_C4',is_c4)
  
      !Config Key   = TMIN_CRIT
      !Config Desc  = minimum temperature limitation, below which the mortality rate will increase, Zhu et al. 2015
      !Config if    = OK_STOMATE
      !Config Def   = undef, 0.0, 0.0, -30.0, -14.0, -30.0, -45.0, -45.0, -60.0, undef, undef, undef, undef
      !Config Help  =
      !Config Units = [degree C]
      CALL getin_p('TMIN_CRIT',tmin_crit) 

    !Config Key   = TYPE_OF_LAI
    !Config Desc  = Type of behaviour of the LAI evolution algorithm 
    !Config if    = OK_SECHIBA
    !Config Def   = inter, inter, inter, inter, inter, inter, inter, inter, inter, inter, inter, inter, inter
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('TYPE_OF_LAI',type_of_lai)

    !Config Key   = NATURAL
    !Config Desc  = natural? 
    !Config if    = OK_SECHIBA, OK_STOMATE
    !Config Def   = y, y, y, y, y, y, y, y, y, y, y, n, n 
    !Config Help  =
    !Config Units = [BOOLEAN]
    CALL getin_p('NATURAL',natural)


    !
    ! Photosynthesis
    !

    !Config Key   = VCMAX_FIX
    !Config Desc  = values used for vcmax when STOMATE is not activated
    !Config if    = OK_SECHIBA and NOT(OK_STOMATE)
    !Config Def   = 0., 40., 50., 30., 35., 40.,30., 40., 35., 60., 60., 70., 70.
    !Config Help  =
    !Config Units = [micromol/m^2/s] 
    CALL getin_p('VCMAX_FIX',vcmax_fix)

    !Config Key   = DOWNREG_CO2
    !Config Desc  = coefficient for CO2 downregulation (unitless)
    !Config if    = OK_CO2
    !Config Def   = 0., 0.38, 0.38, 0.28, 0.28, 0.28, 0.22, 0.22, 0.22, 0.26, 0.26, 0.26, 0.26
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('DOWNREG_CO2',downregulation_co2_coeff)

    !Config Key   = DOWNREG_CO2_NEW
    !Config Desc  = coefficient for CO2 downregulation (unitless)
    !Config if    = OK_CO2 and DOWNREGULATION_CO2_NEW
    !Config Def   = 0., 0.35, 0.35, 0.26, 0.26, 0.26, 0.20, 0.20, 0.20, 0.24,
    !0.03, 0.24, 0.03
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('DOWNREG_CO2_NEW',downregulation_co2_coeff_new)

    !Config Key   = E_KmC
    !Config Desc  = Energy of activation for KmC
    !Config if    = OK_CO2
    !Config Def   = undef,  79430., 79430., 79430., 79430., 79430., 79430., 79430., 79430., 79430., 79430., 79430., 79430.
    !Config Help  = See Medlyn et al. (2002) 
    !Config Units = [J mol-1]
    CALL getin_p('E_KMC',E_KmC)

    !Config Key   = E_KmO
    !Config Desc  = Energy of activation for KmO
    !Config if    = OK_CO2
    !Config Def   = undef, 36380.,  36380.,  36380.,  36380.,  36380., 36380., 36380., 36380., 36380., 36380., 36380., 36380.
    !Config Help  = See Medlyn et al. (2002) 
    !Config Units = [J mol-1]
    CALL getin_p('E_KMO',E_KmO)

    !Config Key   = E_Sco
    !Config Desc  = Energy of activation for Sco
    !Config if    = OK_CO2
    !Config Def   = undef, -24460., -24460., -24460., -24460., -24460., -24460., -24460., -24460., -24460., -24460., -24460., -24460.
    !Config Help  = See Table 2 of Yin et al. (2009) - Value for C4 plants is not mentioned - We use C3 for all plants
    !Config Units = [J mol-1]
    CALL getin_p('E_SCO',E_Sco)
    
    !Config Key   = E_gamma_star
    !Config Desc  = Energy of activation for gamma_star
    !Config if    = OK_CO2
    !Config Def   = undef, 37830.,  37830.,  37830.,  37830.,  37830., 37830., 37830., 37830., 37830., 37830., 37830., 37830.
    !Config Help  = See Medlyn et al. (2002) from Bernacchi al. (2001) 
    !Config Units = [J mol-1]
    CALL getin_p('E_GAMMA_STAR',E_gamma_star)

    !Config Key   = E_Vcmax
    !Config Desc  = Energy of activation for Vcmax
    !Config if    = OK_CO2
    !Config Def   = undef, 71513., 71513., 71513., 71513., 71513., 71513., 71513., 71513., 71513., 67300., 71513., 67300.
    !Config Help  = See Table 2 of Yin et al. (2009) for C4 plants and Kattge & Knorr (2007) for C3 plants (table 3)
    !Config Units = [J mol-1]
    CALL getin_p('E_VCMAX',E_Vcmax)

    !Config Key   = E_Jmax
    !Config Desc  = Energy of activation for Jmax
    !Config if    = OK_CO2
    !Config Def   = undef, 49884., 49884., 49884., 49884., 49884., 49884., 49884., 49884., 49884., 77900., 49884., 77900. 
    !Config Help  = See Table 2 of Yin et al. (2009) for C4 plants and Kattge & Knorr (2007) for C3 plants (table 3)
    !Config Units = [J mol-1]
    CALL getin_p('E_JMAX',E_Jmax)

    !Config Key   = aSV
    !Config Desc  = a coefficient of the linear regression (a+bT) defining the Entropy term for Vcmax
    !Config if    = OK_CO2
    !Config Def   = undef, 668.39, 668.39, 668.39, 668.39, 668.39, 668.39, 668.39, 668.39, 668.39, 641.64, 668.39, 641.64 
    !Config Help  = See Table 3 of Kattge & Knorr (2007) - For C4 plants, we assume that there is no acclimation and that at for a temperature of 25°C, aSV is the same for both C4 and C3 plants (no strong jusitification - need further parametrization)
    !Config Units = [J K-1 mol-1]
    CALL getin_p('ASV',aSV)

    !Config Key   = bSV
    !Config Desc  = b coefficient of the linear regression (a+bT) defining the Entropy term for Vcmax
    !Config if    = OK_CO2
    !Config Def   = undef, -1.07, -1.07, -1.07, -1.07, -1.07, -1.07, -1.07, -1.07, -1.07, 0., -1.07, 0. 
    !Config Help  = See Table 3 of Kattge & Knorr (2007) - For C4 plants, we assume that there is no acclimation
    !Config Units = [J K-1 mol-1 °C-1]
    CALL getin_p('BSV',bSV)

    !Config Key   = TPHOTO_MIN
    !Config Desc  = minimum photosynthesis temperature (deg C)
    !Config if    = OK_STOMATE
    !Config Def   = undef,  -4., -4., -4., -4.,-4.,-4., -4., -4., -4., -4., -4., -4.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('TPHOTO_MIN',tphoto_min)

    !Config Key   = TPHOTO_MAX
    !Config Desc  = maximum photosynthesis temperature (deg C)
    !Config if    = OK_STOMATE
    !Config Def   = undef, 55., 55., 55., 55., 55., 55., 55., 55., 55., 55., 55., 55.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('TPHOTO_MAX',tphoto_max)

    !Config Key   = aSJ
    !Config Desc  = a coefficient of the linear regression (a+bT) defining the Entropy term for Jmax
    !Config if    = OK_CO2
    !Config Def   = undef, 659.70, 659.70, 659.70, 659.70, 659.70, 659.70, 659.70, 659.70, 659.70, 630., 659.70, 630. 
    !Config Help  = See Table 3 of Kattge & Knorr (2007) - and Table 2 of Yin et al. (2009) for C4 plants
    !Config Units = [J K-1 mol-1]
    CALL getin_p('ASJ',aSJ)

    !Config Key   = bSJ
    !Config Desc  = b coefficient of the linear regression (a+bT) defining the Entropy term for Jmax
    !Config if    = OK_CO2
    !Config Def   = undef, -0.75, -0.75, -0.75, -0.75, -0.75, -0.75, -0.75, -0.75, -0.75, 0., -0.75, 0. 
    !Config Help  = See Table 3 of Kattge & Knorr (2007) - For C4 plants, we assume that there is no acclimation
    !Config Units = [J K-1 mol-1 °C-1]
    CALL getin_p('BSJ',bSJ)

    !Config Key   = D_Vcmax
    !Config Desc  = Energy of deactivation for Vcmax
    !Config if    = OK_CO2
    !Config Def   = undef, 200000., 200000., 200000., 200000., 200000., 200000., 200000., 200000., 200000., 192000., 200000., 192000.
    !Config Help  = Medlyn et al. (2002) also uses 200000. for C3 plants (same value than D_Jmax). 'Consequently', we use the value of D_Jmax for C4 plants.
    !Config Units = [J mol-1]
    CALL getin_p('D_VCMAX',D_Vcmax)

    !Config Key   = D_Jmax
    !Config Desc  = Energy of deactivation for Jmax
    !Config if    = OK_CO2
    !Config Def   = undef, 200000., 200000., 200000., 200000., 200000., 200000., 200000., 200000., 200000., 192000., 200000., 192000.
    !Config Help  = See Table 2 of Yin et al. (2009)
    !Config Units = [J mol-1]
    CALL getin_p('D_JMAX',D_Jmax)
    
    !Config Key   = E_gm 
    !Config Desc  = Energy of activation for gm 
    !Config if    = OK_CO2 
    !Config Def   = undef, 49600., 49600., 49600., 49600., 49600., 49600., 49600., 49600., 49600., undef, 49600., undef 
    !Config Help  = See Table 2 of Yin et al. (2009) 
    !Config Units = [J mol-1] 
    CALL getin_p('E_GM',E_gm) 
    
    !Config Key   = S_gm 
    !Config Desc  = Entropy term for gm 
    !Config if    = OK_CO2 
    !Config Def   = undef, 1400., 1400., 1400., 1400., 1400., 1400., 1400., 1400., 1400., undef, 1400., undef 
    !Config Help  = See Table 2 of Yin et al. (2009) 
    !Config Units = [J K-1 mol-1] 
    CALL getin_p('S_GM',S_gm) 
    
    !Config Key   = D_gm 
    !Config Desc  = Energy of deactivation for gm 
    !Config if    = OK_CO2 
    !Config Def   = undef, 437400., 437400., 437400., 437400., 437400., 437400., 437400., 437400., 437400., undef, 437400., undef 
    !Config Help  = See Table 2 of Yin et al. (2009) 
    !Config Units = [J mol-1] 
    CALL getin_p('D_GM',D_gm) 
    
    !Config Key   = E_Rd
    !Config Desc  = Energy of activation for Rd
    !Config if    = OK_CO2
    !Config Def   = undef, 46390., 46390., 46390., 46390., 46390., 46390., 46390., 46390., 46390., 46390., 46390., 46390.
    !Config Help  = See Table 2 of Yin et al. (2009)
    !Config Units = [J mol-1]
    CALL getin_p('E_RD',E_Rd)

    !Config Key   = VCMAX25
    !Config Desc  = Maximum rate of Rubisco activity-limited carboxylation at 25°C
    !Config if    = OK_STOMATE
    !Config Def   = undef, 50., 65., 35., 45., 55., 35., 45., 35., 70., 70., 70., 70.
    !Config Help  =
    !Config Units = [micromol/m^2/s]
    CALL getin_p('VCMAX25',Vcmax25)

    !Config Key   = ARJV
    !Config Desc  = a coefficient of the linear regression (a+bT) defining the Jmax25/Vcmax25 ratio 
    !Config if    = OK_STOMATE
    !Config Def   = undef, 2.59, 2.59, 2.59, 2.59, 2.59, 2.59, 2.59, 2.59, 2.59, 1.715, 2.59, 1.715
    !Config Help  = See Table 3 of Kattge & Knorr (2007) - For C4 plants, we assume that there is no acclimation and that for a temperature of 25°C, aSV is the same for both C4 and C3 plants (no strong jusitification - need further parametrization)
    !Config Units = [mu mol e- (mu mol CO2)-1]
    CALL getin_p('ARJV',arJV)

    !Config Key   = BRJV
    !Config Desc  = b coefficient of the linear regression (a+bT) defining the Jmax25/Vcmax25 ratio 
    !Config if    = OK_STOMATE
    !Config Def   = undef, -0.035, -0.035, -0.035, -0.035, -0.035, -0.035, -0.035, -0.035, -0.035, 0., -0.035, 0.
    !Config Help  = See Table 3 of Kattge & Knorr (2007) -  We assume No acclimation term for C4 plants
    !Config Units = [(mu mol e- (mu mol CO2)-1) (°C)-1]
    CALL getin_p('BRJV',brJV)

    !Config Key   = KmC25
    !Config Desc  = Michaelis–Menten constant of Rubisco for CO2 at 25°C
    !Config if    = OK_CO2
    !Config Def   = undef, 404.9, 404.9, 404.9, 404.9, 404.9, 404.9, 404.9, 404.9, 404.9, 650., 404.9, 650.
    !Config Help  = See Table 2 of Yin et al. (2009) for C4 plants and Medlyn et al. (2002) for C3 plants
    !Config Units = [ubar]
    CALL getin_p('KMC25',KmC25)

    !Config Key   = KmO25
    !Config Desc  = Michaelis–Menten constant of Rubisco for O2 at 25°C
    !Config if    = OK_CO2
    !Config Def   = undef, 278400., 278400., 278400., 278400., 278400., 278400., 278400., 278400., 278400., 450000., 278400., 450000.
    !Config Help  = See Table 2 of Yin et al. (2009) for C4 plants and Medlyn et al. (2002) for C3 plants
    !Config Units = [ubar]
    CALL getin_p('KMO25',KmO25)

    !Config Key   = Sco25
    !Config Desc  = Relative CO2 /O2 specificity factor for Rubisco at 25Â°C
    !Config if    = OK_CO2
    !Config Def   = undef, 2800., 2800., 2800., 2800., 2800., 2800., 2800., 2800., 2800., 2590., 2800., 2590.
    !Config Help  = See Table 2 of Yin et al. (2009)
    !Config Units = [bar bar-1]
    CALL getin_p('SCO25',Sco25)
    
    !Config Key   = gm25 
    !Config Desc  = Mesophyll diffusion conductance at 25ÃÂ°C 
    !Config if    = OK_CO2 
    !Config Def   = undef, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, undef, 0.4, undef 
    !Config Help  = See legend of Figure 6 of Yin et al. (2009) and review by Flexas et al. (2008) - gm is not used for C4 plants 
    !Config Units = [mol m-2 s-1 bar-1] 
    CALL getin_p('GM25',gm25) 
    
    !Config Key   = gamma_star25
    !Config Desc  = Ci-based CO2 compensation point in the absence of Rd at 25°C (ubar)
    !Config if    = OK_CO2
    !Config Def   = undef, 42.75, 42.75, 42.75, 42.75, 42.75, 42.75, 42.75, 42.75, 42.75, 42.75, 42.75, 42.75
    !Config Help  = See Medlyn et al. (2002) for C3 plants - For C4 plants, we use the same value (probably uncorrect)
    !Config Units = [ubar]
    CALL getin_p('gamma_star25',gamma_star25)

    !Config Key   = a1
    !Config Desc  = Empirical factor involved in the calculation of fvpd
    !Config if    = OK_CO2
    !Config Def   = undef, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.72, 0.85, 0.72
    !Config Help  = See Table 2 of Yin et al. (2009)
    !Config Units = [-]
    CALL getin_p('A1',a1)

    !Config Key   = b1
    !Config Desc  = Empirical factor involved in the calculation of fvpd
    !Config if    = OK_CO2
    !Config Def   = undef, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.20, 0.14, 0.20
    !Config Help  = See Table 2 of Yin et al. (2009)
    !Config Units = [-]
    CALL getin_p('B1',b1)

    !Config Key   = g0
    !Config Desc  = Residual stomatal conductance when irradiance approaches zero 
    !Config if    = OK_CO2
    !Config Def   = undef, 0.00625, 0.00625, 0.00625, 0.00625, 0.00625, 0.00625, 0.00625, 0.00625, 0.00625, 0.01875, 0.00625, 0.01875 
    !Config Help  = Value from ORCHIDEE - No other reference.
    !Config Units = [mol m−2 s−1 bar−1]
    CALL getin_p('G0',g0)

    !Config Key   = h_protons
    !Config Desc  = Number of protons required to produce one ATP
    !Config if    = OK_CO2
    !Config Def   = undef, 4., 4., 4., 4., 4., 4., 4., 4., 4., 4., 4., 4. 
    !Config Help  = See Table 2 of Yin et al. (2009) - h parameter
    !Config Units = [mol mol-1]
    CALL getin_p('H_PROTONS',h_protons)

    !Config Key   = fpsir
    !Config Desc  = Fraction of PSII e− transport rate partitioned to the C4 cycle
    !Config if    = OK_CO2
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 0.4, undef, 0.4 
    !Config Help  = See Table 2 of Yin et al. (2009)
    !Config Units = [-]
    CALL getin_p('FPSIR',fpsir)

    !Config Key   = fQ
    !Config Desc  = Fraction of electrons at reduced plastoquinone that follow the Q-cycle
    !Config if    = OK_CO2
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 1., undef, 1.
    !Config Help  = See Table 2 of Yin et al. (2009) - Values for C3 plants are not used
    !Config Units = [-]
    CALL getin_p('FQ',fQ)

    !Config Key   = fpseudo
    !Config Desc  = Fraction of electrons at PSI that follow pseudocyclic transport 
    !Config if    = OK_CO2
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 0.1, undef, 0.1
    !Config Help  = See Table 2 of Yin et al. (2009) - Values for C3 plants are not used
    !Config Units = [-]
    CALL getin_p('FPSEUDO',fpseudo)

    !Config Key   = kp
    !Config Desc  = Initial carboxylation efficiency of the PEP carboxylase
    !Config if    = OK_CO2
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 0.7, undef, 0.7
    !Config Help  = See Table 2 of Yin et al. (2009) 
    !Config Units = [mol m−2 s−1 bar−1]
    CALL getin_p('KP',kp)

    !Config Key   = alpha
    !Config Desc  = Fraction of PSII activity in the bundle sheath
    !Config if    = OK_CO2
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 0.1, undef, 0.1
    !Config Help  = See legend of Figure 6 of Yin et al. (2009)
    !Config Units = [-]
    CALL getin_p('ALPHA',alpha)

    !Config Key   = gbs
    !Config Desc  = Bundle-sheath conductance
    !Config if    = OK_CO2
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 0.003, undef, 0.003
    !Config Help  = See legend of Figure 6 of Yin et al. (2009)
    !Config Units = [mol m−2 s−1 bar−1]
    CALL getin_p('GBS',gbs)

    !Config Key   = theta
    !Config Desc  = Convexity factor for response of J to irradiance
    !Config if    = OK_CO2
    !Config Def   = undef, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7
    !Config Help  = See Table 2 of Yin et al. (2009)   
    !Config Units = [−]
    CALL getin_p('THETA',theta)

    !Config Key   = STRESS_VCMAX
    !Config Desc  = Stress on vcmax
    !Config if    = OK_SECHIBA or OK_STOMATE
    !Config Def   = 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('STRESS_VCMAX', stress_vcmax)
    
    !Config Key   = STRESS_GS
    !Config Desc  = Stress on gs
    !Config if    = OK_SECHIBA or OK_STOMATE
    !Config Def   = 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('STRESS_GS', stress_gs)
    
    !Config Key   = STRESS_GM
    !Config Desc  = Stress on gm
    !Config if    = OK_SECHIBA or OK_STOMATE
    !Config Def   = 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('STRESS_GM', stress_gm)

    !Config Key   = EXT_COEFF
    !Config Desc  = extinction coefficient of the Monsi&Seaki relationship (1953)
    !Config if    = OK_SECHIBA or OK_STOMATE
    !Config Def   = .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('EXT_COEFF',ext_coeff)

!!!!! xuhui for crop rotation
      CALL getin_p('NSTM',nstm)
      IF ( (nstm .LT. 2) .OR. (nstm .GT. 20) ) THEN
        WRITE(numout,*) 'bad value for nstm: ',nstm
        WRITE(numout,*) 'set nstm as 6 (default)'
      ENDIF
      IF ( nstm .LT. 6) THEN
        ! default value did not work properly
        WHERE (pref_soil_veg(:) .GT. nstm)
            pref_soil_veg(:) = nstm
        ENDWHERE
!        Make sure model can work anyway
      ENDIF
!!!!! end crop rotation, xuhui

      !Config Key   = PREF_SOIL_VEG
      !Config Desc  = The soil tile number for each vegetation
      !Config if    = OK_SECHIBA or OK_STOMATE
      !Config Def   = 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 4, 5, 6
      !Config Help  = Gives the number of the soil tile on which we will
      !Config         put each vegetation. This allows to divide the hydrological column
      !Config Units = [-]        
      CALL getin_p('PREF_SOIL_VEG',pref_soil_veg)
      !! we should judge whether the pref_soil_veg is larger than nstm

!!!!! crop parameters      
      !Config Key   = IRRIG_THRESHOLD
      !Config if    = OK_SECHIBA, OK_STOMATE, DO_IRRIGATION
      CALL getin_p('IRRIG_THRESHOLD',irrig_threshold) 

      CALL getin_p('IRRIG_FULFILL',irrig_fulfill) 
      
      
      !
      !
      !Config Key   = OK_LAIDEV
      !Config Desc  = whether or not we open the STICS module
      !Config if    = OK_STOMATE and OK_SECHIBA
      !Config Def   = .false., .false., .false., .false., .false., .false., .false., .false., .false., .false., .false.,  .true., .true.
      !Config Help  =
      !Config Units = [C]
      CALL getin_p('OK_LAIDEV',ok_LAIdev)

!!!qcj++ peatland
      CALL getin_p('IS_PEAT',is_peat)
      CALL getin_p('IS_CROPPEAT',is_croppeat)
      CALL getin_p('IS_SHRUBPEAT',is_shrubpeat)
      CALL getin_p('IS_MOSSPEAT',is_mosspeat)
      CALL getin_p('IS_MINERALWET',is_mineralwet)

      !
      ! Vegetation - Age classes
      !
      !Config Key   = NVMAP
      !Config Desc  = The number of PFTs if we ignore age classes.  
      !               i.e., the number of metaclasses.
      !Config if    = OK_SECHIBA or OK_STOMATE
      !Config Def   = nvm
      !Config Help  = Gives the total number of PFTs ignoring age classes.
      !Config Units = [-]  
      nvmap=nvm
      IF (use_age_class) THEN
         CALL getin_p('GLUC_NVMAP',nvmap)
         IF(nvmap == nvm .AND. .NOT. SingleAgeClass)THEN
            WRITE(numout,*) 'WARNING: The age classes will be used, but'
            WRITE(numout,*) '         the input file indicates that none of the PFTs have age classes.'
            WRITE(numout,*) '         You should change nvmap.'
         ENDIF

         !Config Key   = AGEC_GROUP
         !Config Desc  = The group that each PFT belongs to.  
         !Config if    = OK_SECHIBA or OK_STOMATE
         !Config Def   = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13
         !Config Help  = The group that each PFT belongs to.  If you are not using age classes, this
         !Config         is just equal to the number of the PFT.
         !Config Units = [-]   
         DO ivm=1,nvm
            agec_group(ivm)=ivm
         ENDDO
         CALL getin_p('GLUC_AGEC_GROUP',agec_group)
         IF (.NOT. use_bound_spa) THEN
           CALL getin_p('GLUC_AGE_CLASS_BOUND',age_class_bound)
         ENDIF
      ENDIF

    !Config Key   = EXT_COEFF_VEGETFRAC
    !Config Desc  = extinction coefficient used for the calculation of the bare soil fraction 
    !Config if    = OK_SECHIBA or OK_STOMATE
    !Config Def   = 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('EXT_COEFF_VEGETFRAC',ext_coeff_vegetfrac)

    !
    ! Water-hydrology - sechiba
    !

    !Config Key   = HYDROL_HUMCSTE
    !Config Desc  = Root profile
    !Config Def   = humcste_ref2m or humcste_ref4m depending on zmaxh
    !Config if    = OK_SECHIBA
    !Config Help  = See module constantes_mtc for different default values
    !Config Units = [m]
    CALL getin_p('HYDROL_HUMCSTE',humcste)

    !
    ! Soil - vegetation
    !

    !Config Key   = PREF_SOIL_VEG
    !Config Desc  = The soil tile number for each vegetation
    !Config if    = OK_SECHIBA or OK_STOMATE
    !Config Def   = 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3
    !Config Help  = Gives the number of the soil tile on which we will
    !Config         put each vegetation. This allows to divide the hydrological column
    !Config Units = [-]        
    CALL getin_p('PREF_SOIL_VEG',pref_soil_veg)

  END SUBROUTINE config_pft_parameters


!! ================================================================================================================================
!! SUBROUTINE   : config_sechiba_pft_parameters
!!
!>\BRIEF        This subroutine will read the imposed values for the sechiba pft
!! parameters. It is not called if IMPOSE_PARAM is set to NO. 
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE config_sechiba_pft_parameters()

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.1 Input variables

    !! 0.4 Local variable

    !_ ================================================================================================================================ 

    !
    ! Evapotranspiration -  sechiba
    !

    !Config Key   = RSTRUCT_CONST
    !Config Desc  = Structural resistance 
    !Config if    = OK_SECHIBA
    !Config Def   = 0.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0,  2.5,  2.0,  2.0,  2.0
    !Config Help  =
    !Config Units = [s/m]
    CALL getin_p('RSTRUCT_CONST',rstruct_const)

    !Config Key   = KZERO
    !Config Desc  = A vegetation dependent constant used in the calculation of the surface resistance.
    !Config if    = OK_SECHIBA
    !Config Def   = 0.0, 12.E-5, 12.E-5, 12.e-5, 12.e-5, 25.e-5, 12.e-5,25.e-5, 25.e-5, 30.e-5, 30.e-5, 30.e-5, 30.e-5 
    !Config Help  =
    !Config Units = [kg/m^2/s]
    CALL getin_p('KZERO',kzero)

    !Config Key   = RVEG_PFT
    !Config Desc  = Artificial parameter to increase or decrease canopy resistance.
    !Config if    = OK_SECHIBA
    !Config Def   = 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.
    !Config Help  = This parameter is set by PFT.
    !Config Units = [-]
    CALL getin_p('RVEG_PFT',rveg_pft)    

    !
    ! Water-hydrology - sechiba
    !

    !Config Key   = WMAX_VEG
    !Config Desc  = Maximum field capacity for each of the vegetations (Temporary): max quantity of water
    !Config if    = OK_SECHIBA
    !Config Def   = 150., 150., 150., 150., 150., 150., 150.,150., 150., 150., 150., 150., 150.
    !Config Help  =
    !Config Units = [kg/m^3]
    CALL getin_p('WMAX_VEG',wmax_veg)

    !Config Key   = PERCENT_THROUGHFALL_PFT
    !Config Desc  = Percent by PFT of precip that is not intercepted by the canopy. Default value depend on run mode.
    !Config if    = OK_SECHIBA
    !Config Def   = Case offline+CWRR [0. 0. 0....] else [30. 30. 30.....]
    !Config Help  = During one rainfall event, PERCENT_THROUGHFALL_PFT% of the incident rainfall
    !Config         will get directly to the ground without being intercepted, for each PFT.
    !Config Units = [%]
    CALL getin_p('PERCENT_THROUGHFALL_PFT',throughfall_by_pft)
    throughfall_by_pft(:) = throughfall_by_pft(:) / 100. 


    !
    ! Albedo - sechiba
    !

    !Config Key   = SNOWA_AGED_VIS
    !Config Desc  = Minimum snow albedo value for each vegetation type after aging (dirty old snow), visible albedo
    !Config if    = OK_SECHIBA
    !Config Def   = 0.5, 0., 0., 0.15, 0.14, 0.14, 0.15, 0.14, 0.22, 0.35, 0.35, 0.35, 0.35
    !Config Help  = Values are from the Thesis of S. Chalita (1992), optimized on 04/07/2016
    !Config Units = [-]
    CALL getin_p('SNOWA_AGED_VIS',snowa_aged_vis)

    !Config Key   = SNOWA_AGED_NIR
    !Config Desc  = Minimum snow albedo value for each vegetation type after aging (dirty old snow), near infrared albedo
    !Config if    = OK_SECHIBA
    !Config Def   = 0.35, 0., 0., 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.18, 0.18, 0.18, 0.18
    !Config Help  = Values are from the Thesis of S. Chalita (1992)
    !Config Units = [-]
    CALL getin_p('SNOWA_AGED_NIR',snowa_aged_nir)

    !Config Key   = SNOWA_DEC_VIS
    !Config Desc  = Decay rate of snow albedo value for each vegetation type as it will be used in condveg_snow, visible albedo
    !Config if    = OK_SECHIBA
    !Config Def   = 0.45, 0., 0., 0.1, 0.06, 0.11, 0.10, 0.11, 0.18, 0.60, 0.60, 0.60, 0.60
    !Config Help  = Values are from the Thesis of S. Chalita (1992), optimized on 04/07/2016
    !Config Units = [-]
    CALL getin_p('SNOWA_DEC_VIS',snowa_dec_vis)

    !Config Key   = SNOWA_DEC_NIR
    !Config Desc  = Decay rate of snow albedo value for each vegetation type as it will be used in condveg_snow, near infrared albedo
    !Config if    = OK_SECHIBA
    !Config Def   = 0.45, 0.,  0., 0.06, 0.06, 0.11, 0.06, 0.11, 0.11, 0.52 ,0.52, 0.52, 0.52
    !Config Help  = Values are from the Thesis of S. Chalita (1992)
    !Config Units = [-]
    CALL getin_p('SNOWA_DEC_NIR',snowa_dec_nir)

    !Config Key   = ALB_LEAF_VIS
    !Config Desc  = leaf albedo of vegetation type, visible albedo
    !Config if    = OK_SECHIBA
    !Config Def   = .0, .0397, .0474, .0386, .0484, .0411, .041, .0541, .0435, .0524, .0508, .0509, .0606
    !Config Help  = optimized on 04/07/2016
    !Config Units = [-]
    CALL getin_p('ALB_LEAF_VIS',alb_leaf_vis)

    !Config Key   = ALB_LEAF_NIR
    !Config Desc  = leaf albedo of vegetation type, near infrared albedo
    !Config if    = OK_SECHIBA
    !Config Def   = .0, .227, .214, .193, .208, .244, .177, .218, .213, .252, .265, .272, .244
    !Config Help  = optimized on 04/07/2016
    !Config Units = [-]
    CALL getin_p('ALB_LEAF_NIR',alb_leaf_nir)

    IF ( ok_bvoc ) THEN
       !
       ! BVOC
       !

       !Config Key   = ISO_ACTIVITY
       !Config Desc  = Biogenic activity for each age class : isoprene
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.5, 1.5, 1.5, 0.5
       !Config Help  =
       !Config Units = [-]
       CALL getin_p('ISO_ACTIVITY',iso_activity)

       !Config Key   = METHANOL_ACTIVITY
       !Config Desc  = Isoprene emission factor for each age class : methanol
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 1., 1., 0.5, 0.5
       !Config Help  =
       !Config Units = [-]
       CALL getin_p('METHANOL_ACTIVITY',methanol_activity)

       !Config Key   = EM_FACTOR_ISOPRENE
       !Config Desc  = Isoprene emission factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 24., 24., 8., 16., 45., 8., 18., 0.5, 12., 18., 5., 5.
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_ISOPRENE',em_factor_isoprene)

       !Config Key   = EM_FACTOR_MONOTERPENE
       !Config Desc  = Monoterpene emission factor 
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 2.0, 2.0, 1.8, 1.4, 1.6, 1.8, 1.4, 1.8, 0.8, 0.8,  0.22, 0.22
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_MONOTERPENE',em_factor_monoterpene)

       !Config Key   = C_LDF_MONO 
       !Config Desc  = Monoterpenes fraction dependancy to light
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.6
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_LDF_MONO',LDF_mono)

       !Config Key   = C_LDF_SESQ 
       !Config Desc  = Sesquiterpenes fraction dependancy to light
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.5
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_LDF_SESQ',LDF_sesq)

       !Config Key   = C_LDF_METH 
       !Config Desc  = Methanol fraction dependancy to light
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.8
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_LDF_METH',LDF_meth)

       !Config Key   = C_LDF_ACET 
       !Config Desc  = Acetone fraction dependancy to light
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.2
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_LDF_ACET',LDF_acet)

       !Config Key   = EM_FACTOR_APINENE 
       !Config Desc  = Alfa pinene  emission factor 
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 1.35, 1.35, 0.85, 0.95, 0.75, 0.85, 0.60, 1.98, 0.30, 0.30, 0.09, 0.09
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_APINENE',em_factor_apinene)

       !Config Key   = EM_FACTOR_BPINENE
       !Config Desc  = Beta pinene  emission factor
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 0.30, 0.30, 0.35, 0.25, 0.20, 0.35, 0.12, 0.45, 0.16, 0.12, 0.05, 0.05
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_BPINENE',em_factor_bpinene)

       !Config Key   = EM_FACTOR_LIMONENE
       !Config Desc  = Limonene  emission factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.25, 0.25, 0.20, 0.25, 0.14, 0.20, 0.135, 0.11, 0.19, 0.42, 0.03, 0.03
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_LIMONENE',em_factor_limonene)

       !Config Key   = EM_FACTOR_MYRCENE
       !Config Desc  = Myrcene  emission factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.20, 0.20, 0.12, 0.11, 0.065, 0.12, 0.036, 0.075, 0.08,  0.085, 0.015, 0.015
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_MYRCENE',em_factor_myrcene)

       !Config Key   = EM_FACTOR_SABINENE
       !Config Desc  = Sabinene  emission factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.20, 0.20, 0.12, 0.17, 0.70, 0.12, 0.50, 0.09, 0.085, 0.075, 0.02, 0.02
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_SABINENE',em_factor_sabinene)

       !Config Key   = EM_FACTOR_CAMPHENE 
       !Config Desc  = Camphene  emission factor 
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.15, 0.15, 0.10, 0.10, 0.01, 0.10, 0.01, 0.07, 0.07, 0.08, 0.01, 0.01
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_CAMPHENE',em_factor_camphene)

       !Config Key   = EM_FACTOR_3CARENE 
       !Config Desc  = 3-Carene  emission factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.13, 0.13, 0.42, 0.02, 0.055, 0.42,0.025, 0.125, 0.085, 0.085, 0.065, 0.065
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_3CARENE',em_factor_3carene)

       !Config Key   = EM_FACTOR_TBOCIMENE
       !Config Desc  = T-beta-ocimene  emission factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.25, 0.25, 0.13, 0.09, 0.26, 0.13, 0.20, 0.085, 0.18, 0.18, 0.01, 0.01
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_TBOCIMENE', em_factor_tbocimene)

       !Config Key   = EM_FACTOR_OTHERMONOT
       !Config Desc  = Other monoterpenes  emission factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.17, 0.17, 0.11, 0.11, 0.125, 0.11, 0.274, 0.01, 0.15, 0.155, 0.035, 0.035
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_OTHERMONOT',em_factor_othermonot)

       !Config Key   = EM_FACTOR_SESQUITERP 
       !Config Desc  = Sesquiterpenes  emission factor 
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.45, 0.45, 0.13, 0.3, 0.36, 0.15, 0.3, 0.25, 0.6, 0.6, 0.08, 0.08
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_SESQUITERP',em_factor_sesquiterp)



       !Config Key   = C_BETA_MONO 
       !Config Desc  = Monoterpenes temperature dependency coefficient
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.1
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_BETA_MONO',beta_mono)

       !Config Key   = C_BETA_SESQ 
       !Config Desc  = Sesquiterpenes temperature dependency coefficient
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.17
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_BETA_SESQ',beta_sesq)

       !Config Key   = C_BETA_METH 
       !Config Desc  = Methanol temperature dependency coefficient
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.08
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_BETA_METH',beta_meth)

       !Config Key   = C_BETA_ACET 
       !Config Desc  = Acetone temperature dependency coefficient
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.1
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_BETA_ACET',beta_acet)

       !Config Key   = C_BETA_OXYVOC 
       !Config Desc  = Other oxygenated BVOC temperature dependency coefficient
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.13
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_BETA_OXYVOC',beta_oxyVOC)

       !Config Key   = EM_FACTOR_ORVOC
       !Config Desc  = ORVOC emissions factor 
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5
       !Config Help  =
       !Config Units = [ugC/g/h]  
       CALL getin_p('EM_FACTOR_ORVOC',em_factor_ORVOC)

       !Config Key   = EM_FACTOR_OVOC
       !Config Desc  = OVOC emissions factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5
       !Config Help  =
       !Config Units = [ugC/g/h]        
       CALL getin_p('EM_FACTOR_OVOC',em_factor_OVOC)

       !Config Key   = EM_FACTOR_MBO
       !Config Desc  = MBO emissions factor 
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 2.e-5, 2.e-5, 1.4, 2.e-5, 2.e-5, 0.14, 2.e-5, 2.e-5, 2.e-5, 2.e-5, 2.e-5, 2.e-5
       !Config Help  =
       !Config Units = [ugC/g/h]  
       CALL getin_p('EM_FACTOR_MBO',em_factor_MBO)

       !Config Key   = EM_FACTOR_METHANOL
       !Config Desc  = Methanol emissions factor 
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 0.8, 0.8, 1.8, 0.9, 1.9, 1.8, 1.8, 1.8, 0.7, 0.9, 2., 2.
       !Config Help  =
       !Config Units = [ugC/g/h]  
       CALL getin_p('EM_FACTOR_METHANOL',em_factor_methanol)

       !Config Key   = EM_FACTOR_ACETONE
       !Config Desc  = Acetone emissions factor
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 0.25, 0.25, 0.3, 0.2, 0.33, 0.3, 0.25, 0.25, 0.2, 0.2, 0.08, 0.08
       !Config Help  =
       !Config Units = [ugC/g/h]     
       CALL getin_p('EM_FACTOR_ACETONE',em_factor_acetone)

       !Config Key   = EM_FACTOR_ACETAL
       !Config Desc  = Acetaldehyde emissions factor 
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.2, 0.2, 0.2, 0.2, 0.25, 0.25, 0.16, 0.16, 0.12, 0.12, 0.035, 0.02
       !Config Help  =
       !Config Units = [ugC/g/h]  
       CALL getin_p('EM_FACTOR_ACETAL',em_factor_acetal)

       !Config Key   = EM_FACTOR_FORMAL
       !Config Desc  = Formaldehyde emissions factor
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 0.04, 0.04, 0.08, 0.04, 0.04, 0.04, 0.04, 0.04, 0.025, 0.025, 0.013, 0.013
       !Config Help  = 
       !Config Units = [ugC/g/h]  
       CALL getin_p('EM_FACTOR_FORMAL',em_factor_formal)

       !Config Key   = EM_FACTOR_ACETIC
       !Config Desc  = Acetic Acid emissions factor
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 0.025, 0.025,0.025,0.022,0.08,0.025,0.022,0.013,0.012,0.012,0.008,0.008
       !Config Help  =
       !Config Units = [ugC/g/h]  
       CALL getin_p('EM_FACTOR_ACETIC',em_factor_acetic)

       !Config Key   = EM_FACTOR_FORMIC
       !Config Desc  = Formic Acid emissions factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.015, 0.015, 0.02, 0.02, 0.025, 0.025, 0.015, 0.015,0.010,0.010,0.008,0.008
       !Config Help  =
       !Config Units = [ugC/g/h]  
       CALL getin_p('EM_FACTOR_FORMIC',em_factor_formic)

       !Config Key   = EM_FACTOR_NO_WET
       !Config Desc  = NOx emissions factor wet soil emissions and exponential dependancy factor 
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 2.6, 0.06, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.36, 0.36, 0.36, 0.36
       !Config Help  =
       !Config Units = [ngN/m^2/s]
       CALL getin_p('EM_FACTOR_NO_WET',em_factor_no_wet)

       !Config Key   = EM_FACTOR_NO_DRY
       !Config Desc  = NOx emissions factor dry soil emissions and exponential dependancy factor 
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 8.60, 0.40, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 2.65, 2.65, 2.65, 2.65
       !Config Help  =
       !Config Units = [ngN/m^2/s] 
       CALL getin_p('EM_FACTOR_NO_DRY',em_factor_no_dry)

       !Config Key   = LARCH
       !Config Desc  = Larcher 1991 SAI/LAI ratio
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 0.015, 0.015, 0.003, 0.005, 0.005, 0.003, 0.005, 0.003, 0.005, 0.005, 0.008, 0.008
       !Config Help  =
       !Config Units = [-]  
       CALL getin_p('LARCH',Larch)

    ENDIF ! (ok_bvoc)

  END SUBROUTINE config_sechiba_pft_parameters


!! ================================================================================================================================
!! SUBROUTINE   : config_stomate_pft_parameters 
!!
!>\BRIEF         This subroutine will read the imposed values for the stomate pft
!! parameters. It is not called if IMPOSE_PARAM is set to NO.
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE config_stomate_pft_parameters

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.4 Local variable
   INTEGER(i_std)               :: ivma,ivm,j,k!! index

    !_ ================================================================================================================================

    !
    ! Vegetation structure
    !
!gmjc
     !Config  Key  = IS_GRASSLAND_MANAG
     !Config  Desc = Is the vegetation type a managed grassland ?
     !Config  if  = OK_STOMATE
     !Config  Def  = n, n, n, n, n, n, n, n, n, y, n, n, n
     !Config  Help =
     !Config  Units = NONE
     CALL getin_p('GRM_IS_GRASSLAND_MANAG',is_grassland_manag)
     WRITE(numout,*) 'GRM_IS_GRASSLAND_MANAG',is_grassland_manag
     ! set all managed grassland to unnatural, consistency issue
!mictgmjc set managed grassland as natural too for allowing fire on pasture    
!     DO j = 1,nvm
!       IF (is_grassland_manag(j) == .TRUE.) THEN
!         natural(j) = .FALSE.
!       ENDIF
!     ENDDO
!mictgmjc set managed grassland as pasture to exclude dynamic vegetaion (for DGVM only)
     DO j = 1,nvm
       IF (is_grassland_manag(j) == .TRUE.) THEN
         pasture(j) = .TRUE.
       ENDIF
     ENDDO
     ![chaoyue] we have to warn temporarily if use_age_class is .TRUE.
     !but none of grassland is managed. 
     IF (use_age_class .AND. ALL(.NOT.(is_grassland_manag))) THEN
       CALL ipslerr_p(3, 'config_pft_parameters', &
                         'Age classes are used but none of the ', &
                         'grasslands are managed! looks wierd, please confirm and uncoment this', &
                         'line if this is really what you want.' )
     ENDIF   

     !Config  Key  = IS_GRASSLAND_CUT
     !Config  Desc = Is the vegetation type a cut grassland for management
     !adaptation ?
     !Config  if  = OK_STOMATE
     !Config  Def  = n, n, n, n, n, n, n, n, n, n, n, n, n
     !Config  Help =
     !Config  Units = NONE
     CALL getin_p('GRM_IS_GRASSLAND_CUT',is_grassland_cut)
     WRITE(numout,*) 'GRM_IS_GRASSLAND_CUT',is_grassland_cut
     !Config  Key  = IS_GRASSLAND_GRAZED
     !Config  Desc = Is the vegetation type a grazed grassland for management
     !adaptation ?
     !Config  if  = OK_STOMATE
     !Config  Def  = n, n, n, n, n, n, n, n, n, n, n, n, n
     !Config  Help =
     !Config  Units = NONE
     CALL getin_p('GRM_IS_GRASSLAND_GRAZED',is_grassland_grazed)
     WRITE(numout,*) 'GRM_IS_GRASSLAND_GRAZED',is_grassland_grazed
     !Config  Key  = MANAGEMENT_INTENSITY
     !Config  Desc = management intensity for grassland management
     !adaptation ?
     !Config  if  = OK_STOMATE
     !Config  Def  = n, n, n, n, n, n, n, n, n, n, n, n, n
     !Config  Help =
     !Config  Units = NONE
     CALL getin_p('GRM_MANAGE_INTENSITY',management_intensity)
     WRITE(numout,*) 'GRM_MANAGE_INTENSITY',management_intensity
     !Config  Key  = NB_YEAR_MANAGEMENT
     !Config  Desc = number of years for grassland management
     !adaptation ?
     !Config  if  = OK_STOMATE
     !Config  Def  = n, n, n, n, n, n, n, n, n, n, n, n, n
     !Config  Help =
     !Config  Units = NONE
     CALL getin_p('GRM_NB_YEAR_MANAGEMENT',nb_year_management)
     !Config  Key  = MANAGEMENT_START
     !Config  Desc = start time of grassland management
     !adaptation ?
     !Config  if  = OK_STOMATE
     !Config  Def  = n, n, n, n, n, n, n, n, n, n, n, n, n
     !Config  Help =
     !Config  Units = NONE
     CALL getin_p('GRM_MANAGEMENT_START',management_start)
     !Config  Key  = DEPOSITION_START
     !Config  Desc = start time of N depostion for grassland management
     !adaptation ?
     !Config  if  = OK_STOMATE
     !Config  Def  = n, n, n, n, n, n, n, n, n, n, n, n, n
     !Config  Help =
     !Config  Units = NONE
     CALL getin_p('GRM_DEPOSITION_START',deposition_start)
!end gmjc

     CALL getin_p('GLUC_IS_BIOE1',is_bioe1)

     !!! CROP module
      !
      !Config Key   = SP_CODEPHOT
      !Config Desc  = whether or not sensitive to photoperiod
      !Config if    = OK_STOMATE
      !Config Def   = undef_int,undef_int, undef_int, undef_int, undef_int, undef_int, undef_int, undef_int, undef_int, undef_int, undef_int, 1, 1
      !Config Help  =
      !Config Units = [C]
      CALL getin_p('SP_CODEPHOT',SP_codephot)
     

      !
      !Config Key   = SP_iplt0
      !Config Desc  = sowing date
      !Config if    = OK_STOMATE
      !Config Def   = undef_int,  undef_int,undef_int, undef_int, undef_int, undef_int, undef_int, undef_int, undef_int, undef_int, undef_int, 292, 117
      !Config Help  =
      !Config Units = [C]
      CALL getin_p('SP_IPLT0',SP_iplt0)
      IF (cyc_rot_max .GT. 1) THEN
          CALL getin_p('SP_IPLT1',SP_iplt1)
      ENDIF
      IF (cyc_rot_max .GT. 2) THEN
          CALL getin_p('SP_IPLT2',SP_iplt2)
      ENDIF

      CALL getin_p('CODELAINET',SP_codlainet)
      CALL getin_p('STPLTGER',SP_stpltger)
      CALL getin_p('STADE0',SP_stade0)
      CALL getin_p('DLAIMAX',SP_dlaimax)
!      CALL getin_p('INNSENES',SP_innsenes)
      CALL getin_p('CODEHYPO',SP_codehypo)
      CALL getin_p('LAIPLANTULE',SP_laiplantule)
      CALL getin_p('INNLAI',SP_innlai)
      CALL getin_P('DURVIEF',SP_durvieF)
      CALL getin_p('VLAIMAX',SP_vlaimax)
      CALL getin_p('STLEVAMF',SP_stlevamf)
      CALL getin_p('STDRPMAT',SP_stdrpmat)
      CALL getin_p('STLEVDRP',SP_stlevdrp)
      CALL getin_p('STAMFLAX',SP_stamflax)
      CALL getin_p('NUMAGEBOX',SP_nbox)
      CALL getin_p('LAI0',SP_lai0)
      CALL getin_p('TDMAX',SP_tdmax)
      CALL getin_p('TDMIN',SP_tdmin)
      CALL getin_p('TCXSTOP',SP_tcxstop)
      CALL getin_p('TCMAX',SP_tcmax)
      CALL getin_p('TCMIN',SP_tcmin)
      CALL getin_p('NEFFMAX',SP_neffmax)
      CALL getin_p('NSATRAT',SP_nsatrat)
      CALL getin_p('CODEIR',SP_codeir)
      CALL getin_p('VITIRCARB',SP_vitircarb)
      CALL getin_p('VITIRCARBT',SP_vitircarbT)
      CALL getin_p('SWFACMIN',SP_swfacmin)
      CALL getin_p('IRMAX',SP_irmax)
      CALL getin_p('REPRACMAX',SP_repracpermax)
      CALL getin_p('REPRACMIN',SP_repracpermin)
      CALL getin_p('TMINREMP',SP_tminremp)
      CALL getin_p('TMAXREMP',SP_tmaxremp)
      CALL getin_p('NBJGRAIN',SP_nbjgrain)

      CALL getin_p('DENSITESEM',SP_densitesem)
      CALL getin_p('SLAMAX',SP_slamax)
      CALL getin_p('STLAXSEN',SP_stlaxsen)
      CALL getin_p('STSENLAN',SP_stsenlan)
      CALL getin_p('ZRACPLANTULE',SP_zracplantule)
      CALL getin_p('TGMIN',SP_tgmin)
      CALL getin_p('BDENS',SP_bdens)
       

      !! for dynamic nitrogen processes

      !
      !
      !Config Key   = DY_INN
      !Config Desc  = whether or not we use the dynamic nitrogen processes
      !Config if    = OK_STOMATE
      !Config Def   = .false., .false., .false., .false., .false., .false., .false., .false., .false., .false., .false.,  .false., .false.
      !Config Help  =
      !Config Units = logic
      CALL getin_p('DY_INN',SP_DY_INN)


      !! for dynamic nitrogen processes

      !
      !
      !Config Key   = SP_AVENFERT
      !Config Desc  = the average nitrogen fertilization
      !Config if    = OK_STOMATE
      !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 150.0, 100.0, 100.0
      !Config Help  =
      !Config Units = kg N ha-1
      CALL getin_p('SP_AVENFERT',SP_avenfert)


!!!!! end crop parameters, xuhui

    !Config Key   = SLA
    !Config Desc  = specif leaf area 
    !Config if    = OK_STOMATE
    !Config Def   = 1.5E-2, 1.53E-2, 2.6E-2, 9.26E-3, 2E-2, 2.6E-2, 9.26E-3, 2.6E-2, 1.9E-2, 2.6E-2, 2.6E-2, 2.6E-2, 2.6E-2
    !Config Help  =
    !Config Units = [m^2/gC]
    CALL getin_p('SLA',sla)

    !Config Key   = AVAILABILITY_FACT 
    !Config Desc  = Calculate dynamic mortality in lpj_gap, pft dependent parameter
    !Config If    = OK_STOMATE 
    !Config Def   = undef, 0.14, 0.14, 0.10, 0.10, 0.10, 0.05, 0.05, 0.05, undef, undef, undef, undef 
    !Config Help  = 
    !Config Units = [-]   
    CALL getin_p('AVAILABILITY_FACT',availability_fact)

    !
    ! Allocation - stomate
    !
    !
    !Config Key   = R0 
    !Config Desc  = Standard root allocation 
    !Config If    = OK_STOMATE 
    !Config Def   = undef, .30, .30, .30, .30, .30, .30, .30, .30, .30, .30, .30, .30
    !Config Help  = 
    !Config Units = [-]    
    CALL getin_p('R0',R0)

    !Config Key   = S0 
    !Config Desc  = Standard sapwood allocation 
    !Config If    = OK_STOMATE 
    !Config Def   = undef, .25, .25, .30, .30, .30, .30, .30, .30, .30, .30, .30, .30
    !Config Help  = 
    !Config Units = [-]    
    CALL getin_p('S0',S0)

    !
    ! Respiration - stomate
    !

    !Config Key   = FRAC_GROWTHRESP
    !Config Desc  = fraction of GPP which is lost as growth respiration
    !Config if    = OK_STOMATE
    !Config Def   = undef, .28, .28, .28, .28, .28, .28, .28, .28, .28, .28, .28, .28 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('FRAC_GROWTHRESP',frac_growthresp) 

    !Config Key   = MAINT_RESP_SLOPE_C
    !Config Desc  = slope of maintenance respiration coefficient (1/K), constant c of aT^2+bT+c , tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef, .20, .20, .16, .16, .16, .16, .16, .16, .16, .12, .16, .12 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('MAINT_RESP_SLOPE_C',maint_resp_slope_c) 

    !Config Key   = MAINT_RESP_SLOPE_B
    !Config Desc  = slope of maintenance respiration coefficient (1/K), constant b of aT^2+bT+c , tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef, .0, .0, .0, .0, .0, .0, .0, .0, -.00133, .0, -.00133, .0 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('MAINT_RESP_SLOPE_B',maint_resp_slope_b)

    !Config Key   = MAINT_RESP_SLOPE_A
    !Config Desc  = slope of maintenance respiration coefficient (1/K), constant a of aT^2+bT+c , tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef, .0, .0, .0, .0, .0, .0, .0, .0, .0, .0, .0, .0    
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('MAINT_RESP_SLOPE_A',maint_resp_slope_a)

    !Config Key   = CM_ZERO_LEAF
    !Config Desc  = maintenance respiration coefficient at 0 deg C, for leaves, tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef, 2.35E-3, 2.62E-3, 1.01E-3, 2.35E-3, 2.62E-3, 1.01E-3,2.62E-3, 2.05E-3, 2.62E-3, 2.62E-3, 2.62E-3, 2.62E-3
    !Config Help  =
    !Config Units = [g/g/day]
    CALL getin_p('CM_ZERO_LEAF',cm_zero_leaf)

    !Config Key   = CM_ZERO_SAPABOVE
    !Config Desc  = maintenance respiration coefficient at 0 deg C,for sapwood above, tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4
    !Config Help  =
    !Config Units = [g/g/day]
    CALL getin_p('CM_ZERO_SAPABOVE',cm_zero_sapabove)

    !Config Key   = CM_ZERO_SAPBELOW
    !Config Desc  = maintenance respiration coefficient at 0 deg C, for sapwood below, tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4 
    !Config Help  =
    !Config Units = [g/g/day]
    CALL getin_p('CM_ZERO_SAPBELOW',cm_zero_sapbelow)

    !Config Key   = CM_ZERO_HEARTABOVE
    !Config Desc  = maintenance respiration coefficient at 0 deg C, for heartwood above, tabulated
    !Config if    = OK_STOMATE 
    !Config Def   = undef, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. 
    !Config Help  =
    !Config Units = [g/g/day]
    CALL getin_p('CM_ZERO_HEARTABOVE',cm_zero_heartabove)

    !Config Key   = CM_ZERO_HEARTBELOW
    !Config Desc  = maintenance respiration coefficient at 0 deg C,for heartwood below, tabulated
    !Config if    = OK_STOMATE 
    !Config Def   = undef, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. 
    !Config Help  =
    !Config Units = [g/g/day] 
    CALL getin_p('CM_ZERO_HEARTBELOW',cm_zero_heartbelow)

    !Config Key   = CM_ZERO_ROOT
    !Config Desc  = maintenance respiration coefficient at 0 deg C, for roots, tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef,1.67E-3, 1.67E-3, 1.67E-3, 1.67E-3, 1.67E-3, 1.67E-3,1.67E-3, 1.67E-3, 1.67E-3, 1.67E-3, 1.67E-3, 1.67E-3
    !Config Help  =
    !Config Units = [g/g/day] 
    CALL getin_p('CM_ZERO_ROOT',cm_zero_root)

    !Config Key   = CM_ZERO_FRUIT
    !Config Desc  = maintenance respiration coefficient at 0 deg C, for fruits, tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4,1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4    
    !Config Help  =
    !Config Units = [g/g/day] 
    CALL getin_p('CM_ZERO_FRUIT',cm_zero_fruit)

    !Config Key   = CM_ZERO_CARBRES
    !Config Desc  = maintenance respiration coefficient at 0 deg C, for carbohydrate reserve, tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4,1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4
    !Config Help  =
    !Config Units = [g/g/day] 
    CALL getin_p('CM_ZERO_CARBRES',cm_zero_carbres)

    !
    ! Fire - stomate
    !

    !Config Key   = FLAM
    !Config Desc  = flamability: critical fraction of water holding capacity
    !Config if    = OK_STOMATE
    !Config Def   = undef, .15, .25, .25, .25, .25, .25, .25, .25, .25, .25, .35, .35
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('FLAM',flam)

    !Config Key   = RESIST
    !Config Desc  = fire resistance
    !Config if    = OK_STOMATE
    !Config Def   = undef, .95, .90, .12, .50, .12, .12, .12, .12, .0, .0, .0, .0 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('RESIST',resist)

    !
    ! Flux - LUC
    !

    !Config Key   = COEFF_LCCHANGE_1
    !Config Desc  = Coeff of biomass export for the year
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.897, 0.897, 0.597, 0.597, 0.597, 0.597, 0.597, 0.597, 0.597, 0.597, 0.597, 0.597 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('COEFF_LCCHANGE_1',coeff_lcchange_1)

    !Config Key   = COEFF_LCCHANGE_10
    !Config Desc  = Coeff of biomass export for the decade
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.103, 0.103, 0.299, 0.299, 0.299, 0.299, 0.299, 0.299, 0.299, 0.403, 0.299, 0.403
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('COEFF_LCCHANGE_10',coeff_lcchange_10)

    !Config Key   = COEFF_LCCHANGE_100
    !Config Desc  = Coeff of biomass export for the century
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0., 0., 0.104, 0.104, 0.104, 0.104, 0.104, 0.104, 0.104, 0., 0.104, 0.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('COEFF_LCCHANGE_100',coeff_lcchange_100)
!gmjc
    !Config  Key  =
     !Config  Desc = minimum gdd to allow senescence of crops
     !Config  if  = OK_STOMATE
     !Config  Def  =  ! maximum specific leaf area (m**2/gC)
     !Config  Help =
     !Config  Units = Celsius degrees [C]
     CALL getin_p('SLA_MAX',sla_max)
     !
     !Config  Key  = SLA_MIN
     !Config  Desc = minimum specific leaf area (m**2/gC)
     !Config  if  = OK_STOMATE
     !Config  Def  =
     !Config  Help =
     !Config  Units = Celsius degrees [C]
     CALL getin_p('SLA_MIN',sla_min)
!end gmjc
    !
    ! Phenology
    !

    !Config Key   = LAI_MAX_TO_HAPPY
    !Config Desc  = threshold of LAI below which plant uses carbohydrate reserves
    !Config if    = OK_STOMATE
    !Config Def   = undef, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('LAI_MAX_TO_HAPPY',lai_max_to_happy) 

    !Config Key   = LAI_MAX
    !Config Desc  = maximum LAI, PFT-specific
    !Config if    = OK_STOMATE
    !Config Def   = undef, 7., 7., 5., 5., 5., 4.5, 4.5, 3.0, 2.5, 2.5, 5.,5. 
    !Config Help  =
    !Config Units = [m^2/m^2]
    CALL getin_p('LAI_MAX',lai_max)

    !Config Key   = PHENO_TYPE
    !Config Desc  = type of phenology, 0=bare ground 1=evergreen,  2=summergreen,  3=raingreen,  4=perennial
    !Config if    = OK_STOMATE
    !Config Def   = 0, 1, 3, 1, 1, 2, 1, 2, 2, 4, 4, 2, 3
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('PHENO_TYPE',pheno_type)

    !
    ! Phenology : Leaf Onset
    !

    !Config Key   = PHENO_GDD_CRIT_C
    !Config Desc  = critical gdd, tabulated (C), constant c of aT^2+bT+c
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, 270., 400., 125., 400.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('PHENO_GDD_CRIT_C',pheno_gdd_crit_c)

    !Config Key   = PHENO_GDD_CRIT_B
    !Config Desc  = critical gdd, tabulated (C), constant b of aT^2+bT+c
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef,undef, undef, 6.25, 0., 0., 0.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('PHENO_GDD_CRIT_B',pheno_gdd_crit_b)

    !Config Key   = PHENO_GDD_CRIT_A
    !Config Desc  = critical gdd, tabulated (C), constant a of aT^2+bT+c
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, 0.03125,  0., 0., 0.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('PHENO_GDD_CRIT_A',pheno_gdd_crit_a)

    !Config Key   = PHENO_MOIGDD_T_CRIT
    !Config Desc  = Average temperature threashold for C4 grass used in pheno_moigdd
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 22.0, undef, undef
    !Config Help  =
    !Config Units = [C]
    CALL getin_p('PHENO_MOIGDD_T_CRIT',pheno_moigdd_t_crit)

    !Config Key   = NGD_CRIT
    !Config Desc  = critical ngd, tabulated. Threshold -5 degrees
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, 0., undef, undef, undef, undef, undef
    !Config Help  = NGD : Number of Growing Days.
    !Config Units = [days]
    CALL getin_p('NGD_CRIT',ngd_crit)

    !Config Key   = NCDGDD_TEMP
    !Config Desc  = critical temperature for the ncd vs. gdd function in phenology
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, 5., undef, 0., undef, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [C] 
    CALL getin_p('NCDGDD_TEMP',ncdgdd_temp)

    !Config Key   = HUM_FRAC
    !Config Desc  = critical humidity (relative to min/max) for phenology
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, .5, undef, undef, undef, undef, undef,  undef, .5, .5, .5,.5     
    !Config Help  =
    !Config Units = [%]
    CALL getin_p('HUM_FRAC',hum_frac)

    !Config Key   = HUM_MIN_TIME
    !Config Desc  = minimum time elapsed since moisture minimum
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, 50., undef, undef, undef, undef, undef, undef, 36., 35., 75., 75.
    !Config Help  =
    !Config Units = [days]
    CALL getin_p('HUM_MIN_TIME',hum_min_time)

    !Config Key   = TAU_SAP
    !Config Desc  = sapwood -> heartwood conversion time
    !Config if    = OK_STOMATE
    !Config Def   = undef, 730., 730., 730., 730., 730., 730., 730., 730., undef, undef, undef, undef
    !Config Help  =
    !Config Units = [days]
    CALL getin_p('TAU_SAP',tau_sap)

    !Config Key   = TAU_LEAFINIT
    !Config Desc  = time to attain the initial foliage using the carbohydrate reserve
    !Config if    = OK_STOMATE
    !Config Def   = undef, 10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 10.
    !Config Help  =
    !Config Units = [days]
    CALL getin_p('TAU_LEAFINIT',tau_leafinit) 

    !Config Key   = TAU_FRUIT
    !Config Desc  = fruit lifetime
    !Config if    = OK_STOMATE
    !Config Def   = undef, 90., 90., 90., 90., 90., 90., 90., 90., undef, undef, undef, undef
    !Config Help  =
    !Config Units = [days]
    CALL getin_p('TAU_FRUIT',tau_fruit)

    !Config Key   = ECUREUIL
    !Config Desc  = fraction of primary leaf and root allocation put into reserve
    !Config if    = OK_STOMATE
    !Config Def   = undef, .0, 1., .0, .0, 1., .0, 1., 1., 1., 1., 1., 1.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('ECUREUIL',ecureuil)

    !Config Key   = ALLOC_MIN
    !Config Desc  = minimum allocation above/below = f(age) - 30/01/04 NV/JO/PF
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, undef, undef, undef, undef 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('ALLOC_MIN',alloc_min)

    !Config Key   = ALLOC_MAX
    !Config Desc  = maximum allocation above/below = f(age) - 30/01/04 NV/JO/PF
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('ALLOC_MAX',alloc_max)

    !Config Key   = DEMI_ALLOC 
    !Config Desc  = mean allocation above/below = f(age) - 30/01/04 NV/JO/PF
    !Config if    = OK_STOMATE
    !Config Def   = undef, 5., 5., 5., 5., 5., 5., 5., 5., undef, undef, undef, undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('DEMI_ALLOC',demi_alloc)

    !Config Key   = LEAFLIFE_TAB
    !Config Desc  = leaf longevity
    !Config if    = OK_STOMATE
    !Config Def   = undef, .5, 2., .33, 1., 2., .33, 2., 2., 2., 2., 2., 2. 
    !Config Help  =
    !Config Units = [years]
    CALL getin_p('LEAFLIFE_TAB',leaflife_tab)

    !
    ! Phenology : Senescence
    !
    !
    !Config Key   = LEAFFALL
    !Config Desc  = length of death of leaves, tabulated 
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, 10., undef, undef, 30., undef, 5., 10., 10., 10., 10., 10. 
    !Config Help  =
    !Config Units = [days]
    CALL getin_p('LEAFFALL',leaffall)

    !Config Key   = LEAFAGECRIT
    !Config Desc  = critical leaf age, tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef, 730., 180., 910., 730., 160., 910., 220., 120., 80., 120., 90., 90.  
    !Config Help  =
    !Config Units = [days]
    CALL getin_p('LEAFAGECRIT',leafagecrit) 

    !Config Key   = SENESCENCE_TYPE
    !Config Desc  = type of senescence, tabulated
    !Config if    = OK_STOMATE
    !Config Def   = none, none, dry, none, none, cold, none, cold, cold, mixed, mixed, mixed, mixed 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('SENESCENCE_TYPE',senescence_type) 

    !Config Key   = SENESCENCE_HUM
    !Config Desc  = critical relative moisture availability for senescence
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, .3, undef, undef, undef, undef, undef, undef, .2, .2, .3, .2 
    !Config Help  =
    !Config Units = [-] 
    CALL getin_p('SENESCENCE_HUM',senescence_hum)

    !Config Key   = NOSENESCENCE_HUM
    !Config Desc  = relative moisture availability above which there is no humidity-related senescence
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, .8, undef, undef, undef, undef, undef, undef, 0.6, .3, .3, .3 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('NOSENESCENCE_HUM',nosenescence_hum) 

    !Config Key   = MAX_TURNOVER_TIME
    !Config Desc  = maximum turnover time for grasse
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef,  80.,  80., 80., 80. 
    !Config Help  =
    !Config Units = [days]
    CALL getin_p('MAX_TURNOVER_TIME',max_turnover_time)

    !Config Key   = MIN_TURNOVER_TIME
    !Config Desc  = minimum turnover time for grasse 
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, 10., 10., 10., 10. 
    !Config Help  =
    !Config Units = [days]
    CALL getin_p('MIN_TURNOVER_TIME',min_turnover_time)

    !Config Key   = MIN_LEAF_AGE_FOR_SENESCENCE
    !Config Desc  = minimum leaf age to allow senescence g
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, 90., undef, undef, 90., undef, 60., 60., 30., 30., 30., 30.
    !Config Help  =
    !Config Units = [days] 
    CALL getin_p('MIN_LEAF_AGE_FOR_SENESCENCE',min_leaf_age_for_senescence)

    !Config Key   = SENESCENCE_TEMP_C
    !Config Desc  = critical temperature for senescence (C), constant c of aT^2+bT+c, tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, 16., undef, 14., 10, 5, 5., 5., 10.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('SENESCENCE_TEMP_C',senescence_temp_c)

    !Config Key   = SENESCENCE_TEMP_B
    !Config Desc  = critical temperature for senescence (C), constant b of aT^2+bT+c ,tabulated
    !Config if    = OK_STOMATE 
    !Config Def   = undef, undef, undef, undef, undef, 0., undef, 0., 0., .1, 0., 0., 0.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('SENESCENCE_TEMP_B',senescence_temp_b)

    !Config Key   = SENESCENCE_TEMP_A
    !Config Desc  = critical temperature for senescence (C), constant a of aT^2+bT+c , tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, 0., undef, 0., 0.,.00375, 0., 0., 0. 
    !Config Help  =
    !Config Units = [-] 
    CALL getin_p('SENESCENCE_TEMP_A',senescence_temp_a)

    !Config Key   = GDD_SENESCENCE
    !Config Desc  = minimum gdd to allow senescence of crops  
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 950., 4000.
    !Config Help  =
    !Config Units = [days] 
    CALL getin_p("GDD_SENESCENCE", gdd_senescence)

    !Config Key   = RESIDENCE_TIME
    !Config Desc  = residence time of trees (year
    !Config if    = OK_STOMATE
    !Config Def   = undef, 30, 30, 40, 40, 40, 80, 80, 80, 0, 0, 0, 0, 0
    !Config Help  =
    !Config Units = [years]
    CALL getin_p("RESIDENCE_TIME", residence_time)

    !Config Key   = TCM_CRIT
    !Config Desc  = critical tcm, tabulated 
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, 5.0, 15.5, 15.5, -8.0, -8.0, -8.0, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [C]
    CALL getin_p('TCM_CRIT',tcm_crit)

!qcj++ peatland
    !Config Key   =  WTPWET_CRIT
    !Config Desc  = critical minimum wtp position
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [m]
    CALL getin_p('WTPWET_CRIT',wtpwet_crit)

    !Config Key   =  WTPDRY_CRIT
    !Config Desc  = critical maximum wtp position
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [m]
    CALL getin_p('WTPDRY_CRIT',wtpdry_crit)

    !Config Key   =  WTP_CRIT
    !Config Desc  =  wtp position where plants start to be constrained by WT
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [m]
    CALL getin_p('WTP_CRIT',wtp_crit)

    !Config Key   =  WT_MORTALITY
    !Config Desc  = water table induced mortality
    !Config if    = OK_STOMATE
    !Config Def   = 0,0,0,0,0,0,0,0,0,0,0,0,0,
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('WT_MORTALITY',wt_mortality)


    IF (use_age_class) THEN
      ! Age classes
      ! I want to create a temporary array that indicates which "real" PFT starts
      ! on which index.  This could probably be put somewhere else, but this
      ! routine is only called once a year and this loop is not expensive.

      ! start_index and nagec_pft has length of nvm, only the beginning nvmap elements record
      ! the number of age groups for each PFT. The remaining elements will
      ! be -1 and when we use start_index and nagec_pft in the code, they're always looped
      ! over nvmap.
      start_index(:)=-1
      nagec_pft(:)=-1
      DO ivma=1,nvmap
        ! The start index is just the first place we find this real PFT.
        DO ivm=1,nvm
          IF(agec_group(ivm) .EQ. ivma)THEN
             start_index(ivma)=ivm
             j = 0
             DO WHILE(agec_group(ivm+j) .EQ. ivma)
               j = j+1
               ! we put this condition to handle the case that we have only one single
               ! age class for each PFT, but we still want to use the gross land cover
               ! change module.
               IF ((ivm+j .GT. nvmap .AND. nagec_tree .EQ. 1 .AND. nagec_herb .EQ. 1) .OR. ivm+j .GT. nvm) EXIT
             ENDDO
             nagec_pft(ivma)=j
             EXIT
          ENDIF
        ENDDO
      ENDDO
      ! Check to see if the calculation worked and we found indices for all of them.
      DO ivma=1,nvmap
        IF(start_index(ivma) .LT. 0)THEN
          WRITE(numout,*) 'Could not find a start index for one age class group!'
          WRITE(numout,*) 'Check the input file to make sure the following ivma appears in agec_group'
          WRITE(numout,*) 'ivma,nvmap',ivma,nvmap
          WRITE(numout,*) 'agec_group',agec_group(:)
          STOP
        ENDIF
      ENDDO
      ! Make sure input nagec_tree and nagec_herb is consistent
      DO ivma = 1,nvmap
        IF (nagec_pft(ivma) .GT.1) THEN
         IF (is_tree(start_index(ivma))) THEN
           IF (nagec_pft(ivma) .NE. nagec_tree) THEN
             WRITE(numout,*) 'The real number of age class for trees is not equal to nagec_tree'
             STOP
           ENDIF
         ELSE
           IF (nagec_pft(ivma) .NE. nagec_herb) THEN
             WRITE(numout,*) 'The real number of age class for grass/pasture is not equal to nagec_herb'
             STOP
           ENDIF
         ENDIF
        ENDIF
      ENDDO
    ENDIF

!
!WETLAND CH4 methane
!
!pss+
  !Config Key   = sdepth_v
  !Config Desc  = soil depth for wetland vegetation types
  !Config if    = CH4_CALCUL
  !Config Def   = /0,129,129,129,129,129,129,129,129,79,79,162,162/
  !Config Help  =
  !Config Units = [cm]
  CALL getin_p('SDEPTH_V',sdepth_v)

  !Config Key   = rdepth_v
  !Config Desc  = rooting depth for wetland vegetation types
  !Config if    = CH4_CALCUL
  !Config Def   = /0,64,64,64,64,64,64,64,64,39,39,81,81/
  !Config Help  =
  !Config Units = [cm]
  CALL getin_p('RDEPTH_V',rdepth_v)

  !Config Key   = tveg_v
  !Config Desc  = Plant mediated transport efficiency
  !Config if    = CH4_CALCUL
  !Config Def   = /0,1,1,1,1,1,1,1,1,10,10,15,15/
  !Config Help  =
  !Config Units = [-]
  CALL getin_p('TVEG_V',tveg_v)
!pss-
  
 END SUBROUTINE config_stomate_pft_parameters
!
!=
!


!! ================================================================================================================================
!! SUBROUTINE   : pft_parameters_clear
!!
!>\BRIEF         This subroutine deallocates memory at the end of the simulation. 
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

 SUBROUTINE pft_parameters_clear
   
   l_first_pft_parameters = .TRUE.
   
   IF (ALLOCATED(pft_to_mtc)) DEALLOCATE(pft_to_mtc)
   IF (ALLOCATED(PFT_name)) DEALLOCATE(PFT_name)
   IF (ALLOCATED(veget_ori_fixed_test_1)) DEALLOCATE(veget_ori_fixed_test_1)   
   IF (ALLOCATED(llaimax)) DEALLOCATE(llaimax)
   IF (ALLOCATED(llaimin)) DEALLOCATE(llaimin)
   IF (ALLOCATED(height_presc)) DEALLOCATE(height_presc)   
   IF (ALLOCATED(z0_over_height)) DEALLOCATE(z0_over_height)   
   IF (ALLOCATED(ratio_z0m_z0h)) DEALLOCATE(ratio_z0m_z0h)   
   IF (ALLOCATED(type_of_lai)) DEALLOCATE(type_of_lai)
   IF (ALLOCATED(is_tree)) DEALLOCATE(is_tree)
   IF (ALLOCATED(natural)) DEALLOCATE(natural)
! dgvmjc
    IF (ALLOCATED(pasture)) DEALLOCATE(pasture)
! end dgvmjc
   IF (ALLOCATED(is_deciduous)) DEALLOCATE(is_deciduous)
   IF (ALLOCATED(is_evergreen)) DEALLOCATE(is_evergreen)
   IF (ALLOCATED(is_needleleaf)) DEALLOCATE(is_needleleaf)
   IF (ALLOCATED(is_tropical)) DEALLOCATE(is_tropical)
   IF (ALLOCATED(humcste)) DEALLOCATE(humcste)
   IF (ALLOCATED(pref_soil_veg)) DEALLOCATE(pref_soil_veg)
   IF (ALLOCATED(agec_group)) DEALLOCATE(agec_group)
   IF (ALLOCATED(age_class_bound)) DEALLOCATE(age_class_bound)
   IF (ALLOCATED(start_index)) DEALLOCATE(start_index)
   IF (ALLOCATED(nagec_pft)) DEALLOCATE(nagec_pft)
   IF (ALLOCATED(is_c4)) DEALLOCATE(is_c4)  
!!!qcj++ peatland
   IF (ALLOCATED(is_peat)) DEALLOCATE(is_peat)
   IF (ALLOCATED(is_croppeat)) DEALLOCATE(is_croppeat)
   IF (ALLOCATED(is_shrubpeat)) DEALLOCATE(is_shrubpeat)
   IF (ALLOCATED(is_mosspeat)) DEALLOCATE(is_mosspeat)
   IF (ALLOCATED(is_mineralwet)) DEALLOCATE(is_mineralwet)
   IF (ALLOCATED(vcmax_fix)) DEALLOCATE(vcmax_fix)
   IF (ALLOCATED(downregulation_co2_coeff)) DEALLOCATE(downregulation_co2_coeff) 
   IF (ALLOCATED(downregulation_co2_coeff_new)) DEALLOCATE(downregulation_co2_coeff_new)
   IF (ALLOCATED(E_KmC)) DEALLOCATE(E_KmC)
   IF (ALLOCATED(E_KmO)) DEALLOCATE(E_KmO)
   IF (ALLOCATED(E_Sco)) DEALLOCATE(E_Sco)
   IF (ALLOCATED(E_gamma_star)) DEALLOCATE(E_gamma_star)
   IF (ALLOCATED(E_Vcmax)) DEALLOCATE(E_Vcmax)
   IF (ALLOCATED(E_Jmax)) DEALLOCATE(E_Jmax)
   IF (ALLOCATED(aSV)) DEALLOCATE(aSV)
   IF (ALLOCATED(bSV)) DEALLOCATE(bSV)
   IF (ALLOCATED(tphoto_min)) DEALLOCATE(tphoto_min)
   IF (ALLOCATED(tphoto_max)) DEALLOCATE(tphoto_max)
   IF (ALLOCATED(aSJ)) DEALLOCATE(aSJ)
   IF (ALLOCATED(bSJ)) DEALLOCATE(bSJ)
   IF (ALLOCATED(D_Vcmax)) DEALLOCATE(D_Vcmax)
   IF (ALLOCATED(D_Jmax)) DEALLOCATE(D_Jmax)
   IF (ALLOCATED(E_gm)) DEALLOCATE(E_gm)  
   IF (ALLOCATED(S_gm)) DEALLOCATE(S_gm)  
   IF (ALLOCATED(D_gm)) DEALLOCATE(D_gm)
   IF (ALLOCATED(E_Rd)) DEALLOCATE(E_Rd)
   IF (ALLOCATED(Vcmax25)) DEALLOCATE(Vcmax25)
   IF (ALLOCATED(arJV)) DEALLOCATE(arJV)
   IF (ALLOCATED(brJV)) DEALLOCATE(brJV)
   IF (ALLOCATED(KmC25)) DEALLOCATE(KmC25)
   IF (ALLOCATED(KmO25)) DEALLOCATE(KmO25)
   IF (ALLOCATED(Sco25)) DEALLOCATE(Sco25) 
   IF (ALLOCATED(gm25)) DEALLOCATE(gm25)
   IF (ALLOCATED(gamma_star25)) DEALLOCATE(gamma_star25)
   IF (ALLOCATED(a1)) DEALLOCATE(a1)
   IF (ALLOCATED(b1)) DEALLOCATE(b1)
   IF (ALLOCATED(g0)) DEALLOCATE(g0)
   IF (ALLOCATED(h_protons)) DEALLOCATE(h_protons)
   IF (ALLOCATED(fpsir)) DEALLOCATE(fpsir)
   IF (ALLOCATED(fQ)) DEALLOCATE(fQ)
   IF (ALLOCATED(fpseudo)) DEALLOCATE(fpseudo)
   IF (ALLOCATED(kp)) DEALLOCATE(kp)
   IF (ALLOCATED(alpha)) DEALLOCATE(alpha)
   IF (ALLOCATED(gbs)) DEALLOCATE(gbs)
   IF (ALLOCATED(theta)) DEALLOCATE(theta)
   IF (ALLOCATED(alpha_LL)) DEALLOCATE(alpha_LL)
   IF (ALLOCATED(stress_vcmax)) DEALLOCATE(stress_vcmax) 
   IF (ALLOCATED(stress_gs)) DEALLOCATE(stress_gs) 
   IF (ALLOCATED(stress_gm)) DEALLOCATE(stress_gm)
   IF (ALLOCATED(ext_coeff)) DEALLOCATE(ext_coeff)
   IF (ALLOCATED(ext_coeff_vegetfrac)) DEALLOCATE(ext_coeff_vegetfrac)
   IF (ALLOCATED(rveg_pft)) DEALLOCATE(rveg_pft)
   IF (ALLOCATED(rstruct_const)) DEALLOCATE(rstruct_const)
   IF (ALLOCATED(kzero)) DEALLOCATE(kzero)
   IF (ALLOCATED(wmax_veg)) DEALLOCATE(wmax_veg)
   IF (ALLOCATED(throughfall_by_pft)) DEALLOCATE(throughfall_by_pft)
   IF (ALLOCATED(snowa_aged_vis)) DEALLOCATE(snowa_aged_vis)
   IF (ALLOCATED(snowa_aged_nir)) DEALLOCATE(snowa_aged_nir)
   IF (ALLOCATED(snowa_dec_vis)) DEALLOCATE(snowa_dec_vis)
   IF (ALLOCATED(snowa_dec_nir)) DEALLOCATE(snowa_dec_nir)
   IF (ALLOCATED(alb_leaf_vis)) DEALLOCATE(alb_leaf_vis)
   IF (ALLOCATED(alb_leaf_nir)) DEALLOCATE(alb_leaf_nir)   
   !chaoyue+
   IF (ALLOCATED(permafrost_veg_exists)) DEALLOCATE(permafrost_veg_exists)   
   !chaoyue-
   IF (ALLOCATED(em_factor_isoprene)) DEALLOCATE(em_factor_isoprene)
   IF (ALLOCATED(em_factor_monoterpene)) DEALLOCATE(em_factor_monoterpene)
   IF (ALLOCATED(em_factor_apinene)) DEALLOCATE(em_factor_apinene)
   IF (ALLOCATED(em_factor_bpinene)) DEALLOCATE(em_factor_bpinene)
   IF (ALLOCATED(em_factor_limonene)) DEALLOCATE(em_factor_limonene)
   IF (ALLOCATED(em_factor_myrcene)) DEALLOCATE(em_factor_myrcene)
   IF (ALLOCATED(em_factor_sabinene)) DEALLOCATE(em_factor_sabinene)
   IF (ALLOCATED(em_factor_camphene)) DEALLOCATE(em_factor_camphene)
   IF (ALLOCATED(em_factor_3carene)) DEALLOCATE(em_factor_3carene)
   IF (ALLOCATED(em_factor_tbocimene)) DEALLOCATE(em_factor_tbocimene)
   IF (ALLOCATED(em_factor_othermonot)) DEALLOCATE(em_factor_othermonot)
   IF (ALLOCATED(em_factor_sesquiterp)) DEALLOCATE(em_factor_sesquiterp)
   IF (ALLOCATED(em_factor_ORVOC)) DEALLOCATE(em_factor_ORVOC)
   IF (ALLOCATED(em_factor_OVOC)) DEALLOCATE(em_factor_OVOC)
   IF (ALLOCATED(em_factor_MBO)) DEALLOCATE(em_factor_MBO)
   IF (ALLOCATED(em_factor_methanol)) DEALLOCATE(em_factor_methanol)
   IF (ALLOCATED(em_factor_acetone)) DEALLOCATE(em_factor_acetone)
   IF (ALLOCATED(em_factor_acetal)) DEALLOCATE(em_factor_acetal)
   IF (ALLOCATED(em_factor_formal)) DEALLOCATE(em_factor_formal)
   IF (ALLOCATED(em_factor_acetic)) DEALLOCATE(em_factor_acetic)
   IF (ALLOCATED(em_factor_formic)) DEALLOCATE(em_factor_formic)
   IF (ALLOCATED(em_factor_no_wet)) DEALLOCATE(em_factor_no_wet)
   IF (ALLOCATED(em_factor_no_dry)) DEALLOCATE(em_factor_no_dry)
   IF (ALLOCATED(Larch)) DEALLOCATE(Larch)
   IF (ALLOCATED(leaf_tab)) DEALLOCATE(leaf_tab)
   IF (ALLOCATED(sla)) DEALLOCATE(sla)
   IF (ALLOCATED(availability_fact)) DEALLOCATE(availability_fact)
   IF (ALLOCATED(R0)) DEALLOCATE(R0)
   IF (ALLOCATED(S0)) DEALLOCATE(S0)
   IF (ALLOCATED(L0)) DEALLOCATE(L0)
   IF (ALLOCATED(frac_growthresp)) DEALLOCATE(frac_growthresp)
   IF (ALLOCATED(maint_resp_slope)) DEALLOCATE(maint_resp_slope)
   IF (ALLOCATED(maint_resp_slope_c)) DEALLOCATE(maint_resp_slope_c)
   IF (ALLOCATED(maint_resp_slope_b)) DEALLOCATE(maint_resp_slope_b)
   IF (ALLOCATED(maint_resp_slope_a)) DEALLOCATE(maint_resp_slope_a)
   IF (ALLOCATED(coeff_maint_zero)) DEALLOCATE(coeff_maint_zero)
   IF (ALLOCATED(cm_zero_leaf)) DEALLOCATE(cm_zero_leaf)
   IF (ALLOCATED(cm_zero_sapabove)) DEALLOCATE(cm_zero_sapabove)
   IF (ALLOCATED(cm_zero_sapbelow)) DEALLOCATE(cm_zero_sapbelow)
   IF (ALLOCATED(cm_zero_heartabove)) DEALLOCATE(cm_zero_heartabove)
   IF (ALLOCATED(cm_zero_heartbelow)) DEALLOCATE(cm_zero_heartbelow)
   IF (ALLOCATED(cm_zero_root)) DEALLOCATE(cm_zero_root)
   IF (ALLOCATED(cm_zero_fruit)) DEALLOCATE(cm_zero_fruit)
   IF (ALLOCATED(cm_zero_carbres)) DEALLOCATE(cm_zero_carbres)
   IF (ALLOCATED(flam)) DEALLOCATE(flam)
   IF (ALLOCATED(resist)) DEALLOCATE(resist)
   !spitfire
   IF (ALLOCATED(dens_fuel)) DEALLOCATE(dens_fuel)
   IF (ALLOCATED(f_sh)) DEALLOCATE(f_sh)
   IF (ALLOCATED(crown_length)) DEALLOCATE(crown_length)
   IF (ALLOCATED(BTpar1)) DEALLOCATE(BTpar1)
   IF (ALLOCATED(BTpar2)) DEALLOCATE(BTpar2)
   IF (ALLOCATED(r_ck)) DEALLOCATE(r_ck)
   IF (ALLOCATED(p_ck)) DEALLOCATE(p_ck)
   IF (ALLOCATED(ef_CO2)) DEALLOCATE(ef_CO2)
   IF (ALLOCATED(ef_CO)) DEALLOCATE(ef_CO)
   IF (ALLOCATED(ef_CH4)) DEALLOCATE(ef_CH4)
   IF (ALLOCATED(ef_VOC)) DEALLOCATE(ef_VOC)
   IF (ALLOCATED(ef_TPM)) DEALLOCATE(ef_TPM)
   IF (ALLOCATED(ef_NOx)) DEALLOCATE(ef_NOx)
   IF (ALLOCATED(me)) DEALLOCATE(me)
   IF (ALLOCATED(fire_max_cf_100hr)) DEALLOCATE(fire_max_cf_100hr)
   IF (ALLOCATED(fire_max_cf_1000hr)) DEALLOCATE(fire_max_cf_1000hr)
   !endspit
   IF (ALLOCATED(coeff_lcchange_1)) DEALLOCATE(coeff_lcchange_1)
   IF (ALLOCATED(coeff_lcchange_10)) DEALLOCATE(coeff_lcchange_10)
   IF (ALLOCATED(coeff_lcchange_100)) DEALLOCATE(coeff_lcchange_100)
   IF (ALLOCATED(coeff_indwood_1)) DEALLOCATE(coeff_indwood_1)
   IF (ALLOCATED(coeff_indwood_10)) DEALLOCATE(coeff_indwood_10)
   IF (ALLOCATED(coeff_indwood_100)) DEALLOCATE(coeff_indwood_100)
   IF (ALLOCATED(lai_max_to_happy)) DEALLOCATE(lai_max_to_happy)
   IF (ALLOCATED(lai_max)) DEALLOCATE(lai_max)
!gmjc
   IF (ALLOCATED(is_grassland_manag))DEALLOCATE(is_grassland_manag)
   IF (ALLOCATED(is_grassland_cut))DEALLOCATE(is_grassland_cut)
   IF (ALLOCATED(is_grassland_grazed))DEALLOCATE(is_grassland_grazed)
   IF (ALLOCATED(nb_year_management)) DEALLOCATE(nb_year_management)
   IF (ALLOCATED(management_intensity)) DEALLOCATE(management_intensity)
   IF (ALLOCATED(management_start)) DEALLOCATE(management_start)
   IF (ALLOCATED(deposition_start)) DEALLOCATE(deposition_start)
   IF(ALLOCATED(sla_max))DEALLOCATE(sla_max)
   IF(ALLOCATED(sla_min))DEALLOCATE(sla_min)
!end gmjc
   !bioenergy
   IF (ALLOCATED(is_bioe1))DEALLOCATE(is_bioe1)
   IF (ALLOCATED(pheno_model)) DEALLOCATE(pheno_model)
   IF (ALLOCATED(pheno_type)) DEALLOCATE(pheno_type)
   IF (ALLOCATED(pheno_gdd_crit_c)) DEALLOCATE(pheno_gdd_crit_c)
   IF (ALLOCATED(pheno_gdd_crit_b)) DEALLOCATE(pheno_gdd_crit_b)
   IF (ALLOCATED(pheno_gdd_crit_a)) DEALLOCATE(pheno_gdd_crit_a)
   IF (ALLOCATED(pheno_gdd_crit)) DEALLOCATE(pheno_gdd_crit)
   IF (ALLOCATED(pheno_moigdd_t_crit)) DEALLOCATE(pheno_moigdd_t_crit)
   IF (ALLOCATED(ngd_crit)) DEALLOCATE(ngd_crit)
   IF (ALLOCATED(ncdgdd_temp)) DEALLOCATE(ncdgdd_temp)
   IF (ALLOCATED(hum_frac)) DEALLOCATE(hum_frac)
   IF (ALLOCATED(hum_min_time)) DEALLOCATE(hum_min_time)
   IF (ALLOCATED(tau_sap)) DEALLOCATE(tau_sap)
   IF (ALLOCATED(tau_leafinit)) DEALLOCATE(tau_leafinit)
   IF (ALLOCATED(tau_fruit)) DEALLOCATE(tau_fruit)
   IF (ALLOCATED(ecureuil)) DEALLOCATE(ecureuil)
   IF (ALLOCATED(alloc_min)) DEALLOCATE(alloc_min)
   IF (ALLOCATED(alloc_max)) DEALLOCATE(alloc_max)
   IF (ALLOCATED(demi_alloc)) DEALLOCATE(demi_alloc)
   IF (ALLOCATED(leaflife_tab)) DEALLOCATE(leaflife_tab)
   IF (ALLOCATED(leaffall)) DEALLOCATE(leaffall)
   IF (ALLOCATED(leafagecrit)) DEALLOCATE(leafagecrit)
   IF (ALLOCATED(senescence_type)) DEALLOCATE(senescence_type)
   IF (ALLOCATED(senescence_hum)) DEALLOCATE(senescence_hum)
   IF (ALLOCATED(nosenescence_hum)) DEALLOCATE(nosenescence_hum)
   IF (ALLOCATED(max_turnover_time)) DEALLOCATE(max_turnover_time)
   IF (ALLOCATED(min_turnover_time)) DEALLOCATE(min_turnover_time)
   IF (ALLOCATED(min_leaf_age_for_senescence)) DEALLOCATE(min_leaf_age_for_senescence)
   IF (ALLOCATED(senescence_temp_c)) DEALLOCATE(senescence_temp_c)
   IF (ALLOCATED(senescence_temp_b)) DEALLOCATE(senescence_temp_b)
   IF (ALLOCATED(senescence_temp_a)) DEALLOCATE(senescence_temp_a)
   IF (ALLOCATED(senescence_temp)) DEALLOCATE(senescence_temp)
   IF (ALLOCATED(gdd_senescence)) DEALLOCATE(gdd_senescence)
   IF (ALLOCATED(residence_time)) DEALLOCATE(residence_time)
   IF (ALLOCATED(tmin_crit)) DEALLOCATE(tmin_crit)
   IF (ALLOCATED(tcm_crit)) DEALLOCATE(tcm_crit)
!qcj++ peatland
   IF (ALLOCATED(wtpwet_crit)) DEALLOCATE(wtpwet_crit)
   IF (ALLOCATED(wtpdry_crit)) DEALLOCATE(wtpdry_crit)
   IF (ALLOCATED(wtp_crit)) DEALLOCATE(wtp_crit)
   IF (ALLOCATED(wt_mortality)) DEALLOCATE(wt_mortality)   

   IF (ALLOCATED(lai_initmin)) DEALLOCATE(lai_initmin)
   IF (ALLOCATED(bm_sapl)) DEALLOCATE(bm_sapl)
   IF (ALLOCATED(migrate)) DEALLOCATE(migrate)
   IF (ALLOCATED(maxdia)) DEALLOCATE(maxdia)
   IF (ALLOCATED(cn_sapl)) DEALLOCATE(cn_sapl)
   IF (ALLOCATED(leaf_timecst)) DEALLOCATE(leaf_timecst)
   
!pss+
   IF (ALLOCATED(rdepth_v)) DEALLOCATE(rdepth_v)
   IF (ALLOCATED(sdepth_v)) DEALLOCATE(sdepth_v)
   IF (ALLOCATED(tveg_v)) DEALLOCATE(tveg_v)
!pss-

!!!!! crop parameters
   IF (ALLOCATED(irrig_threshold)) DEALLOCATE(irrig_threshold)
   IF (ALLOCATED(irrig_fulfill)) DEALLOCATE(irrig_fulfill)

   ! DEALLOCATE FOR crop

   IF(ALLOCATED(ok_LAIdev))DEALLOCATE(ok_LAIdev)

   IF(ALLOCATED(SP_codeplante))DEALLOCATE(SP_codeplante)
   IF(ALLOCATED(SP_stade0))DEALLOCATE(SP_stade0)
   IF(ALLOCATED(SP_iplt0))DEALLOCATE(SP_iplt0)
   IF(ALLOCATED(SP_nbox))DEALLOCATE(SP_nbox)
   IF(ALLOCATED(SP_iwater))DEALLOCATE(SP_iwater)
   IF(ALLOCATED(SP_codesimul))DEALLOCATE(SP_codesimul)
   IF(ALLOCATED(SP_codelaitr))DEALLOCATE(SP_codelaitr)
   IF(ALLOCATED(SP_slamax))DEALLOCATE(SP_slamax)
   IF(ALLOCATED(SP_slamin))DEALLOCATE(SP_slamin)
   IF(ALLOCATED(SP_codeperenne))DEALLOCATE(SP_codeperenne)
   IF(ALLOCATED(SP_codcueille))DEALLOCATE(SP_codcueille)
   IF(ALLOCATED(SP_codegdh))DEALLOCATE(SP_codegdh)
   IF(ALLOCATED(SP_codetemp))DEALLOCATE(SP_codetemp)
   IF(ALLOCATED(SP_coderetflo))DEALLOCATE(SP_coderetflo)
   IF(ALLOCATED(SP_codeinnact))DEALLOCATE(SP_codeinnact)
   IF(ALLOCATED(SP_codeh2oact))DEALLOCATE(SP_codeh2oact)
   IF(ALLOCATED(SP_stressdev))DEALLOCATE(SP_stressdev)
   IF(ALLOCATED(SP_innlai))DEALLOCATE(SP_innlai)
   IF(ALLOCATED(SP_innsenes))DEALLOCATE(SP_innsenes)
   IF(ALLOCATED(SP_codebfroid))DEALLOCATE(SP_codebfroid)
   IF(ALLOCATED(SP_codephot))DEALLOCATE(SP_codephot)
   IF(ALLOCATED(SP_codedormance))DEALLOCATE(SP_codedormance)
   IF(ALLOCATED(SP_codefauche))DEALLOCATE(SP_codefauche)
   IF(ALLOCATED(SP_codetempfauche))DEALLOCATE(SP_codetempfauche)
   IF(ALLOCATED(SP_codlainet))DEALLOCATE(SP_codlainet)
   IF(ALLOCATED(SP_codeindetermin))DEALLOCATE(SP_codeindetermin)
   IF(ALLOCATED(SP_codeinitprec))DEALLOCATE(SP_codeinitprec)
   IF(ALLOCATED(SP_culturean))DEALLOCATE(SP_culturean)
   IF(ALLOCATED(SP_jvc))DEALLOCATE(SP_jvc)
   IF(ALLOCATED(SP_tfroid))DEALLOCATE(SP_tfroid)
   IF(ALLOCATED(SP_ampfroid))DEALLOCATE(SP_ampfroid)
   IF(ALLOCATED(SP_jvcmini))DEALLOCATE(SP_jvcmini)
   IF(ALLOCATED(SP_tgmin))DEALLOCATE(SP_tgmin)
   IF(ALLOCATED(SP_stpltger))DEALLOCATE(SP_stpltger)
   IF(ALLOCATED(SP_profsem))DEALLOCATE(SP_profsem)
   IF(ALLOCATED(SP_propjgermin))DEALLOCATE(SP_propjgermin)
   IF(ALLOCATED(SP_tdmax))DEALLOCATE(SP_tdmax)
   IF(ALLOCATED(SP_nbjgerlim))DEALLOCATE(SP_nbjgerlim)
   IF(ALLOCATED(SP_densitesem))DEALLOCATE(SP_densitesem)
   IF(ALLOCATED(SP_vigueurbat))DEALLOCATE(SP_vigueurbat)
   IF(ALLOCATED(SP_codepluiepoquet))DEALLOCATE(SP_codepluiepoquet)
   IF(ALLOCATED(SP_codehypo))DEALLOCATE(SP_codehypo)
   IF(ALLOCATED(SP_elmax))DEALLOCATE(SP_elmax)
   IF(ALLOCATED(SP_belong))DEALLOCATE(SP_belong)
   IF(ALLOCATED(SP_celong))DEALLOCATE(SP_celong)
   IF(ALLOCATED(SP_nlevlim1))DEALLOCATE(SP_nlevlim1)
   IF(ALLOCATED(SP_nlevlim2))DEALLOCATE(SP_nlevlim2)
   IF(ALLOCATED(SP_codrecolte))DEALLOCATE(SP_codrecolte)
   IF(ALLOCATED(SP_variete))DEALLOCATE(SP_variete)
   IF(ALLOCATED(SP_codegermin))DEALLOCATE(SP_codegermin)

   IF(ALLOCATED(S_codeulaivernal))DEALLOCATE(S_codeulaivernal)
   IF(ALLOCATED(SP_swfacmin))DEALLOCATE(SP_swfacmin)
   IF(ALLOCATED(SP_neffmax))DEALLOCATE(SP_neffmax)
   IF(ALLOCATED(SP_nsatrat))DEALLOCATE(SP_nsatrat)

   ! STICS:: LAI CALCULATION
   IF(ALLOCATED(SP_laiplantule))DEALLOCATE(SP_laiplantule)
   IF(ALLOCATED(SP_vlaimax))DEALLOCATE(SP_vlaimax)
   IF(ALLOCATED(SP_stlevamf))DEALLOCATE(SP_stlevamf)
   IF(ALLOCATED(SP_stdrpmat))DEALLOCATE(SP_stdrpmat)
   IF(ALLOCATED(SP_stamflax))DEALLOCATE(SP_stamflax)

   IF(ALLOCATED(SP_udlaimax))DEALLOCATE(SP_udlaimax)
   IF(ALLOCATED(SP_laicomp))DEALLOCATE(SP_laicomp)
   IF(ALLOCATED(SP_adens))DEALLOCATE(SP_adens)
   IF(ALLOCATED(SP_bdens))DEALLOCATE(SP_bdens)
   IF(ALLOCATED(SP_tcxstop))DEALLOCATE(SP_tcxstop)
   IF(ALLOCATED(SP_tcmax))DEALLOCATE(SP_tcmax)
   IF(ALLOCATED(SP_tcmin))DEALLOCATE(SP_tcmin)
   IF(ALLOCATED(SP_dlaimax))DEALLOCATE(SP_dlaimax)
   IF(ALLOCATED(SP_dlaimin))DEALLOCATE(SP_dlaimin)
   IF(ALLOCATED(SP_pentlaimax))DEALLOCATE(SP_pentlaimax)
   IF(ALLOCATED(SP_tigefeuil))DEALLOCATE(SP_tigefeuil)

   IF(ALLOCATED(SP_stlaxsen))DEALLOCATE(SP_stlaxsen)
   IF(ALLOCATED(SP_stsenlan))DEALLOCATE(SP_stsenlan)
   IF(ALLOCATED(SP_stlevdrp))DEALLOCATE(SP_stlevdrp)
   IF(ALLOCATED(SP_stflodrp))DEALLOCATE(SP_stflodrp)
   IF(ALLOCATED(SP_stdrpdes))DEALLOCATE(SP_stdrpdes)
   IF(ALLOCATED(SP_phyllotherme))DEALLOCATE(SP_phyllotherme)

   IF(ALLOCATED(SP_lai0))DEALLOCATE(SP_lai0)
   IF(ALLOCATED(SP_tustressmin))DEALLOCATE(SP_tustressmin)


   ! STICS:: LAI SENESCENCE
   IF(ALLOCATED(SP_nbfgellev))DEALLOCATE(SP_nbfgellev)
   IF(ALLOCATED(SP_ratiodurvieI))DEALLOCATE(SP_ratiodurvieI)
   IF(ALLOCATED(SP_durvieF))DEALLOCATE(SP_durvieF)
   IF(ALLOCATED(SP_ratiosen))DEALLOCATE(SP_ratiosen)
   IF(ALLOCATED(SP_tdmin))DEALLOCATE(SP_tdmin)
   
   ! STICS:: F_humerac
 
   IF(ALLOCATED(SP_sensrsec))DEALLOCATE(SP_sensrsec)
   ! STICS:: gel

   IF(ALLOCATED(SP_codgellev))DEALLOCATE(SP_codgellev)
   IF(ALLOCATED(SP_codgeljuv))DEALLOCATE(SP_codgeljuv)
   IF(ALLOCATED(SP_codgelveg))DEALLOCATE(SP_codgelveg)
   IF(ALLOCATED(SP_tletale))DEALLOCATE(SP_tletale)
   IF(ALLOCATED(SP_tdebgel))DEALLOCATE(SP_tdebgel)
   IF(ALLOCATED(SP_tgellev10))DEALLOCATE(SP_tgellev10)
   IF(ALLOCATED(SP_tgellev90))DEALLOCATE(SP_tgellev90)

   IF(ALLOCATED(SP_tgeljuv10))DEALLOCATE(SP_tgeljuv10)
   IF(ALLOCATED(SP_tgeljuv90))DEALLOCATE(SP_tgeljuv90)
   IF(ALLOCATED(SP_tgelveg10))DEALLOCATE(SP_tgelveg10)
   IF(ALLOCATED(SP_tgelveg90))DEALLOCATE(SP_tgelveg90)




   ! STICS:: Photoperiod
  
   IF(ALLOCATED(SP_sensiphot))DEALLOCATE(SP_sensiphot)
   IF(ALLOCATED(SP_phosat))DEALLOCATE(SP_phosat)
   IF(ALLOCATED(SP_phobase))DEALLOCATE(SP_phobase)
   
   ! STICS:: CARBON ALLOCATION
     
   IF(ALLOCATED(SP_stoprac))DEALLOCATE(SP_stoprac)
   IF(ALLOCATED(SP_zracplantule))DEALLOCATE(SP_zracplantule)
   IF(ALLOCATED(SP_codtrophrac))DEALLOCATE(SP_codtrophrac)
   IF(ALLOCATED(SP_repracpermax))DEALLOCATE(SP_repracpermax)
   IF(ALLOCATED(SP_repracpermin))DEALLOCATE(SP_repracpermin)
   IF(ALLOCATED(SP_krepracperm))DEALLOCATE(SP_krepracperm)
   IF(ALLOCATED(SP_repracseumax))DEALLOCATE(SP_repracseumax)
   IF(ALLOCATED(SP_repracseumin))DEALLOCATE(SP_repracseumin)
   IF(ALLOCATED(SP_krepracseu))DEALLOCATE(SP_krepracseu)
   IF(ALLOCATED(SP_codetemprac))DEALLOCATE(SP_codetemprac)
   IF(ALLOCATED(SP_codedyntalle))DEALLOCATE(SP_codedyntalle)
   IF(ALLOCATED(SP_nbjgrain))DEALLOCATE(SP_nbjgrain)
   IF(ALLOCATED(SP_maxgs))DEALLOCATE(SP_maxgs)
   IF(ALLOCATED(SP_codgelflo))DEALLOCATE(SP_codgelflo)
   IF(ALLOCATED(SP_tgelflo10))DEALLOCATE(SP_tgelflo10)
   IF(ALLOCATED(SP_tgelflo90))DEALLOCATE(SP_tgelflo90)
   IF(ALLOCATED(SP_cgrain))DEALLOCATE(SP_cgrain)
   IF(ALLOCATED(SP_cgrainv0))DEALLOCATE(SP_cgrainv0)
   IF(ALLOCATED(SP_nbgrmax))DEALLOCATE(SP_nbgrmax)
   IF(ALLOCATED(SP_nbgrmin))DEALLOCATE(SP_nbgrmin)
   IF(ALLOCATED(SP_codazofruit))DEALLOCATE(SP_codazofruit)
   IF(ALLOCATED(SP_codeir))DEALLOCATE(SP_codeir)
   IF(ALLOCATED(SP_vitircarb))DEALLOCATE(SP_vitircarb)
   IF(ALLOCATED(SP_irmax))DEALLOCATE(SP_irmax)
   IF(ALLOCATED(SP_vitircarbT))DEALLOCATE(SP_vitircarbT)
   IF(ALLOCATED(SP_codetremp))DEALLOCATE(SP_codetremp)
   IF(ALLOCATED(SP_tminremp))DEALLOCATE(SP_tminremp)
   IF(ALLOCATED(SP_tmaxremp))DEALLOCATE(SP_tmaxremp)
   IF(ALLOCATED(SP_pgrainmaxi))DEALLOCATE(SP_pgrainmaxi)
 
   IF(ALLOCATED(SP_DY_INN))DEALLOCATE(SP_DY_INN)
   IF(ALLOCATED(SP_avenfert))DEALLOCATE(SP_avenfert)
!!!!! end crop parameters

 END SUBROUTINE pft_parameters_clear

END MODULE pft_parameters
