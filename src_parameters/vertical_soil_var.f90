! =================================================================================================================================
! MODULE       : vertical_soil_var
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        
!!
!!\n DESCRIPTION: 
!!                
!! RECENT CHANGE(S):
!!
!! REFERENCE(S)	: 
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_parameters/vertical_soil_var.f90 $
!! $Date: 2017-10-26 15:35:04 +0200 (Thu, 26 Oct 2017) $
!! $Revision: 4719 $
!! \n
!_ ================================================================================================================================

MODULE vertical_soil_var

  USE defprec

  IMPLICIT NONE
  PUBLIC

  !! Dimensioning parameters
  INTEGER(i_std), SAVE      :: ngrnd     !! Number of soil layer for thermo (unitless)
!$OMP THREADPRIVATE(ngrnd)
  INTEGER(i_std), SAVE      :: nslm      !! Number of levels in CWRR (unitless)
!$OMP THREADPRIVATE(nslm)
  REAL(r_std), SAVE         :: zmaxh     !! Maximum depth of soil reservoir in hydrol (m). Old name dpu_max or depth_Wmax
!$OMP THREADPRIVATE(zmaxh)
  REAL(r_std), SAVE         :: zmaxt     !! Maximum depth of the soil thermodynamics (m)
!$OMP THREADPRIVATE(zmaxt)

  !! Variables defining the vertical layering in soil moisture and temperature
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:) :: znt          !! Depth of nodes for thermal (m) 
!$OMP THREADPRIVATE(znt)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:) :: znh          !! Depth of nodes for hydrology (m)
!$OMP THREADPRIVATE(znh)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:) :: dnh          !! Distance between the current node and the one above for hydrology (m)
!$OMP THREADPRIVATE(dnh)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:) :: dlh          !! Soil layer thickness for hydrology (m) 
!$OMP THREADPRIVATE(dlh)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:) :: dlt          !! Soil layer thickness for thermal (m)
!$OMP THREADPRIVATE(dlt)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:) :: zlh          !! Depth of lower layer-interface for hydrology (m)
!$OMP THREADPRIVATE(zlh)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:) :: zlt          !! Depth of lower layer-interface for thermal (m)
!$OMP THREADPRIVATE(zlt)

  REAL(r_std),ALLOCATABLE, DIMENSION(:),SAVE :: diaglev        !! The lower limit of the layer on which soil moisture
                                                               !! (relative) and temperature are going to be diagnosed.
                                                               !! These variables are made for transfering the information
                                                               !! to the biogeophyical processes modelled in STOMATE. 
!$OMP THREADPRIVATE(diaglev)

END MODULE vertical_soil_var
