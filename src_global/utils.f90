!
! Place here all those small routines that do not fit in orchidee logic yet
!
!
!< $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_global/interpol_help.f90 $
!< $Date: 2016-06-17 13:26:43 +0200 (Fri, 17 Jun 2016) $
!< $Author: albert.jornet $
!< $Revision: 3564 $
!
!
MODULE utils 

  ! Modules used :

  USE netcdf
  USE defprec
  USE ioipsl_para

  IMPLICIT NONE

  PRIVATE
  PUBLIC nccheck 
  !
CONTAINS
  !
!! ================================================================================================================================
!! SUBROUTINE 	: nccheck 
!!
!>\BRIEF        Check for netcdf exit status 
!!
!! DESCRIPTION  : Launch an orchidee error message if status variable contains a netcdf error
!!
!! RECENT CHANGE(S) : None
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE nccheck(status)
    INTEGER(i_std), INTENT (IN)         :: status
    CHARACTER(LEN=200)                  :: mesg
    
    IF(status /= nf90_noerr) THEN
      
      WRITE(numout, *) trim(nf90_strerror(status))
      CALL ipslerr_p(3, 'nccheck', 'Netcdf error', 'Check out_orchide_XXXX output files', 'for more information')
    END IF  
  END SUBROUTINE nccheck
!
END MODULE utils 
