! ===============================================================================================================================
! MODULE       : mod_orchidee_omp_data
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        Contains initialization and allocation of variables and functions related to OpenMP parallelization.
!!
!! \n DESCRIPTION : Contains subroutines for initialization and allocation of variables and functions related to 
!!                  OpenMP parallelization.
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCES(S)    : None
!!
!! SVN              :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_parallel/mod_orchidee_omp_data.F90 $
!! $Date: 2017-10-26 15:35:04 +0200 (Thu, 26 Oct 2017) $
!! $Revision: 4719 $
!! \n
!_ ================================================================================================================================
MODULE mod_orchidee_omp_data

!-
  USE defprec
  USE ioipsl
  USE mod_orchidee_para_var

  IMPLICIT NONE

CONTAINS



!!  =============================================================================================================================
!! SUBROUTINE:    barrier2_omp
!!
!>\BRIEF	this routine call two omp barrier to prevent a specific bug when orchidee is coupled to lmdz
!!
!! DESCRIPTION:	  this routine call two omp barrier to prevent a specific bug when orchidee is coupled to lmdz
!!
!!                
!! \n
!_ ==============================================================================================================================
    SUBROUTINE barrier2_omp()

    IMPLICIT NONE

!$OMP BARRIER

  END SUBROUTINE barrier2_omp



  
!!  =============================================================================================================================
!! SUBROUTINE:   Init_orchidee_omp 
!!
!>\BRIEF	define the variables is_ok_omp, is_omp_root, omp_size and omp_rank  in the offline case
!!
!! DESCRIPTION:	  define the variables is_ok_omp, is_omp_root, omp_size and omp_rank  in the offline case
!!
!!                
!! \n
!_ ==============================================================================================================================
  SUBROUTINE Init_orchidee_omp
  IMPLICIT NONE
  
#ifdef CPP_OMP
    IF (is_omp_root) THEN
        is_ok_omp=.TRUE.
    ENDIF
#else    
    is_ok_omp=.FALSE.
#endif


    IF (is_ok_omp) THEN
      STOP 'Open MP is not yet implemented for driver'
    ELSE
      omp_size=1
      omp_rank=0
      is_omp_root=.TRUE.
    ENDIF

  END SUBROUTINE Init_orchidee_omp


!!  =============================================================================================================================
!! SUBROUTINE:  Init_numout_omp  
!!
!>\BRIEF  Define a number for the output file specific to the omp thread. 	  
!!
!! DESCRIPTION:	   Define a number for the output file specific to the omp thread. 	  
!!
!!                
!! \n
!_ ==============================================================================================================================
  SUBROUTINE Init_numout_omp(numout)
    INTEGER, INTENT(in) :: numout
    numout_omp=numout
  END SUBROUTINE Init_numout_omp


!!  =============================================================================================================================
!! SUBROUTINE:  Init_orchidee_omp_data  
!!
!>\BRIEF	  Omp parallelisation in the coupled case. 
!!
!! DESCRIPTION:	   Omp parallelisation in the coupled case. In this routine we will define all omp variables
!!                 is_omp_root, omp_size, omp_rank, nbp_omp_para_nb, nbp_omp_para_begin, nbp_omp_para_end
!!                 nbp_omp_begin, nbp_omp_end, nbp_omp
!!
!!                
!! \n
!_ ==============================================================================================================================
  SUBROUTINE Init_orchidee_omp_data(arg_omp_size,arg_omp_rank,arg_nbp_omp,arg_offset_omp)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: arg_omp_size
    INTEGER, INTENT(IN) :: arg_omp_rank
    INTEGER, INTENT(IN) :: arg_nbp_omp
    INTEGER, INTENT(IN) :: arg_offset_omp
    
    INTEGER    :: i
    
    
    IF (arg_omp_rank==0) THEN
      is_omp_root=.TRUE.
    ELSE
      is_omp_root=.FALSE.
    ENDIF
    
#ifdef CPP_OMP
    IF (is_omp_root) THEN
        is_ok_omp=.TRUE.
    ENDIF
#else    
    is_ok_omp=.FALSE.
#endif

    IF (is_omp_root) omp_size=arg_omp_size

    CALL barrier2_omp()

     omp_rank=arg_omp_rank
    
    IF (is_omp_root) THEN 
      ALLOCATE(nbp_omp_para_nb(0:omp_size-1))
      ALLOCATE(nbp_omp_para_begin(0:omp_size-1))
      ALLOCATE(nbp_omp_para_end(0:omp_size-1))
    ENDIF
    
    CALL barrier2_omp()
    nbp_omp_para_nb(omp_rank)=arg_nbp_omp
    CALL barrier2_omp()
    
    IF (is_omp_root) THEN

      nbp_omp_para_begin(0)=1
      nbp_omp_para_end(0)=nbp_omp_para_nb(0)

      DO i=1,omp_size-1
        nbp_omp_para_begin(i)=nbp_omp_para_end(i-1)+1
        nbp_omp_para_end(i)=nbp_omp_para_begin(i)+nbp_omp_para_nb(i)-1
      ENDDO

    ENDIF

    CALL barrier2_omp()
     
    nbp_omp=nbp_omp_para_nb(omp_rank)
    nbp_omp_begin=nbp_omp_para_begin(omp_rank)
    nbp_omp_end=nbp_omp_para_end(omp_rank)
    
    offset_omp=arg_offset_omp          
    CALL Print_omp_data
    
    CALL Init_synchro_omp()
    
  END SUBROUTINE Init_orchidee_omp_data

!!  =============================================================================================================================
!! SUBROUTINE:    print_omp_data
!!
!>\BRIEF	 print specific omp parallelisation variables 
!!
!! DESCRIPTION:	  	 print specific omp parallelisation variables 
!!
!!                
!! \n
!_ ==============================================================================================================================
  SUBROUTINE print_omp_data
  IMPLICIT NONE

!$OMP CRITICAL  
  PRINT *,'--------> ORCHIDEE TASK ',omp_rank
  PRINT *,'omp_size =',omp_size
  PRINT *,'omp_rank =',omp_rank
  PRINT *,'is_omp_root =',is_omp_root
  PRINT *,'offset_omp',offset_omp
  PRINT *,'nbp_omp_para_nb =',nbp_omp_para_nb
  PRINT *,'nbp_omp_para_begin =',nbp_omp_para_begin
  PRINT *,'nbp_omp_para_end =',nbp_omp_para_end    
  PRINT *,'nbp_omp =',nbp_omp
  PRINT *,'nbp_omp_begin =',nbp_omp_begin
  PRINT *,'nbp_omp_end =',nbp_omp_end    
!$OMP END CRITICAL

  END SUBROUTINE print_omp_data

!!  =============================================================================================================================
!! SUBROUTINE:  Init_synchro_omp  
!!
!>\BRIEF    initialization of  some variables use for the synchronisation of omp threads
!!
!! DESCRIPTION:	  initialization of  some variables use for the synchronisation of omp threads
!!
!!                
!! \n
!_ ==============================================================================================================================
  SUBROUTINE Init_synchro_omp
  IMPLICIT NONE
    
    IF (is_omp_root) THEN
      ALLOCATE(proc_synchro_omp(0:omp_size-1))
      proc_synchro_omp(:)=.FALSE.

      IF ( check_all_transfert ) THEN
         ALLOCATE(omp_function(0:omp_size-1))
         omp_function(:)=-1
      ENDIF
    ENDIF
    CALL barrier2_omp()

  END SUBROUTINE Init_Synchro_omp
  
!!  =============================================================================================================================
!! SUBROUTINE:   Synchro_omp 
!!
!>\BRIEF            routine to make synchronisation of omp threads after a call to a omp routine	  
!!
!! DESCRIPTION:	  routine to make synchronisation of omp threads after a call to a omp routine
!!                add a control to check the time waited for the synchronisation. 
!!                
!! \n
!_ ==============================================================================================================================
  SUBROUTINE Synchro_omp

#ifdef CPP_PARA
    USE mpi
#endif

    IMPLICIT NONE

    INTEGER iter
    LOGICAL, PARAMETER :: check=.TRUE.
    INTEGER, PARAMETER :: iter_max=1
    INTEGER, PARAMETER :: print_iter=1
    INTEGER            :: ierr

    proc_synchro_omp(omp_rank)=.TRUE.
    CALL barrier2_omp()

    iter=0
    DO WHILE (.NOT. ALL(proc_synchro_omp))
       iter=iter+1
       IF ( mod(iter,print_iter) == 0 ) THEN
          IF (numout_omp > 0) THEN
             WRITE(numout_omp,*) "ORCHIDEE SYNCHRO OMP : iter ",iter," rank ",omp_rank," wait for ",proc_synchro_omp
          ELSE
             WRITE(*,*) "ORCHIDEE SYNCHRO OMP : iter ",iter," rank ",omp_rank," wait for ",proc_synchro_omp
          ENDIF
       ENDIF
       IF (check) THEN
          IF (iter > iter_max) THEN
             IF (numout_omp > 0) THEN
                WRITE(numout_omp,*) "TOO MUCH WAIT in Synchro_Omp !! iter ",iter," rank ",omp_rank," wait for ",proc_synchro_omp
                WRITE(numout_omp,*) "We stop here"
                WRITE(numout_omp,*) "omp_function : ",omp_function(:)
             ELSE
                WRITE(*,*) "TOO MUCH WAIT in Synchro_Omp !! iter ",iter," rank ",omp_rank," wait for ",proc_synchro_omp
                WRITE(*,*) "We stop here"
                WRITE(*,*) "omp_function : ",omp_function(:)
             ENDIF
#ifdef CPP_PARA
             CALL MPI_ABORT(MPI_COMM_ORCH, 1, ierr)
#endif     
             STOP 'Fatal error from ORCHIDEE : Synchro_Omp failed'
          ENDIF
       ENDIF
    CALL barrier2_omp()
    ENDDO
    CALL barrier2_omp()
    proc_synchro_omp(omp_rank)=.FALSE.
    CALL barrier2_omp()

   END SUBROUTINE Synchro_omp

!!  =============================================================================================================================
!! SUBROUTINE:    print_omp_function
!!
!>\BRIEF	  
!!
!! DESCRIPTION:	  
!!
!!                
!! \n
!_ ==============================================================================================================================
   SUBROUTINE print_omp_function ()

     IF ( check_all_transfert ) THEN
        CALL barrier2_omp()
        IF (numout_omp > 0) THEN
           WRITE(numout_omp,*) omp_rank,&
                " : ",omp_fct_name(omp_previous),'->',omp_fct_name(omp_function(omp_rank))
           IF (MINVAL(omp_function(:)).LT.MAXVAL(omp_function(:))) &
                WRITE(numout_omp,*) "!!! OMP ERROR : NO MORE SYNCHRO  !!!  ",omp_function(:)
        ELSE
           WRITE(*,*) omp_rank,&
                " : ",omp_fct_name(omp_previous),'->',omp_fct_name(omp_function(omp_rank))
           IF (MINVAL(omp_function(:)).LT.MAXVAL(omp_function(:))) &
                WRITE(*,*) "!!! OMP ERROR : NO MORE SYNCHRO  !!!  ",omp_function(:)
        ENDIF
        CALL barrier2_omp()
     ENDIF

  END SUBROUTINE print_omp_function


END MODULE mod_orchidee_omp_data
