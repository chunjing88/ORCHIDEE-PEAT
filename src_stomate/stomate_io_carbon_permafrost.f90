!< $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_stomate/stomate_io.f90 $ 
!< $Date: 2016-06-17 13:26:43 +0200 (Fri, 17 Jun 2016) $
!< $Author: albert.jornet $
!< $Revision: 3564 $
! IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
MODULE stomate_io_carbon_permafrost
  !---------------------------------------------------------------------
  !-
  !- 
  !-
  !---------------------------------------------------------------------
  USE netcdf
  USE defprec
  USE stomate_data
  USE constantes
  USE constantes_soil
  USE mod_orchidee_para
  USE ioipsl_para 
  USE utils      ! nccheck
#ifdef CPP_PARA
  USE mpi
#endif
  !-
  IMPLICIT NONE
  !-
  PRIVATE
  PUBLIC stomate_io_carbon_permafrost_write, stomate_io_carbon_permafrost_read
  !-
  ! TO CHECK, stomate_finalize also uses this type of var
  INTEGER,PARAMETER                              :: r_typ = NF90_REAL8   !! Specify data format (server dependent)
  !-
CONTAINS
  !-
  !===
  !-
  !! ================================================================================================================================
  !! SUBROUTINE 	: stomate_io_carbon_permafrost_write 
  !!
  !>\BRIEF        Writes stomate permafrost carbon data into a netcdf file
  !!
  !! DESCRIPTION  : It writes into a netcdf files in parallel mode all necessary
  !!                 variables required for spinup (forecesoil)
  !!                
  !!                
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE stomate_io_carbon_permafrost_write (Cforcing_permafrost_name, & 
                        nbp_glo,            start_px,           length_px,      nparan,     nbyear, &
                        index_g,            zz_deep,            zz_coef_deep, &
                        clay,               depth_organic_soil, lalo, &
                        snowdz_2pfcforcing, snowrho_2pfcforcing, soilcarbon_input_2pfcforcing, &
                        tsurf_2pfcforcing,  pb_2pfcforcing,     snow_2pfcforcing, &
                        tprof_2pfcforcing,  fbact_2pfcforcing,  veget_max_2pfcforcing, &
                        rprof_2pfcforcing,  hslong_2pfcforcing, root_exudate_2pfcforcing )
    
    


    CHARACTER(LEN=100), INTENT(in)              :: Cforcing_permafrost_name !! Name of permafrost forcing file
    INTEGER(i_std), INTENT(in)                  :: nbp_glo !nbp_glo is the number of global continental points
    INTEGER(i_std), INTENT(in)                  :: start_px ! Start land point/pixex respect to nbp_glo
    INTEGER(i_std), INTENT(in)                  :: length_px ! Length of lands point/pixel to write
    INTEGER(i_std), INTENT(in)                  :: nparan ! Number of forcesoil timesteps  
    INTEGER(i_std), INTENT(in)                  :: nbyear ! Number of years saved for carbon spinup 
    INTEGER(i_std),DIMENSION(:),INTENT(in)      :: index_g             !! Indices of the terrestrial pixels only (unitless)
    REAL(r_std), DIMENSION(:),   INTENT (in)    :: zz_deep           !! deep vertical profile
    REAL(r_std), DIMENSION(:),   INTENT (in)    :: zz_coef_deep      !! deep vertical profile
    REAL(r_std), DIMENSION(:), INTENT(in)       :: clay                   !! Clay fraction of soil (0-1, unitless), parallel 
    REAL(r_std), DIMENSION(:),   INTENT (in)    :: depth_organic_soil !! how deep is the organic soil?
    REAL(r_std), DIMENSION(:,:),INTENT(in)      :: lalo              !! Geographical coordinates (latitude,longitude) 
    REAL(r_std),DIMENSION(:,:,:), INTENT(in)    :: snowdz_2pfcforcing
    REAL(r_std),DIMENSION(:,:,:), INTENT(in)    :: snowrho_2pfcforcing
    REAL(r_std),DIMENSION(:,:,:,:), INTENT(in)  :: soilcarbon_input_2pfcforcing
    REAL(r_std),DIMENSION(:,:), INTENT(in )     :: tsurf_2pfcforcing
    REAL(r_std),DIMENSION(:,:), INTENT(in)      :: pb_2pfcforcing
    REAL(r_std),DIMENSION(:,:), INTENT(in)      :: snow_2pfcforcing
    REAL(r_std),DIMENSION(:,:,:,:), INTENT(in)  :: tprof_2pfcforcing
    REAL(r_std),DIMENSION(:,:,:,:), INTENT(in)  :: fbact_2pfcforcing
    REAL(r_std),DIMENSION(:,:,:,:), INTENT(in)  :: hslong_2pfcforcing
    REAL(r_std),DIMENSION(:,:,:), INTENT(in)    :: veget_max_2pfcforcing
    REAL(r_std),DIMENSION(:,:,:), INTENT(in)    :: rprof_2pfcforcing
!!!qcj++ peatland
    REAL(r_std),DIMENSION(:,:,:), INTENT(in)    :: root_exudate_2pfcforcing
   
    ! Local Variables
    INTEGER(i_std)                              :: ier, n_directions, i
    INTEGER(i_std)                              :: start(1), ncount(1), start_2d(2), ncount_2d(2), inival, endval 
    INTEGER(i_std)                              :: start_4d(4), ncount_4d(4), start_3d(3), ncount_3d(3)
    INTEGER(i_std),DIMENSION(10)                :: d_id                     !! List each netcdf dimension
    INTEGER(i_std)                              :: vid                      !! Variable identifer of netCDF (unitless)
    INTEGER(i_std)                              :: Cforcing_permafrost_id   !! Permafrost file identifer 
   
    ! Create file 
#ifdef CPP_PARA
    ier = NF90_CREATE (TRIM(Cforcing_permafrost_name),IOR(NF90_NETCDF4,NF90_MPIIO), &
            Cforcing_permafrost_id, comm=MPI_COMM_ORCH, info=MPI_INFO_NULL)
#else
    ier = NF90_CREATE (TRIM(Cforcing_permafrost_name),NF90_NETCDF4, &
            Cforcing_permafrost_id)
#endif
    IF (ier /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'stomate_finalize', &
             &        'PROBLEM creating Cforcing_permafrost file', &
             &        NF90_STRERROR(ier),'')
     END IF


     ! Add variable attribute
     ! Note ::nbp_glo is the number of global continental points
     CALL nccheck( NF90_PUT_ATT (Cforcing_permafrost_id,NF90_GLOBAL, &
          &                           'kjpindex',REAL(nbp_glo,r_std)))
     CALL nccheck( NF90_PUT_ATT (Cforcing_permafrost_id,NF90_GLOBAL, &
          &                           'nparan',REAL(nparan,r_std)))
     CALL nccheck( NF90_PUT_ATT (Cforcing_permafrost_id,NF90_GLOBAL, &
          &                           'nbyear',REAL(nbyear,r_std)))

     ! Add new dimension, variables values from USE
     CALL nccheck( NF90_DEF_DIM (Cforcing_permafrost_id,'points',nbp_glo,d_id(1)))
     CALL nccheck( NF90_DEF_DIM (Cforcing_permafrost_id,'carbtype',ncarb,d_id(2)))
     CALL nccheck( NF90_DEF_DIM (Cforcing_permafrost_id,'vegtype',nvm,d_id(3)))
     CALL nccheck( NF90_DEF_DIM (Cforcing_permafrost_id,'level',ndeep,d_id(4)))
     CALL nccheck( NF90_DEF_DIM (Cforcing_permafrost_id,'time_step',nparan*nbyear,d_id(5)))
     n_directions=2
     CALL nccheck( NF90_DEF_DIM (Cforcing_permafrost_id,'direction',n_directions,d_id(6)))


     ! Add new variable
     CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'points', r_typ,d_id(1),vid))
     CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'carbtype', r_typ,d_id(2),vid))
     CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'vegtype', r_typ,d_id(3),vid))
     CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'level', r_typ,d_id(4),vid))
     CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'time_step',r_typ,d_id(5),vid))
     CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'direction',r_typ,d_id(6),vid))
     CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'index',r_typ,d_id(1),vid))
     CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'clay',r_typ,d_id(1),vid))
     CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'lalo',      r_typ, &
          (/ d_id(1), d_id(6) /),vid))
     !--time-invariant
     CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'zz_deep',r_typ,d_id(4),vid))
     CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'zz_coef_deep',r_typ,d_id(4),vid))
     CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'z_organic',r_typ,d_id(1),vid))
     !--3layers snow 
     CALL nccheck( NF90_DEF_VAR(Cforcing_permafrost_id,'snowdz',r_typ,(/ d_id(1),d_id(2),d_id(5) /),vid))
     CALL nccheck( NF90_DEF_VAR(Cforcing_permafrost_id,'snowrho',r_typ,(/ d_id(1),d_id(2),d_id(5) /),vid))
     !--time-varying 
     CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'soilcarbon_input',r_typ, &
          &                        (/ d_id(1),d_id(2),d_id(3),d_id(5) /),vid))
     CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'pb',r_typ, & 
          &                        (/ d_id(1),d_id(5) /),vid))
     CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'snow',r_typ, &
          &                        (/ d_id(1),d_id(5) /),vid))
     CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'tprof',r_typ, &
          &                        (/ d_id(1),d_id(4),d_id(3),d_id(5) /),vid))
     CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'fbact',r_typ, &
          &                        (/ d_id(1),d_id(4),d_id(3),d_id(5) /),vid))
     CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'hslong',r_typ, &
          &                        (/ d_id(1),d_id(4),d_id(3),d_id(5) /),vid))
     CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'veget_max',r_typ, &
          &                        (/ d_id(1),d_id(3),d_id(5) /),vid))
     CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'rprof',r_typ, &
          &                        (/ d_id(1),d_id(3),d_id(5) /),vid))
     CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'tsurf',r_typ, &
          &                        (/ d_id(1),d_id(5) /),vid))
!!!qcj++ peatland
     CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'root_exudate',r_typ, &
           &                        (/ d_id(1),d_id(3),d_id(5) /),vid))
     CALL nccheck( NF90_ENDDEF (Cforcing_permafrost_id))

     ! Write data
     start=(/ start_px /)
     ncount=(/ length_px /)
     inival=start_px
     endval=start_px + length_px !length_px_end(mpi_rank)
     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'points',vid) )
     CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
     CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, &
          &  (/(REAL(i,r_std),i=inival,endval)/), &
          &  start=start, count=ncount) )

     ! no point to make parallel calls
     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'carbtype',vid))
     CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
     CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, &
          &                        (/(REAL(i,r_std),i=1,ncarb)/)))
     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'vegtype',vid))
     CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
     CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, &
          &                            (/(REAL(i,r_std),i=1,nvm)/)))
     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'level',vid))
     CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
     CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, &
          &                        (/(REAL(i,r_std),i=1,ndeep)/)))
     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'time_step',vid))
     CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
     CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, &
          &                       (/(REAL(i,r_std),i=1,nparan*nbyear)/)))
     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'index',vid))
     CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
     CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, REAL(index_g,r_std) ))

     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'zz_deep',vid))
     CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
     CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, zz_deep ))

     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'zz_coef_deep',vid))
     CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
     CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, zz_coef_deep ))

     ! Parallel writes
     start=(/ start_px /)
     ncount=(/ length_px /)
     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'clay',vid))
     CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
     CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, clay, start=start, count=ncount  ))

     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'z_organic',vid))
     CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
     CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid,depth_organic_soil, start=start, count=ncount))

     start_2d=(/ start_px,1 /)
     ncount_2d=(/ length_px,2 /)
     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'lalo',vid))
     CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
     CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, lalo, start=start_2d, count=ncount_2d ))

     ! putting 3 snow layers
     start_3d=(/ start_px,1,1 /)
     ncount_3d=(/ length_px,nsnow,nparan*nbyear /)
     CALL nccheck( NF90_INQ_VARID(Cforcing_permafrost_id,'snowdz',vid))
     CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
     CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, snowdz_2pfcforcing, start=start_3d ,count=ncount_3d ))

     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'snowrho',vid))
     CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
     CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, snowrho_2pfcforcing, &
                    start=start_3d ,count=ncount_3d ))

     start_4d=(/ start_px,1,1,1 /)
     ncount_4d=(/ length_px,ncarb,nvm,nparan*nbyear /)
     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'soilcarbon_input',vid))
     CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
     CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid,soilcarbon_input_2pfcforcing, &
                    start=start_4d ,count=ncount_4d ))

!!!qcj++ peatland
     start_3d=(/ start_px,1,1 /)
     ncount_3d=(/ length_px,nvm,nparan*nbyear /)
     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'root_exudate',vid))
     CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
     CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid,root_exudate_2pfcforcing, start=start_3d,count=ncount_3d ))

     start_2d=(/ start_px,1 /)
     ncount_2d=(/ length_px,nparan*nbyear /)
     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'tsurf',vid))
     CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
     CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, tsurf_2pfcforcing, start=start_2d, count=ncount_2d ))

     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'pb',vid))
     CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
     CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, pb_2pfcforcing, start=start_2d, count=ncount_2d ))

     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'snow',vid))
     CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
     CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, snow_2pfcforcing, start=start_2d, count=ncount_2d ))

     start_4d=(/ start_px,1,1,1 /)
     ncount_4d=(/ length_px,ndeep,nvm,nparan*nbyear /)
     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'tprof',vid))
     CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
     CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, tprof_2pfcforcing ,start=start_4d, count=ncount_4d))

     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'fbact',vid))
     CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
     CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, fbact_2pfcforcing,start=start_4d, count=ncount_4d ))

     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'hslong',vid))
     CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
     CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, hslong_2pfcforcing,start=start_4d, count=ncount_4d ))

     start_3d=(/ start_px,1,1 /)
     ncount_3d=(/ length_px,nvm,nparan*nbyear /)
     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'veget_max',vid))
     CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
     CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid,veget_max_2pfcforcing, start=start_3d, count=ncount_3d ))

     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'rprof',vid))
     CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
     CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, rprof_2pfcforcing, start=start_3d, count=ncount_3d ))

     ! Finish netcdf file management
     CALL nccheck( NF90_CLOSE (Cforcing_permafrost_id) )

  END SUBROUTINE stomate_io_carbon_permafrost_write
  !-
  !===
  !-
  SUBROUTINE stomate_io_carbon_permafrost_read (Cforcing_permafrost_name,   &
                nparan,             nbyear,     start_px,   length_px,      &
                soilcarbon_input,   pb,         snow,       tsurf,          & 
                tprof,              fbact,      hslong,     rprof,          &
                lalo,               snowdz,     snowrho,    veget_max, root_exudate )

     ! Input Variables
     CHARACTER(LEN=100), INTENT(in)               :: Cforcing_permafrost_name !! Name of permafrost forcing file
     INTEGER(i_std), INTENT(in)                   :: start_px ! Start land point/pixex respect to nbp_glo
     INTEGER(i_std), INTENT(in)                   :: length_px ! Length of lands point/pixel to write
     INTEGER(i_std), INTENT(in)                   :: nparan, nbyear 

     ! Output variables
     REAL(r_std), DIMENSION(:,:,:,:), INTENT(out) :: soilcarbon_input
!!!qcj++ peatland
     REAL(r_std), DIMENSION(:,:,:), INTENT(out) :: root_exudate

     REAL(r_std), DIMENSION(:,:), INTENT(out)     :: pb 
     REAL(r_std), DIMENSION(:,:), INTENT(out)     :: snow
     REAL(r_std), DIMENSION(:,:), INTENT(out)     :: tsurf
     REAL(r_std), DIMENSION(:,:,:,:), INTENT(out) :: tprof
     REAL(r_std), DIMENSION(:,:,:,:), INTENT(out) :: fbact 
     REAL(r_std), DIMENSION(:,:,:,:), INTENT(out) :: hslong
     REAL(r_std), DIMENSION(:,:,:), INTENT(out)   :: rprof
     REAL(r_std), DIMENSION(:,:), INTENT(out)     :: lalo
     REAL(r_std), DIMENSION(:,:,:), INTENT(out)   :: snowdz
     REAL(r_std), DIMENSION(:,:,:), INTENT(out)   :: snowrho 
     REAL(r_std), DIMENSION(:,:,:), INTENT(out)   :: veget_max 

     ! Local Variables
     INTEGER(i_std)                               :: start_2d(2), count_2d(2) 
     INTEGER(i_std)                               :: start_4d(4), count_4d(4), start_3d(3), count_3d(3)
     INTEGER(i_std)                               :: v_id                      !! Variable identifer of netCDF (unitless)
     INTEGER(i_std)                               :: Cforcing_permafrost_id   !! Permafrost file identifer 
     !-
     ! Open FORCESOIL's forcing file to read some basic info (dimensions, variable ID's)
     ! and allocate variables.
     !-
#ifdef CPP_PARA
     CALL nccheck( NF90_OPEN (TRIM(Cforcing_permafrost_name),IOR(NF90_NOWRITE, NF90_MPIIO),Cforcing_permafrost_id, &
                 & comm = MPI_COMM_ORCH, info = MPI_INFO_NULL ))
#else
     CALL nccheck( NF90_OPEN (TRIM(Cforcing_permafrost_name),NF90_NOWRITE,Cforcing_permafrost_id))
#endif
    
     start_4d = (/ start_px, 1, 1, 1 /)
     count_4d = (/ length_px, ncarb, nvm, nparan*nbyear /)
     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'soilcarbon_input',v_id))
     CALL nccheck( NF90_GET_VAR   (Cforcing_permafrost_id,v_id,soilcarbon_input,  &
                    &  start = start_4d, count = count_4d ))

     start_2d=(/ start_px, 1 /)
     count_2d=(/ length_px, nparan*nbyear /)
     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'pb',v_id ))
     CALL nccheck( NF90_GET_VAR (Cforcing_permafrost_id,v_id,pb, &
                & start=start_2d, count=count_2d))

     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'snow',v_id))
     CALL nccheck( NF90_GET_VAR (Cforcing_permafrost_id,v_id,snow, &
                & start=start_2d, count=count_2d))

     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'tsurf',v_id))
     CALL nccheck( NF90_GET_VAR (Cforcing_permafrost_id,v_id,tsurf, &
                & start=start_2d, count=count_2d))

     start_4d=(/ start_px,1,1,1 /)
     count_4d=(/ length_px,ndeep,nvm,nparan*nbyear /)
     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'tprof',v_id))
     CALL nccheck( NF90_GET_VAR (Cforcing_permafrost_id,v_id,tprof, &
                & start=start_4d, count=count_4d))

     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'fbact',v_id))
     CALL nccheck( NF90_GET_VAR (Cforcing_permafrost_id,v_id,fbact, &
                & start=start_4d, count=count_4d))

     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'hslong',v_id))
     CALL nccheck( NF90_GET_VAR (Cforcing_permafrost_id,v_id,hslong, &
                & start=start_4d, count=count_4d))

     start_3d=(/ start_px,1,1 /)
     count_3d=(/ length_px,nvm,nparan*nbyear /)
     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'veget_max',v_id))
     CALL nccheck( NF90_GET_VAR (Cforcing_permafrost_id,v_id,veget_max, &
                & start=start_3d, count=count_3d))

     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'rprof',v_id))
     CALL nccheck( NF90_GET_VAR (Cforcing_permafrost_id,v_id,rprof, &
                & start=start_3d, count=count_3d))

!!!qcj++ peatland
     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'root_exudate',v_id))
     CALL nccheck( NF90_GET_VAR (Cforcing_permafrost_id,v_id,root_exudate, &
                & start=start_3d, count=count_3d))

     start_2d=(/ start_px, 1 /)
     count_2d=(/ length_px, 2 /)
     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'lalo',v_id))
     CALL nccheck( NF90_GET_VAR (Cforcing_permafrost_id,v_id,lalo, &
                & start=start_2d, count=count_2d))

     start_3d=(/ start_px,1,1 /)
     count_3d=(/ length_px,nsnow,nparan*nbyear /)
     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'snowdz',v_id))
     CALL nccheck( NF90_GET_VAR (Cforcing_permafrost_id,v_id,snowdz, &
                & start=start_3d, count=count_3d))

     CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'snowrho',v_id))
     CALL nccheck( NF90_GET_VAR (Cforcing_permafrost_id,v_id,snowrho, &
                & start=start_3d, count=count_3d))
     !- Close Netcdf carbon permafrost file reference
     CALL nccheck( NF90_CLOSE (Cforcing_permafrost_id))

  END SUBROUTINE stomate_io_carbon_permafrost_read

END MODULE stomate_io_carbon_permafrost
