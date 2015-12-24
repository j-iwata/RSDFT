MODULE subspace_diag_sl_module

  use scalapack_module
  use subspace_mate_sl_module
  use subspace_solv_sl_module
  use subspace_rotv_sl_module
  use subspace_diag_variables
  use watch_module

  implicit none

  PRIVATE
  PUBLIC :: subspace_diag_sl

CONTAINS


  SUBROUTINE subspace_diag_sl(k,s,disp_switch)
    implicit none
    logical,intent(IN) :: disp_switch
    integer,intent(IN) :: k,s
    real(8) :: ct(5),et(5),ct0,et0,ct1,et1

    call write_border( 1, " subspace_diag_sl(start)" )

    ct(:)=0.d0
    et(:)=0.d0

    call prep_scalapack( MB_diag )

    allocate( Hsub(LLD_R,LLD_C) )
    Hsub=zero

    call watch(ct0,et0)
    call subspace_mate_sl(k,s)
    call watch(ct1,et1) ; ct(1)=ct1-ct0 ; et(1)=et1-et0

    allocate( Vsub(LLD_R,LLD_C) )
    Vsub=zero

    call watch(ct0,et0)
    call subspace_solv_sl(k,s)
    call watch(ct1,et1) ; ct(2)=ct1-ct0 ; et(2)=et1-et0

    deallocate( Hsub )

    call watch(ct0,et0)
    call subspace_rotv_sl(k,s)
    call watch(ct1,et1) ; ct(3)=ct1-ct0 ; et(3)=et1-et0

    deallocate( Vsub )

!    if ( disp_switch ) then
!       write(*,'(1x,"time(diag_mate1)",2g10.3)') ct(1),et(1)
!!       write(*,'(1x,"time(diag_mate2)",2g10.3)') ct(4),et(4)
!       write(*,'(1x,"time(scalapack)",2g10.3)') ct(2),et(2)
!       write(*,'(1x,"time(diag_rotv1)",2g10.3)') ct(3),et(3)
!!       write(*,'(1x,"time(diag_rotv2)",2g10.3)') ct(5),et(5)
!    end if

    call write_border( 1, " subspace_diag_sl(end)" )

  END SUBROUTINE subspace_diag_sl


END MODULE subspace_diag_sl_module
