MODULE xc_hf_module

  use xc_hybrid_module, only: iflag_hybrid
  use fock_module, only: UpdateWF_fock
  use wf_module, only: unk, occ, hunk, allocate_work_wf
  use parallel_module
  use watch_module

  implicit none

  PRIVATE
  PUBLIC :: init_xc_hf, calc_xc_hf

  real(8) :: E_exchange
  real(8) :: dV
  integer :: SYStype
  integer :: ML_0 ,ML_1
  integer :: MSP_0,MSP_1
  integer :: MBZ_0,MBZ_1
  integer :: MB_0 ,MB_1
  logical :: flag_init=.false.

CONTAINS


  SUBROUTINE init_xc_hf( ML_0_in,ML_1_in, MSP_0_in,MSP_1_in &
       ,MBZ_0_in,MBZ_1_in, MB_0_in,MB_1_in, SYStype_in, dV_in )

    implicit none

    integer,intent(IN) :: SYStype_in
    integer,intent(IN) :: ML_0_in,ML_1_in, MSP_0_in,MSP_1_in
    integer,intent(IN) :: MB_0_in,MB_1_in, MBZ_0_in,MBZ_1_in
    real(8),intent(IN) :: dV_in

    if ( flag_init ) return

    call allocate_work_wf( 2 )

    ML_0       = ML_0_in
    ML_1       = ML_1_in
    MSP_0      = MSP_0_in
    MSP_1      = MSP_1_in
    MBZ_0      = MBZ_0_in
    MBZ_1      = MBZ_1_in
    MB_0       = MB_0_in
    MB_1       = MB_1_in
    SYStype    = SYStype_in
    dV         = dV_in
    E_exchange = 0.0d0
    flag_init  = .true.

  END SUBROUTINE init_xc_hf


  SUBROUTINE calc_xc_hf( Ex )
    implicit none
    real(8),intent(OUT) :: Ex
    real(8) :: sum0,sum1,c
    integer :: s,k,n,i,ierr

    if ( .not. flag_init ) then
       write(*,*) "Please call init_xc_hf first"
       stop "stop@calc_xc_hf"
    end if

    if ( iflag_hybrid == 1 ) then

!       call watcht(disp_switch_parallel,"",0)
       call UpdateWF_fock( SYStype )
!       call watcht(disp_switch_parallel,"UpdateWF_fock",1)

       sum0=0.0d0
       do s=MSP_0,MSP_1
       do k=MBZ_0,MBZ_1
       do n=MB_0 ,MB_1

          if ( abs(occ(n,k,s)) < 1.d-10 ) cycle

          c = occ(n,k,s)

          do i=ML_0,ML_1
#ifdef _DRSDFT_
             sum0 = sum0 + c*unk(i,n,k,s)*hunk(i,n,k,s)
#else
             sum0 = sum0 + c*conjg(unk(i,n,k,s))*hunk(i,n,k,s)
#endif
          end do ! i

       end do ! n
       end do ! k
       end do ! s

       call mpi_allreduce(sum0,sum1,1,mpi_real8,mpi_sum,comm_grid,ierr)
       call mpi_allreduce(sum1,sum0,1,mpi_real8,mpi_sum,comm_band,ierr)
       call mpi_allreduce(sum0,sum1,1,mpi_real8,mpi_sum,comm_bzsm,ierr)
       call mpi_allreduce(sum1,sum0,1,mpi_real8,mpi_sum,comm_spin,ierr)

       E_exchange = 0.5d0*sum0*dV

    end if

    Ex = E_exchange

  END SUBROUTINE calc_xc_hf


END MODULE xc_hf_module
