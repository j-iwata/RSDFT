MODULE sweep_module

  use parallel_module
  use electron_module, only: Nfixed, Ndspin, Nspin, Nband, Nelectron
  use bz_module, only: weight_bz, Nbzsm
  use wf_module
  use cg_module, only: conjugate_gradient
  use array_bound_module, only: ML_0,ML_1,MBZ_0,MBZ_1,MSP_0,MSP_1,MB_0,MB_1
  use gram_schmidt_module
  use io_module
  use total_energy_module, only: calc_with_rhoIN_total_energy
  use fermi_module
  use subspace_diag_module
  use esp_gather_module
  use watch_module
  use hamiltonian_module
  use xc_hybrid_module, only: control_xc_hybrid, get_flag_xc_hybrid
  use io_tools_module
  use eigenvalues_module

  implicit none

  PRIVATE
  PUBLIC :: calc_sweep, init_sweep, read_sweep

  integer :: iconv_check=1
  real(8) :: Echk, Echk0
  real(8) :: tol_Echk=1.d-12
  real(8) :: tol_esp=1.d-7
  real(8) :: max_esperr
  integer :: mb_ref
  integer :: Nsweep

CONTAINS


  SUBROUTINE read_sweep
    implicit none
    call IOTools_readIntegerKeyword( "NSWEEP", Nsweep )
  END SUBROUTINE read_sweep


  SUBROUTINE init_sweep( iconv_check_in, mb_ref_in, tol_in )
    implicit none
    integer,intent(IN) :: iconv_check_in, mb_ref_in
    real(8),intent(IN) :: tol_in
    iconv_check = iconv_check_in
    mb_ref = mb_ref_in
    select case( iconv_check )
    case( 1 )
       tol_Echk = tol_in
    case( 2 )
       tol_esp = tol_in
    case( 3 )
       tol_esp = tol_in
    end select
  END SUBROUTINE init_sweep


  SUBROUTINE calc_sweep( disp_switch, ierr_out, Diter_in, outer_loop_info )
    implicit none
    logical,intent(IN)  :: disp_switch
    integer,intent(OUT) :: ierr_out
    integer,optional,intent(IN)  :: Diter_in
    character(*),optional,intent(IN) :: outer_loop_info
    integer :: iter,s,k,n,m,iflag_hybrid,ierr,Diter
    logical :: flag_exit, flag_conv
    logical :: flag_end, flag_end1, flag_end2
    character(40) :: chr_iter
    character(22) :: add_info
    type(time) :: etime
    type(eigv) :: eval
    logical,external :: exit_program

    Diter = 0
    if ( present(Diter_in) ) then
       Diter = Diter_in
    else
       Diter = Nsweep
    end if

    if ( Diter <= 0 ) return

    call write_border( 0, "" )
    call write_border( 0, " SWEEP START ----------" )

    flag_end  = .false.
    flag_exit = .false.
    flag_conv = .false.
    Echk      = 0.0d0
    ierr_out  = 0
    add_info  = "" ; if ( present(outer_loop_info) ) add_info=outer_loop_info

    allocate( esp0(Nband,Nbzsm,Nspin) ) ; esp0=0.0d0

    call get_flag_xc_hybrid( iflag_hybrid )

    select case( iflag_hybrid )
    case(0,1,3)

       if ( iflag_hunk >= 1 ) then
          do s=MSP_0,MSP_1
          do k=MBZ_0,MBZ_1
             do m=MB_0,MB_1,MB_d
                n=min(m+MB_d-1,MB_1)
                call hamiltonian &
                     (k,s,unk(ML_0,m,k,s),hunk(ML_0,m,k,s),ML_0,ML_1,m,n)
             end do
          end do
          end do
       end if

    case( 2 )

       if ( iflag_hunk >= 1 ) then
          call control_xc_hybrid(0)
          allocate( workwf(ML_0:ML_1,MB_d) ) ; workwf=0.0d0
          do s=MSP_0,MSP_1
          do k=MBZ_0,MBZ_1
             do m=MB_0,MB_1,MB_d
                n=min(m+MB_d-1,MB_1)
                workwf(:,1:n-m+1)=hunk(:,m:n,k,s)
                call hamiltonian &
                     (k,s,unk(ML_0,m,k,s),hunk(ML_0,m,k,s),ML_0,ML_1,m,n)
                hunk(:,m:n,k,s)=hunk(:,m:n,k,s)+workwf(:,1:n-m+1)
             end do ! m
          end do ! k
          end do ! s
          deallocate( workwf )
       end if

       call control_xc_hybrid(1)

    end select

    call calc_with_rhoIN_total_energy( Echk )

    do iter=1,Diter

       write(chr_iter,'(" sweep_iter=",i4,1x,a)') iter, add_info
       call write_border( 0, chr_iter(1:len_trim(chr_iter)) )

       call init_time_watch( etime )

       Echk0=Echk
       esp0 =esp
       do s=MSP_0,MSP_1
       do k=MBZ_0,MBZ_1

          call conjugate_gradient( ML_0,ML_1, Nband, k,s &
                                 ,unk(ML_0,1,k,s), esp(1,k,s), res(1,k,s) )

          call gram_schmidt(1,Nband,k,s)

          call subspace_diag(k,s)

       end do
       end do

       call esp_gather(Nband,Nbzsm,Nspin,esp)

#ifdef _DRSDFT_
       call mpi_bcast( unk, size(unk), MPI_REAL8, 0, comm_fkmb, ierr )
#else
       call mpi_bcast( unk, size(unk), MPI_COMPLEX16, 0, comm_fkmb, ierr )
#endif
       call mpi_bcast( esp, size(esp), MPI_REAL8, 0, comm_fkmb, ierr )

       call calc_fermi(iter,Nfixed,Nband,Nbzsm,Nspin,Nelectron,Ndspin &
                      ,esp,weight_bz,occ,disp_switch)

       call calc_with_rhoIN_total_energy( Echk )

       call conv_check( iter, res, flag_conv )
       call global_watch( .false., flag_end1 )
       flag_end2 = exit_program()
       flag_end  = ( flag_end1 .or. flag_end2 )
       flag_exit = (flag_end.or.flag_conv.or.(iter==Diter))

       if ( disp_switch ) call write_info_sweep
       call construct_eigenvalues( Nband, Nbzsm, Nspin, esp, eval )
       if ( myrank == 0 ) call write_eigenvalues( eval )

       call calc_time_watch( etime )
       if ( disp_switch ) then
          write(*,*)
          write(*,'(1x,"time(sweep)=",f10.3,"(rank0)",f10.3,"(min)" &
               ,f10.3,"(max)")') etime%t0, etime%tmin, etime%tmax
          write(*,*)
       end if

       call write_data( disp_switch, flag_exit )

       if ( flag_exit ) exit

    end do ! iter

    if ( myrank == 0 ) call write_info_esp_wf( 2 )

    deallocate( esp0 )

    ierr_out = iter

    if ( flag_end1 ) then
       ierr_out = -1
       if ( myrank == 0 ) write(*,*) "flag_end=",flag_end
       return
    end if

    if ( flag_end2 ) then
       ierr_out = -3
       if ( myrank == 0 ) write(*,*) "flag_end=",flag_end
       return
    end if

    if ( iter > Diter ) then
       ierr_out = -2
       if ( myrank == 0 ) write(*,*) "sweep not converged"
    end if

    call gather_wf

    call write_border( 0, " SWEEP END ----------" )
    call write_border( 0, "" )

  END SUBROUTINE calc_sweep


  SUBROUTINE conv_check( iter, res, flag_conv )
    implicit none
    integer,intent(IN)  :: iter
    real(8),intent(INOUT) :: res(:,:,:)
    logical,intent(OUT) :: flag_conv
    call write_border( 1, " conv_check(start)" )
    select case( iconv_check )
    case( 1 )
         call conv_check_1( iter, flag_conv )
    case( 2 )
         call conv_check_2( flag_conv )
    case( 3 )
         call conv_check_3( res, flag_conv )
    end select
    call write_border( 1, " conv_check(end)" )
  END SUBROUTINE conv_check

  SUBROUTINE conv_check_1( iter, flag_conv )
    implicit none
    integer,intent(IN)  :: iter
    logical,intent(OUT) :: flag_conv
    flag_conv=.false.
    if ( abs(Echk-Echk0) < tol_Echk ) flag_conv=.true.
  END SUBROUTINE conv_check_1

  SUBROUTINE conv_check_2( flag_conv )
    implicit none
    logical,intent(OUT) :: flag_conv
    integer :: ierr
    real(8) :: err0
    max_esperr = maxval( abs(  esp(1:mb_ref,MBZ_0,MSP_0:MSP_1) &
                             -esp0(1:mb_ref,MBZ_0,MSP_0:MSP_1) ) )
    call mpi_allreduce(max_esperr,err0,1,MPI_REAL8,MPI_MAX,comm_spin,ierr)
    call mpi_allreduce(err0,max_esperr,1,MPI_REAL8,MPI_MAX,comm_bzsm,ierr)
    flag_conv = .false.
    if ( max_esperr < tol_esp ) flag_conv = .true.
  END SUBROUTINE conv_check_2

  SUBROUTINE conv_check_3( res, flag_conv )
    implicit none
    real(8),intent(OUT) :: res(:,:,:)
    logical,intent(OUT) :: flag_conv
    integer :: ierr
    real(8) :: err0
    res=0.0d0
    res=abs(esp-esp0)
    max_esperr = maxval( res(1:mb_ref,:,:) )
    flag_conv = .false.
    if ( max_esperr < tol_esp ) flag_conv = .true.
    where( res < tol_esp ) res=-1.0d0
  END SUBROUTINE conv_check_3


  SUBROUTINE write_info_sweep
    implicit none
    integer :: s,k,n
    call write_border( 1, " write_info_sweep(start)" )
!    write(*,'(a4,a6,a20,2a13,1x)') &
!         "k","n","esp(n,k,s)","esp_err  ","occ(n,k,s)  "
!    do k=1,Nbzsm
!    do n=max(1,nint(Nelectron/2)-5),min(nint(Nelectron/2)+5,Nband)
!       write(*,'(i4,i6,2(f20.15,2g13.5,1x))') k,n &
!            ,(esp(n,k,s),esp(n,k,s)-esp0(n,k,s),occ(n,k,s),s=1,Nspin)
!    end do
!    end do
    if ( iconv_check == 1 ) then
       write(*,'(/,1x,"Echk,dif/tol =",g18.10,2x,g12.5," /",es12.5)') &
            Echk, Echk-Echk0, tol_Echk
    else
       write(*,'(/,1x,"max_esperr/tol, mb_ref =",g12.5," /",es12.5,i7)') &
            max_esperr,tol_esp,mb_ref
    end if
    call write_esp_wf
!    write(*,*) "sum(occ)=",(sum(occ(:,:,s)),s=1,Nspin)
    call write_border( 1, " write_info_sweep(end)" )
  END SUBROUTINE write_info_sweep


END MODULE sweep_module
