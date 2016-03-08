MODULE dos_module

  use bz_module, only: weight_bz
  use parallel_module, only: myrank
  use io_tools_module
  use wf_module, only: esp, occ, res, write_info_esp_wf
  use sweep_module

  implicit none

  PRIVATE
  PUBLIC :: read_param_dos, calc_dos

  character(16),PUBLIC,parameter :: file_kpt_list_dos = "kpt_list_dos"
  character(16),parameter :: file_esp = "dos_eigv"

  integer :: max_loop_dos = 100
  integer :: max_ref_band = 0
  real(8) :: esp_conv_tol = 1.d-5

  integer,parameter :: unit_esp=70
  integer,parameter :: unit_kpt=71

CONTAINS


  SUBROUTINE read_param_dos
    implicit none
    logical :: exist_keyword
    integer :: unit, i
    character(8) :: cbuf
    include 'mpif.h'
    call IOTools_findKeyword( "DOS", exist_keyword, unit_out=unit )
    if ( exist_keyword ) then
       backspace(unit)
       read(unit,*) cbuf, max_ref_band, esp_conv_tol, max_loop_dos
    end if
    call MPI_BCAST( max_loop_dos, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, i )
    call MPI_BCAST( max_ref_band, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, i )
    call MPI_BCAST( esp_conv_tol, 1, MPI_REAL8, 0, MPI_COMM_WORLD, i )
  END SUBROUTINE read_param_dos


  SUBROUTINE calc_dos( ierr )

    implicit none
    integer,intent(OUT) :: ierr
    logical :: disp_sw, flag
    integer :: i,n,k,s,nb,nk,ns
    character(16) :: message

    call write_border( 0, " calc_dos(start)" )

    call read_param_dos

    ierr=0

    call check_disp_switch( disp_sw, 0 )

    nb = size( esp, 1 )
    nk = size( esp, 2 )
    ns = size( esp, 3 )

    if ( max_ref_band == 0 ) then
       max_ref_band = ( nb + nint(sum(occ)/2) )/2
       if ( disp_sw ) write(*,*) "max_ref_band is replaced to ",max_ref_band
    end if

    do s=1,ns
    do k=1,nk
       occ(:,k,s) = weight_bz(k) * (2.0d0/ns)
    end do
    end do

    call init_sweep( 3, max_ref_band, esp_conv_tol )

    call calc_sweep( disp_sw, ierr, max_loop_dos, "(DOS)" )

    if ( myrank == 0 ) then

!       open( unit_esp, FILE=file_esp, POSITION="append" )
!       write(unit_esp,*) "Nband,Nbzsm,Nspin= ",nb,nk,ns
!       do s=1,ns
!       do k=1,nk
!          if ( all( nint(res(1:max_ref_band,k,s)) == -1 ) ) then
!             message="Converged"
!          else
!             message="NotConverged"
!          end if
!          write(unit_esp,'(1x,"spin= ",i1,3x,"kpt= ",3f12.8,3x,a12)') &
!               s,kbb(1:3,k),message
!          do n=1,nb
!             write(unit_esp,*) esp(n,k,s),occ(n,k,s),res(n,k,s)
!          end do
!       end do
!       end do
!       close(unit_esp)

!       open( unit_kpt, FILE=file_kpt_list_dos, POSITION="append" )
!       do i=1,MMBZ
!          k=whole_2_ir(i)
!          if ( any( nint(res(1:max_ref_band,k,:)) /= -1 ) ) cycle
!          write(unit_kpt,*) whole_kpt_grid(1:3,i)
!       end do ! i
!       close(unit_kpt)

    end if

    call write_border( 0, " calc_dos(end)" )

  END SUBROUTINE calc_dos


END MODULE dos_module
