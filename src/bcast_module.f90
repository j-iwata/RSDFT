MODULE bcast_module

  implicit none

  PRIVATE
  PUBLIC :: test_bcast,d_rsdft_bcast,z_rsdft_bcast

  include 'mpif.h'

  integer :: n_opt = 1024
  integer :: n_opt_h = 1024

!  INTERFACE rsdft_bcast
!     MODULE PROCEDURE d_rsdft_bcast,z_rsdft_bcast
!  END INTERFACE

CONTAINS

  SUBROUTINE d_rsdft_bcast(a,n,type,root,comm,ierr)
    implicit none
    integer,intent(IN) :: n,type,root,comm
    real(8),intent(IN) :: a(n)
    integer,intent(OUT) :: ierr
    integer :: m,i
    call MPI_BCAST(a,n,type,root,comm,ierr) ; return
    do i=1,n,n_opt
       m=min(n-i+1,n_opt)
       call MPI_BCAST(a(i),m,type,root,comm,ierr)
    end do
  END SUBROUTINE d_rsdft_bcast

  SUBROUTINE z_rsdft_bcast(a,n,type,root,comm,ierr)
    implicit none
    integer,intent(IN) :: n,type,root,comm
    complex(8),intent(IN) :: a(n)
    integer,intent(OUT) :: ierr
    integer :: m,i
    call MPI_BCAST(a,n,type,root,comm,ierr) ; return
    do i=1,n,n_opt_h
       m=min(n-i+1,n_opt_h)
       call MPI_BCAST(a(i),m,type,root,comm,ierr)
    end do
  END SUBROUTINE z_rsdft_bcast

  SUBROUTINE test_bcast
    implicit none
    integer :: i,j,k,m,mt,n,ierr,myrank
    real(8),allocatable :: a(:)
    real(8) :: ct,ct0,ct1,ctmin,et,et0,et1,etmin
 
    return

    call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
    n=2**20
    allocate( a(n) ) ; a(:)=1.d0
    ctmin=1.d100
    etmin=1.d100
    do j=0,20
       m = 2**j
       mt=0
       do i=1,2**(20-j)
          mt=mt+m
       end do
       call cpu_time(ct0) ; et0=mpi_wtime()
       do k=1,10
          do i=1,2**(20-j)
             call MPI_BCAST(a,m,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          end do
       end do ! k
       call cpu_time(ct1) ; et1=mpi_wtime()
       ct=ct1-ct0
       et=et1-et0
       ctmin=min(ct,ctmin)
       if ( et < etmin ) then
          n_opt = m
          etmin = et
       end if
       if ( myrank == 0 ) then
          write(*,'(1x,3i12,4f10.5)') m,i-1,mt,ct,ctmin,et,etmin
       end if
    end do ! j
    deallocate( a )
    call MPI_BCAST(n_opt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n_opt_h=n_opt/2
    if ( myrank == 0 ) write(*,*) "n_opt, n_opt_h=",n_opt,n_opt_h
  END SUBROUTINE test_bcast

END MODULE bcast_module
