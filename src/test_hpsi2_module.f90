MODULE test_hpsi2_module

  use array_bound_module, only: MB
  use localpot_module
  use hamiltonian_module
  use parallel_module
  use watch_module

  implicit none

  PRIVATE
  PUBLIC :: test_hpsi2

CONTAINS

  SUBROUTINE test_hpsi2(nloop)
    implicit none
    integer,intent(IN) :: nloop
#ifdef _DRSDFT_
    real(8),parameter :: zero=0.0d0, one=1.0d0
    real(8), allocatable :: tpsi(:,:)
    real(8), allocatable :: htpsi(:,:)
#else
    real(8),parameter :: zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0)
    complex(8), allocatable :: tpsi(:,:)
    complex(8), allocatable :: htpsi(:,:)
#endif
    integer :: i, j, n1, n2, nrhs, niter, ierr
    integer :: loop,ii
    real(8), allocatable :: sums(:)
    real(8) :: t0, t1, t2, t3, t4
    logical,save :: flag_allocate=.false.

    n1  = idisp(myrank)+1
    n2  = idisp(myrank)+ircnt(myrank)

    if ( .not.allocated(Vloc) ) then
       allocate( Vloc(n1:n2,1) ) ; Vloc=0.0d0
       flag_allocate=.true.
    end if

    do ii = 0,10

    time_hmlt(:,:)=0.0d0
    time_kine(:,:)=0.0d0
    time_nlpp(:,:)=0.0d0

    nrhs = 2**ii
    if ( nrhs > MB_d ) exit

    niter = MB/nrhs

    allocate( tpsi(n1:n2,nrhs)  ) ; tpsi=one
    allocate( htpsi(n1:n2,nrhs) ) ; htpsi=zero
    allocate( sums(nrhs)        ) ; sums=0.0d0

!    call myrand(tpsi, n2-n1+1, nrhs)
  
    if ( DISP_SWITCH_PARALLEL ) then
       write(*,*) '------------ Start test_hpsi2 ------------'
    end if

    t0 = mpi_wtime()

    do loop=1,nloop
       do i = 1, niter
          call hamiltonian(1, 1, tpsi, htpsi, n1, n2, 1, nrhs) 
       end do
    end do

    t1 = mpi_wtime()

!$OMP parallel private( loop, i )
    do loop=1,nloop
       do i=1,niter
          call op_kinetic(1,tpsi,htpsi,n1,n2,1,nrhs)
       end do
    end do
!$OMP end parallel

    t2 = mpi_wtime()

!$OMP parallel private( loop, i )
    do loop=1,nloop
       do i=1,niter
          call op_localpot(1,n2-n1+1,nrhs,tpsi,htpsi)
       end do
    end do
!$OMP end parallel

    t3 = mpi_wtime()

!$OMP parallel private( loop, i )
    do loop=1,nloop
       do i=1,niter
          call op_nonlocal(1,1,tpsi,htpsi,n1,n2,1,nrhs)
       end do
    end do
!$OMP end parallel

    t4 = mpi_wtime()

    sums(:) = 0d0
    do i = 1, nrhs
       do j = n1, n2
#ifdef _DRSDFT_
          sums(i) = sums(i) + htpsi(j, i)*htpsi(j, i)
#else
          sums(i) = sums(i) + conjg(htpsi(j, i))*htpsi(j, i)
#endif
       end do
    end do
    call mpi_allreduce(MPI_IN_PLACE, sums, nrhs, &
         MPI_DOUBLE_PRECISION, MPI_SUM, comm_grid, ierr)
  
    if ( DISP_SWITCH_PARALLEL ) then
       write(*,*) 'nloop =',nloop
       write(*,*) 'nrhs = ', nrhs
       write(*,*) 'niter = ', niter
       write(*,*) 'time(tot) = ', t1 - t0
       write(*,*) 'time(kin) = ', t2 - t1
       call write_watchb( time_kine, 11, time_kine_indx )
       write(*,*) 'time(loc) = ', t3 - t2
       write(*,*) 'time(nlc) = ', t4 - t3
       call write_watchb( time_nlpp, 3, time_nlpp_indx )
       write(*,*) 'check sum'
       do i = 1, nrhs
          write(*,*) sums(i),sum(sums)
       end do
    end if

    deallocate( sums )
    deallocate( htpsi )
    deallocate( tpsi )

    end do ! ii

    if ( flag_allocate ) then
       deallocate( Vloc )
       flag_allocate=.false.
    end if

  END SUBROUTINE test_hpsi2

  subroutine myrand(V, nrow, ncol)
    implicit none
#ifdef _DRSDFT_
    real(8), intent(out) :: V(nrow,*)
#else
    complex(8), intent(out) :: V(nrow,*)  
#endif
    real(8), allocatable :: tmpV(:,:)
    integer, intent(in) :: nrow, ncol
    integer :: iseed(4)

    allocate(tmpV(nrow,ncol))
    iseed(1) = modulo(myrank_g-2*4096, 4096) ! must be between 0 and 4095
    iseed(2) = modulo(myrank_g-4096, 4096) ! must be between 0 and 4095
    iseed(3) = modulo(myrank_g, 4096) ! must be between 0 and 4095
    iseed(4) = 1 ! must be between 0 and 4095 and odd
  
    call DLARNV(2, iseed, nrow*ncol, tmpV)
    V(:,1:ncol) = tmpV(:,1:ncol)
  end subroutine myrand

END MODULE test_hpsi2_module


