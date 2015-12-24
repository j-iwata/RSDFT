MODULE hartree_sol_module

  use hartree_variables, only: E_hartree,Vh
  use rgrid_module, only: Igrid,dV
  use ggrid_module, only: GG,LLG,MGL,NGgrid,construct_ggrid,destruct_ggrid
  use parallel_module, only: comm_grid, mpi_real8, mpi_sum
  use watch_module
  use hartree_sol_ffte_module, only: calc_hartree_sol_ffte
  use hartree_sol_fftw_module, only: calc_hartree_sol_fftw
  use fft_module

  implicit none

  PRIVATE
  PUBLIC :: calc_hartree_sol

CONTAINS

  SUBROUTINE calc_hartree_sol( n1,n2,n3, rho )

    implicit none
    integer,intent(IN) :: n1,n2,n3
    real(8),intent(IN) :: rho(n1:n2,n3)
    integer :: i,i1,i2,i3
    real(8) :: Eh0,pi4,g2,ctt(0:5),ett(0:5)
    real(8),allocatable :: work(:)
    complex(8),parameter :: z0=(0.d0,0.d0)
    complex(8),allocatable :: zwork0(:,:,:),zwork1(:,:,:)
    logical :: disp_sw

#ifdef _FFTE_
    call calc_hartree_sol_ffte(n1,n2,n3,rho)
    return
#elif _FFTW_
    call calc_hartree_sol_fftw(n1,n2,n3,rho)
    return
#endif

    call write_border( 1, " calc_hartree_sol(start)" )

    ctt(:)=0.d0
    ett(:)=0.d0

    call watch(ctt(0),ett(0))

    allocate( work(n1:n2) ) ; work=0.0d0

    do i=1,n3
       work(n1:n2) = work(n1:n2) + rho(n1:n2,i)
    end do

    call init_fft
    call d1_to_z3_fft( work, zwork0 )

    call watch(ctt(1),ett(1))

    call forward_fft( zwork0, zwork1 )

    call watch(ctt(2),ett(2))

! ---

    call construct_Ggrid(2)

    pi4 = 4.0d0*acos(-1.d0)

    zwork1(:,:,:)=z0
    do i=1,NGgrid(0)
       g2=GG(MGL(i))
       if ( g2 == 0.0d0 ) cycle
       i1=LLG(1,i)
       i2=LLG(2,i)
       i3=LLG(3,i)
       zwork1(i1,i2,i3)=zwork0(i1,i2,i3)*pi4/g2
    end do

    call destruct_Ggrid

! ---

    call watch(ctt(3),ett(3))

    call backward_fft( zwork1, zwork0 )

    call watch(ctt(4),ett(4))

    call finalize_fft

    if ( allocated(zwork0) ) deallocate( zwork0 )

    call z3_to_d1_fft( zwork1, Vh )

    Eh0=0.0d0
    i=n1-1
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       i=i+1
       Eh0 = Eh0 + dble( zwork1(i1,i2,i3) )*work(i)
    end do
    end do
    end do
    Eh0=0.5d0*Eh0*dV
    call mpi_allreduce(Eh0,E_hartree,1,mpi_real8,mpi_sum,comm_grid,i)

    if ( allocated(zwork1) ) deallocate( zwork1 )
    deallocate( work )

    call watch(ctt(5),ett(5))

!    call check_disp_switch( disp_sw, 0 )
!    if ( disp_sw ) then
!       write(*,*) "time(hatree1)=",ctt(1)-ctt(0),ett(1)-ett(0)
!       write(*,*) "time(hatree2)=",ctt(2)-ctt(1),ett(2)-ett(1)
!       write(*,*) "time(hatree3)=",ctt(3)-ctt(2),ett(3)-ett(2)
!       write(*,*) "time(hatree4)=",ctt(4)-ctt(3),ett(4)-ett(3)
!       write(*,*) "time(hatree5)=",ctt(5)-ctt(4),ett(5)-ett(4)
!    end if

    call write_border( 1, " calc_hartree_sol(end)" )

  END SUBROUTINE calc_hartree_sol

END MODULE hartree_sol_module
