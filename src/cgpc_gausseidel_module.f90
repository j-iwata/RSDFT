MODULE cgpc_gausseidel_module

!$  use omp_lib
  use bc_module, only: www, bcset, bcset_1
  use parallel_module
  use omp_variables

  implicit none

  PRIVATE
  PUBLIC :: init_cgpc_gausseidel, cgpc_gausseidel

  integer :: a1b,b1b,a2b,b2b,a3b,b3b
  integer :: ML1,ML2,ML3,ML_0
  integer :: mloop = 6
  real(8) :: c0,c1,c2,c3
  logical :: iswitch_bc = .true.
  real(8),parameter :: tol = 1.d-13

CONTAINS


  SUBROUTINE init_cgpc_gausseidel( Igrid, Ngrid, Hgrid, ggg, mloop_in )
    implicit none
    integer,intent(IN) :: Igrid(2,0:3), Ngrid(0:3), mloop_in
    real(8),intent(IN) :: Hgrid(3), ggg(6)
    ML_0= Igrid(1,0)
    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)
    c1  = 1.0d0/Hgrid(1)**2 * ggg(1)
    c2  = 1.0d0/Hgrid(2)**2 * ggg(2)
    c3  = 1.0d0/Hgrid(3)**2 * ggg(3)
    c0  = 0.5d0/( c1 + c2 + c3 )
    if ( mloop_in > 0 ) mloop = mloop_in
  END SUBROUTINE init_cgpc_gausseidel


  SUBROUTINE cgpc_gausseidel(m,n,Pgk)
    implicit none
    integer,intent(IN) :: m,n
#ifdef _DRSDFT_
    real(8),intent(INOUT) :: Pgk(m,n)
    real(8),parameter :: zero=0.0d0
    real(8),allocatable :: work(:,:,:)
#else
    complex(8),intent(INOUT) :: Pgk(m,n)
    complex(8),parameter :: zero=0.0d0
    complex(8),allocatable :: work(:,:,:)
#endif
    integer :: i,i1,i2,i3,j,loop,j1,j2,j3,ierr
    integer :: m0,m1,n1_omp,n2_omp,a1b_omp,a2b_omp,a3b_omp
    integer :: b1b_omp,b2b_omp,b3b_omp
    real(8),allocatable :: err0(:)
    real(8) :: err1,err

!$OMP parallel private( m0,m1,n1_omp,n2_omp,a1b_omp,a2b_omp,a3b_omp &
!$OMP                  ,b1b_omp,b2b_omp,b3b_omp,i )

    m0 = 0
    m1 = 1
!$  m0 = omp_get_thread_num()
!$  m1 = omp_get_num_threads()

    n1_omp  = Igrid_omp(1,0,m0) - ML_0 + 1
    n2_omp  = Igrid_omp(2,0,m0) - ML_0 + 1
    a1b_omp = Igrid_omp(1,1,m0)
    b1b_omp = Igrid_omp(2,1,m0)
    a2b_omp = Igrid_omp(1,2,m0)
    b2b_omp = Igrid_omp(2,2,m0)
    a3b_omp = Igrid_omp(1,3,m0)
    b3b_omp = Igrid_omp(2,3,m0)

!$omp single
    allocate( err0(0:m1-1) ) ; err0=0.0d0
    allocate( work(a1b:b1b,a2b:b2b,a3b:b3b) )
!$omp end single

    do i3=a3b_omp,b3b_omp
    do i2=a2b_omp,b2b_omp
    do i1=a1b_omp,b1b_omp
       work(i1,i2,i3)=zero
    end do
    end do
    end do

!$omp workshare
    www=zero
!$omp end workshare

    do j=1,n
       do i3=a3b_omp,b3b_omp
       do i2=a2b_omp,b2b_omp
       do i1=a1b_omp,b1b_omp
          i=1+(i1-a1b)+(b1b-a1b+1)*((i2-a2b)+(i3-a3b)*(b2b-a2b+1))
          www(i1,i2,i3,j) = Pgk(i,j)
       end do
       end do
       end do
    end do

!$omp barrier
    if ( iswitch_bc ) call bcset_1(1,n,1,0)
!$omp barrier

    do j=1,n
    do loop=1,mloop
       err0(m0)=0.0d0
       do i3=a3b_omp,b3b_omp
       do i2=a2b_omp,b2b_omp
       do i1=a1b_omp,b1b_omp
          i=1+(i1-a1b)+(b1b-a1b+1)*((i2-a2b)+(i3-a3b)*(b2b-a2b+1))
          www(i1,i2,i3,j) = c0*( 2.0d0*Pgk(i,j) &
               + c1*( www(i1-1,i2,i3,j) + www(i1+1,i2,i3,j) ) &
               + c2*( www(i1,i2-1,i3,j) + www(i1,i2+1,i3,j) ) &
               + c3*( www(i1,i2,i3-1,j) + www(i1,i2,i3+1,j) ) )
          err0(m0) = err0(m0) + abs( www(i1,i2,i3,j) - work(i1,i2,i3) )**2
          work(i1,i2,i3) = www(i1,i2,i3,j)
       end do
       end do
       end do
!$omp barrier
!$omp master
       err1=sum(err0)
       call mpi_allreduce(err1,err,1,MPI_REAL8,MPI_SUM,comm_grid,ierr)
!$omp end master
!$omp barrier
       if ( err < tol ) exit
       if ( iswitch_bc ) call bcset_1(j,j,1,0)
!$omp barrier
    end do ! loop
    end do ! j

    do j=1,n
       do i3=a3b_omp,b3b_omp
       do i2=a2b_omp,b2b_omp
       do i1=a1b_omp,b1b_omp
          i=1+(i1-a1b)+(b1b-a1b+1)*((i2-a2b)+(i3-a3b)*(b2b-a2b+1))
          Pgk(i,j) = www(i1,i2,i3,j)
       end do
       end do
       end do
    end do
!$OMP end parallel

    deallocate( work )

  END SUBROUTINE cgpc_gausseidel


  SUBROUTINE cgpc_gausseidel_0(m,n,Pgk)
    implicit none
    integer,intent(IN) :: m,n
#ifdef _DRSDFT_
    real(8),intent(INOUT) :: Pgk(m,n)
    real(8),parameter :: zero=0.0d0
    real(8),allocatable :: work(:,:,:)
#else
    complex(8),intent(INOUT) :: Pgk(m,n)
    complex(8),parameter :: zero=0.0d0
    complex(8),allocatable :: work(:,:,:)
#endif
    integer :: i,i1,i2,i3,j,loop,j1,j2,j3,ierr
    real(8) :: err0,err

    allocate( work(a1b:b1b,a2b:b2b,a3b:b3b) ) ; work=zero

    www=zero
    do j=1,n
       do i3=a3b,b3b
       do i2=a2b,b2b
       do i1=a1b,b1b
          i=1+(i1-a1b)+(b1b-a1b+1)*((i2-a2b)+(i3-a3b)*(b2b-a2b+1))
          www(i1,i2,i3,j) = Pgk(i,j)
       end do
       end do
       end do
    end do

    if ( iswitch_bc ) call bcset(1,n,1,0)

    do j=1,n
    do loop=1,mloop
       err0=0.0d0
       do i3=a3b,b3b
       do i2=a2b,b2b
       do i1=a1b,b1b
          i=1+(i1-a1b)+(b1b-a1b+1)*((i2-a2b)+(i3-a3b)*(b2b-a2b+1))
          www(i1,i2,i3,j) = c0*( 2.0d0*Pgk(i,j) &
               + c1*( www(i1-1,i2,i3,j) + www(i1+1,i2,i3,j) ) &
               + c2*( www(i1,i2-1,i3,j) + www(i1,i2+1,i3,j) ) &
               + c3*( www(i1,i2,i3-1,j) + www(i1,i2,i3+1,j) ) )
          err0 = err0 + abs( www(i1,i2,i3,j) - work(i1,i2,i3) )**2
          work(i1,i2,i3) = www(i1,i2,i3,j)
       end do
       end do
       end do
       call mpi_allreduce(err0,err,1,MPI_REAL8,MPI_SUM,comm_grid,ierr)
       if ( err < tol ) exit
       if ( iswitch_bc ) call bcset(j,j,1,0)
    end do ! loop
    end do ! j

    do j=1,n
       do i3=a3b,b3b
       do i2=a2b,b2b
       do i1=a1b,b1b
          i=1+(i1-a1b)+(b1b-a1b+1)*((i2-a2b)+(i3-a3b)*(b2b-a2b+1))
          Pgk(i,j) = www(i1,i2,i3,j)
       end do
       end do
       end do
    end do

    deallocate( work )

  END SUBROUTINE cgpc_gausseidel_0


END MODULE cgpc_gausseidel_module
