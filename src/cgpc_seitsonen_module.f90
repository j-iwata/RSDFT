! --- Reference -------------------------------------
! Physical Review B 51, 14057 (1995)
! Ari P. Seitsonen, M. J. Puska, and R. M. Nieminen
!----------------------------------------------------
MODULE cgpc_seitsonen_module

  use localpot_module, only: Vloc
  use kinetic_module, only: op_kinetic

  implicit none

  PRIVATE
  PUBLIC :: init_cgpc_seitsonen, cgpc_seitsonen

  integer :: ML_0, ML_1
  integer :: comm_grid
  real(8) :: dV
  real(8) :: A = 1.d0

CONTAINS


  SUBROUTINE init_cgpc_seitsonen(n1,n2,comm_in,dV_in)
    implicit none
    integer,intent(IN) :: n1,n2,comm_in
    real(8),intent(IN) :: dV_in
    dV = dV_in
    comm_grid = comm_in
    ML_0 = n1
    ML_1 = n2
  END SUBROUTINE init_cgpc_seitsonen


  SUBROUTINE cgpc_seitsonen(k,s,mm,nn,E,wf,Pgk)
    implicit none
    integer,intent(IN) :: k,s,mm,nn
    real(8),intent(IN) :: E(nn)
#ifdef _DRSDFT_
    real(8),intent(IN)    :: wf(mm,nn)
    real(8),intent(INOUT) :: Pgk(mm,nn)
    real(8),allocatable :: Hwf(:,:)
    real(8),parameter :: zero=0.0d0
#else
    complex(8),intent(IN)    :: wf(mm,nn)
    complex(8),intent(INOUT) :: Pgk(mm,nn)
    complex(8),allocatable :: Hwf(:,:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
#endif
    integer :: i,n,ierr
    real(8) :: Kr,x
    real(8),allocatable :: T0(:),T1(:)
    include 'mpif.h'
    real(8),parameter :: c1=2.7d1, c2=1.8d1, c3=1.2d1, c4=8.0d0, c5=1.6d1

    allocate( T0(nn) ) ; T0=0.0d0
    allocate( T1(nn) ) ; T1=0.0d0
    allocate( Hwf(mm,nn) ) ; Hwf=zero

    call op_kinetic(k,wf,Hwf,1,mm,1,nn)

    do n=1,nn
#ifdef _DRSDFT_
       T0(n) = sum( wf(:,n)*Hwf(:,n) )*dV
#else
       T0(n) = sum( conjg(wf(:,n))*Hwf(:,n) )*dV
#endif
    end do

    deallocate( Hwf )

    call MPI_ALLREDUCE(T0,T1,nn,MPI_REAL8,MPI_SUM,comm_grid,ierr)

    T1(:) = 1.0d0/T1(:)

    do n=1,nn
       do i=1,mm
          x  = A*abs( E(n) - Vloc(i-1+ML_0,s) )*T1(n)
          Kr = ( c1 + c2*x + c3*x**2 + c4*x**3 ) &
              /( c1 + c2*x + c3*x**2 + c4*x**3 + c5*x**4 )
          Pgk(i,n) = Kr*Pgk(i,n)
       end do
    end do

    deallocate( T1 )
    deallocate( T0 )

  END SUBROUTINE cgpc_seitsonen


END MODULE cgpc_seitsonen_module
