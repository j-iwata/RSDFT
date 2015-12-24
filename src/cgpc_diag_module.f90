MODULE cgpc_diag_module

!$  use omp_lib
  use kinetic_variables, only: coef_lap0, const_k2
  use localpot_module, only: Vloc
  use ps_nloc2_variables, only: nzlma, Mlma, uVk, JJP, iuV, MJJ

  implicit none

  PRIVATE
  PUBLIC :: init_cgpc_diag, cgpc_diag

  real(8),allocatable :: diag_H(:)

CONTAINS


  SUBROUTINE init_cgpc_diag(n1,n2,k,s,dV)
    implicit none
    integer,intent(IN) :: n1,n2,k,s
    real(8),intent(IN) :: dV
    integer :: i,j,lma
    real(8) :: c

    if ( .not.allocated(diag_H) ) allocate( diag_H(n1:n2) )
    diag_H(:)=0.0d0

    do i=n1,n2
       diag_H(i) = coef_lap0 + const_k2(k) + Vloc(i,s)
    end do

    i=3
    select case(i)
    case( 2 )
       do lma=1,nzlma
          c=iuV(lma)*dV
          do j=1,MJJ(lma)
             i=JJP(j,lma)
             diag_H(i) = diag_H(i) + c*abs(uVk(j,lma,k))**2
          end do
       end do
    case( 3 )
       do lma=1,Mlma
          c=iuV(lma)*dV
          do i=n1,n2
             diag_H(i) = diag_H(i) + c*abs(uVk(i,lma,k))**2
          end do
       end do
    end select

    do i=n1,n2
       if ( diag_H(i) == 0.0d0 ) stop "stop@init_cgpc_diag"
       diag_H(i) = 1.0d0/diag_H(i)
    end do

  END SUBROUTINE init_cgpc_diag


  SUBROUTINE cgpc_diag(m,n,Pgk)
    implicit none
    integer,intent(IN)       :: m,n
#ifdef _DRSDFT_
    real(8),intent(INOUT) :: Pgk(m,n)
#else
    complex(8),intent(INOUT) :: Pgk(m,n)
#endif
    integer :: j

    do j=1,n
       Pgk(:,j) = Pgk(:,j)*diag_H(:)
    end do

  END SUBROUTINE cgpc_diag


END MODULE cgpc_diag_module
