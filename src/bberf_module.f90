MODULE bberf_module ! error function

  PRIVATE
  PUBLIC :: bberf

CONTAINS

!------------------------------------------------------
! written by J. Iwata
!------------------------------------------------------

  FUNCTION bberf( x )
    implicit none
    real(8) :: bberf
    real(8),intent(IN) :: x
    bberf = erf(x)
!   bberf=1.0d0-ccerf(x) ! if buil-in erf is not avaialbe
  END FUNCTION bberf


  FUNCTION ccerf( x )
    implicit none
    real(8) :: ccerf
    real(8),intent(IN) :: x
    real(8),parameter :: pi = 3.141592653589793d0
    real(8) :: y,alpha,tny,eps,delta
    real(8) :: ai,bi,f0,f1,C0,C1,D0,D1,a0
    integer :: nmax,i

    nmax  = 100000
    eps   = epsilon(1.0d0)
    tny   = tiny(1.0d0)
    alpha = 0.5d0

    y = x*x

    f0 = tny
    C0 = f0
    D0 = 0.0d0
    a0 = 1.0d0

    do i=1,nmax
       bi = y + 2.0d0*i - 1.0d0 - alpha
       ai = -(i-1)*(i-1-alpha) + a0
       ai = -(i-1)*(i-1.5d0) + a0
       D1 = bi + ai*D0
       if ( D1 == 0.0d0 ) D1=tny
       C1 = bi + ai/C0
       if ( C1 == 0.0d0 ) C1=tny
       D1 = 1.0d0/D1
       delta = C1*D1
       f1 = f0*delta
       if ( abs(delta-1.0d0) < eps ) exit
       f0 = f1
       C0 = C1
       D0 = D1
       a0 = 0.0d0
    end do

    ccerf = f1*exp(-y)*sqrt(y/pi)

  END FUNCTION ccerf


END MODULE bberf_module
