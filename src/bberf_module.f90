MODULE bberf_module ! error function

  PRIVATE
  PUBLIC :: bberf
  PUBLIC :: bberfc

CONTAINS

!------------------------------------------------------
! written by J. Iwata
!------------------------------------------------------

  FUNCTION bberf( x )
    implicit none
    real(8) :: bberf
    real(8),intent(IN) :: x
#ifdef _NO_ERF_
    bberf=1.0d0-ccerf(x) ! if buil-in erf is not avaialbe
#else
    bberf=erf(x)
#endif
  END FUNCTION bberf


  FUNCTION bberfc( x )
    implicit none
    real(8) :: bberfc
    real(8),intent(IN) :: x
#ifdef _NO_ERF_
    bberfc=ccerf(x) ! if buil-in erfc is not avaialbe
#else
    bberfc=erfc(x)
#endif
  END FUNCTION bberfc


  FUNCTION ccerf( x )
    implicit none
    real(8) :: ccerf
    real(8),intent(IN) :: x
    real(8),parameter :: pi = 3.141592653589793d0
    real(8) :: y,alpha,tny,eps,delta
    real(8) :: ai,bi,f0,f1,C0,C1,D0,D1,a0
    integer :: nmax,i

    if ( abs(x) < 0.1d0 ) then
       ccerf=1.0d0-2.0d0/sqrt(pi)*(x-x**3/3.0d0+0.1d0*x**5 &
                  -x**7/42.0d0+x**9/216.0d0-x**11/1320.0d0 )
       return
    end if

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

    ccerf = sign( f1*exp(-y)*sqrt(y/pi), x ) + 1.0d0-sign(1.0d0,x)

  END FUNCTION ccerf


END MODULE bberf_module
