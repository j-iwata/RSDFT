MODULE expint_module

  implicit none

  PRIVATE
  PUBLIC :: expint

CONTAINS


  FUNCTION expint( n, x )
    implicit none
    real(8) :: expint
    integer,intent(IN) :: n
    real(8),intent(IN) :: x
    expint=expint_2( x )
!    expint=expint_1( n, x )
  END FUNCTION expint


  FUNCTION expint_2( x )

    implicit none
    real(8) :: expint_2
    real(8),intent(IN) :: x
    real(8),parameter :: euler=0.577215664901532860606512d0
    real(8) :: y,tny,eps,delta
    real(8) :: ai,bi,f0,f1,C0,C1,D0,D1,a0
    integer :: nmax,i,j

    nmax = 100000
    eps  = epsilon(1.0d0)
    tny  = tiny(1.0d0)

    if ( x <= 0.0d0 ) then

       stop "bad argument (stop@expint)"

    else if ( x >= 0.5d0 ) then

       f0 = tny
       C0 = f0
       D0 = 0.0d0
       a0 = 1.0d0

       do i=1,nmax
          bi = x + 2.0d0*i - 1.0d0
          ai = -(i-1)*(i-1) + a0
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

       if ( i > nmax ) stop "continued fraction is not converged(stop@expint)"

       expint_2 = f1*exp(-x)

    else if ( x < 0.5d0 ) then

       f0 = 0.0d0

       do i=1,nmax

          f1 = f0

          ai = 1.0d0
          do j=1,i
             ai = ai*( -x/dble(j) )
          end do

          f0 = f0 + ai/dble(i)

          if ( abs(f0-f1) < eps ) exit

       end do

       if ( i > nmax ) stop "power series is not converged(stop@expint)"

       expint_2 = - euler - log(x) - f0

    end if

  END FUNCTION expint_2


  function expint_1(n,x) result(expint)
    implicit none
    integer,intent(in) :: n
    real(8),intent(in) :: x
    real(8) :: expint
    integer,parameter :: maxit=200
    real(8),parameter :: esp=1.d-12, big=huge(x)*esp
    real(8),parameter :: euler=0.577215664901532860606512d0
    integer :: i,nm1,j
    real(8) :: a,b,c,d,del,fact,h,arsum
 
    if ( .not.(n>=0.and.x>=0.d0.and.(x>0.d0.or.n>1)) ) then
       write(*,*) 'Bad arguments in expint.f'
       stop
    end if
 
    if ( n==0 ) then
       expint=exp(-x)/x
       return
    end if
    nm1=n-1
    if ( x==0.d0 ) then
       expint=1.d0/nm1
    else if ( x>1.d0 ) then
       b=x+n
       c=big
       d=1.d0/b
       h=d
       do i=1,maxit
          a=-i*(nm1+i)
          b=b+2.d0
          d=1.d0/(a*d+b)
          c=b+a/c
          del=c*d
          h=h*del
          if ( abs(del-1.d0)<=esp ) exit
       end do
       if ( i>maxit ) then
          write(*,*) 'Continued fraction failed in expint.f'
          stop
       end if
       expint=h*exp(-x)
    else
       if ( nm1/=0 ) then
          expint=1.d0/nm1
       else
          expint=-log(x)-euler
       end if
       fact=1.d0
       do i=1,maxit
          fact=-fact*x/i
          if ( i/=nm1 ) then
             del=-fact/(i-nm1)
          else
             arsum = 0.d0
             do j=1,nm1
                arsum = arsum + 1.d0/j
             end do
             del = fact*(-LOG(x)-euler+arsum)
          end if
          expint=expint+del
          if ( abs(del)<abs(expint)*esp ) exit
       end do
       if ( i>maxit ) then
          write(*,*) 'series failed in expint.f'
          stop
       end if
    end if

  end function expint_1

END MODULE expint_module
