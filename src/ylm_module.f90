MODULE ylm_module ! Spherical Harmonic Function

  implicit none

  PRIVATE
  PUBLIC :: Ylm

CONTAINS

  FUNCTION Ylm( x,y,z, L, M_in )

    implicit none
    real(8) :: Ylm
    real(8),intent(IN) :: x,y,z
    integer,intent(IN) :: L,M_in
    integer :: k1,k2,k,M
    real(8) :: r,r2,Clm,ct,phi,c,pi4
    real(8),parameter :: ep=1.d-10, zero=0.0d0
    real(8),parameter :: half=0.5d0, one=1.0d0, two=2.0d0
    real(8),parameter :: pi=3.141592653589793238d0

    Ylm = zero
    M   = abs(M_in)
    pi4 = 4.0d0*pi

    if ( L < 0 .or. M > L ) then
       write(*,*) "L,M_in=",L,M_in
       stop "Bad arguments (stop at Ylm)"
    end if

    if ( L == 0 ) then
       Ylm=one/sqrt(pi4)
       return
    end if

    if ( x == 0.0d0 .and. y == 0.0d0 .and. z == 0.0d0 ) return

    c=one
    do k=L-M+1,L+M
       c=c*dble(k)
    end do
    Clm=sqrt((two*L+one)/(c*pi4))

    r  = sqrt(x*x+y*y+z*z)
    ct = z/r
    if ( abs(x) < ep ) then
       phi = half*pi*sign(one,y)
    else
       phi = atan(y/x)
       if ( x < 0 ) phi = phi + pi
    end if

    if ( M_in > 0 ) then

       Ylm = sqrt(two)*Clm*Plm(L,M,ct)*cos(phi*M)

    else if ( M_in == 0 ) then

       Ylm = Clm*Plm(L,M,ct)

    else if ( mod(M,2) == 0 ) then

       Ylm =-sqrt(two)*Clm*Plm(L,M,ct)*sin(phi*M)

    else

       Ylm = sqrt(two)*Clm*Plm(L,M,ct)*sin(phi*M)

    end if

    return
  END FUNCTION Ylm

!----------------------------------- Associated Legendre Polynomial Plm(x)

  FUNCTION Plm(L,M,xi)

    implicit none
    real(8) :: Plm
    integer,intent(IN) :: L,M
    real(8),intent(IN) :: xi
    real(8),parameter :: ep=1.d-13,one=1.0d0,two=2.0d0
    real(8) :: pmm,somx2,fact,pmmp1,pll,tmp,x
    integer :: i,ll

    x = xi
    tmp = abs(abs(x)-one) ; if ( tmp < ep ) x = sign(one,x)

    if( m < 0 .or. m > l .or. abs(x) > one ) then
       stop "Bad arguments (stop at Plm)"
    end if

    pmm=one

    if ( m > 0 ) then
       somx2=sqrt((one-x)*(one+x))
       fact=one
       do i=1,M
          pmm=-pmm*fact*somx2
          fact=fact+two
       end do
    end if

    if ( L == M ) then
       Plm = pmm
    else
       pmmp1 = x*(two*m+one)*pmm
       if ( L == M+1 ) then
          Plm=pmmp1
       else
          do ll=M+2,L
             pll=(x*(two*ll-one)*pmmp1-(ll+M-one)*pmm)/(ll-M)
             pmm=pmmp1
             pmmp1=pll
          end do
          Plm=pll
       end if
    end if

    return
  END FUNCTION Plm


END MODULE ylm_module
