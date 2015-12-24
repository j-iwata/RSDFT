MODULE polint_module

  implicit none

  PRIVATE
  PUBLIC :: polint, dpolint

CONTAINS


  SUBROUTINE polint( xa,ya,n,x,y,dy )

    implicit none
    integer,intent(IN)  :: n
    real(8),intent(IN)  :: xa(n),ya(n),x
    real(8),intent(OUT) :: y,dy
    integer :: i,m
    real(8) :: a,b,f0nm,f0nm_0
    real(8),allocatable :: f0(:)

    allocate( f0(n) )
    f0(:) = ya(:)

    f0nm = ya(1)

    do m=1,n-1
       do i=1,n-m
          a = x-xa(i)
          b = x-xa(i+m)
          f0(i) = ( b*f0(i) - a*f0(i+1) )/( b - a )
       end do
       f0nm_0 = f0nm
       f0nm   = sum( f0(1:n-m) )/(n-m)
    end do
    y  = f0(1)
    dy = f0nm - f0nm_0

    deallocate( f0 )

    return
  END SUBROUTINE polint


  SUBROUTINE dpolint( xa,ya,n,x,y,dy )

    implicit none
    integer,intent(IN)  :: n
    real(8),intent(IN)  :: xa(n),ya(n),x
    real(8),intent(OUT) :: y,dy
    integer :: i,m
    real(8) :: a,b,f1nm,f1nm_0
    real(8),allocatable :: f0(:),f1(:)

    allocate( f0(n),f1(n) )
    f0(:) = ya(:)
    f1(:) = 0.0d0

    f1nm = ( ya(2)-ya(1) )/( xa(2)-xa(1) )

    do m=1,n-1
       do i=1,n-m
          a = x-xa(i)
          b = x-xa(i+m)
          f1(i) = ( f0(i) - f0(i+1) + b*f1(i) -a*f1(i+1) )/( b - a )
          f0(i) = ( b*f0(i) - a*f0(i+1) )/( b - a )
       end do
       f1nm_0 = f1nm
       f1nm   = sum( f1(1:n-m) )/(n-m)
    end do

    y  = f1(1)
    dy = f1nm - f1nm_0

    deallocate( f0,f1 )

    return
  END SUBROUTINE dpolint


END MODULE polint_module
