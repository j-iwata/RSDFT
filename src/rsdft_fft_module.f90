MODULE rsdft_fft_module

  implicit none

  PRIVATE
  PUBLIC :: rsdft_fft3d

  integer :: mm
  real(8),parameter :: pi =3.14159265358979323846264338327950d0
 !real(8),parameter :: pi2=6.28318530717958647692528676655901d0
  real(8) :: pi2
  complex(8),parameter :: z0=(0.0d0,0.0d0)
  complex(8),parameter :: zi=(0.0d0,1.0d0)

CONTAINS


  SUBROUTINE rsdft_fft3d( zf, forward_or_backward )
    implicit none
    complex(8),intent(INOUT) :: zf(0:,0:,0:)
    integer,intent(IN) :: forward_or_backward
    complex(8),allocatable :: a(:),b(:),c(:)
    integer :: ml(3),n,i,j,k
    real(8) :: d

    ml(1) = size( zf, 1 )
    ml(2) = size( zf, 2 )
    ml(3) = size( zf, 3 )

    pi2 = 2.0d0*pi
    if ( forward_or_backward < 0 ) pi2=-pi2

    n = maxval( ml )
    allocate( a(0:n-1) ) ; a=z0
    allocate( b(0:n-1) ) ; b=z0
    allocate( c(0:n-1) ) ; c=z0

    do k=0,ml(3)-1
    do j=0,ml(2)-1
       a(:)=z0
       b(:)=z0
       c(:)=z0
       do i=0,ml(1)-1
          a(i) = zf(i,j,k)
       end do
       call rsdft_fft( ml(1), a, b, c )
       do i=0,ml(1)-1
          zf(i,j,k) = b(i)
       end do
    end do ! j
    end do ! k

    do k=0,ml(3)-1
    do i=0,ml(1)-1
       a(:)=z0
       b(:)=z0
       c(:)=z0
       do j=0,ml(2)-1
          a(j) = zf(i,j,k)
       end do
       call rsdft_fft( ml(2), a, b, c )
       do j=0,ml(2)-1
          zf(i,j,k) = b(j)
       end do
    end do ! i
    end do ! k

    do j=0,ml(2)-1
    do i=0,ml(1)-1
       a(:)=z0
       b(:)=z0
       c(:)=z0
       do k=0,ml(3)-1
          a(k) = zf(i,j,k)
       end do
       call rsdft_fft( ml(3), a, b, c )
       do k=0,ml(3)-1
          zf(i,j,k) = b(k)
       end do
    end do ! i
    end do ! j

    if ( forward_or_backward >= 0 ) then
       d = 1.0d0/dble( ml(1)*ml(2)*ml(3) )
       zf(:,:,:) = d*zf(:,:,:)
    end if

    deallocate( c )
    deallocate( b )
    deallocate( a )

  END SUBROUTINE rsdft_fft3d


  SUBROUTINE rsdft_fft( n, a, b, c )
    implicit none
    integer,intent(IN) :: n
    complex(8),intent(INOUT) :: a(0:n-1)
    complex(8),intent(INOUT) :: b(0:)
    complex(8),intent(INOUT) :: c(0:n-1)
    integer :: i,j,k
    integer,allocatable :: info_fact(:),lst_fact(:)
    integer,allocatable :: indx(:)

    allocate( info_fact(2:n) ) ; info_fact=0
    call prime_factor_decomposition( n, info_fact )

    k=sum( info_fact )
    allocate( lst_fact(k) ) ; lst_fact=0

    k=0
    do i=n,2,-1
       do j=1,info_fact(i)
          k=k+1
          lst_fact(k) = i
       end do
    end do

    allocate( indx(0:n-1) ) ; indx=0

    mm=0
    call fft_1_indx( n, lst_fact, indx )
    mm=0
    call fft_1_sub( n, lst_fact, a, b, c )

    c(:)=b(:)
    do i=0,n-1
       j=indx(i)
       b(j)=c(i)
    end do

    deallocate( indx )
    deallocate( lst_fact  )
    deallocate( info_fact )

  END SUBROUTINE rsdft_fft


  SUBROUTINE prime_factor_decomposition( n0, ifact )
    implicit none
    integer,intent(IN)  :: n0
    integer,intent(OUT) :: ifact(2:n0)
    integer :: n,i,j,m,t
    ifact(:)=0
    n=n0
    m=n
    loop_j : do j=1,n0
       loop_i : do i=2,n
          t = mod( m, i )
          if ( t == 0 ) then
             m = m/i
             ifact(i) = ifact(i) + 1
             exit loop_i
          else if ( i == n ) then
             exit loop_j
          end if
       end do loop_i
       n=m
    end do loop_j
  END SUBROUTINE prime_factor_decomposition


  RECURSIVE SUBROUTINE fft_1_sub( n, lst_fact, a, b, c )
    implicit none
    integer,intent(IN) :: n, lst_fact(:)
    complex(8),intent(INOUT) :: a(0:n-1)
    complex(8),intent(INOUT) :: b(0:)
    complex(8),intent(INOUT) :: c(0:)
    integer :: i,j,k,n1,n2

    if ( power_2_check(n) ) then
       call rsdft_fft_2( n, a, b, c )
       return
    end if

    n1 = lst_fact(1)
    n2 = n/n1

    do k=0,n1-1

       c(:)=z0
       do i=0,n2-1
          do j=0,n1-1
             c(i) = c(i) + a(i+j*n2)*exp( -zi*pi2*k/n*(i+j*n2) )
          end do
       end do

       if ( size(lst_fact) == 1 ) then
          do j=0,n2-1
             do i=0,n2-1
                b(mm) = b(mm) + c(i)*exp( -zi*pi2*j/n2*i )
             end do
             mm=mm+1
          end do
       else
          call fft_1_sub( n2, lst_fact(2:), c, b, c(n2:) ) 
       end if

    end do ! k

  END SUBROUTINE fft_1_sub


  RECURSIVE SUBROUTINE fft_1_indx( n, lst_fact, indx )
    implicit none
    integer,intent(IN) :: n, lst_fact(:)
    integer,intent(INOUT) :: indx(0:)
    integer :: j,k,n1,n2,mm0

    if ( power_2_check(n) ) then
       do j=0,n-1
          indx(j+mm) = bit_reverse( j, n )
       end do
       mm=mm+n
       return
    end if

    n1 = lst_fact(1)
    n2 = n/n1

    do k=0,n1-1
       if ( size(lst_fact) == 1 ) then
          do j=0,n2-1
             indx(mm)=j*n1+k
             mm=mm+1
          end do
       else
          mm0=mm
          call fft_1_indx( n2, lst_fact(2:), indx ) 
          do j=mm0,mm-1
             indx(j)=indx(j)*n1+k
          end do
       end if
    end do ! k

  END SUBROUTINE fft_1_indx


  RECURSIVE SUBROUTINE rsdft_fft_2( n, a, b, c )
    implicit none
    integer,intent(IN) :: n
    complex(8),intent(INOUT) :: a(0:n-1)
    complex(8),intent(INOUT) :: b(0:)
    complex(8),intent(INOUT) :: c(0:n-1)
    integer :: i

    do i=0,n/2-1
       c(i) = a(i) + a(i+n/2)
    end do

    if ( n == 2 ) then
       b(mm)=c(0)
       mm=mm+1
    else
       call rsdft_fft_2( n/2, c, b, c(n/2) )
    end if

    do i=0,n/2-1
       c(i) = ( a(i) - a(i+n/2) )*exp( -pi2*zi/n*i )
    end do
    
    if ( n == 2 ) then
       b(mm)=c(0)
       mm=mm+1
    else
       call rsdft_fft_2( n/2, c, b, c(n/2) )
    end if

  END SUBROUTINE rsdft_fft_2


  FUNCTION power_2_check( n )
    implicit none
    logical :: power_2_check
    integer,intent(IN)  :: n
    real(8) :: x
    integer :: i
    x = log(dble(n))/log(2.0d0)
    i = nint(x)
    if ( abs(x-i) < 1.d-10 ) then
       power_2_check = .true.
    else
       power_2_check = .false.
    end if
  END FUNCTION power_2_check


  FUNCTION bit_reverse( i, n )
    implicit none
    integer :: bit_reverse
    integer,intent(IN) :: i, n
    integer :: len_bit,ipos,i_br
    len_bit = nint( log(dble(n))/log(2.0d0) )
    i_br = 0
    do ipos=0,len_bit-1
       call mvbits( i, ipos, 1, i_br, len_bit-1-ipos )
    end do
    bit_reverse = i_br
  END FUNCTION bit_reverse


END MODULE rsdft_fft_module
