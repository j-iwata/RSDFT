SUBROUTINE dot_product(a,b,c,alpha,n,m)
  integer,intent(IN) :: n,m
  integer :: i
  real(8) :: a(*),b(*),c(*),alpha,tmp

  c(1:m)=0.d0

  tmp=0.d0
!$OMP parallel do reduction(+:tmp)
  do i=1,n
     tmp=tmp+a(i)*b(i)
  end do
!$OMP end parallel do
  c(1)=tmp*alpha

  if ( m==2 ) then
     tmp=0.d0
!$OMP parallel do reduction(+:tmp)
     do i=1,n-1,2
        tmp=tmp+a(i)*b(i+1)-a(i+1)*b(i)
     end do
!$OMP end parallel do
     c(m)=tmp*alpha
  end if

  return
END SUBROUTINE dot_product
