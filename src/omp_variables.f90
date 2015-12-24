MODULE omp_variables

!$ use omp_lib

  implicit none

  PRIVATE
  PUBLIC :: init_omp, Igrid_omp,Ngrid_omp
  PUBLIC :: partition_job_omp

  integer,allocatable :: Igrid_omp(:,:,:)
  integer,allocatable :: Ngrid_omp(:,:)
  integer :: nthreads

CONTAINS


  SUBROUTINE init_omp(a1b,b1b,a2b,b2b,a3b,b3b,n1,n2,systype,disp_switch)

    implicit none

    integer,intent(IN) :: a1b,b1b,a2b,b2b,a3b,b3b,n1,n2,systype
    logical,intent(IN) :: disp_switch
    integer :: i1,i2,i3,m,n,k,j,i,ab1,ab2,ab3,nn(3),np(3)
    integer,allocatable :: ntmp(:,:)

    call write_border( 80, " init_omp(start)" )

    ab1 = b1b-a1b+1
    ab2 = b2b-a2b+1
    ab3 = b3b-a3b+1

    m=1
    n=1

!$OMP parallel private(m,n,k,j,i)

!$  n=omp_get_num_threads()

!$OMP single

    nthreads = n
    allocate( Ngrid_omp(0:3,0:n-1)   ) ; Ngrid_omp=0
    allocate( Igrid_omp(2,0:3,0:n-1) ) ; Igrid_omp=0

    m = n

    k     = gcd(ab3,m)
    m     = m/k
    nn(3) = ab3/k
    np(3) = k

    if ( m == 1 ) then

       nn(2) = ab2
       nn(1) = ab1
       np(2) = 1
       np(1) = 1

    else

       k     = gcd(ab2,m)
       m     = m/k
       nn(2) = ab2/k
       np(2) = k

       if ( m == 1 ) then

          nn(1) = ab1
          np(1) = 1

       else

          k     = gcd(ab1,m)
          m     = m/k
          nn(1) = ab1/k
          np(1) = k

!          if ( m == 1 ) then
!          else
!             nn(1) = ab1
!          end if

       end if

    end if

    if ( disp_switch ) then
       write(*,*) "nn(ini)=",ab1,ab2,ab3
       write(*,*) "nn(1:3)=",nn(1:3)
       write(*,*) "np(1:3)=",np(1:3)
       write(*,*) "m      =",m
    end if

    if ( m == 1 ) then

       n=maxval(np)
       allocate( ntmp(3,n) ) ; ntmp=0

       do i=1,ab1
          n=mod(i-1,np(1))+1
          ntmp(1,n)=ntmp(1,n)+1
       end do
       do i=1,ab2
          n=mod(i-1,np(2))+1
          ntmp(2,n)=ntmp(2,n)+1
       end do
       do i=1,ab3
          n=mod(i-1,np(3))+1
          ntmp(3,n)=ntmp(3,n)+1
       end do

       n=-1
       do k=1,np(3)
       do j=1,np(2)
       do i=1,np(1)
          n=n+1
          Ngrid_omp(1,n)=ntmp(1,i)
          Ngrid_omp(2,n)=ntmp(2,j)
          Ngrid_omp(3,n)=ntmp(3,k)
          if ( disp_switch ) write(*,'(1x,"Ngrid_omp",3i5)') Ngrid_omp(1:3,n)
       end do
       end do
       end do

    else

       np( maxloc(nn,1) ) = np( maxloc(nn,1) )*m

       if ( disp_switch ) write(*,*) "np_=",np

       n=maxval(np)
       allocate( ntmp(3,n) ) ; ntmp=0

       do i=1,ab1
          n=mod(i-1,np(1))+1
          ntmp(1,n)=ntmp(1,n)+1
       end do
       do i=1,ab2
          n=mod(i-1,np(2))+1
          ntmp(2,n)=ntmp(2,n)+1
       end do
       do i=1,ab3
          n=mod(i-1,np(3))+1
          ntmp(3,n)=ntmp(3,n)+1
       end do

       n=-1
       do k=1,np(3)
       do j=1,np(2)
       do i=1,np(1)
          n=n+1
          Ngrid_omp(1,n)=ntmp(1,i)
          Ngrid_omp(2,n)=ntmp(2,j)
          Ngrid_omp(3,n)=ntmp(3,k)
          if ( disp_switch ) write(*,'(1x,"Ngrid_omp_",3i5)') Ngrid_omp(1:3,n)
       end do
       end do
       end do

    end if

    if ( systype == 1 ) then

       do j=1,n2-n1+1
          i=mod( j-1, nthreads )
          Ngrid_omp(0,i) = Ngrid_omp(0,i) + 1
       end do

    else

       do i=0,nthreads-1
          Ngrid_omp(0,i) = Ngrid_omp(1,i)*Ngrid_omp(2,i)*Ngrid_omp(3,i)
       end do

    end if

    i=-1
    do i3=1,np(3)
    do i2=1,np(2)
    do i1=1,np(1)
       i=i+1
       Igrid_omp(1,0,i) = sum( Ngrid_omp(0,0:i) ) - Ngrid_omp(0,i) + n1
       Igrid_omp(2,0,i) = Igrid_omp(1,0,i) + Ngrid_omp(0,i) - 1
       Igrid_omp(1,1,i) = sum(ntmp(1,1:i1))-ntmp(1,i1) + a1b
       Igrid_omp(2,1,i) = Igrid_omp(1,1,i) + ntmp(1,i1) - 1
       Igrid_omp(1,2,i) = sum(ntmp(2,1:i2))-ntmp(2,i2) + a2b
       Igrid_omp(2,2,i) = Igrid_omp(1,2,i) + ntmp(2,i2) - 1
       Igrid_omp(1,3,i) = sum(ntmp(3,1:i3))-ntmp(3,i3) + a3b
       Igrid_omp(2,3,i) = Igrid_omp(1,3,i) + ntmp(3,i3) - 1
       if ( disp_switch ) then
          write(*,'(1x,"IGIG",1x,9i8)') Igrid_omp(:,:,i),Ngrid_omp(0,i)
       end if
    end do
    end do
    end do

    if ( disp_switch ) write(*,*) "omp_get_num_threads=",nthreads

    deallocate( ntmp )

!$OMP end single

!$OMP end parallel

    call write_border( 80, " init_omp(end)" )

  END SUBROUTINE init_omp


  FUNCTION gcd(m0,n0)
    implicit none
    integer :: gcd,m0,n0
    integer :: m,n,mtmp,loop

    if ( m0 >= n0 ) then
       m=m0
       n=n0
    else
       m=n0
       n=m0
    end if

    do loop=1,10000
       if ( n == 0 ) exit
       mtmp = n
       n = mod(m,n)
       m = mtmp
    end do

    gcd = m

  END FUNCTION gcd


  SUBROUTINE partition_job_omp( a1b,b1b,a2b,b2b,a3b,b3b,ib1,ib2 )

    implicit none
    integer,intent(INOUT) :: a1b,b1b,a2b,b2b,a3b,b3b,ib1,ib2
    integer :: i1,i2,i3,i4,m,n,k,j,i,ab1,ab2,ab3,nib,nn(4),np(4),myrank_omp
    integer,allocatable :: ntmp(:,:)

    ab1 = b1b-a1b+1
    ab2 = b2b-a2b+1
    ab3 = b3b-a3b+1
    nib = ib2-ib1+1

    m=1
    n=1

!$  n = omp_get_num_threads()
!$  myrank_omp = omp_get_thread_num()

    m = n

    k     = gcd(nib,m)
    m     = m/k
    nn(4) = nib/k
    np(4) = k

    if ( m == 1 ) then

       nn(3) = ab3
       nn(2) = ab2
       nn(1) = ab1
       np(3) = 1
       np(2) = 1
       np(1) = 1

    else

       k     = gcd(ab3,m)
       m     = m/k
       nn(3) = ab3/k
       np(3) = k

       if ( m == 1 ) then

          nn(2) = ab2
          nn(1) = ab1
          np(2) = 1
          np(1) = 1

       else

          k     = gcd(ab2,m)
          m     = m/k
          nn(2) = ab2/k
          np(2) = k

          if ( m == 1 ) then

             nn(1) = ab1
             np(1) = 1

          else

             k     = gcd(ab1,m)
             m     = m/k
             nn(1) = ab1/k
             np(1) = k

          end if

       end if

    end if

!    write(*,*) "nn(ini)=",ab1,ab2,ab3,nib
!    write(*,*) "nn(1:4)=",nn(1:4)
!    write(*,*) "np(1:4)=",np(1:4)
!    write(*,*) "m      =",m

    if ( m == 1 ) then

       i=-1
       loop_a : do i4=0,np(4)-1 
       do i3=0,np(3)-1 
       do i2=0,np(2)-1 
       do i1=0,np(1)-1 
          i=i+1
          if ( i == myrank_omp ) then
             a1b = a1b + nn(1)*i1
             b1b = a1b + nn(1) - 1
             a2b = a2b + nn(2)*i2
             b2b = a2b + nn(2) - 1
             a3b = a3b + nn(3)*i3
             b3b = a3b + nn(3) - 1
             ib1 = ib1 + nn(4)*i4
             ib2 = ib1 + nn(4) - 1
             exit loop_a
          end if
       end do
       end do
       end do
       end do loop_a

    else

       np( maxloc(nn,1) ) = np( maxloc(nn,1) )*m

       n=maxval(np)
       allocate( ntmp(4,0:n-1) ) ; ntmp=0

       do i=1,ab1
          n=mod(i-1,np(1))
          ntmp(1,n)=ntmp(1,n)+1
       end do
       do i=1,ab2
          n=mod(i-1,np(2))
          ntmp(2,n)=ntmp(2,n)+1
       end do
       do i=1,ab3
          n=mod(i-1,np(3))
          ntmp(3,n)=ntmp(3,n)+1
       end do
       do i=1,nib
          n=mod(i-1,np(4))
          ntmp(4,n)=ntmp(4,n)+1
       end do

       i=-1
       loop_b : do i4=0,np(4)-1
       do i3=0,np(3)-1
       do i2=0,np(2)-1
       do i1=0,np(1)-1
          i=i+1
          if ( i == myrank_omp ) then
             a1b = a1b + sum( ntmp(1,1:i1) ) - ntmp(1,i1)
             b1b = a1b + ntmp(1,i1) - 1
             a2b = a2b + sum( ntmp(2,1:i2) ) - ntmp(1,i2)
             b2b = a2b + ntmp(1,i2) - 1
             a3b = a3b + sum( ntmp(3,1:i3) ) - ntmp(1,i3)
             b3b = a3b + ntmp(1,i3) - 1
             ib1 = ib1 + sum( ntmp(4,1:i4) ) - ntmp(1,i4)
             ib2 = ib1 + ntmp(1,i4) - 1
             exit loop_b
          end if
       end do
       end do
       end do
       end do loop_b

       deallocate( ntmp )

    end if

  END SUBROUTINE partition_job_omp


END MODULE omp_variables
