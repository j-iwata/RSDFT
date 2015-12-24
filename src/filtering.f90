MODULE Filtering
!----------------------------------------------------------------------------
! this module is used for filtering radial part of PS
!----------------------------------------------------------------------------
  implicit none
  PRIVATE
  PUBLIC :: opFiltering

CONTAINS

!------------------------------------
  SUBROUTINE opFiltering( qc,L,NRc,NRps,rad,rad1,vrad,tarIN )
    implicit none

    real(8),intent(IN) :: qc
    integer,intent(IN) :: L,NRc,NRps
    real(8),intent(IN) :: rad(NRc),rad1(NRps),vrad(NRc)
    real(8),intent(INOUT) :: tarIN(NRps)
    integer :: i,j
    real(8) :: r,r1
    real(8) :: sb0x,sb1x,sb0y,sb1y,sum0
    real(8),allocatable :: tmp(:)
    real(8),parameter :: const=2.d0/acos(-1.d0) ! HERE

    allocate( tmp(NRc) ) ; tmp(:)=0.d0
    do i=1,NRps
      r=rad1(i)
       select case(L)
       case(0)
          if ( r==0.d0 ) then
             r1=rad(1)             !-- log_mesh points 
             if ( r1==0.d0 ) then
                tmp(1)=qc*qc*qc/3.d0
             else
                tmp(1)=sin(qc*r1)/(r1*r1*r1)-qc*cos(qc*r1)/(r1*r1)
             end if
             do j=2,NRc
                r1=rad(j)         !-- log_mesh points  
                tmp(j)=sin(qc*r1)/(r1*r1*r1)-qc*cos(qc*r1)/(r1*r1)
             end do
          else
             do j=1,NRc
                r1=rad(j)
                if ( r1==0.d0 ) then
                   tmp(j)=sin(qc*r)/(r*r*r)-qc*cos(qc*r)/(r*r)
                else if ( r1==r ) then
                   tmp(j)=(2.d0*qc*r-sin(2.d0*qc*r))/(4.d0*r*r*r)
                else
                   tmp(j)=( sin(qc*(r-r1))/(r-r1)-sin(qc*(r+r1))/(r+r1) ) &
                        /(2.d0*r*r1)
                end if
             end do
          end if
       case(1)
          if ( r==0.d0 ) then
             tarIN(i)=0.d0
             cycle ! i
          else
             do j=1,NRc
                r1=rad(j)
                if ( r1==0.d0 ) then
                   tmp(j)=0.d0
                else if ( r1==r ) then
                   sb0x=sin(qc*r)/(qc*r)
                   sb1x=sb0x/(qc*r)-cos(qc*r)/(qc*r)
                   tmp(j)=(2.d0*qc*r-sin(2.d0*qc*r))/(4.d0*r*r*r) &
                        -qc*qc*sb0x*sb1x/r
                else
                   sb0x=sin(qc*r)/(qc*r)
                   sb0y=sin(qc*r1)/(qc*r1)
                   sb1x=sb0x/(qc*r)-cos(qc*r)/(qc*r)
                   sb1y=sb0y/(qc*r1)-cos(qc*r1)/(qc*r1)
                   tmp(j)=( r1*sb0y*sb1x-r*sb0x*sb1y )*qc*qc/(r*r-r1*r1)
                end if
             end do
          end if
       case(2)
          if ( r==0.d0 ) then
             tarIN(i)=0.d0
             cycle ! i
          else
             do j=1,NRc
                r1=rad(j)
                if ( r1==0.d0 ) then
                   tmp(j)=0.d0
                else if ( r1==r ) then
                   sb1x=sin(qc*r)/(qc*qc*r*r)-cos(qc*r)/(qc*r)
                   tmp(j)=(2.d0*qc*r-sin(2.d0*qc*r))/(4.d0*r*r*r) &
                        -3.d0*qc*sb1x*sb1x/(r*r)
                else
                   sb0x=sin(qc*r)/(qc*r)
                   sb0y=sin(qc*r1)/(qc*r1)
                   sb1x=sb0x/(qc*r)-cos(qc*r)/(qc*r)
                   sb1y=sb0y/(qc*r1)-cos(qc*r1)/(qc*r1)
                   tmp(j)=( r*sb0y*sb1x-r1*sb0x*sb1y )*qc*qc/(r*r-r1*r1) &
                        -3.d0*qc/(r*r1)*sb1x*sb1y
                end if
             end do
          end if
       case(3)
          if ( r==0.d0 ) then
             tarIN(i)=0.d0
             cycle ! i
          else
             do j=1,NRc
                r1=rad(j)
                if ( r1==0.d0 ) then
                   tmp(j)=0.d0
                else if ( r1==r ) then
                   sb0x=sin(qc*r)/(qc*r)
                   sb1x=sin(qc*r)/(qc*qc*r*r)-cos(qc*r)/(qc*r)
                   tmp(j)=(2.d0*qc*r-sin(2.d0*qc*r))/(4.d0*r*r*r) &
                        -3.d0*qc*sb1x*sb1x/(r*r) &
                        -(qc**2)/r*(  3.d0/(qc*r)*(15.d0/((qc*r)**2)-1.d0 ) &
                        *sb1x*sb1x+5.d0/(qc*r)*sb0x*sb0x &
                        +(1.d0-30.d0/(qc*r)**2.d0)*sb0x*sb1x)
                else
                   sb0x=sin(qc*r)/(qc*r)
                   sb0y=sin(qc*r1)/(qc*r1)
                   sb1x=sb0x/(qc*r)-cos(qc*r)/(qc*r)
                   sb1y=sb0y/(qc*r1)-cos(qc*r1)/(qc*r1)
                   tmp(j)= sb0y*sb1x*( 15.d0/(r1*r*r)+(qc**2.d0)*r1/(r**2.d0-r1**2.d0)) &
                        +sb0x*sb1y*( 15.d0/(r*r1*r1)-(qc**2.d0)*r/(r**2.d0-r1**2.d0)) &
                        -45.d0/(qc*r**2.d0*r1**2.d0)*sb1x*sb1y-5.d0*qc/(r*r1)*sb0x*sb0y
                end if
             end do
          end if
       case(4)
          if ( r==0.d0 ) then
             tarIN(i)=0.d0
             cycle ! i
          else
             do j=1,NRc
                r1=rad(j)
                if ( r1==0.d0 ) then
                   tmp(j)=0.d0
                else if ( r1==r ) then
                   sb0x=sin(qc*r)/(qc*r)
                   sb1x=sin(qc*r)/(qc*qc*r*r)-cos(qc*r)/(qc*r)
                   tmp(j)=(2.d0*qc*r-sin(2.d0*qc*r))/(4.d0*r*r*r) &
                        -3.d0*qc*sb1x*sb1x/(r*r) &
                        -(qc**2.d0)/r*(  3.d0/(qc*r)*(15.d0/((qc*r)**2.d0)-1.d0 ) &
                        *sb1x*sb1x+5.d0/(qc*r)*sb0x*sb0x &
                        +(1.d0-30.d0/(qc*r)**2.d0)*sb0x*sb1x) &
                        -(qc**2)/r*( (1575.d0/(qc*r)**5-255.d0/(qc*r)**3.d0 &
                        +10.d0/(qc*r))*sb1x*sb1x+(175.d0/(qc*r)**3.d0 &
                        -5.d0/(qc*r))*sb0x*sb0x+(-1050.d0/(qc*r)**4.d0 &
                        +100.d0/(qc*r)**2.d0-1.d0)*sb0x*sb1x) 
                else
                   sb0x=sin(qc*r)/(qc*r)
                   sb0y=sin(qc*r1)/(qc*r1)
                   sb1x=sb0x/(qc*r)-cos(qc*r)/(qc*r)
                   sb1y=sb0y/(qc*r1)-cos(qc*r1)/(qc*r1)
                   tmp(j)= sb0y*sb1x*(-35.d0/(r1*r1*r) &
                        +(qc**2.d0)*r/(r**2.d0-r1**2.d0) &
                        +525.d0/((qc**2.d0)*(r1**2.d0)*(r**3.d0)) ) &
                        +sb0x*sb1y*(-35.d0/(r*r*r1) &
                        +(qc**2.d0)*r1/(r1**2.d0-r**2.d0) &
                        +525.d0/((qc**2.d0)*(r**2.d0)*(r1**3.d0)) ) &
                        +sb1x*sb1y*(   105.d0*(r**2.d0+r1**2.d0) &
                        /(qc*(r**3.d0)*(r1**3.d0)) &
                        -10.d0*qc/(r*r1) &
                        -1575.d0/((qc**3.d0)*(r**3.d0)*(r1**3.d0))) &
                        +sb0x*sb0y*( -175.d0/(qc*(r**2.d0)*(r1**2.d0)) )
                end if
             end do
          end if
       case default
          stop "opFiltering implementation L"
       end select
       tmp(1:NRc)=tmp(1:NRc)*vrad(1:NRc)
       call simp(tmp(1:NRc),sum0,NRc,2)
       tarIN(i)=sum0*const
    end do	! i

    return
  END SUBROUTINE opFiltering

!----------------------------------

  SUBROUTINE simp(f,s,n,m)
    implicit none
    integer,intent(IN)  :: n,m
    real(8),intent(IN)  :: f(n)
    real(8),intent(OUT) :: s
    real(8),allocatable :: g(:)
    integer :: i,nmax
    nmax=int(n/m)*m
    do i=0,m
       nmax=nmax+i ; if ( nmax>=n ) exit
    end do
    allocate( g(nmax) ) ; g(1:n)=f ; if ( nmax>n ) g(n+1:)=0.d0
    select case(m)
    case default
       s = 0.5d0*(f(1)+f(n)) + sum(f(2:n-1))
    case(2)
       s=0.d0
       do i=1,nmax-2,2
          s = s + g(i) + 4.d0*g(i+1) + g(i+2)
       end do
       s=s/3.d0
    case(4)
       s=0.d0
       do i=1,nmax-4,4
          s=s+7*g(i)+32*g(i+1)+12*g(i+2)+32*g(i+3)+7*g(i+4)
       end do
       s=s*2.d0/45.d0
    case(6)
       s=0.d0
       do i=1,nmax-6,6
          s=s+41*g(i)+216*g(i+1)+27*g(i+2)+272*g(i+3) &
               +27*g(i+4)+216*g(i+5)+41*g(i+6)
       end do
       s=s/140.d0
    end select
    deallocate( g )
    return
  END SUBROUTINE simp

END MODULE Filtering
