MODULE minimal_box_module

  implicit none

  PRIVATE
  PUBLIC :: m_grid_ion,map_grid_ion,mcube_grid_ion,make_minimal_box

  integer :: m_grid_ion
  integer,allocatable :: map_grid_ion(:,:)
  integer,allocatable :: mcube_grid_ion(:,:)

CONTAINS

  SUBROUTINE make_minimal_box(Rc,mm1,mm2,mm3,mm4)
    use aa_module
    use rgrid_module
    real(8),intent(IN)  :: Rc
    integer,intent(OUT) :: mm1,mm2,mm3,mm4
    real(8) :: a1,a2,a3,b1,b2,c1,c2,c3,d1,d2,d3,nv(3)
    real(8) :: x,y,z,r,Rmax,R0
    integer :: m,n,m1,m2,m3,mm,mm0,isize,i1,i2,i3,mmmax,j1,j2,j3
    integer :: mt1,mt2,mt3
    integer :: n1a,n1b,n2a,n2b,n3a,n3b
    integer,allocatable :: jtmp(:,:,:),jtmp0(:,:,:)

    R0 = Rc*Rc

    mm1=0
    mm2=0
    mm3=0
    mm4=0

    a1=Hgrid(1)*Ngrid(1)
    a2=Hgrid(2)*Ngrid(2)
    a3=Hgrid(3)*Ngrid(3)
    b1=sum( aa(1:3,1)*aa(1:3,2) )/(a1*a2)
    b2=sqrt(1.d0-b1**2)
    c1=Hgrid(1)*b2
    c2=Hgrid(2)*b2
    mm1=nint(R0/c1)+1
    mm2=nint(R0/c2)+1

    nv(1)=aa(2,1)*aa(3,2)-aa(3,1)*aa(2,2)
    nv(2)=aa(3,1)*aa(1,2)-aa(1,1)*aa(3,2)
    nv(3)=aa(1,1)*aa(2,2)-aa(2,1)*aa(1,2)
    c3=1.d0/sqrt(sum(nv**2))
    nv(1:3)=nv(1:3)*c3
    x =abs(sum(nv(1:3)*aa(1:3,3)))/a3
    c3=x*Hgrid(3)
    mm3=nint(R0/c3)+1

    mm0=(2*mm1+1)*(2*mm2+1)*(2*mm3+1)
    allocate( jtmp(-mm1:mm1,-mm2:mm2,-mm3:mm3) ) ; jtmp=0

    mmmax = 1000
    Rmax  = 0.d0
    m1    = 0
    m2    = 0
    m3    = 0

    n=0

    c1=1.d0/Ngrid(1)
    c2=1.d0/Ngrid(2)
    c3=1.d0/Ngrid(3)

    do mm=1,mmmax

       m=0
       do i3=-m3,m3
          d3=i3*c3
       do i2=-m2,m2
          d2=i2*c2
       do i1=-m1,m1
          d1=i1*c1
          x=d1*aa(1,1)+d2*aa(1,2)+d3*aa(1,3)
          y=d1*aa(2,1)+d2*aa(2,2)+d3*aa(2,3)
          z=d1*aa(3,1)+d2*aa(3,2)+d3*aa(3,3)
          r=x*x+y*y+z*z
          if ( r <= R0 ) then
             m=m+1
             Rmax=max(Rmax,r)
             jtmp(i1,i2,i3)=jtmp(i1,i2,i3)+1
          end if
       end do
       end do
       end do

       if ( m == n ) exit
       n=m

       m1=m1+1
       m2=m2+1
       m3=m3+1

       if ( m1==mm1 .or. m2==mm2 .or. m3==mm3 ) then
          allocate( jtmp0(-mm1:mm1,-mm2:mm2,-mm3:mm3) )
          jtmp0(:,:,:)=jtmp(:,:,:)
          mt1=mm1 ; if ( m1==mm1 ) mt1=mt1+1
          mt2=mm2 ; if ( m2==mm2 ) mt2=mt2+1
          mt3=mm3 ; if ( m3==mm3 ) mt3=mt3+1
          deallocate( jtmp )
          allocate( jtmp(-mt1:mt1,-mt2:mt2,-mt3:mt3) )
          jtmp=0
          jtmp(-mm1:mm1,-mm2:mm2,-mm3:mm3) &
               = jtmp0(-mm1:mm1,-mm2:mm2,-mm3:mm3)
          mm1=mt1
          mm2=mt2
          mm3=mt3
          deallocate( jtmp0 )
       end if

    end do ! mm

    Rmax=sqrt(Rmax)

    isize = count(jtmp>0)

!- allocate ------------------------------------------------
    if ( allocated(map_grid_ion) ) then
       deallocate( map_grid_ion )
    end if
    if ( allocated(mcube_grid_ion) ) then
       deallocate( mcube_grid_ion )
    end if
    allocate( map_grid_ion(3,isize)  ) ; map_grid_ion=0
    allocate( mcube_grid_ion(2,3)    ) ; mcube_grid_ion=0
!-----------------------------------------------------------

    n=0
    do i3=-mm3,mm3
    do i2=-mm2,mm2
    do i1=-mm1,mm1
       if ( jtmp(i1,i2,i3)>0 ) then
          n=n+1
          map_grid_ion(1,n)=i1
          map_grid_ion(2,n)=i2
          map_grid_ion(3,n)=i3
       end if
    end do
    end do
    end do

    if( n/=isize ) then
       write(*,*) "i/=isize!!!",n,isize
       stop
    end if

    m_grid_ion = n

    mcube_grid_ion(1,1)=minval( map_grid_ion(1,1:m_grid_ion) )
    mcube_grid_ion(1,2)=minval( map_grid_ion(2,1:m_grid_ion) )
    mcube_grid_ion(1,3)=minval( map_grid_ion(3,1:m_grid_ion) )
    mcube_grid_ion(2,1)=maxval( map_grid_ion(1,1:m_grid_ion) )
    mcube_grid_ion(2,2)=maxval( map_grid_ion(2,1:m_grid_ion) )
    mcube_grid_ion(2,3)=maxval( map_grid_ion(3,1:m_grid_ion) )

    deallocate( jtmp )

    mm1 = maxval( abs(mcube_grid_ion(:,1)) ) + 1
    mm2 = maxval( abs(mcube_grid_ion(:,2)) ) + 1
    mm3 = maxval( abs(mcube_grid_ion(:,3)) ) + 1
    mm4 = M_grid_ion

  END SUBROUTINE make_minimal_box

END MODULE minimal_box_module
