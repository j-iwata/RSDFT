MODULE maskf_module

  implicit none

  PRIVATE
  PUBLIC :: nmsk,maskr,xm,dxm,makemaskf

  integer,parameter :: nmsk=201
  real(8) :: maskr(0:nmsk),xm(0:nmsk),dxm

CONTAINS

  SUBROUTINE makemaskf(eta)
    real(8),intent(IN) :: eta
    real(8),allocatable :: X12(:,:),iX12(:,:),lambda(:),rwork(:)
    real(8) :: dx,x1,x2,x,c,pi
    integer,parameter :: matz=1
    integer :: ierr,mxtmp,i,j,LWORK,LIWORK
    integer,allocatable :: iwork(:)
    pi=acos(-1.d0)
    mxtmp=nmsk-1
    dx=1.d0/mxtmp
    allocate( lambda(mxtmp) )
    allocate( X12(mxtmp,mxtmp) )
    do j=1,mxtmp
       x2=j*dx
       do i=j,mxtmp
          x1=i*dx
          if ( i==j ) then
             X12(i,j)=sin((x1+x2)*eta)/(x1+x2)+mxtmp*pi-eta
          else
             X12(i,j)=sin((x1+x2)*eta)/(x1+x2)-sin((x1-x2)*eta)/(x1-x2)
             X12(j,i)=X12(i,j)
          end if
       end do
    end do
    LWORK=1+6*mxtmp+2*mxtmp*mxtmp
    LIWORK=3+5*mxtmp
    allocate( rwork(LWORK),iwork(LIWORK) )
    call DSYEVD('V','L',mxtmp,X12,mxtmp,lambda,rwork,LWORK,iwork,LIWORK,ierr)
    deallocate( iwork,rwork )
    c=X12(1,1)/dx
    do i=1,mxtmp
       x=i*dx
       maskr(i)=X12(i,1)/x/c
    end do
    maskr(0)=1.d0
    maskr(nmsk)=0.d0
    dxm=dx
    do i=0,nmsk
       xm(i)=i*dxm
    end do
    deallocate( lambda )
    deallocate( X12 )
    return
  END SUBROUTINE makemaskf

END MODULE maskf_module
