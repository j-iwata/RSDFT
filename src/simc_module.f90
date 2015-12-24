MODULE simc_module

  use bberf_module

  implicit none

  PRIVATE
  PUBLIC :: simc

  real(8),allocatable :: rads(:),vins(:),wgt(:)
  real(8) :: zvs

CONTAINS


  SUBROUTINE simc(rad,vin,rc,zv,parloc,mesh)
    implicit none
    integer,intent(IN) :: mesh
    real(8),intent(IN) :: rad(mesh),vin(mesh),zv,rc
    real(8),intent(OUT) :: parloc(4)
    call simc_0(rad,vin,rc,zv,parloc,mesh)
    !call simc_1(rad,vin,rc,zv,parloc,mesh)
  END SUBROUTINE simc


  SUBROUTINE simc_1( rad, vin, rc, zv, parloc, mesh )

    implicit none
    integer,intent(IN) :: mesh
    real(8),intent(IN) :: rad(mesh),vin(mesh),zv,rc
    real(8),intent(OUT) :: parloc(4)
    integer,parameter :: lwa0=15, nummin=6
    integer :: info,k,num,ipvt(3),maxnum, lwa
    real(8) :: x(3), nxtsmp, pi
    real(8),allocatable :: wa(:),fvec(:),fjac(:,:)
    real(8),parameter :: lambda=3.5d0
    real(8),parameter :: x1ini=1.0d0, x2ini=0.4d0, x3ini=0.6d0
    real(8),parameter :: tol=1.0d-5, smpstp=0.2d0
    real(8),parameter :: rmax=10.0d0, vmax=100.0d0, vrzmin=3.0d-6

    allocate( rads(mesh)   ) ; rads=0.0d0
    allocate( vins(mesh)   ) ; vins=0.0d0
    allocate( wgt(mesh)    ) ; wgt=0.0d0
    allocate( fvec(mesh)   ) ; fvec=0.0d0
    allocate( fjac(mesh,3) ) ; fjac=0.0d0

    pi = 4.0d0*atan(1.0d0)

    num=0
    nxtsmp=0.0d0

    do k=1,mesh

       if ( rad(k) > nxtsmp ) then

          nxtsmp = nxtsmp + smpstp

          if ( abs(vin(k)) <= vmax ) then
             num = num+1
             rads(num) = rad(k)
             vins(num) = vin(k)
             wgt(num)  = 1.0d0-exp(-(1.2d0*rad(k)/rc)**lambda)
          end if
          if ( (abs(vin(k)*rad(k)+zv) < vrzmin .or. &
                rad(k) > rmax ) .and. num > nummin ) exit

       end if

    end do ! k

    maxnum = num
    lwa = lwa0 + maxnum

    allocate( wa(lwa) ) ; wa=0.0d0

    zvs  = zv
    x(1) = x1ini
    x(2) = x2ini
    x(3) = x3ini

    call levenberg_marquardt(num,3,x,fvec,fjac)

    if ( x(2) < 0.0d0 .or. x(3) < 0.0d0 ) then
       stop "simc: illegally converged."
    end if

    parloc(1) = x(1)
    parloc(2) = x(2)
    parloc(3) = 1.0d0 - x(1)
    parloc(4) = x(3)

    deallocate( wa   )
    deallocate( fjac )
    deallocate( fvec )
    deallocate( wgt  )
    deallocate( vins )
    deallocate( rads )

    return
  END SUBROUTINE simc_1


  SUBROUTINE uscfit_1( m, n, x, fvec, fjac, ldfjac, iflag )
    implicit none
    integer,intent(IN) :: m,n,ldfjac,iflag
    real(8),intent(OUT) :: fvec(m),fjac(ldfjac,n)
    real(8),intent(INOUT) :: x(n)
    real(8) :: pi,xxxx
    integer :: i

    pi = 4.0d0*atan(1.0d0)
    xxxx = 1.0d0 - x(1)

    if ( x(2) < 0.0d0 ) x(2)=0.0d0
    if ( x(3) < 0.0d0 ) x(3)=0.0d0

    if ( iflag == 1 ) then

       do i=1,m
          fvec(i) = (- zvs/rads(i)*( x(1)*bberf(sqrt(x(2))*rads(i)) &
                                   + xxxx*bberf(sqrt(x(3))*rads(i)) &
                                   ) - vins(i) )*sqrt(wgt(i))
       end do

    else if ( iflag == 2 ) then

       do i=1,m
          fjac(i,1) = -zvs/rads(i)*( bberf(sqrt(x(2)*rads(i))) &
                                   - bberf(sqrt(x(3)*rads(i))) &
                                   )*sqrt(wgt(i))
          fjac(i,2) = -zvs/sqrt(pi*x(2))*x(1) &
                      *exp(-x(2)*rads(i)**2)*sqrt(wgt(i))
          fjac(i,3) = -zvs/sqrt(pi*x(3))*(1.0d0-x(1)) &
                      *exp(-x(3)*rads(i)**2)*sqrt(wgt(i))
       end do

    else
       stop "Error in vlfit: iflag must be 1 or 2."
    end if

    return
  END SUBROUTINE uscfit_1


  SUBROUTINE levenberg_marquardt(ndim,nvec,x,fvec,fjac)

    implicit none
    integer,intent(IN) :: ndim,nvec
    real(8),intent(INOUT) :: x(nvec)
    real(8),intent(INOUT) :: fvec(ndim),fjac(ndim,nvec)
    integer,parameter :: max_loop=10000, max_loop0=20
    integer,allocatable :: ipiv(:)
    integer :: m,n,ierr,loop,loop0,num_conv
    real(8),parameter :: delta=1.d-8
    real(8),allocatable :: Hes(:,:),dJ(:),du(:),xtmp(:),ftmp(:),xmin(:)
    real(8) :: c,J,Jtmp,err,Jmin,errmin
    character(80) :: error_message

    allocate( Hes(nvec,nvec) ) ; Hes=0.0d0
    allocate( dJ(nvec)       ) ; dJ=0.0d0
    allocate( du(nvec)       ) ; du=0.0d0
    allocate( ipiv(nvec)     ) ; ipiv=0
    allocate( xtmp(nvec)     ) ; xtmp=0.0d0
    allocate( ftmp(ndim)     ) ; ftmp=0.0d0
    allocate( xmin(ndim)     ) ; xmin=0.0d0

    xmin = x
    Jmin = 1.d100
    num_conv = 0
    errmin = 1.d100

    do loop0=1,max_loop0

       c=0.0001d0

       call uscfit_1( ndim, nvec, x, fvec, fjac, ndim, 1 )
       J=sum( fvec(:)**2 )
       call uscfit_1( ndim, nvec, x, fvec, fjac, ndim, 2 )

       do loop=1,max_loop

          do n=1,nvec
             do m=1,nvec
                Hes(m,n) = sum( fjac(:,m)*fjac(:,n) )
             end do
             Hes(n,n) = Hes(n,n) + c*Hes(n,n)
          end do
          do n=1,nvec
             dJ(n) = sum( fvec(:)*fjac(:,n) )
          end do

          du(:) = -dJ(:)
          call dgesv(nvec,1,Hes,nvec,ipiv,du,nvec,ierr)

          xtmp(:)=x(:)+du(:)
          call uscfit_1( ndim, nvec, xtmp, ftmp, fjac, ndim, 1 )
          Jtmp=sum( ftmp(:)**2 )

          err = sum( du(:)**2 )
          errmin = min( err, errmin )
          if ( err < delta ) then
             num_conv = num_conv + 1
             exit
          end if

          if ( Jtmp > J ) then
             c=10.0d0*c
          else if ( Jtmp <= J ) then
             c=0.1d0*c
             J=Jtmp
             x(:)=xtmp(:)
             fvec(:)=ftmp(:)
             call uscfit_1( ndim, nvec, x, fvec, fjac, ndim, 2 )
          end if

       end do ! loop

!       if ( loop <= max_loop ) then
          if ( Jtmp < Jmin ) then
             Jmin = Jtmp
             xmin = xtmp
          end if
!       end if

       call random_number(x)

    end do ! loop0

    if ( num_conv == 0 ) then
       write(error_message,*) "Jmin,errmin=",Jmin,errmin
       call write_string( error_message )
!       stop "error@levenberg_marquardt"
    end if

    x = xmin

! ---
    deallocate( ftmp )
    deallocate( xtmp )
    deallocate( ipiv )
    deallocate( du   )
    deallocate( dJ   )
    deallocate( Hes  )
  END SUBROUTINE levenberg_marquardt

!#ifdef _TEST_
  SUBROUTINE simc_0(rad,vin,rc,zv,parloc,mesh)
!     $Id: simc.F,v 1.2 1997/06/25 05:07:08 skimu Exp $
!
!     simc: local potential generator
!           MINPACK version
!
    implicit none
    integer :: mesh
    real(8) :: rad(mesh),vin(mesh),parloc(4),zv,rc

    integer,parameter :: maxnum = 2000

    common /parusc/ rads,vins,zvs,wgt
    real*8 rads(maxnum),vins(maxnum),zvs,wgt(maxnum)

    real(8) :: lambda = 3.5d0
    real(8) :: pi

    integer :: k,num
    real(8) :: nxtsmp
    real(8) :: x(3),fvec(maxnum),fjac(maxnum,3)
    integer :: info
    real(8) :: x1ini,x2ini,x3ini
    parameter (x1ini=1.0d0, x2ini=0.4d0, x3ini=0.6d0)
    real(8),parameter :: tol=1.0d-5
    real(8) :: rmax,vmax,vrzmin,smpstp
    integer :: nummin
    parameter(rmax=10.0d0, vmax=100.d0, vrzmin=3.0d-6)
    parameter(smpstp=0.2d0, nummin=6)
    integer :: ipvt(3)
    integer :: lwa
    parameter(lwa=5*3+maxnum)
    real(8) :: wa(lwa)

    pi = 4.0d0*atan(1.0d0)

    num=0
    nxtsmp=0.d0
    do k=1,mesh
       if ( rad(k)>nxtsmp ) then
          nxtsmp = nxtsmp + smpstp
          if ( abs(vin(k)) <= vmax ) then
!**   1.0d0*rc ---- 1.5d0*rc in original
             num=num+1
             if ( num > maxnum ) then
                write(6,*) 'simc: Too many sample points.'
                stop
             end if
             rads(num)=rad(k)
             vins(num)=vin(k)
             wgt(num)=1.d0-dexp(-(1.2d0*rad(k)/rc)**lambda)
          end if
          if ( (abs(vin(k)*rad(k)+zv)<vrzmin .or. rad(k)>rmax) &
               .and. num>nummin) exit
       end if
    end do
    zvs = zv
    x(1) = x1ini
    x(2) = x2ini
    x(3) = x3ini
    call lmder1(uscfit,num,3,x,fvec,fjac,maxnum,tol,info,ipvt,wa,lwa)
!    write(6,*) 'lmder1:',info
    if ( info==0 .or. info==4 .or. info==5 .or. info==6 .or. info==7 ) then
       write(6,*) 'simc: Not converged.'
       write(6,*) 'x(1) = ',x(1)
       write(6,*) 'x(2) = ',x(2)
       write(6,*) 'x(3) = ',x(3)
       write(6,'(A5,3A20)') 'k','rads(k)','vins(k)','wgt(k)'
       do k=1,num
          write(6,'(I5,3g20.12)')k,rads(k),vins(k),wgt(k)
       end do
       stop
    end if
    if ( x(2)<0.0d0 .or. x(3)<0.0d0 ) then
       write(6,*)'simc: illegally converged.'
       stop
    end if
    parloc(1) = x(1)
    parloc(2) = x(2)
    parloc(3) = 1.0d0 - x(1)
    parloc(4) = x(3)
    return
  END SUBROUTINE simc_0
!
!     fitting function for simc
!
  SUBROUTINE uscfit(m,n,x,fvec,fjac,ldfjac,iflag)
    implicit none
    integer :: m,n,ldfjac,iflag
    real(8) :: x(n),fvec(m),fjac(ldfjac,n)
    integer,parameter :: maxnum=2000
    common /parusc/ rad,vin,zv,wgt
    real(8) :: rad(maxnum),vin(maxnum),zv,wgt(maxnum)
    real(8) :: pi
    integer :: i

    pi = 4.0d0*atan(1.0d0)

    if ( x(2) < 0.0d0 ) x(2)=0.0d0
    if ( x(3) < 0.0d0 ) x(3)=0.0d0
    if ( iflag == 1 ) then
       do i=1,m
          fvec(i) = (- zv/rad(i)*(x(1)*bberf(sqrt(x(2))*rad(i)) &
               + (1.0d0-x(1))*bberf(sqrt(x(3))*rad(i))) - vin(i))*sqrt(wgt(i))
       end do
    else if ( iflag == 2 ) then
       do i=1,m
          fjac(i,1) = -zv/rad(i) &
               *(bberf(sqrt(x(2)*rad(i)))-bberf(sqrt(x(3)*rad(i)))) &
               *sqrt(wgt(i))
          fjac(i,2) = -zv/sqrt(pi*x(2))*x(1) &
               *exp(-x(2)*rad(i)**2)*sqrt(wgt(i))
          fjac(i,3) = -zv/sqrt(pi*x(3))*(1.0d0-x(1)) &
               *exp(-x(3)*rad(i)**2)*sqrt(wgt(i))
       end do
    else
       write(6,*) 'Error in vlfit: iflag must be 1 or 2.'
       stop
    end if
    return
  END SUBROUTINE uscfit
!#endif

END MODULE simc_module
