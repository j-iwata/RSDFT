MODULE ps_local_mol_module

  use pseudopot_module, only: Rps,norb,Mr,NRps,rab,parloc,Zps,vql,rad,ippform
  use maskf_module
  use simc_module
  use atom_module, only: Natom,Nelement,aa_atom,ki_atom
  use ps_local_module, only: Vion
  use array_bound_module, only: ML_0,ML_1
  use rgrid_mol_module
  use ps_local_mol_gth_module
  use bberf_module
  use polint_module

  implicit none

  PRIVATE
  PUBLIC :: init_ps_local_mol,construct_ps_local_mol &
       ,vqls,Rcloc,NRcloc

  real(8),allocatable :: vqls(:,:)
  real(8),allocatable :: Rcloc(:)
  integer,allocatable :: NRcloc(:)

CONTAINS


  SUBROUTINE init_ps_local_mol(qcut)
    implicit none
    real(8),intent(IN) :: qcut
    integer :: i,ik,j,m,mm,m0,m1,m2,NRc,MMr,MKI
    real(8) :: pi,Rc,c1,c2,p1,p2,p3,p4,r,vlong,qc,eta
    real(8) :: const,fac,maxerr,x,dy0,y,c,y0,r1,dy,sum0
    real(8),allocatable :: vshort(:),wtmp(:),work(:),tmp(:),tmp1(:)

    if ( any( ippform == 4 ) ) then
       call init_ps_local_mol_gth
       return
    end if

    MKI = Nelement
    qc  = qcut
    MMr = maxval( Mr )
    pi  = acos(-1.d0)
    eta = 8.0d0

    allocate( NRcloc(MKI) ) ; NRcloc=0
    allocate( Rcloc(MKI)  ) ; Rcloc=0.d0
    allocate( vqls(MMr,MKI) ) ; vqls=0.d0

    allocate( vshort(MMr) ) ; vshort=0.d0
    allocate( wtmp(MMr)   ) ; wtmp=0.d0
    allocate( work(MMr)   ) ; work=0.d0
    allocate( tmp(MMr)    ) ; tmp=0.d0
    allocate( tmp1(MMr)   ) ; tmp1=0.d0

    do ik=1,MKI

       MMr = Mr(ik)

       if ( norb(ik) >= 1 ) then

          Rc = maxval( Rps(1:norb(ik),ik) )
          NRc= maxval( NRps(1:norb(ik),ik) )

       else

          Rc=5.d0
          do i=1,MMr
             if ( rad(i,ik) > Rc ) then
                NRc=i
                exit
             end if
          end do

       end if

       call simc(rad(1,ik),vql(1,ik),Rc,Zps(ik),parloc(1,ik),MMr)

       p1=parloc(1,ik) ; p2=sqrt(parloc(2,ik))
       p3=parloc(3,ik) ; p4=sqrt(parloc(4,ik))

       c1=-Zps(ik)*2.d0/sqrt(pi)
       c2=-Zps(ik)
       do i=1,MMr
          r=rad(i,ik)
          if ( r < 1.d-9 ) then
             vlong = c1*(p1*p2+p3*p4)
          else
             vlong = c2/r*( p1*bberf(p2*r)+p3*bberf(p4*r) )
          end if
          vshort(i) = vql(i,ik) - vlong
       end do

!- cut-off radius of the short-range part of the local potential

       NRcloc(ik)=0
       Rcloc(ik)=0.d0
       do i=1,MMr
          if ( abs(vshort(i)) < 1.d-8 ) then
             NRcloc(ik)=i
             Rcloc(ik)=rad(i,ik)
             exit
          end if
       end do

       if ( Rcloc(ik) == 0.0d0 ) then
          do i=1,MMr
             if ( abs(vshort(i)) < 1.d-4 ) then
                NRcloc(ik)=i
                Rcloc(ik)=rad(i,ik)
                exit
             end if
          end do
       end if

       if ( Rcloc(ik) == 0.0d0 ) then
          write(*,*) "log-range part fitting is failed"
          stop "init_local_mol(0)"
       end if

!- Filtering

       const = 2.d0/pi
       fac   = 1.2d0

       do i=1,MMr
          if ( rad(i,ik) >= fac*Rcloc(ik) ) then
             NRcloc(ik)=i
             exit
          end if
       end do
       Rcloc(ik) = rad(NRcloc(ik),ik)

       Rc  = Rcloc(ik)
       NRc = NRcloc(ik)

       call makemaskf(eta)

       wtmp(:)=0.d0
       work(:)=0.d0

       c=1.d0/Rc
       maxerr=0.d0
       do i=1,NRc
          x=rad(i,ik)*c
          if ( x<=dxm ) then
             y0=1.d0 ; dy0=0.d0
          else
             m0=int(x/dxm)
             dy0=1.d10
             do m=1,20
                m1=max(m0-m,1)       ; mm=m1-(m0-m)
                m2=min(m0+m+mm,nmsk) ; mm=(m0+m+mm)-m2
                call polint(xm(m1),maskr(m1),m2-m1+1,x,y,dy)
                if ( abs(dy)<dy0 ) then
                   y0=y ; dy0=abs(dy)
                end if
             end do
          end if
          wtmp(i)=y0
          maxerr=max(maxerr,dy0)
       end do

       do i=1,NRc
          r=rad(i,ik)
          if ( r==0.d0 ) then
             r1=rad(1,ik)
             if ( r1==0.d0 ) then
                tmp1(1)=qc*qc*qc/3.d0
             else
                tmp1(1)=sin(qc*r1)/(r1*r1*r1)-qc*cos(qc*r1)/(r1*r1)
             end if
             do j=2,NRc
                r1=rad(j,ik)
                tmp1(j)=sin(qc*r1)/(r1*r1*r1)-qc*cos(qc*r1)/(r1*r1)
             end do
          else
             do j=1,NRc
                r1=rad(j,ik)
                if ( r1==0.d0 ) then
                   tmp1(j)=sin(qc*r)/(r*r*r)-qc*cos(qc*r)/(r*r)
                else if ( r1==r ) then
                   tmp1(j)=(2*qc*r-sin(2.d0*qc*r))/(4*r*r*r)
                else
                   tmp1(j)=( sin(qc*(r-r1))/(r-r1)-sin(qc*(r+r1))/(r+r1) )/(2.d0*r*r1)
                end if
             end do
          end if
          do j=1,NRc
             tmp(j)=tmp1(j)*vshort(j)*rab(j,ik)*rad(j,ik)*rad(j,ik)/wtmp(j)
          end do
          call simp(tmp(1:NRc),sum0,2)
          vqls(i,ik)=sum0*const*wtmp(i)
       end do

    end do ! ik

    deallocate( tmp1 )
    deallocate( tmp  )
    deallocate( work )
    deallocate( wtmp )
    deallocate( vshort )

  END SUBROUTINE init_ps_local_mol

  SUBROUTINE simp(f,s,m)
    integer,intent(IN)  :: m
    real(8),intent(IN)  :: f(:)
    real(8),intent(OUT) :: s
    real(8),allocatable :: g(:)
    integer :: i,n,nn,nmax
    n=size(f) ; nmax=int(n/m)*m
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


  SUBROUTINE construct_ps_local_mol
    implicit none
    real(8),parameter :: ep=1.d-8
    real(8) :: p1,p2,p3,p4,const1
    real(8) :: Rx,Ry,Rz,r,x,y,z,Rc2,r2
    real(8) :: v,v0,maxerr,err0,err
    integer :: a,ik,i,M_irad,ir,ir0,m,m1,m2,mm,NRc
    integer,allocatable :: irad(:,:)

    if ( any( ippform == 4 ) ) then
       call construct_ps_local_mol_gth
       return
    end if

    if ( .not. allocated(Vion) ) then
       allocate( Vion(ML_0:ML_1) )
    end if
    Vion=0.d0

    const1 = 2.d0/sqrt(acos(-1.d0))

    allocate( irad(0:3000,Nelement) ) ; irad=0
    M_irad=0
    do ik=1,Nelement
       NRc=min( 3000, NRcloc(ik) )
       m=0
       irad(0,ik)=1
       do ir=1,NRc
          m=int(100.d0*rad(ir,ik))+1
          irad( m,ik )=ir
       end do
       ir=irad(0,ik)
       do i=1,m
          if ( irad(i,ik)==0 ) then
             irad(i,ik)=ir
             cycle
          end if
          ir=irad(i,ik)
       end do
       irad(m+1:,ik)=ir
       M_irad=max(M_irad,m)
    end do

    maxerr = 0.d0

    do a=1,Natom

       ik = ki_atom(a)
       Rx = aa_atom(1,a)
       Ry = aa_atom(2,a)
       Rz = aa_atom(3,a)

       p1 =-Zps(ik)*parloc(1,ik)
       p2 = sqrt(parloc(2,ik))
       p3 =-Zps(ik)*parloc(3,ik)
       p4 = sqrt(parloc(4,ik))

       do i=ML_0,ML_1
          x=LL(1,i)*Hsize-Rx
          y=LL(2,i)*Hsize-Ry
          z=LL(3,i)*Hsize-Rz
          r=sqrt(x*x+y*y+z*z)
          if ( r<1.d-9 ) then
             Vion(i)=Vion(i)+const1*(p1*p2+p3*p4)
          else
             Vion(i)=Vion(i)+( p1*bberf(p2*r)+p3*bberf(p4*r) )/r
          end if
       end do

       Rc2 = Rcloc(ik)*Rcloc(ik)
       NRc = NRcloc(ik)

       do i=ML_0,ML_1
          x=LL(1,i)*Hsize-Rx
          y=LL(2,i)*Hsize-Ry
          z=LL(3,i)*Hsize-Rz
          r2=x*x+y*y+z*z
          if ( r2 <= Rc2+1.d-10 ) then
             r=sqrt(r2)
             ir0=irad( int(100.d0*r),ik )
             do ir=ir0,NRc
                if ( r < rad(ir,ik) ) exit
             end do
             err0=1.d10
             v0=0.d0
             do mm=1,20
                m1=max(1,ir-mm)
                m2=min(ir+mm,NRc)
                call polint(rad(m1,ik),vqls(m1,ik),m2-m1+1,r,v,err)
                if ( abs(err)<err0 ) then
                   v0=v
                   err0=abs(err)
                   if ( err0 < ep ) exit
                end if
             end do
             maxerr=max(maxerr,err0)
             Vion(i)=Vion(i)+v0
          end if
       end do

    end do ! a

    deallocate( irad )

  END SUBROUTINE construct_ps_local_mol


END MODULE ps_local_mol_module
