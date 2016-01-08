MODULE xc_vdw_module

  use grid_module, only: grid
  use xc_variables, only: xcene, xcpot
  use ggrid_module, only: NMGL, MGL, GG, NGgrid, LLG, MG_0, MG_1 &
                         ,construct_ggrid, destruct_ggrid
  use gradient_module
  use parallel_module, only: comm_grid, myrank, nprocs, ir_grid, id_grid
  use fd_module, only: fd, construct_nabla_fd, destruct_nabla_fd
  use lattice_module, only: lattice,get_aa_lattice,get_reciprocal_lattice
  use basic_type_factory
  use fft_module

  implicit none

  PRIVATE
  PUBLIC :: init_xc_vdw, calc_xc_vdw

  integer,parameter :: DP=kind(0.0d0)
  integer,parameter :: QP=kind(0.0q0)

  real(8) :: Qmax=5.0d0
  real(8) :: Qmin=0.0d0
  integer :: NumQgrid=20
  real(8) :: SizeQgrid
  real(8),allocatable :: Qgrid(:)
  real(8) :: Dmax=48.6d0
  real(8),parameter :: Zab=-0.8491d0
  integer :: Mmax=12

  real(8),allocatable :: phiG(:,:,:)
  real(8),allocatable :: bmat(:,:),cmat(:,:),dmat(:,:)

  real(8),parameter :: zero_density = 1.d-10
  logical :: flag_init=.false.
  real(8) :: Ex_LDA, Ec_LDA, Ec_vdw

  include 'mpif.h'

CONTAINS


  SUBROUTINE init_xc_vdw
    implicit none
    integer :: i,j,n,np,info

    if ( flag_init ) return

    call write_border( 1, " init_xc_vdw(start)" )

    SizeQgrid = Qmax / NumQgrid
    allocate( Qgrid(0:NumQgrid) ) ; Qgrid=0.0d0
    do i=1,NumQgrid
       Qgrid(i) = Qmin + i*SizeQgrid
    end do

    allocate( phiG(NMGL,0:NumQgrid,0:NumQgrid) )
    phiG=0.0d0

    np=( NumQgrid*(NumQgrid+1) )/2
    n=0
    do j=1,NumQgrid
    do i=j,NumQgrid

       n = n + 1
       if ( np-1 < myrank .or. mod(n-1,nprocs) /= myrank ) cycle

       call calc_phiG( Qgrid(i), Qgrid(j), phiG(1,i,j) )

       if ( i /= j ) phiG(:,j,i) = phiG(:,i,j)

    end do ! i
    end do ! j

    call MPI_ALLREDUCE( MPI_IN_PLACE, phiG, size(phiG), MPI_REAL8 &
         ,MPI_SUM, MPI_COMM_WORLD, info )

    call init1_xc_vdw ! coefficients of 3rd spline

    flag_init = .true.

    call write_border( 1, " init_xc_vdw(end)")

  END SUBROUTINE init_xc_vdw


  SUBROUTINE calc_phiG( qi, qj, phiG )
    implicit none
    real(8),intent(IN)  :: qi,qj
    real(8),intent(OUT) :: phiG(NMGL)
    integer :: nr,ig,i
    real(8) :: rmax,dr,r,pi4,G,sum0,ct0,ct1,sum1
    real(8),allocatable :: phi(:), phiG0(:)

    phiG = 0.0d0

    rmax = Dmax/max( qi, qj )
    nr   = 400
    dr   = rmax/nr

    allocate( phi(0:nr)   ) ; phi=0.0d0
    allocate( phiG0(NMGL) ) ; phiG0=0.0d0

    do i=1,nr

       r = i*dr
       call calc_phi_b( qi*r, qj*r, phi(i) )

       call check_phiG( i, dr, phi(1), phiG )
       sum0=sum((phiG-phiG0)**2)
       if ( sum0 <= 1.d-10 ) exit
       phiG0(:)=phiG(:)

    end do ! i

    deallocate( phiG0 )
    deallocate( phi )

  END SUBROUTINE calc_phiG

  SUBROUTINE check_phiG( nr, dr, phi, phiG )
    implicit none
    integer,intent(IN) :: nr
    real(8),intent(IN) :: dr, phi(nr)
    real(8),intent(OUT) :: phiG(:)
    integer :: i,ig
    real(8) :: r,G,sum0,pi4

    pi4 = 4.0d0*acos(-1.0d0)

    do ig=1,NMGL

       G=sqrt( GG(ig) )

       if ( G <= 1.d-9 ) then

          sum0=0.0d0
          do i=1,nr
             r = i*dr
             sum0 = sum0 + r*r*phi(i)
          end do
          phiG(ig)=pi4*sum0*dr

       else

          sum0=0.0d0
          do i=1,nr
             r = i*dr
             sum0 = sum0 + r*phi(i)*sin(G*r)
          end do
          phiG(ig)=pi4/G*sum0*dr

       end if

    end do ! ig

  END SUBROUTINE check_phiG

  SUBROUTINE calc_phi_b( di, dj, phi )
    implicit none
    real(8),intent(IN)  :: di,dj
    real(8),intent(OUT) :: phi
    integer :: nt,nt0,it,loop
    real(8) :: sum1,phi0,t,dt,dt0,pi
    real(8),parameter :: tol=1.d-5

    pi = acos(-1.0d0)

    nt0 = 8
    dt0 = 0.5d0*pi / nt0

    phi  = 0.0d0
    phi0 = 0.0d0

    nt = nt0/2 ! Interval [0,pi/4] is enough due to symmetry 
    dt = dt0

    call Int_RadialPart( 0.25d0*pi, di, dj, sum1 )

    phi = 0.5d0*( 0.0d0 + sum1 )

    do it=1,nt-1
       t = it*dt
       call Int_RadialPart( t, di, dj, sum1 )
       phi = phi + sum1
    end do ! it

    nt = nt/2

    do loop=1,10

       nt = nt*2
       dt = dt*0.5d0

       do it=1,nt
          t = (2*it-1)*dt
          call Int_RadialPart( t, di, dj, sum1 )
          phi = phi + sum1
       end do ! it

       if ( abs(phi*dt-phi0) <= tol ) exit

       phi0 = phi*dt

    end do ! loop

    phi = 2.0d0*(2.0d0/pi**2)*phi*dt

  END SUBROUTINE calc_phi_b

  SUBROUTINE Int_RadialPart( t, di, dj, sum_rad )
    implicit none
    real(8),intent(IN)  :: t,di,dj
    real(8),intent(OUT) :: sum_rad
    integer :: nx0,nx,ix,loop
    real(8) :: sum0,sum1,const
    real(8) :: xmin,xmax,dx,dx0,x,coshx,sinhx,r
    real(8),parameter :: tol=1.d-8

    xmax = 3.0d0
    xmin =-3.0d0
    nx0  = 50
    dx0  = (xmax-xmin)/nx0

    sum0 = 0.0d0
    sum1 = 0.0d0

    const = 0.5d0*acos(-1.0d0)

    nx = nx0
    dx = dx0
    do ix=1,nx-1
       x=xmin+ix*dx
       coshx = 0.5d0*( exp(x) + exp(-x) )
       sinhx = 0.5d0*( exp(x) - exp(-x) )
       r = exp(const*sinhx)
       sum0 = sum0 + phi_b(r,t,di,dj)*coshx*const*r
!       sum0 = sum0 + phi_c(r,t,di,dj)*coshx*const*r
    end do ! ix
    sum1=sum0*dx

    nx = nx/2

    do loop=1,20

       nx = nx*2
       dx = dx*0.5d0
       do ix=1,nx
          x=xmin+(2*ix-1)*dx
          coshx = 0.5d0*( exp(x) + exp(-x) )
          sinhx = 0.5d0*( exp(x) - exp(-x) )
          r = exp(const*sinhx)
          sum0 = sum0 + phi_b(r,t,di,dj)*coshx*const*r
!          sum0 = sum0 + phi_c(r,t,di,dj)*coshx*const*r
       end do ! ix

       if ( abs(sum0*dx-sum1) <= tol ) exit
       sum1=sum0*dx

    end do ! loop

    sum_rad = sum0*dx

  END SUBROUTINE Int_RadialPart

  FUNCTION phi_b( r, theta, di, dj )
    implicit none
    real(8) :: phi_b
    real(8),intent(IN)  :: r, theta, di, dj
    real(8) :: a,b,WWW,TTT !,gamma,pi
    real(8) :: ct,st,sa,sb,ca,cb,hai,haj,hbi,hbj,vai,vaj,vbi,vbj
    real(8) :: a2,b2,cai,cbi,caj,cbj,gdidi_i,gdjdj_i
    real(8),parameter :: pi=3.141592653589793d0
    real(8),parameter :: gamma=1.396263401595464d0
    real(8),parameter :: thr=3.0d0

!    pi = acos(-1.0d0)
!    gamma = 4.0d0*pi/9.0d0

    a = r*cos(theta)
    b = r*sin(theta)

    a2=a*a
    b2=b*b

    gdidi_i = gamma/(di*di)
    gdjdj_i = gamma/(dj*dj)

    cai = a2*gdidi_i
    cbi = b2*gdidi_i
    caj = a2*gdjdj_i
    cbj = b2*gdjdj_i
!    cai = gamma*a2/(di*di)
!    cbi = gamma*b2/(di*di)
!    caj = gamma*a2/(dj*dj)
!    cbj = gamma*b2/(dj*dj)

!    ct = cos(theta)
!    st = sin(theta)

!    a = r*ct
!    b = r*st

!    ca = cos(a)
!    cb = cos(b)
!    sa = 1.0d0 ; if ( a /= 0.0d0 ) sa=sin(a)/a
!    sb = 1.0d0 ; if ( b /= 0.0d0 ) sb=sin(b)/b
!    sa=sin(a)/a
!    sb=sin(b)/b

!    WWW = ( (thr-a2)*sa-thr*ca )*cos(b) &
!        + ( (thr-b2)*ca + (r*r-thr)*sa )*sin(b)/b
!    WWW = ( (thr-a2)*sa-thr*ca )*cb + ( (thr-b2)*ca + (r*r-thr)*sa )*sb
!    WWW = (thr-a2)*cb*sa + (thr-b2)*ca*sb + (r*r-thr)*sa*sb - thr*ca*cb
!    WWW = (3.0d0-a2)*cb*sa + (3.0d0-b2)*ca*sb &
!        + (r*r-3.0d0)*sa*sb - 3.0d0*ca*cb
!    WWW = 2.0d0*( (3.0d0-a2)*cb*sa + (3.0d0-b2)*ca*sb &
!                + (r*r-3.0d0)*sa*sb - 3.0d0*ca*cb )
!    WWW = 2.0d0*( (3.0d0-a*a)*cb*sa + (3.0d0-b*b)*ca*sb &
!                + (r*r-3.0d0)*sa*sb - 3.0d0*ca*cb )

!    if ( di == 0.0d0 ) then
!       hai = 0.0d0
!       hbi = 0.0d0
!    else
!       hai = 1.0d0 - exp( -gamma*(a/di)**2 )
!       hbi = 1.0d0 - exp( -gamma*(b/di)**2 )
!       hai = 1.0d0 - exp( -cai )
!       hbi = 1.0d0 - exp( -cbi )
!    end if
!    if ( dj == 0.0d0 ) then
!       haj = 0.0d0
!       hbj = 0.0d0
!    else
!       haj = 1.0d0 - exp( -gamma*(a/dj)**2 )
!       hbj = 1.0d0 - exp( -gamma*(b/dj)**2 )
!       haj = 1.0d0 - exp( -caj )
!       hbj = 1.0d0 - exp( -cbj )
!    end if

!    if ( a == 0.0d0 ) then
!       vai = di*di/(2.0d0*gamma)
!       vaj = dj*dj/(2.0d0*gamma)
!    else
!       vai = a2/(2.0d0*hai)
!       vaj = a2/(2.0d0*haj)
!       vai = a2/( 2.0d0-2.0d0*exp(-cai) )
!       vaj = a2/( 2.0d0-2.0d0*exp(-caj) )
       vai = 0.5d0*a2/( 1.0d0-exp(-cai) )
       vaj = 0.5d0*a2/( 1.0d0-exp(-caj) )
!       vai = a2/( 1.0d0-exp(-cai) )
!       vaj = a2/( 1.0d0-exp(-caj) )
!    end if
!    if ( b == 0.0d0 ) then
!       vbi = di*di/(2.0d0*gamma)
!       vbj = dj*dj/(2.0d0*gamma)
!    else
!       vbi = b2/(2.0d0*hbi)
!       vbj = b2/(2.0d0*hbj)
!       vbi = b2/( 2.0d0-2.0d0*exp(-cbi) )
!       vbj = b2/( 2.0d0-2.0d0*exp(-cbj) )
       vbi = 0.5d0*b2/( 1.0d0-exp(-cbi) )
       vbj = 0.5d0*b2/( 1.0d0-exp(-cbj) )
!       vbi = b2/( 1.0d0-exp(-cbi) )
!       vbj = b2/( 1.0d0-exp(-cbj) )
!    end if

    TTT = ( 1.0d0/(vai+vbi) + 1.0d0/(vaj+vbj) ) &
               *( 1.0d0/( (vai+vaj)*(vbi+vbj) ) &
                 +1.0d0/( (vai+vbj)*(vaj+vbi) ) )
!    TTT = ( 1.0d0/(vai+vbi) + 1.0d0/(vaj+vbj) ) &
!               *( 1.0d0/( (vai+vaj)*(vbi+vbj) ) &
!                 +1.0d0/( (vai+vbj)*(vaj+vbi) ) )
!    TTT = 0.5d0*( 1.0d0/(vai+vbi) + 1.0d0/(vaj+vbj) ) &
!               *( 1.0d0/( (vai+vaj)*(vbi+vbj) ) &
!                 +1.0d0/( (vai+vbj)*(vaj+vbi) ) )

    ca = cos(a)
    sa = sin(a)/a
    WWW = ( (thr-a2)*sa-thr*ca )*cos(b) &
        + ( (thr-b2)*ca + (r*r-thr)*sa )*sin(b)/b

    phi_b = r*WWW*TTT

  END FUNCTION phi_b

  FUNCTION phi_c( r, theta, di, dj )
    implicit none
    real(8) :: phi_c
    real(8),intent(IN)  :: r, theta, di, dj
    real(8) :: a,b,gamma,pi,WWW,TTT
    real(8) :: ct,st,sa,sb,ca,cb,hai,haj,hbi,hbj,vai,vaj,vbi,vbj

    pi = acos(-1.0d0)
    gamma = 4.0d0*pi/9.0d0

    ct = cos(theta)
    st = sin(theta)

    a = r*ct
    b = r*st

    ca = cos(a)
    cb = cos(b)
    sa = 1.0d0 ; if ( a /= 0.0d0 ) sa=sin(a)/a
    sb = 1.0d0 ; if ( b /= 0.0d0 ) sb=sin(b)/b

    WWW = 2.0d0*( (3.0d0-a*a)*cb*sa + (3.0d0-b*b)*ca*sb &
                + (r*r-3.0d0)*sa*sb - 3.0d0*ca*cb )

    if ( di == 0.0d0 ) then
       hai = 0.0d0
       hbi = 0.0d0
    else
       hai = 1.0d0 - exp( -gamma*(a/di)**2 )
       hbi = 1.0d0 - exp( -gamma*(b/di)**2 )
    end if
    if ( dj == 0.0d0 ) then
       haj = 0.0d0
       hbj = 0.0d0
    else
       haj = 1.0d0 - exp( -gamma*(a/dj)**2 )
       hbj = 1.0d0 - exp( -gamma*(b/dj)**2 )
    end if

    if ( a == 0.0d0 ) then
       vai = di*di/(2.0d0*gamma)
       vaj = dj*dj/(2.0d0*gamma)
    else
       vai = a*a/(2.0d0*hai)
       vaj = a*a/(2.0d0*haj)
    end if
    if ( b == 0.0d0 ) then
       vbi = di*di/(2.0d0*gamma)
       vbj = dj*dj/(2.0d0*gamma)
    else
       vbi = b*b/(2.0d0*hbi)
       vbj = b*b/(2.0d0*hbj)
    end if

    TTT = 0.5d0*( 1.0d0/(vai+vbi) + 1.0d0/(vaj+vbj) ) &
               *( 1.0d0/( (vai+vaj)*(vbi+vbj) ) &
                 +1.0d0/( (vai+vbj)*(vaj+vbi) ) )

    phi_c = r*WWW*TTT

  END FUNCTION phi_c


  SUBROUTINE init1_xc_vdw
    implicit none
    real(8),allocatable :: mat1(:,:),mat2(:,:),work(:)
    integer,allocatable :: ipiv(:)
    real(8) :: h
    integer :: m,n,i,info

    n = NumQgrid
    h = SizeQgrid
    m = n+1

    allocate( bmat(0:n,0:n) ) ; bmat=0.0d0
    allocate( cmat(0:n,0:n) ) ; cmat=0.0d0
    allocate( dmat(0:n,0:n) ) ; dmat=0.0d0

    allocate( mat1(0:n,0:n) ) ; mat1=0.0d0
    allocate( mat2(0:n,0:n) ) ; mat2=0.0d0
    allocate( ipiv(0:n)     ) ; ipiv=0
    allocate( work(64*m)    ) ; work=0.0d0

    mat1(0,0)=1.0d0
    do i=1,n-1
       mat1(i,i-1)=h
       mat1(i,i  )=4.0d0*h
       mat1(i,i+1)=h
    end do
    mat1(n,n)=1.0d0

    do i=1,n-1
       mat2(i,i-1)= 1.0d0
       mat2(i,i  )=-2.0d0
       mat2(i,i+1)= 1.0d0
    end do
    mat2(0:n,0:n)=mat2(0:n,0:n)*3.0d0/h

    call dgetrf( m,m,mat1,m,ipiv,info)
    call dgetri( m,mat1,m,ipiv,work,size(work),info)

    cmat(0:n,0:n) = matmul( mat1(0:n,0:n), mat2(0:n,0:n) )

    mat1=0.0d0
    do i=0,n-1
       mat1(i,i  )=-1.0d0
       mat1(i,i+1)= 1.0d0
    end do

    dmat(0:n,0:n) = matmul( mat1(0:n,0:n), cmat(0:n,0:n) )/(3.0d0*h)

    bmat(0:n,0:n) = mat1(0:n,0:n)/h - cmat(0:n,0:n)*h - dmat(0:n,0:n)*h**2

    deallocate( work )
    deallocate( ipiv )
    deallocate( mat2 )
    deallocate( mat1 )

  END SUBROUTINE init1_xc_vdw


  SUBROUTINE calc_xc_vdw( rgrid, rho, ene, pot )
    implicit none
    type( grid ),intent(IN) :: rgrid
    type( GSArray ),intent(IN) :: rho
    type( xcene ) :: ene
    type( xcpot ) :: pot
    type( gradient16 ) :: grad16
    logical :: disp_sw
    integer :: i,j,m0,m1,mm,nspin
    real(8) :: rho_tmp(2), trho, edx_lda(2), edc_lda,ss_min,ss_max
    real(QP) :: kf,pi,onethr,fouthr,ss,const,c,q0
    real(8) :: q0_min,q0_max,sbuf(2),rbuf(2)
    real(8) :: vdx_lda(2), vdc_lda(2)
    real(8),allocatable :: qsat(:), pol_spline(:,:)
    real(8),allocatable :: ua(:,:), pol_drv(:,:), f(:), g(:)
    complex(8),allocatable :: theta(:,:)

    call write_border( 1, " calc_xc_vdw(start)" )

    pi     = acos(-1.0q0)
    onethr = 1.0q0/3.0q0
    fouthr = 4.0q0/3.0q0
    const  = Zab/9.0q0
    m0     = rgrid%g1%head
    m1     = rgrid%g1%tail
    nspin  = rho%s_range%size_global

    call construct_gradient16( rgrid, rho, grad16 )

    allocate( qsat(m0:m1) ) ; qsat=0.0d0

    Ex_LDA = 0.0d0
    Ec_LDA = 0.0d0

    ss_min = 1.d100
    ss_max =-1.d100

    pot%c%val(:,:)=0.0d0

    do i = pot%xc%g_range%head, pot%xc%g_range%tail

       rho_tmp(1) = rho%val(i,1)
       rho_tmp(2) = rho%val(i,nspin)
       trho       = sum( rho_tmp(1:nspin) )

       if ( trho < zero_density ) cycle

       call calc_edx_lda( nspin, rho_tmp, edx_lda, vdx_lda )
       call calc_edc_lda( nspin, rho_tmp, edc_lda, vdc_lda )

       Ex_LDA = Ex_LDA + sum( rho_tmp(1:nspin)*edx_lda(1:nspin) )
       Ec_LDA = Ec_LDA + trho*edc_lda

       do j = pot%xc%s_range%head, pot%xc%s_range%tail
          pot%c%val(i,j) = edc_lda + trho*vdc_lda(j)
       end do

       kf = ( 3.0q0*pi*pi*trho )**onethr

       if ( trho <= 0.0d0 ) then
          qsat(i)=Qmax
       else
          ss = grad16%gg(i)/(2.0q0*kf*trho)**2
          q0 = -fouthr*pi*(edx_lda(1)+edc_lda) - const*ss*kf
          c=0.0q0
          do j=1,Mmax
             c=c+(q0/Qmax)**j/j
          end do
          qsat(i) = Qmax*( 1.0q0 - exp(-c) )
          ss_min = min( ss, ss_min )
          ss_max = max( ss, ss_max )
       end if

    end do ! i

    sbuf(1) = minval(qsat)
    sbuf(2) = maxval(qsat)
    call MPI_ALLREDUCE(sbuf(1),q0_min,1,MPI_REAL8,MPI_MIN,comm_grid,i)
    call MPI_ALLREDUCE(sbuf(2),q0_max,1,MPI_REAL8,MPI_MAX,comm_grid,i)

    sbuf(1) = ss_min
    sbuf(2) = ss_max
    call MPI_ALLREDUCE(sbuf(1),ss_min,1,MPI_REAL8,MPI_MIN,comm_grid,i)
    call MPI_ALLREDUCE(sbuf(2),ss_max,1,MPI_REAL8,MPI_MAX,comm_grid,i)

    sbuf(1) = Ex_LDA * rgrid%VolumeElement
    sbuf(2) = Ec_LDA * rgrid%VolumeElement
    call MPI_ALLREDUCE(sbuf,rbuf,2,MPI_REAL8,MPI_SUM,comm_grid,i)
    Ex_LDA = rbuf(1)
    Ec_LDA = rbuf(2)

!    call check_disp_switch( disp_sw, 0 )
!    if ( disp_sw ) then
!       write(*,*) "ss(min,max)=",ss_min,ss_max
!       write(*,*) "q0(min,max)=",q0_min,q0_max
!       write(*,*) "Ex_LDA =",Ex_LDA
!       write(*,*) "Ec_LDA =",Ec_LDA
!       write(*,*) "Exc_LDA=",Ex_LDA+Ec_LDA
!    end if

! ---

    allocate( pol_spline(m0:m1,0:NumQgrid) ) ; pol_spline=0.0d0

    call calc_pol_spline( m0, m1, qsat, pol_spline )

! ---

    mm = max( m1-m0+1, MG_1-MG_0+1 )

    allocate( theta(mm,0:NumQgrid) ) ; theta=(0.0d0,0.0d0)

    do j=0,NumQgrid
       do i=m0,m1
          theta(i-m0+1,j) = sum( rho%val(i,:) )*pol_spline(i,j)
       end do
    end do

    call fft_theta( rgrid, theta )

! ---

    call calc_vdw_energy( rgrid, theta, Ec_vdw )

    ene%Ec = Ec_LDA + Ec_VDW

! ---

    allocate( ua(m0:m1,0:NumQgrid) ) ; ua=0.0d0

    call calc_ua( rgrid, theta, ua )

! ---

    deallocate( theta )

! ---

    allocate( pol_drv(m0:m1,0:NumQgrid) ) ; pol_drv=0.0d0

    call calc_pol_drv( m0, m1, qsat, pol_spline, pol_drv )

! ---

    allocate( f(m0:m1) ) ; f=0.0d0
    allocate( g(m0:m1) ) ; g=0.0d0

    do j=0,NumQgrid
       f(:) = f(:) + ua(:,j)*pol_drv(:,j)
    end do

    call calc_dqdn( rgrid, rho, grad16, f, g )

! ---

    do j=0,NumQgrid
       g(:) = g(:) + ua(:,j)*pol_spline(:,j)
    end do

    c = 1.0d0/rgrid%g1%size_global

    pot%c%val(:,1) = pot%c%val(:,1) + c*g(:)

! ---

    deallocate( g )
    deallocate( f )
    deallocate( pol_drv )
    deallocate( ua )
    deallocate( pol_spline )
    deallocate( qsat )

    call destruct_gradient16( grad16 )

    call write_border( 1, " calc_xc_vdw(end)" )

  END SUBROUTINE calc_xc_vdw


  SUBROUTINE calc_edx_lda( nn, rho, edx_lda, vdx_lda )
    implicit none
    integer,intent(IN)  :: nn
    real(8),intent(IN)  :: rho(nn)
    real(8),optional,intent(OUT) :: edx_lda(nn)
    real(8),optional,intent(OUT) :: vdx_lda(nn)
    real(8) :: onethr, twothr, thrfou, thrpi, pi, factor
    integer :: s

    pi     = acos(-1.0d0)
    onethr = 1.0d0/3.0d0
    twothr = 2.0d0/3.0d0
    thrfou = 3.0d0/4.0d0
    thrpi  = 3.0d0/pi

    factor = 1.0d0
    if ( nn == 2 ) factor = 2.0d0**onethr

    if ( present(edx_lda) ) then
       do s=1,nn
          edx_lda(s) = -factor*thrfou*( thrPi*rho(s) )**onethr
       end do
    end if

    if ( present(vdx_lda) ) then
       do s=1,nn
          if ( rho(s) > 0.0d0 ) then
             vdx_lda(s) = -pi*thrfou/(3.0d0*pi*pi*rho(s))**twothr
          else
             vdx_lda(s)=0.0d0
          end if
       end do
    end if

  END SUBROUTINE calc_edx_lda


  SUBROUTINE calc_edc_lda( nn, rho, edc_lda, vdc_lda )
    implicit none
    integer,intent(IN)  :: nn
    real(8),intent(IN)  :: rho(nn)
    real(8),optional,intent(OUT) :: edc_lda
    real(8),optional,intent(OUT) :: vdc_lda(nn)
    real(8),parameter :: A00  =0.031091d0,A01  =0.015545d0,A02  =0.016887d0
    real(8),parameter :: alp10=0.21370d0 ,alp11=0.20548d0 ,alp12=0.11125d0
    real(8),parameter :: bt10 =7.5957d0  ,bt11 =1.41189d1 ,bt12 =1.0357d1
    real(8),parameter :: bt20 =3.5876d0  ,bt21 =6.1977d0  ,bt22 =3.6231d0
    real(8),parameter :: bt30 =1.6382d0  ,bt31 =3.3662d0  ,bt32 =0.88026d0
    real(8),parameter :: bt40 =0.49294d0 ,bt41 =0.62517d0 ,bt42 =0.49671d0
    real(8) :: one,two,fouthr,pi,const1,const2,onethr,ThrFouPi
    real(8) :: trho,rhoa,rhob,zeta,fz,kf,rs,rssq,ec_U,ec_P,alpc
    real(8) :: deU_dn, deP_dn, dac_dn, dfz_dz, dz_dn(2), dec_dz
    integer :: s

    trho = sum( rho(1:nn) )
    rhoa = rho(1)
    rhob = rho(nn)

    if ( present(edc_lda) ) edc_lda=0.0d0
    if ( present(vdc_lda) ) vdc_lda=0.0d0

    if ( trho <= zero_density ) return

    one      = 1.0d0
    two      = 2.0d0
    fouthr   = 4.0d0/3.0d0
    onethr   = 1.0d0/3.0d0
    pi       = acos(-1.0d0)
    ThrFouPi = 3.0d0/(4.0d0*pi)
    const1   = two**fouthr-two
    const2   = 9.0d0/4.0d0*(two**onethr-one)

    zeta = ( rhoa - rhob )/trho
    fz   = ( (one+zeta)**fouthr + (one-zeta)**fouthr - two )*const1
    kf   = ( 3.0d0*Pi*Pi*trho )**onethr
    rs   = ( ThrFouPi/trho )**onethr
    rssq = sqrt(rs)
    ec_U = -two*A00*( one + alp10*rs ) &
         *log( one + one/( two*A00 &
         *(bt10*rssq + bt20*rs + bt30*rs*rssq + bt40*rs*rs) ) )
    ec_P = -two*A01*( one + alp11*rs ) &
         *log( one + one/( two*A01 &
         *(bt11*rssq + bt21*rs + bt31*rs*rssq + bt41*rs*rs) ) )
    alpc = -two*A02*( one + alp12*rs ) &
         *log( one + one/( two*A02 &
         *(bt12*rssq + bt22*rs + bt32*rs*rssq + bt42*rs*rs) ) )

    if ( present(edc_lda) ) then
       edc_lda = ec_U - alpc*fz*const2*(one-zeta**4) &
            + ( ec_P - ec_U )*fz*zeta**4
    end if

    if ( present(vdc_lda) ) then

       deU_dn = -4.0d0*Pi/9.0d0*rs**4*alp10*ec_U/( 1.0d0 + alp10*rs ) &
            -4.0d0*Pi/9.0d0*rs*rs*( 1.0d0 + alp10*rs ) &
            *( 0.5d0*bt10*sqrt(rs) + bt20*rs + 1.5d0*bt30*rs*sqrt(rs) &
            +2.0d0*bt40*rs*rs ) &
            /( bt10 + bt20*sqrt(rs) + bt30*rs &
            + bt40*rs*sqrt(rs))**2 &
            * exp( ec_U/(2.0d0*A00*(1.0d0+alp10*rs)) )

       deP_dn = -4.0d0*Pi/9.0d0*rs**4*alp11*ec_P/( 1.0d0 + alp11*rs ) &
            -4.0d0*Pi/9.0d0*rs*rs*( 1.0d0 + alp11*rs ) &
            *( 0.5d0*bt11*sqrt(rs) + bt21*rs + 1.5d0*bt31*rs*sqrt(rs) &
            +2.0d0*bt41*rs*rs ) &
            /( bt11 + bt21*sqrt(rs) + bt31*rs &
            + bt41*rs*sqrt(rs))**2 &
            * exp( ec_P/(2.0d0*A01*(1.0d0+alp11*rs)) )

       dac_dn = -4.0d0*Pi/9.0d0*rs**4*alp12*alpc/( 1.0d0 + alp12*rs ) &
            -4.0d0*Pi/9.0d0*rs*rs*( 1.0d0 + alp12*rs ) &
            *( 0.5d0*bt12*sqrt(rs) + bt22*rs + 1.5d0*bt32*rs*sqrt(rs) &
            +2.0d0*bt42*rs*rs ) &
            /( bt12 + bt22*sqrt(rs) + bt32*rs &
            + bt42*rs*sqrt(rs))**2 &
            * exp( alpc/(2.0d0*A02*(1.d0+alp12*rs)) )

       dfz_dz = fouthr*( (one+zeta)**onethr - (one-zeta)**onethr )*const1

       dz_dn(1)  = 2.0d0*rhob/trho
       dz_dn(nn) =-2.0d0*rhoa/trho

       dec_dz = -alpc*dfz_dz*const2*(one-zeta**4) &
               + 4.0d0*alpc*fz*const2*zeta**3 &
               +(ec_P-ec_U)*dfz_dz*zeta**4 &
               +(ec_P-ec_U)*fz*4.0d0*zeta**3

       do s=1,nn
          vdc_lda(s) = deU_dn - dac_dn*fz*const2*(one-zeta**4) &
               +(deP_dn-deU_dn)*fz*zeta**4 + dec_dz*dz_dn(s)
       end do

    end if

  END SUBROUTINE calc_edc_lda


  SUBROUTINE calc_pol_spline( m0, m1, q0, pol )
    implicit none
    integer,intent(IN)  :: m0,m1
    real(8),intent(IN)  :: q0(m0:m1)
    real(8),intent(OUT) :: pol(m0:m1,0:NumQgrid)
    integer :: i,j,k
    real(8) :: q

    pol(:,:)=0.0d0

    do k=m0,m1

       do i=0,NumQgrid
          if ( q0(k) < Qgrid(i) ) exit
       end do
       i=i-1

       q=q0(k)-Qgrid(i)
       do j=0,NumQgrid
          pol(k,j) = bmat(i,j)*q + cmat(i,j)*q**2 + dmat(i,j)*q**3
       end do
       pol(k,i) = pol(k,i) + 1.0d0

! ---------------------------- Wu & Gygi ---
       if ( q0(k) /= 0.0d0 ) then
          do j=0,NumQgrid
             pol(k,j) = Qgrid(j)*pol(k,j)/q0(k)
          end do
       else
          pol(k,:) = 0.0d0
       end if
! ------------------------------------------

    end do ! k

  END SUBROUTINE calc_pol_spline


  SUBROUTINE fft_theta( rgrid, theta )
    implicit none
    type( grid ),intent(IN) :: rgrid
    complex(8),intent(INOUT) :: theta(:,0:)
    integer :: ML,ML1,ML2,ML3,mm,i1,i2,i3,i,j,info
    complex(8),allocatable :: work0(:,:,:), work1(:,:,:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)

    ML  = rgrid%g1%size_global
    ML1 = rgrid%g3%x%size_global
    ML2 = rgrid%g3%y%size_global
    ML3 = rgrid%g3%z%size_global

    allocate( work0(0:ML1-1,0:ML2-1,0:ML3-1) ) ; work0=zero
    allocate( work1(0:ML1-1,0:ML2-1,0:ML3-1) ) ; work1=zero

    call init_fft

    call construct_ggrid(1)

    do j=0,NumQgrid

       work1(:,:,:)=zero
       i=0
       do i3=rgrid%g3%z%head,rgrid%g3%z%tail
       do i2=rgrid%g3%y%head,rgrid%g3%y%tail
       do i1=rgrid%g3%x%head,rgrid%g3%x%tail
          i=i+1
          work1(i1,i2,i3) = theta(i,j)
       end do ! i1
       end do ! i2
       end do ! i3

       call MPI_ALLREDUCE( work1, work0, size(work0) &
            , MPI_COMPLEX16, MPI_SUM, comm_grid, info )

       call forward_fft( work0, work1 )

       theta(:,j) = zero
       do i=MG_0,MG_1
          i1=mod( LLG(1,i)+ML1, ML1 )
          i2=mod( LLG(2,i)+ML2, ML2 )
          i3=mod( LLG(3,i)+ML3, ML3 )
          theta(i-MG_0+1,j) = work0(i1,i2,i3)
       end do

    end do ! j

    call destruct_ggrid

    call finalize_fft

    deallocate( work1 )
    deallocate( work0 )

  END SUBROUTINE fft_theta


  SUBROUTINE calc_vdw_energy( rgrid, theta, Ec )
    implicit none
    complex(8),intent(IN) :: theta(:,0:)
    type( grid ),intent(IN) :: rgrid
    real(8),intent(INOUT) :: Ec
    integer :: a,b,i,j,info
    complex(8) :: z
    real(8) :: c,E0
    logical :: disp_sw

    z=(0.0d0,0.0d0)

    do b=0,NumQgrid
    do a=b,NumQgrid

       c=1.0d0 ; if ( a /= b ) c=2.0d0

       do i=MG_0,MG_1
          j=i-MG_0+1
          z = z + c*theta(j,a)*conjg(theta(j,b))*phiG(MGL(i),a,b)
       end do ! i

    end do ! a
    end do ! b

    E0 = 0.5d0*z * rgrid%g1%size_global * rgrid%VolumeElement
    call MPI_ALLREDUCE(E0,Ec,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,info)

!    call check_disp_switch( disp_sw, 0 )
!    if ( disp_sw ) then
!       write(*,*) "Ec_vdw=",Ec
!       write(*,*) "Ec_vdw+Ec_LDA=",Ec+Ec_LDA
!    end if

  END SUBROUTINE calc_vdw_energy


  SUBROUTINE calc_ua( rgrid, theta, ua )
    implicit none
    complex(8),intent(IN) :: theta(:,0:)
    type( grid ),intent(IN) :: rgrid
    real(8),intent(OUT)   :: ua(:,0:)
    integer :: a,b,i,j,i1,i2,i3
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8),allocatable :: work0(:,:,:),work1(:,:,:)
    integer :: ML1,ML2,ML3,ML,info

    ML  = rgrid%g1%size_global
    ML1 = rgrid%g3%x%size_global
    ML2 = rgrid%g3%y%size_global
    ML3 = rgrid%g3%z%size_global

    allocate( work0(0:ML1-1,0:ML2-1,0:ML3-1) ) ; work0=zero
    allocate( work1(0:ML1-1,0:ML2-1,0:ML3-1) ) ; work1=zero

    call init_fft

    call construct_ggrid(1)

    ua(:,:) = 0.0d0

    do a=1,NumQgrid

       work1(:,:,:)=zero

       do b=1,NumQgrid

          do i=MG_0,MG_1
             j=i-MG_0+1
             i1 = mod( ML1+LLG(1,i), ML1 )
             i2 = mod( ML2+LLG(2,i), ML2 )
             i3 = mod( ML3+LLG(3,i), ML3 )
             work1(i1,i2,i3) = work1(i1,i2,i3) + phiG(MGL(i),a,b)*theta(j,b)
          end do ! i

       end do ! b

       call MPI_ALLREDUCE(work1,work0,size(work0),MPI_COMPLEX16 &
            ,MPI_SUM,MPI_COMM_WORLD,info)

       call backward_fft( work0, work1 )

       i=0
       do i3=rgrid%g3%z%head,rgrid%g3%z%tail
       do i2=rgrid%g3%y%head,rgrid%g3%y%tail
       do i1=rgrid%g3%x%head,rgrid%g3%x%tail
          i=i+1
          ua(i,a) = work0(i1,i2,i3)
       end do
       end do
       end do

    end do ! a

    call destruct_ggrid

    call finalize_fft

    deallocate( work1 )
    deallocate( work0 )

  END SUBROUTINE calc_ua


  SUBROUTINE calc_pol_drv( m0, m1, q0, pol, pol_drv )
    implicit none
    integer,intent(IN)  :: m0,m1
    real(8),intent(IN)  :: q0(m0:m1),pol(m0:m1,0:NumQgrid)
    real(8),intent(OUT) :: pol_drv(m0:m1,0:NumQgrid)
    integer :: i,j,k
    real(8) :: q

    pol_drv(:,:)=0.0d0

    do k=m0,m1

       do i=0,NumQgrid
          if ( q0(k) < Qgrid(i) ) exit
       end do
       i=i-1

       q=q0(k)-Qgrid(i)
       do j=0,NumQgrid
          pol_drv(k,j) = bmat(i,j) + 2.0d0*cmat(i,j)*q + 3.0d0*dmat(i,j)*q*q
       end do

! ---------------------------- Wu & Gygi ---
       if ( q0(k) /= 0.0 ) then
          do j=0,NumQgrid
             pol_drv(k,j)=Qgrid(j)*(-pol(k,j)+q0(k)*pol_drv(k,j))/q0(k)**2
          end do
       else
          pol_drv(k,:)=0.0d0
       end if
! ------------------------------------------

    end do ! k

  END SUBROUTINE calc_pol_drv


  SUBROUTINE calc_dqdn( rgrid, rho, grad16, fin, fou )
    implicit none
    type( grid ) :: rgrid
    type( GSArray ) :: rho
    type( gradient16 ) :: grad16
    real(8),intent(IN)  :: fin(:)
    real(8),intent(OUT) :: fou(:)
    type( fd ) :: nabla
    type( lattice ) :: aa, bb
    real(8) :: trho,rho_tmp(2),vdx_lda(2),vdc_lda(2)
    real(8) :: edx_lda(2),edc_lda
    real(QP) :: pi,FouPiThr,twothr,onethr,sevthr,kf,ss,cm,q0,c,d,const
    real(QP) :: qtrho,b(3,3)
    real(QP),allocatable :: rtmp0(:),rtmp(:,:),ftmp(:)
    integer :: m0,m1,m,ML1,ML2,ML3,Md,i,i1,i2,i3,j,j1,j2,j3,k1,k2,k3,info
    integer :: nspin
    integer,allocatable :: LLL(:,:,:)

    pi       = acos(-1.0q0)
    onethr   = 1.0q0/3.0q0
    twothr   = 2.0q0/3.0q0
    FouPiThr = 4.0q0*pi/3.0q0
    sevthr   = 7.0q0/3.0q0
    const    = Zab/9.0q0

    m0 = rgrid%g1%head
    m1 = rgrid%g1%tail

    ML1 = rgrid%g3%x%size_global
    ML2 = rgrid%g3%y%size_global
    ML3 = rgrid%g3%z%size_global

    nspin = rho%s_range%size_global

    allocate( rtmp(rgrid%g1%size_global,3) ) ; rtmp =0.0q0
    allocate( rtmp0(m0:m1)                 ) ; rtmp0=0.0q0
    allocate( ftmp(m0:m1)                  ) ; ftmp =0.0q0

! ---

    do i=m0,m1

       rho_tmp(1) = rho%val(i,1)
       rho_tmp(2) = rho%val(i,nspin)
       trho       = sum( rho_tmp(1:nspin) )
       qtrho      = trho

       if ( trho <= 0.0d0 ) cycle
!       if ( trho <= zero_density ) cycle

       call calc_edx_lda( nspin, rho_tmp, edx_lda, vdx_lda )
       call calc_edc_lda( nspin, rho_tmp, edc_lda, vdc_lda )

       kf = (3.0q0*pi*pi*qtrho)**onethr
       ss = grad16%gg(i)/(2.0q0*kf*qtrho)**2
       q0 = -FouPithr*(edx_lda(1)+edc_lda) - const*ss*kf

       c=1.0q0
       do j=1,Mmax-1
          c=c+(q0/Qmax)**j
       end do
       d=0.0q0
       do j=1,Mmax
          d=d+(q0/Qmax)**j/j
       end do
       c=c*exp(-d)
       
       ftmp(i) = c*( -FouPiThr*( vdx_lda(1) + vdc_lda(1) )*trho &
                     +const*sevthr*kf*ss )*fin(i-m0+1)

       rtmp0(i) = c*( -const/(2.0q0*qtrho*kf) )*fin(i-m0+1)

    end do ! i

! ---

    rtmp(m0:m1,1) = rtmp0(m0:m1) * grad16%gx(m0:m1)
    rtmp(m0:m1,2) = rtmp0(m0:m1) * grad16%gy(m0:m1)
    rtmp(m0:m1,3) = rtmp0(m0:m1) * grad16%gz(m0:m1)

    deallocate( rtmp0 )

    do i=1,3
       call MPI_ALLGATHERV( rtmp(m0,i), m1-m0+1, MPI_REAL16, rtmp(1,i) &
            , ir_grid, id_grid, MPI_REAL16, comm_grid, info )
    end do

    allocate( LLL(0:ML1-1,0:ML2-1,0:ML3-1) ) ; LLL=0

    i=m0-1
    do i3=rgrid%g3%z%head,rgrid%g3%z%tail
    do i2=rgrid%g3%y%head,rgrid%g3%y%tail
    do i1=rgrid%g3%x%head,rgrid%g3%x%tail
       i=i+1
       LLL(i1,i2,i3) = i
    end do
    end do
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, LLL, size(LLL), MPI_INTEGER &
         ,MPI_SUM, comm_grid, info )

    call construct_nabla_fd( nabla )
    Md = nabla%md

    call get_aa_lattice( aa )
    call get_reciprocal_lattice( aa, bb )
    b(1:3,1)=aa%Length(1)*bb%LatticeVector(1:3,1)/( 2*pi*rgrid%spacing(1) )
    b(1:3,2)=aa%Length(2)*bb%LatticeVector(1:3,2)/( 2*pi*rgrid%spacing(2) )
    b(1:3,3)=aa%Length(3)*bb%LatticeVector(1:3,3)/( 2*pi*rgrid%spacing(3) )

    do i3=0,ML3-1
    do i2=0,ML2-1
    do i1=0,ML1-1
       i=LLL(i1,i2,i3)
       do m=-Md,Md
          cm=nabla%coef(m)*sign(1,m)
          j1=i1+m
          k1=j1/ML1 ; if ( j1 < 0 ) k1=(j1+1)/ML1-1
          j1=j1-k1*ML1
          j =LLL(j1,i2,i3)
          if ( m0 <= j .and. j <= m1 ) then
             ftmp(j) = ftmp(j) + cm*( rtmp(i,1)*b(1,1) &
                                     +rtmp(i,2)*b(2,1) &
                                     +rtmp(i,3)*b(3,1) )
          end if
          j2=i2+m
          k2=j2/ML2 ; if ( j2 < 0 ) k2=(j2+1)/ML2-1
          j2=j2-k2*ML2
          j =LLL(i1,j2,i3)
          if ( m0 <= j .and. j <= m1 ) then
             ftmp(j) = ftmp(j) + cm*( rtmp(i,1)*b(1,2) &
                                     +rtmp(i,2)*b(2,2) &
                                     +rtmp(i,3)*b(3,2) )
          end if
          j3=i3+m
          k3=j3/ML3 ; if ( j3 < 0 ) k3=(j3+1)/ML3-1
          j3=j3-k3*ML3
          j =LLL(i1,i2,j3)
          if ( m0 <= j .and. j <= m1 ) then
             ftmp(j) = ftmp(j) + cm*( rtmp(i,1)*b(1,3) &
                                     +rtmp(i,2)*b(2,3) &
                                     +rtmp(i,3)*b(3,3) )
          end if
       end do ! m
    end do ! i1
    end do ! i2
    end do ! i3

    fou(1:m1-m0+1) = ftmp(m0:m1)

! ---

    deallocate( ftmp )
    deallocate( LLL  )
    deallocate( rtmp )

    call destruct_nabla_fd( nabla )

  END SUBROUTINE calc_dqdn


END MODULE xc_vdw_module
