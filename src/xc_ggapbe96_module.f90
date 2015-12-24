MODULE xc_ggapbe96_module

  use aa_module
  use bb_module
  use bc_module
  use parallel_module
  use fd_module
  use xc_ggapbe96_mol_module

  implicit none

  PRIVATE
  PUBLIC :: init_GGAPBE96, calc_GGAPBE96

  integer :: Igrid(2,0:3)
  integer :: ML_0, ML_1, MSP_0, MSP_1, MSP, comm
  real(8) :: dV
  logical :: flag_init = .true.

  real(8),allocatable :: nab(:)

  real(8),allocatable :: gx(:),gy(:),gz(:)
  real(8),parameter :: zero_density = 1.d-10

  real(8) :: b(3,3)

  integer,allocatable :: LLL2(:,:,:)
  integer :: Md,ML,ML1,ML2,ML3
  real(8) :: Hgrid(3)
  integer :: SYStype

  real(8) :: mu=0.21951d0
  real(8) :: Kp=0.804d0
  real(8) :: beta=0.066725d0

CONTAINS


  SUBROUTINE init_GGAPBE96 &
       ( Igrid_in,MSP_0_in,MSP_1_in,MSP_in,comm_in,dV_in &
        ,Md_in,Hgrid_in,Ngrid_in,SYStype_in, mu_in, Kp_in )
    implicit none
    integer,intent(IN) :: Igrid_in(2,0:3),MSP_0_in,MSP_1_in,MSP_in,comm_in
    real(8),intent(IN) :: dV_in, Hgrid_in(3)
    integer,intent(IN) :: Md_in, Ngrid_in(0:3), SYStype_in
    real(8),optional,intent(IN) :: mu_in, Kp_in
    real(8) :: aaL(3), Pi

    if ( .not.flag_init ) return

    Igrid(:,:) = Igrid_in(:,:)
    ML_0       = Igrid(1,0)
    ML_1       = Igrid(2,0)
    MSP_0      = MSP_0_in
    MSP_1      = MSP_1_in
    MSP        = MSP_in
    comm       = comm_in
    dV         = dV_in
    flag_init  = .false.
    Md         = Md_in
    Hgrid(:)   = Hgrid_in(:)
    ML  = Ngrid_in(0)
    ML1 = Ngrid_in(1)
    ML2 = Ngrid_in(2)
    ML3 = Ngrid_in(3)

    SYStype = SYStype_in

    aaL(1) = sqrt( sum(aa(1:3,1)**2) )
    aaL(2) = sqrt( sum(aa(1:3,2)**2) )
    aaL(3) = sqrt( sum(aa(1:3,3)**2) )
    Pi     = acos(-1.0d0)
    b(:,:) = 0.0d0
    b(1:3,1)=aaL(1)*bb(1:3,1)/(2.0d0*Pi)/Hgrid(1)
    b(1:3,2)=aaL(2)*bb(1:3,2)/(2.0d0*Pi)/Hgrid(2)
    b(1:3,3)=aaL(3)*bb(1:3,3)/(2.0d0*Pi)/Hgrid(3)

    allocate( nab(-Md:Md) ) ; nab=0.d0
    call get_coef_nabla_fd( Md, nab )

    if ( present(mu_in)   ) mu  = mu_in
    if ( present(Kp_in)   ) Kp  = Kp_in
    beta = mu*3.0d0/acos(-1.0d0)**2

    flag_init = .false.

  END SUBROUTINE init_GGAPBE96


  SUBROUTINE calc_GGAPBE96 &
       ( rho, Exc_out, Vxc_out, Ex_out, Ec_out, Vx_out, Vc_out )
    implicit none
    real(8),intent(IN)  :: rho(ML_0:,:)
    real(8),intent(OUT) :: Exc_out
    real(8),optional,intent(OUT) :: Vxc_out(ML_0:,MSP_0:)
    real(8),optional,intent(OUT) :: Ex_out, Ec_out
    real(8),optional,intent(OUT) :: Vx_out(ML_0:,MSP_0:),Vc_out(ML_0:,MSP_0:)
    integer :: i1,i2,i3,i,irank,j1,j2,j3,l1,l2,l3,m1,m2,m3,ierr
    real(8) :: s0(2),s1(2), Ex_part, Ec_part
    real(8),allocatable :: vxc_tmp(:,:)

    if ( flag_init ) then
       write(*,*) "Call INIT_GGAPBE96 first"
       stop "stop@calc_GGAPBE96(in xc_ggapbe96_module)"
    end if

    Ex_part = 0.0d0
    Ec_part = 0.0d0

    Exc_out = 0.0d0
    if ( present(Ex_out)  ) Ex_out =0.0d0
    if ( present(Ec_out)  ) Ec_out =0.0d0
    if ( present(Vxc_out) ) Vxc_out=0.0d0
    if ( present(Vx_out)  ) Vx_out =0.0d0
    if ( present(Vc_out)  ) Vc_out =0.0d0

    select case( SYStype )
    case default

    allocate( LLL2(0:ML1-1,0:ML2-1,0:ML3-1) ) ; LLL2=0
    i=0
    irank=-1
    do i3=1,node_partition(3)
    do i2=1,node_partition(2)
    do i1=1,node_partition(1)
       irank=irank+1
       l1=pinfo_grid(1,irank) ; m1=pinfo_grid(2,irank)+l1-1
       l2=pinfo_grid(3,irank) ; m2=pinfo_grid(4,irank)+l2-1
       l3=pinfo_grid(5,irank) ; m3=pinfo_grid(6,irank)+l3-1
       do j3=l3,m3
       do j2=l2,m2
       do j1=l1,m1
          i=i+1
          LLL2(j1,j2,j3)=i
       end do
       end do
       end do
    end do
    end do
    end do

    call construct_gradient( rho )

    case( 1,2 )

       m1 = ( ML1-1 )/2 + Md
       m2 = ( ML2-1 )/2 + Md
       m3 = ( ML3-1 )/2 + Md

       allocate( LLL2(-m1:m1,-m2:m2,-m3:m3) ) ; LLL2=0

       call get_LLL_mol( m1,m2,m3,LLL2 )

       allocate( gx(ML_0:ML_1) ) ; gx=0.0d0
       allocate( gy(ML_0:ML_1) ) ; gy=0.0d0
       allocate( gz(ML_0:ML_1) ) ; gz=0.0d0

       call construct_gradient_mol( Md,nab,rho,gx,gy,gz )

    end select

! --

    allocate( vxc_tmp(ML_0:ML_1,1:MSP) ) ; vxc_tmp=0.0d0

! -- Exchange --

    call calc_GGAPBE96_x( rho, vxc_tmp, Ex_part )

    if ( present(Vx_out) ) then
       Vx_out(ML_0:ML_1,MSP_0:MSP_1) = vxc_tmp(ML_0:ML_1,MSP_0:MSP_1)
    end if

    if ( present(Vxc_out) ) then
       Vxc_out(ML_0:ML_1,MSP_0:MSP_1) = vxc_tmp(ML_0:ML_1,MSP_0:MSP_1)
    end if

! -- Correlation --

!    call calc_GGAPBE96_c( rho, vxc_tmp, Ec_part )
    call calc_GGAPBE96_c2( rho, vxc_tmp, Ec_part )

    if ( present(Vc_out) ) then
       Vc_out(ML_0:ML_1,MSP_0:MSP_1) = vxc_tmp(ML_0:ML_1,MSP_0:MSP_1)
    end if

    if ( present(Vxc_out) ) then
       Vxc_out(ML_0:ML_1,MSP_0:MSP_1) = Vxc_out(ML_0:ML_1,MSP_0:MSP_1) &
                                      + vxc_tmp(ML_0:ML_1,MSP_0:MSP_1)
    end if

! --

    s0(1) = Ex_part*dV
    s0(2) = Ec_part*dV
    call mpi_allreduce(s0,s1,2,mpi_real8,mpi_sum,comm,ierr)

    s0(1) = s1(1)
    call mpi_allreduce(s0,s1,1,mpi_real8,mpi_sum,comm_spin,ierr)

    Exc_out = s1(1) + s1(2)

    if ( present(Ex_out) ) Ex_out=s1(1)
    if ( present(Ec_out) ) Ec_out=s1(2)

! --

    deallocate( vxc_tmp )

    deallocate( gz, gy, gx )

    deallocate( LLL2 )

  END SUBROUTINE calc_GGAPBE96


  SUBROUTINE construct_gradient(rho)
    implicit none
    real(8),intent(IN) :: rho(ML_0:,:)
    integer :: i,i1,i2,i3,s,m
    real(8) :: g1,g2,g3
    allocate( gx(ML_0:ML_1) ) ; gx=0.0d0
    allocate( gy(ML_0:ML_1) ) ; gy=0.0d0
    allocate( gz(ML_0:ML_1) ) ; gz=0.0d0
    www(:,:,:,:)=0.0d0
    do s=1,MSP
       i=ML_0-1
       do i3=Igrid(1,3),Igrid(2,3)
       do i2=Igrid(1,2),Igrid(2,2)
       do i1=Igrid(1,1),Igrid(2,1)
          i=i+1
          www(i1,i2,i3,1) = www(i1,i2,i3,1) + rho(i,s)
       end do
       end do
       end do
    end do

    call bcset(1,1,Md,0)

    i=ML_0-1
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       g1=0.0d0
       g2=0.0d0
       g3=0.0d0
       do m=1,Md
          g1 = g1 - nab(m)*( www(i1-m,i2,i3,1) - www(i1+m,i2,i3,1) )
          g2 = g2 - nab(m)*( www(i1,i2-m,i3,1) - www(i1,i2+m,i3,1) )
          g3 = g3 - nab(m)*( www(i1,i2,i3-m,1) - www(i1,i2,i3+m,1) )
       end do
       i=i+1
       gx(i) = b(1,1)*g1 + b(1,2)*g2 + b(1,3)*g3
       gy(i) = b(2,1)*g1 + b(2,2)*g2 + b(2,3)*g3
       gz(i) = b(3,1)*g1 + b(3,2)*g2 + b(3,3)*g3
    end do
    end do
    end do
  END SUBROUTINE construct_gradient


  SUBROUTINE calc_GGAPBE96_x( rho, vex, Ex )
    implicit none
    real(8),intent(IN)  :: rho(ML_0:,:)
    real(8),intent(OUT) :: vex(ML_0:,:)
    real(8),intent(OUT) :: Ex
!   real(8),parameter :: mu=0.21951d0, Kp=0.804d0
    integer :: i,ispin,ierr
    real(8) :: rhoa,rhob,trho,cm
    real(8),allocatable :: rrrr(:,:),rtmp(:)
    real(8) :: kf, vx_lda, ex_lda, Fx, Pi, g2, factor
    integer :: i1,i2,i3,j,j1,j2,j3,k1,k2,k3,m,m1,m2,m3

    Pi = acos(-1.0d0)

    factor = 1.0d0
    if ( MSP == 2 ) factor = 2.0d0

    vex = 0.0d0
    Ex  = 0.0d0

    allocate( rtmp(ML_0:ML_1) ) ; rtmp=0.0d0
    allocate( rrrr(ML,3)      ) ; rrrr=0.0d0

    do ispin=MSP_0,MSP_1

       do i=ML_0,ML_1

          trho = factor*rho(i,ispin)

          if ( trho <= zero_density ) cycle

          kf = (3.0d0*Pi*Pi*trho)**(1.0d0/3.0d0)

          ex_lda = -3.0d0/(4.0d0*Pi)*kf
          vx_lda = -1.0d0/Pi*kf

          g2 = gx(i)*gx(i) + gy(i)*gy(i) + gz(i)*gz(i)

          Fx = 1.0d0 + Kp - 4.0d0*Kp*Kp*(trho*kf)**2 &
                            /( 4.0d0*Kp*(trho*kf)**2 + mu*g2 )

          Ex = Ex + trho*ex_lda*Fx

          vex(i,ispin) = vex(i,ispin) &
               + Fx*vx_lda + ( 24.0d0*Pi*Kp*Kp*mu*trho**3*g2 ) &
                            /( 4.0d0*Kp*(trho*kf)**2 + mu*g2 )**2

          rtmp(i) = -18.0d0*Pi*Kp*Kp*mu*trho**4 &
                    /( 4.0d0*Kp*(trho*kf)**2 + mu*g2 )**2

       end do ! i

       rrrr(ML_0:ML_1,1) = rtmp(ML_0:ML_1)*gx(ML_0:ML_1)
       call mpi_allgatherv(rrrr(ML_0,1),ir_grid(myrank_g),mpi_real8 &
            ,rrrr(1,1),ir_grid,id_grid,mpi_real8,comm_grid,ierr)

       rrrr(ML_0:ML_1,2) = rtmp(ML_0:ML_1)*gy(ML_0:ML_1)
       call mpi_allgatherv(rrrr(ML_0,2),ir_grid(myrank_g),mpi_real8 &
            ,rrrr(1,2),ir_grid,id_grid,mpi_real8,comm_grid,ierr)

       rrrr(ML_0:ML_1,3)=rtmp(ML_0:ML_1)*gz(ML_0:ML_1)
       call mpi_allgatherv(rrrr(ML_0,3),ir_grid(myrank_g),mpi_real8 &
            ,rrrr(1,3),ir_grid,id_grid,mpi_real8,comm_grid,ierr)

       select case( SYStype )
       case default

       do i3=0,ML3-1
       do i2=0,ML2-1
       do i1=0,ML1-1
          i=LLL2(i1,i2,i3)
          do m=-Md,Md
             cm=nab(m)*sign(1,m)
             j1=i1+m
             k1=j1/ML1 ; if ( j1<0 ) k1=(j1+1)/ML1-1
             j1=j1-k1*ML1
             j =LLL2(j1,i2,i3)
! The potential vex is calculated at j-th grid point rather than i-th.
! This is because non-transposed nabla matrix Dij is used (See XC.doc).
             if ( ML_0 <= j .and. j <= ML_1 ) then
                vex(j,ispin) = vex(j,ispin) + cm*( rrrr(i,1)*b(1,1) &
                                                  +rrrr(i,2)*b(2,1) &
                                                  +rrrr(i,3)*b(3,1) )
             end if
             j2=i2+m
             k2=j2/ML2 ; if ( j2<0 ) k2=(j2+1)/ML2-1
             j2=j2-k2*ML2
             j =LLL2(i1,j2,i3)
             if ( ML_0 <= j .and. j <= ML_1 ) then
                vex(j,ispin) = vex(j,ispin) + cm*( rrrr(i,1)*b(1,2) &
                                                  +rrrr(i,2)*b(2,2) &
                                                  +rrrr(i,3)*b(3,2) )
             end if
             j3=i3+m
             k3=j3/ML3 ; if ( j3<0 ) k3=(j3+1)/ML3-1
             j3=j3-k3*ML3
             j =LLL2(i1,i2,j3)
             if ( ML_0 <= j .and. j <= ML_1 ) then
                vex(j,ispin) = vex(j,ispin) + cm*( rrrr(i,1)*b(1,3) &
                                                  +rrrr(i,2)*b(2,3) &
                                                  +rrrr(i,3)*b(3,3) )
             end if
          end do ! m
       end do ! i1
       end do ! i2
       end do ! i3

       case( 1,2 )

          m1 = (ML1-1)/2 + Md
          m2 = (ML2-1)/2 + Md
          m3 = (ML3-1)/2 + Md

          call calc_ve_mol &
               ( ML_0, ML_1, m1,m2,m3, Md, nab, vex(:,ispin), LLL2, rrrr )

       end select

    end do ! ispin

    Ex = Ex/factor

    deallocate( rrrr )
    deallocate( rtmp )

  END SUBROUTINE calc_GGAPBE96_x


  SUBROUTINE calc_GGAPBE96_c( rho, vco, Ec )
    implicit none
    real(8),intent(IN)  :: rho(ML_0:,:)
    real(8),intent(OUT) :: vco(ML_0:,:), Ec
    real(8),parameter :: A00  =0.031091d0,A01  =0.015545d0,A02  =0.016887d0
    real(8),parameter :: alp10=0.21370d0 ,alp11=0.20548d0 ,alp12=0.11125d0
    real(8),parameter :: bt10 =7.5957d0  ,bt11 =1.41189d1 ,bt12 =1.0357d1
    real(8),parameter :: bt20 =3.5876d0  ,bt21 =6.1977d0  ,bt22 =3.6231d0
    real(8),parameter :: bt30 =1.6382d0  ,bt31 =3.3662d0  ,bt32 =0.88026d0
    real(8),parameter :: bt40 =0.49294d0 ,bt41 =0.62517d0 ,bt42 =0.49671d0
!   real(8),parameter :: C1=2.14611945579107d0,C2=0.031091d0
    real(8) :: c1,C2
    real(8) :: const1, const2, factor
    integer :: i,j,ispin,m,i1,i2,i3,j1,j2,j3,k1,k2,k3,ierr,m1,m2,m3
    real(8) :: trho, rhoa, rhob
    real(8) :: kf, rs, ec_U, ec_P, ec_lda, phi, g2
    real(8) :: dac_dn, dfz_dz, deU_dn, deP_dn, H, A, T, alpc, fz
    real(8) :: dec_dz, dphi_dz, dH_dA, dH_dT, tmp, dH_dphi, Ai
    real(8) :: dA_dn, dec_dn, Hs, zeta
    real(8) :: dz_dn(2), cm
    real(8),allocatable :: rrrr(:,:), rtmp(:), vtmp(:)
    real(8) :: Pi, one, two, fouthr, onethr
    real(8) :: sevthr, twothr, thrtwo, ThrFouPi

    const1=2.d0**(4.d0/3.d0)-2.d0
    const2=9.d0*(2.d0**(1.d0/3.d0)-1.d0)/4.d0

    factor = 1.0d0
    if ( MSP == 1 ) factor = 0.5d0
 
    Pi       = acos(-1.0d0)
    one      = 1.0d0
    two      = 2.0d0
    fouthr   = 4.0d0/3.0d0
    onethr   = 1.0d0/3.0d0
    ThrFouPi = 3.0d0/(4.0d0*Pi)
    thrtwo   = 3.0d0/2.0d0
    twothr   = 2.0d0/3.0d0
    sevthr   = 7.0d0/3.0d0

    C2 = ( 1.0d0-log(2.0d0) )/Pi**2 ! "gamma" in PBE paper
    C1 = beta/C2

    vco = 0.0d0
    Ec  = 0.0d0

    allocate( rtmp(ML_0:ML_1) ) ; rtmp=0.0d0
    allocate( rrrr(ML,3)      ) ; rrrr=0.0d0

    do i=ML_0,ML_1

       rhoa = rho(i,1  )*factor
       rhob = rho(i,MSP)*factor
       trho = rhoa + rhob

       if ( trho <= zero_density ) cycle

       zeta = ( rhoa - rhob )/trho

       fz = ( (one+zeta)**fouthr + (one-zeta)**fouthr - two )*const1

       kf = ( 3.0d0*Pi*Pi*trho )**onethr

       rs = ( ThrFouPi/trho )**onethr

       ec_U = -2.0d0*A00*( 1.0d0 + alp10*rs ) &
            *log( 1.0d0 + 1.0d0/( 2.0d0*A00 &
            *( bt10*sqrt(rs) + bt20*rs + bt30*rs**thrtwo + bt40*rs*rs) ) )
       ec_P = -2.0d0*A01*( 1.0d0 + alp11*rs ) &
            *log( 1.0d0 + 1.0d0/( 2.0d0*A01 &
            *( bt11*sqrt(rs) + bt21*rs + bt31*rs**thrtwo + bt41*rs*rs) ) )
       alpc = -2.0d0*A02*( 1.0d0 + alp12*rs ) &
            *log( 1.0d0 + 1.0d0/( 2.0d0*A02 &
            *( bt12*sqrt(rs) + bt22*rs + bt32*rs**thrtwo + bt42*rs*rs) ) )

       ec_lda = ec_U - alpc*fz*const2*(one-zeta**4) &
            + ( ec_P - ec_U )*fz*zeta**4

       phi = 0.5d0*( (one+zeta)**twothr + (one-zeta)**twothr )

       T = ( gx(i)*gx(i) + gy(i)*gy(i) + gz(i)*gz(i) )*Pi &
          /( 16.0d0*phi**2*kf*trho**2 )

       Ai = ( exp(-ec_lda/(C2*phi**3)) - 1.0d0 )/C1

       Hs = C2*phi**3*log( one + C1*( Ai*Ai*T + Ai*T*T ) &
                                   /( Ai*Ai + Ai*T + T*T ) )

       Ec = Ec + trho*( ec_lda + Hs )

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

       dec_dz = -alpc*dfz_dz*const2*(one-zeta**4) &
               + 4.0d0*alpc*fz*const2*zeta**3 &
               +(ec_P-ec_U)*dfz_dz*zeta**4 &
               +(ec_P-ec_U)*fz*4.0d0*zeta**3

       if ( zeta == 1.0d0 .or. zeta == -1.0d0 ) then
          dphi_dz = 0.0d0
       else
          dphi_dz = ( (one+zeta)**(-onethr)-(one-zeta)**(-onethr) )/3.d0
       end if

       tmp   = Ai*Ai + Ai*T + T*T
       dH_dA = -phi**3*C1*C2*Ai*T**3*(2.0d0*Ai+T) &
               /( tmp*tmp + C1*T*(Ai*Ai+Ai*T)*tmp )
       dH_dT =  phi**3*C1*C2*Ai**3*(Ai+2.0d0*T) &
               /( tmp*tmp + C1*T*(Ai*Ai+Ai*T)*tmp )

       dH_dphi = 3.0d0*Hs/phi

       dz_dn(1)   = 2.0d0*rhob/trho
       dz_dn(MSP) =-2.0d0*rhoa/trho

       do ispin=MSP_0,MSP_1

          dec_dn = deU_dn - dac_dn*fz*const2*(one-zeta**4) &
               +(deP_dn-deU_dn)*fz*zeta**4 + dec_dz*dz_dn(ispin)

          tmp = exp( -ec_lda/(phi**3*C2) )

          dA_dn = tmp/(C1*C2*phi**3) &
               *( dec_dn - 3.0d0*ec_lda/phi*dphi_dz*dz_dn(ispin) )

          vco(i,ispin) = vco(i,ispin) + ec_lda + Hs + trho*dec_dn &
               + trho*dH_dA*dA_dn - sevthr*T*dH_dT &
               + trho*dH_dphi*dphi_dz*dz_dn(ispin)

       end do ! ispin

       rtmp(i) = dH_dT*Pi/(8.0d0*kf*trho)

    end do ! i

    rrrr(ML_0:ML_1,1) = rtmp(ML_0:ML_1)*gx(ML_0:ML_1)
    call mpi_allgatherv(rrrr(ML_0,1),ir_grid(myrank_g),mpi_real8 &
         ,rrrr(1,1),ir_grid,id_grid,mpi_real8,comm_grid,ierr)

    rrrr(ML_0:ML_1,2) = rtmp(ML_0:ML_1)*gy(ML_0:ML_1)
    call mpi_allgatherv(rrrr(ML_0,2),ir_grid(myrank_g),mpi_real8 &
         ,rrrr(1,2),ir_grid,id_grid,mpi_real8,comm_grid,ierr)

    rrrr(ML_0:ML_1,3) = rtmp(ML_0:ML_1)*gz(ML_0:ML_1)
    call mpi_allgatherv(rrrr(ML_0,3),ir_grid(myrank_g),mpi_real8 &
         ,rrrr(1,3),ir_grid,id_grid,mpi_real8,comm_grid,ierr)

    select case( SYStype )
    case default

    do i3=0,ML3-1
    do i2=0,ML2-1
    do i1=0,ML1-1
       i=LLL2(i1,i2,i3)
       do m=-Md,Md
          cm=nab(m)*sign(1,m)
          j1=i1+m
          k1=j1/ML1 ; if ( j1<0 ) k1=(j1+1)/ML1-1
          j1=j1-k1*ML1
          j =LLL2(j1,i2,i3)
          if ( ML_0 <= j .and. j <= ML_1 ) then
             do ispin=MSP_0,MSP_1
                vco(j,ispin) = vco(j,ispin) + cm*( rrrr(i,1)*b(1,1) &
                                                  +rrrr(i,2)*b(2,1) &
                                                  +rrrr(i,3)*b(3,1) )
             end do
          end if
          j2=i2+m
          k2=j2/ML2 ; if ( j2<0 ) k2=(j2+1)/ML2-1
          j2=j2-k2*ML2
          j =LLL2(i1,j2,i3)
          if ( ML_0 <= j .and. j <= ML_1 ) then
             do ispin=MSP_0,MSP_1
                vco(j,ispin) = vco(j,ispin) + cm*( rrrr(i,1)*b(1,2) &
                                                  +rrrr(i,2)*b(2,2) &
                                                  +rrrr(i,3)*b(3,2) )
             end do
          end if
          j3=i3+m
          k3=j3/ML3 ; if ( j3<0 ) k3=(j3+1)/ML3-1
          j3=j3-k3*ML3
          j =LLL2(i1,i2,j3)
          if ( ML_0 <= j .and. j <= ML_1 ) then
             do ispin=MSP_0,MSP_1
                vco(j,ispin) = vco(j,ispin) + cm*( rrrr(i,1)*b(1,3) &
                                                  +rrrr(i,2)*b(2,3) &
                                                  +rrrr(i,3)*b(3,3) )
             end do
          end if
       end do ! m
    end do ! i1
    end do ! i2
    end do ! i3

    case( 1,2 )

       m1 = (ML1-1)/2 + Md
       m2 = (ML2-1)/2 + Md
       m3 = (ML3-1)/2 + Md

       allocate( vtmp(ML_0:ML_1) ) ; vtmp=0.0d0

       call calc_ve_mol &
            ( ML_0, ML_1, m1,m2,m3, Md, nab, vtmp, LLL2, rrrr )

       do ispin=MSP_0,MSP_1
          vco(:,ispin) = vco(:,ispin) + vtmp(:)
       end do

       deallocate( vtmp )

    end select

    deallocate( rtmp )
    deallocate( rrrr )

    return
  END SUBROUTINE calc_GGAPBE96_c


  SUBROUTINE calc_GGAPBE96_c2( rho, vco, Ec )
    implicit none
    real(8),intent(IN)  :: rho(ML_0:,:)
    real(8),intent(OUT) :: vco(ML_0:,:), Ec
    real(8),parameter :: A00  =0.031091d0,A01  =0.015545d0,A02  =0.016887d0
    real(8),parameter :: alp10=0.21370d0 ,alp11=0.20548d0 ,alp12=0.11125d0
    real(8),parameter :: bt10 =7.5957d0  ,bt11 =1.41189d1 ,bt12 =1.0357d1
    real(8),parameter :: bt20 =3.5876d0  ,bt21 =6.1977d0  ,bt22 =3.6231d0
    real(8),parameter :: bt30 =1.6382d0  ,bt31 =3.3662d0  ,bt32 =0.88026d0
    real(8),parameter :: bt40 =0.49294d0 ,bt41 =0.62517d0 ,bt42 =0.49671d0
!   real(8),parameter :: C1=2.14611945579107d0,C2=0.031091d0
    real(8) :: C1,C2
    real(8) :: const1, const2, factor
    integer :: i,j,ispin,m,i1,i2,i3,j1,j2,j3,k1,k2,k3,ierr,m1,m2,m3
    real(8) :: trho, rhoa, rhob
    real(8) :: kf, rs, ec_U, ec_P, ec_lda, phi, g2
    real(8) :: dac_dn, dfz_dz, deU_dn, deP_dn, H, A, T, alpc, fz
    real(8) :: dec_dz, dphi_dz, dH_dA, dH_dT, tmp, dH_dphi, Ai
    real(8) :: dA_dn, dec_dn, Hs, zeta, dT_dphi
    real(8) :: dz_dn(2), cm, rssq
    real(8) :: dac_drs, deP_drs, deU_drs, drs_dn
    real(8),allocatable :: rrrr(:,:), rtmp(:),vtmp(:)
    real(8) :: Pi, one, two, fouthr, onethr
    real(8) :: sevthr, twothr, thrtwo, ThrFouPi

    const1 = 2.0d0**(4.0d0/3.0d0)-2.0d0
    const2 = 9.0d0*(2.0d0**(1.0d0/3.0d0)-1.0d0)/4.0d0

    factor = 1.0d0
    if ( MSP == 1 ) factor = 0.5d0
 
    Pi       = acos(-1.0d0)
    one      = 1.0d0
    two      = 2.0d0
    fouthr   = 4.0d0/3.0d0
    onethr   = 1.0d0/3.0d0
    ThrFouPi = 3.0d0/(4.0d0*Pi)
    thrtwo   = 3.0d0/2.0d0
    twothr   = 2.0d0/3.0d0
    sevthr   = 7.0d0/3.0d0

    C2 = ( 1.0d0-log(2.0d0) )/Pi**2 ! "gamma" in PBE paper
    C1 = beta/C2

    vco = 0.0d0
    Ec  = 0.0d0

    allocate( rtmp(ML_0:ML_1) ) ; rtmp=0.0d0
    allocate( rrrr(ML,3)      ) ; rrrr=0.0d0

    do i=ML_0,ML_1

       rhoa = rho(i,1  )*factor
       rhob = rho(i,MSP)*factor
       trho = rhoa + rhob

       if ( trho <= zero_density ) cycle

       zeta = ( rhoa - rhob )/trho

       fz = ( (one+zeta)**fouthr + (one-zeta)**fouthr - two )*const1

       kf = ( 3.0d0*Pi*Pi*trho )**onethr

       rs = ( ThrFouPi/trho )**onethr

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

       ec_lda = ec_U - alpc*fz*const2*(one-zeta**4) &
            + ( ec_P - ec_U )*fz*zeta**4

       phi = 0.5d0*( (one+zeta)**twothr + (one-zeta)**twothr )

       T = Pi*( gx(i)*gx(i) + gy(i)*gy(i) + gz(i)*gz(i) ) &
             /( 16.0d0*phi**2*kf*trho**2 )

       A = C1/( exp( -ec_lda/(C2*phi**3) ) - one )

       Hs = C2 * phi**3 * log( one + C1*T*(one+A*T)/(one+A*T+(A*T)**2) )

       Ec = Ec + trho*( ec_lda + Hs )

       drs_dn = -4.0d0*Pi*rs**4/9.0d0

       tmp = bt10*rssq + bt20*rs + bt30*rs*rssq + bt40*rs*rs
       deU_drs = alp10*ec_U/( one + alp10*rs ) &
            +A00*( one + alp10*rs )/rssq &
            *( bt10 + two*bt20*rssq + 3.0d0*bt30*rs + 4.0d0*bt40*rs*rssq ) &
            /( two*A00*tmp*tmp+tmp )

       tmp = bt11*rssq + bt21*rs + bt31*rs*rssq + bt41*rs*rs
       deP_drs = alp11*ec_P/( one + alp11*rs ) &
            +A01*( one + alp11*rs )/rssq &
            *( bt11 + two*bt21*rssq + 3.0d0*bt31*rs + 4.0d0*bt41*rs*rssq ) &
            /( two*A01*tmp*tmp+tmp )

       tmp = bt12*rssq + bt22*rs + bt32*rs*rssq + bt42*rs*rs
       dac_drs = alp12*alpc/( one + alp12*rs ) &
            +A02*( one + alp12*rs )/rssq &
            *( bt12 + two*bt22*rssq + 3.0d0*bt32*rs + 4.0d0*bt42*rs*rssq ) &
            /( two*A02*tmp*tmp+tmp )

       deU_dn = deU_drs * drs_dn
       deP_dn = deP_drs * drs_dn
       dac_dn = dac_drs * drs_dn


       dfz_dz = fouthr*( (one+zeta)**onethr - (one-zeta)**onethr )*const1

       dec_dz = -alpc*dfz_dz*const2*(one-zeta**4) &
               + 4.0d0*alpc*fz*const2*zeta**3 &
               +(ec_P-ec_U)*( dfz_dz*zeta**4 + fz*4.0d0*zeta**3 )

!       if ( zeta == 1.0d0 .or. zeta == -1.0d0 ) then
!          dphi_dz = 0.0d0
!       else
          dphi_dz = ( (one+zeta)**(-onethr)-(one-zeta)**(-onethr) )/3.d0
!       end if

       tmp = one + A*T + (A*T)**2

       dH_dA = -phi**3*C1*C2*A*T**3*(two+A*T)/(tmp**2+C1*T*(one+A*T)*tmp)

       dH_dT =  phi**3*C1*C2*(one+two*A*T)/(tmp**2+C1*T*(one+A*T)*tmp)

       dH_dphi = 3.0d0*Hs/phi

       dz_dn(1)   = 2.0d0*rhob/trho**2
       dz_dn(MSP) =-2.0d0*rhoa/trho**2

       do ispin=MSP_0,MSP_1

          dec_dn = deU_dn - dac_dn*fz*const2*(one-zeta**4) &
               +(deP_dn-deU_dn)*fz*zeta**4 + dec_dz*dz_dn(ispin)

          dA_dn = A*(C1+A)/(C1*C2*phi**3) &
               *( dec_dn - 3.0d0*ec_lda/phi*dphi_dz*dz_dn(ispin) )

          vco(i,ispin) = vco(i,ispin) + ec_lda + Hs + trho*dec_dn &
               + trho*dH_dA*dA_dn &
               + dH_dT*(-sevthr*T - trho*two*T/phi*dphi_dz*dz_dn(ispin) ) &
               + trho*dH_dphi*dphi_dz*dz_dn(ispin)

       end do ! ispin

       rtmp(i) = dH_dT*Pi/(8.0d0*kf*trho*phi**2)

    end do ! i

    rrrr(ML_0:ML_1,1) = rtmp(ML_0:ML_1)*gx(ML_0:ML_1)
    call mpi_allgatherv(rrrr(ML_0,1),ir_grid(myrank_g),mpi_real8 &
         ,rrrr(1,1),ir_grid,id_grid,mpi_real8,comm_grid,ierr)

    rrrr(ML_0:ML_1,2) = rtmp(ML_0:ML_1)*gy(ML_0:ML_1)
    call mpi_allgatherv(rrrr(ML_0,2),ir_grid(myrank_g),mpi_real8 &
         ,rrrr(1,2),ir_grid,id_grid,mpi_real8,comm_grid,ierr)

    rrrr(ML_0:ML_1,3) = rtmp(ML_0:ML_1)*gz(ML_0:ML_1)
    call mpi_allgatherv(rrrr(ML_0,3),ir_grid(myrank_g),mpi_real8 &
         ,rrrr(1,3),ir_grid,id_grid,mpi_real8,comm_grid,ierr)

    select case( SYStype )
    case default

    do i3=0,ML3-1
    do i2=0,ML2-1
    do i1=0,ML1-1
       i=LLL2(i1,i2,i3)
       do m=-Md,Md
          cm=nab(m)*sign(1,m)
          j1=i1+m
          k1=j1/ML1 ; if ( j1<0 ) k1=(j1+1)/ML1-1
          j1=j1-k1*ML1
          j =LLL2(j1,i2,i3)
! The potential vex is calculated at j-th grid point rather than i-th.
! This is because non-transposed nabla matrix Dij is used (See XC.doc).
          if ( ML_0 <= j .and. j <= ML_1 ) then
             do ispin=MSP_0,MSP_1
                vco(j,ispin) = vco(j,ispin) + cm*( rrrr(i,1)*b(1,1) &
                                                  +rrrr(i,2)*b(2,1) &
                                                  +rrrr(i,3)*b(3,1) )
             end do
          end if
          j2=i2+m
          k2=j2/ML2 ; if ( j2<0 ) k2=(j2+1)/ML2-1
          j2=j2-k2*ML2
          j =LLL2(i1,j2,i3)
          if ( ML_0 <= j .and. j <= ML_1 ) then
             do ispin=MSP_0,MSP_1
                vco(j,ispin) = vco(j,ispin) + cm*( rrrr(i,1)*b(1,2) &
                                                  +rrrr(i,2)*b(2,2) &
                                                  +rrrr(i,3)*b(3,2) )
             end do
          end if
          j3=i3+m
          k3=j3/ML3 ; if ( j3<0 ) k3=(j3+1)/ML3-1
          j3=j3-k3*ML3
          j =LLL2(i1,i2,j3)
          if ( ML_0 <= j .and. j <= ML_1 ) then
             do ispin=MSP_0,MSP_1
                vco(j,ispin) = vco(j,ispin) + cm*( rrrr(i,1)*b(1,3) &
                                                  +rrrr(i,2)*b(2,3) &
                                                  +rrrr(i,3)*b(3,3) )
             end do
          end if
       end do ! m
    end do ! i1
    end do ! i2
    end do ! i3

    case( 1,2 )

       m1 = (ML1-1)/2 + Md
       m2 = (ML2-1)/2 + Md
       m3 = (ML3-1)/2 + Md

       allocate( vtmp(ML_0:ML_1) ) ; vtmp=0.0d0

       call calc_ve_mol &
            ( ML_0, ML_1, m1,m2,m3, Md, nab, vtmp, LLL2, rrrr )

       do ispin=MSP_0,MSP_1
          vco(:,ispin) = vco(:,ispin) + vtmp(:)
       end do

       deallocate( vtmp )

    end select

    deallocate( rtmp )
    deallocate( rrrr )

    return
  END SUBROUTINE calc_GGAPBE96_c2


END MODULE xc_ggapbe96_module
