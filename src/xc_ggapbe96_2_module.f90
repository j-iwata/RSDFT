MODULE xc_ggapbe96_2_module

  use gradient_module
  use grid_module, only: grid, get_map_3d_to_1d_grid
  use xc_variables, only: xcpot, xcene
  use fd_module, only: fd, construct_nabla_fd
  use lattice_module, only: lattice, get_aa_lattice, get_reciprocal_lattice
  use parallel_module
  use basic_type_factory

  implicit none

  PRIVATE
  PUBLIC :: calc_GGAPBE96_2

  integer,parameter :: DP=kind(0.0d0)
  integer,parameter :: QP=kind(0.0q0)

  real(QP),parameter :: zero_density = 1.q-10
  real(QP),allocatable :: nab(:)
  real(QP),allocatable :: vx(:,:),vc(:,:)
  real(QP) :: Ex,Ec
  real(QP) :: b(3,3)
  integer,allocatable :: LLL(:,:,:)
  integer :: Md, ML1,ML2,ML3
  integer :: SYStype=0

  real(8) :: mu=0.21951d0
  real(8) :: Kp=0.804d0
  real(8) :: beta=0.066725d0

CONTAINS


  SUBROUTINE calc_GGAPBE96_2( rgrid, rho, ene, pot, mu_in, Kp_in )

    implicit none

    type( grid ),intent(IN) :: rgrid
    type( GSArray ),intent(IN) :: rho
    type( xcene ) :: ene
    type( xcpot ) :: pot
    real(8),optional,intent(IN) :: mu_in, Kp_in

    type(gradient16) :: grad16
    type(fd) :: nabla
    type(lattice) :: aa, bb
    integer :: m1,m2,n1,n2,i,i1,i2,i3
    real(QP) :: Pi
    real(DP) :: sb(2),rb(2)

! ---

    if ( present(mu_in) ) mu=mu_in
    if ( present(Kp_in) ) Kp=Kp_in
    beta = mu*3.0d0/acos(-1.0d0)**2

    call construct_gradient16( rgrid, rho, grad16 )

    call construct_nabla_fd( nabla )

    Md = nabla%md

    if ( .not.allocated(nab) ) allocate( nab(-Md:Md) )
    nab(:) = nabla%coef(:)

    call get_aa_lattice( aa )
    call get_reciprocal_lattice( aa, bb )

    Pi=acos(-1.0q0)
    b(1:3,1)=aa%Length(1)*bb%LatticeVector(1:3,1)/( 2*Pi*rgrid%spacing(1) )
    b(1:3,2)=aa%Length(2)*bb%LatticeVector(1:3,2)/( 2*Pi*rgrid%spacing(2) )
    b(1:3,3)=aa%Length(3)*bb%LatticeVector(1:3,3)/( 2*Pi*rgrid%spacing(3) )

    ML1 = rgrid%g3%x%size_global
    ML2 = rgrid%g3%y%size_global
    ML3 = rgrid%g3%z%size_global

    call get_map_3d_to_1d_grid( rgrid, LLL )

! ---

    m1 = pot%xc%g_range%head
    m2 = pot%xc%g_range%tail
    n1 = pot%xc%s_range%head
    n2 = pot%xc%s_range%tail

    allocate( vx(m1:m2,n1:n2) ) ; vx=0.0q0

    call calc_PBE_x( rho, grad16 )

    allocate( vc(m1:m2,n1:n2) ) ; vc=0.0q0

    call calc_PBE_c( rho, grad16 )

! ---

    sb(1)=Ex*rgrid%VolumeElement
    sb(2)=Ec*rgrid%VolumeElement
    call MPI_ALLREDUCE( sb, rb, 2, MPI_REAL8, MPI_SUM, comm_grid, i )

    ene%Ex  = rb(1)
    ene%Ec  = rb(2)
    ene%Exc = rb(1)+rb(2)

    pot%xc%val(:,:)  = vx(:,:) + vc(:,:)
    if ( allocated(pot%x%val) ) pot%x%val(:,:) = vx(:,:)
    if ( allocated(pot%c%val) ) pot%c%val(:,:) = vc(:,:)

! ---

    deallocate( vc )
    deallocate( vx )
    deallocate( LLL )
    call destruct_gradient16( grad16 )

  END SUBROUTINE calc_GGAPBE96_2


  SUBROUTINE calc_PBE_x( rho, grad16 )
    implicit none
    type( GSArray ) :: rho
    type( gradient16 ) :: grad16
!   real(8),parameter :: mu=0.21951d0, Kp=0.804d0
    integer :: i,ispin,m,ierr
    real(QP) :: trho,cm
    real(QP),allocatable :: rrrr(:,:),rtmp(:)
    real(QP) :: kf, vx_lda, ex_lda, Fx, Pi, g2, factor
    real(QP) :: onethr,const1,const2
    integer :: i1,i2,i3,j,j1,j2,j3,k1,k2,k3
    integer :: mm,m1,m2

    Pi = acos(-1.0q0)

    factor = 1.0q0
    if ( rho%s_range%size_global == 2 ) factor = 2.0q0

    onethr = 1.0q0/3.0q0
    const1 = 3.0q0*Pi*Pi
    const2 = 3.0q0/(4.0q0*Pi)

    m1 = rho%g_range%head
    m2 = rho%g_range%tail
    mm = rho%g_range%size_global
    allocate( rtmp(m1:m2) ) ; rtmp=0.0q0
    allocate( rrrr(mm,3)  ) ; rrrr=0.0q0

    Ex = 0.0q0

    do ispin=rho%s_range%head,rho%s_range%tail

       do i=rho%g_range%head,rho%g_range%tail

          trho = factor*rho%val(i,ispin)

          if ( trho <= zero_density ) cycle

          kf = (const1*trho)**onethr

          ex_lda = -const2*kf
          vx_lda = -kf/Pi

          g2 = grad16%gg(i)

          Fx = 1.0q0 + Kp - 4.0q0*Kp*Kp*(trho*kf)**2 &
                            /( 4.0q0*Kp*(trho*kf)**2 + mu*g2 )

          Ex = Ex + trho*ex_lda*Fx

          vx(i,ispin) = vx(i,ispin) &
               + Fx*vx_lda + ( 24.0q0*Pi*Kp*Kp*mu*trho**3*g2 ) &
                            /( 4.0q0*Kp*(trho*kf)**2 + mu*g2 )**2

          rtmp(i) = -18.0q0*Pi*Kp*Kp*mu*trho**4 &
                    /( 4.0q0*Kp*(trho*kf)**2 + mu*g2 )**2

       end do ! i

       m1 = rho%g_range%head
       m2 = rho%g_range%tail
       rrrr(m1:m2,1) = rtmp(m1:m2)*grad16%gx(m1:m2)
       rrrr(m1:m2,2) = rtmp(m1:m2)*grad16%gy(m1:m2)
       rrrr(m1:m2,3) = rtmp(m1:m2)*grad16%gz(m1:m2)

       do i=1,3
          call mpi_allgatherv(rrrr(m1,i),ir_grid(myrank_g),mpi_real16 &
               ,rrrr(1,i),ir_grid,id_grid,mpi_real16,comm_grid,ierr)
       end do

       select case( SYStype )
       case default

          do i3=0,ML3-1
          do i2=0,ML2-1
          do i1=0,ML1-1
             i=LLL(i1,i2,i3)
             do m=-Md,Md
                cm=nab(m)*sign(1,m)
                j1=i1+m
                k1=j1/ML1 ; if ( j1<0 ) k1=(j1+1)/ML1-1
                j1=j1-k1*ML1
                j =LLL(j1,i2,i3)
! The potential vx is calculated at j-th grid point rather than i-th.
! This is because non-transposed nabla matrix Dij is used (See XC.doc).
                if ( rho%g_range%head <= j .and. j <= rho%g_range%tail ) then
                   vx(j,ispin) = vx(j,ispin) + cm*( rrrr(i,1)*b(1,1) &
                                                   +rrrr(i,2)*b(2,1) &
                                                   +rrrr(i,3)*b(3,1) )
                end if
                j2=i2+m
                k2=j2/ML2 ; if ( j2<0 ) k2=(j2+1)/ML2-1
                j2=j2-k2*ML2
                j =LLL(i1,j2,i3)
                if ( rho%g_range%head <= j .and. j <= rho%g_range%tail ) then
                   vx(j,ispin) = vx(j,ispin) + cm*( rrrr(i,1)*b(1,2) &
                                                   +rrrr(i,2)*b(2,2) &
                                                   +rrrr(i,3)*b(3,2) )
                end if
                j3=i3+m
                k3=j3/ML3 ; if ( j3<0 ) k3=(j3+1)/ML3-1
                j3=j3-k3*ML3
                j =LLL(i1,i2,j3)
                if ( rho%g_range%head <= j .and. j <= rho%g_range%tail ) then
                   vx(j,ispin) = vx(j,ispin) + cm*( rrrr(i,1)*b(1,3) &
                                                   +rrrr(i,2)*b(2,3) &
                                                   +rrrr(i,3)*b(3,3) )
                end if
             end do ! m
          end do ! i1
          end do ! i2
          end do ! i3

       end select

    end do ! ispin

    Ex = Ex/factor

    deallocate( rrrr )
    deallocate( rtmp )

  END SUBROUTINE calc_PBE_x


  SUBROUTINE calc_PBE_c( rho, grad16 )
    implicit none
    type( GSArray ) rho
    type( gradient16 ) :: grad16
    real(8),parameter :: A00  =0.031091d0,A01  =0.015545d0,A02  =0.016887d0
    real(8),parameter :: alp10=0.21370d0 ,alp11=0.20548d0 ,alp12=0.11125d0
    real(8),parameter :: bt10 =7.5957d0  ,bt11 =1.41189d1 ,bt12 =1.0357d1
    real(8),parameter :: bt20 =3.5876d0  ,bt21 =6.1977d0  ,bt22 =3.6231d0
    real(8),parameter :: bt30 =1.6382d0  ,bt31 =3.3662d0  ,bt32 =0.88026d0
    real(8),parameter :: bt40 =0.49294d0 ,bt41 =0.62517d0 ,bt42 =0.49671d0
!   real(8),parameter :: C1=2.14611945579107d0,C2=0.031091d0
    real(QP) :: C1,C2
    real(QP) :: const1, const2, factor
    integer :: i,j,ispin,m,i1,i2,i3,j1,j2,j3,k1,k2,k3,ierr
    integer :: m1,m2,mm,n1,n2
    real(QP) :: trho, rhoa, rhob
    real(QP) :: kf, rs, ec_U, ec_P, ec_lda, phi, g2
    real(QP) :: dac_dn, dfz_dz, deU_dn, deP_dn, H, A, T, alpc, fz
    real(QP) :: dec_dz, dphi_dz, dH_dA, dH_dT, tmp, dH_dphi, Ai
    real(QP) :: dA_dn, dec_dn, Hs, zeta, dT_dphi
    real(QP) :: dz_dn(2), cm, rssq
    real(QP) :: dac_drs, deP_drs, deU_drs, drs_dn
    real(QP),allocatable :: rrrr(:,:), rtmp(:)
    real(QP) :: Pi, one, two, fouthr, onethr
    real(QP) :: sevthr, twothr, thrtwo, ThrFouPi

    const1 = 2.0q0**(4.0q0/3.0q0)-2.0q0
    const2 = 9.0q0*(2.0q0**(1.0q0/3.0q0)-1.0q0)/4.0q0

    factor = 1.0q0
    if ( rho%s_range%size_global == 1 ) factor = 0.5q0
 
    Pi       = acos(-1.0q0)
    one      = 1.0q0
    two      = 2.0q0
    fouthr   = 4.0q0/3.0q0
    onethr   = 1.0q0/3.0q0
    ThrFouPi = 3.0q0/(4.0q0*Pi)
    thrtwo   = 3.0q0/2.0q0
    twothr   = 2.0q0/3.0q0
    sevthr   = 7.0q0/3.0q0

    C2 = ( 1.0q0-log(2.0q0) )/Pi**2 ! "gamma" in PBE paper
    C1 = beta/C2

    Ec = 0.0q0

    mm = rho%g_range%size_global

    allocate( rtmp(rho%g_range%head:rho%g_range%tail) ) ; rtmp=0.0q0
    allocate( rrrr(mm,3)  ) ; rrrr=0.0q0

    do i=rho%g_range%head,rho%g_range%tail

       rhoa = rho%val(i,1)*factor
       rhob = rho%val(i,rho%s_range%size_global)*factor
       trho = rhoa + rhob

       if ( trho <= zero_density ) cycle

       zeta = ( rhoa - rhob )/trho

       fz = ( (one+zeta)**fouthr + (one-zeta)**fouthr - two )*const1

       kf = ( 3.0q0*Pi*Pi*trho )**onethr

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

       phi = 0.5q0*( (one+zeta)**twothr + (one-zeta)**twothr )

       g2 = grad16%gg(i)

       T = Pi*g2/( 16.0q0*phi**2*kf*trho**2 )

       A = C1/( exp( -ec_lda/(C2*phi**3) ) - one )

       Hs = C2 * phi**3 * log( one + C1*T*(one+A*T)/(one+A*T+(A*T)**2) )

       Ec = Ec + trho*( ec_lda + Hs )

       drs_dn = -4.0q0*Pi*rs**4/9.0q0

       tmp = bt10*rssq + bt20*rs + bt30*rs*rssq + bt40*rs*rs
       deU_drs = alp10*ec_U/( one + alp10*rs ) &
            +A00*( one + alp10*rs )/rssq &
            *( bt10 + two*bt20*rssq + 3.0q0*bt30*rs + 4.0q0*bt40*rs*rssq ) &
            /( two*A00*tmp*tmp+tmp )

       tmp = bt11*rssq + bt21*rs + bt31*rs*rssq + bt41*rs*rs
       deP_drs = alp11*ec_P/( one + alp11*rs ) &
            +A01*( one + alp11*rs )/rssq &
            *( bt11 + two*bt21*rssq + 3.0q0*bt31*rs + 4.0q0*bt41*rs*rssq ) &
            /( two*A01*tmp*tmp+tmp )

       tmp = bt12*rssq + bt22*rs + bt32*rs*rssq + bt42*rs*rs
       dac_drs = alp12*alpc/( one + alp12*rs ) &
            +A02*( one + alp12*rs )/rssq &
            *( bt12 + two*bt22*rssq + 3.0q0*bt32*rs + 4.0q0*bt42*rs*rssq ) &
            /( two*A02*tmp*tmp+tmp )

       deU_dn = deU_drs * drs_dn
       deP_dn = deP_drs * drs_dn
       dac_dn = dac_drs * drs_dn


       dfz_dz = fouthr*( (one+zeta)**onethr - (one-zeta)**onethr )*const1

       dec_dz = -alpc*dfz_dz*const2*(one-zeta**4) &
               + 4.0q0*alpc*fz*const2*zeta**3 &
               +(ec_P-ec_U)*( dfz_dz*zeta**4 + fz*4.0q0*zeta**3 )

!       if ( zeta == 1.0q0 .or. zeta == -1.0q0 ) then
!          dphi_dz = 0.0q0
!       else
          dphi_dz = ( (one+zeta)**(-onethr)-(one-zeta)**(-onethr) )/3.0q0
!       end if

       tmp = one + A*T + (A*T)**2

       dH_dA = -phi**3*C1*C2*A*T**3*(two+A*T)/(tmp**2+C1*T*(one+A*T)*tmp)

       dH_dT =  phi**3*C1*C2*(one+two*A*T)/(tmp**2+C1*T*(one+A*T)*tmp)

       dH_dphi = 3.0q0*Hs/phi

       dz_dn(1)   = 2.0q0*rhob/trho**2
       dz_dn(rho%s_range%size_global) =-2.0q0*rhoa/trho**2

       do ispin=rho%s_range%head,rho%s_range%tail

          dec_dn = deU_dn - dac_dn*fz*const2*(one-zeta**4) &
               +(deP_dn-deU_dn)*fz*zeta**4 + dec_dz*dz_dn(ispin)

          dA_dn = A*(C1+A)/(C1*C2*phi**3) &
               *( dec_dn - 3.0q0*ec_lda/phi*dphi_dz*dz_dn(ispin) )

          vc(i,ispin) = vc(i,ispin) + ec_lda + Hs + trho*dec_dn &
               + trho*dH_dA*dA_dn &
               + dH_dT*(-sevthr*T - trho*two*T/phi*dphi_dz*dz_dn(ispin) ) &
               + trho*dH_dphi*dphi_dz*dz_dn(ispin)

       end do ! ispin

       rtmp(i) = dH_dT*Pi/(8.0q0*kf*trho*phi**2)

    end do ! i

    m1 = rho%g_range%head
    m2 = rho%g_range%tail
    rrrr(m1:m2,1) = rtmp(m1:m2)*grad16%gx(m1:m2)
    rrrr(m1:m2,2) = rtmp(m1:m2)*grad16%gy(m1:m2)
    rrrr(m1:m2,3) = rtmp(m1:m2)*grad16%gz(m1:m2)

    do i=1,3
       call mpi_allgatherv(rrrr(m1,i),ir_grid(myrank_g),mpi_real16 &
            ,rrrr(1,i),ir_grid,id_grid,mpi_real16,comm_grid,ierr)
    end do

    select case( SYStype )
    case default

       do i3=0,ML3-1
       do i2=0,ML2-1
       do i1=0,ML1-1
          i=LLL(i1,i2,i3)
          do m=-Md,Md
             cm=nab(m)*sign(1,m)
             j1=i1+m
             k1=j1/ML1 ; if ( j1<0 ) k1=(j1+1)/ML1-1
             j1=j1-k1*ML1
             j =LLL(j1,i2,i3)
! The potential vex is calculated at j-th grid point rather than i-th.
! This is because non-transposed nabla matrix Dij is used (See XC.doc).
             if ( rho%g_range%head <= j .and. j <= rho%g_range%tail ) then
                do ispin=rho%s_range%head,rho%s_range%tail
                   vc(j,ispin) = vc(j,ispin) + cm*( rrrr(i,1)*b(1,1) &
                                                   +rrrr(i,2)*b(2,1) &
                                                   +rrrr(i,3)*b(3,1) )
                end do
             end if
             j2=i2+m
             k2=j2/ML2 ; if ( j2<0 ) k2=(j2+1)/ML2-1
             j2=j2-k2*ML2
             j =LLL(i1,j2,i3)
             if ( rho%g_range%head <= j .and. j <= rho%g_range%tail ) then
                do ispin=rho%s_range%head,rho%s_range%tail
                   vc(j,ispin) = vc(j,ispin) + cm*( rrrr(i,1)*b(1,2) &
                                                   +rrrr(i,2)*b(2,2) &
                                                   +rrrr(i,3)*b(3,2) )
                end do
             end if
             j3=i3+m
             k3=j3/ML3 ; if ( j3<0 ) k3=(j3+1)/ML3-1
             j3=j3-k3*ML3
             j =LLL(i1,i2,j3)
             if ( rho%g_range%head <= j .and. j <= rho%g_range%tail ) then
                do ispin=rho%s_range%head,rho%s_range%tail
                   vc(j,ispin) = vc(j,ispin) + cm*( rrrr(i,1)*b(1,3) &
                                                   +rrrr(i,2)*b(2,3) &
                                                   +rrrr(i,3)*b(3,3) )
                end do
             end if
          end do ! m
       end do ! i1
      end do ! i2
      end do ! i3

    end select

    deallocate( rtmp )
    deallocate( rrrr )

    return
  END SUBROUTINE calc_PBE_c


END MODULE xc_ggapbe96_2_module
