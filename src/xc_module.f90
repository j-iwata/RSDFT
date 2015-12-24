MODULE xc_module

  use rgrid_module, only: Igrid,Ngrid,Hgrid,dV
  use density_module, only: rho, get_range_density, construct_density_v2
  use ps_pcc_module, only: flag_pcc_0, rhoc
  use parallel_module, only: comm_grid, mpi_real8, mpi_integer, mpi_sum
  use array_bound_module, only: MSP,MSP_0,MSP_1,MB_0,MB_1,MBZ_0,MBZ_1
  use kinetic_module, only: SYStype
  use kinetic_variables, only: Md

  use grid_module, only: grid, get_range_rgrid
  use xc_variables, only: xcpot, xcene

  use xc_ldapz81_module
  use xc_pw92_gth_module
  use xc_ggapbe96_module
  use xc_ggapbe96_2_module
  use xc_hybrid_module, only: iflag_hybrid,iflag_hse,alpha_hf
  use xc_hse_module
  use xc_hf_module
  use xc_vdw_module

  use basic_type_factory
  use basic_type_methods
  use io_tools_module

  implicit none

  PRIVATE
  PUBLIC :: calc_xc,Vxc,Exc,E_exchange,E_correlation &
           ,E_exchange_exx
  PUBLIC :: read_xc

  character(8),PUBLIC :: XCtype = 'LDAPZ81'

  real(8),allocatable :: Vxc(:,:)
  real(8) :: Exc,E_exchange,E_correlation
  real(8) :: E_exchange_exx

  real(8),allocatable :: Vx(:,:), Vc(:,:)

  logical :: disp_sw
  logical :: flag_read = .true.

CONTAINS


  SUBROUTINE read_xc
    implicit none
    call write_border( 0, " read_xc(start)" )
    call IOTools_readStringKeyword( "XCTYPE", XCtype )
    flag_read = .false.
    call write_border( 0, " read_xc(end)" )
  END SUBROUTINE read_xc


  SUBROUTINE chk_density(rho,rhoc)
    implicit none
    real(8),optional,intent(IN) :: rho(:,:),rhoc(:)
    real(8),allocatable :: trho(:,:),sb(:,:),rb(:,:)
    integer,allocatable :: ic(:,:),jc(:,:)
    integer :: i,j,m,n,ierr

    if ( present(rho) ) then
       m=size(rho,1)
       n=size(rho,2)
       allocate( trho(m,n) )
       trho=rho
    else if ( present(rhoc) ) then
       m=size(rhoc,1)
       n=1
       allocate( trho(m,n) )
       trho(:,1)=rhoc(:)
    end if

    allocate( sb(3,n) ) ; sb=0.0d0
    allocate( rb(3,n) ) ; rb=0.0d0
    allocate( ic(3,n) ) ; ic=0
    allocate( jc(3,n) ) ; jc=0

    sb(:,:)=0.0d0
    do j=1,n
       do i=1,m
          sb(1,j) = sb(1,j) + trho(i,j)
          if ( trho(i,j) > 0.0d0 ) sb(2,j) = sb(2,j) + trho(i,j)
          if ( trho(i,j) < 0.0d0 ) sb(3,j) = sb(3,j) + trho(i,j)
       end do
       ic(1,j)=count( trho(:,j) == 0.0d0 )
       ic(2,j)=count( trho(:,j) >  0.0d0 )
       ic(3,j)=count( trho(:,j) <  0.0d0 )
    end do

    sb=sb*dV
    call mpi_allreduce(sb,rb,3*n,mpi_real8,mpi_sum,comm_grid,ierr)
    call mpi_allreduce(ic,jc,3*n,mpi_integer,mpi_sum,comm_grid,ierr)

!    if ( disp_sw ) then
!       do j=1,n
!          write(*,'(1x,"(positive)",i8,2x,2g16.8)') jc(2,j),rb(2,j),rb(1,j)
!          write(*,'(1x,"(negative)",i8,2x,2g16.8)') jc(3,j),rb(3,j),rb(1,j)
!          if ( jc(1,j) /= 0 ) write(*,'(1x,"(zero    )",i8)') jc(1,j)
!       end do
!    end if

    deallocate( jc,ic )
    deallocate( rb,sb )

  END SUBROUTINE chk_density


  SUBROUTINE calc_xc

    implicit none

    real(8),allocatable :: rho_tmp(:,:)
    real(8) :: c,mu,kappa
    integer :: s,ML_0,ML_1,MSP_0,MSP_1,MSP
    type( GSarray_v2 ) :: density_v2
    type( GSarray ) :: density
    type( grid ) :: rg
    type( xcpot ) :: pot
    type( xcene ) :: ene

    call write_border( 1, " calc_xc(start)" )
    call check_disp_switch( disp_sw, 0 )

    if ( flag_read ) call read_xc

    call get_range_rgrid( rg )

    call get_range_density( density%g_range, density%s_range )
    call allocateGSArray( density )

    call construct_density_v2( density_v2 )

    call get_range_xc( pot )
    call allocateGSArray( pot%xc )

    ML_0 = pot%xc%g_range%head
    ML_1 = pot%xc%g_range%tail
    MSP_0 = pot%xc%s_range%head
    MSP_1 = pot%xc%s_range%tail
    MSP  = pot%xc%s_range%size_global

    if ( .not.allocated(Vxc) ) then
       allocate( Vxc(ML_0:ML_1,MSP_0:MSP_1) )
       Vxc=0.0d0
    end if

    allocate( rho_tmp(ML_0:ML_1,MSP) )
    rho_tmp(:,:)=rho(:,:)

    if ( flag_pcc_0 ) then
       c=1.0d0
       if ( MSP == 2 ) c=0.5d0
       do s=1,MSP
          rho_tmp(:,s) = rho_tmp(:,s) + c*rhoc(:)
       end do
       call chk_density(rho)
       call chk_density(rhoc=rhoc)
    end if
    call chk_density(rho_tmp)

    density%val(:,:) = rho_tmp(:,:)
    density_v2%val(:,:) = rho_tmp(:,:)

    select case(XCtype)
    case('LDAPZ81')

       call calc_LDAPZ81( density_v2, ene, pot )

       E_exchange    = ene%Ex
       E_correlation = ene%Ec
       Exc           = ene%Exc
       Vxc(:,:)      = pot%xc%val(:,:)

    case('LDAPW92')

       if ( flag_pcc_0 ) stop "PCC is not implemented in LDAPW92" 
       call calc_pw92_gth(ML_0,ML_1,MSP,MSP_0,MSP_1,rho,Exc,Vxc,dV,comm_grid)

    case('GGAPBE96')

       call init_GGAPBE96( Igrid, MSP_0, MSP_1, MSP, comm_grid, dV &
            ,Md, Hgrid, Ngrid, SYStype )

       call calc_GGAPBE96( rho_tmp, Exc, Vxc, E_exchange, E_correlation )

    case('PBE','PBE96')

       call calc_GGAPBE96_2( rg, density, ene, pot )

       E_exchange    = ene%Ex
       E_correlation = ene%Ec
       Exc           = ene%Exc
       Vxc(:,:)      = pot%xc%val(:,:)

    case('PBEsol')

       mu = 10.0d0/81.0d0

       call init_GGAPBE96( Igrid, MSP_0, MSP_1, MSP, comm_grid, dV &
            ,Md, Hgrid, Ngrid, SYStype, mu_in=mu )

       call calc_GGAPBE96( rho_tmp, Exc, Vxc, E_exchange, E_correlation )

    case('revPBE')

       kappa = 1.245d0

       call init_GGAPBE96( Igrid, MSP_0, MSP_1, MSP, comm_grid, dV &
            ,Md, Hgrid, Ngrid, SYStype, Kp_in=kappa )

       call calc_GGAPBE96( rho_tmp, Exc, Vxc, E_exchange, E_correlation )

    case('HF')

       if ( iflag_hybrid == 0 ) then

          if ( disp_sw ) then
             write(*,*) "XCtype, iflag_hybrid =",XCtype, iflag_hybrid
             write(*,*) "LDAPZ81 is called (iflag_hybrid==0)"
          end if

          call calc_LDAPZ81( density_v2, ene, pot )

          E_exchange    = ene%Ex
          E_correlation = ene%Ec
          Exc           = ene%Exc
          Vxc(:,:)      = pot%xc%val(:,:)

       else

          call init_xc_hf( ML_0,ML_1, MSP_0,MSP_1, MBZ_0,MBZ_1 &
                          ,MB_0,MB_1, SYStype, dV )
          call calc_xc_hf( E_exchange_exx )

          Vxc(:,:)      = 0.0d0
          E_exchange    = 0.0d0
          E_correlation = 0.0d0
          Exc           = E_exchange_exx

       end if

    case('PBE0')

       call init_GGAPBE96( Igrid, MSP_0, MSP_1, MSP, comm_grid, dV &
            ,Md, Hgrid, Ngrid, SYStype )

       if ( .not.allocated(Vx) ) then
          allocate( Vx(ML_0:ML_1,MSP_0:MSP_1) ) ; Vx=0.0d0
       end if

       call calc_GGAPBE96( rho_tmp, Exc, Vxc, E_exchange, E_correlation, Vx )

       if ( iflag_hybrid /= 0 ) then

          call init_xc_hf( ML_0,ML_1, MSP_0,MSP_1, MBZ_0,MBZ_1 &
                          ,MB_0,MB_1, SYStype, dV )

          call calc_xc_hf( E_exchange_exx )

          Vxc(:,:) = Vxc(:,:) - alpha_hf*Vx(:,:)

          E_exchange = (1.0d0-alpha_hf)*E_exchange + E_exchange_exx

          Exc = E_exchange + E_correlation

       end if

    case('HSE','HSE06')

       if ( iflag_hybrid == 0 ) then

          if ( disp_sw ) then
             write(*,*) "XCtype, iflag_hybrid =",XCtype, iflag_hybrid
             write(*,*) "GGAPBE96 is called (iflag_hybrid==0)"
          end if

          call init_GGAPBE96( Igrid, MSP_0, MSP_1, MSP, comm_grid &
               , dV ,Md, Hgrid, Ngrid, SYStype )
          call calc_GGAPBE96( rho_tmp, Exc, Vxc, E_exchange, E_correlation )

       else

          call init_xc_hse( ML_0, ML_1, MSP,MSP_0,MSP_1 &
               ,MBZ_0,MBZ_1,MB_0,MB_1,SYStype )

          call calc_xc_hse( iflag_hse, rho, Exc &
               , Vxc, E_exchange, E_correlation )

          call init_xc_hf( ML_0,ML_1, MSP_0,MSP_1, MBZ_0,MBZ_1 &
                          ,MB_0,MB_1, SYStype, dV )

          call calc_xc_hf( E_exchange_exx )

          E_exchange = E_exchange + E_exchange_exx

          Exc = E_exchange + E_correlation

       end if

    case('LCwPBE')

       stop "LCwPBE is not available yet"

    case('VDWDF')

       call allocateGSArray( pot%x )
       call allocateGSArray( pot%c )

       kappa = 1.245d0 ! revPBE

       call calc_GGAPBE96_2( rg, density, ene, pot, Kp_in=kappa )

       E_exchange = ene%Ex

       call init_xc_vdw
       call calc_xc_vdw( rg, density, ene, pot )

       E_correlation = ene%Ec
       Vxc(:,:) = pot%x%val(:,:) + pot%c%val(:,:)

       Exc = E_exchange + E_correlation

    case default

       write(*,*) "Invalid Keyword: XCtype=",XCtype
       stop "stop@calc_xc"

    end select

    deallocate( rho_tmp )

    call write_border( 1, " calc_xc(end)" )

  END SUBROUTINE calc_xc


  SUBROUTINE get_range_xc( pot )
    implicit none
    type(xcpot) :: pot
    pot%xc%g_range%head = Igrid(1,0)
    pot%xc%g_range%tail = Igrid(2,0)
    pot%xc%g_range%size = Igrid(2,0)-Igrid(1,0)+1
    pot%xc%s_range%head = MSP_0
    pot%xc%s_range%tail = MSP_1
    pot%xc%s_range%size = MSP_1-MSP_0+1
    pot%xc%g_range%head_global = 1
    pot%xc%g_range%tail_global = Ngrid(0)
    pot%xc%g_range%size_global = Ngrid(0)
    pot%xc%s_range%head_global = 1
    pot%xc%s_range%tail_global = MSP
    pot%xc%s_range%size_global = MSP
    pot%x%g_range = pot%xc%g_range
    pot%x%s_range = pot%xc%s_range
    pot%c%g_range = pot%xc%g_range
    pot%c%s_range = pot%xc%s_range
  END SUBROUTINE get_range_xc


END MODULE xc_module
