MODULE ps_initrho_module

  use rgrid_module, only: dV,Ngrid,Igrid,Hgrid
  use ggrid_module
  use atom_module, only: Natom,Nelement, ki_atom,aa_atom
  use strfac_module, only: SGK
  use pseudopot_module, only: Mr,rad,rab,cdd,Zps,cdd_coef
  use electron_module, only: Nspin, dspin, Nelectron, Nelectron_spin
  use parallel_module
  use aa_module
  use fft_module
  use polint_module
  use random_initrho_module
  use psv_initrho_module
  use io_tools_module

  implicit none

  PRIVATE
  PUBLIC :: construct_ps_initrho
  PUBLIC :: read_ps_initrho

  logical :: flag_initrho_0
  logical,allocatable :: flag_initrho(:)
  real(8),allocatable :: cddg(:,:)

  complex(8),parameter :: z0=(0.0d0,0.0d0)

  integer :: iswitch_initrho=0
  logical :: analytic=.false.

CONTAINS


  SUBROUTINE read_ps_initrho
    implicit none
    call IOTools_readIntegerKeyword( "INITRHO", iswitch_initrho )
    if ( iswitch_initrho == 2 ) then
       call read_coef_psv_initrho !---> cdd_coef ( pseudopot_module )
    end if
  END SUBROUTINE read_ps_initrho


  SUBROUTINE construct_ps_initrho( rho )
    implicit none
    real(8),intent(OUT) :: rho(:,:)
    integer :: s

    call write_border( 80, " construct_ps_initrho(start)" )

! ---

    if ( allocated(cdd_coef) ) then
       if ( any(cdd_coef/=0.0d0) ) analytic=.true.
    end if

    if ( all(cdd==0.0d0) .and. iswitch_initrho /= 2 ) iswitch_initrho=3

! ---

    if ( disp_switch_parallel ) then
       write(*,*) "iswitch_initrho =",iswitch_initrho
       write(*,*) "analytic        =",analytic
    end if

    select case( iswitch_initrho )
    case default
       call construct_ps_initrho_g( rho )
    case( 1, 2 )
       call construct_ps_initrho_r_spin( rho )
    case( 3 )
       call construct_RandomInitrho( rho )
    end select

! ---

    call check_initrho( rho )

    if ( Nspin == 1 ) then
       call normalize_initrho( Nelectron, rho(:,1) )
    else if ( Nspin == 2 ) then
       do s=1,Nspin
          call normalize_initrho( Nelectron_spin(s), rho(:,s) )
       end do
    else
       stop "stop@construct_ps_initrho"
    end if

    call check_initrho( rho )

    call write_border( 80, " construct_ps_initrho(end)" )

  END SUBROUTINE construct_ps_initrho


  SUBROUTINE init_ps_initrho_g
    implicit none
    integer :: i,ig,ik,MKI
    real(8) :: sum0,G,x,sb,const,Vcell
    real(8),allocatable :: tmp(:)

    allocate( flag_initrho(Nelement) )
    flag_initrho(:)=.false.
    flag_initrho_0 =.false.
    if ( allocated(cdd) ) then
       do ik=1,Nelement
          if ( all(cdd(:,ik)==0.d0) ) cycle
          flag_initrho(ik) = .true.
          flag_initrho_0   = .true.
       end do
    end if

    if ( .not.flag_initrho_0 ) return

    MKI   = Nelement
    Vcell = Ngrid(0)*dV
    const = 1.d0/Vcell
    allocate( cddg(NMGL,MKI) ) ; cddg=0.d0
    do ik=1,MKI
       allocate( tmp(Mr(ik)) )
       do ig=1,NMGL
          G=sqrt(GG(ig))
          if ( G == 0.d0 ) then
             do i=1,Mr(ik)
                tmp(i)=cdd(i,ik)*rab(i,ik)
             end do
          else
             do i=1,Mr(ik)
                x=G*rad(i,ik)
                if ( x<1.d-1 ) then
                   sb=-(1.d0/39916800.d0*x**10-1.d0/362880.d0*x**8 &
                   +1.d0/5040.d0*x**6-1.d0/120.d0*x**4+1.d0/6.d0*x**2-1.d0)
                else
                   sb=sin(x)/x
                end if
                tmp(i)=cdd(i,ik)*sb*rab(i,ik)
             end do
          end if
          call simp(tmp,sum0,2)
          cddg(ig,ik)=sum0*const
       end do
       deallocate( tmp )
    end do ! ik
  END SUBROUTINE init_ps_initrho_g

  SUBROUTINE simp(f,s,m)
    implicit none
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


  SUBROUTINE construct_ps_initrho_g( rho )
    implicit none
    real(8),intent(OUT) :: rho(:,:)
    integer :: i,i1,i2,i3,ik,j,s,MG,ierr
    integer :: ML1,ML2,ML3,m_grid,m_spin
    real(8) :: c,c0
    real(8),allocatable :: rho_tmp(:,:),nelectron_ik(:)
    complex(8),allocatable :: zwork0(:,:,:),zwork1(:,:,:),vg(:)

    call init_ps_initrho_g

    if ( .not. flag_initrho_0 ) return

    MG   = NGgrid(0)
    ML1  = Ngrid(1)
    ML2  = Ngrid(2)
    ML3  = Ngrid(3)

    rho(:,:)=0.0d0

    allocate( zwork0(0:ML1-1,0:ML2-1,0:ML3-1) ) ; zwork0=z0
    allocate( zwork1(0:ML1-1,0:ML2-1,0:ML3-1) ) ; zwork1=z0
    allocate( vg(MG)                          ) ; vg=z0

    m_grid=size( rho, 1 )
    m_spin=size( rho, 2 )

    allocate( rho_tmp(m_grid,Nelement) ) ; rho_tmp=0.0d0

    allocate( nelectron_ik(Nelement) ) ; nelectron_ik=0.0d0

    do i=1,Natom
       ik=ki_atom(i)
       nelectron_ik(ik) = nelectron_ik(ik) + Zps(ik)
    end do

    call construct_Ggrid(2)

    call init_fft

    do ik=1,Nelement

       do i=MG_0,MG_1
          j=MGL(i)
          vg(i)=cddg(j,ik)*SGK(i,ik)
       end do

       call allgatherv_Ggrid( vg )

       zwork0(:,:,:)=z0
       do i=1,NGgrid(0)
          zwork0(LLG(1,i),LLG(2,i),LLG(3,i))=vg(i)
       end do

       call backward_fft( zwork0, zwork1 )

       i=0
       do i3=Igrid(1,3),Igrid(2,3)
       do i2=Igrid(1,2),Igrid(2,2)
       do i1=Igrid(1,1),Igrid(2,1)
          i=i+1
          rho_tmp(i,ik) = rho_tmp(i,ik) + real( zwork0(i1,i2,i3) )
       end do
       end do
       end do

       if ( minval(rho_tmp(:,ik)) < 0.0d0 ) then

          write(*,*) "WARNING: rho is negative at some points" &
               ,minval(rho_tmp(:,ik))

          !rho_tmp(:,ik)=abs( rho_tmp(:,ik) )
          where( rho_tmp(:,ik) < 0.0d0 )
             rho_tmp(:,ik)=0.0d0
          end where

       end if

       c=1.0d0/Nspin
       do s=1,m_spin
          rho(:,s) = rho(:,s) + c*rho_tmp(:,ik)
       end do

    end do ! ik

    call finalize_fft
    call destruct_Ggrid

    deallocate( nelectron_ik )
    deallocate( rho_tmp      )
    deallocate( vg           )
    deallocate( zwork1       )
    deallocate( zwork0       )

  END SUBROUTINE construct_ps_initrho_g


  SUBROUTINE construct_ps_initrho_r( rho )
    implicit none
    real(8),intent(OUT) :: rho(:,:)
    integer :: iatm,ielm,i,i1,i2,i3,j1,j2,j3,m_grid,m_spin,ierr,s
    real(8) :: a,a1,a2,a3,x,y,z,r1,r2,r3,c1,c2,c3,pi4,rr,zi,d
    real(8),allocatable :: rho_tmp(:)

    c1   = 1.0d0/Ngrid(1)
    c2   = 1.0d0/Ngrid(2)
    c3   = 1.0d0/Ngrid(3)
    pi4  = 4.0d0*acos(-1.0d0)

    m_grid = size( rho, 1 )
    m_spin = size( rho, 2 )

    rho(:,:)=0.0d0

    allocate( rho_tmp(m_grid) ) ; rho_tmp=0.0d0

    do ielm=1,Nelement

       rho_tmp(:)=0.0d0

       do iatm=1,Natom

          if ( ki_atom(iatm) /= ielm ) cycle

          a1=aa_atom(1,iatm)
          a2=aa_atom(2,iatm)
          a3=aa_atom(3,iatm)

          do j3=-1,1
          do j2=-1,1
          do j1=-1,1
             i=0
             do i3=Igrid(1,3),Igrid(2,3)
             do i2=Igrid(1,2),Igrid(2,2)
             do i1=Igrid(1,1),Igrid(2,1)
                i=i+1
                r1=i1*c1
                r2=i2*c2
                r3=i3*c3
                x=aa(1,1)*(r1-a1-j1)+aa(1,2)*(r2-a2-j2)+aa(1,3)*(r3-a3-j3)
                y=aa(2,1)*(r1-a1-j1)+aa(2,2)*(r2-a2-j2)+aa(2,3)*(r3-a3-j3)
                z=aa(3,1)*(r1-a1-j1)+aa(3,2)*(r2-a2-j2)+aa(3,3)*(r3-a3-j3)
                rr=x*x+y*y+z*z
                if ( analytic ) then
                   call use_analytic_function( rr,cdd_coef(:,:,ielm),d )
                else
                   call use_interpolation( sqrt(rr),rad(:,ielm),cdd(:,ielm),d )
                end if
                rho_tmp(i) = rho_tmp(i) + d
             end do
             end do
             end do
          end do
          end do
          end do

       end do ! iatm

       d=1.0d0/Nspin
       do s=1,m_spin
          rho(:,s) = rho(:,s) + d*rho_tmp(:)
       end do

    end do ! ielm

    deallocate( rho_tmp )

  END SUBROUTINE construct_ps_initrho_r

  SUBROUTINE use_analytic_function( rr, cdd_coef, d )
    implicit none
    real(8),intent(IN)  :: rr, cdd_coef(:,:)
    real(8),intent(OUT) :: d
    real(8) :: a,b,c
    integer :: i
    d=0.0d0
    do i=1,size(cdd_coef,2)
       a=cdd_coef(1,i)
       b=cdd_coef(2,i)
       c=cdd_coef(3,i)
       d=d+(a+b*rr)*exp(-c*rr)
    end do
  END SUBROUTINE use_analytic_function

  SUBROUTINE use_interpolation( r, rad, cdd, d )
    implicit none
    real(8),intent(IN)  :: r, rad(:), cdd(:)
    real(8),intent(OUT) :: d
    integer :: ir,mr
    real(8) :: rmin,rmax,pi4,err
    pi4  = 4.0d0*acos(-1.0d0)
    rmax = maxval( rad(:) )
    rmin = rad(2)
    mr   = size( rad )
    d    = 0.0d0
    if ( r > rmax ) return
    call get_nearest_index( rad, r, ir )
    if ( ir == 1 .or. r == 0.0d0 ) then
       d = cdd(2)/(pi4*rmin**2)
    else if ( ir < mr-1 ) then
       call polint( rad(ir-1:),cdd(ir-1:),3,r,d,err )
       d = d/(pi4*r*r)
    end if
  END SUBROUTINE use_interpolation

  SUBROUTINE get_nearest_index( rad, r, ir )
    implicit none
    real(8),intent(IN)  :: rad(:), r
    integer,intent(OUT) :: ir
    integer :: ir0,ir1
    ir0=1
    ir1=size( rad )
1   if ( ir1 - ir0 > 1 ) then
       ir = ( ir0 + ir1 )/2
       if ( rad(ir) > r ) then
          ir1=ir
       else
          ir0=ir
       end if
       goto 1
    end if
  END SUBROUTINE get_nearest_index


  SUBROUTINE construct_ps_initrho_r_spin( rho )
    implicit none
    real(8),intent(OUT) :: rho(:,:)
    integer :: iatm,ielm,i,i1,i2,i3,j1,j2,j3,m_grid,m_spin,ierr,s
    real(8) :: c(2),a1,a2,a3,x,y,z,r1,r2,r3,c1,c2,c3,pi4,rr,zi,d
    real(8),allocatable :: rho_tmp(:)

    c1   = 1.0d0/Ngrid(1)
    c2   = 1.0d0/Ngrid(2)
    c3   = 1.0d0/Ngrid(3)
    pi4  = 4.0d0*acos(-1.0d0)

    m_grid = size( rho, 1 )
    m_spin = size( rho, 2 )

    rho(:,:)=0.0d0

    allocate( rho_tmp(m_grid) ) ; rho_tmp=0.0d0

    rho_tmp(:)=0.0d0
    zi=0.0d0

    do iatm=1,Natom

       ielm=ki_atom(iatm)

       if ( allocated(dspin) ) then
          c(1) = 0.5d0*( Zps(ielm) + dspin(iatm) )/Zps(ielm)
          c(2) = 0.5d0*( Zps(ielm) - dspin(iatm) )/Zps(ielm)
       end if

       zi=zi+Zps(ielm)
       a1=aa_atom(1,iatm)
       a2=aa_atom(2,iatm)
       a3=aa_atom(3,iatm)

       do j3=-1,1
       do j2=-1,1
       do j1=-1,1
          i=0
          do i3=Igrid(1,3),Igrid(2,3)
          do i2=Igrid(1,2),Igrid(2,2)
          do i1=Igrid(1,1),Igrid(2,1)
             i=i+1
             r1=i1*c1
             r2=i2*c2
             r3=i3*c3
             x=aa(1,1)*(r1-a1-j1)+aa(1,2)*(r2-a2-j2)+aa(1,3)*(r3-a3-j3)
             y=aa(2,1)*(r1-a1-j1)+aa(2,2)*(r2-a2-j2)+aa(2,3)*(r3-a3-j3)
             z=aa(3,1)*(r1-a1-j1)+aa(3,2)*(r2-a2-j2)+aa(3,3)*(r3-a3-j3)
             rr=x*x+y*y+z*z
             if ( analytic ) then
                call use_analytic_function( rr,cdd_coef(:,:,ielm),d )
             else
                call use_interpolation( sqrt(rr),rad(:,ielm),cdd(:,ielm),d )
             end if
             rho_tmp(i) = rho_tmp(i) + d
             if ( allocated(dspin) ) then
                do s=1,Nspin
                   rho(i,s) = rho(i,s) + c(s)*d
                end do
             end if
          end do
          end do
          end do
       end do
       end do
       end do

    end do ! iatm

    if ( .not.allocated(dspin) ) then
       d=1.0d0/Nspin
       do s=1,Nspin
          rho(:,s) = d*rho_tmp(:)
       end do
    end if

    deallocate( rho_tmp )

  END SUBROUTINE construct_ps_initrho_r_spin


  SUBROUTINE normalize_initrho( N, rho )
    implicit none
    real(8),intent(IN) :: N
    real(8),intent(INOUT) :: rho(:)
    real(8) :: c0,c1
    integer :: i
    c0=sum( rho )*dV
    call MPI_ALLREDUCE( c0, c1, 1, MPI_REAL8, MPI_SUM, comm_grid, i )
    c1=N/c1
    rho(:) = c1*rho(:)
  END SUBROUTINE normalize_initrho


  SUBROUTINE check_initrho( rho )
    implicit none
    real(8),intent(IN) :: rho(:,:)
    integer :: m_spin,s,ierr
    real(8) :: a0,b0,c0,a,b,c
    m_spin=size( rho, 2 )
    if ( disp_switch_parallel ) then
       write(*,'(1x,a12,3a20)') "Nelectron","sum","min","max"
    end if
    do s=1,m_spin
       a0=sum( rho(:,s) )*dV
       b0=minval( rho(:,s) )
       c0=maxval( rho(:,s) )
       call MPI_ALLREDUCE(a0,a,1,MPI_REAL8,MPI_SUM,comm_grid,ierr)
       call MPI_ALLREDUCE(b0,b,1,MPI_REAL8,MPI_MIN,comm_grid,ierr)
       call MPI_ALLREDUCE(c0,c,1,MPI_REAL8,MPI_MAX,comm_grid,ierr)
       if ( disp_switch_parallel ) then
          if ( m_spin == 1 ) then
             write(*,'(1x,f12.5,2x,3f20.10)') Nelectron,a,b,c
          else if ( m_spin == 2 ) then
             write(*,'(1x,f12.5,2x,3f20.10)') Nelectron_spin(s),a,b,c
          end if
       end if
    end do
  END SUBROUTINE check_initrho


END MODULE ps_initrho_module
