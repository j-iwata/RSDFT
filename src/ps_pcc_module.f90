MODULE ps_pcc_module

  use rgrid_module, only: Ngrid,Igrid,dV
  use ggrid_module
  use atom_module, only: Nelement
  use strfac_module, only: SGK
  use pseudopot_module, only: Mr,rad,rab,cdc
  use fft_module

  implicit none

  PRIVATE
  PUBLIC :: rhoc,flag_pcc_0,init_ps_pcc,construct_ps_pcc,cdcg

  real(8),allocatable :: rhoc(:)
  logical :: flag_pcc_0
  logical,allocatable :: flag_pcc(:)
  real(8),allocatable :: cdcg(:,:)

CONTAINS


  SUBROUTINE init_ps_pcc
    implicit none
    integer :: i,ig,ik,MKI
    real(8) :: const,sum0,sb,x,G,Vcell
    real(8),allocatable :: tmp(:)

    call write_border( 80, " init_ps_pcc(start)" )

    allocate( flag_pcc(Nelement) )
    flag_pcc(:) = .false.
    flag_pcc_0  = .false.
    if ( allocated(cdc) ) then
       do ik=1,Nelement
          if ( all(cdc(:,ik)==0.d0) ) cycle
          flag_pcc(ik)=.true.
          flag_pcc_0  =.true.
       end do
    end if

    if ( .not.flag_pcc_0 ) return

    MKI   = Nelement
    Vcell = Ngrid(0)*dV
    const = 4.d0*acos(-1.d0)/Vcell

    if ( any( flag_pcc(:) ) ) then
       allocate( cdcg(NMGL,MKI) ) ; cdcg=0.d0
       do ik=1,MKI
          if ( .not.flag_pcc(ik) ) cycle
          allocate( tmp(Mr(ik)) )
          do ig=1,NMGL
             G=sqrt(GG(ig))
             if ( G == 0.d0 ) then
                do i=1,Mr(ik)
                   tmp(i)=rad(i,ik)*rad(i,ik)*cdc(i,ik)*rab(i,ik)
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
                   tmp(i)=rad(i,ik)*rad(i,ik)*cdc(i,ik)*sb*rab(i,ik)
                end do
             end if
             call simp(tmp,sum0,2)
             cdcg(ig,ik)=sum0*const
          end do
          deallocate( tmp )
       end do !ik
    end if

    call write_border( 80, " init_ps_local(end)" )

  END SUBROUTINE init_ps_pcc

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


  SUBROUTINE construct_ps_pcc
    implicit none
    integer :: a,i,i1,i2,i3,ik,j,MG
    integer :: ML1,ML2,ML3,ML,ML_0,ML_1
    complex(8),allocatable :: zwork(:,:,:),zwork1(:,:,:),vg(:)
    complex(8),parameter :: z0=(0.0d0,0.0d0)

    if ( .not.flag_pcc_0 ) return

    call write_border( 80, " construct_ps_pcc(start)" )

    MG   = NGgrid(0)
    ML   = Ngrid(0)
    ML1  = Ngrid(1)
    ML2  = Ngrid(2)
    ML3  = Ngrid(3)
    ML_0 = Igrid(1,0)
    ML_1 = Igrid(2,0)

    if ( .not.allocated(rhoc) ) then
       allocate( rhoc(ML_0:ML_1) )
       rhoc=0.0d0
    end if

    allocate( zwork(0:ML1-1,0:ML2-1,0:ML3-1) ) ; zwork=z0

    allocate( vg(MG) ) ; vg=z0

    do i=MG_0,MG_1
       j=MGL(i)
       vg(i)=cdcg(j,1)*SGK(i,1)
    end do
    do ik=2,Nelement
       do i=MG_0,MG_1
          j=MGL(i)
          vg(i)=vg(i)+cdcg(j,ik)*SGK(i,ik)
       end do
    end do
    call allgatherv_Ggrid(vg)

    call construct_Ggrid(2)

    zwork(:,:,:)=z0
    do i=1,NGgrid(0)
       zwork(LLG(1,i),LLG(2,i),LLG(3,i))=vg(i)
    end do

    call destruct_Ggrid

    deallocate( vg )

    call init_fft
    call backward_fft( zwork, zwork1 )

!    rhoc(:)=0.d0
!    i=ML_0-1
!    do i3=Igrid(1,3),Igrid(2,3)
!    do i2=Igrid(1,2),Igrid(2,2)
!    do i1=Igrid(1,1),Igrid(2,1)
!       i=i+1
!       rhoc(i)=rhoc(i)+real( zwork(i1,i2,i3) )
!    end do
!    end do
!    end do
    call z3_to_d1_fft( zwork, rhoc ) 

    call finalize_fft

    if ( allocated(zwork1) ) deallocate( zwork1 )
    deallocate( zwork )

    call write_border( 80, " construct_ps_pcc(end)" )

  END SUBROUTINE construct_ps_pcc


END MODULE ps_pcc_module
