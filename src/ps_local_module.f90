MODULE ps_local_module

  use rgrid_module
  use ggrid_module
  use atom_module
  use strfac_module
  use bb_module
  use density_module
  use electron_module
  use pseudopot_module
  use parallel_module
  use watch_module
  use ps_local_gth_module
  use simc_module
  use ffte_sub_module
  use bberf_module
  use ps_local_variables
  use fft_module
  use var_ps_member
  use ps_local_fftw_module, only: construct_ps_local_fftw

  implicit none

  PRIVATE
  PUBLIC :: Vion
  PUBLIC :: init_ps_local
  PUBLIC :: construct_ps_local
  PUBLIC :: construct_ps_local_ffte
  real(8),PUBLIC :: const_ps_local

!  real(8),allocatable :: vqlgl(:,:),vqls(:,:)
  real(8),allocatable :: Vion(:)

  logical :: first_time1=.true.
  integer :: NGHT
  integer,allocatable :: LGHT(:,:),IGHT(:)
  complex(8),allocatable :: fg(:)
  logical :: flag_zero_ave = .false.

CONTAINS


  SUBROUTINE init_ps_local_test( GG, pr, pl )
    implicit none
    real(8),intent(IN) :: GG(:)
    type(ps1d),intent(INOUT) :: pr(:)
    type(pslocal),intent(INOUT) :: pl
    integer :: i,ig,ik,iorb,MMr,NRc,MKI,NMGL
    real(8) :: Rc,p1,p2,p3,p4,vlong,Pi,const,x,r,sb,sum0,G,G2
    real(8) :: Vcell
    real(8),allocatable :: vshort(:),tmp(:)

    MKI   = size(pr)
    NMGL  = size(GG)
    Vcell = Ngrid(0)*dV
    Pi    = acos(-1.0d0)
    const = 4.0d0*Pi/Vcell

    allocate( pl%vqlg(NMGL,MKI) ) ; pl%vqlg=0.0d0

    MMr=0
    do ik=1,MKI
       MMr = max( pr(ik)%Mr, MMr )
    end do

    allocate( vshort(MMr) ) ; vshort=0.0d0
    allocate( tmp(MMr)    ) ; tmp=0.0d0

    do ik=1,MKI

       if ( pr(ik)%ippform == 4 ) then
          call init_ps_local_gth( Vcell, NMGL, ik, GG, pl%vqlg(:,ik) )
          cycle
       end if

       MMr = pr(ik)%Mr

       Rc=0.d0
       NRc=0
       do iorb=1,pr(ik)%norb
          Rc=max( Rc, pr(ik)%Rps(iorb) )
          NRc=max( NRc, pr(ik)%NRps(iorb) )
       end do

       if ( Rc < 1.d-8 ) Rc=5.d0
       if ( NRc == 0 ) then
          do i=1,MMr
             if ( pr(ik)%rad(i) > Rc ) then
                NRc=i
                exit
             end if
          end do
       end if

       call simc(pr(ik)%rad,pr(ik)%vql,Rc,pr(ik)%Zps,pr(ik)%parloc,MMr)

       p1=pr(ik)%parloc(1) ; p2=sqrt(pr(ik)%parloc(2))
       p3=pr(ik)%parloc(3) ; p4=sqrt(pr(ik)%parloc(4))

       do i=1,MMr
          r=pr(ik)%rad(i)
          if ( r < 1.d-9 ) then
             vlong=-2.d0*pr(ik)%Zps/sqrt(Pi)*(p1*p2+p3*p4)
          else
             vlong=-pr(ik)%Zps/r*( p1*bberf(p2*r)+p3*bberf(p4*r) )
          end if
          vshort(i)=pr(ik)%vql(i)-vlong
       end do

       do ig=1,NMGL
          G=sqrt(GG(ig))
          if ( G == 0.d0 ) then
             do i=1,MMr
                tmp(i)=(pr(ik)%rad(i))**2*vshort(i)*pr(ik)%rab(i)
             end do
          else
             do i=1,MMr
                x=G*pr(ik)%rad(i)
                if ( x < 1.d-1 ) then
                   sb=-(1.d0/39916800.d0*x**10-1.d0/362880.d0*x**8 &
                       +1.d0/5040.d0*x**6-1.d0/120.d0*x**4+1.d0/6.d0*x**2-1.d0)
                else
                   sb=sin(x)/x
                end if
                tmp(i)=(pr(ik)%rad(i))**2*vshort(i)*sb*pr(ik)%rab(i)
             end do
          end if
          call simp(tmp(1:MMr),sum0,2)
          pl%vqlg(ig,ik)=sum0*const
       end do

       p1=-pr(ik)%Zps*pr(ik)%parloc(1) ; p2=0.25d0/pr(ik)%parloc(2)
       p3=-pr(ik)%Zps*pr(ik)%parloc(3) ; p4=0.25d0/pr(ik)%parloc(4)
       do ig=1,NMGL
          G2=GG(ig)
          if ( G2 == 0.0d0 ) then
             pl%vqlg(ig,ik) = pl%vqlg(ig,ik)-(p1*p2+p3*p4)*const
          else
             pl%vqlg(ig,ik) = pl%vqlg(ig,ik) &
                  +(p1*exp(-G2*p2)+p3*exp(-G2*p4))/G2*const
          end if
       end do

    end do ! ik

    deallocate( tmp    )
    deallocate( vshort )

! --- const_ps_local

   pl%const_ps_local=0.0d0

   if ( flag_zero_ave ) then

      do ig=1,NMGL
        if ( GG(ig) == 0.0d0 ) exit
      end do
      do i=1,Natom
         ik=ki_atom(i)
         pl%const_ps_local = pl%const_ps_local + pl%vqlg(ig,ik)    
      end do
      pl%vqlg(ig,ik)=0.0d0

   end if

  END SUBROUTINE init_ps_local_test


  SUBROUTINE init_ps_local
    implicit none
    integer :: i,ig,ik,iorb,MMr,NRc,MKI
    real(8) :: Rc,p1,p2,p3,p4,vlong,Pi,const,x,r,sb,sum0,G,G2
    real(8) :: Vcell
    real(8),allocatable :: vshort(:),tmp(:)

    call write_border( 80, " init_ps_local(start)" )

    !call init_ps_local_test( GG, ps, psloc )

    MKI   = Nelement
    Vcell = Ngrid(0)*dV
    Pi    = acos(-1.d0)
    const = 4.d0*Pi/Vcell

    allocate( vqlg(NMGL,MKI)  ) ; vqlg=0.0d0

    MMr=maxval(Mr)
!    allocate( vqls(MMr,MKI)   ) ; vqls=0.d0
!    allocate( vqlgl(NMGL,MKI) ) ; vqlgl=0.d0

    allocate( vshort(MMr) ) ; vshort=0.0d0

    do ik=1,MKI

       if ( ippform(ik) == 4 ) then
          call init_ps_local_gth( Vcell, NMGL, ik, GG, vqlg(1,ik) )
          cycle
       end if

       MMr = Mr(ik)

       Rc=0.d0
       NRc=0
       do iorb=1,norb(ik)
          Rc=max( Rc, Rps(iorb,ik) )
          NRc=max( NRc, NRps(iorb,ik) )
       end do

       if ( Rc<1.d-8 ) Rc=5.d0
       if ( NRc==0 ) then
          do i=1,MMr
             if ( rad(i,ik)>Rc ) then
                NRc=i
                exit
             end if
          end do
       end if

       call simc(rad(1,ik),vql(1,ik),Rc,Zps(ik),parloc(1,ik),MMr)

       p1=parloc(1,ik) ; p2=sqrt(parloc(2,ik))
       p3=parloc(3,ik) ; p4=sqrt(parloc(4,ik))

       do i=1,MMr
          r=rad(i,ik)
          if ( r<1.d-9 ) then
             vlong=-2.d0*Zps(ik)/sqrt(Pi)*(p1*p2+p3*p4)
          else
             vlong=-Zps(ik)/r*( p1*bberf(p2*r)+p3*bberf(p4*r) )
          end if
          vshort(i)=vql(i,ik)-vlong
!          vqls(i,ik)=vql(i,ik)-vlong
       end do

       allocate( tmp(MMr) )

       do ig=1,NMGL
          G=sqrt(GG(ig))
          if ( G == 0.d0 ) then
             do i=1,MMr
                tmp(i)=rad(i,ik)*rad(i,ik)*vshort(i)*rab(i,ik)
             end do
          else
             do i=1,MMr
                x=G*rad(i,ik)
                if ( x<1.d-1 ) then
                   sb=-(1.d0/39916800.d0*x**10-1.d0/362880.d0*x**8 &
                       +1.d0/5040.d0*x**6-1.d0/120.d0*x**4+1.d0/6.d0*x**2-1.d0)
                else
                   sb=sin(x)/x
                end if
                tmp(i)=rad(i,ik)*rad(i,ik)*vshort(i)*sb*rab(i,ik)
             end do
          end if
          call simp(tmp,sum0,2)
          vqlg(ig,ik)=sum0*const
       end do

       p1=-Zps(ik)*parloc(1,ik) ; p2=0.25d0/parloc(2,ik)
       p3=-Zps(ik)*parloc(3,ik) ; p4=0.25d0/parloc(4,ik)
       do ig=1,NMGL
          G2=GG(ig)
          if ( G2 == 0.d0 ) then
!             vqlgl(ig,ik)=-(p1*p2+p3*p4)*const
!             vqlg(ig,ik)=vqlg(ig,ik)+vqlgl(ig,ik)
             vqlg(ig,ik)=vqlg(ig,ik)-(p1*p2+p3*p4)*const
          else
!             vqlgl(ig,ik)=(p1*exp(-G2*p2)+p3*exp(-G2*p4))/G2*const
!             vqlg(ig,ik)=vqlg(ig,ik)+vqlgl(ig,ik)
             vqlg(ig,ik)=vqlg(ig,ik)+(p1*exp(-G2*p2)+p3*exp(-G2*p4))/G2*const
          end if
       end do
       deallocate( tmp )

    end do ! ik

    deallocate( vshort )

! --- const_ps_local

   const_ps_local=0.0d0

   if ( flag_zero_ave ) then

      do ig=1,NMGL
        if ( GG(ig) == 0.0d0 ) exit
      end do
      do i=1,Natom
         ik=ki_atom(i)
         const_ps_local=const_ps_local+vqlg(ig,ik)    
      end do
      vqlg(ig,:)=0.0d0

   end if

    call write_border( 80, " init_ps_local(end)" )

  END SUBROUTINE init_ps_local


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


  SUBROUTINE construct_ps_local_test( SGK, pl )
    implicit none
    complex(8),intent(IN) :: SGK(:,:)
    type(pslocal),intent(INOUT) :: pl
    integer :: i,ik,j,MG,ML1,ML2,ML3,ML_0,ML_1,Nelement
    complex(8),allocatable :: zwork0(:,:,:),zwork1(:,:,:),vg(:)
    complex(8),parameter :: z0=(0.0d0,0.0d0)
    real(8) :: ctt(0:3),ett(0:3)
    logical :: disp_sw

    Nelement = size( pl%vqlg, 2 )

    MG  = NGgrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)
    ML_0= Igrid(1,0)
    ML_1= Igrid(2,0)

    ctt(:)=0.d0
    ett(:)=0.d0

    call watch(ctt(0),ett(0))

    if ( .not.allocated(pl%Vion) ) then
       allocate( pl%Vion(ML_0:ML_1) )
       pl%Vion=0.0d0
    end if

    allocate( vg(MG) )

    do i=MG_0,MG_1
       j=MGL(i)
       vg(i) = pl%vqlg(j,1)*SGK(i,1)
    end do
    do ik=2,Nelement
       do i=MG_0,MG_1
          j=MGL(i)
          vg(i) = vg(i) + pl%vqlg(j,ik)*SGK(i,ik)
       end do
    end do
    call allgatherv_Ggrid(vg)

    call construct_Ggrid(2)

    allocate( zwork0(0:ML1-1,0:ML2-1,0:ML3-1) ) ; zwork0=z0

    do i=1,NGgrid(0)
       zwork0(LLG(1,i),LLG(2,i),LLG(3,i))=vg(i)
    end do

    call destruct_Ggrid

    deallocate( vg )

    call init_fft

    call watch(ctt(1),ett(1))

    call backward_fft( zwork0, zwork1 )

    call watch(ctt(2),ett(2))

    call z3_to_d1_fft( zwork0, pl%Vion )

    call finalize_fft

    if ( allocated(zwork1) ) deallocate( zwork1 )
    if ( allocated(zwork0) ) deallocate( zwork0 )

    call watch(ctt(3),ett(3))

    call check_disp_switch( disp_sw, 0 )
    if ( disp_sw ) then
       write(*,*) "time(const_ps_loc_t1)",ctt(1)-ctt(0),ett(1)-ett(0)
       write(*,*) "time(const_ps_loc_t2)",ctt(2)-ctt(1),ett(2)-ett(1)
       write(*,*) "time(const_ps_loc_t3)",ctt(3)-ctt(2),ett(3)-ett(2)
    end if

  END SUBROUTINE construct_ps_local_test


  SUBROUTINE construct_ps_local
    implicit none
    integer :: a,i,i1,i2,i3,ik,j,MG
    integer :: ML1,ML2,ML3,ML,ML_0,ML_1
    complex(8),allocatable :: zwork0(:,:,:),vg(:)
    complex(8),allocatable :: zwork1(:,:,:)
    real(8) :: ctt(0:3),ett(0:3)
    logical :: disp_sw

    !call construct_ps_local_test( SGK, psloc )

    call write_border( 0, " construct_ps_local(start)" )

    MG  = NGgrid(0)
    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)
    ML_0= Igrid(1,0)
    ML_1= Igrid(2,0)

    ctt(:)=0.d0
    ett(:)=0.d0

    call watch(ctt(0),ett(0))

    if ( .not.allocated(Vion) ) then
       allocate( Vion(ML_0:ML_1) )
       Vion=0.0d0
    end if

#ifdef _FFTE_
    call construct_ps_local_ffte
    return
#elif _FFTW_
    call construct_ps_local_fftw( vqlg, SGK, Vion )
    return
#endif

    allocate( vg(MG) )

    do i=MG_0,MG_1
       j=MGL(i)
       vg(i)=vqlg(j,1)*SGK(i,1)
    end do
    do ik=2,Nelement
       do i=MG_0,MG_1
          j=MGL(i)
          vg(i)=vg(i)+vqlg(j,ik)*SGK(i,ik)
       end do
    end do
    call allgatherv_Ggrid(vg)

    call construct_Ggrid(2)

    allocate( zwork0(0:ML1-1,0:ML2-1,0:ML3-1) )
    zwork0(:,:,:)=(0.d0,0.d0)

    do i=1,NGgrid(0)
       zwork0(LLG(1,i),LLG(2,i),LLG(3,i))=vg(i)
    end do

    call destruct_Ggrid

    deallocate( vg )

    call init_fft

    call watch(ctt(1),ett(1))

    call backward_fft( zwork0, zwork1 )

    call watch(ctt(2),ett(2))

    call z3_to_d1_fft( zwork0, Vion )

    call finalize_fft

    if ( allocated(zwork1) ) deallocate( zwork1 )
    if ( allocated(zwork0) ) deallocate( zwork0 )

    call watch(ctt(3),ett(3))

!    call check_disp_switch( disp_sw, 0 )
!    if ( disp_sw ) then
!       write(*,*) "time(const_ps_loc_1)",ctt(1)-ctt(0),ett(1)-ett(0)
!       write(*,*) "time(const_ps_loc_2)",ctt(2)-ctt(1),ett(2)-ett(1)
!       write(*,*) "time(const_ps_loc_3)",ctt(3)-ctt(2),ett(3)-ett(2)
!    end if

    call write_border( 0, " construct_ps_local(end)" )

  END SUBROUTINE construct_ps_local


  SUBROUTINE construct_ps_local_ffte
    implicit none
    integer :: i,i1,i2,i3,ik,j,n
    integer :: ML,ML1,ML2,ML3,MG,ML_0,ML_1
    real(8) :: ctt(0:3),ett(0:3)
    integer :: a1b,b1b,a2b,b2b,a3b,b3b,ab1,ab12

    call write_border( 80, " construct_ps_local_ffte(start)" )

    MG  = NGgrid(0)
    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)
    ML_0= Igrid(1,0)
    ML_1= Igrid(2,0)
    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)
    ab1 = b1b-a1b+1
    ab12= (b2b-a2b+1)*(b1b-a1b+1)

    if ( first_time1 ) then
       if ( .not.allocated(Vion) ) then
          allocate( Vion(ML_0:ML_1) )
       end if
       call construct_Ggrid(0)
       n=0
       do i=1,NGgrid(0)
          i1=mod( Ngrid(1)+LLG(1,i), Ngrid(1) )
          i2=mod( Ngrid(2)+LLG(2,i), Ngrid(2) )
          i3=mod( Ngrid(3)+LLG(3,i), Ngrid(3) )
          if ( a2b <= i2 .and. i2 <= b2b .and. a3b <= i3 .and. i3 <= b3b ) then
             n=n+1
          end if
       end do
       allocate( LGHT(3,n) ) ; LGHT=0
       allocate( IGHT(n)   ) ; IGHT=0
       n=0
       do i=1,NGgrid(0)
          i1=mod( Ngrid(1)+LLG(1,i), Ngrid(1) )
          i2=mod( Ngrid(2)+LLG(2,i), Ngrid(2) )
          i3=mod( Ngrid(3)+LLG(3,i), Ngrid(3) )
          if ( a2b <= i2 .and. i2 <= b2b .and. a3b <= i3 .and. i3 <= b3b ) then
             n=n+1
             LGHT(1,n)=i1
             LGHT(2,n)=i2
             LGHT(3,n)=i3
             IGHT(n)=i
          end if
       end do
       NGHT=n
       allocate( fg(MG) ) ; fg=(0.0d0,0.0d0)
       call destruct_Ggrid
       first_time1=.false.
    end if

    ctt(:)=0.d0
    ett(:)=0.d0

    call watch(ctt(0),ett(0))

!$OMP parallel private(j,ik)
!$OMP do
    do i=MG_0,MG_1
       j=MGL(i)
       fg(i)=vqlg(j,1)*SGK(i,1)
    end do
!$OMP end do
    do ik=2,Nelement
!$OMP do
       do i=MG_0,MG_1
          j=MGL(i)
          fg(i)=fg(i)+vqlg(j,ik)*SGK(i,ik)
       end do
!$OMP end do
    end do
!$OMP end parallel

    call allgatherv_Ggrid(fg)

!$OMP parallel
!$OMP workshare
    zwork1_ffte(:,:,:)=(0.0d0,0.0d0)
!$OMP end workshare
!$OMP do
    do i=1,NGHT
       zwork1_ffte(LGHT(1,i),LGHT(2,i),LGHT(3,i)) = fg(IGHT(i))
    end do
!$OMP end do
!$OMP end parallel

    call watch(ctt(1),ett(1))

    call pzfft3dv(zwork1_ffte,zwork2_ffte,ML1,ML2,ML3 &
         ,comm_ffty,comm_fftz,npuy,npuz,1)

    call watch(ctt(2),ett(2))

!$OMP parallel do collapse(3) private(i)
    do i3=a3b,b3b
    do i2=a2b,b2b
    do i1=a1b,b1b
       i = ML_0 + i1-a1b + (i2-a2b)*ab1 + (i3-a3b)*ab12
       Vion(i)=real( zwork2_ffte(i1,i2,i3) )*ML
    end do
    end do
    end do
!$OMP end parallel do

    call watch(ctt(3),ett(3))

!    if ( disp_switch_parallel ) then
!       write(*,*) "time(const_ps_loc_ffte1)",ctt(1)-ctt(0),ett(1)-ett(0)
!       write(*,*) "time(const_ps_loc_ffte2)",ctt(2)-ctt(1),ett(2)-ett(1)
!       write(*,*) "time(const_ps_loc_ffte3)",ctt(3)-ctt(2),ett(3)-ett(2)
!    end if

    call write_border( 80, " construct_ps_local_ffte(end)" )

  END SUBROUTINE construct_ps_local_ffte


END MODULE ps_local_module
