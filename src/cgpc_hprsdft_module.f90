MODULE cgpc_hprsdft_module

!$ use omp_lib
  use rgrid_module, only: Igrid,Hgrid
  use rgrid_mol_module, only: Hsize, LL
  use bc_module
  use bcset_0_module, only: bcset_0
  use kinetic_variables, only: ggg
  use parallel_module
  use watch_module, only: watchb_omp, time_cgpc, time_bcfd

  implicit none

  PRIVATE
  PUBLIC :: precond_hprsdft, init_cgpc_hprsdft

#ifdef _DRSDFT_
  real(8),allocatable :: ftmp2(:,:),gtmp2(:,:)
  real(8),allocatable :: rk_pc(:,:), pk_pc(:,:)
  real(8),parameter :: zero=0.0d0
#else
  complex(8),allocatable :: ftmp2(:,:),gtmp2(:,:)
  complex(8),allocatable :: rk_pc(:,:), pk_pc(:,:)
  complex(8),parameter :: zero=(0.0d0,0.0d0)
#endif

  logical :: ompflag1
  integer :: omp1,omp2,ompi,ompnprocs,ompnsw,ompmyrank
  integer :: ompn1,ompn2,ompA3,ompB3
  real(8),allocatable :: ompr1s(:,:)

  integer :: a1b,b1b,a2b,b2b ,mloop,SYStype, ML_0
  real(8) :: H1,H2,H3
  logical :: init_flag=.false.

!$omp threadprivate( omp1,omp2,ompi,ompnsw,ompmyrank,ompnprocs &
!$omp               ,ompn1,ompn2,ompA3,ompB3 )

CONTAINS


  SUBROUTINE init_cgpc_hprsdft( mloop_in, systype_in )
    implicit none
    integer,intent(IN) :: mloop_in, systype_in
    integer :: mm,m1,m2,m3,i,j
    integer :: ompblock,ompblock0

    if ( init_flag ) return
    init_flag = .true.

    mloop = mloop_in
    SYStype = systype_in

    H1=Hgrid(1)
    H2=Hgrid(2)
    H3=Hgrid(3)

    mm = Igrid(2,0) - Igrid(1,0) + 1
    m1 = Igrid(2,1) - Igrid(1,1) + 1
    m2 = Igrid(2,2) - Igrid(1,2) + 1
    m3 = Igrid(2,3) - Igrid(1,3) + 1

    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)

    ML_0 = Igrid(1,0)

!$OMP parallel private( ompblock,ompblock0,i,j )

    ompnprocs = 1
!$  ompnprocs = omp_get_num_threads()
    ompmyrank = 0
!$  ompmyrank = omp_get_thread_num()

    ompblock0=0
    ompblock=0
    do i=1,mm
       j=mod(i-1,ompnprocs)
       if ( j == ompmyrank ) then
          ompblock = ompblock + 1
       else if ( j < ompmyrank ) then
          ompblock0 = ompblock0 + 1
       end if
    end do

    omp1 = ompblock0 + 1
    omp2 = omp1 + ompblock - 1

!    ompn1 = Igrid(1,0) + ompblock0
!    ompn2 = ompn1 + ompblock - 1
    ompn1 = omp1
    ompn2 = omp2

!$OMP single
    allocate( ompr1s(MB_d,0:ompnprocs-1) ) ; ompr1s=0.0d0
!$OMP end single

    ompblock0=0
    ompblock=0
    do i=1,m3
       j=mod(i-1,ompnprocs)
       if ( j == ompmyrank ) then
          ompblock = ompblock + 1
       else if ( j < ompmyrank ) then
          ompblock0 = ompblock0 + 1
       end if
    end do

    ompA3 = Igrid(1,3) + ompblock0
    ompB3 = ompA3 + ompblock - 1

    ompnsw = ( ompA3-Igrid(1,3) )*m1*m2

!$OMP end parallel

    allocate( rk_pc(mm,MB_d) ) ; rk_pc=zero
    allocate( pk_pc(mm,MB_d) ) ; pk_pc=zero
    allocate( gtmp2(mm,MB_d) ) ; gtmp2=zero
    allocate( ftmp2(mm,MB_d) ) ; ftmp2=zero

  END SUBROUTINE init_cgpc_hprsdft


  SUBROUTINE precond_hprsdft( E, k, s, nn, ML0, gk, Pgk )
    implicit none
    integer,intent(IN)       :: k,s,nn,ML0
    real(8),intent(IN)       :: E(nn)
#ifdef _DRSDFT_
    real(8),intent(INOUT) :: gk(ML0,nn), Pgk(ML0,nn)
#else
    complex(8),intent(INOUT) :: gk(ML0,nn), Pgk(ML0,nn)
#endif
    real(8),parameter :: ep=1.d-24
    real(8) :: rr0(nn),rr1(nn),pAp(nn),a(nn),E0(nn),b
    real(8) :: ttmp(2)
    integer :: n,ierr,iloop,mm

    if ( mloop <= 0 ) return

    E0(:) = 0.0d0
    ompflag1 = .FALSE.
    time_bcfd(:,:) = 0.0d0

!$omp parallel private ( n, iloop, b )

    !call watchb_omp( ttmp )

    gtmp2(omp1:omp2,1:nn) = gk(ompn1:ompn2,1:nn)

    !call watchb_omp( ttmp, time_cgpc(1,1) )

    call precond_cg_mat( E0, k, s, ML0, nn )

    !call watchb_omp( ttmp, time_cgpc(1,2) )

    do n=1,nn
       rk_pc(omp1:omp2,n) = gk(ompn1:ompn2,n) - ftmp2(omp1:omp2,n)
       pk_pc(omp1:omp2,n) = rk_pc(omp1:omp2,n)
    end do

    !call watchb_omp( ttmp, time_cgpc(1,1) )

    do n=1,nn
       ompr1s(n,ompmyrank) = sum( abs(rk_pc(omp1:omp2,n))**2 )
    end do
!$omp barrier
!$omp master
    do ompi=1,ompnprocs-1
       ompr1s(1:nn,0) = ompr1s(1:nn,0) + ompr1s(1:nn,ompi)
    end do
    call mpi_allreduce(ompr1s,rr1,nn,mpi_real8,mpi_sum,comm_grid,ierr)
    if ( all( rr1(1:nn)<ep ) ) ompflag1 = .TRUE.
!$omp end master
!$omp barrier

    !call watchb_omp( ttmp, time_cgpc(1,3) )

    if ( ompflag1 ) goto 99

    do iloop=1,mloop

       !call watchb_omp( ttmp )

       gtmp2(omp1:omp2,1:nn) = pk_pc(omp1:omp2,1:nn)

       !call watchb_omp( ttmp, time_cgpc(1,1) )

       call precond_cg_mat( E0, k, s, ML0, nn )

       !call watchb_omp( ttmp, time_cgpc(1,2) )

       do n=1,nn
#ifdef _DRSDFT_
          ompr1s(n,ompmyrank) &
               =sum( pk_pc(omp1:omp2,n)*ftmp2(omp1:omp2,n) )
#else
          ompr1s(n,ompmyrank) &
               =sum( conjg(pk_pc(omp1:omp2,n))*ftmp2(omp1:omp2,n) )
#endif
       end do
!$omp barrier
!$omp master
       do ompi=1,ompnprocs-1
          ompr1s(1:nn,0) = ompr1s(1:nn,0) + ompr1s(1:nn,ompi)
       end do
       call mpi_allreduce(ompr1s,pAp,nn,mpi_real8,mpi_sum,comm_grid,ierr)
       do n=1,nn
          rr0(n)=rr1(n)
          a(n)=rr0(n)/pAp(n)
       end do
!$omp end master
!$omp barrier

       !call watchb_omp( ttmp, time_cgpc(1,3) )

       do n=1,nn
          Pgk(ompn1:ompn2,n) = Pgk(ompn1:ompn2,n) + a(n)*pk_pc(omp1:omp2,n)
       end do

       !call watchb_omp( ttmp, time_cgpc(1,1) )

       if ( iloop == mloop ) exit

       do n=1,nn
          rk_pc(omp1:omp2,n)=rk_pc(omp1:omp2,n)-a(n)*ftmp2(omp1:omp2,n)
       end do

       !call watchb_omp( ttmp, time_cgpc(1,1) )

       do n=1,nn
          ompr1s(n,ompmyrank)=sum( abs(rk_pc(omp1:omp2,n))**2 )
       end do
!$omp barrier
!$omp master
       do ompi=1,ompnprocs-1
          ompr1s(1:nn,0)=ompr1s(1:nn,0)+ompr1s(1:nn,ompi)
       end do
       call mpi_allreduce(ompr1s,rr1,nn,mpi_real8,mpi_sum,comm_grid,ierr)
!$omp end master
!$omp barrier

       !call watchb_omp( ttmp, time_cgpc(1,3) )

       do n=1,nn
          b=rr1(n)/rr0(n)
          pk_pc(omp1:omp2,n) = rk_pc(omp1:omp2,n) + b*pk_pc(omp1:omp2,n)
       end do

       !call watchb_omp( ttmp, time_cgpc(1,1) )

    end do ! iloop

99  continue

!$omp end parallel

    time_cgpc(1:2,8:13) = time_cgpc(1:2,8:13) + time_bcfd(1:2,1:6)

    return

  END SUBROUTINE precond_hprsdft


  SUBROUTINE precond_cg_mat( E, k, s, mm, nn )
    implicit none
    integer,intent(IN) :: k,s,mm,nn
    real(8),intent(INOUT) :: E(nn)
    select case(SYStype)
    case(0)
       call precond_cg_mat2_omp(E,k,s,mm,nn)
    case(1)
       call precond_cg_mat_mol(E,k,s,mm,nn)
!    case(3)
!       call precond_cg_mat_esm(E,k,s,mm,nn)
    end select
  END SUBROUTINE precond_cg_mat


  SUBROUTINE precond_cg_mat2_omp(E,k,s,mm,nn)

    implicit none

    integer,intent(IN) :: k,s,mm,nn
    real(8) :: E(nn)
    real(8) :: c,c1,c2,c3,d,ttmp(2)
    integer :: n,i,i1,i2,i3

    !call watchb_omp( ttmp )

    c =ggg(1)/H1**2+ggg(2)/H2**2+ggg(3)/H3**2
    c1=-0.5d0/H1**2*ggg(1)
    c2=-0.5d0/H2**2*ggg(2)
    c3=-0.5d0/H3**2*ggg(3)

    do n=1,nn
       d=c-E(n)
       do i=omp1,omp2
          ftmp2(i,n)=d*gtmp2(i,n)
       end do
    end do

!$omp barrier
    !call watchb_omp( ttmp, time_cgpc(1,4) )

    do n=1,nn
       i=ompnsw
       do i3=ompA3,ompB3
       do i2=a2b,b2b
       do i1=a1b,b1b
          i=i+1
          www(i1,i2,i3,n)=gtmp2(i,n)
       end do
       end do
       end do
    end do

!$OMP barrier
    !call watchb_omp( ttmp, time_cgpc(1,5) )

    !call bcset_0(1,nn,1,0)
    call bcset_1(1,nn,1,0)

!$OMP barrier
    !call watchb_omp( ttmp, time_cgpc(1,6) )

    do n=1,nn
       i=ompnsw
       do i3=ompA3,ompB3
       do i2=a2b,b2b
       do i1=a1b,b1b
          i=i+1
          ftmp2(i,n)=ftmp2(i,n) &
               +c1*( www(i1+1,i2,i3,n)+www(i1-1,i2,i3,n) ) &
               +c2*( www(i1,i2+1,i3,n)+www(i1,i2-1,i3,n) ) &
               +c3*( www(i1,i2,i3+1,n)+www(i1,i2,i3-1,n) )
       end do
       end do
       end do
    end do

!$omp barrier
    !call watchb_omp( ttmp, time_cgpc(1,7) )

    return

  END SUBROUTINE precond_cg_mat2_omp


  SUBROUTINE precond_cg_mat_mol(E,k,s,mm,nn)
    implicit none
    integer,intent(IN) :: k,s,mm,nn
    real(8),intent(INOUT) :: E(nn)
    real(8) :: c,c1,c2,c3,d
    integer :: n,i,i1,i2,i3,j

    c  =  3.d0/Hsize**2
    c1 = -0.5d0/Hsize**2
    c2 = -0.5d0/Hsize**2
    c3 = -0.5d0/Hsize**2

    do n=1,nn
       d=c-E(n)
       do i=omp1,omp2
          ftmp2(i,n)=d*gtmp2(i,n)
       end do
    end do

!!$OMP workshare
!    www(:,:,:,:)=zero
!!$OMP end workshare
    do n=1,nn
       do i=omp1,omp2
          j=i+ML_0-1
          www( LL(1,j),LL(2,j),LL(3,j),n ) = gtmp2(i,n)
       end do
    end do
!$OMP barrier

    call bcset_1(1,nn,1,0)
!$OMP barrier

    do n=1,nn
       do i=omp1,omp2
          j=i+ML_0-1
          i1=LL(1,j)
          i2=LL(2,j)
          i3=LL(3,j)
          ftmp2(i,n)=ftmp2(i,n)+c1*( www(i1-1,i2,i3,n)+www(i1+1,i2,i3,n) ) &
                               +c2*( www(i1,i2-1,i3,n)+www(i1,i2+1,i3,n) ) &
                               +c3*( www(i1,i2,i3-1,n)+www(i1,i2,i3+1,n) )
       end do
    end do
!$OMP barrier

    return
  END SUBROUTINE precond_cg_mat_mol


END MODULE cgpc_hprsdft_module
