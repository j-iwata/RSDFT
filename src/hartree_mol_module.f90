MODULE hartree_mol_module

!$ use omp_lib
  use rgrid_mol_module, only: LL,KK
  use rgrid_module, only: dV,Hgrid,Ngrid,Igrid
  use parallel_module
  use bc_module
  use fd_module
  use atom_module
  use watch_module
  use io_tools_module
  use ylm_module

  implicit none

  PRIVATE
  PUBLIC :: calc_hartree_mol, init_hartree_mol
  PUBLIC :: timer_reset_hartree_mol, timer_result_hartree_mol

  integer :: MEO = 1
  integer :: Lmax_ME = 8

  logical :: flag_init = .true.
  integer :: M_max, lmmax_ME
  real(8),allocatable :: shf1(:,:),shf2(:,:)
  real(8),allocatable :: lap(:)
  integer :: NMadv
  integer,allocatable :: adv(:),Mdv(:),Ixyz(:,:)
  integer :: Md

  real(8),allocatable :: rholm_0(:,:,:), rholm(:,:,:)
  real(8),allocatable :: vb(:,:)
  real(8),allocatable :: coef(:),Clm(:,:)

  real(8) :: time_ht(2,16)
  character(5) :: indx_ht(16)
  integer :: iter_count_ht(0:3)
  real(8) :: err_chk_ht(2)
  real(8) :: info_ht(2)

  integer :: nthreads
  integer,allocatable :: ir_omp(:),id_omp(:)
  real(8),allocatable :: sum_tmp1(:)

CONTAINS


  SUBROUTINE read_param_hartree_mol
    implicit none
    integer :: itmp(2)
    itmp=-1
    call IOTools_readIntegerKeywords( "MEO", itmp )
    if ( itmp(1) > 0 ) MEO = itmp(1)
    if ( itmp(2) > 0 ) Lmax_ME = itmp(2)
  END SUBROUTINE read_param_hartree_mol


  SUBROUTINE init_hartree_mol( Md_in )
    implicit none
    integer,intent(IN) :: Md_in
    integer :: i,j

! ---

    Md = Md_in

    if ( allocated(lap) ) deallocate(lap)
    allocate( lap(-Md:Md) )
    lap=0.0d0
    call get_coef_lapla_fd(Md,lap)
    lap=lap/Hgrid(1)**2

! ---

    call read_param_hartree_mol

    M_max = Lmax_ME

    lmmax_ME = (Lmax_ME+1)**2

    select case(MEO)
    case(1)
       call prep1_hartree_mol
    case(2)
       call prep2_hartree_mol
       call init_MultiCenterExpansion
    end select

! --- for OpenMP

    nthreads = 1

!$OMP parallel
!$OMP single

!$  nthreads = omp_get_num_threads()

!$OMP end single
!$OMP end parallel

    if ( disp_switch_parallel ) write(*,*) "nthreads=",nthreads

    allocate( ir_omp(0:nthreads-1) ) ; ir_omp=0
    allocate( id_omp(0:nthreads-1) ) ; id_omp=0
    do i=1,ir_grid(myrank_g)
       j=mod(i-1,nthreads)
       ir_omp(j)=ir_omp(j)+1
    end do
    do j=0,nthreads-1
       id_omp(j) = sum( ir_omp(0:j) ) - ir_omp(j) + id_grid(myrank_g)
    end do

    allocate( sum_tmp1(0:nthreads-1) ) ; sum_tmp1=0.0d0

! ---

    flag_init = .false.

    if ( disp_switch_parallel ) then
       write(*,*) "MEO     =",MEO
       write(*,*) "Lmax_ME =",Lmax_ME
       write(*,*) "M_max   =",M_max
       write(*,*) "lmmax_ME=",lmmax_ME
       write(*,*) "Md      =",Md
    end if

  END SUBROUTINE init_hartree_mol


!--------1---------2---------3---------4---------5---------6---------7--
!
! Hartree potential ( with real(8) density )
! ( The reference of this C.-G. Poisson solver is
!                  Phys. Rev. C17, 1682 (1978). )
!
  SUBROUTINE calc_hartree_mol(n1,n2,n3,tn,Vh,Eh,tol,max_iter)
    implicit none
    integer,intent(IN)    :: n1,n2,n3
    real(8),intent(IN)    :: tn(n1:n2,n3)
    real(8),intent(INOUT) :: Vh(n1:n2)
    real(8),intent(OUT)   :: Eh
    real(8),optional,intent(IN) :: tol
    integer,optional,intent(IN) :: max_iter
    integer :: i,ix,iy,iz,iter,ierr,j,m1,m2
    integer :: maxiter,n1_omp,n2_omp,mythread
    integer,parameter :: maxiter_def = 2000
    real(8),parameter :: tol_def = 1.d-20
    real(8),allocatable :: tk(:),zk(:),qk(:)
    real(8) :: sum0,sum1,sum2,b0,ak,ck,const,c,d,pi4,ep
    real(8) :: time_tmp(2),time_tmp_0(2)
    logical :: flag_alloc(2)

    if ( flag_init ) then
       write(*,*) "init_hartree should be called first"
       stop "stop@calc_hartree_mol"
    end if

!    call timer_reset_hartree_mol

    !call watchb( time_tmp ) ; time_tmp_0=time_tmp

! ---
    ep = tol_def
    if ( present(tol) ) ep = tol

    maxiter = maxiter_def
    if ( present(max_iter) ) maxiter = max_iter

    info_ht(1) = ep
    info_ht(2) = maxiter
! ---

    m1  = 1
    m2  = size(KK,2)
    pi4 = 4.0d0*acos(-1.0d0)
    Eh  = 0.0d0

    sum0 = sum( tn(n1:n2,1)*tn(n1:n2,1) )*dV
    call mpi_allreduce(sum0,b0,1,mpi_real8,mpi_sum,comm_grid,ierr)

    if ( b0 == 0.d0 ) then
       Vh(n1:n2) = 0.0d0
       return
    end if

    allocate( tk(n1:n2) ) ; tk=0.0d0
    allocate( zk(n1:n2) ) ; zk=0.0d0
    allocate( qk(n1:n2) ) ; qk=0.0d0

    !call watchb( time_tmp, time_ht(:,1) )

!
! - Boundary condition 1
!

!$OMP parallel
!$OMP do
    do i=n1,n2
       www( LL(1,i),LL(2,i),LL(3,i),1 ) = Vh(i)
    end do
!$OMP end do
    call bcset_1(1,1,Md,0)
!$OMP end parallel

    !call watchb( time_tmp, time_ht(:,2) )

!
! --- Boundary condition 2 (multipole expansion) ---
!

    allocate( vb(m1:m2,n3) )

!$OMP parallel workshare
    vb=0.0d0
!$OMP end parallel workshare

    select case( MEO )
    case(1)
       call SingleCenterExpansion( n1, n2, n3, m1, m2, tn, vb )
    case(2)
       call MultiCenterExpansion( n1, n2, n3, m1, m2, tn, vb )
    end select ! MEO

    !call watchb( time_tmp, time_ht(:,3) )

!$OMP parallel do
    do i=m1,m2
       www( KK(1,i),KK(2,i),KK(3,i),1 ) = vb(i,1)
    end do
!$OMP end parallel do

    deallocate( vb )

    !call watchb( time_tmp, time_ht(:,4) )

!
! --- C.-G. minimization ---
!

    const = 3.d0*lap(0)

!$OMP parallel do
    do i=n1,n2
       zk(i) = -const*Vh(i) - pi4*tn(i,1)
    end do
!$OMP end parallel do

    do j=1,Md
       c=lap(j)
!$OMP parallel do private( ix,iy,iz )
       do i=n1,n2
          ix=LL(1,i) ; iy=LL(2,i) ; iz=LL(3,i)
          zk(i)=zk(i)-c*( www(ix-j,iy,iz,1)+www(ix+j,iy,iz,1) &
                         +www(ix,iy-j,iz,1)+www(ix,iy+j,iz,1) &
                         +www(ix,iy,iz-j,1)+www(ix,iy,iz+j,1) )
       end do
!$OMP end parallel do
    end do

!$OMP parallel workshare
    www(:,:,:,1)=0.d0
!$OMP end parallel workshare

    ak=0.0d0
    ck=0.0d0

    sum0=sum(zk(n1:n2)*zk(n1:n2))*dV
    call mpi_allreduce(sum0,sum1,1,mpi_real8,mpi_sum,comm_grid,ierr)

    !call watchb( time_tmp, time_ht(:,5) )

!
! --- Iteration ---
!
    Iteration : do iter=1,maxiter

       !call watchb( time_tmp )

!$OMP parallel private( i,j,ix,iy,iz,n1_omp,n2_omp,mythread,c,d )

       mythread = 0
!$     mythread = omp_get_thread_num()
       n1_omp = id_omp(mythread) + 1
       n2_omp = id_omp(mythread) + ir_omp(mythread)

       do i=n1_omp,n2_omp
          ix = LL(1,i)
          iy = LL(2,i)
          iz = LL(3,i)
          d = qk(i)
          Vh(i) = Vh(i) + ak*d
          d = d*ck + zk(i)
          www(ix,iy,iz,1) = d
          qk(i) = d
       end do
!$OMP barrier

       call bcset_3(1,1,Md,0)
!$OMP barrier

       !call watchb_omp( time_tmp, time_ht(:,6) )

       do i=n1_omp,n2_omp
          ix = LL(1,i)
          iy = LL(2,i)
          iz = LL(3,i)
          d = const*qk(i)
          do j=1,Md
             c=lap(j)
             d=d+c*( www(ix-j,iy,iz,1)+www(ix+j,iy,iz,1) &
                    +www(ix,iy-j,iz,1)+www(ix,iy+j,iz,1) &
                    +www(ix,iy,iz-j,1)+www(ix,iy,iz+j,1) )
          end do
          tk(i) = d
       end do

       !call watchb_omp( time_tmp, time_ht(:,7) )

       sum_tmp1(mythread) = sum( zk(n1_omp:n2_omp)*tk(n1_omp:n2_omp) )
!$OMP barrier

!$OMP single
       sum0 = sum( sum_tmp1 )*dV
       call mpi_allreduce(sum0,sum2,1,mpi_real8,mpi_sum,comm_grid,ierr)
       ak = sum1/sum2
!$OMP end single

       !call watchb_omp( time_tmp, time_ht(:,8) )

       do i=n1_omp,n2_omp
          zk(i) = zk(i) - ak*tk(i)
       end do

       !call watchb_omp( time_tmp, time_ht(:,9) )

! --- Conv. Check ---

       sum_tmp1(mythread) = sum( zk(n1_omp:n2_omp)*zk(n1_omp:n2_omp) )
!$OMP barrier

!$OMP single
       sum0 = sum( sum_tmp1 )*dV
       call mpi_allreduce(sum0,sum2,1,mpi_real8,mpi_sum,comm_grid,ierr)
       ck   = sum2/sum1
       sum1 = sum2
!$OMP end single

!$OMP end parallel

       if ( sum2/b0 <= ep ) exit

       !call watchb( time_tmp, time_ht(:,10) )

    end do Iteration

    deallocate( qk )
    deallocate( zk )
    deallocate( tk )

!    call timer_result_hartree_mol

    iter_count_ht(0) = iter_count_ht(0) + 1
    iter_count_ht(1) = min( iter, iter_count_ht(1) )
    iter_count_ht(2) = max( iter, iter_count_ht(2) )
    iter_count_ht(3) = iter_count_ht(3) + iter

    !call watchb( time_tmp_0, time_ht(:,11) )

    return

  END SUBROUTINE calc_hartree_mol

!--------1---------2---------3---------4---------5---------6---------7--
!
! Spherical Harmonic Funciton
!
  SUBROUTINE prep1_hartree_mol
    implicit none
    integer :: i,j,k,lm,L,a,m,n,n1,n2,ML0
    integer :: m1,m2
    logical :: flag_alloc(2)
    real(8),parameter :: eps=1.d-20
    real(8) :: r,x,y,z,const,H,pi4

    n1  = idisp(myrank)+1
    n2  = idisp(myrank)+ircnt(myrank)
    ML0 = ircnt(myrank)
    m1  = 1
    m2  = size(KK,2)
    H   = Hgrid(1)
    pi4 = 4.d0*acos(-1.d0)

    if ( allocated(shf1) ) deallocate(shf1)
    if ( allocated(shf2) ) deallocate(shf2)
    allocate( shf1(n1:n2,lmmax_ME) ) ; shf1=0.0d0
    allocate( shf2(lmmax_ME,m1:m2) ) ; shf2=0.0d0

    do i=n1,n2
       x=LL(1,i)*H ; y=LL(2,i)*H ; z=LL(3,i)*H
       r=sqrt(x*x+y*y+z*z)
       lm=0
       do L=0,Lmax_ME
          do m=-L,L
             lm=lm+1
             shf1(i,lm)=Ylm(x,y,z,L,m)*r**L
          end do
       end do
    end do

    do i=m1,m2
       x=KK(1,i)*H ; y=KK(2,i)*H ; z=KK(3,i)*H
       r=sqrt(x*x+y*y+z*z)
       lm=0
       do L=0,Lmax_ME
          const=pi4/(2.d0*L+1.d0)
          do m=-L,L
             lm=lm+1
             shf2(lm,i)=Ylm(x,y,z,L,m)/r**(L+1)*const
          end do
       end do
    end do

    return
  END SUBROUTINE prep1_hartree_mol

!--------1---------2---------3---------4---------5---------6---------7--
!
! Division of Area for Multipole Expansion (MEO=2)
!
  SUBROUTINE prep2_hartree_mol
    implicit none
    real(8),allocatable :: ra(:)
    real(8) :: x,y,z,r2,H
    integer,allocatable :: itmp(:),jtmp(:)
    integer :: i,a,m,n,i1,i2,i3,n1,n2,ierr,maxMdv

    n1 = idisp(myrank)+1
    n2 = idisp(myrank)+ircnt(myrank)
    H  = Hgrid(1)

    allocate( itmp(n1:n2) ) ; itmp=0
    allocate( jtmp(Natom) ) ; jtmp=0

    allocate( ra(n1:n2) )
    ra=1.d10
    do a=1,Natom
       do i=n1,n2
          x=LL(1,i)*H-aa_atom(1,a)
          y=LL(2,i)*H-aa_atom(2,a)
          z=LL(3,i)*H-aa_atom(3,a)
          r2=x*x+y*y+z*z
          if ( r2 < ra(i) ) then
             ra(i)=r2
             itmp(i)=a
          end if
       end do
    end do
    deallocate( ra )

    do a=1,Natom
       jtmp(a)=count(itmp==a) ! # of grid points near the atom a
    end do

    NMadv=count( jtmp>0 )     ! # of atoms (or # of regions) in my rank
    maxMdv=maxval(jtmp)       ! max # of grid points around each atom

!    if ( DISP_SWITCH ) then
       write(*,*) "NMadv,maxMdv,ML,sum(jtmp)=",NMadv,maxMdv,Ngrid(0),sum(jtmp)
!    end if

    if ( allocated(Ixyz) ) deallocate(Ixyz)
    if ( allocated(Mdv)  ) deallocate(Mdv)
    if ( allocated(adv)  ) deallocate(adv)
    allocate( Ixyz(maxMdv,NMadv) )
    allocate( Mdv(NMadv) )
    allocate( adv(NMadv) )

    n=0
    do a=1,Natom
       if ( jtmp(a) > 0 ) then
          n=n+1
          Mdv(n)=jtmp(a)
          adv(n)=a
       end if
    end do

    if ( n /= NMadv ) stop "n/=NMadv!!!"

    do n=1,NMadv
       a=adv(n)
       m=0
       do i=n1,n2
          if ( itmp(i) == a ) then
             m=m+1
             Ixyz(m,n)=i
          end if
       end do
    end do

    deallocate( jtmp,itmp )

    return
  END SUBROUTINE prep2_hartree_mol


  SUBROUTINE SingleCenterExpansion( n1, n2, n3, m1, m2, tn, vb )
    implicit none
    integer,intent(IN)  :: n1,n2,n3,m1,m2
    real(8),intent(IN)  :: tn(n1:n2,n3)
    real(8),intent(INOUT) :: vb(m1:m2,n3)
    integer :: lm,n,i,ierr

    allocate( rholm_0(lmmax_ME,1,n3) ) ; rholm_0=0.0d0
    allocate( rholm(lmmax_ME,1,n3)   ) ; rholm=0.0d0

    do n=1,n3
!$OMP parallel do
       do lm=1,lmmax_ME
          rholm_0(lm,1,n) = sum( tn(n1:n2,n)*shf1(n1:n2,lm) )*dV
       end do
!$OMP end parallel do
    end do

    call mpi_allreduce( rholm_0, rholm, lmmax_ME*n3 &
                       ,mpi_real8, mpi_sum, comm_grid, ierr )

    do n=1,n3
!$OMP parallel do
       do i=m1,m2
          vb(i,n) = sum( shf2(1:lmmax_ME,i)*rholm(1:lmmax_ME,1,n) )
       end do
!$OMP end parallel do
    end do

    deallocate( rholm )
    deallocate( rholm_0 )

  END SUBROUTINE SingleCenterExpansion


  SUBROUTINE init_MultiCenterExpansion
    implicit none
    integer :: m,i,L,k,k1,k2
    real(8) :: pi

    pi = acos(-1.0d0)

    allocate( coef(0:M_max)          ) ; coef=0.0d0
    allocate( Clm(0:Lmax_ME,0:M_max) ) ; Clm=0.0d0

! Pmm=(-1)^m*(2m-1)!!*(1-x^2)^(m/2)

    coef(:)=1.0d0
    do m=1,M_max
       do i=1,2*m-1,2
          coef(m)=coef(m)*i
       end do
       coef(m)=coef(m)*(-1)**m
    end do

    Clm(:,:)=1.0d0
    do L=0,Lmax_ME
       do m=0,L
          k1=L-m
          k2=L+m
          if ( k1 < k2 ) then
             do k=k2,k1+1,-1
                Clm(L,m)=Clm(L,m)*k
             end do
             Clm(L,m)=1.0d0/Clm(L,m)
          else if ( k1 > k2 ) then
             do k=k1,k2+1,-1
                Clm(L,m)=Clm(L,m)*k
             end do
          end if
          Clm(L,m)=sqrt(Clm(L,m)*(2*L+1)/pi)*0.5d0
          if ( m /= 0 ) Clm(L,m)=Clm(L,m)*sqrt(2.0d0)
       end do ! m
    end do ! L

  END SUBROUTINE init_MultiCenterExpansion


  SUBROUTINE MultiCenterExpansion( n1, n2, n3, m1, m2, tn, vb )
    implicit none
    integer,intent(IN)  :: n1,n2,n3,m1,m2
    real(8),intent(IN)  :: tn(n1:n2,n3)
    real(8),intent(INOUT) :: vb(m1:m2,n3)
    integer :: n,a,i,j,k,idv,m,L,lm,ix,iy,iz,ierr
    real(8) :: x,y,z,r,H,pi,pi4
    real(8) :: ak,ck,plm(3)
    real(8),allocatable :: sk(:)

    pi  = acos(-1.0d0)
    pi4 = 4.0d0*pi
    H   = Hgrid(1)

    allocate( rholm_0(lmmax_ME,Natom,n3) ) ; rholm_0=0.0d0
    allocate( rholm(lmmax_ME,Natom,n3)   ) ; rholm=0.0d0

    allocate( sk(0:Lmax_ME) ) ; sk=0.0d0

    n=1

    do idv=1,NMadv

       a=adv(idv)

       do j=1,Mdv(idv)

          i=Ixyz(j,idv)
          x=LL(1,i)*H-aa_atom(1,a)
          y=LL(2,i)*H-aa_atom(2,a)
          z=LL(3,i)*H-aa_atom(3,a)
          r=sqrt(x*x+y*y+z*z)

          if ( r == 0.0d0 ) then
             rholm_0(1,a,n) = rholm_0(1,a,n) + tn(i,n)*Clm(0,0)
             cycle
          end if

          ck=z/r ; if ( abs(ck) > 1.0d0 ) ck=sign(1.0d0,ck)
          if ( abs(x) < 1.d-10 ) then
             ak=0.5d0*pi
             if ( y < 0.0d0 ) ak=ak+pi
          else
             ak=atan(y/x)
             if ( x < 0.0d0 ) ak=ak+pi
          end if

          sk(0)=tn(i,n)
          do L=1,Lmax_ME
             sk(L)=tn(i,n)*r**L
          end do

          lm=0
          do m=0,0

             plm(1) = coef(m)

             lm=lm+1
             rholm_0(lm,a,n) = rholm_0(lm,a,n) + plm(1)*Clm(0,0)*sk(0)

             plm(2) = ck*(2*m+1)*plm(1)

             lm=lm+1
             rholm_0(lm,a,n) = rholm_0(lm,a,n) + plm(2)*Clm(1,0)*sk(1)

             do L=m+2,Lmax_ME

                plm(3) = ( ck*(2*L-1)*plm(2)-(L+m-1)*plm(1) )/dble(L-m)

                lm=lm+1
                rholm_0(lm,a,n) = rholm_0(lm,a,n) + plm(3)*Clm(L,0)*sk(L)

                plm(1) = plm(2)
                plm(2) = plm(3)

             end do ! L

          end do ! m

          do m=1,M_max-1

             plm(1) = coef(m)*(1.0d0-ck*ck)**(0.5d0*m)

             lm=lm+1
             rholm_0(lm,a,n) = rholm_0(lm,a,n) + plm(1)*Clm(m,m)*sk(m)*cos(ak*m)

             lm=lm+1
             rholm_0(lm,a,n) = rholm_0(lm,a,n) - plm(1)*Clm(m,m)*sk(m)*sin(ak*m)*(-1)**m

             plm(2) = ck*(2*m+1)*plm(1)

             lm=lm+1
             rholm_0(lm,a,n) = rholm_0(lm,a,n) + plm(2)*Clm(m+1,m)*sk(m+1)*cos(ak*m)

             lm=lm+1
             rholm_0(lm,a,n) = rholm_0(lm,a,n) - plm(2)*Clm(m+1,m)*sk(m+1)*sin(ak*m)*(-1)**m

             do L=m+2,Lmax_ME

                plm(3) = ( ck*(2*L-1)*plm(2)-(L+m-1)*plm(1) )/dble(L-m)

                lm=lm+1
                rholm_0(lm,a,n) = rholm_0(lm,a,n) + plm(3)*Clm(L,m)*sk(L)*cos(ak*m)

                lm=lm+1
                rholm_0(lm,a,n) = rholm_0(lm,a,n) - plm(3)*Clm(L,m)*sk(L)*sin(ak*m)*(-1)**m

                plm(1) = plm(2)
                plm(2) = plm(3)

             end do ! L

          end do ! m

          do m=M_max,M_max

             plm(1) = coef(m)*(1.0d0-ck*ck)**(0.5d0*m)

             lm=lm+1
             rholm_0(lm,a,n) = rholm_0(lm,a,n) + plm(1)*Clm(m,m)*sk(m)*cos(ak*m)

             lm=lm+1
             rholm_0(lm,a,n) = rholm_0(lm,a,n) - plm(1)*Clm(m,m)*sk(m)*sin(ak*m)*(-1)**m

          end do ! m

       end do ! j

    end do ! n

    rholm_0(:,:,:)=rholm_0(:,:,:)*dV
    call mpi_allreduce( rholm_0, rholm, lmmax_ME*Natom*n3 &
                       ,mpi_real8, mpi_sum, comm_grid, ierr )

! ---

    do a=1,Natom

       do i=m1,m2

          ix=KK(1,i) ; iy=KK(2,i) ; iz=KK(3,i)
          x=ix*H-aa_atom(1,a)
          y=iy*H-aa_atom(2,a)
          z=iz*H-aa_atom(3,a)
          r=sqrt(x*x+y*y+z*z)

          do L=0,Lmax_ME
             sk(L) = pi4/( (2.0d0*L+1.0d0)*r**(L+1) )
          end do

          ck=z/r
          if ( abs(x) < 1.d-10 ) then
             ak=0.5d0*pi
             if ( y < 0.d0 ) ak=ak+pi
          else
             ak=atan(y/x)
             if ( x < 0.d0 ) ak=ak+pi
          end if

          lm=0
          do m=0,0

             plm(1) = coef(m)

             lm=lm+1
             vb(i,n) = vb(i,n) + plm(1)*Clm(0,0)*sk(0)*rholm(lm,a,n)

             plm(2) = ck*(2*m+1)*plm(1)

             lm=lm+1
             vb(i,n) = vb(i,n) + plm(2)*Clm(1,0)*sk(1)*rholm(lm,a,n)

             do L=m+2,Lmax_ME

                plm(3) = ( ck*(2*L-1)*plm(2)-(L+m-1)*plm(1) )/dble(L-m)

                lm=lm+1
                vb(i,n) = vb(i,n) + plm(3)*Clm(L,0)*sk(L)*rholm(lm,a,n)

                plm(1) = plm(2)
                plm(2) = plm(3)

             end do ! L

          end do ! m

          do m=1,M_max-1

             plm(1) = coef(m)*(1.0d0-ck*ck)**(0.5d0*m)

             lm=lm+1
             vb(i,n) = vb(i,n) + plm(1)*Clm(m,m)*sk(m)*cos(ak*m)*rholm(lm,a,n)

             lm=lm+1
             vb(i,n) = vb(i,n) - plm(1)*Clm(m,m)*sk(m)*sin(ak*m)*(-1)**m*rholm(lm,a,n)

             plm(2) = ck*(2*m+1)*plm(1)

             lm=lm+1
             vb(i,n) = vb(i,n) + plm(2)*Clm(m+1,m)*sk(m+1)*cos(ak*m)*rholm(lm,a,n)

             lm=lm+1
             vb(i,n) = vb(i,n) - plm(2)*Clm(m+1,m)*sk(m+1)*sin(ak*m)*(-1)**m*rholm(lm,a,n)

             do L=m+2,Lmax_ME

                plm(3) = ( ck*(2*L-1)*plm(2)-(L+m-1)*plm(1) )/dble(L-m)

                lm=lm+1
                vb(i,n) = vb(i,n) + plm(3)*Clm(L,m)*sk(L)*cos(ak*m)*rholm(lm,a,n)

                lm=lm+1
                vb(i,n) = vb(i,n) - plm(3)*Clm(L,m)*sk(L)*sin(ak*m)*(-1)**m*rholm(lm,a,n)

                plm(1) = plm(2)
                plm(2) = plm(3)

             end do ! L

          end do ! m

          do m=M_max,M_max

             plm(1) = coef(m)*(1.0d0-ck*ck)**(0.5d0*m)

             lm=lm+1
             vb(i,n) = vb(i,n) + plm(1)*Clm(m,m)*sk(m)*cos(ak*m)*rholm(lm,a,n)

             lm=lm+1
             vb(i,n) = vb(i,n) - plm(1)*Clm(m,m)*sk(m)*sin(ak*m)*(-1)**m*rholm(lm,a,n)

          end do ! m

       end do ! i

    end do ! a

    deallocate( sk )
    deallocate( rholm )
    deallocate( rholm_0 )

  END SUBROUTINE MultiCenterExpansion


  SUBROUTINE timer_reset_hartree_mol
    implicit none
    time_ht(:,:)=0.0d0
    indx_ht(1:11) = (/ "init","bc1","bc2","bc2-2" &
                      ,"cg1" ,"cg-bc","cg2","cg3","cg4","cg5","tot" /)
    iter_count_ht=0
    iter_count_ht(1)=1000000
!    err_chk_ht=0.0d0
!    err_chk_ht(1)=1.d100
    info_ht=0.0d0
  END SUBROUTINE timer_reset_hartree_mol

  SUBROUTINE timer_result_hartree_mol
    implicit none
    if ( disp_switch_parallel ) then
       write(*,'(a40," hartree_mol")') repeat("-",40)
       call write_watchb( time_ht, 11, indx_ht )
       write(*,'(1x,"# of iteration (min,max,ave,tot)=",3i6,i8)') &
             iter_count_ht(1:2) &
            ,iter_count_ht(3)/iter_count_ht(0),iter_count_ht(3)
       write(*,*) "ep,maxiter=",info_ht(1),nint(info_ht(2))
       write(*,'(a40)') repeat("-",40)
    end if
  END SUBROUTINE timer_result_hartree_mol


END MODULE hartree_mol_module
