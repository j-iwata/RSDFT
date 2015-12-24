MODULE ps_nloc2_op_module

!$ use omp_lib
  use ps_nloc2_variables
  use watch_module
  use parallel_module
  use rgrid_module, only: dV,Igrid
  use rgrid_mol_module, only: iswitch_eqdiv

  implicit none

  PRIVATE
  PUBLIC :: op_ps_nloc2,    init_op_ps_nloc2
  PUBLIC :: op_ps_nloc2_hp, init_op_ps_nloc2_hp

  integer,allocatable :: ompMJJ1(:,:), ompMJJ2(:,:)
  integer,allocatable :: omplns1(:,:), omplns2(:,:)
  integer :: ompnzlma1,ompnzlma2,ompmyrank,ompnprocs,ompn1,ompn2

  integer,allocatable :: n_i2j(:)
  integer,allocatable :: i2j(:,:,:)

  logical :: init_flag=.false.

!$OMP threadprivate( ompmyrank,ompnprocs,ompnzlma1,ompnzlma2,ompn1,ompn2 )

CONTAINS


  SUBROUTINE init_op_ps_nloc2_hp( init_flag_in )
    implicit none
    logical,intent(IN) :: init_flag_in
    integer :: i,j,k,ompblock,ompblock0,mm

    if ( .not.init_flag_in ) return
    init_flag = .true.

    mm = Igrid(2,0) - Igrid(1,0) + 1

!$OMP parallel private( i,j,k,ompblock,ompblock0 )

    ompnprocs = 1
!$  ompnprocs = omp_get_num_threads()
    ompmyrank = 0
!$  ompmyrank = omp_get_thread_num()

    ompblock=0
    ompblock0=0
    do i=1,nzlma
       j=mod(i-1,ompnprocs)
       if ( j == ompmyrank ) then
          ompblock = ompblock + 1
       else if ( j < ompmyrank ) then
          ompblock0 = ompblock0 + 1
       end if
    end do

    ompnzlma1 = ompblock0 + 1
    ompnzlma2 = ompnzlma1 + ompblock - 1

! ---

!$OMP master
    if ( .not.allocated(omplns1) ) then
       allocate( omplns1(0:np_grid-1,0:ompnprocs-1) )
       allocate( omplns2(0:np_grid-1,0:ompnprocs-1) )
    end if
    omplns1=0
    omplns2=0
    if ( allocated(ompMJJ1) ) deallocate(ompMJJ1)
    if ( allocated(ompMJJ2) ) deallocate(ompMJJ2)
    allocate( ompMJJ1(0:nzlma,0:ompnprocs-1)     ) ; ompMJJ1=0
    allocate( ompMJJ2(0:nzlma,0:ompnprocs-1)     ) ; ompMJJ2=0
!$OMP end master
!$OMP barrier

    omplns1(:,ompmyrank)=-1
    omplns2(:,ompmyrank)=-2
    do i=0,np_grid-1
       do j=1,lma_nsend(i)
          k = sendmap(j,i)
          if ( ompnzlma1 <= k .and. k <= ompnzlma2 ) then
             omplns2(i,ompmyrank) = j
          end if
       end do ! j
       do j=lma_nsend(i),1,-1
          k = sendmap(j,i)
          if ( ompnzlma1 <= k .and. k <= ompnzlma2 ) then
             omplns1(i,ompmyrank) = j
          end if
       end do ! j
    end do ! i

! ---

    ompblock=0
    ompblock0=0
    do i=1,mm
       j=mod(i-1,ompnprocs)
       if ( j == ompmyrank ) then
          ompblock = ompblock + 1
       else if ( j < ompmyrank ) then
          ompblock0 = ompblock0 + 1
       end if
    end do

    ompn1 = Igrid(1,0) + ompblock0
    ompn2 = ompn1 + ompblock - 1

! ---

    ompMJJ1(:,ompmyrank)=-1
    ompMJJ2(:,ompmyrank)=-2
    do k=1,nzlma
       do j=1,MJJ(k)
          i=JJP(j,k)
          if ( ompn1 <= i .and. i<= ompn2 ) then
             ompMJJ1(k,ompmyrank) = j
             exit
          end if
       end do
       do j=MJJ(k),1,-1
          i=JJP(j,k)
          if ( ompn1 <= i .and. i<= ompn2 ) then
             ompMJJ2(k,ompmyrank) = j
             exit
          end if
       end do
    end do

!$OMP end parallel

  END SUBROUTINE init_op_ps_nloc2_hp


  SUBROUTINE op_ps_nloc2_hp(k,tpsi,htpsi,n1,n2,ib1,ib2)
    implicit none
    integer,intent(IN) :: k,n1,n2,ib1,ib2
#ifdef _DRSDFT_
    real(8),intent(IN)  :: tpsi(n1:n2,ib1:ib2)
    real(8),intent(OUT) :: htpsi(n1:n2,ib1:ib2)
    real(8) :: tmp_sum
#else
    complex(8),intent(IN)  :: tpsi(n1:n2,ib1:ib2)
    complex(8),intent(OUT) :: htpsi(n1:n2,ib1:ib2)
    complex(8) :: tmp_sum
#endif
    integer :: i,ib,i1,i2,i3,j,jb,lma,m,ML0,n,nb
    integer :: ierr,nreq
    integer :: irank,jrank,istatus(mpi_status_size,512),ireq(512)
    real(8) :: c, ttmp(2), ttmp1(2)

    if ( Mlma <= 0 ) return

    if ( .not.init_flag ) stop "@op_ps_nloc2_hp(init_flag=.false.)"

    nb = ib2-ib1+1

    !call watchb_omp( ttmp )

    do ib=ib1,ib2
       jb=ib-ib1+1
       do lma=ompnzlma1,ompnzlma2
          tmp_sum=zero
          do j=1,MJJ(lma)
             i=JJP(j,lma)
#ifdef _DRSDFT_
             tmp_sum = tmp_sum + uVk(j,lma,k)*tpsi(i,ib)
#else
             tmp_sum = tmp_sum + conjg(uVk(j,lma,k))*tpsi(i,ib)
#endif
          end do
          uVunk(lma,jb) = iuV(lma)*dV*tmp_sum
       end do ! lma
    end do ! ib

    !call watchb_omp( ttmp, time_nlpp(1,1) )

    select case( iswitch_eqdiv )
    case default

    do i=1,6

!$omp barrier
         
       select case(i)
       case(1,3,5)
          j=i+1
          do ib=1,nb
             uVunk0(ompnzlma1:ompnzlma2,ib)=uVunk(ompnzlma1:ompnzlma2,ib)
          end do
       case(2,4,6)
          j=i-1
       end select
!$omp barrier

       do m=1,nrlma_xyz(i)
          !call watchb_omp( ttmp1 )
          nreq=0
          irank=num_2_rank(m,i)
          jrank=num_2_rank(m,j)
          if ( irank >= 0 ) then
             do ib=1,nb
                i2=omplns1(irank,ompmyrank)-1+lma_nsend(irank)*(ib-1)
                do i1=omplns1(irank,ompmyrank),omplns2(irank,ompmyrank)
                   i2=i2+1
                   sbufnl(i2,irank)=uVunk0(sendmap(i1,irank),ib)
                end do
             end do
!$omp barrier
!$omp master
             nreq=nreq+1
             call mpi_isend(sbufnl(1,irank),lma_nsend(irank)*nb &
                  ,TYPE_MAIN,irank,1,comm_grid,ireq(nreq),ierr)
!$omp end master
          end if
          !call watchb_omp( ttmp1, time_nlpp(1,4) )
!$omp master
          if ( jrank >= 0 ) then
             nreq=nreq+1
             call mpi_irecv(rbufnl(1,jrank),lma_nsend(jrank)*nb &
                  ,TYPE_MAIN,jrank,1,comm_grid,ireq(nreq),ierr)
          end if
          !call watchb_omp( ttmp1, time_nlpp(1,5) )
          call mpi_waitall(nreq,ireq,istatus,ierr)
!$omp end master
!$omp barrier
          !call watchb_omp( ttmp1, time_nlpp(1,6) )
          if ( jrank >= 0 ) then
             do ib=1,nb
                i2=omplns1(jrank,ompmyrank)-1+lma_nsend(jrank)*(ib-1)
                do i1=omplns1(jrank,ompmyrank),omplns2(jrank,ompmyrank)
                   i2=i2+1
                   i3=recvmap(i1,jrank)
                   uVunk(i3,ib) = uVunk(i3,ib) + rbufnl(i2,jrank)
                end do
             end do
          end if
          !call watchb_omp( ttmp1, time_nlpp(1,7) )

       end do ! m

    end do ! i

    case( 2 )

!$OMP barrier
       call comm_eqdiv_ps_nloc2_mol(nzlma,ib1,ib2,uVunk)

    end select

!$omp barrier
    !call watchb_omp( ttmp, time_nlpp(1,2) )

    do ib=ib1,ib2
       jb=ib-ib1+1
       do lma=1,nzlma
          do j=ompMJJ1(lma,ompmyrank),ompMJJ2(lma,ompmyrank)
             i=JJP(j,lma)
             htpsi(i,ib)=htpsi(i,ib)+uVk(j,lma,k)*uVunk(lma,jb)
          end do
       end do
    end do

!$omp barrier
    !call watchb_omp( ttmp, time_nlpp(1,3) )

    return 
      
  END SUBROUTINE op_ps_nloc2_hp


  SUBROUTINE op_ps_nloc2(k,tpsi,htpsi,n1,n2,ib1,ib2)

    implicit none
    integer,intent(IN) :: k,n1,n2,ib1,ib2
#ifdef _DRSDFT_
    real(8),intent(IN)  :: tpsi(n1:n2,ib1:ib2)
    real(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
#else
    complex(8),intent(IN)  :: tpsi(n1:n2,ib1:ib2)
    complex(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
#endif
    integer :: i,ib,j,jb,i1,i2,i3,m,lma,nb,ierr,nreq,lmani,lmanj
    integer :: irank,jrank,istatus(mpi_status_size,512),ireq(512) 
    integer :: nb_0_omp,nb_1_omp,nzlma_0_omp,nzlma_1_omp,m_0,m_1,n_0,n_1
    complex(8) :: zc
    real(8) :: et0,et1,ttmp(2)

    if ( Mlma <= 0 ) return

!$OMP barrier
    !call watchb_omp( ttmp )

    nb = ib2-ib1+1

    call calc_range_omp(nb,nzlma,nb_0_omp,nb_1_omp,nzlma_0_omp,nzlma_1_omp)

    !call watchb_omp( ttmp, time_nlpp(1,4) )

    do ib=nb_0_omp,nb_1_omp
       jb=ib+ib1-1
    do lma=nzlma_0_omp,nzlma_1_omp
       uVunk(lma,ib)=zero
       do j=1,MJJ(lma)
          i=JJP(j,lma)
#ifdef _DRSDFT_
          uVunk(lma,ib)=uVunk(lma,ib)+uVk(j,lma,k)*tpsi(i,jb)
#else
          uVunk(lma,ib)=uVunk(lma,ib)+conjg(uVk(j,lma,k))*tpsi(i,jb)
#endif
       end do
       uVunk(lma,ib)=iuV(lma)*dV*uVunk(lma,ib)
    end do
    end do

    !call watchb_omp( ttmp, time_nlpp(1,1) )

    select case( iswitch_eqdiv )
    case default

       do i=1,6

!$OMP barrier

          select case(i)
          case(1,3,5)
             j=i+1
             do ib=nb_0_omp,nb_1_omp
             do lma=nzlma_0_omp,nzlma_1_omp
                uVunk0(lma,ib)=uVunk(lma,ib)
             end do
             end do
          case(2,4,6)
             j=i-1
          end select
!$OMP barrier

          do m=1,nrlma_xyz(i)
             irank=num_2_rank(m,i)
             jrank=num_2_rank(m,j)
             nreq=0 
             if ( irank >= 0 ) then
                lmani = lma_nsend(irank)
                call calc_range_omp(nb,lmani,m_0,m_1,n_0,n_1)
                do ib=m_0,m_1
                do i1=n_0,n_1
                   i2=i1+(ib-1)*lmani
                   sbufnl(i2,irank)=uVunk0(sendmap(i1,irank),ib)
                end do
                end do
!$OMP barrier
!$OMP master
                nreq=nreq+1
                call mpi_isend(sbufnl(1,irank),lmani*nb &
                     ,TYPE_MAIN,irank,1,comm_grid,ireq(nreq),ierr)
!$OMP end master
             end if
!$OMP master
             if ( jrank >= 0 ) then
                lmanj = lma_nsend(jrank)
                nreq=nreq+1
                call mpi_irecv(rbufnl(1,jrank),lmanj*nb &
                     ,TYPE_MAIN,jrank,1,comm_grid,ireq(nreq),ierr)
             end if
             call mpi_waitall(nreq,ireq,istatus,ierr)
!$OMP end master
!$OMP barrier
             if ( jrank >= 0 ) then
                lmanj = lma_nsend(jrank)
                call calc_range_omp(nb,lmanj,m_0,m_1,n_0,n_1)
                do ib=m_0,m_1
                do i1=n_0,n_1
                   i2=i1+(ib-1)*lmanj
                   i3=recvmap(i1,jrank)
                   uVunk(i3,ib) = uVunk(i3,ib) + rbufnl(i2,jrank)
                end do
                end do
             end if

          end do ! m

       end do ! i

    case( 2 )

!$OMP barrier
       call comm_eqdiv_ps_nloc2_mol(nzlma,ib1,ib2,uVunk)

    end select

!$OMP barrier
    !call watchb_omp( ttmp, time_nlpp(1,2) )

!    do ib=ib1,ib2
!       jb=ib-ib1+1
!!$OMP do
!       do i=n1,n2
!          do m=1,n_i2j(i)
!             j=i2j(1,m,i)
!             lma=i2j(2,m,i)
!             htpsi(i,ib) = htpsi(i,ib) + uVk(j,lma,k)*uVunk(lma,jb)
!          end do
!       end do
!!$OMP end do nowait
!    end do
    do ib=ib1,ib2
       jb=ib-ib1+1
       do lma=1,nzlma
!$OMP do
          do j=1,MJJ(lma)
             i=JJP(j,lma)
             htpsi(i,ib) = htpsi(i,ib) + uVk(j,lma,k)*uVunk(lma,jb)
          end do
!$OMP end do
       end do
    end do

!$OMP barrier
    !call watchb_omp( ttmp, time_nlpp(1,3) )

  END SUBROUTINE op_ps_nloc2

  SUBROUTINE calc_range_omp( nb, nz, nb_0, nb_1, nz_0, nz_1 )

    implicit none
    integer,intent(IN)  :: nb, nz
    integer,intent(OUT) :: nb_0,nb_1,nz_0,nz_1
    integer :: mp,ip,k0,k1,id,a,i,j
    real(8) :: et0,et1
!$  et0=omp_get_wtime()
    mp=1
!$  mp=omp_get_num_threads()
    ip=0
!$  ip=omp_get_thread_num()
    k0=gcd(nb,mp)
    k1=mp/k0
    a=-1
    loop_j : do j=0,k0-1
    do i=0,k1-1
       a=a+1
       if ( a == ip ) then
          nb_0 = j*(nb/k0) + 1
          nb_1 = nb_0 + (nb/k0) - 1
          nz_0 = i*(nz/k1) + 1
          nz_1 = nz_0 + (nz/k1) - 1
          if ( i == k1-1 ) nz_1 = nz
          exit loop_j
       end if
    end do
    end do loop_j
!    et1=omp_get_wtime()
!    write(*,'(1x,8i4)') ip,mp,nb_0,nb_1,nz_0,nz_1
!    write(*,*) "time=",et1-et0

  END SUBROUTINE calc_range_omp

  FUNCTION gcd(m0,n0)
    implicit none
    integer :: gcd,m0,n0
    integer :: m,n,mtmp,loop
    if ( m0 >= n0 ) then
       m=m0
       n=n0
    else
       m=n0
       n=m0
    end if
    do loop=1,10000
       if ( n == 0 ) exit
       mtmp = n
       n = mod(m,n)
       m = mtmp
    end do
    gcd = m
  END FUNCTION gcd

  SUBROUTINE init_op_ps_nloc2
    implicit none
    integer :: ML_0,ML_1,lma,i,j,n

    ML_0 = Igrid(1,0)
    ML_1 = Igrid(2,0)

! --- inverse map ( grid label i --> (j,lma) )

    if ( allocated(i2j) ) deallocate(i2j)
    if ( allocated(n_i2j) ) deallocate(n_i2j)

    allocate( n_i2j(ML_0:ML_1) ) ; n_i2j=0

    do lma=1,nzlma
       do j=1,MJJ(lma)
          i=JJP(j,lma)
          n_i2j(i)=n_i2j(i)+1
       end do
    end do

    n=maxval( n_i2j )
    allocate( i2j(2,n,ML_0:ML_1) ) ; i2j=0

    n_i2j(:)=0
    do lma=1,nzlma
       do j=1,MJJ(lma)
          i=JJP(j,lma)
          n_i2j(i)=n_i2j(i)+1
          i2j(1,n_i2j(i),i)=j
          i2j(2,n_i2j(i),i)=lma
       end do
    end do
  END SUBROUTINE init_op_ps_nloc2


  SUBROUTINE comm_eqdiv_ps_nloc2_mol(nzlma,ib1,ib2,uVunk)
    implicit none
    integer,intent(IN) :: nzlma,ib1,ib2
#ifdef _DRSDFT_
    real(8),intent(INOUT) :: uVunk(nzlma,ib1:ib2)
#else
    complex(8),intent(INOUT) :: uVunk(nzlma,ib1:ib2)
#endif
    integer :: nreq,irank,m,i1,i2,i3,ib,ierr,nb
    integer :: istatus(mpi_status_size,512),ireq(512)
    nb=ib2-ib1+1
    nreq=0
    do irank=0,nprocs_g-1
       m=lma_nsend(irank)
       if ( irank == myrank_g .or. m <= 0 ) cycle
!       i2=0
       do ib=ib1,ib2
!$OMP do
       do i1=1,m
!          i2=i2+1
          i2 = i1 + (ib-ib1)*m
          sbufnl(i2,irank)=uVunk(sendmap(i1,irank),ib)
       end do
!$OMP end do
       end do
!$OMP master
       nreq=nreq+1
       call mpi_isend(sbufnl(1,irank),m*nb,TYPE_MAIN,irank,1 &
            ,comm_grid,ireq(nreq),ierr)
       nreq=nreq+1
       call mpi_irecv(rbufnl(1,irank),m*nb,TYPE_MAIN,irank,1 &
            ,comm_grid,ireq(nreq),ierr)
!$OMP end master
    end do
!$OMP master
    if ( nreq > 0 ) call mpi_waitall(nreq,ireq,istatus,ierr)
!$OMP end master
!$OMP barrier
    do irank=0,nprocs_g-1
       m=lma_nsend(irank)
       if ( irank == myrank_g .or. m <= 0 ) cycle
!       i2=0
       do ib=ib1,ib2
!$OMP do
       do i1=1,m
!          i2=i2+1
          i2 = i1 + (ib-ib1)*m
          i3 = recvmap(i1,irank)
          uVunk(i3,ib) = uVunk(i3,ib) + rbufnl(i2,irank)
       end do
!$OMP end do
       end do
    end do
  END SUBROUTINE comm_eqdiv_ps_nloc2_mol


END MODULE ps_nloc2_op_module

