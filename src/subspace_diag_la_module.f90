MODULE subspace_diag_la_module

  use rgrid_module, only: dV,zdV
  use wf_module, only: unk,esp,hunk,iflag_hunk,MB_0_WF
  use hamiltonian_module
  use parallel_module
  use subspace_diag_variables
  use array_bound_module, only: ML_0,ML_1
  use watch_module

  implicit none

  PRIVATE
  PUBLIC :: subspace_diag_la

CONTAINS

  SUBROUTINE subspace_diag_la(k,s)
    implicit none
    integer,intent(IN) :: k,s
    integer :: ML0,n1,n2,m,n,ierr,MB
    complex(8),allocatable :: work(:)
    integer :: WORK1,WORK2
    integer,save :: LWORK=0,LIWORK,LRWORK
    character(6) :: idiag0
    integer,allocatable :: iwork(:)
    real(8),allocatable :: rwork(:)
#ifdef _DRSDFT_
    real(8),allocatable :: Htmp(:,:)
    real(8),allocatable :: psit(:,:)
    real(8) :: zz
#else
    complex(8),allocatable :: Htmp(:,:)
    complex(8),allocatable :: psit(:,:)
    complex(8) :: zz
#endif
    real(8) :: ct(9),et(9)
    real(8) :: rtmp(1)
    integer :: itmp(1)
    integer,allocatable :: ir(:),id(:)

    call write_border( 1, " subspace_diag_la(start)" )

    n1  = ML_0
    n2  = ML_1
    ML0 = ML_1-ML_0+1
    zz  = 0.5d0*zdV
    MB  = MB_diag

    ct(:)=0.d0
    et(:)=0.d0

#ifdef _DRSDFT_
    idiag0 = "dsyevd"
#else
    idiag0 = "zheevd"
#endif

    allocate( id(0:np_band-1),ir(0:np_band-1) )

    id(0:np_band-1) = id_band(0:np_band-1)*ML0
    ir(0:np_band-1) = ir_band(0:np_band-1)*ML0

    call mpi_allgatherv(unk(n1,MB_0_WF,k,s),ir(myrank_b),TYPE_MAIN &
                       ,unk(n1,1,k,s),ir,id,TYPE_MAIN,comm_band,ierr)
    if ( iflag_hunk >= 1 ) then
       call mpi_allgatherv(hunk(n1,MB_0_WF,k,s),ir(myrank_b),TYPE_MAIN &
                          ,hunk(n1,1,k,s),ir,id,TYPE_MAIN,comm_band,ierr)
    end if

    deallocate( ir,id )

    call watch(ct(1),et(1))

    allocate( Hsub(MB,MB) )
    Hsub=zero

    allocate( psit(ML0,MB) )

    call watch(ct(2),et(2))

    if ( iflag_hunk >= 1 ) then
       psit(:,:)=hunk(:,:,k,s)
    else
       do m=1,MB
          call hamiltonian(k,s,unk(n1,m,k,s),psit(1,m),n1,n2,m,m)
       end do
    end if

    call watch(ct(3),et(3))

#ifdef _DRSDFT_
    call dsyr2k('U','T',MB,ML0,zz,unk(n1,1,k,s),ML0,psit,ML0,zero,Hsub,MB)
!    call dgemm('C','N',MB,MB,ML0, dV,unk(n1,1,k,s),ML0,psit,ML0,zero,Hsub,MB)
#else
    call zher2k('U','C',MB,ML0,zz,unk(n1,1,k,s),ML0,psit,ML0,zero,Hsub,MB)
!    call zgemm('C','N',MB,MB,ML0,zdV,unk(n1,1,k,s),ML0,psit,ML0,zero,Hsub,MB)
!    do n=1,MB
!       do m=1,n
!          Hsub(m,n)=sum(conjg(unk(:,m,k,s))*psit(:,n))*dV
!          write(*,*) m,n,sum(conjg(unk(:,m,k,s))*psit(:,n))*dV
!          write(*,*) m,n,Hsub(m,n)
!          write(*,*) m,n,Hsub(n,m)
!       end do
!    end do
#endif

    call watch(ct(4),et(4))

    deallocate( psit )

! --- allreduce ---

    allocate( Htmp(MB,MB) )
    Htmp=Hsub

    call watch(ct(5),et(5))

    call mpi_allreduce(Htmp,Hsub,MB*MB,TYPE_MAIN,MPI_SUM,comm_grid,ierr)

    call watch(ct(6),et(6))

    deallocate( Htmp )

! --- solve eigenvalue problem ---

    call watch(ct(7),et(7))

    WORK1 = 0
    WORK2 = 0

    select case(idiag0)
    case('zheev')

       LWORK=max(LWORK,2*MB-1)
       LRWORK=3*MB-2
       allocate( work(LWORK),rwork(LRWORK) )

       call zheev('V','U',MB,Hsub,MB,esp(1,k,s),work,2*MB,rwork,ierr)
       if ( ierr==0 ) WORK1=nint( real(work(1)) )

       deallocate( work,rwork )

    case('zheevd')

       LWORK=max(LWORK,2*MB+MB*MB)
       LRWORK=1+5*MB+2*MB*MB ; LIWORK=3+5*MB
       allocate( work(LWORK),rwork(LRWORK),iwork(LIWORK) )

       call zheevd('V','U',MB,Hsub,MB,esp(1,k,s) &
                  ,work,LWORK,rwork,LRWORK,iwork,LIWORK,ierr)
       if ( ierr==0 ) WORK1=nint( real(work(1)) )

       deallocate( work,rwork,iwork )

    case('dsyev')

       if ( LWORK==0 ) LWORK=3*MB-1

       allocate( rwork(LWORK) )

       call DSYEV('V','U',MB,Hsub,MB,esp(1,k,s),rwork,LWORK,ierr)

       WORK1=nint(rwork(1))

       deallocate( rwork )

    case('dsyevd')

!       if ( LWORK==0  ) LWORK=1+6*MB+2*MB*MB
!       if ( LIWORK==0 ) LIWORK=3+5*MB
       if ( LWORK == 0 .and. LIWORK == 0 ) then
          call DSYEVD('V','U',MB,Hsub,MB,esp(1,k,s),rtmp,-1,itmp,-1,ierr)
          LWORK=nint(rtmp(1))
          LIWORK=itmp(1)
       end if

       allocate( rwork(LWORK),iwork(LIWORK) )

       call DSYEVD('V','U',MB,Hsub,MB,esp(1,k,s) &
                   ,rwork,LWORK,iwork,LIWORK,ierr)

       WORK1=nint(rwork(1))
       WORK2=iwork(1)

       deallocate( iwork,rwork )

    end select

    call watch(ct(8),et(8))

! --- Rotation ---

    allocate( psit(ML0,MB) ) ; psit=zero

#ifdef _DRSDFT_
    call dgemm('N','N',ML0,MB,MB,one,unk(n1,1,k,s),ML0 &
               ,Hsub(1,1),MB,zero,psit,ML0)
#else
    call zgemm('N','N',ML0,MB,MB,one,unk(n1,1,k,s),ML0 &
               ,Hsub(1,1),MB,zero,psit,ML0)
#endif
    unk(:,:,k,s)=psit(:,:)

    if ( iflag_hunk >= 1 ) then
#ifdef _DRSDFT_
       call dgemm('N','N',ML0,MB,MB,one,hunk(n1,1,k,s),ML0 &
            ,Hsub(1,1),MB,zero,psit,ML0)
#else
       call zgemm('N','N',ML0,MB,MB,one,hunk(n1,1,k,s),ML0 &
            ,Hsub(1,1),MB,zero,psit,ML0)
#endif
       hunk(:,:,k,s)=psit(:,:)
    end if

    deallocate( psit )

    call watch(ct(9),et(9))

    deallocate( Hsub )

    LWORK=max(LWORK,WORK1)
    LIWORK=max(LIWORK,WORK2)

    call write_border( 1, " subspace_diag_la(end)" )

    return

  END SUBROUTINE subspace_diag_la

END MODULE subspace_diag_la_module
