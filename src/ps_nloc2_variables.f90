MODULE ps_nloc2_variables

  use parallel_module, only: MPI_REAL8,MPI_COMPLEX16,nprocs_g

  integer :: Mlma,nzlma
  integer,allocatable :: JJ_MAP(:,:,:),MJJ_MAP(:),MJJ(:)
  integer :: MMJJ,MAXMJJ,nrlma_xyz(6)
  real(8),allocatable :: uV(:,:)
  integer,allocatable :: amap(:),lmap(:),mmap(:),iorbmap(:)
  integer,allocatable :: lma_nsend(:),iuV(:),nl_rank_map(:)
  integer,allocatable :: num_2_rank(:,:)
  integer,allocatable :: sendmap(:,:),recvmap(:,:)
  integer,allocatable :: JJP(:,:)
  integer :: nl_max_send
#ifdef _DRSDFT_
  real(8),allocatable :: uVk(:,:,:),sbufnl(:,:),rbufnl(:,:)
  real(8),allocatable :: xVk(:,:,:),yVk(:,:,:),zVk(:,:,:)
  real(8),allocatable :: uVunk(:,:),uVunk0(:,:)
  real(8),parameter :: zero=0.d0
  integer,parameter :: TYPE_MAIN=MPI_REAL8
#else
  complex(8),allocatable :: uVk(:,:,:),sbufnl(:,:),rbufnl(:,:)
  complex(8),allocatable :: xVk(:,:,:),yVk(:,:,:),zVk(:,:,:)
  complex(8),allocatable :: uVunk(:,:),uVunk0(:,:)
  complex(8),parameter :: zero=(0.d0,0.d0)
  integer,parameter :: TYPE_MAIN=MPI_COMPLEX16
#endif

CONTAINS

  SUBROUTINE allocate_ps_nloc2(MB_d)
    integer,intent(IN) :: MB_d
    integer :: n
    n=maxval( lma_nsend )*4*MB_d
    if ( allocated(rbufnl) ) deallocate(rbufnl)
    if ( allocated(sbufnl) ) deallocate(sbufnl)
    if ( allocated(uVunk)  ) deallocate(uVunk)
    if ( allocated(uVunk0) ) deallocate(uVunk0)
    allocate( sbufnl(n,0:nprocs_g-1) ) ; sbufnl=zero
    allocate( rbufnl(n,0:nprocs_g-1) ) ; rbufnl=zero
    allocate( uVunk(nzlma,MB_d)  ) ; uVunk=zero
    allocate( uVunk0(nzlma,MB_d) ) ; uVunk0=zero
  END SUBROUTINE allocate_ps_nloc2

  SUBROUTINE checkMapsBeforeForce(myrank)
    implicit none
    integer,intent(IN) :: myrank
    integer :: i
    write(5000+myrank,'(2A7)') 'nzlma','amap'
    write(5100+myrank,'(2A7)') 'nzlma','lmap'
    write(5200+myrank,'(2A7)') 'nzlma','mmap'
    write(5300+myrank,'(2A7)') 'nzlma','iorbmap'
    do i=1,nzlma
      write(5000+myrank,'(2I7)') nzlma,amap(i)
      write(5100+myrank,'(2I7)') nzlma,lmap(i)
      write(5200+myrank,'(2I7)') nzlma,mmap(i)
      write(5300+myrank,'(2I7)') nzlma,iorbmap(i)
    enddo
  END SUBROUTINE checkMapsBeforeForce

END MODULE ps_nloc2_variables
