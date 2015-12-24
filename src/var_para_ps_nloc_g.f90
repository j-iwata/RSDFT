MODULE var_para_ps_nloc_g

  implicit none

  integer :: Mqr
  integer :: c_nzqr
  integer :: MMJJ_Q,MMJJ_t_Q
  integer :: MAXMJJ_Q, MAXMJJ_MAP_Q
  integer :: nrqr_xyz(6)

  integer,allocatable :: JJP_Q(:,:)
  integer,allocatable :: JJ_MAP_Q(:,:,:)
  integer,allocatable :: MJJ_MAP_Q(:)
  integer,allocatable :: MJJ_Q(:)
  integer,allocatable :: nl_rank_map_Q(:)

  integer :: nl_max_send_Q
  integer,allocatable :: qr_nsend(:)
  integer,allocatable :: sendmap_Q(:,:),recvmap_Q(:,:)

  integer,allocatable :: num_2_rank_Q(:,:)

  integer,allocatable :: amap_Q(:),k1map_Q(:),lmamap_Q(:,:)

#ifdef _DRSDFT_
  real(8),allocatable :: sbufnl_Q(:,:),rbufnl_Q(:,:)
  real(8),parameter :: zero=0.d0
#else
  complex(8),allocatable :: sbufnl_Q(:,:),rbufnl_Q(:,:)
  complex(8),parameter :: zero=(0.d0,0.d0)
#endif

CONTAINS

  SUBROUTINE allocateJJMAPQ
    implicit none
    call deallocateJJMAPQ
    allocate(JJ_MAP_Q(6,MMJJ_Q,c_nzqr)) ; JJ_MAP_Q      =0
    allocate(MJJ_MAP_Q(c_nzqr)        ) ; MJJ_MAP_Q     =0
    allocate(MJJ_Q(c_nzqr)            ) ; MJJ_Q         =0
    allocate(nl_rank_map_Q(c_nzqr)    ) ; nl_rank_map_Q =0
    allocate(amap_Q(c_nzqr)           ) ; amap_Q        =0
    allocate(k1map_Q(c_nzqr)          ) ; k1map_Q       =0
    allocate(lmamap_Q(c_nzqr,2)       ) ; lmamap_Q      =0
    return
  END SUBROUTINE allocateJJMAPQ

  SUBROUTINE deallocateJJMAPQ
    implicit none
    if (allocated(JJ_MAP_Q) )     deallocate(JJ_MAP_Q)
    if (allocated(MJJ_MAP_Q))     deallocate(MJJ_MAP_Q)
    if (allocated(MJJ_Q)    )     deallocate(MJJ_Q)
    if (allocated(nl_rank_map_Q)) deallocate(nl_rank_map_Q)
    if (allocated(amap_Q)   ) deallocate(amap_Q)
    if (allocated(k1map_Q)  ) deallocate(k1map_Q)
    if (allocated(lmamap_Q) ) deallocate(lmamap_Q)
    return
  END SUBROUTINE deallocateJJMAPQ

  SUBROUTINE allocateMAPQ(nl_max_send,nprocs_g)
    implicit none
    integer,intent(IN) :: nl_max_send,nprocs_g
    integer :: n
    call deallocateMAPQ
    n=nl_max_send
    allocate(qr_nsend(0:nprocs_g-1)   ) ; qr_nsend  =0
    allocate(sendmap_Q(n,0:nprocs_g-1)) ; sendmap_Q =0
    allocate(recvmap_Q(n,0:nprocs_g-1)) ; recvmap_Q =0
    return
  END SUBROUTINE allocateMAPQ

  SUBROUTINE deallocateMAPQ
    implicit none
    if (allocated(qr_nsend))  deallocate(qr_nsend)
    if (allocated(sendmap_Q)) deallocate(sendmap_Q)
    if (allocated(recvmap_Q)) deallocate(recvmap_Q)
    return
  END SUBROUTINE deallocateMAPQ

  SUBROUTINE allocateBufQ(max_qr_nsend,nprocs_g,MB_d)
    implicit none
    integer,intent(IN) :: max_qr_nsend,nprocs_g,MB_d
    integer :: n
    call deallocateBufQ
    n=max_qr_nsend*MB_d*3
    allocate(rbufnl_Q(n,0:nprocs_g-1)) ; rbufnl_Q=zero
    allocate(sbufnl_Q(n,0:nprocs_g-1)) ; sbufnl_Q=zero
    return
  END SUBROUTINE allocateBufQ

  SUBROUTINE deallocateBufQ
    implicit none
    if (allocated(rbufnl_Q)) deallocate(rbufnl_Q)
    if (allocated(sbufnl_Q)) deallocate(sbufnl_Q)
    return
  END SUBROUTINE deallocateBufQ

END MODULE var_para_ps_nloc_g
