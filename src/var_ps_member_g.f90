MODULE var_ps_member_g

  use parallel_module, only: myrank

  implicit none

  PRIVATE
  PUBLIC :: allocateNzqr
  PUBLIC :: deallocateNzqr
  PUBLIC :: allocateKtoK
  PUBLIC :: deallocateKtoK
  PUBLIC :: allocateQRps
  PUBLIC :: deallocateQRps
  PUBLIC :: sendPSG
  PUBLIC :: allocatePSG

  integer,allocatable,PUBLIC :: nl3v(:,:)
  integer,allocatable,PUBLIC :: l3v(:,:,:)
  
  real(8),allocatable,PUBLIC :: qrL(:,:,:,:)
  integer,PUBLIC :: maxcJ
  real(8),allocatable,PUBLIC :: dqrL(:,:,:,:,:)

  complex(8),allocatable,PUBLIC :: QG(:,:,:)

  real(8),allocatable,PUBLIC :: QRij(:,:)
  
  real(8),allocatable,PUBLIC :: ddi(:,:,:,:)

  real(8),allocatable,PUBLIC :: qqr(:,:,:,:)
  real(8),allocatable,PUBLIC :: qqc(:,:,:,:)

  integer,PUBLIC :: k1max,k2max,k3max,lpsmax
  
  integer,allocatable,PUBLIC :: k1_to_k2(:,:)
  integer,allocatable,PUBLIC :: k1_to_k3(:,:)
  integer,allocatable,PUBLIC :: k1_to_iorb(:,:,:)
  integer,allocatable,PUBLIC :: k1_to_l(:,:,:)
  integer,allocatable,PUBLIC :: k1_to_m(:,:,:)
  integer,allocatable,PUBLIC :: N_k1(:)

  integer,allocatable,PUBLIC :: k2_to_iorb(:,:,:)
  integer,allocatable,PUBLIC :: N_k2(:)

  integer,allocatable,PUBLIC :: Q_NRps(:,:)
  real(8),allocatable,PUBLIC :: Q_Rps(:,:)

  integer,allocatable,PUBLIC :: npq(:)

  integer,PUBLIC :: N_nzqr
  integer,PUBLIC :: N_nlop

  integer,allocatable,PUBLIC :: nzqr_pair(:,:)
  integer,allocatable,PUBLIC :: atommap(:),kk1map(:,:),k1map(:)
  integer,allocatable,PUBLIC :: nlop_pair(:,:)
  real(8),allocatable,PUBLIC :: Dij(:,:),Dij00(:)
  real(8),allocatable,PUBLIC :: qij(:),qij_f(:)

  integer,PUBLIC :: max_Rref=0,max_Lref=0,max_k2=0,max_qgrd=0

#ifdef _DRSDFT_
  real(8),allocatable,PUBLIC :: uVunk(:,:),uVunk0(:,:)
#else
  complex(8),allocatable,PUBLIC :: uVunk(:,:),uVunk0(:,:)
#endif


CONTAINS

!------------------------------------------
  SUBROUTINE allocateNzqr(nspin,Natom,k1max)
    implicit none
    integer,intent(IN) :: nspin,Natom,k1max
    call deallocateNzqr
    allocate(nzqr_pair(N_nzqr,2)) ; nzqr_pair=0
    allocate(atommap(N_nzqr)    ) ; atommap=0
    allocate(k1map(N_nzqr)      ) ; k1map=0
    allocate(kk1map(k1max,Natom)) ; kk1map=0
    allocate(nlop_pair(2,N_nlop)) ; nlop_pair=0
    allocate(Dij(N_nzqr,nspin)  ) ; Dij=0.d0
    allocate(Dij00(N_nzqr)      ) ; Dij00=0.d0
    allocate(qij(N_nzqr)        ) ; qij=0.d0
    allocate(qij_f(N_nzqr)      ) ; qij_f=0.d0
    return
  END SUBROUTINE allocateNzqr

!------------------------------------------
  SUBROUTINE deallocateNzqr
    implicit none
    if (allocated(nzqr_pair)) deallocate(nzqr_pair)
    if (allocated(atommap)  ) deallocate(atommap)
    if (allocated(k1map)    ) deallocate(k1map)
    if (allocated(kk1map)   ) deallocate(kk1map)
    if (allocated(nlop_pair)) deallocate(nlop_pair)
    if (allocated(Dij)      ) deallocate(Dij)
    if (allocated(Dij00)    ) deallocate(Dij00)
    if (allocated(qij)      ) deallocate(qij)
    if (allocated(qij_f)    ) deallocate(qij_f)
    return
  END SUBROUTINE deallocateNzqr
  
!------------------------------------------
  SUBROUTINE allocateKtoK( k1max,k2max,nki,Rrefmax,Lrefmax )
    implicit none
    integer,intent(IN) :: k2max,nki,Rrefmax,Lrefmax
    integer,intent(INOUT) :: k1max
    
    if ( allocated(k1_to_k2) ) then
      call deallocateKtoK
    end if

    k1max  = (Rrefmax*(Lrefmax**2))*(Rrefmax*(Lrefmax**2)+1)/2
    lpsmax = max_Rref*max_Lref
    k3max  = (max_Lref**2)*(max_Lref**2+1)/2

    allocate( k1_to_k2(k1max,nki)     ) ; k1_to_k2(:,:)=0
    allocate( k1_to_k3(k1max,nki)     ) ; k1_to_k3(:,:)=0
    allocate( k1_to_iorb(2,k1max,nki) ) ; k1_to_iorb(:,:,:)=0
    allocate( N_k1(nki)               ) ; N_k1(:)=0

    allocate( k2_to_iorb(2,k2max,nki) ) ; k2_to_iorb(:,:,:)=0
    allocate( N_k2(nki)               ) ; N_k2(:)=0

    allocate( k1_to_l(2,k1max,nki) ) ; k1_to_l=0
    allocate( k1_to_m(2,k1max,nki) ) ; k1_to_m=0

    return
  END SUBROUTINE allocateKtoK

!------------------------------------------
  SUBROUTINE deallocateKtoK
    implicit none
    
    deallocate( k1_to_k2 )
    deallocate( k1_to_k3 )
    deallocate( k1_to_iorb )
    deallocate( N_k1 )

    deallocate( k2_to_iorb )
    deallocate( N_k2 )

    return
  END SUBROUTINE deallocateKtoK

!------------------------------------------

  SUBROUTINE allocateQRps( k2max,nki )
    implicit none
    integer,intent(IN) :: k2max,nki

    if ( allocated(Q_NRps) ) then
      call deallocateQRps
    end if

    allocate( Q_NRps(k2max,nki) ) ; Q_NRps(:,:)=0
    allocate( Q_Rps(k2max,nki)  ) ; Q_Rps(:,:)=0.d0

    return
  END SUBROUTINE allocateQRps

!------------------------------------------
  SUBROUTINE deallocateQRps
    implicit none

    deallocate( Q_NRps )
    deallocate( Q_Rps  )

    return
  END SUBROUTINE deallocateQRps

!------------------------------------------
  SUBROUTINE sendPSG(myrank,Nelement_PP)
    implicit none
    integer,intent(IN) :: myrank,Nelement_PP
    integer :: maxs(1:4),ierr
    include 'mpif.h'
    maxs(1)=max_Lref
    maxs(2)=max_Rref
    maxs(3)=max_k2
    maxs(4)=max_qgrd
    call mpi_bcast(maxs,4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Nelement_PP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if ( myrank /= 0 ) then
       call allocatePSG(maxs(1),maxs(2),maxs(3),maxs(4),Nelement_PP)
    end if
    call mpi_bcast(npq ,size(npq) ,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(ddi ,size(ddi) ,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(qqr ,size(qqr) ,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nl3v,size(nl3v),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(l3v ,size(l3v) ,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(qrL ,size(qrL) ,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    return
  END SUBROUTINE sendPSG

!------------------------------------------
  SUBROUTINE allocatePSG( n_l,n_r,n_k,n_g,nki )
    implicit none
    integer,intent(IN) :: n_l,n_r,n_k,n_g,nki
    real(8),allocatable :: ddi_(:,:,:,:),qqr_(:,:,:,:)
    integer,allocatable :: nl3v_(:,:),l3v_(:,:,:)
    real(8),allocatable :: qrL_(:,:,:,:)
    integer :: ml,mr,mk,mg

    if ( max_Rref==0 .or. max_Lref==0 .or. max_k2==0 ) then
        allocate( npq(nki) ) ; npq(:)=0
        allocate( ddi(n_r,n_r,n_l,nki) ) ; ddi(:,:,:,:)=0.d0
        allocate( qqr(n_r,n_r,n_l,nki) ) ; qqr(:,:,:,:)=0.d0
        allocate( qqc(n_r,n_r,n_l,nki) ) ; qqc(:,:,:,:)=0.d0
        allocate( nl3v(n_k,nki) ) ; nl3v(:,:)=0
        allocate( l3v(n_l,n_k,nki) ) ; l3v(:,:,:)=0
        allocate( qrL(n_g,n_l,n_k,nki) ) ; qrL(:,:,:,:)=0.d0
        max_Rref=n_r
        max_Lref=n_l
        max_k2=n_k
        max_qgrd=n_g
        return
    end if
    ml=max(max_Lref,n_l)
    mr=max(max_Rref,n_r)
    mk=max(max_k2,n_k)
    mg=max(max_qgrd,n_g)
    if ( max_Rref<mr ) then
        allocate( ddi_(mr,mr,ml,nki) ) ; ddi_(:,:,:,:)=0.d0
        allocate( qqr_(mr,mr,ml,nki) ) ; qqr_(:,:,:,:)=0.d0
        ddi_(1:max_Rref,1:max_Rref,1:max_Lref,1:nki)=ddi(1:max_Rref,1:max_Rref,1:max_Lref,1:nki)
        qqr_(1:max_Rref,1:max_Rref,1:max_Lref,1:nki)=qqr(1:max_Rref,1:max_Rref,1:max_Lref,1:nki)
        deallocate( ddi,qqr,qqc )
        allocate( ddi(mr,mr,ml,nki) ) ; ddi(:,:,:,:)=0.d0
        allocate( qqr(mr,mr,ml,nki) ) ; qqr(:,:,:,:)=0.d0
        allocate( qqc(mr,mr,ml,nki) ) ; qqc(:,:,:,:)=0.d0
        ddi(:,:,:,:)=ddi_(:,:,:,:)
        qqr(:,:,:,:)=qqr_(:,:,:,:)
        deallocate( ddi_,qqr_ )
    end if
    if ( max_k2<mk ) then
        allocate( nl3v_(mk,nki) ) ; nl3v_(:,:)=0
        allocate( l3v_(ml,mk,nki) ) ; l3v_(:,:,:)=0
        allocate( qrL_(mg,ml,mk,nki) ) ; qrL_(:,:,:,:)=0.d0
        nl3v_(1:max_k2,1:nki)=nl3v(1:max_k2,1:nki)
        l3v_(1:max_Lref,1:max_k2,1:nki)=l3v(1:max_Lref,1:max_k2,1:nki)
        qrL_(1:max_qgrd,1:max_Lref,1:max_k2,1:nki)=qrL(1:max_qgrd,1:max_Lref,1:max_k2,1:nki)
        deallocate( nl3v,l3v,qrL )
        allocate( nl3v(mk,nki) ) ; nl3v(:,:)=0
        allocate( l3v(ml,mk,nki) ) ; l3v(:,:,:)=0
        allocate( qrL(mg,ml,mk,nki) ) ; qrL(:,:,:,:)=0.d0
        nl3v(:,:)=nl3v_(:,:)
        l3v(:,:,:)=l3v_(:,:,:)
        qrL(:,:,:,:)=qrL_(:,:,:,:)
        deallocate( nl3v_,l3v_,qrL_ )
    end if
    if ( max_Lref<ml ) then
        if ( max_k2>=mk ) then
            allocate( l3v_(ml,mk,nki) ) ; l3v_(:,:,:)=0
            allocate( qrL_(mg,ml,mk,nki) ) ; qrL_(:,:,:,:)=0.d0
            l3v_(1:max_Lref,1:max_k2,1:nki)=l3v(1:max_Lref,1:max_k2,1:nki)
            qrL_(1:max_qgrd,1:max_Lref,1:max_k2,1:nki)=qrL(1:max_qgrd,1:max_Lref,1:max_k2,1:nki)
            deallocate( l3v,qrL )
            allocate( l3v(ml,mk,nki) ) ; l3v(:,:,:)=0
            allocate( qrL(mg,ml,mk,nki) ) ; qrL(:,:,:,:)=0.d0
            l3v(:,:,:)=l3v_(:,:,:)
            qrL(:,:,:,:)=qrL_(:,:,:,:)
            deallocate( l3v_,qrL_ )
        end if
        if ( max_Rref>=mr ) then
            allocate( ddi_(mr,mr,ml,nki) ) ; ddi_(:,:,:,:)=0.d0
            allocate( qqr_(mr,mr,ml,nki) ) ; qqr_(:,:,:,:)=0.d0
            ddi_(1:max_Rref,1:max_Rref,1:max_Lref,1:nki)=ddi(1:max_Rref,1:max_Rref,1:max_Lref,1:nki)
            qqr_(1:max_Rref,1:max_Rref,1:max_Lref,1:nki)=qqr(1:max_Rref,1:max_Rref,1:max_Lref,1:nki)
            deallocate( ddi,qqr,qqc )
            allocate( ddi(mr,mr,ml,nki) ) ; ddi(:,:,:,:)=0.d0
            allocate( qqr(mr,mr,ml,nki) ) ; qqr(:,:,:,:)=0.d0
            allocate( qqc(mr,mr,ml,nki) ) ; qqc(:,:,:,:)=0.d0
            ddi(:,:,:,:)=ddi_(:,:,:,:)
            qqr(:,:,:,:)=qqr_(:,:,:,:)
            deallocate( ddi_,qqr_ )
        end if
    end if
    if ( max_qgrd<mg ) then
        if ( max_Lref>=ml .and. max_k2>=mk ) then
            allocate( qrL_(mg,ml,mk,nki) ) ; qrL_(:,:,:,:)=0.d0
            qrL_(1:max_qgrd,1:max_Lref,1:max_k2,1:nki)=qrL(1:max_qgrd,1:max_Lref,1:max_k2,1:nki)
            deallocate( qrL )
            allocate( qrL(mg,ml,mk,nki) ) ; qrL(:,:,:,:)=0.d0
            qrL(:,:,:,:)=qrL_(:,:,:,:)
            deallocate( qrL_ )
        end if
    end if
    max_Lref=ml
    max_Rref=mr
    max_k2=mk
    max_qgrd=mg

    return
  END SUBROUTINE allocatePSG


END MODULE var_ps_member_g
