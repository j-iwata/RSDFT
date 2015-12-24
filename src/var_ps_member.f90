MODULE var_ps_member

  use parallel_module, only: myrank

  implicit none

  PRIVATE
  PUBLIC :: allocateRps
  PUBLIC :: deallocateRps
  PUBLIC :: allocateRad1
  PUBLIC :: deallocateRad1
  PUBLIC :: send_pseudopot
  PUBLIC :: ps_allocate
  PUBLIC :: ps1d
  PUBLIC :: ps_allocate_ps1d, psg_allocate_ps1d
  PUBLIC :: ps_send_ps1d, psg_send_ps1d

  integer,PUBLIC :: Nelement_PP
  integer,PUBLIC :: Nelement_
  integer,allocatable,PUBLIC :: ippform(:)
  character(30),allocatable,PUBLIC :: file_ps(:)
  real(8),allocatable,PUBLIC :: rad(:,:),rab(:,:),rad1(:,:)
  real(8),allocatable,PUBLIC :: rabr2(:,:)
  real(8),allocatable,PUBLIC :: vql(:,:),viod(:,:,:)
  real(8),allocatable,PUBLIC :: dviod(:,:,:)
  real(8),allocatable,PUBLIC :: cdc(:,:),cdd(:,:)
  real(8),allocatable,PUBLIC :: anorm(:,:)

  real(8),allocatable,PUBLIC :: Rps(:,:)
  integer,allocatable,PUBLIC :: lo(:,:),inorm(:,:),NRps(:,:)

  integer,allocatable,PUBLIC :: Mr(:),norb(:)
  real(8),allocatable,PUBLIC :: parloc(:,:)
  real(8),allocatable,PUBLIC :: Zps(:),Zelement(:)
  real(8),allocatable,PUBLIC :: cdd_coef(:,:,:)

  integer,allocatable,PUBLIC :: NRps0(:,:)
  real(8),allocatable,PUBLIC :: Rps0(:,:)

  integer,PUBLIC :: max_psgrd=0,max_psorb=0,max_ngauss=0

  integer,allocatable,PUBLIC :: nlf(:)
  integer,allocatable,PUBLIC :: nrf(:,:)
  integer,allocatable,PUBLIC :: no(:,:)
  
  real(8),allocatable,PUBLIC :: hnml(:,:,:,:),knml(:,:,:,:)
  real(8),allocatable,PUBLIC :: hnl(:,:,:),knl(:,:,:)
  real(8),allocatable,PUBLIC :: Rcloc(:)

  integer,PUBLIC :: ps_type = 0

  type ps1d
     integer :: allocation_status = 0
     integer :: ps_type
     integer :: ippform
     character(30) :: file_ps
     real(8),allocatable :: rad(:),rab(:),rad1(:),rabr2(:)
     real(8),allocatable :: vql(:)
     real(8),allocatable :: viod(:,:), dviod(:,:)
     real(8),allocatable :: cdc(:), cdd(:)
     real(8),allocatable :: anorm(:)
     real(8),allocatable :: Rps(:)
     integer,allocatable :: lo(:), inorm(:), NRps(:)
     integer :: Mr, norb
     real(8) :: parloc(4)
     real(8) :: Zps, Zelement
     integer,allocatable :: NRps0(:)
     real(8),allocatable :: Rps0(:)
     integer :: nlf, nrf_max
     integer,allocatable :: nrf(:)
     integer,allocatable :: no(:)
     real(8),allocatable :: hnml(:,:,:),knml(:,:,:)
     real(8),allocatable :: hnl(:,:),knl(:,:)
     real(8) :: Rcloc
     integer :: ngauss
     real(8),allocatable :: cdd_coef(:,:)
! atomic pseudo wave function
     real(8),allocatable :: ups(:,:)
     real(8),allocatable :: Dij(:,:)
! uspp
     integer ::npq
     integer,allocatable :: nl3v(:)
     integer,allocatable :: l3v(:,:)
     real(8),allocatable :: ddi(:,:,:)
     real(8),allocatable :: qqr(:,:,:)
     real(8),allocatable :: qqc(:,:,:)
     real(8),allocatable :: qrL(:,:,:)
  end type ps1d

  type(ps1d),allocatable,PUBLIC :: ps(:)

  include 'mpif.h'

CONTAINS

!------------------------------------------
  SUBROUTINE allocateRps
    implicit none
    integer :: m
    if (allocated( NRps0 )) then
      call deallocateRps
    end if
    m=max_psorb
    allocate( NRps0(1:m,1:Nelement_) ) ; NRps0=0
    allocate( Rps0(1:m,1:Nelement_)  ) ; Rps0=0.d0
    NRps0(:,:)=NRps(:,:)
    Rps0(:,:)=Rps(:,:)
    return
  END SUBROUTINE allocateRps
!------------------------------------------
  SUBROUTINE deallocateRps
    implicit none
    deallocate( NRps0 )
    deallocate( Rps0  )
    return
  END SUBROUTINE deallocateRps
!------------------------------------------
  SUBROUTINE allocateRad1(m)
    implicit none
    integer,intent(IN) :: m
    if (allocated( rad1 )) then
      if (m>size(rad1(:,1))) then
        call deallocateRad1
      else
        return
      endif
    end if
    allocate( rad1(m,Nelement_)  ) ; rad1=0.d0
    return
  END SUBROUTINE allocateRad1
!------------------------------------------
  SUBROUTINE deallocateRad1
    implicit none
    deallocate( Rad1  )
    return
  END SUBROUTINE deallocateRad1

!------------------------------------------
  SUBROUTINE send_pseudopot(myrank)
    implicit none
    integer,intent(IN) :: myrank
    integer :: m,n,ierr
    m=max_psgrd
    n=max_psorb
    call mpi_bcast(m,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Nelement_PP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if ( myrank /= 0 ) then
       call ps_allocate(m,n)
    end if
    call mpi_bcast(Mr    ,size(Mr),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(norb  ,size(norb),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Zps   ,size(Zps),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(parloc,size(parloc),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(anorm ,size(anorm),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(inorm ,size(inorm),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Rps   ,size(Rps),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(NRps  ,size(NRps),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(lo    ,size(lo),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(no    ,size(no),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(vql   ,size(vql),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(cdd   ,size(cdd),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(cdc   ,size(cdc),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(rad   ,size(rad),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(rab   ,size(rab),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(viod  ,size(viod),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Rcloc ,size(Rcloc),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
! uspp
    call mpi_bcast(nlf   ,size(nlf),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nrf   ,size(nrf),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(rabr2 ,size(rabr2),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!
    call mpi_bcast(max_ngauss,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if ( max_ngauss /= 0 ) then
       if ( myrank /= 0 ) then
          allocate( cdd_coef(3,max_ngauss,Nelement_PP) ) ; cdd_coef(:,:,:)=0.0d0
       end if
       call mpi_bcast(cdd_coef,size(cdd_coef),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    end if
!
    call mpi_bcast(hnl ,size(hnl) ,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(knl ,size(knl) ,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(hnml,size(hnml),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(knml,size(knml),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(ps_type,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  END SUBROUTINE send_pseudopot

!-------------------------------------------------------
  SUBROUTINE ps_allocate(n_grd,n_orb)
    implicit none
    integer,intent(IN) :: n_grd,n_orb
    integer :: mg,mo
    real(8),allocatable :: vql_tmp(:,:),cdd_tmp(:,:),rad_tmp(:,:)
    real(8),allocatable :: cdc_tmp(:,:),viod_tmp(:,:,:)
    real(8),allocatable :: anorm_tmp(:,:),Rps_tmp(:,:),rab_tmp(:,:)
    real(8),allocatable :: rabr2_tmp(:,:)
    integer,allocatable :: inorm_tmp(:,:),lo_tmp(:,:),NRps_tmp(:,:)
    integer,allocatable :: nrf_tmp(:,:),no_tmp(:,:)
    if ( .not.allocated(hnl) ) then
       allocate( hnl(3,0:2,Nelement_PP) ) ; hnl=0.0d0
       allocate( knl(3,1:2,Nelement_PP) ) ; knl=0.0d0
    end if
    if ( .not.allocated(hnml) ) then
       allocate( hnml(3,3,0:2,Nelement_PP) ) ; hnml=0.0d0
       allocate( knml(3,3,1:2,Nelement_PP) ) ; knml=0.0d0
    end if
    if ( max_psgrd==0 .or. max_psorb==0 ) then
       allocate( Mr(Nelement_PP)   ) ; Mr(:)=0
       allocate( norb(Nelement_PP) ) ; norb(:)=0
       allocate( nlf(Nelement_PP)   ) ; nlf(:)=0
       allocate( nrf(n_orb,Nelement_PP) ) ; nrf(:,:)=0
       allocate( Zps(Nelement_PP)      ) ; Zps(:)=0.d0
       allocate( Zelement(Nelement_PP) ) ; Zelement(:)=0.d0
       allocate( parloc(4,Nelement_PP)    ) ; parloc(:,:)=0.d0
       allocate( anorm(n_orb,Nelement_PP) ) ; anorm(:,:)=0.d0
       allocate( inorm(n_orb,Nelement_PP) ) ; inorm(:,:)=0
       allocate( Rps(n_orb,Nelement_PP)   ) ; Rps(:,:)=0.d0
       allocate( NRps(n_orb,Nelement_PP)  ) ; NRps(:,:)=0
       allocate( lo(n_orb,Nelement_PP)    ) ; lo(:,:)=0
       allocate( no(n_orb,Nelement_PP) ) ; no(:,:)=0
       allocate( vql(n_grd,Nelement_PP)   ) ; vql(:,:)=0.d0
       allocate( cdd(n_grd,Nelement_PP)   ) ; cdd(:,:)=0.d0
       allocate( cdc(n_grd,Nelement_PP)   ) ; cdc(:,:)=0.d0
       allocate( rad(n_grd,Nelement_PP)   ) ; rad(:,:)=0.d0
       allocate( rab(n_grd,Nelement_PP)   ) ; rab(:,:)=0.d0
       allocate( viod(n_grd,n_orb,Nelement_PP) ) ; viod(:,:,:)=0.d0
       allocate( Rcloc(Nelement_PP) ) ; Rcloc(:)=0.0d0
       allocate( rabr2(n_grd,Nelement_PP) ) ; rabr2(:,:)=0.d0
       max_psgrd=n_grd
       max_psorb=n_orb
       return
    end if
    mg = max( max_psgrd, n_grd )
    mo = max( max_psorb, n_orb )
    if ( max_psgrd < mg ) then
       allocate( vql_tmp(mg,Nelement_PP) ) ; vql_tmp(:,:)=0.d0
       allocate( cdd_tmp(mg,Nelement_PP) ) ; cdd_tmp(:,:)=0.d0
       allocate( rad_tmp(mg,Nelement_PP) ) ; rad_tmp(:,:)=0.d0
       allocate( rab_tmp(mg,Nelement_PP) ) ; rab_tmp(:,:)=0.d0
       allocate( cdc_tmp(mg,Nelement_PP) ) ; cdc_tmp(:,:)=0.d0
       allocate( rabr2_tmp(mg,Nelement_PP) ) ; rabr2_tmp(:,:)=0.d0
       vql_tmp(1:max_psgrd,1:Nelement_PP) = vql(1:max_psgrd,1:Nelement_PP)
       cdd_tmp(1:max_psgrd,1:Nelement_PP) = cdd(1:max_psgrd,1:Nelement_PP)
       rad_tmp(1:max_psgrd,1:Nelement_PP) = rad(1:max_psgrd,1:Nelement_PP)
       rab_tmp(1:max_psgrd,1:Nelement_PP) = rab(1:max_psgrd,1:Nelement_PP)
       cdc_tmp(1:max_psgrd,1:Nelement_PP) = cdc(1:max_psgrd,1:Nelement_PP)
       rabr2_tmp(1:max_psgrd,1:Nelement_PP)=rabr2_tmp(1:max_psgrd,1:Nelement_PP)
       deallocate( rabr2 )
       deallocate( cdc )
       deallocate( rab )
       deallocate( rad )
       deallocate( cdd )
       deallocate( vql )
       allocate( vql(mg,Nelement_PP) ) ; vql(:,:)=0.d0
       allocate( cdd(mg,Nelement_PP) ) ; cdd(:,:)=0.d0
       allocate( rad(mg,Nelement_PP) ) ; rad(:,:)=0.d0
       allocate( rab(mg,Nelement_PP) ) ; rab(:,:)=0.d0
       allocate( cdc(mg,Nelement_PP) ) ; cdc(:,:)=0.d0
       allocate( rabr2(mg,Nelement_PP) ) ; rabr2(:,:)=0.d0
       vql(:,:)=vql_tmp(:,:)
       cdd(:,:)=cdd_tmp(:,:)
       rad(:,:)=rad_tmp(:,:)
       rab(:,:)=rab_tmp(:,:)
       cdc(:,:)=cdc_tmp(:,:)
       rabr2(:,:)=rabr2_tmp(:,:)
       deallocate( rabr2_tmp )
       deallocate( cdc_tmp )
       deallocate( rab_tmp )
       deallocate( rad_tmp )
       deallocate( cdd_tmp )
       deallocate( vql_tmp )
       allocate( viod_tmp(mg,mo,Nelement_PP) ) ; viod_tmp(:,:,:)=0.d0
       viod_tmp(1:max_psgrd,1:max_psorb,1:Nelement_PP) &
            = viod(1:max_psgrd,1:max_psorb,1:Nelement_PP)
       deallocate( viod )
       allocate( viod(mg,mo,Nelement_PP) ) ; viod(:,:,:)=0.d0
       viod(:,:,:)=viod_tmp(:,:,:)
       deallocate( viod_tmp )
    end if
    if ( max_psorb < mo ) then
       allocate( anorm_tmp(mo,Nelement_PP) ) ; anorm_tmp(:,:)=0.d0
       allocate( inorm_tmp(mo,Nelement_PP) ) ; inorm_tmp(:,:)=0
       allocate( lo_tmp(mo,Nelement_PP) ) ; lo_tmp(:,:)=0
       allocate( no_tmp(mo,Nelement_PP) ) ; no_tmp(:,:)=0
       allocate( Rps_tmp(mo,Nelement_PP) ) ; Rps_tmp(:,:)=0.d0
       allocate( NRps_tmp(mo,Nelement_PP) ) ; NRps_tmp(:,:)=0
       allocate( nrf_tmp(mo,Nelement_PP) ) ; nrf_tmp(:,:)=0
       anorm_tmp(1:max_psorb,1:Nelement_PP) = anorm(1:max_psorb,1:Nelement_PP)
       inorm_tmp(1:max_psorb,1:Nelement_PP) = inorm(1:max_psorb,1:Nelement_PP)
       lo_tmp(1:max_psorb,1:Nelement_PP) = lo(1:max_psorb,1:Nelement_PP)
       no_tmp(1:max_psorb,1:Nelement_PP) = no(1:max_psorb,1:Nelement_PP)
       Rps_tmp(1:max_psorb,1:Nelement_PP) = Rps(1:max_psorb,1:Nelement_PP)
       NRps_tmp(1:max_psorb,1:Nelement_PP) = NRps(1:max_psorb,1:Nelement_PP)
       nrf_tmp(1:max_psorb,1:Nelement_PP) = nrf(1:max_psorb,1:Nelement_PP)
       deallocate( nrf )
       deallocate( NRps )
       deallocate( Rps )
       deallocate( no )
       deallocate( lo )
       deallocate( inorm )
       deallocate( anorm )
       allocate( anorm(mo,Nelement_PP) ) ; anorm(:,:)=0.d0
       allocate( inorm(mo,Nelement_PP) ) ; inorm(:,:)=0
       allocate( lo(mo,Nelement_PP)    ) ; lo(:,:)=0
       allocate( no(mo,Nelement_PP)    ) ; no(:,:)=0
       allocate( Rps(mo,Nelement_PP)   ) ; Rps(:,:)=0.d0
       allocate( NRps(mo,Nelement_PP)  ) ; NRps(:,:)=0
       allocate( nrf(mo,Nelement_PP)   ) ; nrf(:,:)=0
       anorm(:,:) = anorm_tmp(:,:)
       inorm(:,:) = inorm_tmp(:,:)
       lo(:,:) = lo_tmp(:,:)
       no(:,:) = no_tmp(:,:)
       Rps(:,:) = Rps_tmp(:,:)
       NRps(:,:) = NRps_tmp(:,:)
       nrf(:,:) = nrf_tmp(:,:)
       deallocate( nrf_tmp )
       deallocate( NRps_tmp )
       deallocate( Rps_tmp )
       deallocate( lo_tmp )
       deallocate( no_tmp )
       deallocate( inorm_tmp )
       deallocate( anorm_tmp )
       if ( max_psgrd >= mg ) then
          allocate( viod_tmp(mg,mo,Nelement_PP) ) ; viod_tmp(:,:,:)=0.d0
          viod_tmp(1:max_psgrd,1:max_psorb,1:Nelement_PP) &
               = viod(1:max_psgrd,1:max_psorb,1:Nelement_PP)
          deallocate( viod )
          allocate( viod(mg,mo,Nelement_PP) ) ; viod(:,:,:)=0.d0
          viod(:,:,:)=viod_tmp(:,:,:)
          deallocate( viod_tmp )
       end if
    end if
    max_psgrd = mg
    max_psorb = mo
  END SUBROUTINE ps_allocate


  SUBROUTINE ps_allocate_ps1d( ps )
    implicit none
    type(ps1d),intent(INOUT) :: ps
    integer :: n_grd, n_orb

    if ( ps%allocation_status > 0 ) return
    ps%allocation_status = 1

    allocate( ps%hnl(3,0:2)    ) ; ps%hnl =0.0d0
    allocate( ps%knl(3,1:2)    ) ; ps%knl =0.0d0
    allocate( ps%hnml(3,3,0:2) ) ; ps%hnml=0.0d0
    allocate( ps%knml(3,3,1:2) ) ; ps%knml=0.0d0

    n_grd = ps%Mr
    n_orb = ps%norb

    ps%nlf      = 0
    ps%Zps      = 0.0d0
    ps%Zelement = 0.0d0
    ps%parloc   = 0.0d0
    ps%Rcloc    = 0.0d0  

    if ( n_orb /= 0 ) then
       allocate( ps%nrf(n_orb)        ) ; ps%nrf=0
       allocate( ps%anorm(n_orb)      ) ; ps%anorm=0.0d0
       allocate( ps%inorm(n_orb)      ) ; ps%inorm=0
       allocate( ps%Rps(n_orb)        ) ; ps%Rps=0.0d0
       allocate( ps%NRps(n_orb)       ) ; ps%NRps=0
       allocate( ps%lo(n_orb)         ) ; ps%lo=0
       allocate( ps%no(n_orb)         ) ; ps%no=0
       allocate( ps%Dij(n_orb,n_orb)  ) ; ps%Dij=0.0d0
    end if

    if ( n_grd /= 0 ) then
       allocate( ps%vql(n_grd)        ) ; ps%vql=0.0d0
       allocate( ps%cdc(n_grd)        ) ; ps%cdc=0.0d0
       allocate( ps%cdd(n_grd)        ) ; ps%cdd=0.0d0
       allocate( ps%rad(n_grd)        ) ; ps%rad=0.0d0
       allocate( ps%rab(n_grd)        ) ; ps%rab=0.0d0
       allocate( ps%rabr2(n_grd)      ) ; ps%rabr2=0.0d0
    end if

    if ( n_grd /=0 .and. n_orb /= 0 ) then
       allocate( ps%viod(n_grd,n_orb) ) ; ps%viod=0.0d0
       allocate( ps%ups(n_grd,n_orb)  ) ; ps%ups=0.0d0
    end if

  END SUBROUTINE ps_allocate_ps1d


  SUBROUTINE psg_allocate_ps1d( ps )
    implicit none
    type(ps1d),intent(INOUT) :: ps
    integer :: n_g,n_k,n_l,n_r
    if ( ps%allocation_status > 1 ) return
    ps%allocation_status = 2
    n_g=ps%Mr
    n_k=ps%npq
    n_l=ps%nlf
    n_r=ps%nrf_max
    allocate( ps%ddi(n_r,n_r,n_l) ) ; ps%ddi(:,:,:)=0.0d0
    allocate( ps%qqr(n_r,n_r,n_l) ) ; ps%qqr(:,:,:)=0.0d0
    allocate( ps%qqc(n_r,n_r,n_l) ) ; ps%qqc(:,:,:)=0.0d0
    allocate( ps%nl3v(n_k)        ) ; ps%nl3v(:)=0
    allocate( ps%l3v(n_l,n_k)     ) ; ps%l3v(:,:)=0
    allocate( ps%qrL(n_g,n_l,n_k) ) ; ps%qrL(:,:,:)=0.0d0
  END SUBROUTINE psg_allocate_ps1d


  SUBROUTINE ps_send_ps1d( ps )
    implicit none
    type(ps1d),intent(INOUT) :: ps
    integer :: i

    call mpi_bcast(ps%Mr,1,mpi_integer,0,mpi_comm_world,i)
    call mpi_bcast(ps%norb,1,mpi_integer,0,mpi_comm_world,i)

    call ps_allocate_ps1d( ps )

    call mpi_bcast(ps%hnl,size(ps%hnl),mpi_real8,0,mpi_comm_world,i)
    call mpi_bcast(ps%knl,size(ps%knl),mpi_real8,0,mpi_comm_world,i)
    call mpi_bcast(ps%hnml,size(ps%hnml),mpi_real8,0,mpi_comm_world,i)
    call mpi_bcast(ps%knml,size(ps%knml),mpi_real8,0,mpi_comm_world,i)

    call mpi_bcast(ps%nlf,1,mpi_integer,0,mpi_comm_world,i)
    call mpi_bcast(ps%Zps,1,mpi_real8,0,mpi_comm_world,i)
    call mpi_bcast(ps%Zelement,1,mpi_real8,0,mpi_comm_world,i)
    call mpi_bcast(ps%parloc,size(ps%parloc),mpi_real8,0,mpi_comm_world,i)
    call mpi_bcast(ps%Rcloc,1,mpi_real8,0,mpi_comm_world,i)

    if ( ps%norb > 0 ) then
       call mpi_bcast(ps%nrf,size(ps%nrf),mpi_integer,0,mpi_comm_world,i)
       call mpi_bcast(ps%anorm,size(ps%anorm),mpi_real8,0,mpi_comm_world,i)
       call mpi_bcast(ps%inorm,size(ps%inorm),mpi_integer,0,mpi_comm_world,i)
       call mpi_bcast(ps%Rps,size(ps%Rps),mpi_real8,0,mpi_comm_world,i)
       call mpi_bcast(ps%NRps,size(ps%NRps),mpi_integer,0,mpi_comm_world,i)
       call mpi_bcast(ps%lo,size(ps%lo),mpi_integer,0,mpi_comm_world,i)
       call mpi_bcast(ps%no,size(ps%no),mpi_integer,0,mpi_comm_world,i)
       call mpi_bcast(ps%Dij,size(ps%Dij),mpi_real8,0,mpi_comm_world,i)
    end if

    if ( ps%Mr > 0 ) then
       call mpi_bcast(ps%vql,size(ps%vql),mpi_real8,0,mpi_comm_world,i)
       call mpi_bcast(ps%cdc,size(ps%cdc),mpi_real8,0,mpi_comm_world,i)
       call mpi_bcast(ps%cdd,size(ps%cdd),mpi_real8,0,mpi_comm_world,i)
       call mpi_bcast(ps%rad,size(ps%rad),mpi_real8,0,mpi_comm_world,i)
       call mpi_bcast(ps%rab,size(ps%rab),mpi_real8,0,mpi_comm_world,i)
    end if

    if ( ps%Mr > 0 .and. ps%norb > 0 ) then
       call mpi_bcast(ps%viod,size(ps%viod),mpi_real8,0,mpi_comm_world,i)
       call mpi_bcast(ps%ups,size(ps%ups),mpi_real8,0,mpi_comm_world,i)
    end if

  END SUBROUTINE ps_send_ps1d


  SUBROUTINE psg_send_ps1d( ps )
    implicit none
    type(ps1d),intent(INOUT) :: ps
    integer :: i

    call mpi_bcast(ps%npq,1,mpi_integer,0,mpi_comm_world,i)
    call mpi_bcast(ps%nrf_max,1,mpi_integer,0,mpi_comm_world,i)

    call psg_allocate_ps1d( ps )

    call mpi_bcast(ps%ddi,size(ps%ddi),mpi_real8,0,mpi_comm_world,i)
    call mpi_bcast(ps%qqr,size(ps%qqr),mpi_real8,0,mpi_comm_world,i)
    call mpi_bcast(ps%qqc,size(ps%qqc),mpi_real8,0,mpi_comm_world,i)
    call mpi_bcast(ps%nl3v,size(ps%nl3v),mpi_integer,0,mpi_comm_world,i)
    call mpi_bcast(ps%l3v,size(ps%l3v),mpi_integer,0,mpi_comm_world,i)
    call mpi_bcast(ps%qrL,size(ps%qrL),mpi_real8,0,mpi_comm_world,i)    

  END SUBROUTINE psg_send_ps1d


END MODULE var_ps_member
