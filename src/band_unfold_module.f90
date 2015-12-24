MODULE band_unfold_module

  use aa_module, only: aa
  use bb_module, only: bb, calc_bb
  use rgrid_module, only: Ngrid, Igrid
  use band_variables, only: nbk, ak, nfki
  use parallel_module
  use wf_module, only: unk, esp
  use fft_module

  implicit none

  PRIVATE
  PUBLIC :: read_band_unfold, init_band_unfold, band_unfold &
           ,finalize_band_unfold, iswitch_banduf

  logical :: iswitch_banduf=.false.

  integer :: unit_uf
  integer :: nktrj,nktrj_0
  integer :: nbk_pc
  integer :: mg_pc, mg_sc
  integer,allocatable :: LG_pc(:,:)
  integer,allocatable :: map_ktrj(:,:)

  real(8) :: ax_pc, aa_pc(3,3), bb_pc(3,3)
  real(8),allocatable :: kbb_pc(:,:),kxyz_pc(:,:)
  real(8),allocatable :: kbb_0(:,:),kxyz_0(:,:),kbb(:,:)
  real(8),allocatable :: weight_uf(:,:,:)

#ifdef _DRSDFT_
  integer,parameter :: TYPE_MAIN=MPI_REAL8
  real(8),parameter :: zero=0.0d0
#else
  integer,parameter :: TYPE_MAIN=MPI_COMPLEX16
  complex(8),parameter :: zero=(0.0d0,0.0d0)
#endif

CONTAINS


  SUBROUTINE read_band_unfold( rank, unit )
    implicit none
    integer,intent(IN) :: rank,unit
    character(6) :: cbuf,ckey
    integer :: i
    call write_border( 0, " read_band_unfold(start)" )
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:6) == "BANDUF" ) then
             backspace(unit)
             read(unit,*) cbuf
             read(unit,*) ax_pc
             read(unit,*) aa_pc(1:3,1)
             read(unit,*) aa_pc(1:3,2)
             read(unit,*) aa_pc(1:3,3)
             iswitch_banduf=.true.
          end if
       end do
999    continue
       write(*,*) "iswitch_banduf=",iswitch_banduf
    end if
    call mpi_bcast(iswitch_banduf,1,MPI_LOGICAL,0,MPI_COMM_WORLD,i)
    call mpi_bcast(ax_pc,1,MPI_REAL8,0,MPI_COMM_WORLD,i)
    call mpi_bcast(aa_pc,9,MPI_REAL8,0,MPI_COMM_WORLD,i)
    call write_border( 0, " read_band_unfold(end)" )
  END SUBROUTINE read_band_unfold


  SUBROUTINE init_band_unfold( nktrj_io, ktrj_out, unit_in, disp_switch )

    implicit none
    integer,intent(IN)  :: unit_in
    logical,intent(IN)  :: disp_switch
    integer,intent(INOUT) :: nktrj_io
    real(8),intent(OUT) :: ktrj_out(6,nktrj_io)
    integer :: iktrj,MB,MS

    call write_border( 0, " init_band_unfold(start)" )

    unit_uf = unit_in
    if ( myrank == 0 ) open(unit_uf,file="band_ufld")

    aa_pc(:,:) = ax_pc*aa_pc(:,:)

    call calc_bb( aa_pc, bb_pc )

    nktrj_0 = nktrj_io

    call init0_band_unfold

    ktrj_out=0.0d0
    ktrj_out(1:3,1:nktrj) = kbb(1:3,1:nktrj)

    nktrj_io = nktrj

    call init1_band_unfold

    MB=sum(ir_band)
    MS=sum(ir_spin)
    allocate( weight_uf(nktrj_0,MB,MS) ) ; weight_uf=0.0d0

    if ( disp_switch ) then
       write(*,*) "size[weight_uf](MB)=",size(weight_uf)*8.0d0/1024.d0**2
       write(*,*) "mg_pc=",mg_pc
       write(*,*) "(aa_pc)"
       write(*,'(1x, f20.15)') ax_pc
       write(*,'(1x,3f20.15)') aa_pc(1:3,1) 
       write(*,'(1x,3f20.15)') aa_pc(1:3,2) 
       write(*,'(1x,3f20.15)') aa_pc(1:3,3) 
       write(*,*) "(bb_pc)"
       write(*,'(1x,3f20.15)') bb_pc(1:3,1) 
       write(*,'(1x,3f20.15)') bb_pc(1:3,2) 
       write(*,'(1x,3f20.15)') bb_pc(1:3,3) 
       write(*,'(1x,a4,1x,2(10x,a10,10x))') "iktrj","k(SC)","k(PC)"
       do iktrj=1,nktrj_0
          write(*,'(1x,i4,1x,2(3f10.6,1x))') &
               iktrj,kbb_0(1:3,iktrj),kbb_pc(1:3,iktrj)
       end do
       write(*,'(1x,a4,1x,10x,a10,10x,a4)') "iktrj","k'(SC)"," "
       do iktrj=1,nktrj
          write(*,'(1x,i4,1x,3f10.6,1x,i4)') &
               iktrj,kbb(1:3,iktrj),count(map_ktrj(:,1)==iktrj)
       end do
    end if

    call write_border( 0, " init_band_unfold(end)" )

  END SUBROUTINE init_band_unfold


  SUBROUTINE init0_band_unfold

    implicit none
    real(8) :: pi2
    integer :: ibk,j,i,iktrj,loop

! --- k points in PC ---

    allocate( kbb_pc(3,nktrj_0)  ) ; kbb_pc=0.0d0
    allocate( kxyz_pc(3,nktrj_0) ) ; kxyz_pc=0.0d0

    j=0
    do ibk=1,nbk
       do i=0,nfki(ibk)-1
          j=j+1
          kbb_pc(1,j) = ak(1,ibk) + i*( ak(1,ibk+1) - ak(1,ibk) )/nfki(ibk)
          kbb_pc(2,j) = ak(2,ibk) + i*( ak(2,ibk+1) - ak(2,ibk) )/nfki(ibk)
          kbb_pc(3,j) = ak(3,ibk) + i*( ak(3,ibk+1) - ak(3,ibk) )/nfki(ibk)
       end do ! i
    end do ! ibk
    kbb_pc(1,nktrj_0) = ak(1,nbk+1)
    kbb_pc(2,nktrj_0) = ak(2,nbk+1)
    kbb_pc(3,nktrj_0) = ak(3,nbk+1)

    do iktrj=1,nktrj_0
       kxyz_pc(1,iktrj) = sum( bb_pc(1,:)*kbb_pc(:,iktrj) )
       kxyz_pc(2,iktrj) = sum( bb_pc(2,:)*kbb_pc(:,iktrj) )
       kxyz_pc(3,iktrj) = sum( bb_pc(3,:)*kbb_pc(:,iktrj) )
    end do

! --- corresponding k-points in SC ---

    pi2 = 2.0d0*acos(-1.0d0)

    allocate( kbb_0(3,nktrj_0) ) ; kbb_0=0.0d0
    allocate( kxyz_0(3,nktrj_0)  ) ; kxyz_0=0.0d0

    do iktrj=1,nktrj_0
       kbb_0(1,iktrj) = sum( aa(:,1)*kxyz_pc(:,iktrj) )/pi2
       kbb_0(2,iktrj) = sum( aa(:,2)*kxyz_pc(:,iktrj) )/pi2
       kbb_0(3,iktrj) = sum( aa(:,3)*kxyz_pc(:,iktrj) )/pi2
    end do

    do iktrj=1,nktrj_0
       do i=1,3
          do loop=1,1000
             if ( kbb_0(i,iktrj) > 0.5d0 ) then
                kbb_0(i,iktrj) = kbb_0(i,iktrj) - 1.0d0
             else if ( kbb_0(i,iktrj) < -0.5d0 ) then
                kbb_0(i,iktrj) = kbb_0(i,iktrj) + 1.0d0
             else
                exit
             end if
          end do ! loop
       end do ! i
       kxyz_0(1:3,iktrj) = matmul( bb(1:3,1:3), kbb_0(1:3,iktrj) )
    end do ! iktrj

! ---

    allocate( kbb(3,nktrj_0)      ) ; kbb=0.0d0
    allocate( map_ktrj(nktrj_0,2) ) ; map_ktrj=0

    nktrj=1
    kbb(1:3,1)=kbb_0(1:3,1)
    map_ktrj(1,1)=1
    map_ktrj(1,2)=0

    do i=2,nktrj_0
       do j=1,nktrj
          if ( all( abs(kbb_0(:,i)-kbb(:,j)) < 1.d-12 ) ) then
             map_ktrj(i,1)=j
             map_ktrj(i,2)=0
             exit
          else if ( all( abs(kbb_0(:,i)+kbb(:,j)) < 1.d-12 ) ) then
             map_ktrj(i,1)=j
             map_ktrj(i,2)=1
             exit
          end if
       end do
       if ( j > nktrj ) then
          nktrj=nktrj+1
          kbb(:,nktrj)=kbb_0(:,i)
          map_ktrj(i,1)=j
          map_ktrj(i,2)=0
       end if
    end do

  END SUBROUTINE init0_band_unfold


  SUBROUTINE init1_band_unfold

    implicit none
    integer :: i,j,i1,i2,i3
    integer,allocatable :: LG_tmp(:,:)
    real(8) :: pi2
    real(8) :: matab(3,3),vtmp(3),utmp(3)

    pi2 = 2.0d0*acos(-1.0d0)

    do j=1,3
    do i=1,3
       matab(i,j)=sum( aa_pc(1:3,i)*bb(1:3,j) )/pi2
    end do
    end do
!    do i=1,3
!       write(*,'(1x,3f15.10)') matab(i,1:3)
!    end do

    mg_sc = Ngrid(1)*Ngrid(2)*Ngrid(3)

    allocate( LG_tmp(3,mg_sc) ) ; LG_tmp=0

    mg_pc=0
    do i3=-(Ngrid(3)-1)/2,Ngrid(3)/2
    do i2=-(Ngrid(2)-1)/2,Ngrid(2)/2
    do i1=-(Ngrid(1)-1)/2,Ngrid(1)/2

       vtmp(1) = i1
       vtmp(2) = i2
       vtmp(3) = i3

       utmp(:) = matmul( matab(:,:), vtmp(:) )

       vtmp(:) = nint( utmp(:) )

       if ( all( abs(vtmp-utmp) < 1.d-13 ) ) then
          mg_pc=mg_pc+1
          LG_tmp(1,mg_pc)=i1
          LG_tmp(2,mg_pc)=i2
          LG_tmp(3,mg_pc)=i3
       end if

    end do
    end do
    end do

    allocate( LG_pc(3,mg_pc) ) ; LG_pc=0
    LG_pc(:,:) = LG_tmp(:,1:mg_pc)

    deallocate( LG_tmp )

  END SUBROUTINE init1_band_unfold


  SUBROUTINE band_unfold( jktrj, disp_switch )

    implicit none
    integer,intent(IN) :: jktrj
    logical,intent(IN) :: disp_switch
    complex(8),allocatable :: zwork0(:,:,:),zwork1(:,:,:)
    integer :: ML,ML1,ML2,ML3,MSP_0,MSP_1,MBZ_0,MBZ_1,MB_0,MB_1,ML_0,ML_1
    integer :: s,k,n,i,i1,i2,i3,ierr,LG_sc(3),iktrj
    real(8) :: vtmp(3),utmp(3),pi2,sum0

    if ( .not.iswitch_banduf ) return

    call write_border( 0, " band_unfold(start)" )

    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    ML_0  = id_grid(myrank_g) + 1
    ML_1  = id_grid(myrank_g) + ir_grid(myrank_g)
    MB_0  = id_band(myrank_b) + 1
    MB_1  = id_band(myrank_b) + ir_band(myrank_b)
    MBZ_0 = id_bzsm(myrank_k) + 1
    MBZ_1 = id_bzsm(myrank_k) + ir_bzsm(myrank_k)
    MSP_0 = id_spin(myrank_s) + 1
    MSP_1 = id_spin(myrank_s) + ir_spin(myrank_s)

    call init_fft

! ---

    pi2 = 2.0d0*acos(-1.0d0)

    allocate( zwork0(0:ML1-1,0:ML2-1,0:ML3-1) ) ; zwork0=(0.0d0,0.0d0)
    allocate( zwork1(0:ML1-1,0:ML2-1,0:ML3-1) ) ; zwork1=(0.0d0,0.0d0)

    do iktrj=1,nktrj_0

       if ( map_ktrj(iktrj,1) /= jktrj ) cycle

       vtmp(:) = kxyz_0(:,iktrj) - kxyz_pc(:,iktrj)
       utmp(:) = matmul( vtmp, aa(:,:) )/pi2
       vtmp(:) = nint( utmp(:) )

       if ( any( abs(vtmp-utmp) > 1.d-13 ) ) stop "stop@band_unfold"

       LG_sc(:) = nint( utmp(:) )

! ---

       do s=MSP_0,MSP_1
       do k=MBZ_0,MBZ_0
       do n=MB_0 ,MB_1

          zwork0(:,:,:)=zero
          i=ML_0-1
          do i3=Igrid(1,3),Igrid(2,3)
          do i2=Igrid(1,2),Igrid(2,2)
          do i1=Igrid(1,1),Igrid(2,1)
             i=i+1
             zwork0(i1,i2,i3) = unk(i,n,k,s)
          end do
          end do
          end do
          if ( map_ktrj(iktrj,2) == 1 ) zwork0=conjg(zwork0)

          call mpi_allreduce( MPI_IN_PLACE, zwork0, ML, MPI_COMPLEX16 &
               ,MPI_SUM, comm_grid, ierr )

          call forward_fft( zwork0, zwork1 )

          sum0=0.0d0
          do i=1,mg_pc
             i1 = mod( LG_pc(1,i)-LG_sc(1)+ML1, ML1 )
             i2 = mod( LG_pc(2,i)-LG_sc(2)+ML2, ML2 )
             i3 = mod( LG_pc(3,i)-LG_sc(3)+ML3, ML3 )
             sum0=sum0+abs( zwork0(i1,i2,i3) )**2
          end do
          weight_uf(iktrj,n,s)=sum0

       end do ! n
       end do ! k
       end do ! s

! --- write unfolding data ---

       if ( myrank == 0 ) then
          write(unit_uf,'(1x,"iktrj=",i6,2x,2(3f10.6,2x))') &
               iktrj,kxyz_pc(1:3,iktrj),kxyz_0(1:3,iktrj)
          do s=MSP_0,MSP_1
          do k=MBZ_0,MBZ_0
          do n=MB_0 ,MB_1
             write(unit_uf,'(1x,2i6,f20.15,2x,g20.10)') &
                  n,s,esp(n,k,s),weight_uf(iktrj,n,s)
          end do
          end do
          end do
       end if

! ---

    end do ! iktrj

    deallocate( zwork1 )
    deallocate( zwork0 )

    call finalize_fft

    call write_border( 0, " band_unfold(end)" )

  END SUBROUTINE band_unfold


  SUBROUTINE finalize_band_unfold
    if ( iswitch_banduf .and. myrank == 0 ) close(unit_uf)
  END SUBROUTINE finalize_band_unfold


END MODULE band_unfold_module
