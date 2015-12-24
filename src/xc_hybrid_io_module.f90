MODULE xc_hybrid_io_module

  use parallel_module
  use wf_module
  use rgrid_module
  use rgrid_mol_module, only: LL

  implicit none

  PRIVATE
  PUBLIC :: read_xc_hybrid_io

CONTAINS

  SUBROUTINE read_xc_hybrid_io( file_wf2, SYStype, IO_ctrl, DISP_SWITCH &
       ,n1,n2, MSP_0,MSP_1, unk_hf, occ_hf, kbb_hf &
       ,FKMB,FKMB_0,FKMB_1, FKBZ,FKBZ_0,FKBZ_1,FKMMBZ )

    implicit none

    character(*),intent(IN) :: file_wf2
    integer,intent(IN) :: SYStype, IO_ctrl
    logical,intent(IN) :: DISP_SWITCH
    integer,intent(IN) :: n1,n2,MSP_0,MSP_1
#ifdef _DRSDFT_
    real(8),allocatable,intent(INOUT) :: unk_hf(:,:,:,:)
#else
    complex(8),allocatable,intent(INOUT) :: unk_hf(:,:,:,:)
#endif
    real(8),allocatable,intent(INOUT) :: occ_hf(:,:,:), kbb_hf(:,:)
    integer,intent(INOUT) :: FKMB,FKMB_0,FKMB_1
    integer,intent(INOUT) :: FKBZ,FKBZ_0,FKBZ_1,FKMMBZ

    integer :: itmp(13),n,k,s,i,i1,i2,i3,j1,j2,j3,mx,my,mz,m,ierr
    integer :: ML_tmp,ML1_tmp,ML2_tmp,ML3_tmp,MB_tmp,MB1_tmp,MB2_tmp
    integer :: ML,ML1,ML2,ML3,MB,ML0,MBZ_tmp,MSP_tmp
    integer :: MB_0,MB_1,MBZ,MBZ_0,MBZ_1,MSP
    integer :: IO_ctrl0, OC, TYPE_WF
    integer,allocatable :: LL2(:,:),LL_tmp(:,:),ir(:),id(:)
    logical :: flag_related, flag_newformat
    character(5) :: cmyrank
    character(32) :: file_wf_split
    real(8) :: aa(3,3),bb(3,3)
#ifdef _DRSDFT_
    integer,parameter :: TYPE_MAIN=MPI_REAL8
    real(8),parameter :: zero=0.0d0
    real(8),allocatable :: utmp(:), utmp3(:,:,:)
#else
    integer,parameter :: TYPE_MAIN=MPI_COMPLEX16
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8),allocatable :: utmp(:), utmp3(:,:,:)
#endif

    if ( DISP_SWITCH ) then
       write(*,'(a40," read_xc_hybrid_io")') repeat("-",40)
    end if

! ---

    ML    = Ngrid(0)
    ML1   = Ngrid(1)
    ML2   = Ngrid(2)
    ML3   = Ngrid(3)
    MB    = sum(ir_band)
    ML0   = n2 - n1 + 1
    MBZ   = sum(ir_bzsm)
    MSP   = sum(ir_spin)
    MB_0  = id_band(myrank_b)+1
    MB_1  = id_band(myrank_b)+ir_band(myrank_b)
    MBZ_0 = id_bzsm(myrank_k)+1
    MBZ_1 = id_bzsm(myrank_k)+ir_bzsm(myrank_k)

! ---

    allocate( LL2(3,ML)    ) ; LL2=0
    allocate( LL_tmp(3,ML) ) ; LL_tmp=0

    if ( SYStype == 0 ) then

       i=n1-1
       do i3=Igrid(1,3),Igrid(2,3)
       do i2=Igrid(1,2),Igrid(2,2)
       do i1=Igrid(1,1),Igrid(2,1)
          i=i+1
          LL2(1,i)=i1
          LL2(2,i)=i2
          LL2(3,i)=i3
       end do
       end do
       end do

    else if ( SYStype == 1 ) then

       LL2(1:3,n1:n2) = LL(1:3,n1:n2)

    end if

    allocate( ir(0:np_grid-1), id(0:np_grid-1) )

    ir(0:np_grid-1)=3*ir_grid(0:np_grid-1)
    id(0:np_grid-1)=3*id_grid(0:np_grid-1)
    call mpi_allgatherv(LL2(1,n1),ir(myrank_g),mpi_integer &
         ,LL2,ir,id,mpi_integer,comm_grid,ierr)

    deallocate( id,ir )

! ---

    if ( myrank == 0 ) then

       open(3,file=file_wf2,form='unformatted')
       read(3) itmp(1)

       flag_newformat = .false.
       if ( itmp(1) < 0 ) flag_newformat = .true.
       write(*,*) "flag_newformat=",flag_newformat

       if ( flag_newformat ) then
          read(3) itmp(1:4)
          read(3) itmp(5:7)
          read(3) itmp(8:10)
          read(3) itmp(11:13)
       else
          rewind 3
          read(3) itmp(1:4)
          read(3) itmp(5:7)
       end if

    end if

    call mpi_bcast(flag_newformat,1,mpi_logical,0,mpi_comm_world,ierr)
    call mpi_bcast(itmp,13,mpi_integer,0,mpi_comm_world,ierr)
    ML_tmp  = itmp(1)
    ML1_tmp = itmp(2)
    ML2_tmp = itmp(3)
    ML3_tmp = itmp(4)
    FKMB    = itmp(5)
    MB1_tmp = itmp(6)
    MB2_tmp = itmp(7)
    FKBZ    = itmp(8)
    MSP_tmp = itmp(9)
    FKMMBZ  = itmp(10)
    IO_ctrl0= itmp(11)
    OC      = itmp(12)
    TYPE_WF = itmp(13)

    if ( .not.flag_newformat ) then
       FKBZ = MBZ
       MSP_tmp = MSP
    end if

    if ( DISP_SWITCH ) then
       write(*,*) "ML ,ML_tmp              =",ML ,ML_tmp
       write(*,*) "ML1    ,ML2    ,ML3     =",ML1,ML2,ML3
       write(*,*) "ML1_tmp,ML2_tmp,ML3_tmp =",ML1_tmp,ML2_tmp,ML3_tmp
       write(*,*) "FKMB ,MB                =",FKMB ,MB
       write(*,*) "MB1_tmp, MB2_tmp        =",MB1_tmp,MB2_tmp
       if ( flag_newformat ) then
          write(*,*) "FKBZ, MSP_tmp, FKMMBZ   =",FKBZ,MSP_tmp,FKMMBZ
          write(*,*) "IO_ctrl0, OC, TYPE_WF   =",IO_ctrl0, OC, TYPE_WF
          write(*,*) "(new format)"
       else
          write(*,*) "FKBZ = MBZ              =",FKBZ
          write(*,*) "MSP_tmp = MSP           =",MSP
          write(*,*) "(old format)"
       end if
    end if
    if ( ML_tmp /= ML .or. ML1_tmp /= ML1 .or. &
         ML2_tmp /= ML2 .or. ML3_tmp /= ML3 ) stop "stop@simple_wf_io_read"

    if ( MB1_tmp /= 1 .or. MB2_tmp /= FKMB ) then
       if ( DISP_SWITCH ) then
          write(*,*) "******** WARNING! ********"
          write(*,*) "MB1_tmp,MB2_tmp=",MB1_tmp,MB2_tmp
          write(*,*) "stop!"
       end if
       stop
    end if

! ---

    if ( myrank == 0 ) then
       read(3) LL_tmp(:,:)
    end if
    call mpi_bcast(LL_tmp,3*ML,mpi_integer,0,mpi_comm_world,ierr)

    if ( DISP_SWITCH ) then
       write(*,*) "minval(LL_tmp),maxval(LL_tmp)"
       write(*,*) minval(LL_tmp(1,1:ML) ),maxval( LL_tmp(1,1:ML))
       write(*,*) minval(LL_tmp(2,1:ML) ),maxval( LL_tmp(2,1:ML))
       write(*,*) minval(LL_tmp(3,1:ML) ),maxval( LL_tmp(3,1:ML))
    end if

! ---

    FKBZ_0 = 1
    FKBZ_1 = FKBZ
    !FKMB_0 = 1
    !FKMB_1 = FKMB

! ---

    allocate( occ_hf(FKMB_0:FKMB_1,FKBZ_0:FKBZ_1,MSP) )
    occ_hf=0.0d0

! ---

    if ( myrank == 0 ) then
       read(3) occ_hf(:,:,:)
    end if
    call mpi_bcast(occ_hf,size(occ_hf),mpi_real8,0,mpi_comm_world,ierr)

    if ( myrank == 0 ) then
       write(*,*) "sum(occ_hf)=",sum(occ_hf)
    end if

! ---

    if ( flag_newformat ) then

       allocate( kbb_hf(3,FKBZ_0:FKBZ_1) ) ; kbb_hf=0.0d0

       if ( myrank == 0 ) read(3) aa,bb,kbb_hf

       call mpi_bcast(kbb_hf,size(kbb_hf),mpi_real8,0,mpi_comm_world,ierr)

    end if

! ---

    allocate( unk_hf(n1:n2,FKMB_0:FKMB_1,FKBZ_0:FKBZ_1,MSP_0:MSP_1) )
    unk_hf=zero

! ---

    if ( IO_ctrl == 0 ) then

       if ( SYStype == 0 ) then
          allocate( utmp3(0:ML1-1,0:ML2-1,0:ML3-1) ) ; utmp3=zero
       else if ( SYStype == 1 ) then
          mx=(Ngrid(1)-1)/2
          my=(Ngrid(2)-1)/2
          mz=(Ngrid(3)-1)/2
          allocate( utmp3(-mx:mx,-my:my,-mz:mz) ) ; utmp3=0.d0
       end if

       allocate( utmp(ML) ) ; utmp=zero

       do s=1,MSP
       do k=1,FKBZ
       do n=1,FKMB

          if ( myrank == 0 ) then

             read(3) utmp

             do i=1,ML
                i1=LL_tmp(1,i) ; i2=LL_tmp(2,i) ; i3=LL_tmp(3,i)
                utmp3(i1,i2,i3)=utmp(i)
             end do
             do i=1,ML
                i1=LL2(1,i) ; i2=LL2(2,i) ; i3=LL2(3,i)
                utmp(i)=utmp3(i1,i2,i3)
             end do

          end if

          call mpi_barrier(mpi_comm_world,ierr)

          call mpi_bcast(utmp,ML,TYPE_MAIN,0,mpi_comm_world,ierr)

          if ( n < FKMB_0 .or. FKMB_1 < n ) cycle
          if ( s < MSP_0 .or. MSP_1 < s ) cycle
          unk_hf(n1:n2,n,k,s) = utmp(n1:n2)

       end do ! n
       end do ! k
       end do ! s

       deallocate( utmp  )
       deallocate( utmp3 )

    else if ( IO_ctrl == 3 ) then

       if ( myrank == 0 ) close(3)

       write(cmyrank,'(i5.5)') myrank
       file_wf_split = trim(file_wf2)//"."//trim(adjustl(cmyrank))

       open(3,file=file_wf_split,form="unformatted")

       do s=MSP_0,MSP_1
       do k=MBZ_0,MBZ_1
       do n=MB_0 ,MB_1

          if ( n < FKMB_0 .or. FKMB_1 < n ) cycle
          read(3) unk_hf(n1:n2,n,k,s)

       end do ! n
       end do ! k
       end do ! s

       m = size( unk_hf,1 )*size( unk_hf,2 )
       do s=MSP_0,MSP_1
       do k=MBZ_0,MBZ_1
          call MPI_ALLREDUCE( MPI_IN_PLACE, unk_hf(:,:,k,s), m &
               , TYPE_MAIN, MPI_SUM, comm_band, ierr )
       end do ! k
       end do ! s

       m = size( unk_hf,1 )*size( unk_hf,2 )*size( unk_hf,3 )
       do s=MSP_0,MSP_1
          call MPI_ALLREDUCE( MPI_IN_PLACE, unk_hf(:,:,:,s), m &
               , TYPE_MAIN, MPI_SUM, comm_bzsm, ierr )
       end do ! s

    end if

    if ( IO_ctrl == 0 .and. myrank == 0 .or. IO_ctrl == 3 ) close(3)

    if ( DISP_SWITCH ) then
       write(*,*) "read from ",file_wf2
    end if

    deallocate( LL_tmp )
    deallocate( LL2 )

    return

  END SUBROUTINE read_xc_hybrid_io

END MODULE xc_hybrid_io_module
