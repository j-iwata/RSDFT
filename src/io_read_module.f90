MODULE io_read_module

  use parallel_module
  use wf_module
  use rgrid_module
  use rgrid_mol_module, only: LL

  implicit none

  PRIVATE
  PUBLIC :: simple_wf_io_read

CONTAINS

  SUBROUTINE simple_wf_io_read( file_wf2, SYStype, IO_ctrl, DISP_SWITCH )

    implicit none

    character(*),intent(IN) :: file_wf2
    integer,intent(IN) :: SYStype, IO_ctrl
    logical,intent(IN) :: DISP_SWITCH

    integer :: ierr,istatus(MPI_STATUS_SIZE,123),irank
    integer :: itmp(13),n,k,s,i,i1,i2,i3,j1,j2,j3,mx,my,mz
    integer :: ML_tmp,ML1_tmp,ML2_tmp,ML3_tmp,MB_tmp,MB1_tmp,MB2_tmp
    integer :: ML,ML1,ML2,ML3,MB,n1,n2,ML0,MBZ_tmp,MSP_tmp,MMBZ_tmp
    integer :: MB_0,MB_1,MBZ,MBZ_0,MBZ_1,MSP,MSP_0,MSP_1
    integer :: IO_ctrl0, OC, TYPE_WF
    integer,allocatable :: LL2(:,:),LL_tmp(:,:),ir(:),id(:)
    real(8) :: aa0(3,3),bb0(3,3)
    real(8),allocatable :: kbb0(:,:)
    logical :: flag_related, flag_newformat
    character(5) :: cmyrank
    character(32) :: file_wf_split
#ifdef _DRSDFT_
    integer,parameter :: TYPE_MAIN=MPI_REAL8
    integer,parameter :: type_wf_0=1
    real(8),parameter :: zero=0.0d0
    real(8),allocatable :: utmp(:), utmp3(:,:,:)
    real(4),allocatable :: utmpSP(:)
#else
    integer,parameter :: TYPE_MAIN=MPI_COMPLEX16
    integer,parameter :: type_wf_0=0
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8),allocatable :: utmp(:), utmp3(:,:,:)
    complex(4),allocatable :: utmpSP(:)
#endif
    real(8),allocatable :: dtmp(:)
    real(4),allocatable :: dtmpSP(:)

    call write_border( 0, " simple_wf_io_read(start)" )

! ---

    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)
    MB  = MB_WF
    n1  = ML_0_WF
    n2  = ML_1_WF
    ML0 = n2 - n1 + 1

    MB_0  = MB_0_WF
    MB_1  = MB_1_WF
    MBZ   = MK_WF
    MBZ_0 = MK_0_WF
    MBZ_1 = MK_1_WF
    MSP   = MS_WF
    MSP_0 = MS_0_WF
    MSP_1 = MS_1_WF

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
    MB_tmp  = itmp(5)
    MB1_tmp = itmp(6)
    MB2_tmp = itmp(7)
    MBZ_tmp = itmp(8)
    MSP_tmp = itmp(9)
    MMBZ_tmp= itmp(10)
    IO_ctrl0= itmp(11)
    OC      = itmp(12)
    TYPE_WF = itmp(13)

    if ( DISP_SWITCH ) then
       write(*,*) "ML ,ML_tmp                =",ML ,ML_tmp
       write(*,*) "ML1    ,ML2    ,ML3       =",ML1,ML2,ML3
       write(*,*) "ML1_tmp,ML2_tmp,ML3_tmp   =",ML1_tmp,ML2_tmp,ML3_tmp
       write(*,*) "MB ,MB_tmp                =",MB ,MB_tmp
       write(*,*) "MB1_tmp, MB2_tmp          =",MB1_tmp,MB2_tmp
       if ( flag_newformat ) then
          write(*,*) "MBZ_tmp, MSP_tmp,MMBZ_tmp =",MBZ_tmp,MSP_tmp,MMBZ_tmp
          write(*,*) "IO_ctrl0, OC, TYPE_WF     =",IO_ctrl0, OC, TYPE_WF
          write(*,*) "(new format)"
       else
          write(*,*) "(old format)"
       end if
    end if
    if ( ML_tmp /= ML .or. ML1_tmp /= ML1 .or. &
         ML2_tmp /= ML2 .or. ML3_tmp /= ML3 ) stop "stop@simple_wf_io_read"

    if ( MB1_tmp /= 1 .or. MB2_tmp < MB ) then
       if ( DISP_SWITCH ) then
          write(*,*) "******** WARNING! ********"
          write(*,*) "MB1_tmp,MB2_tmp=",MB1_tmp,MB2_tmp
          write(*,*) "stop!"
       end if
       stop
    end if

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

    if ( myrank == 0 ) then
       read(3) occ(:,:,:)
    end if
    call mpi_bcast(occ,MB*MBZ*MSP,mpi_real8,0,mpi_comm_world,ierr)

    if ( myrank == 0 ) then
       write(*,*) "sum(occ)=",sum(occ)
    end if

! ---

    if ( flag_newformat ) then

       if ( myrank == 0 ) read(3)

    end if

! ---

    if ( SYStype == 0 ) then

       allocate( utmp3(0:ML1-1,0:ML2-1,0:ML3-1) ) ; utmp3=zero

    else if ( SYStype == 1 ) then

       mx=(Ngrid(1)-1)/2
       my=(Ngrid(2)-1)/2
       mz=(Ngrid(3)-1)/2
       allocate( utmp3(-mx:mx,-my:my,-mz:mz) ) ; utmp3=0.d0

    end if

    if ( IO_ctrl == 3 ) then

       if ( myrank == 0 ) close(3)

       write(cmyrank,'(i5.5)') myrank
       file_wf_split = trim(file_wf2)//"."//trim(adjustl(cmyrank))

       open(3,file=file_wf_split,form="unformatted")

    end if

    allocate( utmp(ML) ) ; utmp=zero
    if ( OC == 14 .or. OC == 15 ) then
       allocate( utmpSP(ML) ) ; utmpSP=zero
    end if

    if ( type_wf == 1 ) then
       allocate( dtmp(ML) ) ; dtmp=0.0d0
       if ( OC == 14 .or. OC == 15 ) then
          allocate( dtmpSP(ML) ) ; dtmpSP=0.0d0
       end if
    end if

    do s=1,MSP
    do k=1,MBZ
    do n=MB1_tmp,min(MB2_tmp,MB)

       do irank=0,nprocs-1
          i1=id_band(id_class(irank,4))+1
          j1=id_band(id_class(irank,4))+ir_band(id_class(irank,4))
          i2=id_bzsm(id_class(irank,5))+1
          j2=id_bzsm(id_class(irank,5))+ir_bzsm(id_class(irank,5))
          i3=id_spin(id_class(irank,6))+1
          j3=id_spin(id_class(irank,6))+ir_spin(id_class(irank,6))
          if ( id_grid(id_class(irank,0))==0 .and. &
               i1<=n .and. n<=j1 .and. &
               i2<=k .and. k<=j2 .and. &
               i3<=s .and. s<=j3 ) exit
       end do
       if ( irank >= nprocs ) then
          write(*,*) "ERROR!(stpo@simple_wf_io_read): myrank=",myrank
          stop
       end if
            
       flag_related = .false.
       if ( MBZ_0<=k .and. k<=MBZ_1 .and. MB_0 <=n .and. n<=MB_1 &
            .and. MSP_0<=s .and. s<=MSP_1 ) then
          flag_related = .true.
       end if

       if ( IO_ctrl == 0 ) then

          if ( myrank == 0 ) then

             select case( OC )
             case default
                if ( type_wf == 1 ) then
                   read(3) dtmp
                   utmp=dtmp
                else
                   read(3) utmp
                end if
             case(14,15)
                if ( type_wf == 1 ) then
                   read(3) dtmpSP
                   utmp=dtmpSP
                else
                   read(3) utmpSP
                   utmp=utmpSP
                end if
             end select

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

          if ( irank /= 0 .and. myrank_f == 0 ) then

             if ( myrank == 0 ) then
                call mpi_send(utmp,ML,TYPE_MAIN,irank,0, &
                     mpi_comm_world,ierr)
             end if
             if ( myrank == irank ) then
                call mpi_recv(utmp,ML,TYPE_MAIN,0,0, &
                     mpi_comm_world,istatus,ierr)
             end if

          end if

          call mpi_barrier(mpi_comm_world,ierr)

          if ( flag_related ) then
             call mpi_scatterv(utmp,ir_grid,id_grid,TYPE_MAIN &
                  ,unk(n1,n,k,s),ML0,TYPE_MAIN,0,comm_grid,ierr)
          end if

       else if ( IO_ctrl == 3 ) then

          if ( flag_related ) then
             select case(OC)
             case default
                if ( type_wf == 1 ) then
                   read(3) dtmp(n1:n2)
                   if ( flag_related ) unk(n1:n2,n,k,s)=dtmp(n1:n2)
                else
                   read(3) utmp(n1:n2)
                   if ( flag_related ) unk(n1:n2,n,k,s)=utmp(n1:n2)
                end if
             case(14,15)
                if ( type_wf == 1 ) then
                   read(3) dtmpSP(n1:n2)
                   if ( flag_related ) unk(n1:n2,n,k,s)=dtmpSP(n1:n2)
                else
                   read(3) utmpSP(n1:n2)
                   if ( flag_related ) unk(n1:n2,n,k,s)=utmpSP(n1:n2)
                end if
             end select
          end if

       end if

    end do ! n
    end do ! k
    end do ! s

    if ( IO_ctrl == 0 .and. myrank == 0 .or. IO_ctrl == 3 ) close(3)

    if ( allocated(dtmpSP) ) deallocate( dtmpSP )
    if ( allocated(dtmp) ) deallocate( dtmp )
    if ( allocated(utmpSP) ) deallocate( utmpSP )
    deallocate( utmp )
    deallocate( utmp3 )
    deallocate( LL_tmp )
    deallocate( LL2 )

#ifdef _DRSDFT_
    call mpi_bcast( unk,size(unk),MPI_REAL8,0,comm_fkmb,ierr )
#else
    call mpi_bcast( unk,size(unk),MPI_COMPLEX16,0,comm_fkmb,ierr )
#endif

    if ( DISP_SWITCH ) then
       write(*,*) "read from ",file_wf2
    end if

    call write_border( 0, " simple_wf_io_read(start)" )

    return

  END SUBROUTINE simple_wf_io_read

END MODULE io_read_module
