MODULE io_write_module

  use parallel_module
  use rgrid_mol_module, only: LL
  use rgrid_module, only: Igrid, Ngrid
  use wf_module
  use aa_module
  use bb_module
  use bz_module, only: kbb, MMBZ

  implicit none

  PRIVATE
  PUBLIC :: simple_wf_io_write

CONTAINS


  SUBROUTINE simple_wf_io_write &
       ( file_wf, IO_ctrl, OC, SYStype, MBwr1, MBwr2, disp_switch )

    implicit none
    character(*),intent(IN) :: file_wf
    integer,intent(IN) :: IO_ctrl, OC, SYStype, MBwr1,MBwr2
    logical,intent(IN) :: disp_switch
    integer,parameter :: unit = 1
    integer :: i,i1,i2,i3,j1,j2,j3,k,n,s,n1,n2
    integer :: istatus(MPI_STATUS_SIZE,123),irank,ierr
    integer,allocatable :: irc(:),ids(:)
    integer,allocatable :: LL2(:,:)
    logical :: flag_related
    integer :: ML,ML1,ML2,ML3,ML0
    integer :: MB,MB_0,MB_1
    integer :: MK,MK_0,MK_1,MS,MS_0,MS_1
    character(len=5) :: cmyrank
    character(len=32) :: file_wf_split
    logical :: type_d, type_z
#ifdef _DRSDFT_
    real(8),allocatable :: utmp(:)
    real(4),allocatable :: utmpSP(:)
    real(8),parameter :: zero=0.0d0
    integer,parameter :: TYPE_MAIN=MPI_REAL8
    integer,parameter :: TYPE_WF=1
#else
    complex(8),allocatable :: utmp(:)
    complex(4),allocatable :: utmpSP(:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    integer,parameter :: TYPE_MAIN=MPI_COMPLEX16
    integer,parameter :: TYPE_WF=0
#endif

    call write_border( 1, " simple_wf_io_write(start)" )
    if ( DISP_SWITCH ) write(*,*) "OC =",OC

    n1  = idisp(myrank)+1
    n2  = idisp(myrank)+ircnt(myrank)
    ML0 = ircnt(myrank)

    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    MB   = MB_WF
    MB_0 = MB_0_WF
    MB_1 = MB_1_WF
    MK   = MK_WF
    MK_0 = MK_0_WF
    MK_1 = MK_1_WF
    MS   = MS_WF
    MS_0 = MS_0_WF
    MS_1 = MS_1_WF

! ---

    allocate( LL2(3,ML) ) ; LL2=0

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

    allocate( irc(0:np_grid-1),ids(0:np_grid-1) )
    irc(0:np_grid-1)=3*ir_grid(0:np_grid-1)
    ids(0:np_grid-1)=3*id_grid(0:np_grid-1)

    call mpi_allgatherv(LL2(1,n1),irc(myrank_g),mpi_integer, &
         LL2,irc,ids,mpi_integer,comm_grid,ierr)

    deallocate( ids, irc )

! ---

    allocate( utmp(ML) ) ; utmp=zero

    select case(OC)
    case default

       if ( DISP_SWITCH ) then
          if ( TYPE_WF == 0 ) print*,'WF: complex(8) -> complex(8)'
          if ( TYPE_WF == 1 ) print*,'WF: real(8) -> real(8)'
       end if

    case(4,5,14,15)

       if ( DISP_SWITCH ) then
          if ( TYPE_WF == 0 ) print*,'WF: complex(8) -> complex(4)'
          if ( TYPE_WF == 1 ) print*,'WF: real(8) -> real(4)'
       end if
       allocate( utmpSP(ML) ) ; utmpSP=zero

    end select

    if ( myrank == 0 ) then

       open( unit, file=file_wf, form="unformatted" )

       if ( OC < 10 ) then

          write(unit) -1
          write(unit) ML,ML1,ML2,ML3
          write(unit) MB,MBwr1,MBwr2
          write(unit) MK, MS, MMBZ
          write(unit) IO_ctrl, OC, TYPE_WF
          write(unit) LL2(:,:)
          write(unit) occ(:,:,:)
          write(unit) aa,bb,kbb

       else if ( OC >= 10 ) then   !--- old format ---

          write(unit) ML,ML1,ML2,ML3
          write(unit) MB,MBwr1,MBwr2
          write(unit) LL2(:,:)
          write(unit) occ(:,:,:)

       end if

    end if

    if ( IO_ctrl == 3 ) then

       if ( myrank == 0 ) close(unit)

       write(cmyrank,'(i5.5)') myrank
       file_wf_split = trim(file_wf)//"."//trim(adjustl(cmyrank))
       open( unit, file=file_wf_split, form="unformatted" )

    end if

    do s=1,MS
    do k=1,MK
    do n=MBwr1,MBwr2

       do irank=0,nprocs-1
          i1=id_band(id_class(irank,4))+1
          j1=id_band(id_class(irank,4))+ir_band(id_class(irank,4))
          i2=id_bzsm(id_class(irank,5))+1
          j2=id_bzsm(id_class(irank,5))+ir_bzsm(id_class(irank,5))
          i3=id_spin(id_class(irank,6))+1
          j3=id_spin(id_class(irank,6))+ir_spin(id_class(irank,6))
          if ( id_grid(id_class(irank,0))==0 .and. i1<=n .and. &
               n<=j1 .and. i2<=k .and. k<=j2 .and. i3<=s .and. &
               s<=j3 ) exit
       end do
       if ( irank >= nprocs ) then
          write(*,*) "ERROR(write_data)",myrank
          stop
       end if

       flag_related = ( MK_0 <= k .and. k <= MK_1 .and. &
                        MB_0 <= n .and. n <= MB_1 .and. &
                        MS_0 <= s .and. s <= MS_1 )

       if ( IO_ctrl == 0 ) then

          if ( flag_related ) then
             call mpi_gatherv(unk(n1,n,k,s),ML0,TYPE_MAIN, &
                  utmp,ir_grid,id_grid,TYPE_MAIN,0,comm_grid,ierr)
          end if
          call mpi_barrier(mpi_comm_world,ierr)
          if ( irank /= 0 .and. myrank_f == 0 ) then
             if ( irank == myrank ) then
                call mpi_send(utmp,ML,TYPE_MAIN,0,0,mpi_comm_world,ierr)
             end if
             if ( myrank == 0 ) then
                call mpi_recv(utmp,ML,TYPE_MAIN,irank,0, &
                     mpi_comm_world,istatus,ierr)
             end if
          end if

          if ( myrank == 0 ) then
             select case(OC)
             case default
                write(unit) utmp(:)
             case(4,5,14,15)
                utmpSP(:)=utmp(:)
                write(unit) utmpSP(:)
             end select
          end if

       else if ( IO_ctrl == 3 ) then

          if ( flag_related ) then
             select case(OC)
             case default
                write(unit) unk(n1:n2,n,k,s)
             case(4,5,14,15)
                utmpSP(n1:n2)=unk(n1:n2,n,k,s)
                write(unit) utmpSP(n1:n2)
             end select
          end if

       end if

    end do ! n
    end do ! k
    end do ! s

    if ( IO_ctrl == 0 .and. myrank == 0 .or. IO_ctrl == 3 ) then
       close(1)
    end if

    if ( allocated(utmpSP) ) deallocate( utmpSP )
    deallocate( utmp )

    deallocate( LL2 )

    if ( DISP_SWITCH ) then
       write(*,*) "write to ",file_wf
       write(*,*) "MBwr1,MBwr2 =",MBwr1,MBwr2
    end if

    call write_border( 1, " simple_wf_io_write(end)" )

  END SUBROUTINE simple_wf_io_write


END MODULE io_write_module
