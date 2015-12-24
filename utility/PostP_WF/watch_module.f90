MODULE watch_module

  implicit none

  PRIVATE
  PUBLIC :: watch,watcht,read_watch,global_watch,read_oldformat_watch

  real(8) :: ct0=0.d0, ctt=0.d0
  real(8) :: ett=0.d0
  integer :: count0=0

  real(8) :: etime_limit

  real(8) :: global_ctime0,global_etime0
  logical :: flag_count_start=.true.

CONTAINS

  SUBROUTINE read_watch(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i
    character(7) :: cbuf,ckey
    etime_limit = 1.d100
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:7) == "ETLIMIT" ) then
             backspace(unit)
             read(unit,*) cbuf,etime_limit
          end if
       end do
999    continue
       write(*,*) "etime_limit=",etime_limit
    end if
    call send_watch(0)
  END SUBROUTINE read_watch

  SUBROUTINE read_oldformat_watch(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    if ( rank == 0 ) then
       read(unit,*) etime_limit
       write(*,*) "etime_limit=",etime_limit
    end if
    call send_watch(0)
  END SUBROUTINE read_oldformat_watch

  SUBROUTINE send_watch(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(etime_limit,1,MPI_REAL8,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_watch

  SUBROUTINE watch(ctime,etime)
    real(8),intent(OUT) :: ctime,etime
    include 'mpif.h'
    call cpu_time(ctime)
    etime=mpi_wtime()
  END SUBROUTINE watch

  SUBROUTINE watcht(disp_switch,indx,icnt)
    logical,intent(IN) :: disp_switch
    character(*),intent(IN) :: indx 
    integer,intent(IN) :: icnt
    integer :: count,count_rate
    real(8) :: ct
    call cpu_time(ct)
    call system_clock(count,count_rate)
    if ( icnt == 0 ) then
       ctt=0.d0
       ett=0.d0
    else if ( icnt == 1 ) then
       ctt=ct-ct0
       ett=real(count-count0)/real(count_rate)
       if ( disp_switch ) write(*,*) "timet(",indx,")=",ctt,ett
    end if
    ct0=ct
    count0=count
  END SUBROUTINE watcht

  SUBROUTINE global_watch(flag_timelimit)
    logical,intent(OUT) :: flag_timelimit
    integer :: ierr
    real(8) :: ct,et
    include 'mpif.h'
    flag_timelimit=.false.
    call cpu_time(ct)
    et=mpi_wtime()
    if ( flag_count_start ) then
       global_ctime0=ct
       global_etime0=et
       flag_count_start=.false.
       return
    end if
    ct=et-global_etime0
    call mpi_allreduce(ct,et,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
    if ( et > etime_limit ) flag_timelimit=.true.
  END SUBROUTINE global_watch

END MODULE watch_module
