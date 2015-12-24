MODULE watch_module

  implicit none

  PRIVATE
  PUBLIC :: watch,watchs,watcht,read_watch,global_watch,read_oldformat_watch
  PUBLIC :: watcha, write_watcha
  PUBLIC :: watchb, write_watchb, watchb_omp
  PUBLIC :: time_cgpc, time_hmlt, time_kine, time_nlpp, time_bcfd
  PUBLIC :: time_cgpc_indx, time_hmlt_indx, time_kine_indx, time_nlpp_indx
  PUBLIC :: time
  PUBLIC :: init_time_watch, calc_time_watch

  include 'mpif.h'

  real(8) :: ct0=0.d0, ctt=0.d0
  real(8) :: ett=0.d0
  integer :: count0=0

  real(8) :: time_cgpc(2,16)
  real(8) :: time_hmlt(2,4)
  real(8) :: time_kine(2,16)
  real(8) :: time_nlpp(2,8)
  real(8) :: time_bcfd(2,8)
  character(5) :: time_hmlt_indx(4)
  character(5) :: time_kine_indx(16)
  character(5) :: time_nlpp_indx(8)
  character(5) :: time_cgpc_indx(16)
  data time_hmlt_indx(1:4)/"kine","loc","nlc","exx"/
  data time_kine_indx(1:11)/"kine1","kine2","bc","kine3","kine4" &
                     ,"recv","send","spack","waita","final","totbc"/
  data time_nlpp_indx(1:7)/"nlc1","nlcom","nlc2","com1","com2","com3","com4"/
  data time_cgpc_indx(8:13)/"recv","send","spack","waita","final","totbc"/

  real(8) :: etime_limit

  real(8) :: global_ctime0,global_etime0
  logical :: flag_count_start=.true.

  integer,parameter :: max_timea_counter=64
  real(8) :: timea(0:max_timea_counter,2)

  type time
     real(8) :: t0, t1
     real(8) :: tmin, tmax
  end type time

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
    call mpi_bcast(etime_limit,1,MPI_REAL8,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_watch

  SUBROUTINE watch(ctime,etime)
    real(8),intent(OUT) :: ctime,etime
    call cpu_time(ctime)
    etime=mpi_wtime()
  END SUBROUTINE watch

  SUBROUTINE watchs(ctime,etime,icnt)
    real(8),intent(INOUT) :: ctime,etime
    integer,intent(IN) :: icnt
    integer :: count,count_rate
    real(8) :: ct
    call cpu_time(ct)
    call system_clock(count,count_rate)
    if ( icnt == 1 ) then
       ctime=ctime+ct-ct0
       etime=etime+real(count-count0)/real(count_rate)
    end if
    ct0=ct
    count0=count
  END SUBROUTINE watchs

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

  SUBROUTINE global_watch(disp_switch,flag_timelimit)
    implicit none
    logical,intent(IN) :: disp_switch
    logical,optional,intent(OUT) :: flag_timelimit
    integer :: ierr
    real(8) :: ct,et,s(2),r(2)
    call cpu_time(ct)
    et=mpi_wtime()
    if ( flag_count_start ) then
       global_ctime0=ct
       global_etime0=et
       flag_count_start=.false.
       return
    end if
    s(1)=ct-global_ctime0
    s(2)=et-global_etime0
    call mpi_allreduce(s,r,2,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
    if ( present(flag_timelimit) ) then
       flag_timelimit=.false.
       if ( r(2) > etime_limit ) flag_timelimit=.true.
    end if
    if ( disp_switch ) write(*,'(1x,"TIME(END)",3f12.5)') r(1:2)
  END SUBROUTINE global_watch


  SUBROUTINE watcha( icounter )
    implicit none
    integer,intent(INOUT) :: icounter
    integer :: n
    n = min( max(icounter+1,0), max_timea_counter )
    call cpu_time( timea(n,1) )
    timea(n,2) = mpi_wtime()
    icounter = n
  END SUBROUTINE watcha

  SUBROUTINE write_watcha( n, indx, unit )
    implicit none
    integer,intent(IN) :: n
    character(*),intent(IN) :: indx
    integer,optional,intent(IN) :: unit
    integer :: i,j,u
    u=6
    if ( present(unit) ) u=unit
    do i=1,n
       write(u,'(1x,i3,2x,a,3x,2x,2f15.5)') &
            i, "timea("//indx//")", (timea(i,j)-timea(i-1,j),j=1,2)
    end do
  END SUBROUTINE write_watcha


  SUBROUTINE watchb( t_tmp, t_out )
    implicit none
    real(8),intent(INOUT) :: t_tmp(2)
    real(8),optional,intent(INOUT) :: t_out(2)
    real(8) :: tnow(2)
    call cpu_time( tnow(1) )
    tnow(2) = mpi_wtime()
    if ( present(t_out) ) t_out(:) = t_out(:) + tnow(:) - t_tmp(:)
    t_tmp(:) = tnow(:)
  END SUBROUTINE watchb

  SUBROUTINE write_watchb( t_results, n, indx, unit )
    implicit none
    integer,intent(IN) :: n
    real(8),intent(IN) :: t_results(2,n)
    character(*),intent(IN) :: indx(n)
    integer,optional,intent(IN) :: unit
    integer :: i,j,u
    u=6
    if ( present(unit) ) u=unit
    do i=1,n
       write(u,'(1x,i3,2x,a,3x,2x,2f15.5)') &
            i, "timeb("//indx(i)//")", (t_results(j,i),j=1,2)
    end do
  END SUBROUTINE write_watchb


  SUBROUTINE watchb_omp( t_tmp, t_out )
    implicit none
    real(8),intent(INOUT) :: t_tmp(2)
    real(8),optional,intent(INOUT) :: t_out(2)
    real(8) :: tnow(2)
!$OMP master
    call cpu_time( tnow(1) )
    tnow(2) = mpi_wtime()
    if ( present(t_out) ) t_out(:) = t_out(:) + tnow(:) - t_tmp(:)
    t_tmp(:) = tnow(:)
!$OMP end master
  END SUBROUTINE watchb_omp


  SUBROUTINE init_time_watch( t )
    implicit none
    type(time),intent(INOUT) :: t
    t%t0 = mpi_wtime()
  END SUBROUTINE init_time_watch

  SUBROUTINE calc_time_watch( t, comm_in )
    implicit none
    type(time),intent(INOUT) :: t
    integer,optional,intent(IN) :: comm_in
    integer :: i, comm
    t%t1 = mpi_wtime()
    t%t0 = t%t1 - t%t0
    comm=MPI_COMM_WORLD ; if ( present(comm_in) ) comm=comm_in
    call mpi_allreduce( t%t0, t%tmin, 1, MPI_REAL8, MPI_MIN, comm, i )
    call mpi_allreduce( t%t0, t%tmax, 1, MPI_REAL8, MPI_MAX, comm, i )
  END SUBROUTINE calc_time_watch


END MODULE watch_module
