MODULE info_module

  implicit none

  PRIVATE
  PUBLIC :: open_info, write_info, close_info

  integer,parameter :: unit_info = 100, unit_stdout = 6
  character(30) :: file_name = "RSDFT_INFO"
  integer :: myrank
  character(8)  :: date_start, date_end
  character(10) :: time_start, time_end

CONTAINS

  SUBROUTINE open_info(rank)
    implicit none
    integer,intent(IN) :: rank
    myrank = rank
!    if ( rank == 0 ) then
!       open(unit_info,file=file_name)
!    end if
    if ( myrank == 0 ) call header_info( unit_stdout )
  END SUBROUTINE open_info

  SUBROUTINE write_info(info)
    implicit none
    character(*),intent(IN) :: info
    if ( info == "" ) return
    if ( myrank == 0 ) then
       write(*,*) "INFO: ",info
       write(unit_info,*) "INFO: ",info
    end if
  END SUBROUTINE write_info

  SUBROUTINE close_info
    if ( myrank == 0 ) call footer_info( unit_stdout )
!    if ( myrank == 0 ) close(unit_info)
  END SUBROUTINE close_info

  SUBROUTINE header_info( unit )
    implicit none
    integer,intent(IN) :: unit
!    character(8)  :: date
!    character(10) :: time

    write(unit,'(a70)') repeat("-",70)
    write(unit,*) "RSDFT ver.1.2.2"

    write(unit,'(a60," header_info")') repeat("-",60)

    call date_and_time(DATE=date_start,TIME=time_start)
    call write_date_and_time( unit, "Start time  ", date_start, time_start )
    write(unit,*) "preprocessor option list ( -cpp -Dxxxx )"
#ifdef _DRSDFT_
       write(unit,*) "_DRSDFT_"," (wave functions are REAL8)"
#else
       write(unit,*) "        "," (wave functions are COMPLEX16)"
#endif
#ifdef _SPLINE_
       write(unit,*) "_SPLINE_"," (spline interpolation is use in ps_nloc2)"
#endif
#ifdef _FFTE_
       write(unit,*) "_FFTE_"," (FFTE is used)"
#endif
#ifdef _FFTW_
       write(unit,*) "_FFTW_"," (FFTW is used)"
#endif
#ifdef _LAPACK_
       write(unit,*) "_LAPACK_"," (LAPACK is used instead of ScaLAPACK"
#endif

  END SUBROUTINE header_info


  SUBROUTINE write_date_and_time( unit, indx, date, time )
    implicit none
    integer,intent(IN) :: unit
    character(*),intent(IN) :: indx, date, time
    write(unit,*) indx//date(1:4)//"/"//date(5:6)//"/"//date(7:8) &
                      //" "//time(1:2)//":"//time(3:4)//":"//time(5:10)
  END SUBROUTINE write_date_and_time


  SUBROUTINE footer_info( unit )
    implicit none
    integer,intent(IN) :: unit
    call date_and_time(DATE=date_end,TIME=time_end)
    call write_date_and_time( unit, "Start time  ", date_start, time_start )
    call write_date_and_time( unit, "End   time  ", date_end, time_end )
  END SUBROUTINE footer_info


END MODULE info_module
