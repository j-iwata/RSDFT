SUBROUTINE stop_program( indx )
  implicit none
  character(*),intent(IN) :: indx
  integer :: ierr
  include 'mpif.h'
  write(*,'(/,a30 " stop_program is called !!! ",a30)') &
       repeat("-",30),repeat("-",30)
  write(*,*) indx
  call MPI_FINALIZE( ierr )
  stop "stop@stop_program"
END SUBROUTINE stop_program


SUBROUTINE stop_program_f( indx )
  implicit none
  character(*),intent(IN) :: indx
  integer :: errcode, ierr
  include 'mpif.h'
  write(*,'(/,a30 " stop_program_f is called !!! ",a30)') &
       repeat("-",30),repeat("-",30)
  write(*,*) indx
  call MPI_ABORT( MPI_COMM_WORLD, errcode, ierr )
  stop "stop@stop_program_f"
END SUBROUTINE stop_program_f
