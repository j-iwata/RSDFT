FUNCTION exit_program()
  implicit none
  logical :: exit_program
  integer :: ierr,myrank
  include 'mpif.h'
  call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )
  exit_program=.false.
  if ( myrank == 0 ) then
     if ( .not.exit_program ) inquire( FILE="exit", EXIST=exit_program )
     if ( .not.exit_program ) inquire( FILE="EXIT", EXIST=exit_program )
  end if
  call MPI_BCAST( exit_program,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr )
END FUNCTION exit_program
