MODULE io_tools_module

  implicit none

  PRIVATE
  PUBLIC :: init_io_tools
  PUBLIC :: IOTools_readStringKeyword
  PUBLIC :: IOTools_readIntegerKeyword
  PUBLIC :: IOTools_readIntegerKeywords
  PUBLIC :: IOTools_readReal8Keyword
  PUBLIC :: IOTools_readReal8Keywords
  PUBLIC :: IOTools_readIntegerString
  PUBLIC :: IOTools_findKeyword
!  PUBLIC :: IOTools_readRealVectorKeyword
!  PUBLIC :: IOTools_readLogicalKeyword
!  PUBLIC :: IOTools_readIntegerVectorKeyword
!  PUBLIC :: IOTools_bcastIntegerParameter

  integer,parameter :: max_trial_read = 10000

  include 'mpif.h'

  integer :: myrank
  integer :: unit_default
  integer :: unit_input_result = 110

CONTAINS


  SUBROUTINE init_io_tools( myrank_in, unit_in )
    implicit none
    integer,intent(IN) :: myrank_in, unit_in
    myrank = myrank_in
    unit_default = unit_in
!    if ( myrank == 0 ) open(unit_input_result,file="input_result")
  END SUBROUTINE init_io_tools


  SUBROUTINE IOTools_readStringKeyword( keyword, variable, unit_in )
    implicit none
    character(*),intent(IN) :: keyword
    character(*),intent(INOUT) :: variable
    integer,optional,intent(IN) :: unit_in
    character(10) :: cbuf,ckey
    integer :: i,unit
    unit=unit_default ; if ( present(unit_in) ) unit=unit_in
    if ( myrank == 0 ) then
       rewind unit
       do i=1,max_trial_read
          read(unit,*,END=999) cbuf
          call convertToCapital(cbuf,ckey)
          if ( ckey == keyword ) then
             backspace(unit)
             read(unit,*) cbuf, variable
             write(*,'(1x,A10," : ",A10)') keyword, variable
             exit
          end if
       end do ! i
    end if
999 call mpi_bcast( variable,len(variable),MPI_CHARACTER,0,MPI_COMM_WORLD,i )
  END SUBROUTINE IOTools_readStringKeyword


  SUBROUTINE IOTools_readIntegerKeyword( keyword, variable, unit_in )
    implicit none
    character(*),intent(IN) :: keyword
    integer,intent(INOUT) :: variable
    integer,optional,intent(IN) :: unit_in
    character(10) :: cbuf,ckey
    integer :: i,unit
    unit=unit_default ; if ( present(unit_in) ) unit=unit_in
    if ( myrank == 0 ) then
       rewind unit
       do i=1,max_trial_read
          read(unit,*,END=999) cbuf
          call convertToCapital(cbuf,ckey)
          if ( ckey == keyword ) then
             backspace(unit)
             read(unit,*) cbuf,variable
             write(*,'(1x,A10," : ",3I10)') keyword,variable
             exit
          end if
       end do ! i
999    continue
       !write(unit_input_result,'(a10,3x,i10)') keyword,variable
    end if
    call MPI_BCAST(variable,1,MPI_INTEGER,0,MPI_COMM_WORLD,i)
  END SUBROUTINE IOTools_readIntegerKeyword


  SUBROUTINE IOTools_readIntegerKeywords( keyword, variables, unit_in )
    implicit none
    character(*),intent(IN) :: keyword
    integer,intent(INOUT) :: variables(:)
    integer,optional,intent(IN) :: unit_in
    character(10) :: cbuf,ckey
    integer :: i,unit
    unit=unit_default ; if ( present(unit_in) ) unit=unit_in
    if ( myrank == 0 ) then
       rewind unit
       do i=1,max_trial_read
          read(unit,*,END=999) cbuf
          call convertToCapital(cbuf,ckey)
          if ( ckey == keyword ) then
             backspace(unit)
             read(unit,*) cbuf,variables(:)
             write(*,'(1x,A10," : ",3I10)') keyword,variables(:)
             exit
          end if
       end do ! i
    end if
999 call MPI_BCAST(variables,size(variables),MPI_INTEGER,0,MPI_COMM_WORLD,i)
  END SUBROUTINE IOTools_readIntegerKeywords


  SUBROUTINE IOTools_readReal8Keyword( keyword, variable, unit_in )
    implicit none
    character(*),intent(IN) :: keyword
    real(8),intent(INOUT) :: variable
    integer,optional,intent(IN) :: unit_in
    character(10) :: cbuf,ckey
    integer :: i,unit
    unit=unit_default ; if ( present(unit_in) ) unit=unit_in
    if ( myrank == 0 ) then
       rewind unit
       do i=1,max_trial_read
          read(unit,*,END=999) cbuf
          call convertToCapital(cbuf,ckey)
          if ( ckey == keyword ) then
             backspace(unit)
             read(unit,*) cbuf,variable
             if ( abs(variable) < 1.d-2 ) then
                write(*,'(1x,A10," : ",ES15.7)') keyword,variable
             else
                write(*,'(1x,A10," : ",F15.10)') keyword,variable
             end if
             exit
          end if
       end do ! i
    end if
999 call MPI_BCAST(variable,1,MPI_REAL8,0,MPI_COMM_WORLD,i)
  END SUBROUTINE IOTools_readReal8Keyword


  SUBROUTINE IOTools_readReal8Keywords( keyword, variables, unit_in )
    implicit none
    character(*),intent(IN) :: keyword
    real(8),intent(INOUT) :: variables(:)
    integer,optional,intent(IN) :: unit_in
    character(10) :: cbuf,ckey
    integer :: i,unit
    unit=unit_default ; if ( present(unit_in) ) unit=unit_in
    if ( myrank == 0 ) then
       rewind unit
       do i=1,max_trial_read
          read(unit,*,END=999) cbuf
          call convertToCapital(cbuf,ckey)
          if ( ckey == keyword ) then
             backspace(unit)
             read(unit,*) cbuf,variables(:)
             if ( any(abs(variables) < 1.d-2) ) then
                write(*,'(1x,A10," : ",3ES15.7)') keyword,variables(:)
             else
                write(*,'(1x,A10," : ",3F15.10)') keyword,variables(:)
             end if
             exit
          end if
       end do ! i
    end if
999 call MPI_BCAST(variables,size(variables),MPI_REAL8,0,MPI_COMM_WORLD,i)
  END SUBROUTINE IOTools_readReal8Keywords


  SUBROUTINE IOTools_readIntegerString( keyword, variable1, variable2, u, norewind )
    implicit none
    character(*),intent(IN) :: keyword
    integer,intent(INOUT) :: variable1
    character(*),intent(INOUT) :: variable2
    integer,optional,intent(IN) :: u
    logical,optional,intent(IN) :: norewind
    character(10) :: cbuf,ckey
    integer :: i,unit
    unit=unit_default ; if ( present(u) ) unit=u
    if ( myrank == 0 ) then
       if ( .not.present(norewind) ) rewind unit
       do i=1,max_trial_read
          read(unit,*,END=999) cbuf
          call convertToCapital(cbuf,ckey)
          if ( ckey == keyword ) then
             backspace(unit)
             read(unit,*) cbuf,variable1, variable2
             write(*,'(1x,A10," : ",i4,2x,a20)') keyword,variable1,variable2
             exit
          end if
       end do ! i
    end if
999 call MPI_BCAST(variable1,1,MPI_INTEGER,0,MPI_COMM_WORLD,i)
    call MPI_BCAST(variable2,len(variable2),MPI_CHARACTER,0,MPI_COMM_WORLD,i)
  END SUBROUTINE IOTools_readIntegerString


  SUBROUTINE IOTools_findKeyword( keyword, hasKeyword, unit_out, flag_bcast )
    implicit none
    character(*),intent(IN) :: keyword
    logical,intent(OUT) :: hasKeyword
    integer,optional,intent(OUT) :: unit_out
    logical,optional,intent(IN) :: flag_bcast
    integer :: i,unit
    character(10) :: cbuf,ckey
    unit=unit_default
    if ( present(unit_out) ) unit_out=unit
    hasKeyword = .false.
    if ( myrank == 0 ) then
       rewind unit
       do i=1,max_trial_read
          read(unit,*,END=999) cbuf
          call convertToCapital(cbuf,ckey)
          if ( ckey == keyword ) then
             hasKeyword = .true.
             exit
          end if
       end do ! i
       if ( hasKeyword ) write(*,'(1x,A10)') keyword
    end if
999 continue
    if ( present(flag_bcast) ) then
       call MPI_BCAST( hasKeyword, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, i )
    end if
  END SUBROUTINE IOTools_findKeyword

#ifdef TEST
  SUBROUTINE IOTools_readRealKeyword(keyword,unit_number,keyword_variable)
    implicit none
    character(*),intent(IN) :: keyword
    integer,intent(IN) :: unit_number
    real(8),intent(OUT) :: keyword_variable
    logical :: hasKeyword=.false.
    integer :: i
    character(10) :: cbuf,ckey
    integer :: keyword_length
    keyword_length=len_trim(keyword)
    rewind unit_number
    do i=1,10000
      read(unit_number,*,END=999) cbuf
      call convertToCapital(cbuf,ckey)
      if ( ckey==keyword ) then
        backspace(unit_number)
        read(unit_number,*) cbuf,keyword_variable
        hasKeyword=.true.
        exit
      endif
    enddo
999 continue
    if ( hasKeyword ) write(*,'(1x,A10," : ",F20.12)') keyword,keyword_variable
  END SUBROUTINE IOTools_readRealKeyword

  SUBROUTINE IOTools_readRealVectorKeyword(keyword,unit_number,keyword_variable)
    implicit none
    character(*),intent(IN) :: keyword
    integer,intent(IN) :: unit_number
    real(8),intent(OUT) :: keyword_variable(1:3)
    logical :: hasKeyword=.false.
    integer :: i
    character(10) :: cbuf,ckey
    integer :: keyword_length
    keyword_length=len_trim(keyword)
    rewind unit_number
    do i=1,10000
      read(unit_number,*,END=999) cbuf
      call convertToCapital(cbuf,ckey)
      if ( ckey==keyword ) then
        backspace(unit_number)
        read(unit_number,*) cbuf,keyword_variable
        hasKeyword=.true.
        exit
      endif
    enddo
999 continue
    if ( hasKeyword ) write(*,'(1x,A10," : ",3F20.12)') keyword,keyword_variable
  END SUBROUTINE IOTools_readRealVectorKeyword

  SUBROUTINE IOTools_readIntegerVectorKeyword(keyword,unit_number,keyword_variable)
    implicit none
    character(*),intent(IN) :: keyword
    integer,intent(IN) :: unit_number
    integer,intent(OUT) :: keyword_variable(1:3)
    logical :: hasKeyword=.false.
    integer :: i
    character(10) :: cbuf,ckey
    integer :: keyword_length
    keyword_length=len_trim(keyword)
    rewind unit_number
    do i=1,10000
      read(unit_number,*,END=999) cbuf
      call convertToCapital(cbuf,ckey)
      if ( ckey==keyword ) then
        backspace(unit_number)
        read(unit_number,*) cbuf,keyword_variable
        hasKeyword=.true.
        exit
      endif
    enddo
999 continue
    if ( hasKeyword ) write(*,'(1x,A10," : ",3I10)') keyword,keyword_variable
  END SUBROUTINE IOTools_readIntegerVectorKeyword

  SUBROUTINE IOTools_readLogicalKeyword(keyword,unit_number,keyword_variable)
    implicit none
    character(*),intent(IN) :: keyword
    integer,intent(IN) :: unit_number
    logical,intent(OUT) :: keyword_variable
    logical :: hasKeyword=.false.
    integer :: i
    character(10) :: cbuf,ckey
    integer :: keyword_length
    keyword_length=len_trim(keyword)
    rewind unit_number
    do i=1,10000
      read(unit_number,*,END=999) cbuf
      call convertToCapital(cbuf,ckey)
      if ( ckey==keyword ) then
        backspace(unit_number)
        read(unit_number,*) cbuf,keyword_variable
        hasKeyword=.true.
        exit
      endif
    enddo
999 continue
    if ( hasKeyword ) write(*,'(1x,A10," : ",L10)') keyword,keyword_variable
  END SUBROUTINE IOTools_readLogicalKeyword
#endif


  SUBROUTINE IOTools_bcastIntegerParameter( i )
    implicit none
    integer,intent(INOUT) :: i(:)
    integer :: ierr
    call MPI_BCAST( i, size(i), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
  END SUBROUTINE IOTools_bcastIntegerParameter


  SUBROUTINE convertToCapital(cbuf,CKEY)
    implicit none
    character(*),intent(IN)  :: cbuf
    character(*),intent(OUT) :: CKEY
    integer :: j,k,n
    n=len_trim(cbuf)
    CKEY=""
    do j=1,n
      k=iachar( cbuf(j:j) )
      if ( k >= 97 ) k=k-32
      CKEY(j:j) = achar(k)
    end do
  END SUBROUTINE convertToCapital


END MODULE io_tools_module
