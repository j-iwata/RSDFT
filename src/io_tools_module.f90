MODULE io_tools_module

  implicit none

  PRIVATE
  PUBLIC :: init_io_tools
  PUBLIC :: IOTools_readStringKeyword
  PUBLIC :: IOTools_readIntegerKeyword
  PUBLIC :: IOTools_readReal8Keyword
  PUBLIC :: IOTools_readReal8Keywords
  PUBLIC :: IOTools_readIntegerString
  PUBLIC :: IOTools_findKeyword

  integer,parameter :: max_trial_read = 10000

#ifndef _NOMPI_
  include 'mpif.h'
#endif

  integer :: myrank = 0
  integer :: unit_default = 1
  integer :: unit_input_result = 110
  logical :: flag_init = .false.

  INTERFACE IOTools_readIntegerKeyword
     MODULE PROCEDURE IOTools_readIntegerKeyword_sca &
                     ,IOTools_readIntegerKeyword_vec
  END INTERFACE

CONTAINS


  SUBROUTINE init_io_tools( myrank_in, unit_in )
    implicit none
    integer,intent(IN) :: myrank_in, unit_in
    myrank = myrank_in
    unit_default = unit_in
!    if ( myrank == 0 ) open(unit_input_result,file="input_result")
    flag_init = .true.
  END SUBROUTINE init_io_tools

  SUBROUTINE check_init
    implicit none
    integer :: ierr
#ifdef _NOMPI_
    flag_init = .true.
    return
#else
    include 'mpif.h'
    if ( .not.flag_init ) then
       call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )
    end if
    flag_init = .true.
#endif
  END SUBROUTINE check_init


  SUBROUTINE IOTools_readStringKeyword( keyword, variable, unit_in )
    implicit none
    character(*),intent(IN) :: keyword
    character(*),intent(INOUT) :: variable
    integer,optional,intent(IN) :: unit_in
    character(10) :: cbuf,ckey
    integer :: i,unit
    call check_init
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
999    continue
    end if
#ifndef _NOMPI_
    call MPI_BCAST( variable,len(variable),MPI_CHARACTER,0,MPI_COMM_WORLD,i )
#endif
  END SUBROUTINE IOTools_readStringKeyword


  SUBROUTINE IOTools_readIntegerKeyword_sca( keyword, variable, unit_in )
    implicit none
    character(*),intent(IN) :: keyword
    integer,intent(INOUT) :: variable
    integer,optional,intent(IN) :: unit_in
    character(10) :: cbuf,ckey
    integer :: i,unit
    call check_init
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
#ifndef _NOMPI_
    call MPI_BCAST(variable,1,MPI_INTEGER,0,MPI_COMM_WORLD,i)
#endif
  END SUBROUTINE IOTools_readIntegerKeyword_sca


  SUBROUTINE IOTools_readIntegerKeyword_vec( keyword, variables, unit_in )
    implicit none
    character(*),intent(IN) :: keyword
    integer,intent(INOUT) :: variables(:)
    integer,optional,intent(IN) :: unit_in
    character(10) :: cbuf,ckey
    integer :: i,unit
    call check_init
    unit=unit_default ; if ( present(unit_in) ) unit=unit_in
    if ( myrank == 0 ) then
       rewind unit
       do i=1,max_trial_read
          read(unit,*,END=999) cbuf
          call convertToCapital(cbuf,ckey)
          if ( ckey == keyword ) then
             backspace(unit)
             read(unit,*) cbuf,variables(:)
             write(*,'(1x,A10," : ",7I10)') keyword,variables(:)
             exit
          end if
       end do ! i
999    continue
    end if
#ifndef _NOMPI_
    call MPI_BCAST(variables,size(variables),MPI_INTEGER,0,MPI_COMM_WORLD,i)
#endif
  END SUBROUTINE IOTools_readIntegerKeyword_vec


  SUBROUTINE IOTools_readReal8Keyword( keyword, variable, unit_in )
    implicit none
    character(*),intent(IN) :: keyword
    real(8),intent(INOUT) :: variable
    integer,optional,intent(IN) :: unit_in
    character(10) :: cbuf,ckey
    integer :: i,unit
    call check_init
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
999    continue
    end if
#ifndef _NOMPI_
    call MPI_BCAST(variable,1,MPI_REAL8,0,MPI_COMM_WORLD,i)
#endif
  END SUBROUTINE IOTools_readReal8Keyword


  SUBROUTINE IOTools_readReal8Keywords( keyword, variables, unit_in )
    implicit none
    character(*),intent(IN) :: keyword
    real(8),intent(INOUT) :: variables(:)
    integer,optional,intent(IN) :: unit_in
    character(10) :: cbuf,ckey
    integer :: i,unit
    call check_init
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
999    continue
    end if
#ifndef _NOMPI_
    call MPI_BCAST(variables,size(variables),MPI_REAL8,0,MPI_COMM_WORLD,i)
#endif
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
    call check_init
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
999    continue
    end if
#ifndef _NOMPI_
    call MPI_BCAST(variable1,1,MPI_INTEGER,0,MPI_COMM_WORLD,i)
    call MPI_BCAST(variable2,len(variable2),MPI_CHARACTER,0,MPI_COMM_WORLD,i)
#endif
  END SUBROUTINE IOTools_readIntegerString


  SUBROUTINE IOTools_findKeyword( keyword, hasKeyword, unit_out, flag_bcast )
    implicit none
    character(*),intent(IN) :: keyword
    logical,intent(OUT) :: hasKeyword
    integer,optional,intent(OUT) :: unit_out
    logical,optional,intent(IN) :: flag_bcast
    integer :: i,unit
    character(10) :: cbuf,ckey
    call check_init
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
999    continue
    end if
    if ( present(flag_bcast) ) then
#ifndef _NOMPI_
       call MPI_BCAST( hasKeyword, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, i )
#endif
    end if
  END SUBROUTINE IOTools_findKeyword


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
