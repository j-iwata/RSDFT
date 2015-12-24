MODULE kinetic_variables

  implicit none

  PRIVATE
  PUBLIC :: Md, coef_lap0, coef_lap, coef_nab &
           ,coef_nabk, const_k2, zcoef_kin, ggg &
           ,flag_nab, flag_n12, flag_n23, flag_n31
  PUBLIC :: wk
  PUBLIC :: SYStype
  PUBLIC :: read_kinetic, read_oldformat_kinetic
  PUBLIC :: coef_kin

  integer :: SYStype=0
  integer :: Md
  logical :: flag_nab, flag_n12, flag_n23, flag_n31
  real(8) :: coef_lap0, ggg(6)
  real(8),allocatable :: coef_lap(:,:), coef_nab(:,:)
  real(8),allocatable :: coef_nabk(:,:,:), const_k2(:)
  real(8),allocatable :: coef_kin(:)
  complex(8),allocatable :: zcoef_kin(:,:,:)
#ifdef _DRSDFT_
  real(8),allocatable :: wk(:,:,:,:)
#else
  complex(8),allocatable :: wk(:,:,:,:)
#endif

CONTAINS

  SUBROUTINE read_kinetic(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i
    character(7) :: cbuf,ckey
    Md = 6
    SYStype = 0
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:2) == "MD" ) then
             backspace(unit)
             read(unit,*) cbuf,Md
          else if ( ckey(1:7) == "SYSTYPE" ) then
             backspace(unit)
             read(unit,*) cbuf,SYStype
          end if
       end do
999    continue
       write(*,*) "Md =",Md
       write(*,*) "SYStype =",SYStype
    end if
    call send_kinetic(0)
  END SUBROUTINE read_kinetic

  SUBROUTINE read_oldformat_kinetic(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    if ( rank == 0 ) then
       read(unit,*) Md, SYStype
       write(*,*) "Md =",Md
       write(*,*) "SYStype =",SYStype
    end if
    call send_kinetic(0)
  END SUBROUTINE read_oldformat_kinetic

  SUBROUTINE send_kinetic(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(Md,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(SYStype,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_kinetic

END MODULE kinetic_variables
