MODULE band_variables

  use io_tools_module

  implicit none

  integer :: nbk,mb_band,mb2_band,maxiter_band,nskip_band
  real(8) :: esp_conv_tol=-1.d0
  integer,allocatable :: nfki(:)
  real(8),allocatable :: ak(:,:)

  integer,parameter :: unit_band_eigv=100
  integer,parameter :: unit_band_ovlp=110
  integer,parameter :: unit_band_dedk=120
  integer,parameter :: unit_band_ufld=130

CONTAINS

  SUBROUTINE read_band
    implicit none
    integer :: i,ierr,unit
    logical :: exist_keyword
    character(8) :: cbuf
    include 'mpif.h'
    call write_border( 0, " read_band(start)" )
    esp_conv_tol = 1.d-5
    maxiter_band = 50
    nskip_band   = 0
    nbk          = 0
    mb_band      = 0
    mb2_band     = 0
    call IOTools_findKeyword( "BAND", exist_keyword, unit_out=unit )
    if ( exist_keyword ) then
       backspace(unit)
       read(unit,*) cbuf,nbk,mb_band,mb2_band,esp_conv_tol,maxiter_band
       allocate( ak(3,nbk+1) ) ; ak=0.0d0
       allocate( nfki(nbk)   ) ; nfki=0
       read(unit,*) ak(1,1:nbk+1)
       read(unit,*) ak(2,1:nbk+1)
       read(unit,*) ak(3,1:nbk+1)
       read(unit,*) nfki(1:nbk)
       write(*,*) "nbk          =",nbk
       write(*,*) "mb_band      =",mb_band
       write(*,*) "mb2_band     =",mb2_band
       write(*,*) "esp_von_tol  =",esp_conv_tol
       write(*,*) "maxiter_band =",maxiter_band
    end if
    call mpi_bcast(nbk,1,mpi_integer,0,mpi_comm_world,ierr)
    if ( .not.allocated(ak) ) then
       allocate( ak(3,nbk+1) ) ; ak=0.0d0
       allocate( nfki(nbk)   ) ; nfki=0
    end if
    call mpi_bcast(ak,size(ak),mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(nfki,size(nfki),mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(mb_band,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(mb2_band,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(esp_conv_tol,1,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(maxiter_band,1,mpi_integer,0,mpi_comm_world,ierr)
    call IOTools_readIntegerKeyword( "SKIPBAND", nskip_band )
    call write_border( 0, " read_band(end)" )
  END SUBROUTINE read_band

END MODULE band_variables
