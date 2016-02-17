MODULE electron_module

  use atom_module, only: Natom,Nelement,ki_atom
  use pseudopot_module, only: Zps
  use bz_module, only: Nbzsm,weight_bz

  implicit none

  PRIVATE
  PUBLIC :: Nband,Nspin,Nelectron,Next_electron &
           ,Ndspin,Nfixed &
           ,read_electron,count_electron &
           ,Nelectron_spin, dspin
  PUBLIC :: check_Nband_electron

  integer :: Nband, Nspin, Nfixed
  real(8) :: Nelectron,Next_electron,Ndspin
  integer :: MB_0,MB_1,MSP_0,MSP_1

  integer,parameter :: unit_sc=980
  real(8) :: Nelectron_spin(2)
  real(8),allocatable :: dspin(:)

CONTAINS

  SUBROUTINE read_electron(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i
    character(6) :: cbuf,ckey
    Nband=0
    Nspin=1
    Nfixed=0
    Next_electron=0.d0
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:5) == "NBAND" ) then
             backspace(unit)
             read(unit,*) cbuf,Nband
          else if ( ckey(1:5) == "NSPIN" ) then
             backspace(unit)
             read(unit,*) cbuf,Nspin
          else if ( ckey(1:6) == "NDSPIN" ) then
             backspace(unit)
             read(unit,*) cbuf,Ndspin
          else if ( ckey(1:6) == "NFIXED" ) then
             backspace(unit)
             read(unit,*) cbuf,Nfixed
          else if ( ckey(1:6) == "NEXTE" ) then
             backspace(unit)
             read(unit,*) cbuf,Next_electron
          end if
       end do
999    continue
       write(*,*) "Nband=",Nband
       write(*,*) "Nspin=",Nspin
       write(*,*) "Ndspin=",Ndspin
       write(*,*) "Nfixed=",Nfixed
       write(*,*) "Next_electron=",Next_electron
       if ( Nspin == 1 .and. Ndspin /= 0.d0 ) then
          Ndspin=0.d0
          write(*,*) "Ndspin is replaced to 0.0: Ndspin=",Ndspin
       end if
    end if
    call send_electron(0)
    if ( Ndspin < 0.0d0 ) then
       allocate( dspin(Natom) ) ; dspin=0.0d0
       call Read_SpinConf(rank)
    end if
  END SUBROUTINE read_electron


  SUBROUTINE Read_SpinConf(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: i1,i2,i,ierr
    real(8) :: diff_ele
    include 'mpif.h'
    if ( rank == 0 ) then
       rewind unit_sc
       do i=1,Natom
          read(unit_sc,*,END=90) i1,i2,diff_ele
          dspin(i1:i2)=diff_ele
       end do
90     continue
    end if
    call MPI_BCAST(dspin,Natom,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    Ndspin=sum(dspin)
    if ( rank == 0 ) write(*,*) "Ndspin is replaced: Ndspin=",Ndspin
  END SUBROUTINE Read_SpinConf


  SUBROUTINE send_electron(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(Nband,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Nspin,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Next_electron,1,MPI_REAL8,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Ndspin,1,MPI_REAL8,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Nfixed,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_electron


  SUBROUTINE count_electron
    implicit none
    integer :: i
    call write_border( 0, " count_electron(start)" )
    Nelectron=0.0d0
    do i=1,Natom
       Nelectron = Nelectron + Zps(ki_atom(i))
    end do
    Nelectron = Nelectron + Next_electron
    Nelectron_spin(1) = 0.5d0*Nelectron + 0.5d0*Ndspin
    Nelectron_spin(2) = 0.5d0*Nelectron - 0.5d0*Ndspin
    call write_border( 0, " count_electron(end)" )
  END SUBROUTINE count_electron


  SUBROUTINE check_Nband_electron
    implicit none
    real(8),parameter :: factor=1.5d0
    character(30) :: mesg
    call write_border( 0, " check_Nband_electron(start)" )
    if ( dble(Nband) < 0.5d0*Nelectron ) then
       Nband = nint( 0.5d0*Nelectron * factor )
       Nband = max( Nband, 8 )
       write(mesg,'(1x,"Nband is replaced to ",i8)') Nband
       call write_string( mesg )
    end if
    call write_border( 0, " check_Nband_electron(end)" )
  END SUBROUTINE check_Nband_electron


END MODULE electron_module
