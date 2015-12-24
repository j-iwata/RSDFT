MODULE atom_module

  implicit none

  PRIVATE
  PUBLIC :: Natom,Nelement,aa_atom,ki_atom,read_atom

  integer :: Natom,Nelement
  integer,allocatable :: ki_atom(:)
  real(8),allocatable :: aa_atom(:,:)

CONTAINS

  SUBROUTINE read_atom(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i
    if ( rank == 0 ) then
       rewind unit
       read(unit,*) Nelement,Natom
       write(*,*) "Nelment,Natom=",Nelement,Natom
    end if
    call send_atom_1(0)
    allocate( aa_atom(3,Natom) ) ; aa_atom=0.d0
    allocate( ki_atom(Natom)   ) ; ki_atom=0
    if ( rank == 0 ) then
       do i=1,Natom
          read(unit,*) ki_atom(i),aa_atom(1:3,i)
       end do
       write(*,'(8x,a7,3a18)') "ki_atom","aa_atom1","aa_atom2","aa_atom3"
       if ( Natom <= 11 ) then
          do i=1,Natom
             write(*,'(1x,i5,2x,i7,3f18.12)') i,ki_atom(i),aa_atom(:,i)
          end do
       else
          do i=1,min(5,Natom)
             write(*,'(1x,i5,2x,i7,3f18.12)') i,ki_atom(i),aa_atom(:,i)
          end do
          write(*,'(1x,10x,".")')
          write(*,'(1x,10x,".")')
          write(*,'(1x,10x,".")')
          do i=Natom-5,Natom
             write(*,'(1x,i5,2x,i7,3f18.12)') i,ki_atom(i),aa_atom(:,i)
          end do
       end if
    end if
    call send_atom_2(0)
  END SUBROUTINE read_atom

  SUBROUTINE send_atom_1(myrank)
    implicit none
    integer,intent(IN) :: myrank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(Natom,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Nelement,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_atom_1

  SUBROUTINE send_atom_2(myrank)
    implicit none
    integer,intent(IN) :: myrank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(ki_atom,Natom,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(aa_atom,3*Natom,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_atom_2

END MODULE atom_module
