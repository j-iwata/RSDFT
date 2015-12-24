MODULE eion_mol_module

  use pseudopot_module, only: Zps
  use atom_module, only: Natom,ki_atom,aa_atom
  use parallel_module

  implicit none

  PRIVATE
  PUBLIC :: calc_eion_mol

  integer,allocatable :: icnt(:),idis(:)
  logical :: flag_init=.true.

CONTAINS

  SUBROUTINE calc_eion_mol(Eion)
    implicit none
    real(8),intent(OUT) :: Eion
    real(8) :: r
    integer :: i,j,ik,jk,ierr

    Eion=0.0d0

    if ( flag_init ) then
       allocate( icnt(0:nprocs-1) ) ; icnt=0
       allocate( idis(0:nprocs-1) ) ; idis=0
       do i=1,Natom
          j=mod(i-1,nprocs)
          icnt(j)=icnt(j)+1
       end do
       do j=0,nprocs-1
          idis(j)=sum(icnt(0:j))-icnt(j)
       end do
       flag_init=.false.
    end if

    do i=idis(myrank)+1,idis(myrank)+icnt(myrank)

       ik=ki_atom(i)

       do j=1,i-1

          jk=ki_atom(j)

          r=sqrt( ( aa_atom(1,i) - aa_atom(1,j) )**2 &
                 +( aa_atom(2,i) - aa_atom(2,j) )**2 &
                 +( aa_atom(3,i) - aa_atom(3,j) )**2 )

          Eion = Eion + Zps(ik)*Zps(jk)/r

       end do ! j

    end do ! i

    call mpi_allreduce &
         (MPI_IN_PLACE,Eion,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

  END SUBROUTINE calc_eion_mol

END MODULE eion_mol_module
