MODULE force_ion_mol_module

  use atom_module
  use pseudopot_module

  implicit none

  PRIVATE
  PUBLIC :: calc_force_ion_mol

CONTAINS

  SUBROUTINE calc_force_ion_mol(force)
    implicit none
    real(8),intent(INOUT) :: force(:,:)
    integer :: a,b
    real(8) :: r3,c1,ax,ay,az,bx,by,bz

    force(1:3,1:Natom) = 0.0d0

    do a=1,Natom

       ax = aa_atom(1,a)
       ay = aa_atom(2,a)
       az = aa_atom(3,a)

       do b=1,Natom

          if ( a == b ) cycle

          bx = aa_atom(1,b)
          by = aa_atom(2,b)
          bz = aa_atom(3,b)

          r3 = sqrt( (ax-bx)**2 + (ay-by)**2 + (az-bz)**2 )**3
          c1 = Zps(ki_atom(a))*Zps(ki_atom(b))/r3

          force(1,a) = force(1,a) + c1*(ax-bx)
          force(2,a) = force(2,a) + c1*(ay-by)
          force(3,a) = force(3,a) + c1*(az-bz)

       end do ! b

    end do ! a

    return
  END SUBROUTINE calc_force_ion_mol

END MODULE force_ion_mol_module
