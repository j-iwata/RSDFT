MODULE force_local_mol_gth_module

  use atom_module, only: Natom, aa_atom, ki_atom
  use rgrid_module, only: Igrid, dV
  use rgrid_mol_module, only: LL, Hsize
  use parallel_module
  use bberf_module
  use pseudopot_module, only: parloc, Rcloc, Zps

  implicit none

  PRIVATE
  PUBLIC :: calc_force_local_mol_gth

CONTAINS

  SUBROUTINE calc_force_local_mol_gth( trho, force )
    implicit none
    real(8),intent(IN) :: trho(:)
    real(8),intent(INOUT) :: force(:,:)
    integer :: a,ik,i,ML_0,ML_1,ierr
    real(8) :: Rx,Ry,Rz,Rc,C1,C2,C3,C4,f
    real(8) :: x,y,z,r,cnst0,cnst1,cnst2
    real(8),allocatable :: ftmp(:,:)

    ML_0 = Igrid(1,0)
    ML_1 = Igrid(2,0)

    allocate( ftmp(3,Natom) ) ; ftmp=0.0d0

    cnst0 = 2.0d0/sqrt( acos(-1.0d0) )

    do a=1,Natom

       ik = ki_atom(a)
       Rx = aa_atom(1,a)
       Ry = aa_atom(2,a)
       Rz = aa_atom(3,a)

       Rc = Rcloc(ik)

       C1 = parloc(1,ik)
       C2 = parloc(2,ik)
       C3 = parloc(3,ik)
       C4 = parloc(4,ik)

       cnst1 = 1.0d0/Rc
       cnst2 = 1.0d0/(sqrt(2.0d0)*Rc)

       do i=ML_0,ML_1

          x = LL(1,i)*Hsize - Rx
          y = LL(2,i)*Hsize - Ry
          z = LL(3,i)*Hsize - Rz
          r = sqrt( x*x + y*y + z*z )

          if ( r < 1.d-9 ) then

             f = Zps(ik)*cnst0*cnst2**3/6.0d0 &
                  - cnst1**2*C1 + 2.0d0*C2*cnst1**2

          else

             f = Zps(ik)*bberf(cnst2*r)/r**2 &
               - Zps(ik)*cnst0*cnst2*exp(-(cnst2*r)**2)/r &
               - r*cnst1**2*exp(-0.5d0*(cnst1*r)**2)*( C1 &
               + C2*(cnst1*r)**2 + C3*(cnst1*r)**4 + C4*(cnst1*r)**6 ) &
               + exp(-0.5d0*(cnst1*r)**2)*( 2.0d0*C2*cnst1**2*r &
               + 4.0d0*C3*cnst1**4*r**3 + 6.0d0*C4*cnst1**6*r**5 )

             f = f/r

          end if

          ftmp(1,a) = ftmp(1,a) + f*x*trho(i-ML_0+1)
          ftmp(2,a) = ftmp(2,a) + f*y*trho(i-ML_0+1)
          ftmp(3,a) = ftmp(3,a) + f*z*trho(i-ML_0+1)

       end do ! i

    end do ! a

    call mpi_allreduce( MPI_IN_PLACE, ftmp, 3*Natom, MPI_REAL8, MPI_SUM &
                       ,comm_grid, ierr)

    force(:,:) = force(:,:) + ftmp(:,:)*dV

    deallocate( ftmp )

  END SUBROUTINE calc_force_local_mol_gth


END MODULE force_local_mol_gth_module
