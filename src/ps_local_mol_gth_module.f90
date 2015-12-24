MODULE ps_local_mol_gth_module

  use pseudopot_module, only: parloc, Zps, Rcloc
  use ps_local_module, only: Vion
  use atom_module, only: Natom, aa_atom, ki_atom
  use rgrid_mol_module, only: LL, Hsize
  use rgrid_module, only: Igrid
  use bberf_module

  implicit none

  PRIVATE
  PUBLIC :: init_ps_local_mol_gth, construct_ps_local_mol_gth

CONTAINS


  SUBROUTINE init_ps_local_mol_gth
    implicit none
  END SUBROUTINE init_ps_local_mol_gth


  SUBROUTINE construct_ps_local_mol_gth( disp_sw )
    implicit none
    logical,optional,intent(IN) :: disp_sw
    integer :: a,ik,i,ML_0,ML_1
    real(8) :: Rx,Ry,Rz,x,y,z,r,cnst0,cnst1,cnst2
    real(8) :: C1,C2,C3,C4,Rc

    if ( present(disp_sw) ) then
       if ( disp_sw ) then
          write(*,'(a50," construct_ps_local_mol_gth(start)")') repeat("-",50)
       end if
    end if

    ML_0 = Igrid(1,0)
    ML_1 = Igrid(2,0)

    if ( .not. allocated(Vion) ) then
       allocate( Vion(ML_0:ML_1) )
    end if
    Vion=0.0d0

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

       cnst1=1.0d0/Rc
       cnst2=1.0d0/(sqrt(2.0d0)*Rc)

       do i=ML_0,ML_1
          x=LL(1,i)*Hsize-Rx
          y=LL(2,i)*Hsize-Ry
          z=LL(3,i)*Hsize-Rz
          r=sqrt(x*x+y*y+z*z)
          if ( r < 1.d-9 ) then
             Vion(i) = Vion(i) - Zps(ik)*cnst2*cnst0 + C1
          else
             Vion(i) = Vion(i) - Zps(ik)*bberf(cnst2*r)/r &
                  + exp( -0.5d0*(r*cnst1)**2 )*( C1 &
                  + C2*(r*cnst1)**2 + C3*(r*cnst1)**4 + C4*(r*cnst1)**6 )
          end if
       end do ! i

    end do ! a

    if ( present(disp_sw) ) then
       if ( disp_sw ) then
          write(*,'(a50," construct_ps_local_mol_gth(end)")') repeat("-",50)
       end if
    end if

  END SUBROUTINE construct_ps_local_mol_gth


END MODULE ps_local_mol_gth_module
