MODULE random_initrho_module

  use aa_module, only: aa
  use rgrid_module, only: Ngrid, Igrid
  use pseudopot_module, only: Rps, Zps
  use atom_module, only: aa_atom, ki_atom, Natom
  use parallel_module, only: disp_switch_parallel

  implicit none

  PRIVATE
  PUBLIC :: construct_RandomInitrho

CONTAINS

  SUBROUTINE construct_RandomInitrho( rho )
    implicit none
    real(8),intent(OUT) :: rho(:,:)
    integer :: iatm,ielm,i1,i2,i3,j1,j2,j3,i,msp,s
    real(8) :: a1,a2,a3,x,y,z,r2,rc2,c1,c2,c3,d,d0,Zcharge

    if ( disp_switch_parallel ) write(*,*) "random initial density"

    rho = 0.0d0
    msp = size( rho, 2 )

    c1 = 1.0d0/Ngrid(1)
    c2 = 1.0d0/Ngrid(2)
    c3 = 1.0d0/Ngrid(3)

    Zcharge=0.0d0

    do iatm=1,Natom

       ielm    = ki_atom(iatm)
       rc2     = maxval( Rps(:,ki_atom(ielm)) )
       Zcharge = Zcharge + Zps(ielm)

       do j3=-1,1
       do j2=-1,1
       do j1=-1,1

          a1 = aa_atom(1,iatm) + j1
          a2 = aa_atom(2,iatm) + j2
          a3 = aa_atom(3,iatm) + j3

          i=0
          do i3=Igrid(1,3),Igrid(2,3)
          do i2=Igrid(1,2),Igrid(2,2)
          do i1=Igrid(1,1),Igrid(2,1)

             i = i+1

             x = (i1*c1-a1)*aa(1,1)+(i2*c2-a2)*aa(1,2)+(i3*c3-a3)*aa(1,3)
             y = (i1*c1-a1)*aa(1,1)+(i2*c2-a2)*aa(1,2)+(i3*c3-a3)*aa(1,3)
             z = (i1*c1-a1)*aa(1,1)+(i2*c2-a2)*aa(1,2)+(i3*c3-a3)*aa(1,3)

             r2 = x*x + y*y + z*z

             if ( r2 <= rc2 ) then

                call random_number(d)

                rho(i,1) = rho(i,1) + Zps(ielm)*d

             end if

          end do
          end do
          end do

       end do
       end do
       end do

    end do ! iatm

    do s=1,msp
       rho(:,s) = rho(:,1)
    end do

  END SUBROUTINE construct_RandomInitrho

END MODULE random_initrho_module
