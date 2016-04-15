MODULE efield_module

  use aa_module, only: aa
  use rgrid_variables, only: Igrid, Ngrid
  use io_tools_module

  implicit none

  PRIVATE
  PUBLIC :: sawtooth_efield

  real(8) :: Evec(3)
  integer :: icontrol = 0

CONTAINS


  SUBROUTINE init_efield
    implicit none
    Evec(:)=0.0d0
    call IOTools_readReal8Keywords( "EFIELD", Evec )
    icontrol = -1
    if ( any( Evec /= 0.0d0 ) ) icontrol = 1
  END SUBROUTINE init_efield


  SUBROUTINE sawtooth_efield( v )
    implicit none
    real(8),intent(INOUT) :: v(:)
    real(8) :: ee(3,3),x,y,z
    integer :: i1,i2,i3,i
    if ( icontrol ==  0 ) call init_efield
    if ( icontrol == -1 ) return
    ee(:,1)=aa(:,1)/Ngrid(1)
    ee(:,2)=aa(:,2)/Ngrid(2)
    ee(:,3)=aa(:,3)/Ngrid(3)
    i=0
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       i=i+1
       x=ee(1,1)*i1+ee(1,2)*i2+ee(1,3)*i3
       y=ee(2,1)*i1+ee(2,2)*i2+ee(2,3)*i3
       z=ee(3,1)*i1+ee(3,2)*i2+ee(3,3)*i3
       v(i) = v(i) + x*Evec(1) + y*Evec(2) + z*Evec(3)
    end do
    end do
    end do
  END SUBROUTINE sawtooth_efield


END MODULE efield_module
