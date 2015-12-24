MODULE aa_module

  use lattice_module, only: lattice, construct_aa_lattice
  use io_tools_module

  implicit none

  PRIVATE
  PUBLIC :: init_aa, read_aa
  PUBLIC :: construct_aa
  PUBLIC :: get_org_aa
  PUBLIC :: set_org_aa

  real(8),PUBLIC :: ax
  real(8),PUBLIC :: Va
  real(8),PUBLIC :: aa(3,3)

  real(8) :: ax_org
  real(8) :: aa_org(3,3)

CONTAINS

  SUBROUTINE init_aa( aa_obj )
    implicit none
    type(lattice),intent(INOUT) :: aa_obj
    if ( ax == 0.0d0 ) then
       ax = aa_obj%LatticeConstant
       aa = aa_obj%LatticeVector/ax
    else
       aa_obj%LatticeConstant = ax
       aa_obj%LatticeVector = aa
       call construct_aa_lattice( aa_obj )
    end if
    call set_org_aa( ax, aa )
    call construct_aa
  END SUBROUTINE init_aa


  SUBROUTINE read_aa
    implicit none
    ax=0.0d0
    aa(:,:)=0.0d0
    aa(1,1)=1.0d0
    aa(2,2)=1.0d0
    aa(3,3)=1.0d0
    call IOTools_readReal8Keyword( "AX", ax )
    call IOTools_readReal8Keywords( "A1", aa(:,1) )
    call IOTools_readReal8Keywords( "A2", aa(:,2) )
    call IOTools_readReal8Keywords( "A3", aa(:,3) )
  END SUBROUTINE read_aa


  SUBROUTINE construct_aa
    implicit none
    aa(:,:)=ax*aa(:,:)
    Va = aa(1,1)*aa(2,2)*aa(3,3)+aa(1,2)*aa(2,3)*aa(3,1) &
        +aa(1,3)*aa(2,1)*aa(3,2)-aa(1,3)*aa(2,2)*aa(3,1) &
        -aa(1,2)*aa(2,1)*aa(3,3)-aa(1,1)*aa(2,3)*aa(3,2)
  END SUBROUTINE construct_aa


  SUBROUTINE set_org_aa(ax_in,aa_in)
    implicit none
    real(8),intent(IN) :: ax_in,aa_in(3,3)
    ax_org=ax_in
    aa_org(:,:)=aa_in(:,:)
  END SUBROUTINE set_org_aa


  SUBROUTINE get_org_aa(ax_out,aa_out)
    implicit none
    real(8),intent(OUT) :: ax_out,aa_out(3,3)
    ax_out=ax_org
    aa_out(:,:)=aa_org(:,:)
  END SUBROUTINE get_org_aa


END MODULE aa_module
