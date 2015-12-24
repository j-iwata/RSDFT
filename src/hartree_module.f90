MODULE hartree_module

  use hartree_variables
  use hartree_sol_module
  use hartree_sol_ffte_module
  use hartree_mol_module
 !use esm_hartree_module
  use hartree_ene_module

  implicit none

  PRIVATE
  PUBLIC :: calc_hartree, init_hartree

  integer :: SYStype=0
  integer :: Nspin
  integer :: Md
  integer :: ML0, ML1

CONTAINS


  SUBROUTINE init_hartree( Igrid, Nspin_in, Md_in, SYStype_in )
    implicit none
    integer,intent(IN) :: Igrid(2,0:3), Nspin_in, Md_in, SYStype_in

    call write_border( 0, " init_hartree(start)" )

    ML0     = Igrid(1,0)
    ML1     = Igrid(2,0)
    Nspin   = Nspin_in
    Md      = Md_in
    SYStype = SYStype_in

    allocate( Vh(ML0:ML1) )
    Vh=0.0d0

    select case( SYStype )
    case(1)
       call init_hartree_mol(Md)
    !case(2)
    !   call init_esm_hartree(Md)
    end select

    call write_border( 0, " init_hartree(end)" )

  END SUBROUTINE init_hartree


  SUBROUTINE calc_hartree(n1,n2,n3,rho)
    implicit none
    integer,intent(IN) :: n1,n2,n3
    real(8),intent(IN) :: rho(n1:n2,n3)
    real(8),allocatable :: trho(:)
    integer :: i,s

    call write_border( 1, " calc_hartree(start)" )

    select case(SYStype)
    case default

       call calc_hartree_sol(n1,n2,n3,rho)

    case(1)

       allocate( trho(n1:n2) )
!$OMP parallel do
       do i=n1,n2
          trho(i) = rho(i,1)
       end do
!$OMP end parallel do
       do s=2,n3
!$OMP parallel do
          do i=n1,n2
             trho(i) = trho(i) + rho(i,s)
          end do
!$OMP end parallel do
       end do
       call calc_hartree_mol(n1,n2,1,trho,Vh,E_hartree)
       call calc_hartree_ene( trho, Vh, E_hartree )
       deallocate( trho )

!    case(3)

!       call calc_esm_hartree(n1,n2,n3,rho,Vh,E_hartree)

    end select

    call write_border( 1, " calc_hartree(end)" )

  END SUBROUTINE calc_hartree


END MODULE hartree_module
