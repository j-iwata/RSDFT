MODULE fock_cg_module

!  use parallel_module
!  use rgrid_module, only: Ngrid, Igrid
!  use ggrid_module, only: NGgrid, Ecut
!  use bb_module, only: bb
!  use xc_hybrid_module, only: q_fock, R_hf, omega &
!                             ,iflag_lcwpbe, iflag_hse, iflag_hf, iflag_pbe0
!  use watch_module
!  use bz_module, only: kbb
!  use fock_ffte_module

  use hartree_mol_module, only: calc_hartree_mol

  implicit none

  PRIVATE
  PUBLIC :: ct_fock_cg, et_focK_cg, fock_cg

  real(8) :: ct_fock_cg(10),et_fock_cg(10)

CONTAINS


  SUBROUTINE Fock_cg( n1, n2, k, q, trho, tVh, t )

    implicit none

    integer,intent(IN) :: n1,n2,k,q,t
#ifdef _DRSDFT_
    real(8),intent(IN)    :: trho(n1:n2)
    real(8),intent(INOUT) :: tVh(n1:n2)
#else
    complex(8),intent(IN)    :: trho(n1:n2)
    complex(8),intent(INOUT) :: tVh(n1:n2)
#endif
    real(8) :: Edummy

#ifdef _DRSDFT_
    call calc_hartree_mol( n1, n2, 1, trho, tVh, Edummy, tol=1.d-8 )
#endif

    return
  END SUBROUTINE Fock_cg


END MODULE fock_cg_module
