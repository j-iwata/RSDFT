MODULE ps_pcc_mol_module

  use ps_pcc_module

  implicit none

  PRIVATE
  PUBLIC :: init_ps_pcc_mol,construct_ps_pcc_mol

CONTAINS


  SUBROUTINE init_ps_pcc_mol
    implicit none
    flag_pcc_0 = .false.
  END SUBROUTINE init_ps_pcc_mol


  SUBROUTINE construct_ps_pcc_mol
    implicit none
    return
  END SUBROUTINE construct_ps_pcc_mol


END MODULE ps_pcc_mol_module
