MODULE global_variables

  use aa_module
  use bb_module
  use atom_module
  use rgrid_sol_module
  use rgrid_module
  use ggrid_module
  use bz_module
  use fd_module
  use strfac_module
  use pseudopot_module
  use ps_local_module
  use ps_pcc_module
  use ps_initrho_module
  use ps_nloc2_init_module
  use ps_nloc2_module
  use bc_module
  use electron_module
  use density_module
  use parallel_module
  use wf_module
  use localpot_module
  use nonlocal_module
  use eion_module
  use gram_schmidt_module
  use gram_schmidt_t_module
  use hartree_variables
  use hartree_module
  use xc_hybrid_module
  use xc_module
  use scalapack_module
  use subspace_diag_module
  use subspace_diag_la_module
  use subspace_diag_sl_module
  use subspace_mate_sl_module
  use subspace_solv_sl_module
  use subspace_rotv_sl_module
  use kinetic_variables
  use kinetic_module
  use hamiltonian_module
  use cgpc_module
  use cg_module
  use total_energy_module
  use mixing_module
  use esp_gather_module
  use fermi_module
  use watch_module
  use io_module
  use array_bound_module
  use atomopt_module

!  use esm_rgrid_module
!  use esm_rshell_module
!  use esm_cylindrical_test
!  use ps_local_rs_module
!  use esm_genpot_module
!  use esm_kinetic_module

  use rgrid_mol_module
  use ps_local_mol_module, only: init_ps_local_mol,construct_ps_local_mol
  use ps_pcc_mol_module
  use ps_initrho_mol_module
  use ps_nloc2_mol_module
  use bc_mol_module
  use kinetic_mol_module

  use ps_gth_module
  use ps_nloc_mr_module

  use bcast_module

!  use localpot2_variables
!  use localpot2_ion_module
!  use localpot2_density_module
!  use localpot2_vh_module
!  use localpot2_xc_module
!  use localpot2_module
!  use localpot2_te_module
!  use localpot2_Smatrix_module

  use ps_nloc3_module

  use test_hpsi2_module
  use test_force_module

  use info_module

  use init_occ_electron_module

  use esp_calc_module

  use symmetry_module

  use force_module

  use sweep_module
  use scf_module
  use scf_chefsi_module

  implicit none

  integer :: iswitch_scf,iswitch_opt,iswitch_band
  integer :: iswitch_test,iswitch_tddft
  real(8) :: etime_limit
  logical :: disp_switch

END MODULE global_variables
