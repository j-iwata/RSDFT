PROGRAM Real_Space_DFT

  use global_variables
  use parameters_module
  use func2gp_module
  use band_module
  use dos_module, only: calc_dos
  use hamiltonian_matrix_module
  use rtddft_mol_module
  use omp_variables, only: init_omp
  use test_rtsol_module
  use ps_nloc_initiate_module
  use ps_getDij_module
  use io_tools_module, only: init_io_tools, IOTools_readIntegerKeyword
  use lattice_module
  use ffte_sub_module, only: init_ffte_sub
  use fftw_module
  use vdw_grimme_module
  use efield_module

  implicit none
  integer,parameter :: unit_input_parameters = 1
  integer,parameter :: unit_atomic_coordinates = 970
  real(8) :: ct0,ct1,et0,et1,exc_tmp,eh_tmp,eion_tmp,tmp,shift_factor
  integer :: i,n,k,s,iter,m,ierr,i1,i2,i3,m1,m2,m3,j,mm1,mm2,mm3,info
  real(8),allocatable :: force(:,:),forcet(:,:),vtmp(:)
  type(lattice) :: aa_obj, bb_obj
  logical,parameter :: recalc_esp=.true.
  real(8) :: Etot, Ehwf
  integer :: info_level=0

! --- start MPI ---

  call start_mpi_parallel

! --- global time counter start ---

  call global_watch(.false.)

! --- info ---

  call open_info(myrank)

! --- init_io_tools ---

  call init_io_tools( myrank, unit_input_parameters )

! --- DISP_SWITCH ---

  DISP_SWITCH = (myrank==0)
  disp_switch_parallel = (myrank==0)

  call check_disp_switch( DISP_SWITCH, 1 )

  call IOTools_readIntegerKeyword( "INFOLEVEL", info_level )

  call check_disp_length( info_level, 1 )
     
! --- input parameters ---

  call read_atom( myrank, unit_atomic_coordinates, aa_obj )

  call read_parameters

! ---  Type of System ( RS-SOL or RS-MOL ) ---

  call IOTools_readIntegerKeyword( "SYSTYPE", Systype )

! --- Lattice ---

  call write_border( 0, " main( aa & bb )(start)" )

  call init_aa( aa_obj )

  if ( SYStype == 0 ) then
     call convert_to_aa_coordinates_atom( aa_obj, aa_atom )
  else if ( SYStype == 1 ) then
     call convert_to_xyz_coordinates_atom( aa_obj, aa_atom )
  end if

  call backup_aa_lattice( aa_obj )

  call write_border( 0, " main( aa & bb )(end)" )

! --- Pseudopotential ---

  call read_pseudopot( Nelement, myrank )

! --- info atoms ---

  call write_info_atom( Zps, file_ps )

! --- count the total # of electrons & check # of bands ---

  call count_electron

  call check_Nband_electron

! --- init_force ---

  call init_force

! --- R-space Grid ---

  call Init_Rgrid( SYStype, Md, unit=2 )

! --- Reciprocal Lattice ---

  call construct_bb(aa)
  call get_reciprocal_lattice( aa_obj, bb_obj )

! --- G-space Grid ---

  call Init_Ggrid( Ngrid, bb, Hgrid, disp_switch )

! --- Test ( Egg-Box Effect ) ---

  if ( iswitch_test == 2 ) then
     aa_atom(1,1) = aa_atom(1,1) + Hgrid(1)*0.5d0/ax
     if ( disp_switch ) then
        write(*,*) "--- EGG BOX TEST !!! ---"
        write(*,*) aa_atom(1,1),aa_atom(1,1)*ax,Hgrid(1)*0.5d0
     end if
  end if

! --- Symmetry ---

  call init_symmetry( Ngrid,dV,aa,bb, Natom,ki_atom,aa_atom )

! --- Brillouin Zone sampling ---

  call generate_bz

  if ( myrank == 0 ) call write_info_bz( bb )

! --- initial set up for parallel computation ---

!  call test_bcast

  call init_scalapack( Nband )

  call init_parallel( Ngrid, Nband, Nbzsm, Nspin )

  call InitParallel_Rgrid

  call InitParallel_Ggrid( nprocs, myrank )

  call prep_symmetry( Igrid )

! --- initialization for thread-parallel computation ---

  call init_omp( Igrid(1,1),Igrid(2,1),Igrid(1,2),Igrid(2,2) &
                ,Igrid(1,3),Igrid(2,3),Igrid(1,0),Igrid(2,0) &
                ,SYStype, disp_switch )

! --- Initialization for FFT ---

  if ( SYStype == 0 ) then

     call init_ffte_sub(Igrid(1,1:3),Ngrid(1:3),node_partition(1:3),comm_grid)

     call init_fftw( Ngrid(1:3), node_partition(1:3), comm_grid, myrank_g )

  end if

!- FD boundary set -

  call init_bcset( Md, SYStype )

! --- kinetic energy oprator coefficients ---

  call init_kinetic( aa, bb, Nbzsm, kbb, Hgrid, Igrid, MB_d, DISP_SWITCH )

! --- ??? ---

  call set_array_bound

! --- Pseudopotential, initial density, and partial core correction ---

!  call read_pseudopot( Nelement, myrank )


!-------- init density 

  call init_density(Nelectron,dV)

!----------------------- SOL sol -----

  if ( SYStype == 0 ) then

     call init_ps_local
     call init_ps_pcc

     call construct_strfac  !----- structure factor

     call construct_ps_local
     call construct_ps_pcc
     call construct_ps_initrho( rho )

     call destruct_strfac   !----- structure factor

     call ps_nloc_initiate( Gcut )

!----------------------- MOL mol -----

  else if ( SYStype == 1 ) then

     call init_ps_local_mol(Gcut)
     call init_ps_pcc_mol
     call init_ps_initrho_mol

     call Construct_RgridMol(Igrid)

     call construct_ps_local_mol
     call construct_ps_pcc_mol
     call construct_ps_initrho_mol
     call normalize_density

     call ps_nloc2_init(Gcut)
     call prep_ps_nloc2_mol

     call ConstructBoundary_RgridMol(Md,Igrid)

  end if

! --- External electric field (included in the local ionic potential) ---

  call sawtooth_efield( Vion )

! --- symmetrize density ---

  call sym_rho( ML_0, ML_1, Nspin, MSP_0, MSP_1, rho )

!-------------------- Hamiltonian Test

  if ( iswitch_test == 1 ) then
     call test_hpsi2( 10 )
     goto 900
  end if

!--------------------

! --- Ion-Ion ---

  call init_eion( SYStype, disp_switch )

! --- Initialization of subspace diagonalization ---

  call init_subspace_diag( Nband )

! --- Initial wave functions ---

  call init_wf( SYStype )

  do s=MSP_0,MSP_1
  do k=MBZ_0,MBZ_1
     call gram_schmidt(1,Nband,k,s)
  end do
  end do

! --- Initial occupation ---

  call init_occ_electron(Nelectron,Ndspin,Nbzsm,weight_bz,occ)

  if ( DISP_SWITCH ) then
     write(*,'(a60," main")') repeat("-",60)
     write(*,*) "Natom    =",Natom
     write(*,*) "Nelement =",Nelement
     write(*,*) "Nband =",Nband
     write(*,*) "Nspin =",Nspin
     write(*,*) "Nbzsm =",Nbzsm,MBZ
     write(*,*) "Nelectron =",Nelectron
     write(*,*) "Next_electron =",Next_electron
     write(*,*) "Ndspin,Nfixed =",Ndspin,Nfixed
     write(*,*) "Zps   =",Zps(1:Nelement)
     write(*,*) "sum(occ)=",sum(occ)
     if ( Nspin == 2 ) then
        write(*,*) "sum(occ(up))  =",sum(occ(:,:,1))
        write(*,*) "sum(occ(down))=",sum(occ(:,:,Nspin))
     endif
     do n=max(1,nint(Nelectron/2)-10),min(nint(Nelectron/2)+10,Nband)
        do k=1,Nbzsm
           write(*,*) n,k,(occ(n,k,s),s=1,Nspin)
        end do
     end do
  end if

! --- Initial setup for Hybrid XC functional ---

  call init_xc_hybrid( ML_0, ML_1, Nelectron, Nspin, Nband &
       , MMBZ, Nbzsm, MBZ_0, MBZ_1, MSP, MSP_0, MSP_1, MB_0, MB_1 &
       , kbb, bb, Va, SYStype, np_fkmb, disp_switch )

! --- Initial Potential ---

  call init_hartree( Igrid, Nspin, Md, SYStype )
  call calc_hartree( ML_0, ML_1, MSP, rho )

  call calc_xc

  call init_localpot

  do s=MSP_0,MSP_1
     Vloc(:,s) = Vion(:) + Vh(:) + Vxc(:,s)
  end do

!-------------------- Real-Time Test

  if ( iswitch_test == 3 ) then
     call test_rtsol
     goto 900
  end if

! --- Read previous w.f. , density , potentials ---

  call read_data(disp_switch)

  call getDij

! the following GS should be performed when MB1_tmp is smaller than Nband,
! otherwise not necessary

  do s=MSP_0,MSP_1
  do k=MBZ_0,MBZ_1
     call gram_schmidt(1,Nband,k,s)
  end do
  end do

! ---

  if ( GetParam_IO(1) == 1 .or. GetParam_IO(1) >= 3 ) then
     call control_xc_hybrid(1)
     if ( disp_switch ) write(*,*) "iflag_hybrid=",iflag_hybrid
  end if

! ---
! The followings are just to get H and XC energies,
! but potentials are also recalculated with new rho.

  call calc_hartree(ML_0,ML_1,MSP,rho)
  call calc_xc
  do s=MSP_0,MSP_1
     Vloc(:,s) = Vion(:) + Vh(:) + Vxc(:,s)
  end do

! ---

  call getDij

! --- init_vdW_Grimme ---

  call init_vdw_grimme( aa, ki_atom, zn_atom )
  call calc_E_vdw_grimme( aa_atom )

! --- total energy ---

  call calc_with_rhoIN_total_energy( Ehwf )
  call calc_total_energy( recalc_esp, Etot )

! ---

  if ( iswitch_dos == 1 ) then

     call control_xc_hybrid(1)
     call calc_dos( ierr )
     goto 900

  else

     call Init_IO( "sweep" )
     call calc_sweep( disp_switch, ierr )
     call Init_IO( "" )
     if ( ierr == -1 ) goto 900

  end if

! ---

  select case( iswitch_scf )
  case( 1 )
     call calc_scf( disp_switch, ierr, tol_force_in=feps )
     if ( ierr < 0 ) goto 900
     call calc_total_energy( recalc_esp, Etot, 6 )
  case( 2 )
     call calc_scf_chefsi( Diter_scf_chefsi, ierr, disp_switch )
     if ( ierr < 0 ) goto 900
     call calc_total_energy( recalc_esp, Etot, 6 )
  case( -1 )
     if ( nprocs == 1 ) then
        call construct_hamiltonian_matrix( Ngrid(0) )
     end if
     goto 900
  end select

! ---

!
! --- BAND ---
!

  if ( iswitch_band == 1 ) then
     call control_xc_hybrid(1)
     call Init_IO( "band" )
     call band(nint(Nelectron*0.5d0),disp_switch)
     call Init_IO( "" )
  end if

!
! --- Force test, atomopt, CPMD ---
!

  select case( iswitch_opt )
  case( -1 )

     call test_force(SYStype)

  case( 1,2 ) ! --- atomopt ---

     call atomopt(iswitch_opt,disp_switch)

     call calc_total_energy( recalc_esp, Etot, 6 )

  end select

!
! --- TDDFT ---
!
  select case( iswitch_tddft )
  case( 0 )

  case( 1,2 )

     select case( SYStype )
     case( 1 )
        call init_rtddft_mol( 1, myrank )
        call rtddft_mol( iswitch_tddft )
        goto 900
     case default
        write(*,*) "real-time tddft is available only for rsmol"
        goto 900
     end select

  case default

     write(*,'(1x,"iswitch_tddft=",i2," is not available")') iswitch_tddft

  end select

! --- finalize ---

  if ( DISP_SWITCH ) write(*,*) "END_PROGRAM : MAIN" 

900 continue

  call global_watch(disp_switch)
  if ( SYStype == 0 ) call finalize_fftw
  call close_info
  call end_mpi_parallel

END PROGRAM Real_Space_DFT
