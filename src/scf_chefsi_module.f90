MODULE scf_chefsi_module

  use aa_module
  use parallel_module
  use electron_module
  use localpot_module
  use mixing_module
  use xc_hybrid_module
  use xc_module
  use hartree_variables, only: Vh
  use hartree_module, only: calc_hartree
  use ps_local_module
  use bz_module
  use wf_module
  use cg_module
  use array_bound_module
  use gram_schmidt_module
  use io_module
  use total_energy_module
  use fermi_module
  use subspace_diag_module
  use esp_gather_module
  use density_module
  use watch_module
  use ggrid_module, only: Ecut
  use rgrid_module, only: dV, Ngrid
  use esp_calc_module

  !use localpot2_variables, only: vloc_dense,vloc_dense_old,rho_nl &
  !                              ,vxc_nl,vh_nl,vion_nl
  !use localpot2_module, only: flag_localpot2, test2_localpot2
  !use localpot2_density_module, only: localpot2_density
  !use localpot2_vh_module, only: localpot2_vh
  !use localpot2_xc_module, only: localpot2_xc
  !use localpot2_ion_module, only: localpot2_calc_eion
  !use localpot2_te_module, only: localpot2_te, diff_Etot_lpot2

  use ChebyshevFilter_module

  implicit none

  PRIVATE
  PUBLIC :: calc_scf_chefsi, read_scf_chefsi, Diter_scf_chefsi

  integer :: Diter_scf_chefsi = 100
  integer :: Ndiag            = 1
  logical :: second_diag      =.false.
  real(8) :: scf_conv(2)      = 0.0d0
  real(8) :: fmax_conv        = 0.0d0
  real(8) :: etot_conv        = 0.0d0

CONTAINS


  SUBROUTINE read_scf_chefsi( rank, unit )
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i,ierr
    character(8) :: cbuf,ckey
    scf_conv(1)=1.d-15
    scf_conv(2)=0.0d0
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:7) == "SCFCONV" ) then
             backspace(unit)
             read(unit,*) cbuf,scf_conv
          else if ( ckey(1:8) == "FMAXCONV" ) then
             backspace(unit)
             read(unit,*) cbuf,fmax_conv
          else if ( ckey(1:8) == "ETOTCONV" ) then
             backspace(unit)
             read(unit,*) cbuf,etot_conv
          else if ( ckey(1:5) == "DITER" ) then
             backspace(unit)
             read(unit,*) cbuf,Diter_scf_chefsi
          else if ( ckey(1:5) == "NDIAG" ) then
             backspace(unit)
             read(unit,*) cbuf,Ndiag,second_diag
          end if
       end do
999    continue
       write(*,*) "scf_conv         =",scf_conv
       write(*,*) "fmax_conv        =",fmax_conv
       write(*,*) "etot_conv        =",etot_conv
       write(*,*) "Diter_scf_chefsi =",Diter_scf_chefsi
       write(*,*) "Ndiag            =",Ndiag
       write(*,*) "second_diag      =",second_diag
    end if
    call mpi_bcast( scf_conv  ,2,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(fmax_conv  ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(etot_conv  ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Diter_scf_chefsi,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Ndiag      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(second_diag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
!
    call Init_ChebyshevFilter( rank, unit )
!
  END SUBROUTINE read_scf_chefsi


  SUBROUTINE calc_scf_chefsi( Diter, ierr_out, disp_switch )
    implicit none
    integer,intent(IN)  :: Diter
    integer,intent(OUT) :: ierr_out
    logical,intent(IN) :: disp_switch
    integer :: iter,s,k,n,m,ierr,idiag
    integer :: ML01,MSP01,ib1,ib2
    real(8) :: ct0,et0,ct1,et1
    logical :: flag_exit,flag_end,flag_conv
    real(8) :: Etot, Ehwf

    flag_end  = .false.
    flag_exit = .false.
    flag_conv = .false.
    ierr_out  = 0

    ML01      = ML_1-ML_0+1
    MSP01     = MSP_1-MSP_0+1
    ib1       = max(1,nint(Nelectron/2)-20)
    ib2       = min(nint(Nelectron/2)+80,Nband)

    call init_mixing(ML01,MSP,MSP_0,MSP_1,comm_grid,comm_spin &
                    ,dV,rho(ML_0,MSP_0),Vloc(ML_0,MSP_0) &
                    ,ir_grid,id_grid,myrank)

    allocate( esp0(Nband,Nbzsm,Nspin) ) ; esp0=0.0d0

    do iter=1,Diter

       if ( disp_switch ) then
          write(*,'(a40," scf_chefsi_iter=",i4)') repeat("-",40),iter
       end if

       call watch(ct0,et0)

       esp0=esp

       do s=MSP_0,MSP_1
       do k=MBZ_0,MBZ_1

          if ( iter == 1 ) then
             call watcht(disp_switch,"",0)
             call subspace_diag(k,s)
             call watcht(disp_switch,"diag",1)
          end if

          do idiag=1,Ndiag

             if ( disp_switch ) then
                write(*,'(a5," idiag=",i4)') repeat("-",5),idiag
             end if

             call watcht(disp_switch,"",0)

             call ChebyshevFilter( k,s,MB_0,MB_1 )

             call watcht(disp_switch,"chef",1)

             call gram_schmidt(1,Nband,k,s)

             call watcht(disp_switch,"gs  ",1)

             if ( second_diag .or. idiag < Ndiag ) then
                call subspace_diag(k,s)
                call watcht(disp_switch,"diag",1)
             else if ( idiag == Ndiag ) then
                call esp_calc(k,s,unk(ML_0,MB_0,k,s) &
                             ,ML_0,ML_1,MB_0,MB_1,esp(MB_0,k,s))
                call watcht(disp_switch,"esp_calc",1)
             end if

          end do ! idiag

       end do ! k
       end do ! s

       call watcht(disp_switch,"    ",0)

       call esp_gather(Nband,Nbzsm,Nspin,esp)

       call watcht(disp_switch,"esp_gather",1)

       call calc_fermi(iter,Nfixed,Nband,Nbzsm,Nspin,Nelectron,Ndspin &
                      ,esp,weight_bz,occ,disp_switch)

       call watcht(disp_switch,"fermi",1)

       call calc_with_rhoIN_total_energy( Ehwf )

       call watcht(disp_switch,"harris",1)

! ---
       call calc_density ! n_out
       call calc_hartree(ML_0,ML_1,MSP,rho)
       call control_xc_hybrid(1)
       call calc_xc
       call control_xc_hybrid(2)
       call calc_total_energy( .false., Etot )
! ---

       call watcht(disp_switch,"etot",1)

! ---
       do s=MSP_0,MSP_1
          Vloc(:,s) = Vion(:) + Vh(:) + Vxc(:,s)
       end do

       call perform_mixing( ML01, MSP_0, MSP_1, rho(ML_0,MSP_0) &
            ,Vloc(ML_0,MSP_0), disp_switch )

       if ( mod(imix,2) == 0 ) then
          call normalize_density
          m=(ML_1-ML_0+1)*(MSP_1-MSP_0+1)
          call mpi_allgather &
               (rho(ML_0,MSP_0),m,mpi_real8,rho,m,mpi_real8,comm_spin,ierr)
          call calc_hartree(ML_0,ML_1,MSP,rho)
          call calc_xc
          do s=MSP_0,MSP_1
             Vloc(:,s) = Vion(:) + Vh(:) + Vxc(:,s)
          end do
       end if
       call control_xc_hybrid(1)
! ---

       call watcht(disp_switch,"mixing",1)

!       if ( flag_localpot2 ) then
!          call sub_localpot2_scf( disp_switch )
!          flag_conv=.false.
!          if ( abs(diff_etot_lpot2) < 1.d-10 ) flag_conv=.true.
!       end if

       call write_info_scf( ib1, ib2, iter, disp_switch )

       call watch(ct1,et1)
       if ( disp_switch ) write(*,*) "time(scf)",ct1-ct0,et1-et0

       call global_watch(.false.,flag_end)

       flag_exit = (flag_conv.or.flag_end)

       call watcht(disp_switch,"",0)
       call write_data(disp_switch,flag_exit)
       call watcht(disp_switch,"io",1)

       if ( flag_exit ) exit

    end do ! iter

    if ( myrank == 0 ) write(*,*) "------------ SCF result ----------"
    call write_info_scf( 1, Nband, iter, myrank==0 )

    deallocate( esp0 )

    ierr_out = iter

    if ( flag_end ) then
       ierr_out = -1
       if ( myrank == 0 ) write(*,*) "flag_end=",flag_end
       return
    end if

    if ( iter > Diter ) then
       ierr_out = -2
       if ( myrank == 0 ) write(*,*) "scf not converged"
    end if

    call gather_wf

  END SUBROUTINE calc_scf_chefsi


  SUBROUTINE write_info_scf( ib1, ib2, iter, disp_switch )
    implicit none
    integer,intent(IN) :: ib1, ib2, iter
    logical,intent(IN) :: disp_switch
    integer :: s,k,n,nb1,nb2,i,u(3)
    u(:) = (/ 6, 98, 99 /)
    do i=1,3
       if ( u(i) == 6 ) then
          nb1 = ib1
          nb2 = ib2
       else
          nb1 = 1
          nb2 = Nband
       end if
       if ( u(i) == 6  .and. .not.disp_switch ) cycle
       if ( u(i) /= 6  .and. myrank /= 0 ) cycle
       if ( u(i) == 98 .and. myrank == 0 ) rewind 98
       if ( u(i) == 99 .and. myrank == 0 ) then
          write(u(i),'("AX", f20.15)') ax
          write(u(i),'("A1",3f20.15)') aa(1:3,1)/ax
          write(u(i),'("A2",3f20.15)') aa(1:3,2)/ax
          write(u(i),'("A3",3f20.15)') aa(1:3,3)/ax
          write(u(i),'("VA", f30.15)') Va
          write(u(i),'("NGRID",3i5,i10)') Ngrid(1:3),Ngrid(0)
          write(u(i),'("ECUT ",f10.5)') Ecut
          write(u(i),'("XC",a10)') XCtype
          write(u(i),*) "sum(occ)=",(sum(occ(:,:,s)),s=1,Nspin)
!          write(u(i),*) "iter,sqerr=",iter,sqerr_out(1:Nspin)
       end if
       write(u(i),'(a4,a6,a20,2a13,1x)') &
            "k","n","esp(n,k,s)","esp_err","occ(n,k,s)"
       do k=1,Nbzsm
       do n=nb1,nb2
          write(u(i),'(i4,i6,2(f20.15,2g13.5,1x))') k,n &
               ,(esp(n,k,s),esp(n,k,s)-esp0(n,k,s),occ(n,k,s),s=1,Nspin)
       end do
       end do
    end do
  END SUBROUTINE write_info_scf


!  SUBROUTINE sub_localpot2_scf( disp_switch )
!    implicit none
!    logical,intent(IN) :: disp_switch
!    integer :: mm1,mm2,mm3
!    real(8) :: eion_tmp,eh_tmp,exc_tmp
!    call localpot2_density( rho_nl )
!    call localpot2_calc_eion( vion_nl, rho_nl, eion_tmp )
!    call localpot2_vh( Ecut, rho_nl, vh_nl, eh_tmp )
!    call localpot2_xc( rho_nl, vxc_nl, exc_tmp )
!    vloc_dense=vion_nl+vh_nl+vxc_nl
!    vloc_dense=beta*vloc_dense+(1.d0-beta)*vloc_dense_old
!    vloc_dense_old=vloc_dense
!    call test2_localpot2( vloc_dense )
!    call localpot2_te( eion_tmp, eh_tmp, exc_tmp, disp_switch )
!  END SUBROUTINE sub_localpot2_scf


END MODULE scf_chefsi_module
