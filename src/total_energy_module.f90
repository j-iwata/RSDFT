MODULE total_energy_module

  use rgrid_module, only: dV
  use hamiltonian_module
  use hartree_variables, only: Vh, E_hartree
  use xc_module, only: Vxc,E_exchange,E_correlation,Exc,E_exchange_exx
  use eion_module, only: Eewald
  use wf_module, only: unk,esp,occ
  use localpot_module, only: Vloc
  use ps_local_module, only: Vion, const_ps_local
  use density_module, only: rho
  use parallel_module
  use fermi_module, only: efermi,Eentropy
  use array_bound_module, only: ML_0,ML_1,MB,MB_0,MB_1 &
                               ,MBZ,MBZ_0,MBZ_1,MSP,MSP_0,MSP_1
  use fock_module
  use var_sys_parameter, only: pp_kind
  use vdw_grimme_module

  implicit none

  PRIVATE
  PUBLIC :: calc_total_energy
  PUBLIC :: calc_with_rhoIN_total_energy

  integer :: scf_iter_

  real(8) :: Ekin,Eloc,Enlc,Eeig,Eion,Fene,Evdw
  real(8) :: Etot_0=0.d0
  real(8) :: Ekin_0=0.d0
  real(8) :: Eloc_0=0.d0
  real(8) :: Enlc_0=0.d0
  real(8) :: Eeig_0=0.d0
  real(8) :: Eion_0=0.d0
  real(8) :: Ehat_0=0.d0
  real(8) :: Exc_0 =0.d0
  real(8) :: Ex_0  =0.d0
  real(8) :: Ec_0  =0.d0

  real(8) :: efermi_0  =0.d0
  real(8) :: Eentropy_0=0.d0
  real(8) :: Fene_0    =0.d0

  real(8) :: Ehwf
  real(8) :: Ehwf_0  = 0.d0
  real(8) :: Eloc_in = 0.d0
  real(8) :: Ehat_in = 0.d0
  real(8) :: Exc_in  = 0.d0
  real(8) :: Eion_in = 0.d0

  real(8) :: diff_etot = 0.d0

CONTAINS


  SUBROUTINE calc_total_energy( flag_recalc_esp, Etot, unit_in )
    implicit none
    logical,intent(IN) :: flag_recalc_esp
    real(8),intent(INOUT) :: Etot
    integer,optional,intent(IN) :: unit_in
    integer :: i,n,k,s,n1,n2,ierr,nb1,nb2,unit
    real(8) :: s0(4),s1(4),uu,cnst
    real(8),allocatable :: esp0(:,:,:,:),esp1(:,:,:,:)
    real(8),allocatable :: esp0_Q(:,:,:),esp1_Q(:,:,:)
#ifdef _DRSDFT_
    real(8),parameter :: zero=0.d0
    real(8),allocatable :: work(:,:)
    real(8),allocatable :: work00(:,:)
#else
    complex(8),parameter :: zero=(0.d0,0.d0)
    complex(8),allocatable :: work(:,:)
    complex(8),allocatable :: work00(:,:)
#endif
    include 'mpif.h'

    call write_border( 1, " calc_total_energy(start)" )

    Etot_0 = Etot

    Etot = 0.d0
    Ekin = 0.d0
    Eloc = 0.d0
    Enlc = 0.d0
    Eeig = 0.d0
    Eion = 0.d0
    Fene = 0.d0
    Evdw = 0.d0

    n1 = ML_0
    n2 = ML_1

    if ( flag_recalc_esp ) then

       allocate( esp0(MB,MBZ,MSP,4) ) ; esp0=0.d0
       allocate( esp0_Q(MB,MBZ,MSP) ) ; esp0_Q=0.d0
       allocate( work(n1:n2,MB_d)   ) ; work=zero
       allocate( work00(n1:n2,MB_d) ) ; work00=zero

       do s=MSP_0,MSP_1
       do k=MBZ_0,MBZ_1
       do n=MB_0,MB_1,MB_d

          nb1=n
          nb2=min(nb1+MB_d-1,MB_1)

!---------------------------------------------------- kinetic

          work=zero
!$OMP parallel
          call op_kinetic(k,unk(n1,n,k,s),work,n1,n2,nb1,nb2)
!$OMP end parallel
          do i=nb1,nb2
#ifdef _DRSDFT_
          esp0(i,k,s,1)=sum( unk(:,i,k,s)*work(:,i-nb1+1) )*dV
#else
          esp0(i,k,s,1)=sum( conjg(unk(:,i,k,s))*work(:,i-nb1+1) )*dV
#endif
          end do

!---------------------------------------------------- local

          work=zero
!$OMP parallel
          call op_localpot(s,n2-n1+1,nb2-nb1+1,unk(n1,n,k,s),work)
!$OMP end parallel
          do i=nb1,nb2
#ifdef _DRSDFT_
          esp0(i,k,s,2)=sum( unk(:,i,k,s)*work(:,i-nb1+1) )*dV
#else
          esp0(i,k,s,2)=sum( conjg(unk(:,i,k,s))*work(:,i-nb1+1) )*dV
#endif
          end do

!---------------------------------------------------- nonlocal

          if ( pp_kind == 'USPP' ) then

             work=zero
             work00=zero
!$OMP parallel
             call op_nonlocal(k,s,unk(n1,n,k,s),work,n1,n2,nb1,nb2,work00)
!$OMP end parallel
             do i=nb1,nb2
#ifdef _DRSDFT_
                esp0(i,k,s,3)=sum( unk(:,i,k,s)*work(:,i-nb1+1) )*dV
                esp0_Q(i,k,s)=sum( unk(:,i,k,s)*work00(:,i-nb1+1) )*dV
#else
                esp0(i,k,s,3)=sum( conjg(unk(:,i,k,s))*work(:,i-nb1+1) )*dV
                esp0_Q(i,k,s)=sum( conjg(unk(:,i,k,s))*work00(:,i-nb1+1) )*dV
#endif
             end do

          else if ( pp_kind == 'NCPP' ) then

             work=zero
!$OMP parallel
             call op_nonlocal(k,s,unk(n1,n,k,s),work,n1,n2,nb1,nb2)
!$OMP end parallel
             do i=nb1,nb2
#ifdef _DRSDFT_
                esp0(i,k,s,3)=sum( unk(:,i,k,s)*work(:,i-nb1+1) )*dV
#else
                esp0(i,k,s,3)=sum( conjg(unk(:,i,k,s))*work(:,i-nb1+1) )*dV
#endif
             end do

          end if

!---------------------------------------------------- fock

          work=zero
          call op_fock(k,s,n1,n2,n,n,unk(n1,n,k,s),work)
#ifdef _DRSDFT_
          esp0(n,k,s,4)=sum( unk(:,n,k,s)*work(:,1) )*dV
#else
          esp0(n,k,s,4)=sum( conjg(unk(:,n,k,s))*work(:,1) )*dV
#endif

       end do ! n
       end do ! k
       end do ! s

       deallocate( work   )
       deallocate( work00 )

       allocate( esp1(MB,MBZ,MSP,4) )
       allocate( esp1_Q(MB,MBZ,MSP) )

       n=MB*MBZ*MSP*4
       call mpi_allreduce(esp0,esp1,n,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
       esp1=esp1/np_fkmb
       n=MB*MBZ*MSP
       call mpi_allreduce(esp0_Q,esp1_Q,n,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
       esp1_Q=esp1_Q/np_fkmb

       Ekin = sum( occ(:,:,:)*esp1(:,:,:,1) )
       Eloc = sum( occ(:,:,:)*esp1(:,:,:,2) )
        
       if ( pp_kind == 'USPP' ) then

          Enlc = sum( occ(:,:,:)*esp1_Q(:,:,:) )

       else if ( pp_kind == 'NCPP' ) then

          Enlc = sum( occ(:,:,:)*esp1(:,:,:,3) )

       endif

       esp(:,:,:)=esp1(:,:,:,1)+esp1(:,:,:,2)+esp1(:,:,:,3)+esp1(:,:,:,4)

       deallocate( esp1 )
       deallocate( esp0 )
       deallocate( esp1_Q )
       deallocate( esp0_Q )

    end if ! flag_recalc_esp

    Eeig = sum( occ(:,:,:)*esp(:,:,:) )
    cnst = sum( occ(:,:,:) )*const_ps_local

    select case( pp_kind )
    case( "USPP" )

    s0(:)=0.d0
    s1(:)=0.d0
    do s=MSP_0,MSP_1
       do i=n1,n2
          s0(1) = s0(1) + rho(i,s)*Vloc(i,s) ! rho inlculdes Qij(r) terms
          s0(2) = s0(2) + rho(i,s)*Vion(i)
       end do
    end do
    s1(:)=s0(:)*dV

    call mpi_allreduce(s1,s0,2,mpi_real8,mpi_sum,comm_grid,ierr)
    call mpi_allreduce(s0,s1,2,mpi_real8,mpi_sum,comm_spin,ierr)

    case( "NCPP" )

    s0(:)=0.d0
    do s=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1
    do n=MB_0,MB_1
       s1(:)=0.d0
       do i=n1,n2
          uu=abs(unk(i,n,k,s))**2
          s1(1) = s1(1) + uu*Vloc(i,s)
          s1(2) = s1(2) + uu*Vion(i)
          s1(3) = s1(3) + uu*Vh(i)          ! not used?
          s1(4) = s1(4) + uu*Vxc(i,s)       ! not used?
       end do
       s0(:)=s0(:)+occ(n,k,s)*s1(:)
    end do
    end do
    end do
    s1(:)=s0(:)*dV
    call mpi_allreduce(s1,s0,4,mpi_real8,mpi_sum,comm_grid,ierr)
    call mpi_allreduce(s0,s1,4,mpi_real8,mpi_sum,comm_band,ierr)
    call mpi_allreduce(s1,s0,4,mpi_real8,mpi_sum,comm_bzsm,ierr)
    call mpi_allreduce(s0,s1,4,mpi_real8,mpi_sum,comm_spin,ierr)
!    s0(:)=s0(:)*dV/np_fkmb
!    call mpi_allreduce(s0,s1,4,mpi_real8,mpi_sum,MPI_COMM_WORLD,ierr)

    end select

    Eloc = s1(1)
    Eion = s1(2)

    call get_E_vdw_grimme( Evdw )

    Etot = Eeig - Eloc + E_hartree + Exc + Eion + Eewald &
         - 2*E_exchange_exx + Evdw + cnst

    Ehwf = Eeig - Eloc_in + Ehat_in + Exc_in + Eion_in + Eewald &
         - 2*E_exchange_exx + Evdw + cnst

    Fene = Etot - Eentropy

    diff_etot = Etot - Etot_0
!    diff_etot = Etot - Ehwf

    unit=99 ; if ( present(unit_in) ) unit=unit_in
    call write_info_total_energy( Etot, (myrank==0), unit )

!    if ( present(flag_rewind) ) then
!       call write_info_total_energy( disp_switch, flag_rewind )
!    else
!       call write_info_total_energy( disp_switch, .false. )
!    end if
!    if ( disp_switch ) then
!       write(*,'(1x,"Total Energy   =",f16.8,2x,"(Hartree)")') Etot
!       write(*,'(1x,"Harris Energy  =",f16.8,2x,"(Hartree)")') Ehwf
!       write(*,'(1x,"difference    =",g13.5)') Etot-Ehwf
!       write(*,'(1x,"Total (Harris) Energy =",f16.8,2x,"(",f16.8,")" &
!            ,2x,"(Hartree)")') Etot, Ehwf
!       write(*,'(1x,"difference =",g13.5)') Etot-Ehwf
!    end if

    Ekin_0 = Ekin
    Eloc_0 = Eloc
    Enlc_0 = Enlc
    Eion_0 = Eion
    Ehat_0 = E_hartree
    Exc_0  = Exc
    Ex_0   = E_exchange
    Ec_0   = E_correlation
    Eeig_0 = Eeig

    efermi_0   = efermi
    Eentropy_0 = Eentropy
    Fene_0     = Fene

    call write_border( 1, " calc_total_energy(end)" )

  END SUBROUTINE calc_total_energy


  SUBROUTINE calc_with_rhoIN_total_energy( Etot )
    implicit none
    real(8),optional,intent(OUT) :: Etot
    real(8) :: sb(2),rb(2),Eeig_tmp
    integer :: s,ierr
    call write_border( 1, " calc_with_rhoIN_total_energy(start)" )
    sb(:)=0.d0
    do s=MSP_0,MSP_1
       sb(1) = sb(1) + sum(Vloc(:,s)*rho(:,s))
       sb(2) = sb(2) + sum(Vion(:)*rho(:,s))
    end do
    call mpi_allreduce(sb,rb,2,MPI_REAL8,MPI_SUM,comm_grid,ierr)
    call mpi_allreduce(rb,sb,2,MPI_REAL8,MPI_SUM,comm_spin,ierr)
    Eloc_in = sb(1)*dV
    Eion_in = sb(2)*dV
    Ehat_in = E_hartree
    Exc_in  = Exc
    Eeig_tmp=sum( occ(:,:,:)*esp(:,:,:) )
    call get_E_vdw_grimme( Evdw )
    Etot = Eeig_tmp - Eloc_in + Ehat_in + Exc_in + Eion_in + Eewald &
         - 2*E_exchange_exx + const_ps_local*sum(occ) + Evdw
    call write_border( 1, " calc_with_rhoIN_total_energy(end)" )
  END SUBROUTINE calc_with_rhoIN_total_energy


  SUBROUTINE write_info_total_energy( Etot, flag_write, u )
    implicit none
    real(8),intent(IN) :: Etot
    logical,intent(IN) :: flag_write
    integer,intent(IN) :: u
    if ( flag_write ) then
       if ( u == 6 ) then
          call write_border( 0, " total energy" )
       else
          rewind u
       end if
       write(u,*) "Total Energy ",Etot
       write(u,*) "Harris Energy",Ehwf
       write(u,*) "Ion-Ion                    ",Eewald
       write(u,*) "Local Potential            ",Eloc
       write(u,*) "Ion Local Potential        ",Eion
       if ( Enlc /= 0.0d0 ) write(u,*) "Ion Nonlocal Potential     ",Enlc
       if ( Ekin /= 0.0d0 ) write(u,*) "Kinetic Energy             ",Ekin
       write(u,*) "Hartree Energy             ",E_hartree
       write(u,*) "Exchange-Correlation Energy",Exc
       write(u,*) "Exchange Energy            ",E_exchange
       write(u,*) "Correlation Energy         ",E_correlation
       write(u,*) "Sum of eigenvalues         ",Eeig
       write(u,*) "Fermi energy               ",efermi
       if ( u == 6 ) call write_border( 0, "" )
    end if
!    u(:) = (/ 6, 99 /)
!    do i=1,2
!       if ( u(i) == 6 .and. .not.disp_switch ) cycle
!       if ( u(i) /= 6 .and. myrank /= 0 ) cycle
!       if ( u(i) /= 6 .and. myrank == 0 .and. flag_rewind ) rewind u(i)
!       if ( Evdw /= 0.0d0 ) write(u(i),*) '(VDW) ',Evdw
!       if ( const_ps_local /= 0.0d0 ) write(u(i),*) '(cnst)',const_ps_local*sum(occ)
!       write(u(i),*) '(EII) ',Eewald
!       write(u(i),*) '(KIN) ',Ekin, Ekin-Ekin_0
!       write(u(i),*) '(LOC) ',Eloc, Eloc-Eloc_0
!       write(u(i),*) '(NLC) ',Enlc, Enlc-Enlc_0
!       write(u(i),*) '(ION) ',Eion, Eion-Eion_0
!       write(u(i),*) '(HTR) ',E_hartree, E_hartree-Ehat_0
!       write(u(i),*) '(EXC) ',Exc,  Exc-Exc_0
!       write(u(i),*) '(EXG) ',E_exchange, E_exchange-Ex_0
!       write(u(i),*) '(COR) ',E_correlation, E_correlation-Ec_0
!       write(u(i),*) '(EIG) ',Eeig, Eeig-Eeig_0
!!       write(u(i),*) '(HWF) ',Ehwf, Ehwf-Etot
!!       write(u(i),*) '(TOT) ',Etot, Etot_0-Etot
!       write(u(i),*) '(efermi)  ',efermi, efermi-efermi_0
!       !write(u(i),*) '(entropy) ',Eentropy,Eentropy-Eentropy_0
!       !write(u(i),*) '(FreeEne) ',Fene,Fene-Fene_0
!       call flush(u(i))
!    end do
  END SUBROUTINE write_info_total_energy


END MODULE total_energy_module
