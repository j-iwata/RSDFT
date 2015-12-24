MODULE density_module

  use wf_module
  use parallel_module, only: comm_grid,comm_band,comm_bzsm,comm_spin,ir_grid,id_grid,ir_spin,id_spin,myrank_g,myrank_s,myrank
  use symmetry_module, only: sym_rho
  use basic_type_factory
  use basic_type_methods
  use rhonks_g_module, only: get_rhonks
  use var_sys_parameter, only: pp_kind
  use array_bound_module, only: get_grid_range_local, get_grid_range_globl &
                               ,get_spin_range_local, get_spin_range_globl

  implicit none

  PRIVATE
  PUBLIC :: init_density
  PUBLIC :: normalize_density
  PUBLIC :: calc_density
  PUBLIC :: get_range_density, construct_density_v2
  PUBLIC :: writeDensity

  real(8),allocatable,PUBLIC :: rho(:,:)
  real(8),PUBLIC :: sum_dspin(2)

  integer :: ML_RHO,ML_0_RHO,ML_1_RHO
  integer :: MS_RHO,MS_0_RHO,MS_1_RHO

  real(8) :: Nelectron_RHO,dV_RHO

  include 'mpif.h'

CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE get_range_density( g_range, s_range )
    implicit none
    type( ArrayRange1D ),intent(OUT) :: g_range,s_range
    g_range%head = ML_0_RHO
    g_range%tail = ML_1_RHO
    g_range%size = ML_1_RHO - ML_0_RHO + 1
    s_range%head = 1
    s_range%tail = MS_RHO
    s_range%size = MS_RHO
    g_range%head_global = 1
    g_range%tail_global = ML_RHO
    g_range%size_global = ML_RHO
    s_range%head_global = 1
    s_range%tail_global = MS_RHO
    s_range%size_global = MS_RHO
  END SUBROUTINE get_range_density

  SUBROUTINE get_range_density_v2( g_range, s_range )
    implicit none
    type( pArrayRange1D ),intent(OUT) :: g_range,s_range
    call get_grid_range_local( g_range%local )
    call get_grid_range_globl( g_range%globl )
    call get_spin_range_local( s_range%local )
    call get_spin_range_globl( s_range%globl )
    g_range%alloc = g_range%local
    s_range%alloc = s_range%globl
  END SUBROUTINE get_range_density_v2

  SUBROUTINE construct_density_v2( density )
    implicit none
    type( GSArray_v2 ),intent(OUT) :: density
    call get_range_density_v2( density%g_range, density%s_range )
    call allocateGSArray_v2( density )
  END SUBROUTINE construct_density_v2

  SUBROUTINE init_density(Nelectron,dV)
    implicit none
    real(8),intent(IN) :: Nelectron,dV
    integer :: i,s

    call write_border( 1, " init_density(start)" )

    Nelectron_RHO = Nelectron
    dV_RHO        = dV

    ML_RHO   = sum( ir_grid )
    ML_0_RHO = id_grid(myrank_g) + 1
    ML_1_RHO = id_grid(myrank_g) + ir_grid(myrank_g)
    MS_RHO   = sum( ir_spin )
    MS_0_RHO = id_spin(myrank_s) + 1
    MS_1_RHO = id_spin(myrank_s) + ir_spin(myrank_s)

    if ( .not. allocated(rho) ) then
       allocate( rho(ML_0_RHO:ML_1_RHO,MS_RHO) )
       call random_number(rho)
       call normalize_density
    end if

    call write_border( 1, " init_density(end)" )

  END SUBROUTINE init_density

!-----------------------------------------------------------------------

  SUBROUTINE normalize_density
    implicit none
    real(8) :: c,d
    integer :: ierr
    include 'mpif.h'
    call write_border( 1, " normalize_density(start)" )
    c=sum(rho)*dV_RHO
    call mpi_allreduce(c,d,1,MPI_REAL8,MPI_SUM,comm_grid,ierr)
    c=Nelectron_RHO/d
    rho=c*rho
    call write_border( 1, " normalize_density(end)" )
  END SUBROUTINE normalize_density

!-----------------------------------------------------------------------

  SUBROUTINE calc_density( Ntot )
    ! IN:	unk(:.n.k.s),ML_0,ML_1,MS_0,MS_1
    ! OUT:	rho(n1:n2,s)
    implicit none
    real(8),optional,intent(OUT) :: Ntot(:)
    integer :: n,k,s,i
    integer :: n1,n2,n0
    real(8),allocatable :: rhonks(:)

    call write_border( 1, " calc_density(start)" )

    select case ( pp_kind )
    case ( 'NCPP' )

       rho(:,:)=0.0d0
       do s=MS_0_WF,MS_1_WF
          do k=MK_0_WF,MK_1_WF
             do n=MB_0_WF ,MB_1_WF
                rho(:,s)=rho(:,s)+occ(n,k,s)*abs( unk(:,n,k,s) )**2
             end do
          end do
       end do
       call sym_rho( ML_0_RHO, ML_1_RHO, MS_RHO, MS_0_RHO, MS_1_RHO, rho )
       call reduce_and_gather

    case ( 'USPP' )

       n1 = ML_0_WF
       n2 = ML_1_WF
       n0 = ML_1_WF - ML_0_WF + 1

       allocate( rhonks(n1:n2) ) ; rhonks(:)=0.0d0

       rho(:,:)=0.0d0
       do s=MS_0_WF,MS_1_WF
          do k=MK_0_WF,MK_1_WF
             do n=MB_0_WF,MB_1_WF
                rhonks(:)=0.d0
                call get_rhonks( rhonks,n1,n2,n,k,s )
                rho(:,s) = rho(:,s) + occ(n,k,s)*rhonks(:)
             end do
          end do
       end do

       deallocate( rhonks )

       call reduce_and_gather

    end select

    if ( present(Ntot) ) call calc_spin_density( Ntot )

    call write_border( 1, " calc_density(end)" )

  END SUBROUTINE calc_density

!-----------------------------------------------------------------------

  SUBROUTINE reduce_and_gather
    implicit none
    integer :: n,k,s,m,ierr
    include 'mpif.h'
    m=ML_1_RHO-ML_0_RHO+1
    do s=MS_0_RHO,MS_1_RHO
       call mpi_allreduce(MPI_IN_PLACE,rho(ML_0_RHO,s),m,mpi_real8,mpi_sum,comm_band,ierr)
       call mpi_allreduce(MPI_IN_PLACE,rho(ML_0_RHO,s),m,mpi_real8,mpi_sum,comm_bzsm,ierr)
    end do
    ! The following assumes all 'MS_1-MS_0+1' are the same
    m=m*(MS_1_RHO-MS_0_RHO+1)
    call mpi_allgather(rho(ML_0_RHO,MS_0_RHO),m,mpi_real8,rho,m,mpi_real8 &
         ,comm_spin,ierr)
  END SUBROUTINE reduce_and_gather


  SUBROUTINE calc_spin_density( Ntot )
    implicit none
    real(8),intent(OUT) :: Ntot(:)
    integer :: ierr
    real(8) :: tmp(4)
    tmp(:) = 0.0d0
    tmp(1) = sum( rho(:,1) )*dV_RHO
    if ( MS_RHO == 2 ) then
       tmp(2) = sum( rho(:,MS_RHO) )*dV_RHO
       tmp(3) = sum(     rho(:,1)-rho(:,MS_RHO)  )*dV_RHO
       tmp(4) = sum( abs(rho(:,1)-rho(:,MS_RHO)) )*dV_RHO
    end if
    call mpi_allreduce( tmp,Ntot,4,mpi_real8,mpi_sum,comm_grid,ierr)
  END SUBROUTINE calc_spin_density


  SUBROUTINE writeDensity( iter )
    implicit none
    integer,intent(IN) :: iter
    integer :: s,i
    call reduce_and_gather
    do s=MS_0_RHO,MS_1_RHO
       do i=ML_0_RHO,ML_1_RHO
          write(5000+iter,'(g20.12)') rho(i,s)
       end do
    end do
  END SUBROUTINE writeDensity


END MODULE density_module
