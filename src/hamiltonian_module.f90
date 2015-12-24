MODULE hamiltonian_module

  use parallel_module, only: myrank
  use kinetic_module
  use localpot_module
  use nonlocal_module
  use fock_module
  use watch_module

  implicit none

  PRIVATE
  PUBLIC :: hamiltonian,op_kinetic,op_localpot,op_nonlocal &
           ,ctt_hamil,ett_hamil

  real(8) :: ctt_hamil(4),ett_hamil(4)

CONTAINS

  SUBROUTINE hamiltonian(k,s,tpsi,htpsi,n1,n2,ib1,ib2)
    implicit none
    integer,intent(IN) :: k,s,n1,n2,ib1,ib2
#ifdef _DRSDFT_
    real(8),intent(IN)  :: tpsi(n1:n2,ib1:ib2)
    real(8),intent(OUT) :: htpsi(n1:n2,ib1:ib2)
#else
    complex(8),intent(IN)  :: tpsi(n1:n2,ib1:ib2)
    complex(8),intent(OUT) :: htpsi(n1:n2,ib1:ib2)
#endif
    real(8) :: ttmp(2)

!$OMP parallel

!$OMP workshare
    htpsi=(0.d0,0.d0)
!$OMP end workshare

    call watchb_omp( ttmp )

! --- Kinetic energy ---

    call op_kinetic(k,tpsi,htpsi,n1,n2,ib1,ib2)

!$OMP barrier

    call watchb_omp( ttmp, time_hmlt(1,1) )

! --- local potential ---

    call op_localpot(s,n2-n1+1,ib2-ib1+1,tpsi,htpsi)

!$OMP barrier

    call watchb_omp( ttmp, time_hmlt(1,2) )

! --- nonlocal potential ---

    call op_nonlocal(k,s,tpsi,htpsi,n1,n2,ib1,ib2)

!$OMP barrier

    call watchb_omp( ttmp, time_hmlt(1,3) )

!$OMP end parallel

    call watchb( ttmp )

    call op_fock(k,s,n1,n2,ib1,ib2,tpsi,htpsi)

    call watchb( ttmp, time_hmlt(1,4) )

  END SUBROUTINE hamiltonian

END MODULE hamiltonian_module
