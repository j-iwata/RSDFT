MODULE ps_nloc_mol_gth_module

!  use aa_module
!  use rgrid_module
!  use atom_module
  use pseudopot_module

  implicit none

  PRIVATE
  PUBLIC :: init_ps_nloc_mol_gth

  logical :: flag_init = .false.
  real(8),allocatable :: Rps0(:,:)

CONTAINS


  SUBROUTINE init_ps_nloc_mol_gth( disp_switch )
    implicit none
    logical,intent(IN) :: disp_switch
    real(8) :: rmax,dr,pi,ep,const1,const2,gamma,r,v
    integer :: MMr,iorb,m,n,ielm,L,i

    if ( flag_init ) return
    flag_init = .true.

    rmax = 30.d0
    MMr  = 3000
    dr   = rmax/MMr
    pi   = acos(-1.0d0)
    ep   = 1.d-10
    m=size( Rps,1 )
    n=size( Rps,2 )
    allocate( Rps0(m,n) )
    Rps0(:,:) = 0.0d0
    Rps0(:,:) = Rps(:,:)
    Rps(:,:) = 0.0d0

    do ielm=1,n
    do iorb=1,norb(ielm)
       n=no(iorb,ielm)
       L=lo(iorb,ielm)
       gamma=sqrt(pi)
       do i=1,L+2*n-1
          gamma=gamma*(i-0.5d0)
       end do
       const1=sqrt(2.0d0)/(Rps0(iorb,ielm)**(L+2*n-0.5d0)*sqrt(gamma))
       const2=0.5d0/(Rps0(iorb,ielm)**2)
       do i=1,MMr
          r=i*dr
          v=const1*r**(L+2*n-2)*exp(-r*r*const2)
          if ( abs(v) < ep ) then
             Rps(iorb,ielm)=max( Rps(iorb,ielm),r )
             exit
          end if
       end do
    end do ! iorb
    end do ! ielm

  END SUBROUTINE init_ps_nloc_mol_gth


END MODULE ps_nloc_mol_gth_module
