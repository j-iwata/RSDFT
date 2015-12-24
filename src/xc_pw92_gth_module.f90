MODULE xc_pw92_gth_module

  implicit none

  PRIVATE
  PUBLIC :: calc_pw92_gth

  real(8),parameter :: a0=0.4581652932831429d0 ,da0=0.119086804055547d0
  real(8),parameter :: a1=2.217058676663745d0  ,da1=0.6157402568883345d0
  real(8),parameter :: a2=0.7405551735357053d0 ,da2=0.1574201515892867d0
  real(8),parameter :: a3=0.01968227878617998d0,da3=0.003532336663397157d0
  real(8),parameter :: b1=1.0d0                ,db1=0.0d0
  real(8),parameter :: b2=4.504130959426697d0  ,db2=0.2673612973836267d0
  real(8),parameter :: b3=1.110667363742916d0  ,db3=0.2052004607777787d0
  real(8),parameter :: b4=0.02359291751427506d0,db4=0.004200005045691381d0

  include 'mpif.h'

CONTAINS

  SUBROUTINE calc_pw92_gth(n1,n2,nspin,s1,s2,rho,Exc,Vxc,dV,comm_grid)

    integer,intent(IN)  :: n1,n2,nspin,s1,s2,comm_grid
    real(8),intent(IN)  :: rho(n1:n2,nspin),dV
    real(8),intent(OUT) :: Exc,Vxc(n1:n2,s1:s2)
    real(8) :: trho,zeta,d0,d1,d2,d3,d4,d5,fx,rs,epsilon_xc
    real(8) :: a0z,a1z,a2z,a3z,b1z,b2z,b3z,b4z,factor(2)
    real(8) :: de_dr,dr_dn,de_df,df_dz
    integer :: i,ierr

    Exc = 0.d0
    Vxc = 0.d0

    d0 = 0.5d0*nspin
    d1 = 4.d0/3.d0
    d2 = 1.d0/( 2.d0*(2.d0**(1.d0/3.d0)-1.d0) )
    d3 = 3.d0/(4.d0*acos(-1.d0))
    d4 = 1.d0/3.d0
    d5 = 1.d0/(4.d0*acos(-1.d0))

    factor(1) = 2.d0
    factor(2) =-2.d0

    do i=n1,n2

       trho = d0*( rho(i,1)+rho(i,nspin) )

       if ( trho <= 0.d0 ) cycle

       zeta = ( rho(i,1)-rho(i,nspin) )/trho

       fx = ( abs(1.d0+zeta)**d1 + abs(1.d0-zeta)**d1 - 2.d0 )*d2

       rs = (d3/trho)**d4

       a0z = a0 + da0*fx
       a1z = a1 + da1*fx
       a2z = a2 + da2*fx
       a3z = a3 + da3*fx

       b1z = b1 + db1*fx
       b2z = b2 + db2*fx
       b3z = b3 + db3*fx
       b4z = b4 + db4*fx

       epsilon_xc = -( a0z + a1z*rs + a2z*rs**2 + a3z*rs**3 ) &
                    /( b1z*rs + b2z*rs**2 + b3z*rs**3 + b4z*rs**4 )

       Exc = Exc + trho*epsilon_xc

       de_df = ( -( da0 + da1*rs + da2*rs**2 + da3*rs**3 ) &
                 *( b1z*rs + b2z*rs**2 + b3z*rs**3 + b4z*rs**4 ) &
                 +( a0z + a1z*rs + a2z*rs**2 + a3z*rs**3 ) &
                 *( db1*rs + db2*rs**2 + db3*rs**3 + db4*rs**4 ) &
               )/( b1z*rs + b2z*rs**2 + b3z*rs**3 + b4z*rs**4 )**2

       df_dz = d1*d2*( abs(1.d0+zeta)**d4 - abs(1.d0-zeta)**d4 )

       de_dr = ( -( a1z + 2*a2z*rs + 3*a3z*rs**2 )*( b1z*rs + b2z*rs**2 + b3z*rs**3 + b4z*rs**4 ) &
                 +( a0z + a1z*rs + a2z*rs**2 + a3z*rs**3 ) &
                 *( b1z + 2*b2z*rs + 3*b3z*rs**2 + 4*b4z*rs**3 ) &
               )/( b1z*rs + b2z*rs**2 + b3z*rs**3 + b4z*rs**4 )**2

       dr_dn = -d5/(trho*rs)**2

       Vxc(i,s1) = epsilon_xc + trho*( de_dr*dr_dn + de_df*df_dz*factor(s1)*rho(i,s2)/trho**2 )
       Vxc(i,s2) = epsilon_xc + trho*( de_dr*dr_dn + de_df*df_dz*factor(s2)*rho(i,s1)/trho**2 )

    end do

    d0=Exc*dV
    call mpi_allreduce(d0,Exc,1,MPI_REAL8,MPI_SUM,comm_grid,ierr)

  END SUBROUTINE calc_pw92_gth

END MODULE xc_pw92_gth_module
