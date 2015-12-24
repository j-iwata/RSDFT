MODULE force_ewald_module

  use ewald_variables, only: eta, mg, mr, LG, LR, ipair, mpair
  use aa_module, only: aa, Va
  use bb_module, only: bb
  use atom_module, only: Natom, aa_atom, ki_atom
  use pseudopot_module, only: Zps
  use parallel_module

  implicit none

  PRIVATE
  PUBLIC :: calc_force_ewald

CONTAINS

  SUBROUTINE calc_force_ewald(MI,force3)
    implicit none
    integer,intent(IN) :: MI
    real(8),intent(OUT) :: force3(3,MI)
    real(8) :: sqeta,x,y,z,r,rr,ss,const,fewldg(3),fewldr(3)
    real(8) :: const1,const2,const3
    real(8) :: GR,c,c1,c2,sum0
    real(8) :: Vcell,pi,pi2
    real(8) :: sum1,sum2,sum3,sum_tmp(3)
    integer :: i,j,k,n,a,b,ab,ierr
    real(8) :: gg,gx,gy,gz
    real(8),allocatable :: work(:,:)

    force3(:,:) = 0.d0

    pi    = acos(-1.d0)
    pi2   = 2.d0*pi
    Vcell = abs(Va)
    sqeta = sqrt(eta)

    const1 = 2.d0*sqrt(eta/pi)
    const2 = 4.d0*pi/Vcell
    const3 = 1.d0/(4.d0*eta)
!$OMP parallel do private( sum_tmp,GR,sum0,gx,gy,gz,gg,c,c1 )
    do b=1,MI
       sum_tmp(:)=0.d0
       do i=1,mg
          if ( all(LG(1:3,i)==0) ) cycle
          sum0=0.d0
          do a=1,Natom
             GR=pi2*( LG(1,i)*(aa_atom(1,a)-aa_atom(1,b)) &
                     +LG(2,i)*(aa_atom(2,a)-aa_atom(2,b)) &
                     +LG(3,i)*(aa_atom(3,a)-aa_atom(3,b)) )
             sum0=sum0+Zps(ki_atom(a))*sin(GR)
          end do
          gx=bb(1,1)*LG(1,i)+bb(1,2)*LG(2,i)+bb(1,3)*LG(3,i)
          gy=bb(2,1)*LG(1,i)+bb(2,2)*LG(2,i)+bb(2,3)*LG(3,i)
          gz=bb(3,1)*LG(1,i)+bb(3,2)*LG(2,i)+bb(3,3)*LG(3,i)
          gg=gx*gx+gy*gy+gz*gz
          c=exp(-const3*gg)/gg
          sum_tmp(1)=sum_tmp(1)+gx*sum0*c
          sum_tmp(2)=sum_tmp(2)+gy*sum0*c
          sum_tmp(3)=sum_tmp(3)+gz*sum0*c
       end do
       c1=Zps(ki_atom(b))*const2
       force3(1,b) = force3(1,b) - sum_tmp(1)*c1
       force3(2,b) = force3(2,b) - sum_tmp(2)*c1
       force3(3,b) = force3(3,b) - sum_tmp(3)*c1
    end do ! b
!$OMP end parallel do

#ifdef TEST
    do b=1,MI
       do a=1,MI
          if ( a==b ) cycle
          sum1=0.d0
          sum2=0.d0
          sum3=0.d0
          do i=1,mr
             x=sum( aa(1,1:3)*(aa_atom(1:3,b)-aa_atom(1:3,a)+LR(1:3,i)) )
             y=sum( aa(2,1:3)*(aa_atom(1:3,b)-aa_atom(1:3,a)+LR(1:3,i)) )
             z=sum( aa(3,1:3)*(aa_atom(1:3,b)-aa_atom(1:3,a)+LR(1:3,i)) )
             rr=x*x+y*y+z*z
             r=sqrt(rr)
             c=-( erfc(sqeta*r) + const1*exp(-eta*rr)*r )/(r*rr)
             sum1=sum1+c*x
             sum2=sum2+c*y
             sum3=sum3+c*z
          end do
          c2=2.d0*Zps(ki_atom(a))*Zps(ki_atom(b))
          force3(1,b) = force3(1,b) - 0.5d0*sum1*c2
          force3(2,b) = force3(2,b) - 0.5d0*sum2*c2
          force3(3,b) = force3(3,b) - 0.5d0*sum3*c2
       end do
    end do ! b
#endif

    do ab=1,mpair
       a=ipair(1,ab)
       b=ipair(2,ab)
       if ( a == b ) cycle
       c2=2.d0*Zps(ki_atom(a))*Zps(ki_atom(b))
       sum1=0.d0
       sum2=0.d0
       sum3=0.d0
!$OMP parallel do private( x,y,z,rr,r,c ) reduction(+:sum1,sum2,sum3)
       do i=1,mr
          x=sum( aa(1,1:3)*(aa_atom(1:3,b)-aa_atom(1:3,a)+LR(1:3,i)) )
          y=sum( aa(2,1:3)*(aa_atom(1:3,b)-aa_atom(1:3,a)+LR(1:3,i)) )
          z=sum( aa(3,1:3)*(aa_atom(1:3,b)-aa_atom(1:3,a)+LR(1:3,i)) )
          rr=x*x+y*y+z*z
          r=sqrt(rr)
          c=-( erfc(sqeta*r) + const1*exp(-eta*rr)*r )/(r*rr)
          sum1=sum1+c*x
          sum2=sum2+c*y
          sum3=sum3+c*z
       end do
!$OMP end parallel do
       force3(1,b) = force3(1,b) - 0.5d0*sum1*c2
       force3(2,b) = force3(2,b) - 0.5d0*sum2*c2
       force3(3,b) = force3(3,b) - 0.5d0*sum3*c2
       sum1=0.d0
       sum2=0.d0
       sum3=0.d0
!$OMP parallel do private( x,y,z,rr,r,c ) reduction(+:sum1,sum2,sum3)
       do i=1,mr
          x=sum( aa(1,1:3)*(aa_atom(1:3,a)-aa_atom(1:3,b)+LR(1:3,i)) )
          y=sum( aa(2,1:3)*(aa_atom(1:3,a)-aa_atom(1:3,b)+LR(1:3,i)) )
          z=sum( aa(3,1:3)*(aa_atom(1:3,a)-aa_atom(1:3,b)+LR(1:3,i)) )
          rr=x*x+y*y+z*z
          r=sqrt(rr)
          c=-( erfc(sqeta*r) + const1*exp(-eta*rr)*r )/(r*rr)
          sum1=sum1+c*x
          sum2=sum2+c*y
          sum3=sum3+c*z
       end do
!$OMP end parallel do
       force3(1,a) = force3(1,a) - 0.5d0*sum1*c2
       force3(2,a) = force3(2,a) - 0.5d0*sum2*c2
       force3(3,a) = force3(3,a) - 0.5d0*sum3*c2
    end do

    call mpi_allreduce(MPI_IN_PLACE,force3,3*MI,mpi_real8,mpi_sum,mpi_comm_world,ierr)

  END SUBROUTINE calc_force_ewald

END MODULE force_ewald_module
