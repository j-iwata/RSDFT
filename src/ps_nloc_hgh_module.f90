MODULE ps_nloc_hgh_module

  use aa_module
  use rgrid_module
  use atom_module
  use pseudopot_module
  use ylm_module

  implicit none

  PRIVATE
  PUBLIC :: prep_ps_nloc_hgh, init_force_ps_nloc_hgh, init_ps_nloc_hgh

  logical :: flag_init = .false.
  real(8),allocatable :: Rps0(:,:)

CONTAINS


  SUBROUTINE prep_ps_nloc_hgh(natm,n_0,L_0,MMJJ_0,M_grid_ion,map_grid_ion &
                             ,icheck_tmp3,JJ_tmp,MJJ_tmp,uV_tmp,nzlma,MMJJ)
    implicit none
    integer,intent(IN) :: natm,n_0,L_0,MMJJ_0
    integer,intent(IN) :: M_grid_ion
    integer,intent(IN) :: map_grid_ion(3,M_grid_ion)
    integer,intent(INOUT) :: icheck_tmp3(natm,n_0,2*L_0+1)
    integer,intent(INOUT) :: JJ_tmp(6,MMJJ_0,n_0,natm)
    integer,intent(INOUT) :: MJJ_tmp(n_0,natm)
    real(8),intent(INOUT) :: uV_tmp(MMJJ_0,n_0,natm)
    integer,intent(OUT) :: nzlma,MMJJ
    integer :: i,j,i1,i2,i3,id1,id2,id3,iorb
    integer :: m,n,L,ic1,ic2,ic3,ik,a,ML1,ML2,ML3
    integer :: i1_0,i2_0,i3_0,k1,k2,k3,lma,lma0
    real(8) :: c1,c2,c3,pi,gamma,Rx,Ry,Rz,v0,r2,r,x,y,z,d1,d2,d3
    real(8) :: Rps2,const1,const2

    pi                 = acos(-1.d0)
    ML1                = Ngrid(1)
    ML2                = Ngrid(2)
    ML3                = Ngrid(3)
    c1                 = 1.d0/Ngrid(1)
    c2                 = 1.d0/Ngrid(2)
    c3                 = 1.d0/Ngrid(3)
    icheck_tmp3(:,:,:) = 0
    MMJJ               = 0
    nzlma              = 0
    lma                = 0
    lma0               = 0

    do a=1,Natom

       Rx = aa(1,1)*aa_atom(1,a)+aa(1,2)*aa_atom(2,a)+aa(1,3)*aa_atom(3,a)
       Ry = aa(2,1)*aa_atom(1,a)+aa(2,2)*aa_atom(2,a)+aa(2,3)*aa_atom(3,a)
       Rz = aa(3,1)*aa_atom(1,a)+aa(3,2)*aa_atom(2,a)+aa(3,3)*aa_atom(3,a)

       ic1 = nint( aa_atom(1,a)*Ngrid(1) )
       ic2 = nint( aa_atom(2,a)*Ngrid(2) )
       ic3 = nint( aa_atom(3,a)*Ngrid(3) )

       ik = ki_atom(a)

       do iorb=1,norb(ik)

          Rps2 = Rps(iorb,ik)**2
          L    = lo(iorb,ik)
          n    = no(iorb,ik)

          gamma=sqrt(Pi)
          do i=1,L+2*n-1
             gamma=gamma*(i-0.5d0)
          end do
          const1 = sqrt(2.d0)/(Rps0(iorb,ik)**(L+2*n-0.5d0)*sqrt(gamma))
          const2 = 0.5d0/(Rps0(iorb,ik)*Rps0(iorb,ik))

          j=0
          do i=1,M_grid_ion

             i1 = map_grid_ion(1,i)
             i2 = map_grid_ion(2,i)
             i3 = map_grid_ion(3,i)

             id1 = ic1 + i1
             id2 = ic2 + i2
             id3 = ic3 + i3

             k1=id1/ML1 ; if ( id1<0 ) k1=(id1+1)/ML1-1 ; i1_0=id1-k1*ML1
             k2=id2/ML2 ; if ( id2<0 ) k2=(id2+1)/ML2-1 ; i2_0=id2-k2*ML2
             k3=id3/ML3 ; if ( id3<0 ) k3=(id3+1)/ML3-1 ; i3_0=id3-k3*ML3

             if ( Igrid(1,1) <= i1_0 .and. i1_0 <= Igrid(2,1) .and. &
                  Igrid(1,2) <= i2_0 .and. i2_0 <= Igrid(2,2) .and. &
                  Igrid(1,3) <= i3_0 .and. i3_0 <= Igrid(2,3) ) then

                d1 = id1*c1
                d2 = id2*c2
                d3 = id3*c3

                x  = aa(1,1)*d1+aa(1,2)*d2+aa(1,3)*d3-Rx
                y  = aa(2,1)*d1+aa(2,2)*d2+aa(2,3)*d3-Ry
                z  = aa(3,1)*d1+aa(3,2)*d2+aa(3,3)*d3-Rz
                r2 = x*x+y*y+z*z

                if ( r2 > Rps2 + 1.d-10 ) cycle

                r  = sqrt(r2)
                v0 = const1*r**(L+2*n-2)*exp(-r2*const2)

                nzlma = lma+2*L+1
                j=j+1
                JJ_tmp(1,j,iorb,a) = i1_0
                JJ_tmp(2,j,iorb,a) = i2_0
                JJ_tmp(3,j,iorb,a) = i3_0
                JJ_tmp(4,j,iorb,a) = k1
                JJ_tmp(5,j,iorb,a) = k2
                JJ_tmp(6,j,iorb,a) = k3
                uV_tmp(j,iorb,a)   = v0

             end if

          end do ! i ( 1 - M_grid_ion )

          MJJ_tmp(iorb,a)=j

       end do ! iorb
    end do ! a

  END SUBROUTINE prep_ps_nloc_hgh


  SUBROUTINE init_ps_nloc_hgh( disp_switch )
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
             if ( disp_switch ) write(*,*) ielm,iorb,n,L,Rps(iorb,ielm)
             exit
          end if
       end do
    end do ! iorb
    end do ! ielm

  END SUBROUTINE init_ps_nloc_hgh


  SUBROUTINE init_force_ps_nloc_hgh &
       (MMJJ,nzlma,amap,lmap,mmap,iorbmap,MJJ_MAP,JJ_MAP,Y1,Y2,Y3,duVdR)
    implicit none
    integer,intent(IN) :: MMJJ,nzlma
    integer,intent(IN) :: amap(nzlma),lmap(nzlma),mmap(nzlma),iorbmap(nzlma)
    integer,intent(IN) :: MJJ_MAP(nzlma),JJ_MAP(6,MMJJ,nzlma)
    real(8),intent(IN) :: Y1(0:3,-3:3,0:4,-4:4)
    real(8),intent(IN) :: Y2(0:3,-3:3,0:4,-4:4)
    real(8),intent(IN) :: Y3(0:3,-3:3,0:4,-4:4)
    real(8),intent(OUT) :: duVdR(3,MMJJ,nzlma)
    integer :: lma
    integer :: a,L,m,iorb,ik,n,j,L1,L1z
    real(8) :: c1,c2,c3,pi,cnst,cnst1,cnst2,cnst3,tmp0,tmp1,r2,e,rp
    real(8) :: gamma,Rx,Ry,Rz,x,y,z,r,d1,d2,d3,v0,v1,yy1,yy2,yy3

!$OMP workshare
    duVdR=0.d0
!$OMP end workshare

!$OMP single
    c1    = 1.d0/Ngrid(1)
    c2    = 1.d0/Ngrid(2)
    c3    = 1.d0/Ngrid(3)
    pi    = acos(-1.d0)
    cnst  = sqrt(4.d0*pi/3.d0)
!$OMP end single

!$OMP do private( a,L,m,n,iorb,ik,Rx,Ry,Rz,d1,d2,d3,x,y,z,r,v0,v1,gamma,r2 &
!$OMP            ,yy1,yy2,yy3,tmp0,tmp1,lma,j,L1,L1z,cnst1,cnst2,cnst3,e,rp )
    do lma=1,nzlma

       a    = amap(lma) ; if ( a <= 0 ) cycle
       L    = lmap(lma)
       m    = mmap(lma)
       iorb = iorbmap(lma)
       ik   = ki_atom(a)
       n    = no(iorb,ik)
       gamma=sqrt(pi)
       do j=1,L+2*n-1
          gamma=gamma*(j-0.5d0)
       end do
       cnst1=cnst*sqrt(2.d0)/(Rps0(iorb,ik)**(L+2*n-0.5d0)*sqrt(gamma))
       cnst2=0.5d0/(Rps0(iorb,ik)*Rps0(iorb,ik))
       cnst3=2.d0*cnst2
       Rx=aa(1,1)*aa_atom(1,a)+aa(1,2)*aa_atom(2,a)+aa(1,3)*aa_atom(3,a)
       Ry=aa(2,1)*aa_atom(1,a)+aa(2,2)*aa_atom(2,a)+aa(2,3)*aa_atom(3,a)
       Rz=aa(3,1)*aa_atom(1,a)+aa(3,2)*aa_atom(2,a)+aa(3,3)*aa_atom(3,a)

       if ( L == 0 .and. n == 1 ) then

          cnst1=cnst1*cnst3
          do j=1,MJJ_MAP(lma)
             d1 = c1*JJ_MAP(1,j,lma)+JJ_MAP(4,j,lma)
             d2 = c2*JJ_MAP(2,j,lma)+JJ_MAP(5,j,lma)
             d3 = c3*JJ_MAP(3,j,lma)+JJ_MAP(6,j,lma)
             x  = aa(1,1)*d1+aa(1,2)*d2+aa(1,3)*d3-Rx
             y  = aa(2,1)*d1+aa(2,2)*d2+aa(2,3)*d3-Ry
             z  = aa(3,1)*d1+aa(3,2)*d2+aa(3,3)*d3-Rz
             r2 = x*x+y*y+z*z
             r  = sqrt(r2)
             e  = exp(-r2*cnst2)
             v1 =-cnst1*r*e
             yy1=0.d0
             yy2=0.d0
             yy3=0.d0
             do L1z=-1,1
                tmp1 = v1*Ylm(x,y,z,1,L1z)
                yy1=yy1+tmp1*Y1(L,m,1,L1z)
                yy2=yy2+tmp1*Y2(L,m,1,L1z)
                yy3=yy3-tmp1*Y3(L,m,1,L1z)
             end do
             duVdR(1,j,lma)=yy1
             duVdR(2,j,lma)=yy2
             duVdR(3,j,lma)=yy3
          end do ! j

       else !( L == 0 .and. n == 1 )

          do j=1,MJJ_MAP(lma)
             d1 = c1*JJ_MAP(1,j,lma)+JJ_MAP(4,j,lma)
             d2 = c2*JJ_MAP(2,j,lma)+JJ_MAP(5,j,lma)
             d3 = c3*JJ_MAP(3,j,lma)+JJ_MAP(6,j,lma)
             x  = aa(1,1)*d1+aa(1,2)*d2+aa(1,3)*d3-Rx
             y  = aa(2,1)*d1+aa(2,2)*d2+aa(2,3)*d3-Ry
             z  = aa(3,1)*d1+aa(3,2)*d2+aa(3,3)*d3-Rz
             r2 = x*x+y*y+z*z
             r  = sqrt(r2)
             e  = exp(-r2*cnst2)
             rp = r**(L+2*n-3)
             v0 = rp*e
             v1 = (L+2*n-2-cnst3*r2)*rp*e
             yy1=0.d0
             yy2=0.d0
             yy3=0.d0
             do L1=abs(L-1),L+1
                tmp0 = cnst1*( v1 + 0.5d0*(2+L*(L+1)-L1*(L1+1))*v0 )
                do L1z=-L1,L1
                   tmp1=tmp0*Ylm(x,y,z,L1,L1z)
                   yy1=yy1+tmp1*Y1(L,m,L1,L1z)
                   yy2=yy2+tmp1*Y2(L,m,L1,L1z)
                   yy3=yy3-tmp1*Y3(L,m,L1,L1z)
                end do
             end do ! L1
             duVdR(1,j,lma)=yy1
             duVdR(2,j,lma)=yy2
             duVdR(3,j,lma)=yy3
          end do ! j

       end if !( L == 0 .and. n == 1 )

    end do ! lma
!$OMP end do

  END SUBROUTINE init_force_ps_nloc_hgh


END MODULE ps_nloc_hgh_module
