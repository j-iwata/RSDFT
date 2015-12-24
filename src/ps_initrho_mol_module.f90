MODULE ps_initrho_mol_module

  use pseudopot_module
  use atom_module
  use array_bound_module
  use density_module
  use rgrid_mol_module
  use polint_module

  PRIVATE
  PUBLIC :: init_ps_initrho_mol,construct_ps_initrho_mol

  logical :: flag_initrho_0
  logical,allocatable :: flag_initrho(:)

CONTAINS


  SUBROUTINE init_ps_initrho_mol
    implicit none
    integer :: ir,ik,MKI
    real(8) :: pi4

    MKI=Nelement
    pi4=4.d0*acos(-1.d0)

    allocate( flag_initrho(MKI) )
    flag_initrho(:)=.false.
    flag_initrho_0 =.false.
    if ( allocated(cdd) ) then
       do ik=1,MKI
          if ( all(cdd(:,ik)==0.d0) ) cycle
          flag_initrho(ik) = .true.
          flag_initrho_0   = .true.
       end do
    end if

    if ( .not.flag_initrho_0 ) return

    do ik=1,MKI
       do ir=2,Mr(ik)
          cdd(ir,ik)=cdd(ir,ik)/(pi4*rad(ir,ik)**2)
       end do
       if ( rad(1,ik) == 0.d0 ) then
          cdd(1,ik)=cdd(2,ik)
       else
          cdd(1,ik)=cdd(1,ik)/(pi4*rad(1,ik)**2)
       end if
    end do
  END SUBROUTINE init_ps_initrho_mol


  SUBROUTINE construct_ps_initrho_mol
    implicit none
    integer :: i,ir,ir0,ik,MKI,Mr0,Mr1,nr,a,M_irad,m,mm,m1,m2
    integer,allocatable ::NRc(:),irad(:,:)
    real(8) :: c,x,y,z,r,v,v0,err,err0,maxerr

    if ( .not.flag_initrho_0 ) return

    if ( .not. allocated(rho) ) then
       allocate( rho(ML_0:ML_1,MSP) )
       rho=0.d0
    end if

    MKI=Nelement

    allocate( NRc(MKI) )

    do ik=1,MKI
       Mr1=Mr(ik)
       Mr0=maxloc( cdd(1:Mr1,ik) ,1 )
       do ir=Mr0,Mr1
          if ( cdd(ir,ik) < 1.d-9 ) then
             NRc(ik)=ir
             exit
          end if
       end do
    end do ! ik

    allocate( irad(0:3000,MKI) ) ; irad=0

    M_irad=0
    do ik=1,MKI
       nr=min( 3000, NRc(ik) )
       m=0
       irad(0,ik)=1
       do ir=1,nr
          m=int(100.d0*rad(ir,ik))+1
          irad( m,ik )=ir
       end do
       ir=irad(0,ik)
       do i=1,m
          if ( irad(i,ik)==0 ) then
             irad(i,ik)=ir
             cycle
          end if
          ir=irad(i,ik)
       end do
       irad(m+1:,ik)=ir
       M_irad=max(M_irad,m)
    end do

    rho(:,:)=0.d0

    maxerr=0.d0
    do a=1,Natom
       ik = ki_atom(a)
       nr = NRc(ik)
       do i=ML_0,ML_1
          x=LL(1,i)*Hsize-aa_atom(1,a)
          y=LL(2,i)*Hsize-aa_atom(2,a)
          z=LL(3,i)*Hsize-aa_atom(3,a)
          r=sqrt(x*x+y*y+z*z)
          if ( r <= rad(nr,ik) ) then
             ir0=irad( int(100.d0*r),ik )
             do ir=ir0,nr
                if ( r < rad(ir,ik) ) exit
             end do
             err0=1.d10
             v0=0.d0
             do mm=1,20
                m1=max(1,ir-mm)
                m2=min(ir+mm,nr)
                call polint(rad(m1,ik),cdd(m1,ik),m2-m1+1,r,v,err)
                if ( abs(err) < err0 ) then
                   v0=v
                   err0=abs(err)
                   if ( err0 < 1.d-9 ) exit
                end if
             end do
             maxerr=max(maxerr,err0)
             rho(i,1)=rho(i,1)+v0
          end if
       end do ! i
    end do ! a

    deallocate( irad )
    deallocate( NRc )

    if ( MSP > 1 ) then
       c=1.d0/MSP
       rho(:,1)=c*rho(:,1)
       do i=2,MSP
          rho(:,i)=rho(:,1)
       end do
    end if

  END SUBROUTINE construct_ps_initrho_mol


END MODULE ps_initrho_mol_module
