MODULE ps_q_init_module

  use var_ps_member
  use var_ps_member_g
  use parallel_module, only: myrank
  use polint_module
  use maskf_module
  use Filtering, only: opFiltering

  implicit none
  
  PRIVATE
  PUBLIC :: initKtoKPSQ_ps_q_init
  PUBLIC :: ps_q_init
  PUBLIC :: ps_q_init_derivative

CONTAINS

!-----------------------------------------------

  SUBROUTINE initKtoKPSQ_ps_q_init

    implicit none
    integer :: ik,l1,l2,i1,i2,m1,m2
    integer :: k1,k2,k3,nr1,nr2
    integer :: mm1,mm2,mmin,mmax
    logical :: disp_switch_local
    logical :: registered_k2(1:max_k2)

    call write_border( 0," initKtoKPSQ_pq_q_init" )

    disp_switch_local=(myrank==0)

    call allocateKtoK( k1max, max_k2, Nelement_, max_Rref, max_Lref )

    do ik=1,Nelement_

       k1=0
       k2=0
       nr1=0
       do l1=1,nlf(ik)
       do i1=1,nrf(l1,ik)
          nr1=nr1+1
          nr2=0
          do l2=1,nlf(ik)
          do i2=1,nrf(l2,ik)
             nr2=nr2+1
             if ( l2 > l1 ) cycle
             if ( .not.( l1 == l2 .and. i2 > i1 ) ) then
                k2=k2+1
                do m1=1,2*l1-1
                do m2=1,2*l2-1
                   if ( .not.( l1 == l2 .and. i1 == i2 .and. m2 > m1 ) ) then
                      k1=k1+1
                      mm1=(l1-1)**2+m1
                      mm2=(l2-1)**2+m2
                      mmin=min(mm1,mm2)
                      mmax=max(mm1,mm2)
                      k3=(mmax-1)*mmax/2 + mmin
                      k1_to_k2(k1,ik)=k2
                      k1_to_k3(k1,ik)=k3
                      k1_to_iorb(1,k1,ik)=nr1
                      k1_to_iorb(2,k1,ik)=nr2
                      k1_to_l(1,k1,ik)=l1-1
                      k1_to_l(2,k1,ik)=l2-1
                      k1_to_m(1,k1,ik)=m1-l1
                      k1_to_m(2,k1,ik)=m2-l2
                   end if
                end do ! m2
                end do ! m1
             end if
          end do ! i2
          end do ! l2
       end do ! i1
       end do ! l1
       N_k1(ik)=k1
    end do ! ik

    N_k2(:)=0
    k2_to_iorb(:,:,:)=0

    do ik=1,Nelement_
       registered_k2(:)=.false.
       do k1=1,N_k1(ik)
          k2=k1_to_k2(k1,ik)
          if ( registered_k2(k2) ) cycle
          registered_k2(k2)=.true.
          N_k2(ik)=N_k2(ik)+1
          k2_to_iorb(1,k2,ik)=k1_to_iorb(1,k1,ik)
          k2_to_iorb(2,k2,ik)=k1_to_iorb(2,k1,ik)
       end do
    end do

    return
  END SUBROUTINE initKtoKPSQ_ps_q_init

!-----------------------------------------------

  SUBROUTINE ps_q_init( qcut, rcfac, qcfac, etafac )

    implicit none
    real(8),intent(IN) :: qcut
    real(8),intent(IN) :: rcfac,qcfac,etafac
    integer :: MMr,NRc,Q_NRc
    integer :: ik,k2,i,m,m0,m1,m2,m3,ll3,L,l1,l2,i1,i2
    integer :: iorb1,iorb2
    real(8),parameter :: dr=2.d-3, ep=1.d-14
    real(8) :: Rc,x,y,dy,y0,dy0
    real(8) :: maxerr
    real(8) :: qc,Q_rcfac,Q_etafac
    real(8),allocatable :: vrad(:),Q_wm(:,:,:)

    call write_border( 0, " ps_q_init" )

    Q_rcfac  = rcfac
    Q_etafac = etafac
    k2max    = max_k2

    qc = qcut*qcfac
    if ( qc <= 0.d0 ) qc=qcut

    call allocateQRps( k2max,Nelement_ )

! ---

    do ik=1,Nelement_

       MMr = Mr(ik)

       do k2=1,N_k2(ik)

          iorb1 = k2_to_iorb(1,k2,ik)
          iorb2 = k2_to_iorb(2,k2,ik)

          Rc  = max( Rps0(iorb1,ik), Rps0(iorb2,ik) )*Q_rcfac
          NRc = minloc( abs(rad(1:MMr,ik)-Rc), 1 )
          if ( rad(NRc,ik) < Rc ) NRc = NRc + 1
          if ( NRc > MMr ) then
             write(*,*) "NRc,MMr= ",NRc,MMr
             stop "rcfac is too large:stop@ps_q_init(1)"
          end if

          Q_NRps(k2,ik) = NRc
          Q_Rps(k2,ik)  = rad(NRc,ik)

       end do ! k2

    end do ! ik

! ---

    NRc = maxval( Q_NRps )
    allocate( Q_wm(NRc,k2max,Nelement_) ) ; Q_wm=0.0d0

    do ik=1,Nelement_

       do k2=1,N_k2(ik)

          NRc = Q_NRps(k2,ik)
          Rc  = Q_Rps(k2,ik)
          call makemaskf( Q_etafac )

          maxerr=0.d0
          do i=1,NRc
             x=rad(i,ik)/Rc
             if ( x<=dxm ) then
                y0=1.d0 ; dy0=0.d0
             else
                m0=int(x/dxm)
                dy0=1.d10
                do m=1,20
                   m1=max(m0-m,1) ; m2=min(m0+m,nmsk)
                   call polint(xm(m1),maskr(m1),m2-m1+1,x,y,dy)
                   if ( abs(dy)<dy0 ) then
                      y0=y ; dy0=abs(dy)
                   end if
                end do
             end if
             Q_wm(i,k2,ik)=y0
             maxerr=max(maxerr,dy0)
          end do

       end do ! k2

    end do ! ik

! ---

    do ik=1,Nelement_
       do k2=1,N_k2(ik)
          Q_NRps(k2,ik)=Q_Rps(k2,ik)/dr+1
          if ( (Q_NRps(k2,ik)-1)*dr < Q_Rps(k2,ik) ) then
             Q_NRps(k2,ik)=Q_NRps(k2,ik)+1
          endif
       end do
    end do

    MMr=max( maxval(Mr), maxval(Q_NRps) )

    call allocateRad1( MMr )

    do ik=1,Nelement_
       do i=1,MMr
          rad1(i,ik)=(i-1)*dr
       end do
    end do

    NRc   = maxval( NRps0 )
    Q_NRc = maxval( Q_NRps )
    Q_NRc = max( Q_NRc, NRc )
    call resize_qrL( Q_NRc )

! ---

    allocate( vrad(NRc) ) ; vrad=0.0d0

    do ik=1,Nelement_
       do k2=1,N_k2(ik)
          iorb1=k2_to_iorb(1,k2,ik)
          iorb2=k2_to_iorb(2,k2,ik)
          NRc=max( NRps0(iorb1,ik),NRps0(iorb2,ik) )
          do ll3=1,nl3v(k2,ik)
             L=l3v(ll3,k2,ik)-1
             vrad(1:NRc)=qrL(1:NRc,ll3,k2,ik)*rab(1:NRc,ik)/Q_wm(1:NRc,k2,ik)
             qrL(:,ll3,k2,ik)=0.0d0
             call opFiltering( qc,L,NRc,Q_NRps(k2,ik),rad(1,ik),rad1(1,ik) &
                              ,vrad,qrL(1,ll3,k2,ik) )
          end do ! ll3
       end do ! k2
    end do ! ik

    deallocate( vrad )

! ---

    do ik=1,Nelement_

       do k2=1,N_k2(ik)

          NRc=Q_NRps(k2,ik)
          Rc =Q_Rps(k2,ik)
          call makemaskf( Q_etafac )

          maxerr=0.d0
          do i=1,NRc
             x=(i-1)*dr/Rc
             if ( x<dxm ) then
                y0=1.d0 ; dy0=0.d0
             else
                m0=int(x/dxm)
                dy0=1.d10
                do m=1,20
                   m1=max(m0-m,1) ; m2=min(m0+m,nmsk)
                   call polint(xm(m1),maskr(m1),m2-m1+1,x,y,dy)
                   if ( abs(dy)<dy0 ) then
                      y0=y ; dy0=abs(dy)
                   end if
                end do
             end if
             if ( maxerr < dy0 ) maxerr=dy0
             do ll3=1,nl3v(k2,ik)
                qrL(i,ll3,k2,ik)=y0*qrL(i,ll3,k2,ik)
             end do
          end do ! i
       end do ! k2

    end do ! ik

    deallocate( Q_wm )

    return
  END SUBROUTINE ps_q_init


  SUBROUTINE resize_qrL( Q_NRc )
    implicit none
    integer,intent(IN) :: Q_NRc
    real(8),allocatable :: qrLtmp(:,:,:,:)
    integer :: m0,m1,m2,m3
    if ( size(qrL,1) < Q_NRc ) then
       m0=size(qrL,1)
       m1=size(qrL,2)
       m2=size(qrL,3)
       m3=size(qrL,4)
       allocate( qrLtmp(m0,m1,m2,m3) ) ; qrLtmp=0.0d0
       qrLtmp(:,:,:,:)=qrL(:,:,:,:)
       deallocate( qrL )
       allocate( qrL(Q_NRc,m1,m2,m3) ) ; qrL=0.0d0
       qrL(1:m0,1:m1,1:m2,1:m3)=qrLtmp(1:m0,1:m1,1:m2,1:m3)
       deallocate( qrLtmp )
    end if
  END SUBROUTINE resize_qrL


  SUBROUTINE ps_q_init_derivative
    implicit none
    integer :: ik,L,NRc,J,i,m,m1,m2,lm,ll3,ik1,ik2,cJ
    real(8) :: maxerr,y,dy,y0,dy0,r,const0,const1
    real(8),allocatable :: dwork(:,:,:)
    real(8),parameter :: sqrt_4pi_3=sqrt(4.d0*acos(-1.d0)/3.d0)

    call write_border( 0, " ps_q_init_derivative(start)" )

    maxcJ=0
    do ik=1,Nelement_
       do ik1=1,N_k1(ik)
          ik2=k1_to_k2(ik1,ik)
          do ll3=1,nl3v(ik2,ik)
             L=l3v(ll3,ik2,ik)-1
             cJ=0
             do J=abs(L-1),L+1
                cJ=cJ+1
             end do
             maxcJ=max(cJ,maxcJ)
          end do
       end do
    end do

    NRc = size( qrL, 1 ) !max_psgrd
    l   = size( qrL, 2 ) !max_Lref
    ik2 = size( qrL, 3 ) !max_k2
    ik  = Nelement_

    allocate( dqrL(NRc,l,ik2,ik,maxcJ) ) ; dqrL=0.0d0

    allocate( dwork(NRc,l,ik2) ) ; dwork=0.0d0

    const0 = sqrt( 4.0d0*acos(-1.0d0)/3.0d0 )

    do ik=1,Nelement_

       do ik2=1,N_k2(ik)
          NRc=Q_NRps(ik2,ik)
          do ll3=1,nl3v(ik2,ik)
             L=l3v(ll3,ik2,ik)-1
             maxerr=0.0d0
             do i=1,NRc
                dy0=1.d10
                do m=1,20
                   m1=max(i-m,1) ; m2=min(i+m,NRc)
                   call dpolint( rad1(m1,ik),qrL(m1,ll3,ik2,ik) &
                                ,m2-m1+1,rad1(i,ik),y,dy )
                   if ( abs(dy)<dy0 ) then
                      y0=y ; dy0=abs(dy)
                   end if
                end do ! m
                maxerr=max(maxerr,dy0)
                dwork(i,ll3,ik2)=y0
             end do ! i
          enddo ! ll3
       end do ! ik2

       do ik2=1,N_k2(ik)
          NRc=Q_NRps(ik2,ik)
          do ll3=1,nl3v(ik2,ik)
             L=l3v(ll3,ik2,ik)-1
             cJ=0
             do J=abs(L-1),L+1
                cJ=cJ+1
                const1=0.5d0*(2.d0+L*(L+1)-J*(J+1))
                do i=1,NRc
                   r = rad1(i,ik)
                   dqrL(i,ll3,ik2,ik,cJ) = const0*r*r*( r*dwork(i,ll3,ik2) &
                        + const1*qrL(i,ll3,ik2,ik) )
                end do
             end do ! J
          end do ! ll3
       enddo ! ik2

    end do ! ik

    deallocate( dwork )

    call write_border( 0, " ps_q_init_derivative(end)" )

  END SUBROUTINE ps_q_init_derivative


END MODULE ps_q_init_module
