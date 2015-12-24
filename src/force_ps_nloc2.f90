MODULE force_ps_nloc2

!----------------------------------------------------------------------------
! this module calculates force for ultrasoft PS
! this module use subroutines from force_sub_module.f90
!----------------------------------------------------------------------------

  use atom_module, only: ki_atom,aa_atom,Natom
  use aa_module, only: aa
  use var_ps_member
  use var_ps_member_g
  use var_para_ps_nloc_g, only: MAXMJJ_MAP_Q, MJJ_MAP_Q, JJ_MAP_Q
  use ps_nloc2_variables, only: amap,iorbmap,mmap,lmap,nzlma,Mlma &
                               ,sendmap,sbufnl,recvmap,rbufnl &
                               ,lma_nsend,num_2_rank, nrlma_xyz&
                               ,MJJ,MJJ_MAP,JJ_MAP,JJP,MMJJ,uVk
  use rgrid_module
  use array_bound_module
  use parallel_module, only: MB_d, disp_switch_parallel, comm_grid
  use localpot_module, only: Vloc
  use para_rgrid_comm, only: do3StepComm_F
  use force_sub_module
  use ylm_module
  use ps_q_init_module, only: ps_q_init_derivative
  use bz_module
  use wf_module
  use watch_module
  use force_sub_sub_module

  implicit none

  PRIVATE
  PUBLIC :: calcForcePSnonLoc2

  include 'mpif.h'

  type(gaunt) :: GCL1
  type(gaunt_ll) :: GCLL

CONTAINS

  SUBROUTINE calcForcePSnonLoc2(MI,force2)
    implicit none
    integer,intent(IN) :: MI
    real(8),intent(OUT) :: force2(3,MI)
    integer :: ib1,ib2
    integer :: i,j,k,s,n,ir,L1,L1z,NRc,irank,jrank
    integer :: nreq,max_nreq
    integer :: a0,lm0,lm1,lma
    integer :: a,ik,m,L
    integer :: ierr,ir0
    integer,allocatable :: ireq(:),istatus(:,:),ilm1(:,:,:)
    real(8) :: a1,a2,a3
    real(8) :: kr,c
    real(8) :: tmp
    real(8) :: ctt(0:4),ett(0:4)
    real(8) :: yy1,yy2,yy3
    real(8),allocatable :: work2(:,:),duVdR(:,:,:)
    real(8) :: Dij_f(MB_0:MB_1),const_f(MB_0:MB_1)
#ifdef _DRSDFT_
    real(8) :: ztmp
    real(8),allocatable :: wtmp5(:,:,:,:,:)
    real(8),allocatable :: uVunk_tmp(:,:,:,:)
    real(8),parameter :: zero=0.0d0
#else
    complex(8) :: ztmp
    complex(8),allocatable :: wtmp5(:,:,:,:,:)
    complex(8),allocatable :: uVunk_tmp(:,:,:,:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
#endif
    integer :: i0,iorb0
    real(8) :: forceQ(3,MI)

    real(8),parameter :: pi2=2.d0*acos(-1.d0)
! d1,d2 need to be removed from this module, define a different variable name
    real(8) :: d1,d2,d3
    integer :: i1,i2,i3

    L = max(  maxval(l3v), maxval(lo) )
    call construct_gaunt_coef_L1( L, GCL1 )
    L  = maxval( lo )
    call construct_gaunt_coef_LL( L, L, GCLL )

    force2(:,:) = 0.0d0

    if ( Mlma <= 0 ) return

    maxerr=0.d0
    ctt(:)=0.d0
    ett(:)=0.d0

    call setLocalIndexForBoundary(Igrid)
    !OUT: a1b,a2b,a3b,ab1,ab2,ab3

    call setConstGridWidth(Ngrid)
    !OUT: ML1,ML2,ML3,c1,c2,c3

    if ( .not.allocated(ilm1) ) then
       L1=maxval(lo)+1
       n=maxval(norb)
       allocate( ilm1(0:L1,n,Nelement_) ) ; ilm1=0
       do ik=1,Nelement_
          lm1=0
          do iorb=1,norb(ik)
             L=lo(iorb,ik)
             do L1=abs(L-1),L+1
                lm1=lm1+1
                ilm1(L1,iorb,ik)=lm1
             end do
          end do
       end do
    end if

    allocate( wtmp5(0:3,nzlma,MB_0:MB_1,MBZ_0:MBZ_1,MSP_0:MSP_1) )
    allocate( duVdR(3,MMJJ,nzlma) )

    call setLocalAtoms(MI,aa_atom,Igrid)
    !OUT: isAtomInThisNode
!$OMP parallel

!$OMP workshare
    wtmp5=zero
    duVdR=0.0d0
!$OMP end workshare

!$OMP master
    call watch(ctt(0),ett(0))
!$OMP end master

!!$OMP do schedule(dynamic) firstprivate( maxerr ) &
!!$OMP    private( a,L,m,iorb,ik,Rx,Ry,Rz,NRc,d1,d2,d3,x,y,z,r  &
!!$OMP            ,ir,ir0,yy1,yy2,yy3,err0,err,tmp0,tmp1,m1,m2  &
!!$OMP            ,lma,j,L1,L1z,lm1,im )
    do lma=1,nzlma

       call getAtomInfoFrom_lma(lma)
       !OUT: iatom,iL,im,iorb,ikind,notAtom

       if ( noAtomHere ) cycle

       call getAtomPosition_real
       !OUT: Rx,Ry,Rz

       NRc=NRps(iorb,ikind)

!!$OMP parallel do firstprivate( maxerr ) &
!!$OMP             private( d1,d2,d3,x,y,z,r,ir0,yy1,yy2,yy3,lm1,err,err0 &
!!$OMP                     ,tmp0,tmp1,m1,m2,j,L1,im,L1z )
       do j=1,MJJ_MAP(lma)

          call getAtomCenteredPositionFrom_lma(lma,j)
          !OUT: x,y,z,r

          call get_nearest_index( rad1(:,ikind), r, ir )

          yy1=0.0d0
          yy2=0.0d0
          yy3=0.0d0
          do L1=abs(iL-1),iL+1
             lm1=ilm1(L1,iorb,ikind)
             if ( abs(x)>1.d-14 .or. abs(y)>1.d-14 .or. &
                  abs(z)>1.d-14 .or. L1==0 ) then
                call interpolate_dviod(lm1,ir,NRc) !OUT: tmp0
                do L1z=-L1,L1
                   tmp=tmp0*Ylm(x,y,z,L1,L1z)
                   yy1=yy1+tmp*GCL1%x(iL,im,L1,L1z)
                   yy2=yy2+tmp*GCL1%y(iL,im,L1,L1z)
                   yy3=yy3+tmp*GCL1%z(iL,im,L1,L1z)
                end do
             end if
          end do ! L1

          duVdR(1,j,lma)= yy1
          duVdR(2,j,lma)= yy2
          duVdR(3,j,lma)=-yy3

       end do ! j
!!$OMP end parallel do
    end do ! lma
!!$OMP end do

!$OMP master
    call watch(ctt(1),ett(1))
!$OMP end master

    do s=MSP_0,MSP_1
       do k=MBZ_0,MBZ_1
!$OMP do schedule(dynamic) private( c,i,d1,d2,d3,kr,ztmp,i1,i2,i3 )
          do n=MB_0,MB_1
             if ( occ(n,k,s) < 1.d-10 ) cycle
             do lma=1,nzlma
                if ( MJJ_MAP(lma) == MJJ(lma) ) then
#ifdef _DRSDFT_
                   do j=1,MJJ(lma)
                      i=JJP(j,lma)
                      wtmp5(0,lma,n,k,s)=wtmp5(0,lma,n,k,s) &
                           +uVk(j,lma,k)*unk(i,n,k,s)
                      wtmp5(1,lma,n,k,s)=wtmp5(1,lma,n,k,s) &
                           +duVdR(1,j,lma)*unk(i,n,k,s)
                      wtmp5(2,lma,n,k,s)=wtmp5(2,lma,n,k,s) &
                           +duVdR(2,j,lma)*unk(i,n,k,s)
                      wtmp5(3,lma,n,k,s)=wtmp5(3,lma,n,k,s) &
                           +duVdR(3,j,lma)*unk(i,n,k,s)
                   end do
#else
                   do j=1,MJJ(lma)
                      i=JJP(j,lma)
                      wtmp5(0,lma,n,k,s)=wtmp5(0,lma,n,k,s) &
                           +uVk(j,lma,k)*conjg(unk(i,n,k,s))
                   end do
                   do j=1,MJJ_MAP(lma)
                      i=JJP(j,lma)
                      d1=c1*JJ_MAP(1,j,lma)+JJ_MAP(4,j,lma)
                      d2=c2*JJ_MAP(2,j,lma)+JJ_MAP(5,j,lma)
                      d3=c3*JJ_MAP(3,j,lma)+JJ_MAP(6,j,lma)
                      kr=pi2*(kbb(1,k)*d1+kbb(2,k)*d2+kbb(3,k)*d3)
                      ztmp=dcmplx(cos(kr),sin(kr))*unk(i,n,k,s)
                      wtmp5(1,lma,n,k,s)=wtmp5(1,lma,n,k,s)+duVdR(1,j,lma)*ztmp
                      wtmp5(2,lma,n,k,s)=wtmp5(2,lma,n,k,s)+duVdR(2,j,lma)*ztmp
                      wtmp5(3,lma,n,k,s)=wtmp5(3,lma,n,k,s)+duVdR(3,j,lma)*ztmp
                   end do
#endif
                else ! --- MJJ(lma) /= MJJ_MAP(lma) ---
#ifdef _DRSDFT_
                   do j=1,MJJ(lma)
                      i=JJP(j,lma)
                      wtmp5(0,lma,n,k,s)=wtmp5(0,lma,n,k,s) &
                           +uVk(j,lma,k)*unk(i,n,k,s)
                   end do
                   do j=1,MJJ_MAP(lma)
                      i1=JJ_MAP(1,j,lma)
                      i2=JJ_MAP(2,j,lma)
                      i3=JJ_MAP(3,j,lma)
                      i = i1-a1b + (i2-a2b)*ab1 + (i3-a3b)*ab1*ab2 + ML_0
                      wtmp5(1,lma,n,k,s)=wtmp5(1,lma,n,k,s) &
                           +duVdR(1,j,lma)*unk(i,n,k,s)
                      wtmp5(2,lma,n,k,s)=wtmp5(2,lma,n,k,s) &
                           +duVdR(2,j,lma)*unk(i,n,k,s)
                      wtmp5(3,lma,n,k,s)=wtmp5(3,lma,n,k,s) &
                           +duVdR(3,j,lma)*unk(i,n,k,s)
                   end do
#else
                   do j=1,MJJ(lma)
                      i=JJP(j,lma)
                      wtmp5(0,lma,n,k,s)=wtmp5(0,lma,n,k,s) &
                           +uVk(j,lma,k)*conjg(unk(i,n,k,s))
                   end do
                   do j=1,MJJ_MAP(lma)
                      i1=JJ_MAP(1,j,lma)
                      i2=JJ_MAP(2,j,lma)
                      i3=JJ_MAP(3,j,lma)
                      i = i1-a1b + (i2-a2b)*ab1 + (i3-a3b)*ab1*ab2 + ML_0
                      d1=c1*i1+JJ_MAP(4,j,lma)
                      d2=c2*i2+JJ_MAP(5,j,lma)
                      d3=c3*i3+JJ_MAP(6,j,lma)
                      kr=pi2*(kbb(1,k)*d1+kbb(2,k)*d2+kbb(3,k)*d3)
                      ztmp=dcmplx(cos(kr),sin(kr))*unk(i,n,k,s)
                      wtmp5(1,lma,n,k,s)=wtmp5(1,lma,n,k,s)+duVdR(1,j,lma)*ztmp
                      wtmp5(2,lma,n,k,s)=wtmp5(2,lma,n,k,s)+duVdR(2,j,lma)*ztmp
                      wtmp5(3,lma,n,k,s)=wtmp5(3,lma,n,k,s)+duVdR(3,j,lma)*ztmp
                   end do
#endif
                end if
             end do ! lma
          end do ! n
!$OMP end do
       end do ! k
    end do ! s

!$OMP master
    call watch(ctt(2),ett(2))
!$OMP end master

!$OMP single
    deallocate( duVdR )
!$OMP end single

!$OMP workshare
    force2(:,:)=0.0d0
!$OMP end workshare

!$OMP single
    if ( nzlma > 0 ) then
       do s=MSP_0,MSP_1
       do k=MBZ_0,MBZ_1
       do n=MB_0,MB_1,MB_d
          ib1=n
          ib2=min(ib1+MB_d-1,MB_1)
          if ( occ(n,k,s) < 1.d-10 ) cycle
          call do3StepComm_F(nrlma_xyz,num_2_rank,sendmap,recvmap &
               ,lma_nsend,sbufnl,rbufnl,nzlma,ib1,ib2,wtmp5(0,1,ib1,k,s))
       end do ! n
       end do ! k
       end do ! s
    end if
    do s=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1
       do m=1,N_nzqr
          lma1=nzqr_pair(m,1)
          lma2=nzqr_pair(m,2)
          a=amap(lma1)
          if ( a <= 0 ) cycle
          if ( isAtomInThisNode(a) ) then
             do n=MB_0,MB_1
                Dij_f(n)=Dij(m,s)-esp(n,k,s)*qij_f(m)
                const_f(n)=Dij_f(n)*(-2.d0)*occ(n,k,s)*dV*dV
             end do
             if ( lma1 < lma2  ) stop 'Nzqr_pair is strange'
             if ( lma1 == lma2 ) then
                force2(1,a)=force2(1,a) &
                     +sum( const_f(MB_0:MB_1) &
                             *real(wtmp5(0,lma1,MB_0:MB_1,k,s) &
                                  *wtmp5(1,lma2,MB_0:MB_1,k,s),8) )
                force2(2,a)=force2(2,a) &
                     +sum( const_f(MB_0:MB_1) &
                             *real(wtmp5(0,lma1,MB_0:MB_1,k,s) &
                                  *wtmp5(2,lma2,MB_0:MB_1,k,s),8) )
                force2(3,a)=force2(3,a) &
                     +sum( const_f(MB_0:MB_1) &
                             *real( wtmp5(0,lma1,MB_0:MB_1,k,s) &
                                   *wtmp5(3,lma2,MB_0:MB_1,k,s),8) )
             else
                force2(1,a)=force2(1,a) &
                     +sum( const_f(MB_0:MB_1) &
                             *real( wtmp5(0,lma1,MB_0:MB_1,k,s) &
                                   *wtmp5(1,lma2,MB_0:MB_1,k,s),8) )
                force2(2,a)=force2(2,a) &
                     +sum( const_f(MB_0:MB_1) &
                             *real( wtmp5(0,lma1,MB_0:MB_1,k,s) &
                                   *wtmp5(2,lma2,MB_0:MB_1,k,s),8) )
                force2(3,a)=force2(3,a) &
                     +sum( const_f(MB_0:MB_1) &
                             *real( wtmp5(0,lma1,MB_0:MB_1,k,s) &
                                   *wtmp5(3,lma2,MB_0:MB_1,k,s),8) )
                force2(1,a)=force2(1,a) &
                     +sum( const_f(MB_0:MB_1) &
                             *real( wtmp5(0,lma2,MB_0:MB_1,k,s) &
                                   *wtmp5(1,lma1,MB_0:MB_1,k,s),8) )
                force2(2,a)=force2(2,a) &
                     +sum( const_f(MB_0:MB_1) &
                             *real( wtmp5(0,lma2,MB_0:MB_1,k,s) &
                                   *wtmp5(2,lma1,MB_0:MB_1,k,s),8) )
                force2(3,a)=force2(3,a) &
                     +sum( const_f(MB_0:MB_1) &
                             *real( wtmp5(0,lma2,MB_0:MB_1,k,s) &
                                   *wtmp5(3,lma1,MB_0:MB_1,k,s),8) )
             end if
          endif
       end do ! m
    end do ! k
    end do ! s

    call mpi_allreduce &
         (MPI_IN_PLACE,force2,3*MI,mpi_real8,mpi_sum,mpi_comm_world,ierr)

    allocate( uVunk_tmp(nzlma,MB_0:MB_1,MBZ_0:MBZ_1,MSP_0:MSP_1) )
    uVunk_tmp=zero
#ifdef _DRSDFT_
    uVunk_tmp(:,:,:,:)=wtmp5(0,:,:,:,:)*dV
#else
    uVunk_tmp(:,:,:,:)=conjg(wtmp5(0,:,:,:,:)*dV)
#endif

    deallocate( wtmp5 )

!$OMP end single

!$OMP master
    call watch(ctt(3),ett(3))
!$OMP end master

! -------------------- Q_part -----

    call calcForceQ(uVunk_tmp,MI,forceQ)

    do a=1,MI
       force2(1,a)=force2(1,a)+forceQ(1,a)
       force2(2,a)=force2(2,a)+forceQ(2,a)
       force2(3,a)=force2(3,a)+forceQ(3,a)
    enddo

!$OMP master
    call watch(ctt(4),ett(4))
!$OMP end master

!$OMP end parallel

!    if ( disp_switch_parallel ) then
!       write(*,*) "time(force_nloc_uspp_1)",ctt(1)-ctt(0),ett(1)-ett(0)
!       write(*,*) "time(force_nloc_uspp_2)",ctt(2)-ctt(1),ett(2)-ett(1)
!       write(*,*) "time(force_nloc_uspp_3)",ctt(3)-ctt(2),ett(3)-ett(2)
!       write(*,*) "time(force_nloc_uspp_4)",ctt(4)-ctt(3),ett(4)-ett(3)
!    end if

    return
  END SUBROUTINE calcForcePSnonLoc2


  SUBROUTINE calcForceQ(uVunk_tmp,MI,forceQ)
    implicit none
    integer,intent(IN) :: MI
#ifdef _DRSDFT_
    real(8),intent(INOUT) :: uVunk_tmp(nzlma,MB_0:MB_1,MBZ_0:MBZ_1,MSP_0:MSP_1)
#else
    complex(8),intent(IN) :: uVunk_tmp(nzlma,MB_0:MB_1,MBZ_0:MBZ_1,MSP_0:MSP_1)
#endif
    real(8),intent(INOUT) :: forceQ(3,MI)
    integer :: i1,i2,i3
    integer :: i,j,k,s,n,ir,L,L1,L1z,NRc,irank,jrank
    integer :: iqr,ll3,cJ
    integer :: ib1,ib2
    real(8) :: yq1,yq2,yq3,tmp,tmp1,tmp2,tmp3
    integer :: nreq,max_nreq
    integer :: a,a0,ik,m,lm0,lm1,lma,m1,m2
    integer :: ierr,ir0
    integer,allocatable :: ireq(:),istatus(:,:)
    real(8),parameter :: ep=1.d-8
    real(8) :: err,err0
    real(8) :: a1,a2,a3
    real(8) :: kr,c
    real(8),parameter :: pi2=2.d0*acos(-1.d0)
    real(8) :: yy1,yy2,yy3
    real(8),allocatable :: work2(:,:),dQY(:,:,:)
    real(8) :: dQY_tmp(1:3)
    real(8) :: ztmp
    real(8),allocatable :: rtmp5(:,:,:)
    integer :: i0,iorb0
    integer :: k1,k2,k3
    integer :: d1,d2,d3
    logical,save :: flag_init=.true.

    if ( flag_init ) then
       call ps_q_init_derivative
       flag_init = .false.
    end if

    forceQ(:,:) = 0.0d0

    maxerr=0.0d0

    allocate( rtmp5(0:2,N_nzqr,MSP_0:MSP_1) )
    allocate( dQY(3,MAXMJJ_MAP_Q,N_nzqr) )

!$OMP parallel

!$OMP workshare
    rtmp5=0.0d0
    dQY=0.0d0
!$OMP end workshare

!!$OMP do schedule(dynamic) firstprivate( maxerr ) &
!!$OMP    private( a,L,m,ik,Rx,Ry,Rz,NRc,d1,d2,d3,x,y,z,r  &
!!$OMP            ,ir,ir0,yy1,yy2,yy3,tmp1,m1,m2  &
!!$OMP            ,lma,j,L1,L1z,lm1 )

    do iqr=1,N_nzqr

       call getAtomInfoFrom_iqr( iqr )
      !OUT: iatom,ikind,ik1,ikk1,ik2
      !     lma1,lma2,l_1,l_2,m_1,m_2,iorb1,iorb2,noAtomHere

      if ( noAtomHere ) cycle

      call getAtomPosition_real
      !OUT: Rx,Ry,Rz
      NRc=max(NRps(iorb1,ikind),NRps(iorb2,ikind))

      do j=1,MJJ_MAP_Q(iqr)

         call getAtomCenteredPositionFrom_iqr(iqr,j)
         !OUT: x,y,z,r

         call get_nearest_index( rad1(:,ikind), r, ir )

         dQY_tmp(:)=0.0d0
         do ll3=1,nl3v(ik2,ikind)  ! ll3: info of l_1,m_1,l_2,m_2
            L=l3v(ll3,ik2,ikind)-1 ! L: {s:0,p:1,d:2,...}
            cJ=0 ! cJ is just the index to run L1 which do not start from 1
            do L1=abs(L-1),L+1
               cJ=cJ+1
               tmp0=0.0d0
               err0=0.0d0
               if ( abs(x)>1.d-14 .or. abs(y)>1.d-14 .or. &
                    abs(z)>1.d-14 .or. L1==0 ) then
                  call interpolate_dqrL(ll3,cJ,ir,NRc) !OUT: tmp0,maxerr
                  yy1=0.0d0
                  yy2=0.0d0
                  yy3=0.0d0
                  do L1z=-L1,L1
                  ! L1z: M
                  ! M summation can be done, for Q depend on L, {l1,m1,l2,m2}
                     yq1=0.0d0
                     yq2=0.0d0
                     yq3=0.0d0
                     do M=-L,L
                        yq1=yq1+GCLL%yyy(l_1,m_1,l_2,m_2,L,M)*GCL1%x(L,M,L1,L1z)
                        yq2=yq2+GCLL%yyy(l_1,m_1,l_2,m_2,L,M)*GCL1%y(L,M,L1,L1z)
                        yq3=yq3+GCLL%yyy(l_1,m_1,l_2,m_2,L,M)*GCL1%z(L,M,L1,L1z)
                     end do
                     tmp3=Ylm(x,y,z,L1,L1z)
                     yy1=yy1+yq1*tmp3
                     yy2=yy2+yq2*tmp3
                     yy3=yy3+yq3*tmp3
                  end do ! L1z
                  dQY_tmp(1)=dQY_tmp(1)+tmp0*yy1
                  dQY_tmp(2)=dQY_tmp(2)+tmp0*yy2
                  dQY_tmp(3)=dQY_tmp(3)+tmp0*yy3
               end if
            end do ! L1
         enddo ! ll3
         dQY(1,j,iqr)=dQY_tmp(1)
         dQY(2,j,iqr)=dQY_tmp(2)
         dQY(3,j,iqr)=dQY_tmp(3)
      end do ! j

   end do ! iqr
!!$OMP end do

   do s=MSP_0,MSP_1
!!$OMP do schedule(dynamic) private( i,kr,ztmp )
      do iqr=1,N_nzqr
         do j=1,MJJ_MAP_Q(iqr)
            i1=JJ_MAP_Q(1,j,iqr)
            i2=JJ_MAP_Q(2,j,iqr)
            i3=JJ_MAP_Q(3,j,iqr)
            i=i1-a1b+(i2-a2b)*ab1+(i3-a3b)*ab1*ab2+ML_0
            ztmp=Vloc(i,s)
            rtmp5(0,iqr,s)=rtmp5(0,iqr,s)+dQY(1,j,iqr)*ztmp
            rtmp5(1,iqr,s)=rtmp5(1,iqr,s)+dQY(2,j,iqr)*ztmp
            rtmp5(2,iqr,s)=rtmp5(2,iqr,s)-dQY(3,j,iqr)*ztmp
         end do ! j
      end do ! iqr
!!$OMP end do
   end do ! s

!$OMP single
   deallocate( dQY )
!$OMP end single

!$OMP single
    do s=MSP_0,MSP_1
       call comm_dQRij( rtmp5(:,:,s) )
    enddo

    do s=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1

       do m=1,N_nzqr

          lma1 = nzqr_pair(m,1)
          lma2 = nzqr_pair(m,2)

          a = amap(lma1)
          if ( a <= 0 ) cycle
          if ( .not.isAtomInThisNode(a) ) cycle
          if ( lma1 < lma2 ) stop 'Nzqr_pair is strange'

          tmp1 = -rtmp5(0,m,s)*dV
          tmp2 = -rtmp5(1,m,s)*dV
          tmp3 = -rtmp5(2,m,s)*dV

          if ( lma1 == lma2 ) then
             forceQ(1,a)=forceQ(1,a) &
                  +sum( occ(MB_0:MB_1,k,s) &
                       *abs(uVunk_tmp(lma1,MB_0:MB_1,k,s))**2 )*tmp1
             forceQ(2,a)=forceQ(2,a) &
                  +sum( occ(MB_0:MB_1,k,s) &
                       *abs(uVunk_tmp(lma1,MB_0:MB_1,k,s))**2 )*tmp2
             forceQ(3,a)=forceQ(3,a) &
                  +sum( occ(MB_0:MB_1,k,s) &
                       *abs(uVunk_tmp(lma1,MB_0:MB_1,k,s))**2 )*tmp3
          else
#ifdef _DRSDFT_
             forceQ(1,a)=forceQ(1,a) &
                  +2.0d0*sum(occ(MB_0:MB_1,k,s) &
                  *real(uVunk_tmp(lma1,MB_0:MB_1,k,s) &
                       *uVunk_tmp(lma2,MB_0:MB_1,k,s),8))*tmp1
             forceQ(2,a)=forceQ(2,a) &
                  +2.0d0*sum(occ(MB_0:MB_1,k,s) &
                  *real(uVunk_tmp(lma1,MB_0:MB_1,k,s) &
                       *uVunk_tmp(lma2,MB_0:MB_1,k,s),8))*tmp2
             forceQ(3,a)=forceQ(3,a) &
                  +2.0d0*sum(occ(MB_0:MB_1,k,s) &
                  *real(uVunk_tmp(lma1,MB_0:MB_1,k,s) &
                       *uVunk_tmp(lma2,MB_0:MB_1,k,s),8))*tmp3
#else
             forceQ(1,a)=forceQ(1,a)+sum(occ(MB_0:MB_1,k,s)* &
                  & real(conjg(uVunk_tmp(lma1,MB_0:MB_1,k,s)) &
                              *uVunk_tmp(lma2,MB_0:MB_1,k,s)))*tmp1
             forceQ(2,a)=forceQ(2,a)+sum(occ(MB_0:MB_1,k,s)* &
                  & real(conjg(uVunk_tmp(lma1,MB_0:MB_1,k,s)) &
                              *uVunk_tmp(lma2,MB_0:MB_1,k,s)))*tmp2
             forceQ(3,a)=forceQ(3,a)+sum(occ(MB_0:MB_1,k,s)* &
                  & real(conjg(uVunk_tmp(lma1,MB_0:MB_1,k,s)) &
                              *uVunk_tmp(lma2,MB_0:MB_1,k,s)))*tmp3
             forceQ(1,a)=forceQ(1,a)+sum(occ(MB_0:MB_1,k,s)* &
                  & real(conjg(uVunk_tmp(lma2,MB_0:MB_1,k,s)) &
                              *uVunk_tmp(lma1,MB_0:MB_1,k,s)))*tmp1
             forceQ(2,a)=forceQ(2,a)+sum(occ(MB_0:MB_1,k,s)* &
                  & real(conjg(uVunk_tmp(lma2,MB_0:MB_1,k,s)) &
                              *uVunk_tmp(lma1,MB_0:MB_1,k,s)))*tmp2
             forceQ(3,a)=forceQ(3,a)+sum(occ(MB_0:MB_1,k,s)* &
                  & real(conjg(uVunk_tmp(lma2,MB_0:MB_1,k,s)) &
                              *uVunk_tmp(lma1,MB_0:MB_1,k,s)))*tmp3
#endif
          end if

       end do ! m

    end do ! k
    end do ! s

!$OMP end single

!$OMP single
    call mpi_allreduce( MPI_IN_PLACE, forceQ, size(forceQ) &
                      , mpi_real8,mpi_sum,mpi_comm_world,ierr )
!$OMP end single

!$OMP single
    deallocate( rtmp5 )
!$OMP end single

!$OMP end parallel

    return
  END SUBROUTINE calcForceQ


  SUBROUTINE get_nearest_index( rad, r, ir )
    implicit none
    real(8),intent(IN)  :: rad(:), r
    integer,intent(OUT) :: ir
    integer :: ir0,ir1
    ir0=1
    ir1=size( rad )
1   if ( ir1 - ir0 > 1 ) then
       ir = ( ir0 + ir1 )/2
       if ( rad(ir) > r ) then
          ir1=ir
       else
          ir0=ir
       end if
       goto 1
    end if
  END SUBROUTINE get_nearest_index


  SUBROUTINE comm_dQRij( f )
    implicit none
    real(8),intent(INOUT) :: f(:,:)
    integer :: iqr,lma1,lma2,a1,i1,i2,m1,m2,ierr
    real(8),allocatable :: ftmp(:,:,:,:,:,:)

    i1 = maxval( norb )
    i2 = maxval( lo )
    allocate( ftmp(3,Natom,i1,i1,-i2:i2,-i2:i2) ) ; ftmp=0.0d0

    do iqr=1,N_nzqr

       lma1 = nzqr_pair(iqr,1)
       lma2 = nzqr_pair(iqr,2)

       a1 = amap(lma1)

       i1 = iorbmap(lma1)
       i2 = iorbmap(lma2)

       m1 = mmap(lma1)
       m2 = mmap(lma2)

       ftmp(:,a1,i1,i2,m1,m2) = ftmp(:,a1,i1,i2,m1,m2) + f(:,iqr)

    end do ! iqr

    call MPI_ALLREDUCE( MPI_IN_PLACE, ftmp, size(ftmp), MPI_REAL8 &
                       ,MPI_SUM, comm_grid, ierr )

    do iqr=1,N_nzqr

       lma1 = nzqr_pair(iqr,1)
       lma2 = nzqr_pair(iqr,2)

       a1 = amap(lma1)

       i1 = iorbmap(lma1)
       i2 = iorbmap(lma2)

       m1 = mmap(lma1)
       m2 = mmap(lma2)

       f(:,iqr) = ftmp(:,a1,i1,i2,m1,m2)

    end do ! iqr

    deallocate( ftmp )

  END SUBROUTINE comm_dQRij


END MODULE force_ps_nloc2
