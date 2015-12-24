MODULE force_nloc2_mol_module

  use parallel_module
  use rgrid_module
  use rgrid_mol_module
  use wf_module
  use electron_module
  use array_bound_module
  use atom_module
  use pseudopot_module
  use ps_nloc2_variables
  use ps_nloc2_init_module
  use ylm_module
  use polint_module

  implicit none

  PRIVATE
  PUBLIC :: calc_force_nloc2_mol

  logical :: flag_init = .true.
  integer,allocatable :: ilm1(:,:,:),irad(:,:)
  real(8) :: Y1(0:3,-3:3,0:4,-4:4)
  real(8) :: Y2(0:3,-3:3,0:4,-4:4)
  real(8) :: Y3(0:3,-3:3,0:4,-4:4)

CONTAINS


  SUBROUTINE calc_force_nloc2_mol(force)
    implicit none
    real(8),intent(INOUT) :: force(:,:)
    integer :: lma,ir,ir0,L1,L1z,j,lm1,s,n,i,i1,i2,i3,i4
    integer :: NRc,a,L,m,iorb,ik,m1,m2,mm,im,ierr
    integer :: irank,jrank,nreq
    integer,allocatable :: istatus(:,:),ireq(:)
    real(8) :: maxerr,err,err0,x,y,z,r
    real(8) :: Rx,Ry,Rz,Rc,c
    real(8) :: v,v0,yy1,yy2,yy3
    real(8),allocatable :: duVdR(:,:,:)
#ifdef _DRSDFT_
    real(8),allocatable :: w(:,:,:,:),wtmp(:,:)
    real(8),parameter :: zero=0.0d0
#else
    complex(8),allocatable :: w(:,:,:,:),wtmp(:,:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
#endif
    real(8),allocatable :: ftmp(:,:)
    logical,allocatable :: a_rank(:)

    if ( Mlma <= 0 ) return

    if ( flag_init ) then
       call init_force_nloc2_mol
       call init_force_nloc2_derivative_mol
       flag_init = .false.
    end if

    allocate( ftmp(3,Natom) ) ; ftmp=0.0d0

    allocate( duVdR(3,MMJJ,nzlma) ) ; duVdR=0.0d0

    maxerr = 0.0d0

    do lma=1,nzlma

       a    = amap(lma) ; if ( a <= 0 ) cycle
       L    = lmap(lma)
       m    = mmap(lma)
       iorb = iorbmap(lma)
       ik   = ki_atom(a)

       Rx = aa_atom(1,a)
       Ry = aa_atom(2,a)
       Rz = aa_atom(3,a)

       Rc  = Rps(iorb,ik)
       NRc = NRps(iorb,ik)

       do j=1,MJJ_MAP(lma)

          x = JJ_MAP(1,j,lma)*Hsize - Rx
          y = JJ_MAP(2,j,lma)*Hsize - Ry
          z = JJ_MAP(3,j,lma)*Hsize - Rz
          r = sqrt(x*x+y*y+z*z)

          ir0 = irad( int(100.d0*r), ik )
          do ir=ir0,NRc
             if ( r < rad1(ir,ik) ) exit
          end do
          if ( ir < 1 ) stop "stop@calc_force_nloc2_mol(1)"

          yy1 = 0.0d0
          yy2 = 0.0d0
          yy3 = 0.0d0

          do L1=abs(L-1),L+1

             lm1 = ilm1(L1,iorb,ik)

             if ( abs(x) > 1.d-14 .or. abs(y) > 1.d-14 .or. &
                  abs(z) > 1.d-14 .or. L1 == 0 ) then

                if ( r < rad1(2,ik) ) then

                   v0 = dviod(2,lm1,ik)/rad1(2,ik)**2
                   err0 = 0.0d0

                else if ( ir <= NRc ) then

                   err0 = 1.d10
                   v0 = 0.0d0
                   do im=1,20
                      m1=max(1,ir-im)
                      m2=min(ir+im,NRc)
                      mm=m2-m1+1
                      call polint(rad1(m1,ik),dviod(m1,lm1,ik),mm,r,v,err)
                      if ( abs(err) < err0 ) then
                         v0 = v
                         err0 = abs(err)
                         if ( err0 < 1.d-8 ) exit
                      end if
                   end do
                   v0 = v0/r**2

                else

                   stop "stop@calc_force_nloc2_mol(2)"

                end if
                maxerr = max(maxerr,err0)

             end if

             do L1z=-L1,L1
                v = v0*Ylm(x,y,z,L1,L1z)
                yy1 = yy1 + v*Y1(L,m,L1,L1z)
                yy2 = yy2 + v*Y2(L,m,L1,L1z)
                yy3 = yy3 - v*Y3(L,m,L1,L1z)
             end do

          end do ! L1

          duVdR(1,j,lma) = yy1
          duVdR(2,j,lma) = yy2
          duVdR(3,j,lma) = yy3

       end do ! j

    end do ! lma

!--

    allocate( w(0:3,nzlma,MB_0:MB_1,MSP_0:MSP_1) ) ; w=zero

    do s=MSP_0,MSP_1
    do n=MB_0,MB_1

       if ( occ(n,1,s) == 0.0d0 ) cycle

       c = -2.0d0*occ(n,1,s)*dV*dV

       do lma=1,nzlma

          do j=1,MJJ(lma)
             i=JJP(j,lma)
#ifdef _DRSDFT_
             w(0,lma,n,s) = w(0,lma,n,s) + uVk(j,lma,1)*unk(i,n,1,s)
#else             
             w(0,lma,n,s) = w(0,lma,n,s) + uVk(j,lma,1)*conjg(unk(i,n,1,s))
#endif
             w(1,lma,n,s) = w(1,lma,n,s) + duVdR(1,j,lma)*unk(i,n,1,s)
             w(2,lma,n,s) = w(2,lma,n,s) + duVdR(2,j,lma)*unk(i,n,1,s)
             w(3,lma,n,s) = w(3,lma,n,s) + duVdR(3,j,lma)*unk(i,n,1,s)
          end do ! j
          w(0,lma,n,s) = iuV(lma)*c*w(0,lma,n,s)

       end do ! lma

    end do ! n
    end do ! s

    deallocate( duVdR )

!--

    m = 2*maxval( nrlma_xyz )
    allocate( ireq(m) ) ; ireq=0
    allocate( istatus(MPI_STATUS_SIZE,m) ) ; istatus=0
    allocate( wtmp(0:3,nzlma) ) ; wtmp=zero
    allocate( a_rank(Natom) ) ; a_rank=.false.

    do a=1,Natom
       i1=nint(aa_atom(1,a)/Hsize)
       i2=nint(aa_atom(2,a)/Hsize)
       i3=nint(aa_atom(3,a)/Hsize)
       if ( Igrid(1,1) <= i1 .and. i1 <= Igrid(2,1) .and. &
            Igrid(1,2) <= i2 .and. i2 <= Igrid(2,2) .and. &
            Igrid(1,3) <= i3 .and. i3 <= Igrid(2,3) ) then
          a_rank(a) = .true.
       end if
    end do

    do s=MSP_0,MSP_1
    do n=MB_0,MB_1

       if ( occ(n,1,s) == 0.d0 ) cycle

       if ( iswitch_eqdiv == 1 ) then

       do i=1,6
          select case(i)
          case(1,3,5)
             j=i+1
             wtmp(:,:)=w(:,:,n,s)
          case(2,4,6)
             j=i-1
          end select
          do m=1,nrlma_xyz(i)
             nreq=0
             irank=num_2_rank(m,i)
             jrank=num_2_rank(m,j)
             if ( irank >= 0 ) then
                i1=0
                do i2=1,lma_nsend(irank)
                do i3=0,3
                   i1=i1+1
                   sbufnl(i1,irank)=wtmp(i3,sendmap(i2,irank))
                end do
                end do
                nreq=nreq+1
                call mpi_isend(sbufnl(1,irank),4*lma_nsend(irank) &
                     ,TYPE_MAIN,irank,1,comm_grid,ireq(nreq),ierr)
             end if
             if ( jrank >= 0 ) then
                nreq=nreq+1
                call mpi_irecv(rbufnl(1,jrank),4*lma_nsend(jrank) &
                     ,TYPE_MAIN,jrank,1,comm_grid,ireq(nreq),ierr)
             end if
             call mpi_waitall(nreq,ireq,istatus,ierr)
             if ( jrank >= 0 ) then
                i1=0
                do i2=1,lma_nsend(jrank)
                   i4=recvmap(i2,jrank)
                do i3=0,3
                   i1=i1+1
                   w(i3,i4,n,s) = w(i3,i4,n,s) + rbufnl(i1,jrank)
                end do
                end do
             end if
          end do ! m
       end do ! i

       else if ( iswitch_eqdiv == 2 ) then

          call comm_eqdiv_ps_nloc2_mol(nzlma,1,1,w(0,1,n,s))

       end if ! iswitch_eqdiv

       do lma=1,nzlma
          a=amap(lma)
          if ( a <= 0 ) cycle
          if ( .not.a_rank(a) ) cycle
          ftmp(1,a) = ftmp(1,a) + w(0,lma,n,s)*w(1,lma,n,s)
          ftmp(2,a) = ftmp(2,a) + w(0,lma,n,s)*w(2,lma,n,s)
          ftmp(3,a) = ftmp(3,a) + w(0,lma,n,s)*w(3,lma,n,s)
       end do

    end do ! n
    end do ! s

    deallocate( a_rank )
    deallocate( wtmp )
    deallocate( istatus )
    deallocate( ireq )

    deallocate( w )

    call mpi_allreduce &
         (MPI_IN_PLACE,ftmp,3*Natom,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

    force(:,:) = ftmp(:,:)

    deallocate( ftmp )

  END SUBROUTINE calc_force_nloc2_mol


  SUBROUTINE init_force_nloc2_mol
    implicit none
    integer :: L,L1,n,lm1,ik,iorb
    integer :: M_irad,m,NRc,i,ir

    Y1=0.d0
    Y1( 0, 0, 1, 1) =  0.282094791773878d0
    Y1( 1,-1, 2,-2) = -0.218509686118416d0
    Y1( 1, 0, 2, 1) =  0.218509686118416d0
    Y1( 1, 1, 0, 0) =  0.282094791773878d0
    Y1( 1, 1, 2, 0) = -0.126156626101008d0
    Y1( 1, 1, 2, 2) =  0.218509686118416d0
    Y1( 2,-2, 1,-1) = -0.218509686118416d0
    Y1( 2,-2, 3,-3) = -0.226179013159540d0
    Y1( 2,-2, 3,-1) =  0.058399170081902d0
    Y1( 2,-1, 3,-2) = -0.184674390922372d0
    Y1( 2, 0, 1, 1) = -0.126156626101008d0
    Y1( 2, 0, 3, 1) =  0.202300659403421d0
    Y1( 2, 1, 1, 0) =  0.218509686118416d0
    Y1( 2, 1, 3, 0) = -0.143048168102669d0
    Y1( 2, 1, 3, 2) =  0.184674390922372d0
    Y1( 2, 2, 1, 1) =  0.218509686118416d0
    Y1( 2, 2, 3, 1) = -0.058399170081902d0
    Y1( 2, 2, 3, 3) =  0.226179013159540d0
    Y1( 3,-3, 2,-2) = -0.226179013159540d0
    Y1( 3,-1, 2,-2) =  0.058399170081901d0
    Y1( 3,-2, 2,-1) = -0.184674390922371d0
    Y1( 3, 1, 2, 0) =  0.202300659403420d0
    Y1( 3, 0, 2, 1) = -0.143048168102668d0
    Y1( 3, 2, 2, 1) =  0.184674390922371d0
    Y1( 3, 1, 2, 2) = -0.058399170081901d0
    Y1( 3, 3, 2, 2) =  0.226179013159540d0
    Y1( 3,-3, 4,-4) = -0.230329432980890d0
    Y1( 3,-2, 4,-3) = -0.199471140200716d0
    Y1( 3,-3, 4,-2) =  0.043528171377568d0
    Y1( 3,-1, 4,-2) = -0.168583882836183d0
    Y1( 3,-2, 4,-1) =  0.075393004386513d0
    Y1( 3, 1, 4, 0) = -0.150786008773026d0
    Y1( 3, 0, 4, 1) =  0.194663900273006d0
    Y1( 3, 2, 4, 1) = -0.075393004386513d0
    Y1( 3, 1, 4, 2) =  0.168583882836183d0
    Y1( 3, 3, 4, 2) = -0.043528171377568d0
    Y1( 3, 2, 4, 3) =  0.199471140200716d0
    Y1( 3, 3, 4, 4) =  0.230329432980890d0
    Y2=0.d0
    Y2( 0, 0, 1,-1) =  0.282094791773878d0
    Y2( 1,-1, 0, 0) =  0.282094791773878d0
    Y2( 1,-1, 2, 0) = -0.126156626101008d0
    Y2( 1,-1, 2, 2) = -0.218509686118416d0
    Y2( 1, 0, 2,-1) =  0.218509686118416d0
    Y2( 1, 1, 2,-2) = -0.218509686118416d0
    Y2( 2,-2, 1, 1) = -0.218509686118416d0
    Y2( 2,-2, 3, 1) =  0.058399170081902d0
    Y2( 2,-2, 3, 3) =  0.226179013159540d0
    Y2( 2,-1, 1, 0) =  0.218509686118416d0
    Y2( 2,-1, 3, 0) = -0.143048168102669d0
    Y2( 2,-1, 3, 2) = -0.184674390922372d0
    Y2( 2, 0, 1,-1) = -0.126156626101008d0
    Y2( 2, 0, 3,-1) =  0.202300659403421d0
    Y2( 2, 1, 3,-2) = -0.184674390922372d0
    Y2( 2, 2, 1,-1) = -0.218509686118416d0
    Y2( 2, 2, 3,-3) =  0.226179013159540d0
    Y2( 2, 2, 3,-1) =  0.058399170081902d0
    Y2( 3, 1, 2,-2) =  0.058399170081901d0
    Y2( 3, 3, 2,-2) =  0.226179013159540d0
    Y2( 3, 0, 2,-1) = -0.143048168102668d0
    Y2( 3, 2, 2,-1) = -0.184674390922371d0
    Y2( 3,-1, 2, 0) =  0.202300659403420d0
    Y2( 3,-2, 2, 1) = -0.184674390922371d0
    Y2( 3,-3, 2, 2) =  0.226179013159540d0
    Y2( 3,-1, 2, 2) =  0.058399170081901d0
    Y2( 3, 3, 4,-4) = -0.230329432980890d0
    Y2( 3, 2, 4,-3) =  0.199471140200716d0
    Y2( 3, 1, 4,-2) = -0.168583882836183d0
    Y2( 3, 3, 4,-2) = -0.043528171377568d0
    Y2( 3, 0, 4,-1) =  0.194663900273006d0
    Y2( 3, 2, 4,-1) =  0.075393004386513d0
    Y2( 3,-1, 4, 0) = -0.150786008773026d0
    Y2( 3,-2, 4, 1) =  0.075393004386513d0
    Y2( 3,-3, 4, 2) = -0.043528171377568d0
    Y2( 3,-1, 4, 2) = -0.168583882836183d0
    Y2( 3,-2, 4, 3) =  0.199471140200716d0
    Y2( 3,-3, 4, 4) = -0.230329432980890d0
    Y3=0.d0
    Y3( 0, 0, 1, 0) =  0.282094791773878d0
    Y3( 1,-1, 2,-1) =  0.218509686118416d0
    Y3( 1, 0, 0, 0) =  0.282094791773878d0
    Y3( 1, 0, 2, 0) =  0.252313252202016d0
    Y3( 1, 1, 2, 1) =  0.218509686118416d0
    Y3( 2,-2, 3,-2) =  0.184674390922372d0
    Y3( 2,-1, 1,-1) =  0.218509686118416d0
    Y3( 2,-1, 3,-1) =  0.233596680327607d0
    Y3( 2, 0, 1, 0) =  0.252313252202016d0
    Y3( 2, 0, 3, 0) =  0.247766695083476d0
    Y3( 2, 1, 1, 1) =  0.218509686118416d0
    Y3( 2, 1, 3, 1) =  0.233596680327607d0
    Y3( 2, 2, 3, 2) =  0.184674390922372d0
    Y3( 3,-2, 2,-2) =  0.184674390922371d0
    Y3( 3,-1, 2,-1) =  0.233596680327607d0
    Y3( 3, 0, 2, 0) =  0.247766695083476d0
    Y3( 3, 1, 2, 1) =  0.233596680327607d0
    Y3( 3, 2, 2, 2) =  0.184674390922371d0
    Y3( 3,-3, 4,-3) =  0.162867503967639d0
    Y3( 3,-2, 4,-2) =  0.213243618622923d0
    Y3( 3,-1, 4,-1) =  0.238413613504448d0
    Y3( 3, 0, 4, 0) =  0.246232521229829d0
    Y3( 3, 1, 4, 1) =  0.238413613504448d0
    Y3( 3, 2, 4, 2) =  0.213243618622923d0
    Y3( 3, 3, 4, 3) =  0.162867503967639d0

    L1=maxval(lo)+1
    n =maxval(norb)
    allocate( ilm1(0:L1,n,Nelement) ) ; ilm1=0
    do ik=1,Nelement
       lm1=0
       do iorb=1,norb(ik)
          L=lo(iorb,ik)
          do L1=abs(L-1),L+1
             lm1=lm1+1
             ilm1(L1,iorb,ik)=lm1
          end do
       end do
    end do

    allocate( irad(0:3000,Nelement) ) ; irad=0
    M_irad=0
    do ik=1,Nelement
       NRc=maxval( NRps(:,ik) )
       NRc=min( 3000, NRc )
       m=0
       irad(0,ik)=1
       do ir=1,NRc
          m=int(100.d0*rad1(ir,ik))+1
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

  END SUBROUTINE init_force_nloc2_mol


  SUBROUTINE init_force_nloc2_derivative_mol
    implicit none
    integer :: lm,ik,m,iorb,i,L,NRc,J,im,m1,m2,mm
    real(8) :: c,err,err0,v,v0,maxerr,r
    real(8),allocatable :: dvrad(:,:,:) 

    lm=0
    do ik=1,Nelement
       m=0
       do iorb=1,norb(ik)
          if ( lo(iorb,ik)==0 ) then
             m=m+1
          else
             m=m+3
          end if
       end do
       lm=max(m,lm)
    end do

    NRc=maxval(NRps)
    allocate( dviod(NRc,lm,Nelement) ) ; dviod=0.0d0
    allocate( dvrad(NRc,lm,Nelement) ) ; dvrad=0.0d0

    maxerr=0.0d0
    do ik=1,Nelement
       do iorb=1,norb(ik)
          L=lo(iorb,ik)
          NRc=NRps(iorb,ik)
          do i=1,NRc
             err0=1.d10
             v0=0.0d0
             do im=1,20
                m1=max(i-im,1) 
                m2=min(i+im,NRc)
                mm=m2-m1+1
                call dpolint(rad1(m1,ik),viod(m1,iorb,ik),mm,rad1(i,ik),v,err)
                if ( abs(err) < err0 ) then
                   v0=v
                   err0=abs(err)
                end if
             end do ! im
             dviod(i,iorb,ik)=v0
             maxerr=max(maxerr,err0)
          end do ! i
       end do ! iorb
    end do ! ik

    do ik=1,Nelement
       lm=0
       do iorb=1,norb(ik)
          L=lo(iorb,ik)
          NRc=NRps(iorb,ik)
          do J=abs(L-1),L+1
             lm=lm+1
             c=0.5d0*(2.0d0+L*(L+1)-J*(J+1))
             do i=1,NRc
                r=rad1(i,ik)
                dvrad(i,lm,ik) = r*r*dviod(i,iorb,ik) + c*r*viod(i,iorb,ik)
             end do ! i
          end do ! J
       end do ! iorb
    end do ! ik

    c = sqrt( 4.0d0*acos(-1.d0)/3.0d0 )
    do ik=1,Nelement
       lm=0
       do iorb=1,norb(ik)
          L=lo(iorb,ik)
          NRc=NRps(iorb,ik)
          do J=abs(L-1),L+1
             lm=lm+1
             do i=1,NRc
                dviod(i,lm,ik) = c*dvrad(i,lm,ik)
             end do ! i
          end do ! J
       end do ! iorb
    end do ! ik

    deallocate( dvrad )

  END SUBROUTINE init_force_nloc2_derivative_mol


  SUBROUTINE comm_eqdiv_ps_nloc2_mol(nzlma,ib1,ib2,w)
    implicit none
    integer,intent(IN) :: nzlma,ib1,ib2
#ifdef _DRSDFT_
    real(8),intent(INOUT) :: w(4,nzlma,ib1:ib2)
#else
    complex(8),intent(INOUT) :: w(4,nzlma,ib1:ib2)
#endif
    integer :: nreq,irank,m,i1,i2,i3,i4,ib,ierr,nb
    integer :: istatus(mpi_status_size,512),ireq(512)
    nb=ib2-ib1+1
    nreq=0
    do irank=0,nprocs_g-1
       m=lma_nsend(irank)*4
       if ( irank == myrank_g .or. m <= 0 ) cycle
       i2=0
       do ib=ib1,ib2
       do i1=1,lma_nsend(irank)
       do i3=1,4
          i2=i2+1
          sbufnl(i2,irank)=w(i3,sendmap(i1,irank),ib)
       end do
       end do
       end do
       nreq=nreq+1
       call mpi_isend(sbufnl(1,irank),m*nb,TYPE_MAIN,irank,1 &
            ,comm_grid,ireq(nreq),ierr)
       nreq=nreq+1
       call mpi_irecv(rbufnl(1,irank),m*nb,TYPE_MAIN,irank,1 &
            ,comm_grid,ireq(nreq),ierr)
    end do
    if ( nreq > 0 ) call mpi_waitall(nreq,ireq,istatus,ierr)
    do irank=0,nprocs_g-1
       m=lma_nsend(irank)*4
       if ( irank == myrank_g .or. m <= 0 ) cycle
       i2=0
       do ib=ib1,ib2
       do i1=1,lma_nsend(irank)
          i4=recvmap(i1,irank)
       do i3=1,4
          i2=i2+1
          w(i3,i4,ib)=w(i3,i4,ib)+rbufnl(i2,irank)
       end do
       end do
       end do
    end do
  END SUBROUTINE comm_eqdiv_ps_nloc2_mol


END MODULE force_nloc2_mol_module
