MODULE ps_nloc_mr_module

  use aa_module
  use atom_module
  use rgrid_module
  use parallel_module
  use array_bound_module
  use pseudopot_module
  use ps_nloc2_init_module, only: rad1,dviod
  use ps_nloc2_variables
  use ps_nloc_hgh_module
  use minimal_box_module
  use bz_module
  use watch_module
  use wf_module
  use ylm_module
  use hsort_module
  use polint_module
  use spline_module

  implicit none

  PRIVATE
  PUBLIC :: prep_ps_nloc_mr &
           ,op_ps_nloc_mr,calc_force_ps_nloc_mr &
           ,prep_uvk_ps_nloc_mr,prep_rvk_ps_nloc_mr

  real(8),allocatable :: y2a(:,:,:),y2b(:,:,:)
  integer,allocatable :: ilm1(:,:,:)

  integer :: N_nzqr
  integer,allocatable :: nzqr_pair(:,:)
  real(8),allocatable :: Dij00(:)

CONTAINS


  SUBROUTINE prep_ps_nloc_mr
    implicit none
    complex(8) :: ztmp0
    integer,allocatable :: icheck_tmp1(:),icheck_tmp2(:),itmp(:,:)
    integer,allocatable :: icheck_tmp3(:,:,:),icheck_tmp4(:,:,:)
    integer,allocatable :: sendmap_tmp(:,:),recvmap_tmp(:,:),ireq(:)
    integer,allocatable :: lma_nsend_tmp(:),maps_tmp(:,:),itmp1(:)
    integer,allocatable :: irad(:,:),nl_rank_map_tmp(:),itmp3(:,:)
    integer,allocatable :: itmp2(:),LLp(:,:)
    integer,allocatable :: JJ_tmp(:,:,:,:),MJJ_tmp(:,:)
    integer,allocatable :: jtmp3(:,:,:),mtmp3(:),istatus(:,:)
    integer :: a,i,j,k,L,m,n,mm1,mm2,mm3,m1,m2,ML0,k1,k2,k3
    integer :: i1,i2,i3,j1,j2,j3,ik,ir,m0,iorb,mm,ierr,ir0,irlma
    integer :: ic1,ic2,ic3,id1,id2,id3,ii1,ii2,ii3,iii1,iii2,iii3,lma1,lma2
    integer :: Nintp_0,nzlma_0,M_irad,NRc,MMJJ_0,lma,lma0,i1_0,i2_0,i3_0
    integer :: nreq,ibuf(3,3),irank
    real(8),parameter :: ep=1.d-8
    real(8) :: x,y,z,r,Rx,Ry,Rz,Rps2,v,v0,d1,d2,d3,r2,kr,pi2
    real(8) :: tmp0,tmp1,tmp2,tmp3,c1,c2,c3,maxerr,err0,err
    real(8),allocatable :: uV_tmp(:,:,:),work(:)
    real(8) :: ctt(0:5),ett(0:5)
    integer :: ML1,ML2,ML3,a1b,b1b,a2b,b2b,a3b,b3b
    integer :: ab1,ab2,ab3,a1,a2,l1,l2
    integer :: np1,np2,np3,nrlma
    logical,allocatable :: lcheck_tmp1(:,:)

    Mlma=0
    do i=1,Natom
       ik=ki_atom(i)
       do iorb=1,norb(ik)
          Mlma=Mlma+2*lo(iorb,ik)+1
       end do
    end do

    if ( .not.allocated(y2a) .and. all(ippform /= 4) ) then
       NRc=maxval(NRps)
       n=maxval(norb)
       allocate( y2a(NRc,n,Nelement) )
       y2a=0.d0
       do ik=1,Nelement
       do iorb=1,norb(ik)
          d1=0.d0
          d2=0.d0
          call spline(rad1(1,ik),viod(1,iorb,ik),NRps(iorb,ik),d1,d2,y2a(1,iorb,ik))
       end do
       end do
    end if

    if ( all(ippform == 4) ) call init_ps_nloc_hgh( disp_switch_parallel ) 

    if ( Mlma < nprocs_g ) then
       nzlma_0 = Mlma
    else
       nzlma_0 = min(Mlma*125/nprocs_g,Mlma)
    end if

    ctt(:)=0.d0
    ett(:)=0.d0

    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)
    ab1 = Igrid(2,1)-Igrid(1,1)+1
    ab2 = Igrid(2,2)-Igrid(1,2)+1
    ab3 = Igrid(2,3)-Igrid(1,3)+1

    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    r=maxval(Rps)+maxval(Hgrid(1:3))+1.d-8
    call make_minimal_box(r,mm1,mm2,mm3,MMJJ_0)
    mm1 = maxval( abs(mcube_grid_ion(:,1)) ) + 1
    mm2 = maxval( abs(mcube_grid_ion(:,2)) ) + 1
    mm3 = maxval( abs(mcube_grid_ion(:,3)) ) + 1

    MMJJ_0 = M_grid_ion

    L=maxval(lo)
    n=maxval(norb)
    allocate( icheck_tmp3(Natom,n,2*L+1) ) ; icheck_tmp3=0
    allocate( JJ_tmp(6,MMJJ_0,n,Natom)   ) ; JJ_tmp=0
    allocate( MJJ_tmp(n,Natom)           ) ; MJJ_tmp=0
    allocate( uV_tmp(MMJJ_0,n,Natom)     ) ; uV_tmp=0.d0

    call watch(ctt(0),ett(0))

    if ( any(ippform == 4) ) then

       call prep_ps_nloc_hgh(Natom,n,L,MMJJ_0,M_grid_ion,map_grid_ion &
                            ,icheck_tmp3,JJ_tmp,MJJ_tmp,uV_tmp,nzlma,MMJJ)

       if ( .not.all( ippform == 4 ) ) then
          write(*,*) "Mixed use of different pseudopotenial is forbidden"
          stop "stop@prep_ps_nloc_mr"
       end if

    else

#ifndef _SPLINE_

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

#endif

       c1                 = 1.d0/Ngrid(1)
       c2                 = 1.d0/Ngrid(2)
       c3                 = 1.d0/Ngrid(3)
       maxerr             = 0
       icheck_tmp3(:,:,:) = 0
       MMJJ               = 0
       nzlma              = 0
       lma                = 0
       lma0               = 0

!$OMP parallel do schedule(dynamic) firstprivate( maxerr ) &
!$OMP    private( Rx,Ry,Rz,ic1,ic2,ic3,ik,iorb,Rps2,NRc,L,j,i,i1,i2,i3 &
!$OMP            ,id1,id2,id3,k1,k2,k3,i1_0,i2_0,i3_0,d1,d2,d3,x,y,z,r2,r &
!$OMP            ,v0,err0,ir0,ir,mm,m1,m2,v,err )
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
             NRc  = NRps(iorb,ik)
             L    = lo(iorb,ik)
             j    = 0

             do i=1,M_grid_ion

                i1 = map_grid_ion(1,i)
                i2 = map_grid_ion(2,i)
                i3 = map_grid_ion(3,i)

                id1 = ic1 + i1
                id2 = ic2 + i2
                id3 = ic3 + i3

                k1=id1/ML1 ; if ( id1<0 ) k1=(id1+1)/ML1-1
                k2=id2/ML2 ; if ( id2<0 ) k2=(id2+1)/ML2-1
                k3=id3/ML3 ; if ( id3<0 ) k3=(id3+1)/ML3-1
                i1_0=id1-k1*ML1
                i2_0=id2-k2*ML2
                i3_0=id3-k3*ML3

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

                   if ( r2 > Rps2+1.d-10 ) cycle

                   r    = sqrt(r2)
                   v0   = 0.d0
                   err0 = 0.d0

                   if ( abs(x)>1.d-14 .or. abs(y)>1.d-14 .or. &
                        abs(z)>1.d-14 .or. L==0 ) then
#ifdef _SPLINE_
                      call splint(rad1(1,ik),viod(1,iorb,ik),y2a(1,iorb,ik),NRc,r,v0)
#else
                      ir0=irad( int(100.d0*r),ik )
                      do ir=ir0,NRc
                         if ( r<rad1(ir,ik) ) exit
                      end do
                      if ( ir <= 2 ) then
                         v0=viod(2,iorb,ik)
                         if ( ir < 1 ) stop "ps_nloc_mr(0)"
                      else if ( ir <= NRc ) then
                         err0=1.d10
                         do mm=1,20
                            m1=max(1,ir-mm)
                            m2=min(ir+mm,NRc)
                            call polint &
                                 (rad1(m1,ik),viod(m1,iorb,ik),m2-m1+1,r,v,err)
                            if ( abs(err)<err0 ) then
                               v0=v
                               err0=abs(err)
                               if ( err0<ep ) exit
                            end if
                         end do
                      else
                         write(*,*) "ps_nloc_mr(1)",ir,NRc,viod(NRc,iorb,ik)
                         write(*,*) viod(NRc+1,iorb,ik),r,rad1(ir,ik)
                         stop
                      end if
                      maxerr=max(maxerr,err0)
#endif
                   end if

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
!$OMP end parallel do

#ifndef _SPLINE_
       deallocate( irad )
#endif

    end if ! pselect

    lma=0
    do a=1,Natom
       ik=ki_atom(a)
       do iorb=1,norb(ik)
          j=MJJ_tmp(iorb,a)
          if ( j > 0 ) then
             L=lo(iorb,ik)
             nzlma=nzlma+2*L+1
             do m=1,2*L+1
                lma=lma+1
                icheck_tmp3(a,iorb,m)=lma
             end do
          end if
       end do
    end do
    MMJJ = maxval( MJJ_tmp )

    allocate( lcheck_tmp1(Mlma,0:np_grid-1) )
    lcheck_tmp1(:,:)=.false.
    lma=0
    do a=1,Natom
       ik=ki_atom(a)
       do iorb=1,norb(ik)
          L=lo(iorb,ik)
          j=MJJ_tmp(iorb,a)
          do m=1,2*L+1
             lma=lma+1
             if ( j > 0 ) then
                lcheck_tmp1(lma,myrank_g)=.true.
             end if
          end do
       end do
    end do
    call mpi_allgather(lcheck_tmp1(1,myrank_g),Mlma,mpi_logical &
                      ,lcheck_tmp1,Mlma,mpi_logical,comm_grid,ierr)

    call watch(ctt(1),ett(1))

! for grid-parallel computation

    nzlma_0 = min(nzlma_0*2,Mlma)
    nrlma   = 0

    n=maxval( node_partition(1:3) )
    allocate( itmp(n,3) ) ; itmp=0
    allocate( lma_nsend_tmp(0:nprocs_g-1) )
    allocate( icheck_tmp1(0:nprocs_g-1)   )
    allocate( icheck_tmp2(0:nprocs_g-1)   )
    allocate( nl_rank_map_tmp(0:nprocs_g) )
    allocate( maps_tmp(nzlma_0,6) )
    allocate( sendmap_tmp(nzlma_0,0:nprocs_g-1) )
    allocate( recvmap_tmp(nzlma_0,0:nprocs_g-1) )

    maps_tmp(:,:)      = 0
    sendmap_tmp(:,:)   = 0
    recvmap_tmp(:,:)   = 0
    icheck_tmp1(:)     = 0
    icheck_tmp2(:)     = 0
    lma_nsend_tmp(:)   = 0
    nl_rank_map_tmp(:) =-1

    np1 = node_partition(1)
    np2 = node_partition(2)
    np3 = node_partition(3)

    lma=0
    do a=1,Natom
       ik=ki_atom(a)
    do iorb=1,norb(ik)
       L=lo(iorb,ik)
    do m=-L,L
       lma=lma+1
       icheck_tmp1(:)=0
       do n=0,np_grid-1
          if ( lcheck_tmp1(lma,n) ) icheck_tmp1(n) = 1
       end do
       icheck_tmp1(myrank_g) = icheck_tmp3(a,iorb,m+L+1)

       itmp(:,:)=0
       n=-1
       do i3=1,node_partition(3)
       do i2=1,node_partition(2)
       do i1=1,node_partition(1)
          n=n+1
          if ( icheck_tmp1(n) == 0 ) cycle
          itmp(i1,1) = i1
          itmp(i2,2) = i2
          itmp(i3,3) = i3
       end do
       end do
       end do
       k1=count( itmp(:,1)>0 )
       k2=count( itmp(:,2)>0 )
       k3=count( itmp(:,3)>0 )
       ic1=0
       id1=np1
       do i=1,np1
          if ( ic1==0 .and. itmp(i,1)/=0 ) then
             ic1=i
          else if ( ic1/=0 .and. itmp(i,1)==0 ) then
             id1=i-1
             exit
          end if
       end do
       if ( id1-ic1+1/=k1 ) then
          i1=0
          j1=np1
          do i=id1+1,np1
             if ( i1==0 .and. itmp(i,1)/=0 ) then
                i1=i
             else if ( i1/=0 .and. itmp(i,1)==0 ) then
                j1=i-1
                exit
             end if
          end do
          i1=i1-np1
          j1=j1-np1
          ic1=i1
       end if
       ic2=0
       id2=np2
       do i=1,np2
          if ( ic2==0 .and. itmp(i,2)/=0 ) then
             ic2=i
          else if ( ic2/=0 .and. itmp(i,2)==0 ) then
             id2=i-1
             exit
          end if
       end do
       if ( id2-ic2+1/=k2 ) then
          i2=0
          j2=np2
          do i=id2+1,np2
             if ( i2==0 .and. itmp(i,2)/=0 ) then
                i2=i
             else if ( i2/=0 .and. itmp(i,2)==0 ) then
                j2=i-1
                exit
             end if
          end do
          i2=i2-np2
          j2=j2-np2
          ic2=i2
       end if
       ic3=0
       id3=np3
       do i=1,np3
          if ( ic3==0 .and. itmp(i,3)/=0 ) then
             ic3=i
          else if ( ic3/=0 .and. itmp(i,3)==0 ) then
             id3=i-1
             exit
          end if
       end do
       if ( id3-ic3+1/=k3 ) then
          i3=0
          j3=np3
          do i=id3+1,np3
             if ( i3==0 .and. itmp(i,3)/=0 ) then
                i3=i
             else if ( i3/=0 .and. itmp(i,3)==0 ) then
                j3=i-1
                exit
             end if
          end do
          i3=i3-np3
          j3=j3-np3
          ic3=i3
       end if
       do j3=ic3,id3
       do j2=ic2,id2
       do j1=ic1,id1
          k1=mod(j1+np1-1,np1)+1
          k2=mod(j2+np2-1,np2)+1
          k3=mod(j3+np3-1,np3)+1
          k = k1-1 + (k2-1)*np1 + (k3-1)*np1*np2
          if ( icheck_tmp1(k)==0 ) icheck_tmp1(k)=-1
       end do
       end do
       end do
       do n=0,nprocs_g-1
          if ( icheck_tmp1(n)/=0 ) then
             icheck_tmp2(n)=icheck_tmp2(n)+1
          end if
       end do
       if ( icheck_tmp1(myrank_g)/=0 ) then
          if ( icheck_tmp1(myrank_g)>0 ) then
             maps_tmp(icheck_tmp2(myrank_g),1)=icheck_tmp1(myrank_g)
          end if
          maps_tmp(icheck_tmp2(myrank_g),2)=inorm(iorb,ik)
          maps_tmp(icheck_tmp2(myrank_g),3)=a
          maps_tmp(icheck_tmp2(myrank_g),4)=L
          maps_tmp(icheck_tmp2(myrank_g),5)=m
          maps_tmp(icheck_tmp2(myrank_g),6)=iorb

          do n=0,nprocs_g-1
             if ( n==myrank_g .or. icheck_tmp1(n)==0 ) cycle
             lma_nsend_tmp(n)=lma_nsend_tmp(n)+1  
             sendmap_tmp(lma_nsend_tmp(n),n)=icheck_tmp2(myrank_g)
             recvmap_tmp(lma_nsend_tmp(n),n)=icheck_tmp2(n)
             if ( any(nl_rank_map_tmp(0:nrlma)==n) ) cycle
             nrlma=nrlma+1
             nl_rank_map_tmp(nrlma)=n
          end do
       end if

    end do ! m
    end do ! iorb
    end do ! a

    call watch(ctt(2),ett(2))

    nzlma = icheck_tmp2(myrank_g)

    deallocate( itmp )
    deallocate( icheck_tmp2 )
    deallocate( icheck_tmp1 )
    deallocate( icheck_tmp3 )
    deallocate( lcheck_tmp1 )

    if ( allocated(uV) ) then
       deallocate( uV )
       deallocate( JJ_MAP )
       deallocate( MJJ_MAP )
       deallocate( MJJ )
       deallocate( iuV )
       deallocate( amap )
       deallocate( lmap )
       deallocate( mmap )
       deallocate( iorbmap )
       deallocate( nl_rank_map )
    end if
    allocate( uV(MMJJ,nzlma)       ) ; uV=0.d0
    allocate( JJ_MAP(6,MMJJ,nzlma) ) ; JJ_MAP=0
    allocate( MJJ_MAP(nzlma)       ) ; MJJ_MAP=0
    allocate( MJJ(nzlma)           ) ; MJJ=0
    allocate( iuV(nzlma)           ) ; iuV=0
    allocate( amap(nzlma)          ) ; amap=0
    allocate( lmap(nzlma)          ) ; lmap=0
    allocate( mmap(nzlma)          ) ; mmap=0
    allocate( iorbmap(nzlma)       ) ; iorbmap=0
    allocate( nl_rank_map(nrlma)   ) ; nl_rank_map=-1

    do i=1,nrlma
       nl_rank_map(i)=nl_rank_map_tmp(i)
    end do

    deallocate( nl_rank_map_tmp )

    do lma=1,nzlma
       if ( maps_tmp(lma,1) == 0 ) cycle
       iuV(lma)     = maps_tmp(lma,2)
       amap(lma)    = maps_tmp(lma,3)
       lmap(lma)    = maps_tmp(lma,4)
       mmap(lma)    = maps_tmp(lma,5)
       iorbmap(lma) = maps_tmp(lma,6)
    end do

    c1=1.d0/ML1
    c2=1.d0/ML2
    c3=1.d0/ML3

!$OMP parallel do private( a,l,m,iorb,Rx,Ry,Rz,j,i1,i2,i3,k1,k2,k3,d1,d2,d3,x,y,z )
    do lma=1,nzlma
       if ( maps_tmp(lma,1) == 0 ) cycle
       a    = amap(lma)
       l    = lmap(lma)
       m    = mmap(lma)
       iorb = iorbmap(lma)
       MJJ_MAP(lma) = MJJ_tmp(iorb,a)
       Rx=aa(1,1)*aa_atom(1,a)+aa(1,2)*aa_atom(2,a)+aa(1,3)*aa_atom(3,a)
       Ry=aa(2,1)*aa_atom(1,a)+aa(2,2)*aa_atom(2,a)+aa(2,3)*aa_atom(3,a)
       Rz=aa(3,1)*aa_atom(1,a)+aa(3,2)*aa_atom(2,a)+aa(3,3)*aa_atom(3,a)
       do j=1,MJJ_MAP(lma)
          i1=JJ_tmp(1,j,iorb,a)
          i2=JJ_tmp(2,j,iorb,a)
          i3=JJ_tmp(3,j,iorb,a)
          k1=JJ_tmp(4,j,iorb,a)
          k2=JJ_tmp(5,j,iorb,a)
          k3=JJ_tmp(6,j,iorb,a)
          d1=c1*i1+k1
          d2=c2*i2+k2
          d3=c3*i3+k3
          x = aa(1,1)*d1+aa(1,2)*d2+aa(1,3)*d3-Rx
          y = aa(2,1)*d1+aa(2,2)*d2+aa(2,3)*d3-Ry
          z = aa(3,1)*d1+aa(3,2)*d2+aa(3,3)*d3-Rz
          uV(j,lma) = uV_tmp(j,iorb,a)*Ylm(x,y,z,l,m)
          JJ_MAP(1:6,j,lma) = JJ_tmp(1:6,j,iorb,a)
       end do
    end do
!$OMP end parallel do

    deallocate( MJJ_tmp )
    deallocate( JJ_tmp )
    deallocate( uV_tmp )
    deallocate( maps_tmp )

    call watch(ctt(3),ett(3))

    allocate( icheck_tmp4(a1b:b1b,a2b:b2b,a3b:b3b) )
    icheck_tmp4=0
    do lma=1,nzlma
       j=0
       icheck_tmp4=0
       do i=1,MJJ_MAP(lma)
          i1=JJ_MAP(1,i,lma)
          i2=JJ_MAP(2,i,lma)
          i3=JJ_MAP(3,i,lma)
          if ( icheck_tmp4(i1,i2,i3)==0 ) then
             j=j+1
             icheck_tmp4(i1,i2,i3)=j
          end if
       end do
       MJJ(lma)=j
    end do
    MAXMJJ = maxval( MJJ(1:nzlma) )
    deallocate( icheck_tmp4 )

    nl_max_send = maxval( lma_nsend_tmp )

    if ( allocated(lma_nsend) ) then
       deallocate( lma_nsend )
       deallocate( sendmap )
       deallocate( recvmap )
    end if
    allocate( lma_nsend(0:nprocs_g-1) ) ; lma_nsend=0
    allocate( sendmap(nl_max_send,0:nprocs_g-1) ) ; sendmap=0
    allocate( recvmap(nl_max_send,0:nprocs_g-1) ) ; recvmap=0

    do n=0,nprocs_g-1
       sendmap(1:nl_max_send,n) = sendmap_tmp(1:nl_max_send,n)
       lma_nsend(n) = lma_nsend_tmp(n)
    end do


    allocate( ireq(2*nprocs_g) )
    allocate( istatus(MPI_STATUS_SIZE,2*nprocs_g) )
    nreq=0
    do n=0,nprocs_g-1
       if ( lma_nsend(n)<=0 .or. n==myrank_g ) cycle
       nreq=nreq+1
       call mpi_isend(recvmap_tmp(1,n),lma_nsend(n),mpi_integer,n,1 &
            ,comm_grid,ireq(nreq),ierr)
       nreq=nreq+1
       call mpi_irecv(recvmap(1,n) ,lma_nsend(n),mpi_integer,n,1 &
            ,comm_grid,ireq(nreq),ierr)
    end do
    call mpi_waitall(nreq,ireq,istatus,ierr)
    deallocate( istatus )
    deallocate( ireq )

    deallocate( recvmap_tmp,sendmap_tmp,lma_nsend_tmp )

    allocate( LLp(3,0:nprocs_g-1) )
    n=-1
    do i3=0,node_partition(3)-1
    do i2=0,node_partition(2)-1
    do i1=0,node_partition(1)-1
       n=n+1
       LLp(1,n)=i1
       LLp(2,n)=i2
       LLp(3,n)=i3
    end do
    end do
    end do

    allocate( itmp(3,nrlma) ) ; itmp=0
    allocate( itmp1(nrlma), work(nrlma) )
    allocate( itmp2(nrlma),itmp3(3,nrlma) )

    do irlma=1,nrlma
       n=nl_rank_map(irlma)
       itmp(1,irlma)=LLp(1,n)-LLp(1,myrank_g)
       itmp(2,irlma)=LLp(2,n)-LLp(2,myrank_g)
       itmp(3,irlma)=LLp(3,n)-LLp(3,myrank_g)
    end do

    nrlma_xyz(1:6)=0

    m=0
    n=0
    do i=1,nrlma
       if( itmp(2,i)==0 .and. itmp(3,i)==0 .and. itmp(1,i)>0 )then
          n=n+1
          work(n)=itmp(1,i)
          itmp2(n)=i
       end if
    end do
    if ( n>0 ) then
       call indexx(n,work,itmp1)
       do i=1,n
          j=itmp2( itmp1(i) )
          itmp3(:,m+i)=itmp(:,j)
       end do
    end if
    m=m+n
    nrlma_xyz(1)=nrlma_xyz(1)+n
    n=0
    do i=1,nrlma
       if( itmp(2,i)==0 .and. itmp(3,i)==0 .and. itmp(1,i)<0 )then
          n=n+1
          work(n)=itmp(1,i)
          itmp2(n)=i
       end if
    end do
    if ( n>0 ) then
       call indexx(n,work,itmp1)
       do i=1,n
          j=itmp2(itmp1(i))
          itmp3(:,m+n-i+1)=itmp(:,j)
       end do
    end if
    m=m+n
    nrlma_xyz(2)=nrlma_xyz(2)+n

    n=0
    do i=1,nrlma
       if( itmp(1,i)==0 .and. itmp(3,i)==0 .and. itmp(2,i)>0 )then
          n=n+1
          work(n)=itmp(2,i)
          itmp2(n)=i
       end if
    end do
    if ( n>0 ) then
       call indexx(n,work,itmp1)
       do i=1,n
          j=itmp2( itmp1(i) )
          itmp3(:,m+i)=itmp(:,j)
       end do
    end if
    m=m+n
    nrlma_xyz(3)=nrlma_xyz(3)+n
    n=0
    do i=1,nrlma
       if( itmp(1,i)==0 .and. itmp(3,i)==0 .and. itmp(2,i)<0 )then
          n=n+1
          work(n)=itmp(2,i)
          itmp2(n)=i
       end if
    end do
    if ( n>0 ) then
       call indexx(n,work,itmp1)
       do i=1,n
          j=itmp2(itmp1(i))
          itmp3(:,m+n-i+1)=itmp(:,j)
       end do
    end if
    m=m+n
    nrlma_xyz(4)=nrlma_xyz(4)+n

    n=0
    do i=1,nrlma
       if( itmp(1,i)==0 .and. itmp(2,i)==0 .and. itmp(3,i)>0 )then
          n=n+1
          work(n)=itmp(3,i)
          itmp2(n)=i
       end if
    end do
    if ( n>0 ) then
       call indexx(n,work,itmp1)
       do i=1,n
          j=itmp2( itmp1(i) )
          itmp3(:,m+i)=itmp(:,j)
       end do
    end if
    m=m+n
    nrlma_xyz(5)=nrlma_xyz(5)+n
    n=0
    do i=1,nrlma
       if( itmp(1,i)==0 .and. itmp(2,i)==0 .and. itmp(3,i)<0 )then
          n=n+1
          work(n)=itmp(3,i)
          itmp2(n)=i
       end if
    end do
    if ( n>0 ) then
       call indexx(n,work,itmp1)
       do i=1,n
          j=itmp2(itmp1(i))
          itmp3(:,m+n-i+1)=itmp(:,j)
       end do
    end if
    m=m+n
    nrlma_xyz(6)=nrlma_xyz(6)+n


    n=maxval( nrlma_xyz )
    if ( allocated(num_2_rank) ) then
       deallocate( num_2_rank )
    end if
    allocate( num_2_rank(n,6) )
    num_2_rank(:,:)=MPI_PROC_NULL


    m=0
    do i=1,6
       do j=1,nrlma_xyz(i)
          m=m+1
          i1=itmp3(1,m)+LLp(1,myrank_g)
          i2=itmp3(2,m)+LLp(2,myrank_g)
          i3=itmp3(3,m)+LLp(3,myrank_g)
          k = i1 + i2*np1 + i3*np1*np2
          num_2_rank(j,i)=k
       end do
    end do

    deallocate( itmp,itmp1,itmp2,itmp3,work )
    deallocate( LLp )

    do i=1,5,2
       n=max( nrlma_xyz(i),nrlma_xyz(i+1) )
       nrlma_xyz(i)=n
       nrlma_xyz(i+1)=n
    end do

    call watch(ctt(4),ett(4))

    call allocate_ps_nloc2(MB_d)

    if ( allocated(JJP) ) then
       deallocate( JJP )
       deallocate( uVk )
    end if
    allocate( JJP(MAXMJJ,nzlma) ) ; JJP=0
    allocate( uVk(MAXMJJ,nzlma,MBZ_0:MBZ_1) ) ; uVk=0.d0

    call prep_uvk_ps_nloc_mr(MBZ_0,MBZ_1,kbb(1,MBZ_0))

    call watch(ctt(5),ett(5))

!--- multi-reference cofficients

    N_nzqr=0
    do lma1=1,nzlma
       a1 = amap(lma1) ; if( a1==0 ) cycle
       l1 = lmap(lma1)
       m1 = mmap(lma1)
       i1 = no( iorbmap(lma1), ki_atom(a1) )
    do lma2=1,nzlma
       a2 = amap(lma2) ; if ( a2==0 ) cycle
       l2 = lmap(lma2)
       m2 = mmap(lma2)
       i2 = no( iorbmap(lma2), ki_atom(a2) )
       if ( a1==a2 .and. l1==l2 .and. m1==m2 ) then
          N_nzqr=N_nzqr+1
       end if
    end do ! lma2
    end do ! lma1

    if ( disp_switch_parallel ) write(*,*) "N_nzqr=",N_nzqr

    if ( allocated(Dij00)     ) deallocate(Dij00)
    if ( allocated(nzqr_pair) ) deallocate(nzqr_pair)
    allocate( Dij00(N_nzqr)       ) ; Dij00=0.0d0
    allocate( nzqr_pair(2,N_nzqr) ) ; nzqr_pair=0

    n=0
    do lma1=1,nzlma
       a1 = amap(lma1) ; if ( a1==0 ) cycle
       l1 = lmap(lma1)
       m1 = mmap(lma1)
       i1 = no( iorbmap(lma1), ki_atom(a1) )
    do lma2=1,nzlma
       a2 = amap(lma2) ; if ( a2==0 ) cycle
       l2 = lmap(lma2)
       m2 = mmap(lma2)
       i2 = no( iorbmap(lma2), ki_atom(a2) )
       if ( a1==a2 .and. l1==l2 .and. m1==m2 ) then
          n=n+1
          nzqr_pair(1,n)=lma1
          nzqr_pair(2,n)=lma2
          Dij00(n)=hnml(i1,i2,l1,ki_atom(a1))
       end if
    end do ! lma2
    end do ! lma1


    if ( disp_switch_parallel ) then
       write(*,*) "time(ps_nloc_mr_1)",ctt(1)-ctt(0),ett(1)-ett(0)
       write(*,*) "time(ps_nloc_mr_2)",ctt(2)-ctt(1),ett(2)-ett(1)
       write(*,*) "time(ps_nloc_mr_3)",ctt(3)-ctt(2),ett(3)-ett(2)
       write(*,*) "time(ps_nloc_mr_4)",ctt(4)-ctt(3),ett(4)-ett(3)
       write(*,*) "time(ps_nloc_mr_5)",ctt(5)-ctt(4),ett(5)-ett(4)
    end if

  END SUBROUTINE prep_ps_nloc_mr


  SUBROUTINE prep_uvk_ps_nloc_mr(k0,k1,kbb)
    implicit none
    integer,intent(IN) :: k0,k1
    real(8),intent(IN) :: kbb(3,k0:k1)
    integer :: a1b,b1b,a2b,b2b,a3b,b3b,ab1,ab2,ab3
    integer :: i,j,k,j3,lma,i0,i1,i2,i3,m1,m2,m3
    integer,allocatable :: icheck_tmp4(:,:,:)
    real(8) :: c1,c2,c3,d1,d2,d3,pi2,kr
    complex(8) :: ztmp0

    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)

    c1=1.d0/Ngrid(1)
    c2=1.d0/Ngrid(2)
    c3=1.d0/Ngrid(3)

    ab1=b1b-a1b+1
    ab2=b2b-a2b+1
    ab3=b3b-a3b+1

    allocate( icheck_tmp4(a1b:b1b,a2b:b2b,a3b:b3b) )
    icheck_tmp4=0

    pi2 = 2.d0*acos(-1.d0)

    do k=k0,k1
       d1=pi2*kbb(1,k)
       d2=pi2*kbb(2,k)
       d3=pi2*kbb(3,k)
       do lma=1,nzlma
          j=0
          icheck_tmp4=0
          do i=1,MJJ_MAP(lma)
             i1=JJ_MAP(1,i,lma)
             i2=JJ_MAP(2,i,lma)
             i3=JJ_MAP(3,i,lma)
             m1=JJ_MAP(4,i,lma)
             m2=JJ_MAP(5,i,lma)
             m3=JJ_MAP(6,i,lma)
             j3=icheck_tmp4(i1,i2,i3)
             kr=d1*(c1*i1+m1)+d2*(c2*i2+m2)+d3*(c3*i3+m3)
             ztmp0=dcmplx(cos(kr),-sin(kr))*uV(i,lma)
             if ( j3==0 ) then
                j=j+1
                icheck_tmp4(i1,i2,i3)=j
                uVk(j,lma,k)=ztmp0
                JJP(j,lma) = i1-a1b + (i2-a2b)*ab1 + (i3-a3b)*ab1*ab2 + ML_0
             else
                uVk(j3,lma,k)=uVk(j3,lma,k)+ztmp0
             end if
          end do
       end do ! lma
    end do ! k

    deallocate( icheck_tmp4 )

  END SUBROUTINE prep_uvk_ps_nloc_mr


  SUBROUTINE prep_rvk_ps_nloc_mr(k0,k1,kbb)
    implicit none
    integer,intent(IN) :: k0,k1
    real(8),intent(IN) :: kbb(3,k0:k1)
    integer :: a1b,b1b,a2b,b2b,a3b,b3b,ab1,ab2,ab3
    integer :: i,j,k,j3,lma,i1,i2,i3,m1,m2,m3
    integer,allocatable :: icheck_tmp4(:,:,:)
    real(8) :: c1,c2,c3,d1,d2,d3,pi2,kr,x,y,z
    complex(8) :: ztmp0
    integer,save :: array_size(4)=0

    if ( array_size(1) /= MAXMJJ .or. array_size(2) /= nzlma .or. &
         array_size(3) /= MBZ_0  .or. array_size(4) /= MBZ_1 ) then
       if ( allocated(xVk) ) deallocate(xVk)
       if ( allocated(yVk) ) deallocate(yVk)
       if ( allocated(zVk) ) deallocate(zVk)
       allocate( xVk(MAXMJJ,nzlma,MBZ_0:MBZ_1) )
       allocate( yVk(MAXMJJ,nzlma,MBZ_0:MBZ_1) )
       allocate( zVk(MAXMJJ,nzlma,MBZ_0:MBZ_1) )
       array_size(1)=MAXMJJ
       array_size(2)=nzlma
       array_size(3)=MBZ_0
       array_size(4)=MBZ_1
    end if
    xVk=0.d0
    yVk=0.d0
    zVk=0.d0

    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)

    c1=1.d0/Ngrid(1)
    c2=1.d0/Ngrid(2)
    c3=1.d0/Ngrid(3)

    ab1=b1b-a1b+1
    ab2=b2b-a2b+1
    ab3=b3b-a3b+1

    allocate( icheck_tmp4(a1b:b1b,a2b:b2b,a3b:b3b) )
    icheck_tmp4=0

    pi2 = 2.d0*acos(-1.d0)

    do k=k0,k1
       d1=pi2*kbb(1,k)
       d2=pi2*kbb(2,k)
       d3=pi2*kbb(3,k)
       do lma=1,nzlma
          j=0
          icheck_tmp4=0
          do i=1,MJJ_MAP(lma)
             i1=JJ_MAP(1,i,lma)
             i2=JJ_MAP(2,i,lma)
             i3=JJ_MAP(3,i,lma)
             m1=JJ_MAP(4,i,lma)
             m2=JJ_MAP(5,i,lma)
             m3=JJ_MAP(6,i,lma)
             j3=icheck_tmp4(i1,i2,i3)
             kr=d1*(c1*i1+m1)+d2*(c2*i2+m2)+d3*(c3*i3+m3)
             ztmp0=dcmplx(cos(kr),-sin(kr))*uV(i,lma)
             x=(c1*i1+m1)*aa(1,1)+(c2*i2+m2)*aa(1,2)+(c3*i3+m3)*aa(1,3)
             y=(c1*i1+m1)*aa(2,1)+(c2*i2+m2)*aa(2,2)+(c3*i3+m3)*aa(2,3)
             z=(c1*i1+m1)*aa(3,1)+(c2*i2+m2)*aa(3,2)+(c3*i3+m3)*aa(3,3)
             if ( j3==0 ) then
                j=j+1
                icheck_tmp4(i1,i2,i3)=j
                xVk(j,lma,k)=x*ztmp0
                yVk(j,lma,k)=y*ztmp0
                zVk(j,lma,k)=z*ztmp0
             else
                xVk(j3,lma,k)=xVk(j3,lma,k)+x*ztmp0
                yVk(j3,lma,k)=yVk(j3,lma,k)+y*ztmp0
                zVk(j3,lma,k)=zVk(j3,lma,k)+z*ztmp0
             end if
          end do
       end do ! lma
    end do ! k

    deallocate( icheck_tmp4 )

  END SUBROUTINE prep_rvk_ps_nloc_mr


  SUBROUTINE op_ps_nloc_mr(k,tpsi,htpsi,n1,n2,ib1,ib2)
    implicit none
    integer,intent(IN) :: k,n1,n2,ib1,ib2
#ifdef _DRSDFT_
    real(8),intent(IN)  :: tpsi(n1:n2,ib1:ib2)
    real(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
    real(8),allocatable :: uVunk(:,:),uVunk0(:,:)
#else
    complex(8),intent(IN)  :: tpsi(n1:n2,ib1:ib2)
    complex(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
    complex(8),allocatable :: uVunk(:,:),uVunk0(:,:)
#endif
    integer :: i,ib,j,i1,i2,m,lma,nb,ierr,nreq,lma1,lma2
    integer :: irank,jrank,istatus(mpi_status_size,512),ireq(512)
    complex(8) :: zc

    nb = ib2-ib1+1

    if ( Mlma <= 0 ) return

    allocate( uVunk(nzlma,ib1:ib2),uVunk0(nzlma,ib1:ib2) )

!$OMP parallel

    do ib=ib1,ib2
!$OMP do
       do lma=1,nzlma
          uVunk(lma,ib)=zero
          do j=1,MJJ(lma)
             i=JJP(j,lma)
#ifdef _DRSDFT_
             uVunk(lma,ib)=uVunk(lma,ib)+uVk(j,lma,k)*tpsi(i,ib)
#else
             uVunk(lma,ib)=uVunk(lma,ib)+conjg(uVk(j,lma,k))*tpsi(i,ib)
#endif
          end do
          uVunk(lma,ib)=uVunk(lma,ib)*dV
       end do
!$OMP end do
    end do

!$OMP single
    do i=1,6
       select case(i)
       case(1,3,5)
!!$OMP single
          j=i+1
!!$OMP end single
!!$OMP workshare
          uVunk0(:,:)=uVunk(:,:)
!!$OMP end workshare
       case(2,4,6)
!!$OMP single
          j=i-1
!!$OMP end single
       end select
!!$OMP single
       do m=1,nrlma_xyz(i)
          nreq=0
          irank=num_2_rank(m,i)
          jrank=num_2_rank(m,j)
          if( irank>=0 )then
             i2=0
             do ib=ib1,ib2
                do i1=1,lma_nsend(irank)
                   i2=i2+1
                   sbufnl(i2,irank)=uVunk0(sendmap(i1,irank),ib)
                end do
             end do
             nreq=nreq+1
             call mpi_isend(sbufnl(1,irank),lma_nsend(irank)*nb &
                  ,TYPE_MAIN,irank,1,comm_grid,ireq(nreq),ierr)
          end if
          if( jrank>=0 )then
             nreq=nreq+1
             call mpi_irecv(rbufnl(1,jrank),lma_nsend(jrank)*nb &
                  ,TYPE_MAIN,jrank,1,comm_grid,ireq(nreq),ierr)
          end if
          call mpi_waitall(nreq,ireq,istatus,ierr)
          if( jrank>=0 )then
             i2=0
             do ib=ib1,ib2
                do i1=1,lma_nsend(jrank)
                   i2=i2+1
                   uVunk(recvmap(i1,jrank),ib) &
                        =uVunk(recvmap(i1,jrank),ib)+rbufnl(i2,jrank)
                end do
             end do
          end if
       end do
!!$OMP end single
    end do
!$OMP end single

    do ib=ib1,ib2
       do m=1,N_nzqr
          lma1=nzqr_pair(1,m)
          lma2=nzqr_pair(2,m)
!$OMP do
          do j=1,MJJ(lma1)
             htpsi(JJP(j,lma1),ib)=htpsi(JJP(j,lma1),ib) &
                  +Dij00(m)*uVk(j,lma1,k)*uVunk(lma2,ib)
          end do
!$OMP end do
       end do
    end do

!$OMP end parallel

    deallocate( uVunk0,uVunk )

  END SUBROUTINE op_ps_nloc_mr


  SUBROUTINE calc_force_ps_nloc_mr(MI,force2)
    implicit none
    integer,intent(IN) :: MI
    real(8),intent(OUT) :: force2(3,MI)
    integer :: i1,i2,i3,lma1,lma2
    integer :: i,j,k,s,n,ir,iorb,L,L1,L1z,NRc,irank,jrank
    integer :: nreq,max_nreq
    integer :: a,a0,ik,m,lm0,lm1,lma,im,m1,m2
    integer :: ierr,M_irad,ir0
    integer,allocatable :: ireq(:),istatus(:,:),irad(:,:),ilm1(:,:,:)
    real(8),parameter :: ep=1.d-8
    real(8),save :: Y1(0:3,-3:3,0:4,-4:4)
    real(8),save :: Y2(0:3,-3:3,0:4,-4:4)
    real(8),save :: Y3(0:3,-3:3,0:4,-4:4)
    real(8) :: err,err0,maxerr,Rx,Ry,Rz
    real(8) :: a1,a2,a3,c1,c2,c3,d1,d2,d3
    real(8) :: x,y,z,r,kr,pi2,c
    real(8) :: tmp,tmp0,tmp1
    real(8) :: ctt(0:3),ett(0:3)
    real(8) :: yy1,yy2,yy3
    real(8),allocatable :: work2(:,:),duVdR(:,:,:)
#ifdef _DRSDFT_
    real(8) :: ztmp
    real(8),allocatable :: wtmp5(:,:,:,:,:),vtmp2(:,:)
#else
    complex(8) :: ztmp
    complex(8),allocatable :: wtmp5(:,:,:,:,:),vtmp2(:,:)
#endif
    logical,save :: flag_Y = .true.
    logical,allocatable :: a_rank(:)
    integer :: ML1,ML2,ML3,i0,iorb0
    integer :: k1,k2,k3,a1b,a2b,a3b,ab1,ab2,ab3

    force2(:,:) = 0.0d0

    if ( Mlma <= 0 ) return

    if ( flag_Y ) then
       Y1=0.0d0
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
       Y2=0.0d0
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
       Y3=0.0d0
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
       flag_Y = .false.
    end if

    pi2 = 2.d0*acos(-1.d0)

    maxerr=0.d0
    ctt(:)=0.d0
    ett(:)=0.d0

    a1b=Igrid(1,1)
    a2b=Igrid(1,2)
    a3b=Igrid(1,3)
    ab1=Igrid(2,1)-Igrid(1,1)+1
    ab2=Igrid(2,2)-Igrid(1,2)+1
    ab3=Igrid(2,3)-Igrid(1,3)+1

    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    c1=1.d0/ML1
    c2=1.d0/ML2
    c3=1.d0/ML3

    if ( .not.allocated(ilm1) ) then
       L1=maxval(lo)+1
       n=maxval(norb)
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
    end if

    if ( .not.allocated(y2b) .and. all(ippform /= 4) ) then
       lm1=maxval(ilm1)
       NRc=maxval(NRps)
       allocate( y2b(NRc,lm1,Nelement) )
       y2b=0.d0
       do ik=1,Nelement
       do iorb=1,norb(ik)
          L=lo(iorb,ik)
          do L1=abs(L-1),L+1
             lm1=ilm1(L1,iorb,ik)
             d1=0.d0
             d2=0.d0
             call spline(rad1(1,ik),dviod(1,lm1,ik),NRps(iorb,ik),d1,d2,y2b(1,lm1,ik))
          end do
       end do
       end do
    end if

    allocate( wtmp5(0:3,nzlma,MB_0:MB_1,MBZ_0:MBZ_1,MSP_0:MSP_1) )
    allocate( vtmp2(0:3,nzlma) )
    allocate( a_rank(Natom) )
    allocate( duVdR(3,MMJJ,nzlma) )

!$OMP parallel

!$OMP workshare
    wtmp5=zero
    a_rank(:)=.false.
!$OMP end workshare

!$OMP do private( i1,i2,i3,k1,k2,k3 )
    do a=1,Natom
       i1 = nint( aa_atom(1,a)*ML1 )
       i2 = nint( aa_atom(2,a)*ML2 )
       i3 = nint( aa_atom(3,a)*ML3 )
       k1 = i1/ML1 ; if ( i1<0 ) k1=(i1+1)/ML1-1
       k2 = i2/ML2 ; if ( i2<0 ) k2=(i2+1)/ML2-1
       k3 = i3/ML3 ; if ( i3<0 ) k3=(i3+1)/ML3-1
       i1 = i1 - k1*ML1
       i2 = i2 - k2*ML2
       i3 = i3 - k3*ML3
       if ( Igrid(1,1) <= i1 .and. i1 <= Igrid(2,1) .and. &
            Igrid(1,2) <= i2 .and. i2 <= Igrid(2,2) .and. &
            Igrid(1,3) <= i3 .and. i3 <= Igrid(2,3) ) then
          a_rank(a)=.true.
       end if
    end do
!$OMP end do

    if ( any(ippform == 4) ) then

!$OMP master
       call watch(ctt(0),ett(0))
!$OMP end master

!$OMP single
       call init_force_ps_nloc_hgh &
            (MMJJ,nzlma,amap,lmap,mmap,iorbmap,MJJ_MAP,JJ_MAP,Y1,Y2,Y3,duVdR)
!$OMP end single

!$OMP master
       call watch(ctt(1),ett(1))
!$OMP end master

    else

#ifndef _SPLINE_

!$OMP single
       allocate( irad(0:3000,Nelement) )
       irad=0
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
!$OMP end single

#endif

!$OMP workshare
       duVdR=0.d0
!$OMP end workshare

!$OMP master
       call watch(ctt(0),ett(0))
!$OMP end master

!$OMP do schedule(dynamic) firstprivate( maxerr ) &
!$OMP    private( a,L,m,iorb,ik,Rx,Ry,Rz,NRc,d1,d2,d3,x,y,z,r  &
!$OMP            ,ir,ir0,yy1,yy2,yy3,err0,err,tmp0,tmp1,m1,m2  &
!$OMP            ,lma,j,L1,L1z,lm1,im )
       do lma=1,nzlma
          a    = amap(lma)
          if ( a <= 0 ) cycle
          L    = lmap(lma)
          m    = mmap(lma)
          iorb = iorbmap(lma)
          ik   = ki_atom(a)
          Rx=aa(1,1)*aa_atom(1,a)+aa(1,2)*aa_atom(2,a)+aa(1,3)*aa_atom(3,a)
          Ry=aa(2,1)*aa_atom(1,a)+aa(2,2)*aa_atom(2,a)+aa(2,3)*aa_atom(3,a)
          Rz=aa(3,1)*aa_atom(1,a)+aa(3,2)*aa_atom(2,a)+aa(3,3)*aa_atom(3,a)
          NRc=NRps(iorb,ik)
!!$OMP parallel do firstprivate( maxerr ) &
!!$OMP             private( d1,d2,d3,x,y,z,r,ir0,yy1,yy2,yy3,lm1,err,err0 &
!!$OMP                     ,tmp0,tmp1,m1,m2,j,L1,im,L1z )
          do j=1,MJJ_MAP(lma)
             d1=c1*JJ_MAP(1,j,lma)+JJ_MAP(4,j,lma)
             d2=c2*JJ_MAP(2,j,lma)+JJ_MAP(5,j,lma)
             d3=c3*JJ_MAP(3,j,lma)+JJ_MAP(6,j,lma)
             x = aa(1,1)*d1+aa(1,2)*d2+aa(1,3)*d3-Rx
             y = aa(2,1)*d1+aa(2,2)*d2+aa(2,3)*d3-Ry
             z = aa(3,1)*d1+aa(3,2)*d2+aa(3,3)*d3-Rz
             r = sqrt(x*x+y*y+z*z)
#ifndef _SPLINE_
             ir0=irad( int(100.d0*r),ik )
             do ir=ir0,NRc
                if ( r<rad1(ir,ik) ) exit
             end do
#endif
             yy1=0.d0
             yy2=0.d0
             yy3=0.d0
             do L1=abs(L-1),L+1
                lm1=ilm1(L1,iorb,ik)
                if ( abs(x)>1.d-14 .or. abs(y)>1.d-14 &
                .or. abs(z)>1.d-14 .or. L1==0 ) then
#ifdef _SPLINE_
                   if ( r < rad1(2,ik) ) then
                      tmp0=dviod(2,lm1,ik)/(rad1(2,ik)**2)
                   else
                      call splint(rad1(1,ik),dviod(1,lm1,ik),y2b(1,lm1,ik),NRc,r,tmp0)
                      tmp0=tmp0/(r*r)
                   end if
#else
                   if ( ir <= 2 ) then
                      err0=0.d0
                      tmp0=dviod(2,lm1,ik)/(rad1(2,ik)**2)
                      if ( ir < 1 ) stop "calc_force_ps_nloc_mr"
                   else if ( ir <= NRc ) then
                      err0=1.d10
                      do im=1,20
                         m1=max(1,ir-im)
                         m2=min(ir+im,NRc)
                         call polint(rad1(m1,ik),dviod(m1,lm1,ik) &
                              ,m2-m1+1,r,tmp1,err)
                         if ( abs(err)<err0 ) then
                            tmp0=tmp1
                            err0=abs(err)
                            if ( err0<ep ) exit
                         end if
                      end do
                      tmp0=tmp0/(r*r)
                   else
                      write(*,*) "force_ps_nloc_mr",ir,NRc
                      stop
                   end if
                   maxerr=max(maxerr,err0)
#endif
                   do L1z=-L1,L1
                      tmp1=tmp0*Ylm(x,y,z,L1,L1z)
                      yy1=yy1+tmp1*Y1(L,m,L1,L1z)
                      yy2=yy2+tmp1*Y2(L,m,L1,L1z)
                      yy3=yy3-tmp1*Y3(L,m,L1,L1z)
                   end do
                end if
             end do ! L1

             duVdR(1,j,lma)=yy1
             duVdR(2,j,lma)=yy2
             duVdR(3,j,lma)=yy3

          end do ! j
!!$OMP end parallel do
       end do ! lma
!$OMP end do

!$OMP master
       call watch(ctt(1),ett(1))
!$OMP end master

#ifndef _SPLINE_
!$OMP single
       deallocate( irad )
!$OMP end single
#endif

    end if

    do s=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1
!$OMP do schedule(dynamic) private( c,i,d1,d2,d3,kr,ztmp,i1,i2,i3 )
    do n=MB_0,MB_1
       if ( occ(n,k,s) == 0.d0 ) cycle
       c=-2.d0*occ(n,k,s)*dV*dV
       do lma=1,nzlma
          if ( MJJ_MAP(lma) == MJJ(lma) ) then
#ifdef _DRSDFT_
             do j=1,MJJ(lma)
                i=JJP(j,lma)
                wtmp5(0,lma,n,k,s)=wtmp5(0,lma,n,k,s)+uVk(j,lma,k)*unk(i,n,k,s)
                wtmp5(1,lma,n,k,s)=wtmp5(1,lma,n,k,s)+duVdR(1,j,lma)*unk(i,n,k,s)
                wtmp5(2,lma,n,k,s)=wtmp5(2,lma,n,k,s)+duVdR(2,j,lma)*unk(i,n,k,s)
                wtmp5(3,lma,n,k,s)=wtmp5(3,lma,n,k,s)+duVdR(3,j,lma)*unk(i,n,k,s)
             end do
             wtmp5(0,lma,n,k,s)=c*wtmp5(0,lma,n,k,s)
#else
             do j=1,MJJ(lma)
                i=JJP(j,lma)
                wtmp5(0,lma,n,k,s)=wtmp5(0,lma,n,k,s) &
                     +uVk(j,lma,k)*conjg(unk(i,n,k,s))
             end do
             wtmp5(0,lma,n,k,s)=c*wtmp5(0,lma,n,k,s)
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
                wtmp5(0,lma,n,k,s)=wtmp5(0,lma,n,k,s)+uVk(j,lma,k)*unk(i,n,k,s)
             end do
             wtmp5(0,lma,n,k,s)=c*wtmp5(0,lma,n,k,s)
             do j=1,MJJ_MAP(lma)
                i1=JJ_MAP(1,j,lma)
                i2=JJ_MAP(2,j,lma)
                i3=JJ_MAP(3,j,lma)
                i = i1-a1b + (i2-a2b)*ab1 + (i3-a3b)*ab1*ab2 + ML_0
                wtmp5(1,lma,n,k,s)=wtmp5(1,lma,n,k,s)+duVdR(1,j,lma)*unk(i,n,k,s)
                wtmp5(2,lma,n,k,s)=wtmp5(2,lma,n,k,s)+duVdR(2,j,lma)*unk(i,n,k,s)
                wtmp5(3,lma,n,k,s)=wtmp5(3,lma,n,k,s)+duVdR(3,j,lma)*unk(i,n,k,s)
             end do
#else
             do j=1,MJJ(lma)
                i=JJP(j,lma)
                wtmp5(0,lma,n,k,s)=wtmp5(0,lma,n,k,s) &
                     +uVk(j,lma,k)*conjg(unk(i,n,k,s))
             end do
             wtmp5(0,lma,n,k,s)=c*wtmp5(0,lma,n,k,s)
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
    max_nreq=2*maxval( nrlma_xyz )
    allocate( ireq(max_nreq) )
    allocate( istatus(MPI_STATUS_SIZE,max_nreq) )
!$OMP end single

!$OMP workshare
    force2(:,:)=0.d0
!$OMP end workshare

!$OMP single
    do s=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1
    do n=MB_0,MB_1

       if ( occ(n,k,s) == 0.d0 ) cycle

       do i=1,6
          select case(i)
          case(1,3,5)
             j=i+1
             vtmp2(:,:)=wtmp5(:,:,n,k,s)
          case(2,4,6)
             j=i-1
          end select
          do m=1,nrlma_xyz(i)
             nreq=0
             irank=num_2_rank(m,i)
             jrank=num_2_rank(m,j)
             if( irank >= 0 )then
                i1=0
                do i2=1,lma_nsend(irank)
                do i3=0,3
                   i1=i1+1
                   sbufnl(i1,irank)=vtmp2(i3,sendmap(i2,irank))
                end do
                end do
                nreq=nreq+1
                call mpi_isend(sbufnl(1,irank),4*lma_nsend(irank) &
                     ,TYPE_MAIN,irank,1,comm_grid,ireq(nreq),ierr)
             end if
             if( jrank >= 0 )then
                nreq=nreq+1
                call mpi_irecv(rbufnl(1,jrank),4*lma_nsend(jrank) &
                     ,TYPE_MAIN,jrank,1,comm_grid,ireq(nreq),ierr)
             end if
             call mpi_waitall(nreq,ireq,istatus,ierr)
             if( jrank >= 0 )then
                i1=0
                do i2=1,lma_nsend(jrank)
                do i3=0,3
                   i1=i1+1
                   wtmp5(i3,recvmap(i2,jrank),n,k,s) &
                        =wtmp5(i3,recvmap(i2,jrank),n,k,s)+rbufnl(i1,jrank)
                end do
                end do
             end if
          end do ! m
       end do ! i

       do m=1,N_nzqr
          lma1=nzqr_pair(1,m)
          lma2=nzqr_pair(2,m)
          a=amap(lma1)
          if ( a <= 0 ) cycle
          if ( a_rank(a) ) then
             force2(1,a)=force2(1,a) &
                  +Dij00(m)*real(wtmp5(0,lma1,n,k,s)*wtmp5(1,lma2,n,k,s),8)
             force2(2,a)=force2(2,a) &
                  +Dij00(m)*real(wtmp5(0,lma1,n,k,s)*wtmp5(2,lma2,n,k,s),8)
             force2(3,a)=force2(3,a) &
                  +Dij00(m)*real(wtmp5(0,lma1,n,k,s)*wtmp5(3,lma2,n,k,s),8)
          end if
       end do

    end do ! n
    end do ! k
    end do ! s

    deallocate( istatus )
    deallocate( ireq )

    allocate( work2(3,Natom) )
!$OMP end single

!$OMP workshare
    work2(:,:)=force2(:,:)
!$OMP end workshare

!$OMP single
    call mpi_allreduce(work2,force2,3*Natom &
         ,mpi_real8,mpi_sum,mpi_comm_world,ierr)
    deallocate( work2 )
    deallocate( a_rank )
    deallocate( vtmp2 )
    deallocate( wtmp5 )
!$OMP end single

!$OMP master
    call watch(ctt(3),ett(3))
!$OMP end master

!$OMP end parallel

    if ( disp_switch_parallel ) then
       write(*,*) "time(force_nloc_mr_1)",ctt(1)-ctt(0),ett(1)-ett(0)
       write(*,*) "time(force_nloc_mr_2)",ctt(2)-ctt(1),ett(2)-ett(1)
       write(*,*) "time(force_nloc_mr_3)",ctt(3)-ctt(2),ett(3)-ett(2)
    end if

  END SUBROUTINE calc_force_ps_nloc_mr


END MODULE ps_nloc_mr_module
