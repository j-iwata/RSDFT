MODULE ps_nloc2_mol_module

  use rgrid_module
  use rgrid_mol_module
  use atom_module
  use array_bound_module, only: ML_0,ML_1
  use parallel_module
  use pseudopot_module
  use ps_nloc2_init_module
  use ps_nloc2_variables
  use minimal_box_module
  use ps_nloc_mol_gth_module
  use ps_nloc2_op_module, only: init_op_ps_nloc2_hp
  use ylm_module
  use hsort_module
  use polint_module

  implicit none

  PRIVATE
  PUBLIC :: prep_ps_nloc2_mol

CONTAINS


  SUBROUTINE prep_ps_nloc2_mol
    implicit none
    complex(8) :: ztmp0
    integer,allocatable :: icheck_tmp1(:),icheck_tmp2(:),itmp(:,:)
    integer,allocatable :: icheck_tmp3(:,:,:),icheck_tmp4(:,:,:)
    integer,allocatable :: sendmap_tmp(:,:),recvmap_tmp(:,:),ireq(:)
    integer,allocatable :: lma_nsend_tmp(:),maps_tmp(:,:),itmp1(:)
    integer,allocatable :: irad(:,:),nl_rank_map_tmp(:),itmp3(:,:)
    integer,allocatable :: itmp2(:),LLp(:,:),LLL(:,:,:)
    integer,allocatable :: JJ_tmp(:,:,:,:),MJJ_tmp(:,:)
    integer,allocatable :: jtmp3(:,:,:),mtmp3(:),istatus(:,:)
    integer :: a,i,j,k,L,m,n,mm1,mm2,mm3,m1,m2,ML0,k1,k2,k3
    integer :: i1,i2,i3,j1,j2,j3,ik,ir,m0,iorb,mm,ierr,ir0,irlma
    integer :: ic1,ic2,ic3,id1,id2,id3,ii1,ii2,ii3,iii1,iii2,iii3
    integer :: Nintp_0,nzlma_0,M_irad,NRc,MMJJ_0,lma,lma0,i1_0,i2_0,i3_0
    integer :: nreq,ibuf(3,3),irank
    real(8),parameter :: ep=1.d-8
    real(8) :: x,y,z,r,Rx,Ry,Rz,Rps2,v,v0,d1,d2,d3,r2,kr,pi2
    real(8) :: tmp0,tmp1,tmp2,tmp3,c1,c2,c3,maxerr,err0,err
    real(8),allocatable :: uV_tmp(:,:,:),work(:)
    real(8) :: rmin(3),rmax(0:3),rsize
    integer :: ML1,ML2,ML3,a1b,b1b,a2b,b2b,a3b,b3b
    integer :: MMJJ_t,ab1,ab2,ab3
    integer :: np1,np2,np3,nrlma

    if ( DISP_SWITCH_PARALLEL ) then
       write(*,'(a60," prep_ps_nloc2_mol")') repeat("-",60)
    end if

    Mlma=0
    do i=1,Natom
       ik=ki_atom(i)
       do iorb=1,norb(ik)
          Mlma=Mlma+2*lo(iorb,ik)+1
       end do
    end do

    if ( Mlma <= 0 ) return

    if ( any( ippform == 4 ) ) then
       call init_ps_nloc_mol_gth( disp_switch_parallel )
    end if

    call GetGridSize_RgridMol( Rsize_out=rsize )

    rmin(:)= 1.d100
    rmax(:)=-1.d100
    do a=1,Natom
       r2 = sum( aa_atom(1:3,a)**2 )
       rmax(0) = max( sqrt(r2), rmax(0) )
       rmin(1) = min( aa_atom(1,a), rmin(1) )
       rmin(2) = min( aa_atom(2,a), rmin(2) )
       rmin(3) = min( aa_atom(3,a), rmin(3) )
       rmax(1) = max( aa_atom(1,a), rmax(1) )
       rmax(2) = max( aa_atom(2,a), rmax(2) )
       rmax(3) = max( aa_atom(3,a), rmax(3) )
    end do
    if ( DISP_SWITCH_PARALLEL ) then
       write(*,*) "rmin(1),rmin(1)+max(Rps)=",rmin(1),rmin(1)-maxval(Rps)
       write(*,*) "rmin(2),rmin(2)+max(Rps)=",rmin(2),rmin(2)-maxval(Rps)
       write(*,*) "rmin(3),rmin(3)+max(Rps)=",rmin(3),rmin(3)-maxval(Rps)
       write(*,*) "rmax(0),rmax(0)+max(Rps)=",rmax(0),rmax(0)+maxval(Rps)
       write(*,*) "rmax(1),rmax(1)+max(Rps)=",rmax(1),rmax(1)+maxval(Rps)
       write(*,*) "rmax(2),rmax(2)+max(Rps)=",rmax(2),rmax(2)+maxval(Rps)
       write(*,*) "rmax(3),rmax(3)+max(Rps)=",rmax(3),rmax(3)+maxval(Rps)
    end if
    rmin(1)=rmin(1)-maxval(Rps)
    rmin(2)=rmin(2)-maxval(Rps)
    rmin(3)=rmin(3)-maxval(Rps)
    rmax(0)=rmax(0)+maxval(Rps)
    rmax(1)=rmax(1)+maxval(Rps)
    rmax(2)=rmax(2)+maxval(Rps)
    rmax(3)=rmax(3)+maxval(Rps)
    if ( abs(rmin(1))>(Hgrid(1)*Ngrid(1)/2) .or. &
         abs(rmax(1))>(Hgrid(1)*Ngrid(1)/2) .or. &
         abs(rmin(2))>(Hgrid(2)*Ngrid(2)/2) .or. &
         abs(rmax(2))>(Hgrid(2)*Ngrid(2)/2) .or. &
         abs(rmin(3))>(Hgrid(3)*Ngrid(3)/2) .or. &
         abs(rmax(3))>(Hgrid(3)*Ngrid(3)/2) .or. &
             rmax(0) > Rsize     ) then
       if ( myrank == 0 ) write(*,*) "Rsize is too small!!"
       call mpi_finalize(ierr)
       stop "stop@ps_nloc2_mol"
    end if

    if ( Mlma < nprocs_g ) then
       nzlma_0 = Mlma
    else
       nzlma_0 = min(Mlma*125/nprocs_g,Mlma)
    end if

    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)

    r=maxval(Rps)+Hsize+1.d-8
    call make_minimal_box(r,mm1,mm2,mm3,MMJJ_0)
    mm1 = maxval( abs(mcube_grid_ion(:,1)) ) + 1
    mm2 = maxval( abs(mcube_grid_ion(:,2)) ) + 1
    mm3 = maxval( abs(mcube_grid_ion(:,3)) ) + 1

    MMJJ_0 = M_grid_ion

    L=maxval(lo)
    n=maxval(norb)
    allocate( icheck_tmp3(Natom,n,2*L+1) ) ; icheck_tmp3=0
    allocate( JJ_tmp(3,MMJJ_0,n,Natom)   ) ; JJ_tmp=0
    allocate( MJJ_tmp(n,Natom)           ) ; MJJ_tmp=0
    allocate( uV_tmp(MMJJ_0,n,Natom)     ) ; uV_tmp=0.d0

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

    maxerr             = 0
    icheck_tmp3(:,:,:) = 0
    MMJJ               = 0
    MMJJ_t             = 0
    nzlma              = 0
    lma                = 0
    lma0               = 0

    do a=1,Natom

       Rx = aa_atom(1,a)
       Ry = aa_atom(2,a)
       Rz = aa_atom(3,a)

       ic1 = nint( aa_atom(1,a)/Hsize )
       ic2 = nint( aa_atom(2,a)/Hsize )
       ic3 = nint( aa_atom(3,a)/Hsize )

       ik = ki_atom(a)

       do iorb=1,norb(ik)

          Rps2 = Rps(iorb,ik)**2
          NRc  = NRps(iorb,ik)
          L    = lo(iorb,ik)
          j    = 0

          do i=1,M_grid_ion

             id1 = ic1 + map_grid_ion(1,i)
             id2 = ic2 + map_grid_ion(2,i)
             id3 = ic3 + map_grid_ion(3,i)

             if ( Igrid(1,1) <= id1 .and. id1 <= Igrid(2,1) .and. &
                  Igrid(1,2) <= id2 .and. id2 <= Igrid(2,2) .and. &
                  Igrid(1,3) <= id3 .and. id3 <= Igrid(2,3) ) then

                x = id1*Hsize - Rx
                y = id2*Hsize - Ry
                z = id3*Hsize - Rz

                r2 = x*x + y*y + z*z

                if ( r2 > Rps2+1.d-10 ) cycle

                r    = sqrt(r2)
                v0   = 0.d0
                err0 = 0.d0

                if ( abs(x)>1.d-14 .or. abs(y)>1.d-14 .or. &
                     abs(z)>1.d-14 .or. L==0 ) then
                   ir0=irad( int(100.d0*r),ik )
                   do ir=ir0,NRc
                      if ( r<rad1(ir,ik) ) exit
                   end do
                   if ( ir <= 2 ) then
                      v0=viod(2,iorb,ik)
                      if ( ir < 1 ) stop "ps_nloc2_mol(0)"
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
                      write(*,*) "ps_nloc2_mol(1)",ir,NRc,viod(NRc,iorb,ik)
                      write(*,*) viod(NRc+1,iorb,ik),r,rad1(ir,ik)
                      stop
                   end if
                   maxerr=max(maxerr,err0)
                end if

!                if ( v0 /= 0.d0 ) then
                   nzlma = lma + 2*L+1
                   j=j+1
                   JJ_tmp(1,j,iorb,a) = id1
                   JJ_tmp(2,j,iorb,a) = id2
                   JJ_tmp(3,j,iorb,a) = id3
                   uV_tmp(j,iorb,a)   = v0
!                end if

             end if

          end do ! i ( 1 - M_grid_ion )

          MMJJ = max( MMJJ, j )
          if ( nzlma == lma+2*L+1 ) then
             do m=1,2*L+1
                lma=lma+1
                icheck_tmp3(a,iorb,m)=lma
             end do
             MJJ_tmp(iorb,a)=j
          end if

       end do ! iorb
    end do ! a

    deallocate( irad )


    if ( iswitch_eqdiv == 2 ) then

       if ( allocated(MJJ_MAP) ) then
          deallocate(MJJ_MAP)
          deallocate(JJ_MAP)
          deallocate(iorbmap)
          deallocate(mmap)
          deallocate(lmap)
          deallocate(amap)
          deallocate(iuV)
          deallocate(uVk)
          deallocate(JJP)
          deallocate(MJJ)
          deallocate(lma_nsend)
       end if

       allocate( lma_nsend(0:nprocs_g-1) ) ; lma_nsend=0
       allocate( MJJ(nzlma) ) ; MJJ=0
       allocate( JJP(MMJJ,nzlma) ) ; JJP=0
       allocate( uVk(MMJJ,nzlma,1) ) ; uVk=0.0d0
       allocate( iuV(nzlma) ) ; iuV=0

       allocate( amap(nzlma)    ) ; amap=0.0d0
       allocate( lmap(nzlma)    ) ; lmap=0.0d0
       allocate( mmap(nzlma)    ) ; mmap=0.0d0
       allocate( iorbmap(nzlma) ) ; iorbmap=0.0d0

       allocate( JJ_MAP(3,MMJJ,nzlma) ) ; JJ_MAP=0
       allocate( MJJ_MAP(nzlma)       ) ; MJJ_MAP=0

       allocate( icheck_tmp1(0:nprocs_g-1) ) ; icheck_tmp1=0
       allocate( icheck_tmp2(0:nprocs_g-1) ) ; icheck_tmp2=0
       allocate( LLL(a1b:b1b,a2b:b2b,a3b:b3b) ) ; LLL=0

       do i=ML_0,ML_1
          LLL( LL(1,i),LL(2,i),LL(3,i) ) = i
       end do

       do a=1,Natom
          ik=ki_atom(a)
          Rx=aa_atom(1,a)
          Ry=aa_atom(2,a)
          Rz=aa_atom(3,a)
          do iorb=1,norb(ik)
             L=lo(iorb,ik)
             do m=-L,L

                icheck_tmp1(:)=0
                call mpi_allgather(icheck_tmp3(a,iorb,m+L+1),1,mpi_integer &
                     ,icheck_tmp1,1,mpi_integer,comm_grid,ierr)

                do irank=0,nprocs_g-1

                   if ( icheck_tmp1(irank) /= 0 ) then

                      icheck_tmp2(irank)=icheck_tmp2(irank)+1
                      lma=icheck_tmp2(irank)

                      if ( irank == myrank_g ) then

                         amap(lma) = a
                         lmap(lma) = L
                         mmap(lma) = m
                         iorbmap(lma) = iorb

                         iuV(lma)=inorm(iorb,ik)
                         MJJ(lma)=MJJ_tmp(iorb,a)
                         MJJ_MAP(lma)=MJJ_tmp(iorb,a)
                         do j=1,MJJ(lma)
                            i1=JJ_tmp(1,j,iorb,a)
                            i2=JJ_tmp(2,j,iorb,a)
                            i3=JJ_tmp(3,j,iorb,a)
                            JJ_MAP(1:3,j,lma)=JJ_tmp(1:3,j,iorb,a)
                            x=i1*Hsize-Rx
                            y=i2*Hsize-Ry
                            z=i3*Hsize-Rz
                            JJP(j,lma)=LLL(i1,i2,i3)
                            uVk(j,lma,1)=uV_tmp(j,iorb,a)*Ylm(x,y,z,L,m)
                         end do ! j

                      else

                         if ( icheck_tmp1(myrank_g) /= 0 ) then
                            lma_nsend(irank) = lma_nsend(irank) + 1
                         end if

                      end if

                   end if

                end do ! irank

             end do ! m
          end do ! iorb
       end do ! a

       deallocate( LLL )

       if ( allocated(sendmap) ) then
          deallocate( recvmap )
          deallocate( sendmap )
       end if

       nl_max_send=maxval( lma_nsend )
       allocate( sendmap(nl_max_send,0:nprocs_g-1) ) ; sendmap=0
       allocate( recvmap(nl_max_send,0:nprocs_g-1) ) ; recvmap=0

       lma_nsend(:)=0
       icheck_tmp2(:)=0
       do a=1,Natom
          ik=ki_atom(a)
          do iorb=1,norb(ik)
             L=lo(iorb,ik)
             do m=-L,L

                icheck_tmp1(:)=0
                call mpi_allgather(icheck_tmp3(a,iorb,m+L+1),1,mpi_integer &
                     ,icheck_tmp1,1,mpi_integer,comm_grid,ierr)

                do irank=0,nprocs_g-1
                   if ( icheck_tmp1(irank) == 0 ) cycle
                   icheck_tmp2(irank) = icheck_tmp2(irank) + 1
                end do

                do irank=0,nprocs_g-1

                   if ( icheck_tmp1(irank) == 0 .or. icheck_tmp1(myrank_g) == 0 .or. irank == myrank_g ) cycle

                   lma_nsend(irank) = lma_nsend(irank) + 1

                   sendmap( lma_nsend(irank),irank ) = icheck_tmp2(myrank_g)

                   recvmap( lma_nsend(irank),irank ) = icheck_tmp2(myrank_g)

                end do

             end do ! m
          end do ! iorb
       end do ! a

       deallocate( icheck_tmp2 )
       deallocate( icheck_tmp1 )

       call allocate_ps_nloc2(MB_d)

    else


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

    do a=1,Natom
       ik=ki_atom(a)
    do iorb=1,norb(ik)
       L=lo(iorb,ik)
    do m=-L,L
       icheck_tmp1(:)=0
       call mpi_allgather(icheck_tmp3(a,iorb,m+L+1),1,mpi_integer &
                         ,icheck_tmp1,1,mpi_integer,comm_grid,ierr)
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
          stop "ps_nloc2_mol(2)"
!          i1=0
!          j1=np1
!          do i=id1+1,np1
!             if ( i1==0 .and. itmp(i,1)/=0 ) then
!                i1=i
!             else if ( i1/=0 .and. itmp(i,1)==0 ) then
!                j1=i-1
!                exit
!             end if
!          end do
!          i1=i1-np1
!          j1=j1-np1
!          ic1=i1
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
          stop "ps_nloc2_mol(3)"
!          i2=0
!          j2=np2
!          do i=id2+1,np2
!             if ( i2==0 .and. itmp(i,2)/=0 ) then
!                i2=i
!             else if ( i2/=0 .and. itmp(i,2)==0 ) then
!                j2=i-1
!                exit
!             end if
!          end do
!          i2=i2-np2
!          j2=j2-np2
!          ic2=i2
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
          stop "ps_nloc2_mol(4)"
!          i3=0
!          j3=np3
!          do i=id3+1,np3
!             if ( i3==0 .and. itmp(i,3)/=0 ) then
!                i3=i
!             else if ( i3/=0 .and. itmp(i,3)==0 ) then
!                j3=i-1
!                exit
!             end if
!          end do
!          i3=i3-np3
!          j3=j3-np3
!          ic3=i3
       end if
       do k3=ic3,id3
       do k2=ic2,id2
       do k1=ic1,id1
!          k1=mod(j1+np1-1,np1)+1
!          k2=mod(j2+np2-1,np2)+1
!          k3=mod(j3+np3-1,np3)+1
          k = k1-1 + (k2-1)*np1 + (k3-1)*np1*np2
          if ( icheck_tmp1(k) == 0 ) icheck_tmp1(k) = -1
       end do
       end do
       end do
       do n=0,nprocs_g-1
          if ( icheck_tmp1(n) /= 0 ) icheck_tmp2(n) = icheck_tmp2(n) + 1
       end do
       if ( icheck_tmp1(myrank_g) /= 0 ) then
          if ( icheck_tmp1(myrank_g) > 0 ) then
             maps_tmp( icheck_tmp2(myrank_g),1 ) = icheck_tmp1(myrank_g)
          end if
          maps_tmp( icheck_tmp2(myrank_g),2 ) = inorm(iorb,ik)
          maps_tmp( icheck_tmp2(myrank_g),3 ) = a
          maps_tmp( icheck_tmp2(myrank_g),4 ) = L
          maps_tmp( icheck_tmp2(myrank_g),5 ) = m
          maps_tmp( icheck_tmp2(myrank_g),6 ) = iorb

          do n=0,nprocs_g-1
             if ( n == myrank_g .or. icheck_tmp1(n) == 0 ) cycle
             lma_nsend_tmp(n) = lma_nsend_tmp(n) + 1  
             sendmap_tmp( lma_nsend_tmp(n),n ) = icheck_tmp2(myrank_g)
             recvmap_tmp( lma_nsend_tmp(n),n ) = icheck_tmp2(n)
             if ( any(nl_rank_map_tmp(0:nrlma)==n) ) cycle
             nrlma = nrlma + 1
             nl_rank_map_tmp(nrlma) = n
          end do
       end if

    end do ! m
    end do ! iorb
    end do ! a

    nzlma = icheck_tmp2(myrank_g)

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
    allocate( JJ_MAP(3,MMJJ,nzlma) ) ; JJ_MAP=0
    allocate( MJJ_MAP(nzlma)       ) ; MJJ_MAP=0
    allocate( MJJ(nzlma)     ) ; MJJ=0
    allocate( iuV(nzlma)     ) ; iuV=0
    allocate( amap(nzlma)    ) ; amap=0
    allocate( lmap(nzlma)    ) ; lmap=0
    allocate( mmap(nzlma)    ) ; mmap=0
    allocate( iorbmap(nzlma) ) ; iorbmap=0
    allocate( nl_rank_map(nrlma) ) ; nl_rank_map=-1

    do lma=1,nzlma
       if ( maps_tmp(lma,1) == 0 ) cycle
       iuV(lma)     = maps_tmp(lma,2)
       amap(lma)    = maps_tmp(lma,3)
       lmap(lma)    = maps_tmp(lma,4)
       mmap(lma)    = maps_tmp(lma,5)
       iorbmap(lma) = maps_tmp(lma,6)
    end do

    do lma=1,nzlma
       if ( maps_tmp(lma,1) == 0 ) cycle
       a    = amap(lma)
       l    = lmap(lma)
       m    = mmap(lma)
       iorb = iorbmap(lma)
       MJJ_MAP(lma) = MJJ_tmp(iorb,a)
       Rx=aa_atom(1,a)
       Ry=aa_atom(2,a)
       Rz=aa_atom(3,a)
       do j=1,MJJ_MAP(lma)
          x = JJ_tmp(1,j,iorb,a)*Hsize - Rx
          y = JJ_tmp(2,j,iorb,a)*Hsize - Ry
          z = JJ_tmp(3,j,iorb,a)*Hsize - Rz
          uV(j,lma) = uV_tmp(j,iorb,a)*Ylm(x,y,z,l,m)
          JJ_MAP(1:3,j,lma) = JJ_tmp(1:3,j,iorb,a)
       end do
    end do
     
    deallocate( MJJ_tmp )
    deallocate( JJ_tmp )
    deallocate( uV_tmp )


    allocate( icheck_tmp4(a1b:b1b,a2b:b2b,a3b:b3b) )
    icheck_tmp4=0
    do lma=1,nzlma
       j=0
       icheck_tmp4=0
       do i=1,MJJ_MAP(lma)
          i1=JJ_MAP(1,i,lma)
          i2=JJ_MAP(2,i,lma)
          i3=JJ_MAP(3,i,lma)
          if ( icheck_tmp4(i1,i2,i3) == 0 ) then
             j=j+1
             icheck_tmp4(i1,i2,i3) = j
          end if
       end do
       MJJ(lma)=j
    end do
    MAXMJJ = maxval( MJJ(1:nzlma) )
    deallocate( icheck_tmp4 )


    do i=1,nrlma
       nl_rank_map(i)=nl_rank_map_tmp(i)
    end do

    deallocate( itmp )
    deallocate( nl_rank_map_tmp )
    deallocate( icheck_tmp2 )
    deallocate( icheck_tmp1 )
    deallocate( maps_tmp )
    deallocate( icheck_tmp3 )

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
       itmp(1,irlma) = LLp(1,n) - LLp(1,myrank_g)
       itmp(2,irlma) = LLp(2,n) - LLp(2,myrank_g)
       itmp(3,irlma) = LLp(3,n) - LLp(3,myrank_g)
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

    call allocate_ps_nloc2(MB_d)

    if ( allocated(JJP) ) then
       deallocate( JJP )
       deallocate( uVk )
    end if
    allocate( JJP(MAXMJJ,nzlma) ) ; JJP=0
    allocate( uVk(MAXMJJ,nzlma,1) )


    allocate( LLL(a1b:b1b,a2b:b2b,a3b:b3b) ) ; LLL=0
    do i=ML_0,ML_1
       LLL( LL(1,i),LL(2,i),LL(3,i) ) = i
    end do
    allocate( icheck_tmp4(a1b:b1b,a2b:b2b,a3b:b3b) )
    icheck_tmp4=0

    do lma=1,nzlma
       j=0
       icheck_tmp4=0
       do i=1,MJJ_MAP(lma)
          i1 = JJ_MAP(1,i,lma)
          i2 = JJ_MAP(2,i,lma)
          i3 = JJ_MAP(3,i,lma)
          j3 = icheck_tmp4(i1,i2,i3)
          if ( j3 == 0 ) then
             j=j+1
             icheck_tmp4(i1,i2,i3)=j
             uVk(j,lma,1) = uV(i,lma)
             JJP(j,lma) = LLL(i1,i2,i3)
          else
             uVk(j3,lma,1) = uVk(j3,lma,1) + uV(i,lma)
          end if
       end do
    end do ! lma

    deallocate( icheck_tmp4 )
    deallocate( LLL )


    end if ![ iswitch_eqdiv ]

    call init_op_ps_nloc2_hp( .true. )

  END SUBROUTINE prep_ps_nloc2_mol


END MODULE ps_nloc2_mol_module

