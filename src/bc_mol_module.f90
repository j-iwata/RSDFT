MODULE bc_mol_module

  use bc_variables, only: fdinfo_send,fdinfo_recv,n_neighbor
  use rgrid_mol_module, only: map_g2p_rgrid_mol

  implicit none

  PRIVATE
  PUBLIC :: init_bcset_mol

  include 'mpif.h'

CONTAINS

  SUBROUTINE init_bcset_mol(Md,Ngrid,nproc,myrnk,comm,pinfo_grid)
    implicit none
    integer,intent(IN) :: Ngrid(3),Md,nproc,myrnk,comm
    integer,intent(IN) :: pinfo_grid(8,0:nproc-1)
    real(8) :: H,Rc2,r2,z
    integer :: n,a1,a2,a3,b1,b2,b3,ix,iy,iz,i,j,jx,jy,jz
    integer :: ix0,ix1,iy0,iy1,iz0,iz1,irank,i1,i2,i3
    integer :: mx,my,mz,nreq,itags,itagr,ierr
    integer,allocatable :: map_grid_2_pinfo(:,:,:)
    integer,allocatable :: neighbor_info(:,:,:)
    integer,allocatable :: itmp(:),jtmp(:)
    integer,allocatable :: ireq(:),istatus(:,:)
    include 'mpif.h'

    if ( myrnk == 0 ) write(*,'(a60," init_bcset_mol(start)")') repeat("-",60)

    mx  = (Ngrid(1)-1)/2
    my  = (Ngrid(2)-1)/2
    mz  = (Ngrid(3)-1)/2

    a1=-mx-Md ; a2=-my-Md ; a3=-mz-Md
    b1= mx+Md ; b2= my+Md ; b3= mz+Md
    allocate( map_grid_2_pinfo(a1:b1,a2:b2,a3:b3) )
    map_grid_2_pinfo(:,:,:)=-1

    call map_g2p_rgrid_mol(a1,b1,a2,b2,a3,b3,map_grid_2_pinfo,nproc,pinfo_grid)

!
! --- FD Boundary conditions ( eqdiv ) ---
!

    allocate( neighbor_info(8,6,100) )
    allocate( itmp(0:nproc-1) )
    allocate( jtmp(0:nproc-1) )

    n_neighbor(:)=0
    neighbor_info(:,:,:)=0
    neighbor_info(1:2,:,:)=MPI_PROC_NULL
    itmp(:)=-1
    jtmp(:)=0

    iy0=    pinfo_grid(3,myrnk)+1
    iy1=iy0+pinfo_grid(4,myrnk)-1
    iz0=    pinfo_grid(5,myrnk)+1
    iz1=iz0+pinfo_grid(6,myrnk)-1

    ix0=pinfo_grid(1,myrnk)+1-Md
    ix1=pinfo_grid(1,myrnk)
    itmp(:)=-1
    jtmp(:)= 0
    do ix=ix0,ix1
       if ( ix < -mx .or. mx < ix ) cycle
       do iz=iz0,iz1
       do iy=iy0,iy1
          irank = map_grid_2_pinfo(ix,iy,iz)
          if ( 0 <= irank .and. irank <= nproc-1 ) jtmp(irank)=jtmp(irank)+1
       end do
       end do
    end do
    n_neighbor(1)=0
    do i=0,nproc-1
       if ( jtmp(i)>0 ) then
          n_neighbor(1)=n_neighbor(1)+1
          itmp( n_neighbor(1) )=i
       end if
    end do
    do i=1,n_neighbor(1)
       irank=itmp(i)
       ix=max( ix0, pinfo_grid(1,irank)+1 )
       jx=min( ix1, pinfo_grid(1,irank)+pinfo_grid(2,irank) )
       iy=max( iy0, pinfo_grid(3,irank)+1 )
       jy=min( iy1, pinfo_grid(3,irank)+pinfo_grid(4,irank) )
       iz=max( iz0, pinfo_grid(5,irank)+1 )
       jz=min( iz1, pinfo_grid(5,irank)+pinfo_grid(6,irank) )
       neighbor_info( 1, 1, i )=myrnk
       neighbor_info( 2, 1, i )=irank
       neighbor_info( 3, 1, i )=ix
       neighbor_info( 4, 1, i )=jx
       neighbor_info( 5, 1, i )=iy
       neighbor_info( 6, 1, i )=jy
       neighbor_info( 7, 1, i )=iz
       neighbor_info( 8, 1, i )=jz
    end do

    ix0=pinfo_grid(1,myrnk)+pinfo_grid(2,myrnk)+1
    ix1=pinfo_grid(1,myrnk)+pinfo_grid(2,myrnk)+Md
    itmp(:)=-1
    jtmp(:)= 0
    do ix=ix0,ix1
       if ( ix < -mx .or. mx < ix ) cycle
       do iz=iz0,iz1
       do iy=iy0,iy1
          irank = map_grid_2_pinfo(ix,iy,iz)
          if ( 0 <= irank .and. irank <= nproc-1 ) jtmp(irank)=jtmp(irank)+1
       end do
       end do
    end do
    n_neighbor(2)=0
    do i=0,nproc-1
       if ( jtmp(i)>0 ) then
          n_neighbor(2)=n_neighbor(2)+1
          itmp( n_neighbor(2) )=i
       end if
    end do
    do i=1,n_neighbor(2)
       irank=itmp(i)
       ix=max( ix0, pinfo_grid(1,irank)+1 )
       jx=min( ix1, pinfo_grid(1,irank)+pinfo_grid(2,irank) )
       iy=max( iy0, pinfo_grid(3,irank)+1 )
       jy=min( iy1, pinfo_grid(3,irank)+pinfo_grid(4,irank) )
       iz=max( iz0, pinfo_grid(5,irank)+1 )
       jz=min( iz1, pinfo_grid(5,irank)+pinfo_grid(6,irank) )
       neighbor_info( 1, 2, i )=myrnk
       neighbor_info( 2, 2, i )=irank
       neighbor_info( 3, 2, i )=ix
       neighbor_info( 4, 2, i )=jx
       neighbor_info( 5, 2, i )=iy
       neighbor_info( 6, 2, i )=jy
       neighbor_info( 7, 2, i )=iz
       neighbor_info( 8, 2, i )=jz
    end do

    ix0=    pinfo_grid(1,myrnk)+1
    ix1=ix0+pinfo_grid(2,myrnk)-1
    iz0=    pinfo_grid(5,myrnk)+1
    iz1=iz0+pinfo_grid(6,myrnk)-1

    iy0=pinfo_grid(3,myrnk)+1-Md
    iy1=pinfo_grid(3,myrnk)
    itmp(:)=-1
    jtmp(:)= 0
    do iy=iy1,iy0,-1
       if ( iy < -my .or. my < iy ) cycle
       do iz=iz0,iz1
       do ix=ix0,ix1
          irank = map_grid_2_pinfo(ix,iy,iz)
          if ( 0 <= irank .and. irank <= nproc-1 ) jtmp(irank)=jtmp(irank)+1
       end do
       end do
    end do
    n_neighbor(3)=0
    do i=0,nproc-1
       if ( jtmp(i)>0 ) then
          n_neighbor(3)=n_neighbor(3)+1
          itmp( n_neighbor(3) )=i
       end if
    end do
    do i=1,n_neighbor(3)
       irank=itmp(i)
       ix=max( ix0, pinfo_grid(1,irank)+1 )
       jx=min( ix1, pinfo_grid(1,irank)+pinfo_grid(2,irank) )
       iy=max( iy0, pinfo_grid(3,irank)+1 )
       jy=min( iy1, pinfo_grid(3,irank)+pinfo_grid(4,irank) )
       iz=max( iz0, pinfo_grid(5,irank)+1 )
       jz=min( iz1, pinfo_grid(5,irank)+pinfo_grid(6,irank) )
       neighbor_info( 1, 3, i )=myrnk
       neighbor_info( 2, 3, i )=irank
       neighbor_info( 3, 3, i )=ix
       neighbor_info( 4, 3, i )=jx
       neighbor_info( 5, 3, i )=iy
       neighbor_info( 6, 3, i )=jy
       neighbor_info( 7, 3, i )=iz
       neighbor_info( 8, 3, i )=jz
    end do

    iy0=pinfo_grid(3,myrnk)+pinfo_grid(4,myrnk)+1
    iy1=pinfo_grid(3,myrnk)+pinfo_grid(4,myrnk)+Md
    itmp(:)=-1
    jtmp(:)= 0
    do iy=iy0,iy1
       if ( iy < -my .or. my < iy ) cycle
       do iz=iz0,iz1
       do ix=ix0,ix1
          irank = map_grid_2_pinfo(ix,iy,iz)
          if ( 0 <= irank .and. irank <= nproc-1 ) jtmp(irank)=jtmp(irank)+1
       end do
       end do
    end do
    n_neighbor(4)=0
    do i=0,nproc-1
       if ( jtmp(i)>0 ) then
          n_neighbor(4)=n_neighbor(4)+1
          itmp( n_neighbor(4) )=i
       end if
    end do
    do i=1,n_neighbor(4)
       irank=itmp(i)
       ix=max( ix0, pinfo_grid(1,irank)+1 )
       jx=min( ix1, pinfo_grid(1,irank)+pinfo_grid(2,irank) )
       iy=max( iy0, pinfo_grid(3,irank)+1 )
       jy=min( iy1, pinfo_grid(3,irank)+pinfo_grid(4,irank) )
       iz=max( iz0, pinfo_grid(5,irank)+1 )
       jz=min( iz1, pinfo_grid(5,irank)+pinfo_grid(6,irank) )
       neighbor_info( 1, 4, i )=myrnk
       neighbor_info( 2, 4, i )=irank
       neighbor_info( 3, 4, i )=ix
       neighbor_info( 4, 4, i )=jx
       neighbor_info( 5, 4, i )=iy
       neighbor_info( 6, 4, i )=jy
       neighbor_info( 7, 4, i )=iz
       neighbor_info( 8, 4, i )=jz
    end do

    ix0=    pinfo_grid(1,myrnk)+1
    ix1=ix0+pinfo_grid(2,myrnk)-1
    iy0=    pinfo_grid(3,myrnk)+1
    iy1=iy0+pinfo_grid(4,myrnk)-1

    iz0=pinfo_grid(5,myrnk)+1-Md
    iz1=pinfo_grid(5,myrnk)
    itmp(:)=-1
    jtmp(:)= 0
    do iz=iz1,iz0,-1
       if ( iz < -mz .or. mz < iz ) cycle
       do iy=iy0,iy1
       do ix=ix0,ix1
          irank = map_grid_2_pinfo(ix,iy,iz)
          if ( 0 <= irank .and. irank <= nproc-1 ) jtmp(irank)=jtmp(irank)+1
       end do
       end do
    end do
    n_neighbor(5)=0
    do i=0,nproc-1
       if ( jtmp(i)>0 ) then
          n_neighbor(5)=n_neighbor(5)+1
          itmp( n_neighbor(5) )=i
       end if
    end do
    do i=1,n_neighbor(5)
       irank=itmp(i)
       ix=max( ix0, pinfo_grid(1,irank)+1 )
       jx=min( ix1, pinfo_grid(1,irank)+pinfo_grid(2,irank) )
       iy=max( iy0, pinfo_grid(3,irank)+1 )
       jy=min( iy1, pinfo_grid(3,irank)+pinfo_grid(4,irank) )
       iz=max( iz0, pinfo_grid(5,irank)+1 )
       jz=min( iz1, pinfo_grid(5,irank)+pinfo_grid(6,irank) )
       neighbor_info( 1, 5, i )=myrnk
       neighbor_info( 2, 5, i )=irank
       neighbor_info( 3, 5, i )=ix
       neighbor_info( 4, 5, i )=jx
       neighbor_info( 5, 5, i )=iy
       neighbor_info( 6, 5, i )=jy
       neighbor_info( 7, 5, i )=iz
       neighbor_info( 8, 5, i )=jz
    end do

    iz0=pinfo_grid(5,myrnk)+pinfo_grid(6,myrnk)+1
    iz1=pinfo_grid(5,myrnk)+pinfo_grid(6,myrnk)+Md
    itmp(:)=-1
    jtmp(:)= 0
    do iz=iz0,iz1
       if ( iz < -mz .or. mz < iz ) cycle
       do iy=iy0,iy1
       do ix=ix0,ix1
          irank = map_grid_2_pinfo(ix,iy,iz)
          if ( 0 <= irank .and. irank <= nproc-1 ) jtmp(irank)=jtmp(irank)+1
       end do
       end do
    end do
    n_neighbor(6)=0
    do i=0,nproc-1
       if ( jtmp(i)>0 ) then
          n_neighbor(6)=n_neighbor(6)+1
          itmp( n_neighbor(6) )=i
       end if
    end do
    do i=1,n_neighbor(6)
       irank=itmp(i)
       ix=max( ix0, pinfo_grid(1,irank)+1 )
       jx=min( ix1, pinfo_grid(1,irank)+pinfo_grid(2,irank) )
       iy=max( iy0, pinfo_grid(3,irank)+1 )
       jy=min( iy1, pinfo_grid(3,irank)+pinfo_grid(4,irank) )
       iz=max( iz0, pinfo_grid(5,irank)+1 )
       jz=min( iz1, pinfo_grid(5,irank)+pinfo_grid(6,irank) )
       neighbor_info( 1, 6, i )=myrnk
       neighbor_info( 2, 6, i )=irank
       neighbor_info( 3, 6, i )=ix
       neighbor_info( 4, 6, i )=jx
       neighbor_info( 5, 6, i )=iy
       neighbor_info( 6, 6, i )=jy
       neighbor_info( 7, 6, i )=iz
       neighbor_info( 8, 6, i )=jz
    end do

    deallocate( jtmp )
    deallocate( itmp )
    deallocate( map_grid_2_pinfo )

    if ( allocated(fdinfo_recv) ) deallocate(fdinfo_recv)
    if ( allocated(fdinfo_send) ) deallocate(fdinfo_send)

    n=maxval( n_neighbor(1:6) )
    allocate( fdinfo_send(10,n,6) ) ; fdinfo_send=0
    allocate( fdinfo_recv(10,n,6) ) ; fdinfo_recv=0

    fdinfo_send(7,:,:)=MPI_PROC_NULL
    fdinfo_recv(8,:,:)=MPI_PROC_NULL

    do i=1,6
       do n=1,n_neighbor(i)
          ix=neighbor_info( 3, i, n )
          jx=neighbor_info( 4, i, n )
          iy=neighbor_info( 5, i, n )
          jy=neighbor_info( 6, i, n )
          iz=neighbor_info( 7, i, n )
          jz=neighbor_info( 8, i, n )
          fdinfo_recv(1,n,i)=ix
          fdinfo_recv(2,n,i)=jx
          fdinfo_recv(3,n,i)=iy
          fdinfo_recv(4,n,i)=jy
          fdinfo_recv(5,n,i)=iz
          fdinfo_recv(6,n,i)=jz
          fdinfo_recv(7,n,i)=neighbor_info( 1, i, n ) !=myrnk
          fdinfo_recv(8,n,i)=neighbor_info( 2, i, n ) !=irank
          fdinfo_recv(9,n,i)=(jx-ix+1)*(jy-iy+1)*(jz-iz+1)
          select case(i)
          case(1,2)
             fdinfo_recv(10,n,i)=(jy-iy+1)*(jz-iz+1)
          case(3,4)
             fdinfo_recv(10,n,i)=(jx-ix+1)*(jz-iz+1)
          case(5,6)
             fdinfo_recv(10,n,i)=(jx-ix+1)*(jy-iy+1)
          end select
       end do
    end do

    deallocate( neighbor_info )

    n=maxval(n_neighbor)*6*2
    allocate( ireq(n) ) ; ireq(:)=0
    allocate( istatus(MPI_STATUS_SIZE,n) ) ; istatus=0

    nreq=0
    do j=1,6
       do i=1,n_neighbor(j)
          irank=fdinfo_recv(8,i,j)
          itags=1
          itagr=1
          nreq=nreq+1
          call mpi_irecv(fdinfo_send(1,i,j),10,mpi_integer,irank &
                        ,itagr,comm,ireq(nreq),ierr)
          nreq=nreq+1
          call mpi_isend(fdinfo_recv(1,i,j),10,mpi_integer,irank &
                        ,itags,comm,ireq(nreq),ierr)
       end do
    end do

    call mpi_waitall(nreq,ireq,istatus,ierr)

    deallocate( istatus )
    deallocate( ireq )

    if ( myrnk == 0 ) write(*,'(a60," init_bcset_mol(end)")') repeat("-",60)

    return

  END SUBROUTINE init_bcset_mol

END MODULE bc_mol_module
