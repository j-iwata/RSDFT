MODULE bc_module

  use rgrid_module, only: Ngrid,Igrid
  use parallel_module
  use bc_mol_module, only: init_bcset_mol
  use bc_variables, only: fdinfo_send,fdinfo_recv,n_neighbor,www,Md &
                         ,sbuf,rbuf,TYPE_MAIN,zero
  use watch_module
  use bcset_1_module, only: bcset_1
  use bcset_3_module, only: bcset_3

  implicit none

  PRIVATE
  PUBLIC :: init_bcset
  PUBLIC :: allocate_bcset
  PUBLIC :: bcset
  PUBLIC :: bcset_1
  PUBLIC :: bcset_3
  PUBLIC :: n_neighbor, fdinfo_send, fdinfo_recv
  PUBLIC :: www

CONTAINS


  SUBROUTINE bcset(ib1,ib2,ndepth,idir)
    implicit none
    integer,intent(IN) :: ib1,ib2,ndepth,idir
    integer :: a1,a2,a3,b1,b2,b3,nb,ns,ms,mt
    integer :: m,n,ndata,i1,i2,i3,ib,ierr
    integer :: c1,c2,c3,d1,d2,d3,irank,nreq,itags,itagr,ireq(36)
    integer :: i,j,l
    integer :: istatus(mpi_status_size,123)

    a1=Igrid(1,1)
    b1=Igrid(2,1)
    a2=Igrid(1,2)
    b2=Igrid(2,2)
    a3=Igrid(1,3)
    b3=Igrid(2,3)

    nb=ib2-ib1+1

    nreq=0

!(1) [sender][2]-->[1][receiver]

    if ( idir==0 .or. idir==1 .or. idir==2 ) then
       do n=1,n_neighbor(1)
          if ( ndepth==1 ) then
             ndata = fdinfo_recv(10,n,1)*nb
          else
             ndata = fdinfo_recv(9,n,1)*nb
          end if
          if ( ndata<1 ) cycle
          irank = fdinfo_recv(8,n,1)
          itagr = 10
          nreq  = nreq + 1
          call mpi_irecv(rbuf(1,n,1),ndata,TYPE_MAIN,irank,itagr &
                        ,comm_grid,ireq(nreq),ierr)
       end do
       do n=1,n_neighbor(2)
          if ( ndepth==1 ) then
             c1=b1
             d1=b1
             j=fdinfo_send(10,n,2)*nb
          else
             c1=fdinfo_send(1,n,2)
             d1=fdinfo_send(2,n,2)
             j=fdinfo_send(9,n,2)*nb
          end if
          c2=fdinfo_send(3,n,2) ; d2=fdinfo_send(4,n,2)
          c3=fdinfo_send(5,n,2) ; d3=fdinfo_send(6,n,2)
          if ( j<1 ) cycle
          ndata=0
          do ib=ib1,ib2
             do i3=c3,d3
             do i2=c2,d2
             do i1=c1,d1
                ndata=ndata+1
                sbuf(ndata,n,2)=www(i1,i2,i3,ib)
             end do
             end do
             end do
          end do
          if ( ndata/=j ) stop
          irank = fdinfo_send(7,n,2)
          itags = 10
          nreq  = nreq + 1
          call mpi_isend(sbuf(1,n,2),ndata,TYPE_MAIN,irank,itags &
                        ,comm_grid,ireq(nreq),ierr)
       end do
    end if

!    call mpi_waitall(nreq,ireq,istatus,ierr) ; nreq=0

!(3) [sender][4]-->[3][receiver]

    if ( idir==0 .or. idir==3 .or. idir==4 ) then
       do n=1,n_neighbor(3)
          if ( ndepth==1 ) then
             ndata = fdinfo_recv(10,n,3)*nb
          else
             ndata = fdinfo_recv(9,n,3)*nb
          end if
          if ( ndata<1 ) cycle
          irank = fdinfo_recv(8,n,3)
          itagr = 30
          nreq  = nreq + 1
          call mpi_irecv(rbuf(1,n,3),ndata,TYPE_MAIN,irank,itagr &
                        ,comm_grid,ireq(nreq),ierr)
       end do
       do n=1,n_neighbor(4)
          if ( ndepth==1 ) then
             c2=b2
             d2=b2
             j=fdinfo_send(10,n,4)*nb
          else
             c2=fdinfo_send(3,n,4)
             d2=fdinfo_send(4,n,4)
             j=fdinfo_send(9,n,4)*nb
          end if
          c1=fdinfo_send(1,n,4) ; d1=fdinfo_send(2,n,4)
          c3=fdinfo_send(5,n,4) ; d3=fdinfo_send(6,n,4)
          if ( j<1 ) cycle
          ndata=0
          do ib=ib1,ib2
             do i3=c3,d3
             do i2=c2,d2
             do i1=c1,d1
                ndata=ndata+1
                sbuf(ndata,n,4)=www(i1,i2,i3,ib)
             end do
             end do
             end do
          end do
          if ( ndata/=j ) stop
          irank = fdinfo_send(7,n,4)
          itags = 30
          nreq  = nreq + 1
          call mpi_isend(sbuf(1,n,4),ndata,TYPE_MAIN,irank,itags &
                        ,comm_grid,ireq(nreq),ierr)
       end do
    end if

!    call mpi_waitall(nreq,ireq,istatus,ierr) ; nreq=0

!(5) [sender][6]-->[5][receiver]

    if ( idir==0 .or. idir==5 .or. idir==6 ) then
       do n=1,n_neighbor(5)
          if ( ndepth==1 ) then
             ndata = fdinfo_recv(10,n,5)*nb
          else
             ndata = fdinfo_recv(9,n,5)*nb
          end if
          if ( ndata<1 ) cycle
          irank = fdinfo_recv(8,n,5)
          itagr = 50
          nreq  = nreq + 1
          call mpi_irecv(rbuf(1,n,5),ndata,TYPE_MAIN,irank,itagr &
                        ,comm_grid,ireq(nreq),ierr)
       end do
       do n=1,n_neighbor(6)
          if ( ndepth==1 ) then
             c3=b3
             d3=b3
             j=fdinfo_send(10,n,6)*nb
          else
             c3=fdinfo_send(5,n,6)
             d3=fdinfo_send(6,n,6)
             j=fdinfo_send(9,n,6)*nb
          end if
          c1=fdinfo_send(1,n,6) ; d1=fdinfo_send(2,n,6)
          c2=fdinfo_send(3,n,6) ; d2=fdinfo_send(4,n,6)
          if ( j<1 ) cycle
          ndata=0
          do ib=ib1,ib2
             do i3=c3,d3
             do i2=c2,d2
             do i1=c1,d1
                ndata=ndata+1
                sbuf(ndata,n,6)=www(i1,i2,i3,ib)
             end do
             end do
             end do
          end do
          if ( ndata/=j ) stop
          irank = fdinfo_send(7,n,6)
          itags = 50
          nreq  = nreq + 1
          call mpi_isend(sbuf(1,n,6),ndata,TYPE_MAIN,irank,itags &
                        ,comm_grid,ireq(nreq),ierr)
       end do
    end if

!    call mpi_waitall(nreq,ireq,istatus,ierr) ; nreq=0

!(2) [receiver][2]<--[1][sender]

    if ( idir==0 .or. idir==1 .or. idir==2 ) then
       do n=1,n_neighbor(2)
          if ( ndepth==1 ) then
             ndata = fdinfo_recv(10,n,2)*nb
          else
             ndata = fdinfo_recv(9,n,2)*nb
          end if
          if ( ndata<1 ) cycle
          irank = fdinfo_recv(8,n,2)
          itagr = 20
          nreq  = nreq + 1
          call mpi_irecv(rbuf(1,n,2),ndata,TYPE_MAIN,irank,itagr &
                        ,comm_grid,ireq(nreq),ierr)
       end do
       do n=1,n_neighbor(1)
          if ( ndepth==1 ) then
             c1=a1
             d1=a1
             j=fdinfo_send(10,n,1)*nb
          else
             c1=fdinfo_send(1,n,1)
             d1=fdinfo_send(2,n,1)
             j=fdinfo_send(9,n,1)*nb
          end if
          c2=fdinfo_send(3,n,1) ; d2=fdinfo_send(4,n,1)
          c3=fdinfo_send(5,n,1) ; d3=fdinfo_send(6,n,1)
          if ( j<1 ) cycle
          ndata=0
          do ib=ib1,ib2
             do i3=c3,d3
             do i2=c2,d2
             do i1=c1,d1
                ndata=ndata+1
                sbuf(ndata,n,1)=www(i1,i2,i3,ib)
             end do
             end do
             end do
          end do
          if ( ndata/=j ) stop
          irank = fdinfo_send(7,n,1)
          itags = 20
          nreq  = nreq + 1
          call mpi_isend(sbuf(1,n,1),ndata,TYPE_MAIN,irank,itags &
                        ,comm_grid,ireq(nreq),ierr)
       end do
    end if

!    call mpi_waitall(nreq,ireq,istatus,ierr) ; nreq=0

!(4) [receiver][4]<--[3][sender]

    if ( idir==0 .or. idir==3 .or. idir==4 ) then
       do n=1,n_neighbor(4)
          if ( ndepth==1 ) then
             ndata = fdinfo_recv(10,n,4)*nb
          else
             ndata = fdinfo_recv(9,n,4)*nb
          end if
          if ( ndata<1 ) cycle
          irank = fdinfo_recv(8,n,4)
          itagr = 40
          nreq  = nreq + 1
          call mpi_irecv(rbuf(1,n,4),ndata,TYPE_MAIN,irank,itagr &
                        ,comm_grid,ireq(nreq),ierr)
       end do
       do n=1,n_neighbor(3)
          if ( ndepth==1 ) then
             c2=a2
             d2=a2
             j=fdinfo_send(10,n,3)*nb
          else
             c2=fdinfo_send(3,n,3)
             d2=fdinfo_send(4,n,3)
             j=fdinfo_send(9,n,3)*nb
          end if
          c1=fdinfo_send(1,n,3) ; d1=fdinfo_send(2,n,3)
          c3=fdinfo_send(5,n,3) ; d3=fdinfo_send(6,n,3)
          if ( j<1 ) cycle
          ndata=0
          do ib=ib1,ib2
             do i3=c3,d3
             do i2=c2,d2
             do i1=c1,d1
                ndata=ndata+1
                sbuf(ndata,n,3)=www(i1,i2,i3,ib)
             end do
             end do
             end do
          end do
          if ( ndata/=j ) stop
          irank = fdinfo_send(7,n,3)
          itags = 40
          nreq  = nreq + 1
          call mpi_isend(sbuf(1,n,3),ndata,TYPE_MAIN,irank,itags &
                        ,comm_grid,ireq(nreq),ierr)
       end do
    end if

!    call mpi_waitall(nreq,ireq,istatus,ierr) ; nreq=0

!(6) [receiver][6]<--[5][sender]

    if ( idir==0 .or. idir==5 .or. idir==6 ) then
       do n=1,n_neighbor(6)
          if ( ndepth==1 ) then
             ndata = fdinfo_recv(10,n,6)*nb
          else
             ndata = fdinfo_recv(9,n,6)*nb
          end if
          if ( ndata<1 ) cycle
          irank = fdinfo_recv(8,n,6)
          itagr = 60
          nreq  = nreq + 1
          call mpi_irecv(rbuf(1,n,6),ndata,TYPE_MAIN,irank,itagr &
                        ,comm_grid,ireq(nreq),ierr)
       end do
       do n=1,n_neighbor(5)
          if ( ndepth==1 ) then
             c3=a3
             d3=a3
             j=fdinfo_send(10,n,5)*nb
          else
             c3=fdinfo_send(5,n,5)
             d3=fdinfo_send(6,n,5)
             j=fdinfo_send(9,n,5)*nb
          end if
          c1=fdinfo_send(1,n,5) ; d1=fdinfo_send(2,n,5)
          c2=fdinfo_send(3,n,5) ; d2=fdinfo_send(4,n,5)
          if ( j<1 ) cycle
          ndata=0
          do ib=ib1,ib2
             do i3=c3,d3
             do i2=c2,d2
             do i1=c1,d1
                ndata=ndata+1
                sbuf(ndata,n,5)=www(i1,i2,i3,ib)
             end do
             end do
             end do
          end do
          if ( ndata/=j ) stop
          irank = fdinfo_send(7,n,5)
          itags = 60
          nreq  = nreq + 1
          call mpi_isend(sbuf(1,n,5),ndata,TYPE_MAIN,irank,itags &
                        ,comm_grid,ireq(nreq),ierr)
       end do
    end if

    call mpi_waitall(nreq,ireq,istatus,ierr) ; nreq=0

    do m=1,6
       do n=1,n_neighbor(m)

          c1=fdinfo_recv(1,n,m)
          d1=fdinfo_recv(2,n,m)
          c2=fdinfo_recv(3,n,m)
          d2=fdinfo_recv(4,n,m)
          c3=fdinfo_recv(5,n,m)
          d3=fdinfo_recv(6,n,m)

          if ( Md>ndepth .and. ndepth==1 ) then
             if ( fdinfo_recv(10,n,m)<1 ) cycle
             select case(m)
             case(1)
                c1=a1-1
                d1=c1
             case(2)
                c1=b1+1
                d1=c1
             case(3)
                c2=a2-1
                d2=c2
             case(4)
                c2=b2+1
                d2=c2
             case(5)
                c3=a3-1
                d3=c3
             case(6)
                c3=b3+1
                d3=c3
             end select
          else
             if ( fdinfo_recv(9,n,m)<1 ) cycle
          end if
          i=0
          do ib=ib1,ib2
             do i3=c3,d3
             do i2=c2,d2
             do i1=c1,d1
                i=i+1
                www(i1,i2,i3,ib)=rbuf(i,n,m)
             end do
             end do
             end do
          end do

       end do
    end do

    return

  END SUBROUTINE bcset


  SUBROUTINE init_bcset( Md_in, SYStype )
    implicit none
    integer,intent(IN) :: Md_in, SYStype
    Md = Md_in
    select case( SYStype )
    case default
       call init_bcset_sol
    case( 1 )
       call init_bcset_mol(Md,Ngrid(1),np_grid,myrank_g,comm_grid,pinfo_grid)
    end select
    call allocate_bcset
  END SUBROUTINE init_bcset


  SUBROUTINE init_bcset_sol
    implicit none
    integer :: a1,a2,a3,b1,b2,b3,a1b,b1b,a2b,b2b,a3b,b3b
    integer,allocatable :: map_grid_2_pinfo(:,:,:,:)
    integer,allocatable :: ireq(:)
    integer :: i,i1,i2,i3,m,n,j,j1,j2,j3,nc,m1,m2,m3
    integer :: jrank,ip,fp,ns,irank,nreq,itags,itagr,ierr
    integer :: ML1,ML2,ML3
    integer :: istatus(MPI_STATUS_SIZE,123)

    call write_border( 80, " init_bcset_sol(start)" )

    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    a1b=Igrid(1,1) ; b1b=Igrid(2,1)
    a2b=Igrid(1,2) ; b2b=Igrid(2,2)
    a3b=Igrid(1,3) ; b3b=Igrid(2,3)

    a1=-Md ; b1=Ngrid(1)-1+Md
    a2=-Md ; b2=Ngrid(2)-1+Md
    a3=-Md ; b3=Ngrid(3)-1+Md
    
    allocate( map_grid_2_pinfo(a1:b1,a2:b2,a3:b3,2) )
    map_grid_2_pinfo(:,:,:,:)=0

    do i=0,np_grid-1
       i1=pinfo_grid(1,i)
       m1=pinfo_grid(2,i)
       i2=pinfo_grid(3,i)
       m2=pinfo_grid(4,i)
       i3=pinfo_grid(5,i)
       m3=pinfo_grid(6,i)
       do j3=i3,i3+m3-1
       do j2=i2,i2+m2-1
       do j1=i1,i1+m1-1
          map_grid_2_pinfo(j1,j2,j3,1)=i
          map_grid_2_pinfo(j1,j2,j3,2)=pinfo_grid(7,i)
       end do
       end do
       end do
    end do

    do i3=a3,b3
       j3=mod(i3+ML3,ML3)
    do i2=a2,b2
       j2=mod(i2+ML2,ML2)
    do i1=a1,b1
       j1=mod(i1+ML1,ML1)
       map_grid_2_pinfo(i1,i2,i3,1)=map_grid_2_pinfo(j1,j2,j3,1)
       map_grid_2_pinfo(i1,i2,i3,2)=map_grid_2_pinfo(j1,j2,j3,2)
    end do
    end do
    end do

    do i=2,6,2
       n=minval( pinfo_grid(i,0:np_grid-1) )
       do nc=1,Md
          if ( nc*n >= Md ) exit
       end do
       if ( nc<1 .or. Md<nc ) stop
       n_neighbor(i-1)=nc
       n_neighbor(i)  =nc
    end do

    n=maxval(n_neighbor(:))
    allocate( fdinfo_send(10,n,6) ) ; fdinfo_send=0
    allocate( fdinfo_recv(10,n,6) ) ; fdinfo_recv=0

    fdinfo_send(7,:,:)=MPI_PROC_NULL
    fdinfo_recv(8,:,:)=MPI_PROC_NULL

    jrank=map_grid_2_pinfo(a1b-1,a2b,a3b,1)
    ip=a1b-1
    fp=ip
    nc=1
    ns=(b2b-a2b+1)*(b3b-a3b+1)
    fdinfo_recv(1,nc,1)=ip
    fdinfo_recv(2,nc,1)=fp
    fdinfo_recv(3,nc,1)=a2b
    fdinfo_recv(4,nc,1)=b2b
    fdinfo_recv(5,nc,1)=a3b
    fdinfo_recv(6,nc,1)=b3b
    fdinfo_recv(7,nc,1)=myrank_g
    fdinfo_recv(8,nc,1)=jrank
    fdinfo_recv(9,nc,1)=ns
    fdinfo_recv(10,nc,1)=ns
    do i=a1b-2,a1b-Md,-1
       irank=map_grid_2_pinfo(i,a2b,a3b,1)
       if ( irank/=jrank ) then
          jrank=irank
          nc=nc+1
          fp=i
       end if
       fdinfo_recv(1,nc,1)=i
       fdinfo_recv(2,nc,1)=fp
       fdinfo_recv(3,nc,1)=a2b
       fdinfo_recv(4,nc,1)=b2b
       fdinfo_recv(5,nc,1)=a3b
       fdinfo_recv(6,nc,1)=b3b
       fdinfo_recv(7,nc,1)=myrank_g
       fdinfo_recv(8,nc,1)=jrank
       fdinfo_recv(9,nc,1)=(fp-i+1)*ns
    end do

    jrank=map_grid_2_pinfo(b1b+1,a2b,a3b,1)
    ip=b1b+1
    fp=ip
    nc=1
    ns=(b2b-a2b+1)*(b3b-a3b+1)
    fdinfo_recv(1,nc,2)=ip
    fdinfo_recv(2,nc,2)=fp
    fdinfo_recv(3,nc,2)=a2b
    fdinfo_recv(4,nc,2)=b2b
    fdinfo_recv(5,nc,2)=a3b
    fdinfo_recv(6,nc,2)=b3b
    fdinfo_recv(7,nc,2)=myrank_g
    fdinfo_recv(8,nc,2)=jrank
    fdinfo_recv(9,nc,2)=ns
    fdinfo_recv(10,nc,2)=ns
    do i=b1b+1,b1b+Md
       irank=map_grid_2_pinfo(i,a2b,a3b,1)
       if ( irank/=jrank ) then
          jrank=irank
          nc=nc+1
          ip=i
       end if
       fdinfo_recv(1,nc,2)=ip
       fdinfo_recv(2,nc,2)=i
       fdinfo_recv(3,nc,2)=a2b
       fdinfo_recv(4,nc,2)=b2b
       fdinfo_recv(5,nc,2)=a3b
       fdinfo_recv(6,nc,2)=b3b
       fdinfo_recv(7,nc,2)=myrank_g
       fdinfo_recv(8,nc,2)=jrank
       fdinfo_recv(9,nc,2)=(i-ip+1)*ns
    end do

    jrank=map_grid_2_pinfo(a1b,a2b-1,a3b,1)
    ip=a2b-1
    fp=ip
    nc=1
    ns=(b1b-a1b+1)*(b3b-a3b+1)
    fdinfo_recv(1,nc,3)=a1b
    fdinfo_recv(2,nc,3)=b1b
    fdinfo_recv(3,nc,3)=ip
    fdinfo_recv(4,nc,3)=fp
    fdinfo_recv(5,nc,3)=a3b
    fdinfo_recv(6,nc,3)=b3b
    fdinfo_recv(7,nc,3)=myrank_g
    fdinfo_recv(8,nc,3)=jrank
    fdinfo_recv(9,nc,3)=ns
    fdinfo_recv(10,nc,3)=ns
    do i=a2b-2,a2b-Md,-1
       irank=map_grid_2_pinfo(a1b,i,a3b,1)
       if ( irank/=jrank ) then
          jrank=irank
          nc=nc+1
          fp=i
       end if
       fdinfo_recv(1,nc,3)=a1b
       fdinfo_recv(2,nc,3)=b1b
       fdinfo_recv(3,nc,3)=i
       fdinfo_recv(4,nc,3)=fp
       fdinfo_recv(5,nc,3)=a3b
       fdinfo_recv(6,nc,3)=b3b
       fdinfo_recv(7,nc,3)=myrank_g
       fdinfo_recv(8,nc,3)=jrank
       fdinfo_recv(9,nc,3)=(fp-i+1)*ns
    end do

    jrank=map_grid_2_pinfo(a1b,b2b+1,a3b,1)
    ip=b2b+1
    fp=ip
    nc=1
    ns=(b1b-a1b+1)*(b3b-a3b+1)
    fdinfo_recv(1,nc,4)=a1b
    fdinfo_recv(2,nc,4)=b1b
    fdinfo_recv(3,nc,4)=ip
    fdinfo_recv(4,nc,4)=fp
    fdinfo_recv(5,nc,4)=a3b
    fdinfo_recv(6,nc,4)=b3b
    fdinfo_recv(7,nc,4)=myrank_g
    fdinfo_recv(8,nc,4)=jrank
    fdinfo_recv(9,nc,4)=ns
    fdinfo_recv(10,nc,4)=ns
    do i=b2b+1,b2b+Md
       irank=map_grid_2_pinfo(a1b,i,a3b,1)
       if ( irank/=jrank ) then
          jrank=irank
          nc=nc+1
          ip=i
       end if
       fdinfo_recv(1,nc,4)=a1b
       fdinfo_recv(2,nc,4)=b1b
       fdinfo_recv(3,nc,4)=ip
       fdinfo_recv(4,nc,4)=i
       fdinfo_recv(5,nc,4)=a3b
       fdinfo_recv(6,nc,4)=b3b
       fdinfo_recv(7,nc,4)=myrank_g
       fdinfo_recv(8,nc,4)=jrank
       fdinfo_recv(9,nc,4)=(i-ip+1)*ns
    end do

    jrank=map_grid_2_pinfo(a1b,a2b,a3b-1,1)
    ip=a3b-1
    fp=ip
    nc=1
    ns=(b2b-a2b+1)*(b1b-a1b+1)
    fdinfo_recv(1,nc,5)=a1b
    fdinfo_recv(2,nc,5)=b1b
    fdinfo_recv(3,nc,5)=a2b
    fdinfo_recv(4,nc,5)=b2b
    fdinfo_recv(5,nc,5)=ip
    fdinfo_recv(6,nc,5)=fp
    fdinfo_recv(7,nc,5)=myrank_g
    fdinfo_recv(8,nc,5)=jrank
    fdinfo_recv(9,nc,5)=ns
    fdinfo_recv(10,nc,5)=ns
    do i=a3b-2,a3b-Md,-1
       irank=map_grid_2_pinfo(a1b,a2b,i,1)
       if ( irank/=jrank ) then
          jrank=irank
          nc=nc+1
          fp=i
       end if
       fdinfo_recv(1,nc,5)=a1b
       fdinfo_recv(2,nc,5)=b1b
       fdinfo_recv(3,nc,5)=a2b
       fdinfo_recv(4,nc,5)=b2b
       fdinfo_recv(5,nc,5)=i
       fdinfo_recv(6,nc,5)=fp
       fdinfo_recv(7,nc,5)=myrank_g
       fdinfo_recv(8,nc,5)=jrank
       fdinfo_recv(9,nc,5)=(fp-i+1)*ns
    end do

    jrank=map_grid_2_pinfo(a1b,a2b,b3b+1,1)
    ip=b3b+1
    fp=ip
    nc=1
    ns=(b2b-a2b+1)*(b1b-a1b+1)
    fdinfo_recv(1,nc,6)=a1b
    fdinfo_recv(2,nc,6)=b1b
    fdinfo_recv(3,nc,6)=a2b
    fdinfo_recv(4,nc,6)=b2b
    fdinfo_recv(5,nc,6)=ip
    fdinfo_recv(6,nc,6)=fp
    fdinfo_recv(7,nc,6)=myrank_g
    fdinfo_recv(8,nc,6)=jrank
    fdinfo_recv(9,nc,6)=ns
    fdinfo_recv(10,nc,6)=ns
    do i=b3b+1,b3b+Md
       irank=map_grid_2_pinfo(a1b,a2b,i,1)
       if ( irank/=jrank ) then
          jrank=irank
          nc=nc+1
          ip=i
       end if
       fdinfo_recv(1,nc,6)=a1b
       fdinfo_recv(2,nc,6)=b1b
       fdinfo_recv(3,nc,6)=a2b
       fdinfo_recv(4,nc,6)=b2b
       fdinfo_recv(5,nc,6)=ip
       fdinfo_recv(6,nc,6)=i
       fdinfo_recv(7,nc,6)=myrank_g
       fdinfo_recv(8,nc,6)=jrank
       fdinfo_recv(9,nc,6)=(i-ip+1)*ns
    end do

    deallocate( map_grid_2_pinfo )

    n=maxval(n_neighbor)*6*2
    allocate( ireq(n) ) ; ireq(:)=0

    nreq=0
    do j=1,6
       do i=1,n_neighbor(j)
          irank=fdinfo_recv(8,i,j)
          select case(j)
          case(1,3,5)
             itags=10*(j+1)+i
             itagr=10*j+i
          case(2,4,6)
             itags=10*(j-1)+i
             itagr=10*j+i
          end select
          nreq=nreq+1
          call mpi_irecv(fdinfo_send(1,i,j),10,mpi_integer,irank &
                        ,itagr,comm_grid,ireq(nreq),ierr)
          nreq=nreq+1
          call mpi_isend(fdinfo_recv(1,i,j),10,mpi_integer,irank &
                        ,itags,comm_grid,ireq(nreq),ierr)
       end do
    end do

    call mpi_waitall(nreq,ireq,istatus,ierr)

    deallocate( ireq )

    do i=1,6
       do n=1,n_neighbor(i)
          select case(i)
          case(1,2)
             fdinfo_send(1,n,i)=mod(fdinfo_send(1,n,i)+ML1,ML1)
             fdinfo_send(2,n,i)=mod(fdinfo_send(2,n,i)+ML1,ML1)
          case(3,4)
             fdinfo_send(3,n,i)=mod(fdinfo_send(3,n,i)+ML2,ML2)
             fdinfo_send(4,n,i)=mod(fdinfo_send(4,n,i)+ML2,ML2)
          case(5,6)
             fdinfo_send(5,n,i)=mod(fdinfo_send(5,n,i)+ML3,ML3)
             fdinfo_send(6,n,i)=mod(fdinfo_send(6,n,i)+ML3,ML3)
          end select
       end do
    end do

    call allocate_bcset

    call write_border( 80, " init_bcset_sol(end)" )

  END SUBROUTINE init_bcset_sol


  SUBROUTINE allocate_bcset
    implicit none
    integer :: a1,a2,a3,b1,b2,b3
    a1=Igrid(1,1)-Md ; a2=Igrid(1,2)-Md ; a3=Igrid(1,3)-Md
    b1=Igrid(2,1)+Md ; b2=Igrid(2,2)+Md ; b3=Igrid(2,3)+Md
    if ( allocated(www) ) deallocate(www)
    allocate( www(a1:b1,a2:b2,a3:b3,MB_d) )
    www(:,:,:,:)=zero
    a1=maxval( n_neighbor(1:6) )
    a2=maxval( fdinfo_send(9,1:a1,1:6) )*MB_d
    a3=maxval( fdinfo_recv(9,1:a1,1:6) )*MB_d
    if ( allocated(sbuf) ) deallocate(sbuf)
    if ( allocated(rbuf) ) deallocate(rbuf)
    allocate( sbuf(a2,a1,6) )
    allocate( rbuf(a3,a1,6) )
    sbuf(:,:,:)=zero
    rbuf(:,:,:)=zero
  END SUBROUTINE allocate_bcset


END MODULE bc_module
