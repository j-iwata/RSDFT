MODULE parallel_module

  implicit none

  include 'mpif.h'

!  PRIVATE
  PRIVATE :: send_parallel
  PUBLIC :: node_partition &
           ,comm_grid, comm_band, comm_bzsm, comm_spin &
           ,myrank,myrank_g,myrank_b,myrank_k,myrank_s &
           ,nprocs,nprocs_g,nprocs_k,nprocs_s &
           ,np_band,np_spin,np_bzsm,np_grid &
           ,id_class,ircnt,idisp,ir_spin,id_spin &
           ,ir_band,id_band,ir_grid,id_grid,ir_bzsm,id_bzsm &
           ,pinfo_grid,MB_d &
           ,read_parallel,init_parallel &
           ,start_mpi_parallel,end_mpi_parallel &
           ,disp_switch_parallel
  PUBLIC :: comm_fkmb, myrank_f, np_fkmb, ir_fkmb, id_fkmb

  integer,parameter :: max_parallel=7
  integer :: node_partition(max_parallel)
  integer :: nprocs, myrank
  integer :: comm_grid, myrank_g, nprocs_g, np_grid
  integer :: comm_spin, myrank_s, nprocs_s, np_spin
  integer :: comm_band, myrank_b, nprocs_b, np_band
  integer :: comm_bzsm, myrank_k, nprocs_k, np_bzsm
  integer :: comm_fkmb, myrank_f, nprocs_f, np_fkmb
  integer :: MB_d
  integer,allocatable :: id_class(:,:),ircnt(:),idisp(:)
  integer,allocatable :: ir_grid(:),id_grid(:)
  integer,allocatable :: ir_band(:),id_band(:)
  integer,allocatable :: ir_bzsm(:),id_bzsm(:)
  integer,allocatable :: ir_spin(:),id_spin(:)
  integer,allocatable :: ir_fkmb(:),id_fkmb(:)
  integer,allocatable :: pinfo_grid(:,:)
  logical :: disp_switch_parallel=.false.

CONTAINS


  SUBROUTINE start_mpi_parallel
    integer :: ierr,iprovided
    call mpi_init(ierr)
    !call mpi_init_thread( MPI_THREAD_MULTIPLE, iprovided, ierr)
    call mpi_comm_size(mpi_comm_world,nprocs,ierr)
    call mpi_comm_rank(mpi_comm_world,myrank,ierr)
  END SUBROUTINE start_mpi_parallel


  SUBROUTINE end_mpi_parallel
    integer :: ierr
    call mpi_finalize(ierr)
  END SUBROUTINE end_mpi_parallel


  SUBROUTINE read_parallel(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i
    character(5) :: cbuf,ckey
    node_partition(:)=1
    MB_d = 0
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:5) == "PROCS" ) then
             backspace(unit)
             read(unit,*) cbuf,node_partition(1:max_parallel)
          else if ( ckey(1:3) == "MBD" ) then
             backspace(unit)
             read(unit,*) cbuf,MB_d
          end if
       end do
999    continue
       write(*,'(1x,"node_partition(1:6)=",6i4)') node_partition(1:6)
       write(*,*) "MB_d=",MB_d
    end if
    call send_parallel(0)
  END SUBROUTINE read_parallel


  SUBROUTINE read_oldformat_parallel(rank,unit)
    integer,intent(IN) :: rank,unit
    node_partition(:)=1
    if ( rank == 0 ) then
       read(unit,*) node_partition(1:max_parallel)
       read(unit,*) MB_d
       write(*,'(1x,"node_partition(1:6)=",9i4)') &
            node_partition(1:max_parallel)
       write(*,*) "MB_d=",MB_d
    end if
    call send_parallel(0)
  END SUBROUTINE read_oldformat_parallel


  SUBROUTINE send_parallel(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    call mpi_bcast(node_partition,max_parallel,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(MB_d,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_parallel


  SUBROUTINE init_parallel(Ngrid,Nband,Nbzsm,Nspin)
    implicit none
    integer,intent(IN) :: Ngrid(0:3),Nband,Nspin,Nbzsm
    integer :: m,n,ierr,irank,jrank,itags,itagr,nreq
    integer :: istatus(MPI_STATUS_SIZE,123)
    integer :: ip,fp,nc,ns,a1,a2,a3,a4,a5,a6,a7,b1,b2,b3,b4,b5,b6
    integer :: i1,i2,i3,i4,ib,j0,j1,j2,j3,m0,m1,m2,m3,i,j,k,id,ie,ir
    integer :: ML1,ML2,ML3,np1,np2,np3,MLI(2,3)
    integer,allocatable :: mtmp(:),ntmp(:,:),ireq(:)
    integer,allocatable :: map_grid_2_pinfo(:,:,:,:)

    call write_border( 80, " init_parallel(start)" )

    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    n=1
    do i=1,max_parallel
       n=n*node_partition(i)
    end do

    if ( n /= nprocs ) then
       write(*,*) "n,nprocs=",n,nprocs
       stop
    end if

    np1=node_partition(1)
    np2=node_partition(2)
    np3=node_partition(3)
    if ( disp_switch_parallel ) then
       write(*,*) "np1,np2,np3=",np1,np2,np3
       write(*,*) "np_grid    =",np1*np2*np3
       write(*,*) "np_band    =",node_partition(4)
       write(*,*) "np_bz      =",node_partition(5)
       write(*,*) "np_spin    =",node_partition(6)
       write(*,*) "np_fkmb    =",node_partition(7)
       write(*,*) "nprocs     =",nprocs
    end if
    if ( np1>ML1 .or. np2>ML2 .or. np3>ML3 .or. &
         node_partition(4)>Nband .or. node_partition(6)>nspin ) then
       write(*,'(1x,12i6)') node_partition(1:6),ML1,ML2,ML3,Nband,nspin
       stop "stop@init_parallel"
    end if
    if ( node_partition(5)>Nbzsm ) then
       write(*,'(1x,i6)') Nbzsm
       stop
    end if

! --- class ---

    allocate( id_class(0:nprocs-1,0:max_parallel) )
    id_class=-1

    j=-1
    do a7=0,node_partition(7)-1
    do a6=0,node_partition(6)-1
    do a5=0,node_partition(5)-1
    do a4=0,node_partition(4)-1
       i=-1
       do a3=0,node_partition(3)-1
       do a2=0,node_partition(2)-1
       do a1=0,node_partition(1)-1
          i=i+1
          j=j+1
          id_class(j,0)=i
          id_class(j,1)=a1
          id_class(j,2)=a2
          id_class(j,3)=a3
          id_class(j,4)=a4
          id_class(j,5)=a5
          id_class(j,6)=a6
          id_class(j,7)=a7
       end do ! a1
       end do ! a2
       end do ! a3
    end do ! a4
    end do ! a5
    end do ! a6
    end do ! a7

    if ( disp_switch_parallel ) then
       write(*,'(1x,5a5)') "irank","grid","band","bz","spin"
       do i=0,nprocs-1
          write(*,'(1x,6(2x,i3))') i,id_class(i,0),id_class(i,4:7)
       end do
    end if

! --- fock-band parallel ---

    if ( disp_switch_parallel ) write(*,'("--- fock-band-parallel ---")')

    np_fkmb = node_partition(7)

    if ( np_fkmb > Nband ) then
       write(*,*) "np_fkmb is too large!"
       write(*,*) ",np_fkmb, Nband=",np_fkmb,Nband
       stop
    end if

    allocate( id_fkmb(0:np_fkmb-1) ) ; id_fkmb=0
    allocate( ir_fkmb(0:np_fkmb-1) ) ; ir_fkmb=0

    do i=0,Nband-1
       j=mod(i,np_fkmb)
       ir_fkmb(j)=ir_fkmb(j)+1
    end do

    do i=0,np_fkmb-1
       id_fkmb(i)=sum( ir_fkmb(0:i) )-ir_fkmb(i)
    end do

    if ( disp_switch_parallel ) then
       write(*,*) "Nband  =",Nband
       write(*,*) "np_fkmb=",np_fkmb
       write(*,'(1x,3a8)') "i","id_fkmb","ir_fkmb"
       do i=0,np_fkmb-1
          write(*,'(1x,3i8)') i,id_fkmb(i),ir_fkmb(i)
       end do
    end if

    nprocs_f = count( id_class(:,7)==id_class(myrank,7) )

    if ( disp_switch_parallel ) write(*,*) "nprocs_f =",nprocs_f

! --- spin parallel ---

    if ( disp_switch_parallel ) write(*,'("--- spin-parallel ---")')

    np_spin = node_partition(6)

    if ( np_spin>nspin ) then
       write(*,*) "np_spin is too large!"
       write(*,*) ",np_spin,nspin=",np_spin,nspin
       stop
    end if

    allocate( id_spin(0:np_spin-1) ) ; id_spin=0
    allocate( ir_spin(0:np_spin-1) ) ; ir_spin=0

    do i=0,nspin-1
       j=mod(i,np_spin)
       ir_spin(j)=ir_spin(j)+1
    end do

    do i=0,np_spin-1
       id_spin(i)=sum( ir_spin(0:i) )-ir_spin(i)
    end do

    if ( disp_switch_parallel ) then
       write(*,*) "nspin  =",nspin
       write(*,*) "np_spin=",np_spin
       write(*,'(1x,3a8)') "i","id_spin","ir_spin"
       do i=0,np_spin-1
          write(*,'(1x,3i8)') i,id_spin(i),ir_spin(i)
       end do
    end if

    nprocs_s = count( id_class(:,6)==id_class(myrank,6) .and. &
                      id_class(:,7)==id_class(myrank,7) )

    if ( disp_switch_parallel ) write(*,*) "nprocs_s =",nprocs_s

! --- k parallel ---

    if ( disp_switch_parallel ) write(*,'("--- k-parallel ---")')

    np_bzsm = node_partition(5)

    allocate( id_bzsm(0:np_bzsm-1) ) ; id_bzsm=0
    allocate( ir_bzsm(0:np_bzsm-1) ) ; ir_bzsm=0

    do i=0,max(Nbzsm,np_bzsm)-1
       j=mod(i,np_bzsm)
       ir_bzsm(j)=ir_bzsm(j)+1
    end do

    do i=0,np_bzsm-1
       id_bzsm(i)=sum( ir_bzsm(0:i) )-ir_bzsm(i)
    end do

    if ( disp_switch_parallel ) then
       write(*,*) "Nbzsm    =",Nbzsm
       write(*,*) "np_bzsm=",np_bzsm
       write(*,'(1x,3a8)') "i","id_bzsm","ir_bzsm"
       do i=0,np_bzsm-1
          write(*,'(1x,3i8)') i,id_bzsm(i),ir_bzsm(i)
       end do
    end if

    nprocs_k = count( id_class(:,5)==id_class(myrank,5) .and. &
                      id_class(:,6)==id_class(myrank,6) .and. &
                      id_class(:,7)==id_class(myrank,7) )

    if ( disp_switch_parallel ) write(*,*) "nprocs_k =",nprocs_k

! --- band parallel ---

    if ( disp_switch_parallel ) write(*,'("--- band-parallel ---")')

    np_band = node_partition(4)

    allocate( id_band(0:np_band-1) ) ; id_band=0
    allocate( ir_band(0:np_band-1) ) ; ir_band=0

    do i=0,Nband-1
       j=mod(i,np_band)
       ir_band(j)=ir_band(j)+1
    end do

    do i=0,np_band-1
       id_band(i)=sum( ir_band(0:i) )-ir_band(i)
    end do

    if ( disp_switch_parallel ) then
       write(*,*) "Nband  =",Nband
       write(*,*) "np_band=",np_band
       write(*,'(1x,3a8)') "i","id_band","ir_band"
       do i=0,np_band-1
          write(*,'(1x,3i8)') i,id_band(i),ir_band(i)
       end do
    end if

    nprocs_b = count( id_class(:,4)==id_class(myrank,4) .and.  &
                      id_class(:,5)==id_class(myrank,5) .and.  &
                      id_class(:,6)==id_class(myrank,6) .and.  &
                      id_class(:,7)==id_class(myrank,7)       )

    if ( disp_switch_parallel ) write(*,*) "nprocs_b =",nprocs_b

! --- grid parallel ---

    if ( disp_switch_parallel ) write(*,'("--- grid-parallel ---")')

    np_grid = node_partition(1)*node_partition(2)*node_partition(3)
    if ( np_spin*np_bzsm*np_band*np_grid*np_fkmb /= nprocs ) stop "init_parallel"

    allocate( id_grid(0:np_grid-1)      ) ; id_grid=0
    allocate( ir_grid(0:np_grid-1)      ) ; ir_grid=0
    allocate( pinfo_grid(8,0:np_grid-1) ) ; pinfo_grid=0

    nprocs_g = count( id_class(:,0)==id_class(myrank,0) .and. &
                      id_class(:,4)==id_class(myrank,4) .and. &
                      id_class(:,5)==id_class(myrank,5) .and. &
                      id_class(:,6)==id_class(myrank,6) .and. &
                      id_class(:,7)==id_class(myrank,7)       )

! ---

    if ( disp_switch_parallel ) then
       write(*,*) "nprocs_g =",nprocs_g
       write(*,*) "nprocs_b =",nprocs_b
       write(*,*) "nprocs_k =",nprocs_k
       write(*,*) "nprocs_s =",nprocs_s
       write(*,*) "nprocs_f =",nprocs_f
       write(*,*) "nprocs   =",nprocs
    end if

! --- communicators ---

    m=0
    i=0
    do i4=0,np_fkmb-1
    do i3=0,np_spin-1
    do i2=0,np_bzsm-1
    do i1=0,np_band-1
       i=i+1
       if ( id_class(myrank,4)==i1 .and. &
            id_class(myrank,5)==i2 .and. &
            id_class(myrank,6)==i3 .and. &
            id_class(myrank,7)==i4       ) m=i
    end do
    end do
    end do
    end do
    call mpi_comm_split(mpi_comm_world,m,myrank,comm_grid,ierr)
    call mpi_comm_rank(comm_grid,myrank_g,ierr)
    call mpi_comm_size(comm_grid,nprocs_g,ierr)

    m=0
    i=0
    do i4=0,np_fkmb-1
    do i3=0,np_spin-1
    do i2=0,np_bzsm-1
    do i1=0,np_grid-1
       i=i+1
       if ( id_class(myrank,0)==i1 .and. &
            id_class(myrank,5)==i2 .and. &
            id_class(myrank,6)==i3 .and. &
            id_class(myrank,7)==i4       ) m=i
    end do
    end do
    end do
    end do
    call mpi_comm_split(mpi_comm_world,m,myrank,comm_band,ierr)
    call mpi_comm_rank(comm_band,myrank_b,ierr)
    call mpi_comm_size(comm_band,nprocs_b,ierr)

    m=0
    i=0
    do i4=0,np_fkmb-1
    do i3=0,np_spin-1
    do i2=0,np_band-1
    do i1=0,np_grid-1
       i=i+1
       if ( id_class(myrank,0)==i1 .and. &
            id_class(myrank,4)==i2 .and. &
            id_class(myrank,6)==i3 .and. &
            id_class(myrank,7)==i4 ) m=i
    end do
    end do
    end do
    end do
    call mpi_comm_split(mpi_comm_world,m,myrank,comm_bzsm,ierr)
    call mpi_comm_rank(comm_bzsm,myrank_k,ierr)
    call mpi_comm_size(comm_bzsm,nprocs_k,ierr)

    m=0
    i=0
    do i4=0,np_fkmb-1
    do i3=0,np_bzsm-1
    do i2=0,np_band-1
    do i1=0,np_grid-1
       i=i+1
       if ( id_class(myrank,0)==i1 .and. &
            id_class(myrank,4)==i2 .and. &
            id_class(myrank,5)==i3 .and. &
            id_class(myrank,7)==i4 ) m=i
    end do
    end do
    end do
    end do
    call mpi_comm_split(mpi_comm_world,m,myrank,comm_spin,ierr)
    call mpi_comm_rank(comm_spin,myrank_s,ierr)
    call mpi_comm_size(comm_spin,nprocs_s,ierr)

    m=0
    i=0
    do i4=0,np_spin-1
    do i3=0,np_bzsm-1
    do i2=0,np_band-1
    do i1=0,np_grid-1
       i=i+1
       if ( id_class(myrank,0)==i1 .and. &
            id_class(myrank,4)==i2 .and. &
            id_class(myrank,5)==i3 .and. &
            id_class(myrank,6)==i4 ) m=i
    end do
    end do
    end do
    end do
    call mpi_comm_split(mpi_comm_world,m,myrank,comm_fkmb,ierr)
    call mpi_comm_rank(comm_fkmb,myrank_f,ierr)
    call mpi_comm_size(comm_fkmb,nprocs_f,ierr)

    if ( disp_switch_parallel ) then
       write(*,*) "comm_world, nprocs   =",MPI_COMM_WORLD,nprocs
       write(*,*) "comm_grid,  nprocs_g =",comm_grid,nprocs_g
       write(*,*) "comm_band,  nprocs_b =",comm_band,nprocs_b
       write(*,*) "comm_bzsm,  nprocs_k =",comm_bzsm,nprocs_k
       write(*,*) "comm_spin,  nprocs_s =",comm_spin,nprocs_s
       write(*,*) "comm_fkmb,  nprocs_f =",comm_fkmb,nprocs_f
    end if

    if ( myrank_g /= id_class(myrank,0) ) then
       write(*,*) "myrank,myrank_g=",myrank,myrank_g
       stop
    end if

! --- band bundle ---

    n = min( MB_d,ir_band(id_class(myrank,4)) )
    MB_d = max( n,1 )

! ---

    allocate( idisp(0:nprocs-1) ) ; idisp=-1
    allocate( ircnt(0:nprocs-1) ) ; ircnt=0

    call write_border( 80, " init_parallel(end)" )

  END SUBROUTINE init_parallel


END MODULE parallel_module
