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

  integer,parameter :: max_parallel=6
  integer :: node_partition(max_parallel)
  integer :: comm_grid, comm_band, comm_bzsm, comm_spin
  integer :: myrank,myrank_g,myrank_b,myrank_k,myrank_s
  integer :: nprocs,nprocs_g,nprocs_b,nprocs_k,nprocs_s
  integer :: np_band,np_grid,np_spin,np_bzsm
  integer :: MB_d
  integer,allocatable :: id_class(:,:),ircnt(:),idisp(:)
  integer,allocatable :: ir_band(:),id_band(:)
  integer,allocatable :: ir_bzsm(:),id_bzsm(:)
  integer,allocatable :: ir_spin(:),id_spin(:)
  integer,allocatable :: ir_grid(:),id_grid(:)
  integer,allocatable :: pinfo_grid(:,:)
  logical :: disp_switch_parallel=.false.

CONTAINS


  SUBROUTINE start_mpi_parallel
    integer :: ierr,i
    !call mpi_init(ierr)
    call mpi_init_thread( MPI_THREAD_MULTIPLE, i, ierr)
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
       write(*,'(1x,"node_partition(1:6)=",6i4)') node_partition(1:6)
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
    integer :: ip,fp,nc,ns,a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6
    integer :: i1,i2,i3,ib,j0,j1,j2,j3,m0,m1,m2,m3,i,j,k,id,ie,ir
    integer :: ML1,ML2,ML3,np1,np2,np3,MB,MLI(2,3)
    integer,allocatable :: mtmp(:),ntmp(:,:),ireq(:)
    integer,allocatable :: map_grid_2_pinfo(:,:,:,:)

    if ( disp_switch_parallel ) write(*,'(a60," init_parallel")') repeat("-",60)

    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)
    MB  = Nband

    n=1
    do i=1,6
       n=n*node_partition(i)
    end do

    if ( n/=nprocs ) then
       write(*,*) "n,nprocs=",n,nprocs
       stop
    end if

    np1=node_partition(1)
    np2=node_partition(2)
    np3=node_partition(3)
    if ( disp_switch_parallel ) then
       write(*,*) "np1,np2,np3=",np1,np2,np3
       write(*,*) "np_band    =",node_partition(4)
       write(*,*) "np_bz      =",node_partition(5)
       write(*,*) "np_spin    =",node_partition(6)
       write(*,*) "np_grid    =",np1*np2*np3
       write(*,*) "nprocs     =",nprocs
    end if
    if ( np1>ML1 .or. np2>ML2 .or. np3>ML3 .or. &
         node_partition(4)>Nband .or. node_partition(6)>nspin ) then
       write(*,'(1x,12i6)') node_partition(1:6),ML1,ML2,ML3,MB,nspin
       stop "stop@init_parallel"
    end if
    if ( node_partition(5)>Nbzsm ) then
       write(*,'(1x,i6)') Nbzsm
       stop
    end if

! --- class ---

    allocate( id_class(0:nprocs-1,0:6) )
    id_class=-1

    j=-1
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
       end do
       end do
       end do
    end do
    end do
    end do

!    if ( disp_switch_parallel ) then
!       write(*,'(1x,5a5)') "irank","grid","band","bz","spin"
!       do i=0,nprocs-1
!          write(*,'(1x,5(2x,i3))') i,id_class(i,0),id_class(i,4:6)
!       end do
!    end if

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

    nprocs_s = count( id_class(:,6)==id_class(myrank,6) )

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
                      id_class(:,6)==id_class(myrank,6) )

    if ( disp_switch_parallel ) write(*,*) "nprocs_k =",nprocs_k

! --- band parallel ---

    if ( disp_switch_parallel ) write(*,'("--- band-parallel ---")')

    np_band = node_partition(4)

    allocate( id_band(0:np_band-1) ) ; id_band=0
    allocate( ir_band(0:np_band-1) ) ; ir_band=0

    do i=0,MB-1
       j=mod(i,np_band)
       ir_band(j)=ir_band(j)+1
    end do

    do i=0,np_band-1
       id_band(i)=sum( ir_band(0:i) )-ir_band(i)
    end do

    if ( disp_switch_parallel ) then
       write(*,*) "MB     =",MB
       write(*,*) "np_band=",np_band
       write(*,'(1x,3a8)') "i","id_band","ir_band"
       do i=0,np_band-1
          write(*,'(1x,3i8)') i,id_band(i),ir_band(i)
       end do
    end if

    nprocs_b = count( id_class(:,4)==id_class(myrank,4) .and.  &
                      id_class(:,5)==id_class(myrank,5) .and.  &
                      id_class(:,6)==id_class(myrank,6)       )

    if ( disp_switch_parallel ) write(*,*) "nprocs_b =",nprocs_b

! --- grid parallel ---

    if ( disp_switch_parallel ) write(*,'("--- grid-parallel ---")')

    np_grid = node_partition(1)*node_partition(2)*node_partition(3)
    if ( np_spin*np_bzsm*np_band*np_grid/=nprocs ) stop

    allocate( id_grid(0:np_grid-1)      ) ; id_grid=0
    allocate( ir_grid(0:np_grid-1)      ) ; ir_grid=0
    allocate( pinfo_grid(8,0:np_grid-1) ) ; pinfo_grid=0

    nprocs_g = count( id_class(:,0)==id_class(myrank,0) .and. &
                      id_class(:,4)==id_class(myrank,4) .and. &
                      id_class(:,5)==id_class(myrank,5) .and. &
                      id_class(:,6)==id_class(myrank,6)       )

    if ( disp_switch_parallel ) then
       write(*,*) "nprocs_g =",nprocs_g
       write(*,*) "nprocs_b =",nprocs_b
       write(*,*) "nprocs_k =",nprocs_k
       write(*,*) "nprocs_s =",nprocs_s
       write(*,*) "nprocs   =",nprocs
    end if

! --- communicators ---

    allocate( ntmp(0:nprocs-1,4) ) ; ntmp=0

    m=0
    i=0
    do i3=0,np_spin-1
    do i2=0,np_bzsm-1
    do i1=0,np_band-1
       i=i+1
       if ( id_class(myrank,4)==i1 .and. &
            id_class(myrank,5)==i2 .and. &
            id_class(myrank,6)==i3       ) m=i
    end do
    end do
    end do
    call mpi_comm_split(mpi_comm_world,m,myrank,comm_grid,ierr)
    call mpi_comm_rank(comm_grid,myrank_g,ierr)
    call mpi_comm_size(comm_grid,nprocs_g,ierr)

    m=0
    i=0
    do i3=0,np_spin-1
    do i2=0,np_bzsm-1
    do i1=0,np_grid-1
       i=i+1
       if ( id_class(myrank,0)==i1 .and. &
            id_class(myrank,5)==i2 .and. &
            id_class(myrank,6)==i3       ) m=i
    end do
    end do
    end do
    call mpi_comm_split(mpi_comm_world,m,myrank,comm_band,ierr)
    call mpi_comm_rank(comm_band,myrank_b,ierr)
    call mpi_comm_size(comm_band,nprocs_b,ierr)

    m=0
    i=0
    do i3=0,np_spin-1
    do i2=0,np_band-1
    do i1=0,np_grid-1
       i=i+1
       if ( id_class(myrank,0)==i1 .and. &
            id_class(myrank,4)==i2 .and. id_class(myrank,6)==i3 ) m=i
    end do
    end do
    end do
    call mpi_comm_split(mpi_comm_world,m,myrank,comm_bzsm,ierr)
    call mpi_comm_rank(comm_bzsm,myrank_k,ierr)
    call mpi_comm_size(comm_bzsm,nprocs_k,ierr)

    m=0
    i=0
    do i3=0,np_bzsm-1
    do i2=0,np_band-1
    do i1=0,np_grid-1
       i=i+1
       if ( id_class(myrank,0)==i1 .and. &
            id_class(myrank,4)==i2 .and. id_class(myrank,5)==i3 ) m=i
    end do
    end do
    end do
    call mpi_comm_split(mpi_comm_world,m,myrank,comm_spin,ierr)
    call mpi_comm_rank(comm_spin,myrank_s,ierr)
    call mpi_comm_size(comm_spin,nprocs_s,ierr)

    deallocate( ntmp )

    if ( disp_switch_parallel ) then
       write(*,*) "comm_world, nprocs   =",MPI_COMM_WORLD,nprocs
       write(*,*) "comm_grid,  nprocs_g =",comm_grid,nprocs_g
       write(*,*) "comm_band,  nprocs_b =",comm_band,nprocs_b
       write(*,*) "comm_bzsm,  nprocs_k =",comm_bzsm,nprocs_k
       write(*,*) "comm_spin,  nprocs_s =",comm_spin,nprocs_s
    end if

    if ( myrank_g/=id_class(myrank,0) ) then
       write(*,*) "myrank,myrank_g=",myrank,myrank_g
       stop
    end if

! --- band bundle ---

    n = min( MB_d,ir_band(id_class(myrank,4)) )
    MB_d = max( n,1 )

! ---

    allocate( idisp(0:nprocs-1) ) ; idisp=-1
    allocate( ircnt(0:nprocs-1) ) ; ircnt=0

  END SUBROUTINE init_parallel


END MODULE parallel_module



SUBROUTINE convert_capital(cbuf,CKEY)
  implicit none
  character(*),intent(IN)  :: cbuf
  character(*),intent(OUT) :: CKEY
  integer :: j,k,n
  n=len_trim(cbuf)
  CKEY=""
  do j=1,n
     k=iachar( cbuf(j:j) )
     if ( k >= 97 ) k=k-32
     CKEY(j:j) = achar(k)
  end do
END SUBROUTINE convert_capital



PROGRAM cubegen_wf3

  use parallel_module

  implicit none

  integer,parameter :: u1=10,u2=970,u3=3,iu=2
  integer :: ML,ML1,ML2,ML3,zatom(9),Ngrid(0:3)
  integer :: i1,i2,i3,i,j,n,k,s,MI,MKI,ierr
  integer :: n1,n2,k1,k2,s1,s2,b1,b2
  integer,allocatable :: LL(:,:),Kion(:)
  real(8),allocatable :: urtmp(:,:,:,:),rho(:,:,:),occ(:,:,:)
  real(8),allocatable :: asi(:,:),rsi(:,:),rtmp(:)
  real(8) :: ax,aa(3,3)
  character(30) :: cbuf,file_name,file_wf_split
  character(5) :: cmyrank
  complex(8),allocatable :: uztmp(:,:,:,:),ztmp(:)

  integer :: MB,MB1,MB2
  integer :: MBZ,MSP,itmp(6)
  integer :: ndata, iflag_abs
  integer,allocatable :: idata(:,:)
  logical :: flag_real8

! ---

  call start_mpi_parallel

  disp_switch_parallel = ( myrank == 0 )

! ---

  if ( myrank == 0 ) then
  rewind u3
  read(u3,*) cbuf, iflag_abs
  read(u3,*) file_name
  read(u3,*) MBZ
  read(u3,*) MSP
  read(u3,*) ndata
  end if

  call mpi_bcast( cbuf,30,mpi_character,0,mpi_comm_world,ierr )
  call mpi_bcast( iflag_abs,1,mpi_integer,0,mpi_comm_world,ierr )
  call mpi_bcast( file_name,30,mpi_character,0,mpi_comm_world,ierr )
  call mpi_bcast( iflag_abs,1,mpi_integer,0,mpi_comm_world,ierr )
  call mpi_bcast( MBZ,1,mpi_integer,0,mpi_comm_world,ierr )
  call mpi_bcast( MSP,1,mpi_integer,0,mpi_comm_world,ierr )
  call mpi_bcast( ndata,1,mpi_integer,0,mpi_comm_world,ierr )

! ---

  if ( cbuf == "r" .or. cbuf == "R" ) then
     flag_real8 = .true.
  else if ( cbuf == "c" .or. cbuf == "C" ) then
     flag_real8 = .false.
     if ( iflag_abs == 0 ) then
        stop "complex wave function should be converted to a real quantity"
     end if
  end if

  if ( myrank == 0 ) then
  write(*,*) "flag_real8=",flag_real8, cbuf
  write(*,*) "iflag_abs =",iflag_abs
  write(*,*) "file_name =",file_name
  write(*,*) "MBZ, MSP  =",MBZ, MSP
  write(*,*) "ndata     =",ndata
  end if

  allocate( idata(3,ndata) ) ; idata=0

  if ( myrank == 0 ) then
  do i=1,ndata
     read(u3,*) idata(1:3,i)
  end do
  end if

  call mpi_bcast(idata,3*ndata,mpi_integer,0,mpi_comm_world,ierr)

! ---

  call read_parallel( myrank, u3 )

! ---

  if ( myrank == 0 ) then
  read(u2,*) cbuf,ax
  read(u2,*) cbuf,aa(1:3,1)
  read(u2,*) cbuf,aa(1:3,2)
  read(u2,*) cbuf,aa(1:3,3)
  read(u2,*)
  read(u2,*) MKI,MI,zatom(1:MKI)
  end if

  call mpi_bcast( ax,1,mpi_real8,0,mpi_comm_world,ierr )
  call mpi_bcast( aa,9,mpi_real8,0,mpi_comm_world,ierr )
  call mpi_bcast( MKI,1,mpi_integer,0,mpi_comm_world,ierr )
  call mpi_bcast( MI,1,mpi_integer,0,mpi_comm_world,ierr )
  call mpi_bcast( zatom,9,mpi_integer,0,mpi_comm_world,ierr )

  allocate( asi(3,MI) ) ; asi=0.0d0
  allocate( rsi(3,MI) ) ; asi=0.0d0
  allocate( Kion(MI)  ) ; Kion=0

  if ( myrank == 0 ) then
  do i=1,MI
     read(u2,*) Kion(i),asi(1:3,i)
  end do
  end if

  call mpi_bcast( asi,3*MI,mpi_real8,0,mpi_comm_world,ierr )
  call mpi_bcast( Kion,MI,mpi_integer,0,mpi_comm_world,ierr )

! ---

  aa(:,:)=ax*aa(:,:)

  rsi(:,:)=matmul( aa,asi )

! ---

  if ( myrank == 0 ) then
     open(u1,file=file_name,status='old',form='unformatted')
     read(u1) i
     if ( i <= 0 ) then
        write(*,*) "--- new format ---"
        read(u1) ML,ML1,ML2,ML3
        read(u1) MB,MB1,MB2
        read(u1) itmp(1:3)
        read(u1) itmp(4:6)
     else
        rewind u1
        read(u1) ML,ML1,ML2,ML3
        read(u1) MB,MB1,MB2
     end if
  endif

  call mpi_bcast(ML,1,mpi_integer,0,mpi_comm_world,ierr)
  call mpi_bcast(ML1,1,mpi_integer,0,mpi_comm_world,ierr)
  call mpi_bcast(ML2,1,mpi_integer,0,mpi_comm_world,ierr)
  call mpi_bcast(ML3,1,mpi_integer,0,mpi_comm_world,ierr)
  call mpi_bcast(MB,1,mpi_integer,0,mpi_comm_world,ierr)
  call mpi_bcast(MB1,1,mpi_integer,0,mpi_comm_world,ierr)
  call mpi_bcast(MB2,1,mpi_integer,0,mpi_comm_world,ierr)

  allocate( LL(3,ML)        ) ; LL=0.0d0
  allocate( occ(MB,MBZ,MSP) ) ; occ=0.0d0

  if ( myrank == 0 ) then
  read(u1) LL
  read(u1) occ
  close(u1)
  end if

  call mpi_bcast(LL,3*ML,mpi_integer,0,mpi_comm_world,ierr)
  call mpi_bcast(occ,size(occ),mpi_integer,0,mpi_comm_world,ierr)

  if ( myrank == 0 ) then
     write(*,*) "ML,ML1,ML2,ML3=",ML,ML1,ML2,ML3
     write(*,*) "MB=",MB
  end if

! ---

  Ngrid(0) = ML
  Ngrid(1) = ML1
  Ngrid(2) = ML2
  Ngrid(3) = ML3

  call init_parallel( Ngrid, MB, MBZ, MSP )

! ---

!  ir_grid(:)=Ngrid(0)/np_grid
!    n=Ngrid(0)-sum(ir_grid)
!  do i=1,n
!     j=mod(i-1,np_grid)
!     ir_grid(j)=ir_grid(j)+1
!  end do
!  do j=0,np_grid-1
!     id_grid(j)=sum(ir_grid(0:j))-ir_grid(j)
!  end do

  call InitParallel_RgridSol( node_partition, np_grid, pinfo_grid )

  id_grid(:) = pinfo_grid(7,:)
  ir_grid(:) = pinfo_grid(8,:)

  if ( myrank == 0 ) then
     do i=0,np_grid-1
        write(*,*) i,id_grid(i)+1,id_grid(i)+ir_grid(i),ir_grid(i)
     end do
     write(*,*) "sum(ir_grid),ML=",sum(ir_grid),ML
  end if

! ---

  n1 = id_grid(myrank_g)+1
  n2 = id_grid(myrank_g)+ir_grid(myrank_g)
  b1 = id_band(myrank_b)+1
  b2 = id_band(myrank_b)+ir_band(myrank_b)
  k1 = id_bzsm(myrank_k)+1
  k2 = id_bzsm(myrank_k)+ir_bzsm(myrank_k)
  s1 = id_spin(myrank_s)+1
  s2 = id_spin(myrank_s)+ir_spin(myrank_s)

  if ( flag_real8 ) then
     allocate( urtmp(n1:n2,b1:b2,k1:k2,s1:s2) ) ; urtmp=0.0d0
     allocate( rtmp(ML) ) ; rtmp=0.0d0
  else
     allocate( uztmp(n1:n2,b1:b2,k1:k2,s1:s2) ) ; uztmp=(0.0d0,0.0d0)
     allocate( ztmp(ML) ) ; ztmp=(0.0d0,0.0d0)
  end if

  allocate( rho(0:ML1-1,0:ML2-1,0:ML3-1) ) ; rho=0.0d0

! ---

  write(cmyrank,'(i5.5)') myrank
  file_wf_split = trim(file_name)//"."//trim(adjustl(cmyrank))

  open(u1,file=file_wf_split,status='old',form='unformatted')

  do s=id_spin(myrank_s)+1,id_spin(myrank_s)+ir_spin(myrank_s)
  do k=id_bzsm(myrank_k)+1,id_bzsm(myrank_k)+ir_bzsm(myrank_k)
  do n=id_band(myrank_b)+1,id_band(myrank_b)+ir_band(myrank_b)

     if ( flag_real8 ) then
        read(u1) urtmp(:,n,k,s)
     else
        read(u1) uztmp(:,n,k,s)
     end if

  end do ! n
  end do ! k
  end do ! s

  close(u1)


! ---

  do i=1,ndata

     n = idata(1,i)
     k = idata(2,i)
     s = idata(3,i)

     if ( b1 <= n .and. n <= b2 .and. &
          k1 <= k .and. k <= k2 .and. &
          s1 <= s .and. s <= s2 ) then

        if ( flag_real8 ) then
           call mpi_allgatherv( urtmp(n1,n,k,s),ir_grid(myrank_g) &
           ,mpi_real8,rtmp,ir_grid,id_grid,mpi_real8,comm_grid,ierr)
        else
           call mpi_allgatherv( uztmp(n1,n,k,s),ir_grid(myrank_g) &
          ,ztmp,mpi_complex16,ir_grid,id_grid,mpi_complex16,comm_grid,ierr)
        end if

        call gen_cube(n,k,s)

     end if

     call mpi_barrier( mpi_comm_world, ierr )

  end do ! i

! ---

  call end_mpi_parallel

! ---

CONTAINS


  SUBROUTINE gen_cube(n,k,s)

    implicit none

    integer,intent(IN) :: n,k,s
    integer :: i
    real(8) :: aa_del(3,3), r0(3)
    character(30) :: name, file_out
    character(1) :: cs
    character(3) :: ck
    character(5) :: cn

    aa_del(:,1)=aa(:,1)/ML1
    aa_del(:,2)=aa(:,2)/ML2
    aa_del(:,3)=aa(:,3)/ML3

    name  = 'SYS1'
    r0(:) = 0.0d0

    write(cn,'(i5.5)') n
    write(ck,'(i3.3)') k
    write(cs,'(i1.1)') s

    file_out = "wf_"//cn//"_"//ck//"_"//cs//".cub"

    select case( iflag_abs )
    case( 0 )
       if ( flag_real8 ) then
          do i=1,ML
             rho(LL(1,i),LL(2,i),LL(3,i)) = rtmp(i)
          end do
       else
          stop "complex wf can not be plot as it is"
       end if
    case( 1 )
       if ( flag_real8 ) then
          do i=1,ML
             rho(LL(1,i),LL(2,i),LL(3,i)) = abs( rtmp(i) )
          end do
       else
          do i=1,ML
             rho(LL(1,i),LL(2,i),LL(3,i)) = abs( ztmp(i) )
          end do
       end if
    case( 2 )
       if ( flag_real8 ) then
          do i=1,ML
             rho(LL(1,i),LL(2,i),LL(3,i)) = rtmp(i)**2
          end do
       else
          do i=1,ML
             rho(LL(1,i),LL(2,i),LL(3,i)) = abs( ztmp(i) )**2
          end do
       end if
    end select

    if ( myrank_g == 0 ) then

    open(iu,file=file_out)
    write(iu,100) name
    write(iu,100) name
    write(iu,110) MI,r0
    write(iu,110) ML1,aa_del(:,1)
    write(iu,110) ML2,aa_del(:,2)
    write(iu,110) ML3,aa_del(:,3)
    do i=1,MI
       write(iu,110) zatom(Kion(i)),real(zatom(Kion(i)),8),rsi(:,i)
    end do
    do i1=0,ML1-1
       do i2=0,ML2-1
          write(iu,*) rho(i1,i2,:)
       end do
    end do
    close(iu)

    end if

100 format(a4)
110 format(i5,4f12.6)

  END SUBROUTINE gen_cube


  SUBROUTINE InitParallel_RgridSol( np, np_grid, pinfo_grid )
    implicit none
    integer,intent(IN)  :: np(3), np_grid
    integer,intent(OUT) :: pinfo_grid(8,0:np_grid-1)
    integer :: i,j,n,i1,i2,i3
    integer,allocatable :: ntmp(:,:)
    n=maxval(np)
    allocate( ntmp(0:n-1,3) ) ; ntmp=0
    do j=1,3
       do i=0,Ngrid(j)-1
          n=mod(i,np(j))
          ntmp(n,j)=ntmp(n,j)+1
       end do
    end do
    n=-1
    i= 0
    do i3=0,np(3)-1
    do i2=0,np(2)-1
    do i1=0,np(1)-1
       n=n+1
       pinfo_grid(1,n)=sum( ntmp(0:i1,1) )-ntmp(i1,1)
       pinfo_grid(2,n)=ntmp(i1,1)
       pinfo_grid(3,n)=sum( ntmp(0:i2,2) )-ntmp(i2,2)
       pinfo_grid(4,n)=ntmp(i2,2)
       pinfo_grid(5,n)=sum( ntmp(0:i3,3) )-ntmp(i3,3)
       pinfo_grid(6,n)=ntmp(i3,3)
       pinfo_grid(7,n)=i
       pinfo_grid(8,n)=ntmp(i1,1)*ntmp(i2,2)*ntmp(i3,3)
       i=i+pinfo_grid(8,n)
    end do
    end do
    end do
    deallocate( ntmp )
  END SUBROUTINE InitParallel_RgridSol


END PROGRAM cubegen_wf3
