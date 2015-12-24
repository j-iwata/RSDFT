MODULE io2_module

  use parallel_module
  use rgrid_module, only: Igrid, dV, Hgrid
  use array_bound_module
  use wf_module
  use watch_module

  implicit none

  PRIVATE
  PUBLIC :: read_io2, read_data_io2

  integer :: nprocs_new, nmax
  integer,allocatable :: nmap(:),iomap(:,:,:)
  integer :: MB_0_IO, MB_1_IO

CONTAINS


  SUBROUTINE read_io2( rank )
    implicit none
    integer,intent(IN) :: rank
    integer,parameter  :: unit0=10
    integer :: i,j,ierr,idummy
    integer :: ML_tmp,ML1_tmp,ML2_tmp,ML3_tmp
    integer :: MB_tmp,MB1_tmp,MB2_tmp
    character(10) :: cbuf
!---
    if ( rank == 0 ) then
       open(unit0,file='input_wfiodir',status='old')
       read(unit0,*) nprocs_new
       allocate( nmap(nprocs_new) ) ; nmap=0
       do i=1,nprocs_new
          read(unit0,*) nmap(i)
          do j=1,nmap(i)
             read(unit0,*)
          end do
       end do
       nmax=maxval( nmap )
       allocate( iomap(0:6,nmax,nprocs_new) ) ; iomap=0
       rewind unit0
       read(unit0,*) nprocs_new
       do i=1,nprocs_new
          read(unit0,*) nmap(i)
          do j=1,nmap(i)
             read(unit0,*) iomap(0,j,i)
          end do
       end do
       do i=1,nprocs_new
          read(unit0,*)
          do j=1,nmap(i)
             read(unit0,*) cbuf,idummy,iomap(1:6,j,i)
          end do
       end do
       close(unit0)
    end if
!---
    call mpi_bcast(nmax,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(nprocs_new,1,mpi_integer,0,mpi_comm_world,ierr)
    if ( rank /= 0 ) then
       allocate( nmap(nprocs_new) ) ; nmap=0
       allocate( iomap(0:6,nmax,nprocs_new) ) ; iomap=0
    end if
    call mpi_bcast(nmap,nprocs_new,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(iomap,size(iomap),mpi_integer,0,mpi_comm_world,ierr)
!---
    if ( myrank == 0 ) then
       open(3,file="wf.dat1",status="old",form="unformatted")
       read(3) ML_tmp,ML1_tmp,ML2_tmp,ML3_tmp
       read(3) MB_tmp,MB1_tmp,MB2_tmp
       close(3)
       write(*,*) "MB_0_IO,MB_1_IO=",1,MB_tmp
    end if
    call mpi_bcast(MB_tmp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    MB_0_IO=1
    MB_1_IO=MB_tmp
!---
  END SUBROUTINE read_io2


  SUBROUTINE read_data_io2( disp_switch )
    implicit none
    logical,intent(IN) :: disp_switch
    integer,parameter :: u1=10
    integer :: i,n,k,s,i1,i2,i3,i0,a1a,b1a,a2a,b2a,a3a,b3a
    integer :: a1c,b1c,a2c,b2c,a3c,b3c,a1b,b1b,a2b,b2b,a3b,b3b
    integer,allocatable :: LLL(:,:,:)
    character(5) :: crank_old
    character(32) :: filename
    real(8) :: ct0,ct1,et0,et1
#ifdef _DRSDFT_
    real(8),allocatable :: w(:,:,:)
    real(8),parameter :: z0=0.0d0
#else
    complex(8),allocatable :: w(:,:,:)
    complex(8),parameter :: z0=(0.0d0,0.0d0)
#endif
    if ( disp_switch_parallel ) then
       write(*,'(a40," read_data_io2")') repeat("-",60)
    end if
    call watch(ct0,et0)
    call read_io2(myrank)
    a1b=Igrid(1,1)
    b1b=Igrid(2,1)
    a2b=Igrid(1,2)
    b2b=Igrid(2,2)
    a3b=Igrid(1,3)
    b3b=Igrid(2,3)
    allocate( LLL(a1b:b1b,a2b:b2b,a3b:b3b) ) ; LLL=0
    i=ML_0_WF-1
    do i3=a3b,b3b
    do i2=a2b,b2b
    do i1=a1b,b1b
       i=i+1
       LLL(i1,i2,i3)=i
    end do
    end do
    end do
    do i=1,nmap(myrank+1)
       write(crank_old,'(i5.5)') iomap(0,i,myrank+1)
       filename="wf.dat1."//crank_old
       if ( disp_switch_parallel ) write(*,*) filename
       open(u1,file=filename,form='unformatted')
       a1a=iomap(1,i,myrank+1)
       b1a=iomap(2,i,myrank+1)
       a2a=iomap(3,i,myrank+1)
       b2a=iomap(4,i,myrank+1)
       a3a=iomap(5,i,myrank+1)
       b3a=iomap(6,i,myrank+1)
       allocate( w(a1a:b1a,a2a:b2a,a3a:b3a) ) ; w=z0
       a1c=max(a1b,a1a)
       b1c=min(b1b,b1a)
       a2c=max(a2b,a2a)
       b2c=min(b2b,b2a)
       a3c=max(a3b,a3a)
       b3c=min(b3b,b3a)
       do s=MS_0_WF,MS_1_WF
       do k=MK_0_WF,MK_1_WF
       do n=MB_0_IO,MB_1_IO
          w=z0
          read(u1) w
          do i3=a3c,b3c
          do i2=a2c,b2c
          do i1=a1c,b1c
             i0=LLL(i1,i2,i3)
             unk(i0,n,k,s)=w(i1,i2,i3)
          end do
          end do
          end do
       end do ! n
       end do ! k
       end do ! s
       deallocate( w )
       close(u1)
    end do ! i
    deallocate( LLL )
    call watch(ct1,et1)
    if ( disp_switch_parallel ) then
       write(*,*) "time(read_data_io2)=",ct1-ct0,et1-et0
       write(*,'(a40," read_data_io2(end)")') repeat("-",60)
    end if
  END SUBROUTINE read_data_io2


END MODULE io2_module
