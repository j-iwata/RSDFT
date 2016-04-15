PROGRAM dos2gp

  !use libtetrabz, only: libtetrabz_dos
  use libtetrabz_mpi, only: libtetrabz_mpi_dos

  implicit none

  integer,parameter :: u1=1, u2=2, u5=5, u99=99, u10=10
  integer,parameter :: max_loop=10000000
  integer :: ltetra
  integer :: k,i,k1,k2,k3,i1,i2,i3,n,ne,s
  integer :: nk, n_irwedge, n_whole, indx_range(2,3)
  integer :: nge(3),ngw(3)
  integer :: nband,nspin,nkpnt
  integer :: myrank
  integer,allocatable :: kbz(:,:),ktmp(:,:,:),kgrd(:,:,:)
  integer,allocatable :: itmp(:,:)
  real(8),allocatable :: eval(:,:,:),occp(:,:,:)
  real(8),allocatable :: wbz(:),kpnt(:,:),weight(:)
  real(8),allocatable :: eig(:,:,:,:),wght(:,:,:,:,:),e0(:)
  real(8),parameter :: HT=27.2116d0
  real(8) :: emin,emax,de,e1,e2
  real(8) :: bb(3,3), ktry(3), err, etime(0:1)
  character(80) :: cbuf
  character(16) :: file_bz, file_eigv
  character(8) :: cbuf3(3)
  include 'mpif.h'

  call MPI_INIT(i)
  call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, i )

  etime(0)=mpi_wtime()

  file_bz   = "bz_info"
  file_eigv = "eigenvalues"

!------------------------------

  open(u1,file=file_bz,status='old')
  read(u1,*)
  read(u1,*) bb(1:3,1)
  read(u1,*) bb(1:3,2)
  read(u1,*) bb(1:3,3)
  read(u1,'(a)') cbuf
  call get_num_from_chr( cbuf, nk )
  read(u1,'(a)') cbuf
  call get_num_from_chr( cbuf, n_irwedge )
  read(u1,'(a)') cbuf
  call get_num_from_chr( cbuf, n_whole )

  read(u1,*)
  read(u1,*)

  allocate( kbz(3,n_irwedge) ) ; kbz=0
  allocate( wbz(n_irwedge)   ) ; wbz=0.0d0

  do i=1,n_irwedge
     read(u1,*) k, kbz(1:3,i), wbz(i)
  end do

  read(u1,'(a)') cbuf
  call get_nums_from_chr( cbuf, indx_range(:,1) )
  read(u1,'(a)') cbuf
  call get_nums_from_chr( cbuf, indx_range(:,2) )
  read(u1,'(a)') cbuf
  call get_nums_from_chr( cbuf, indx_range(:,3) )

  allocate( ktmp(indx_range(1,1):indx_range(2,1) &
                ,indx_range(1,2):indx_range(2,2) &
                ,indx_range(1,3):indx_range(2,3)) ) ; ktmp=0

  read(u1,*)

  do i=1,n_whole
     read(u1,*) k,k1,k2,k3,ktmp(k1,k2,k3)
  end do

  close(u1)

! ---

  k=max( size(ktmp,1), size(ktmp,2), size(ktmp,3) )
  allocate( itmp(k,3) ) ; itmp=0 

  loop_3 : do k3=indx_range(1,3),indx_range(2,3)
  do k2=indx_range(1,2),indx_range(2,2)
     if ( .not.( all(ktmp(:,k2,k3)/=0) ) ) exit loop_3
  end do
  end do loop_3
  nge(1) = count( ktmp(:,k2,k3) /= 0 )
  i=0
  do k1=indx_range(1,1),indx_range(2,1)
     if ( ktmp(k1,k2,k3) /= 0 ) then
        i=i+1
        itmp(i,1)=k1
     end if
  end do

  loop_1 : do k1=indx_range(1,1),indx_range(2,1)
  do k3=indx_range(1,3),indx_range(2,3)
     if ( .not.( all(ktmp(k1,:,k3)/=0) ) ) exit loop_1
  end do
  end do loop_1
  nge(2) = count( ktmp(k1,:,k3) /= 0 )
  i=0
  do k2=indx_range(1,2),indx_range(2,2)
     if ( ktmp(k1,k2,k3) /= 0 ) then
        i=i+1
        itmp(i,2)=k2
     end if
  end do

  loop_2 : do k2=indx_range(1,2),indx_range(2,2)
  do k1=indx_range(1,1),indx_range(2,1)
     if ( .not.( all(ktmp(k1,k2,:)/=0) ) ) exit loop_2
  end do
  end do loop_2
  nge(3) = count( ktmp(k1,k2,:) /= 0 )
  i=0
  do k3=indx_range(1,3),indx_range(2,3)
     if ( ktmp(k1,k2,k3) /= 0 ) then
        i=i+1
        itmp(i,3)=k3
     end if
  end do

! ---

  allocate( kgrd(nge(1),nge(2),nge(3)) ) ; kgrd=0

  do i3=1,nge(3)
  do i2=1,nge(2)
  do i1=1,nge(1)
     kgrd(i1,i2,i3) = ktmp(itmp(i1,1),itmp(i2,2),itmp(i3,3))
  end do
  end do
  end do

! ---

  open(u2,file=file_eigv)

  nkpnt=0

  do i=1,max_loop
     read(u2,*,END=10) cbuf
     if ( cbuf == "Nband" ) then
        backspace(u2)
        read(u2,*) cbuf3, nband, nspin, k
        nkpnt=nkpnt+1
     end if
  end do
10 continue

  allocate( eval(nband,nkpnt,nspin) ) ; eval=0.0d0
  allocate( kpnt(3,nkpnt) ) ; kpnt=0.0d0
  allocate( weight(nkpnt) ) ; weight=0.0d0

  rewind u2

  do k=1,nkpnt
     read(u2,*)
     read(u2,*) kpnt(1:3,k),weight(k)
     do n=1,nband
        read(u2,*) i,( eval(n,k,s), s=1,nspin )
     end do
  end do

  close(u2)

! ---

  allocate( eig(nband,nge(1),nge(2),nge(3)) ) ; eig=0.0d0

  do i3=1,nge(3)
  do i2=1,nge(2)
  do i1=1,nge(1)

     i = kgrd(i1,i2,i3)

     ktry(1:3) = dble(kbz(1:3,i))/dble(nk)

     do k=1,nkpnt
        err = sum( (ktry-kpnt(:,k))**2 )
        if ( err < 1.d-10 ) then
           do n=1,nband
              eig(n,i1,i2,i3) = eval(n,k,1)
           end do
           exit
        end if
     end do ! k

  end do ! i1
  end do ! i2
  end do ! i3

!------------------------------

  emin = minval( eval )
  emax = maxval( eval )

  if ( myrank == 0 ) then
     write(*,*) "emin     (HT,eV)=",emin,emin*HT
     write(*,*) "emax     (HT,eV)=",emax,emax*HT
     write(*,*) "emax-emin(HT,eV)=",emax-emin,(emax-emin)*HT
  end if

!------------------------------

  if ( myrank == 0 ) then
     write(*,*) "# of energy grid: ne="
     read(u5,*) ne
     write(*,*) ne
  end if
  call MPI_BCAST( ne, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, i )

  if ( myrank == 0 ) then
     write(*,*) "energy range: e1,e2= ( 0,0: default values are set )"
     read(u5,*) e1,e2
  end if
  call MPI_BCAST( e1, 1, MPI_REAL8, 0, MPI_COMM_WORLD, i )
  call MPI_BCAST( e2, 1, MPI_REAL8, 0, MPI_COMM_WORLD, i )
  if ( e1 == 0.0d0 ) e1=emin
  if ( e2 == 0.0d0 ) e2=emax
  if ( myrank == 0 ) write(*,*) e1,e2

  if ( myrank == 0 ) then
     write(*,*) "ltetra [1:linear tetrahedron, 2:opt-linear tetrahedron]"
     read(u5,*) ltetra
     write(*,*) ltetra
  end if
  call MPI_BCAST( ltetra, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, i )

! ---

  ngw(1:3) = nge(1:3)

  allocate( wght(ne,nband,ngw(1),ngw(2),ngw(3)) ) ; wght=0.0d0
  allocate( e0(ne) ) ; e0=0.0d0


  de=(emax-emin)/ne
  do i=1,ne
     e0(i) = emin + (i-1)*de
  end do

  !call libtetrabz_dos( ltetra, bb, nband, nge, eig, ngw, wght, ne, e0 )
  call libtetrabz_mpi_dos &
       ( ltetra, MPI_COMM_WORLD, bb, nband, nge, eig, ngw, wght, ne, e0 )

  if ( myrank == 0 ) then
     rewind 10
     do i=1,ne
        write(10,*) e0(i), sum( wght(i,:,:,:,:) )
     end do
  end if

  etime(1)=mpi_wtime()
  if ( myrank == 0 ) write(*,*) "TIME=",etime(1)-etime(0)

  call MPI_FINALIZE( i )

CONTAINS

  SUBROUTINE get_num_from_chr( cbuf, n )
    implicit none
    character(*),intent(IN) :: cbuf
    integer,intent(OUT) :: n
    integer :: i
    do i=1,len_trim(cbuf)
       if ( cbuf(i:i) == ":" ) exit
    end do
    read(cbuf(i+1:),*) n
  END SUBROUTINE get_num_from_chr

  SUBROUTINE get_nums_from_chr( cbuf, n )
    implicit none
    character(*),intent(IN) :: cbuf
    integer,intent(OUT) :: n(:)
    integer :: i
    do i=1,len_trim(cbuf)
       if ( cbuf(i:i) == ":" ) exit
    end do
    read(cbuf(i+1:),*) n(:)
  END SUBROUTINE get_nums_from_chr


END PROGRAM dos2gp
