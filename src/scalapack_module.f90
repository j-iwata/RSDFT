MODULE scalapack_module

  use parallel_module
  use io_tools_module

  implicit none

  PRIVATE
  PUBLIC :: prep_scalapack, init_scalapack, UPLO &
           ,NPROW,NPCOL,MBSIZE,NBSIZE,LLD_R,LLD_C,DESCA,DESCB,DESCZ &
           ,NP0,NQ0,NPX,NQX,usermap
  PUBLIC :: read_scalapack

  integer :: NPROW=0
  integer :: NPCOL=0
  integer :: MBSIZE=0
  integer :: NBSIZE=0
  integer :: LLD_R,LLD_C
  integer :: NP0,NQ0,NPX,NQX
  integer :: DESCA(9),DESCB(9),DESCZ(9)
  integer,allocatable :: usermap(:,:,:)
  character(1) :: UPLO

  logical :: iblacs = .false.
  logical :: flag_read = .true.

CONTAINS


  SUBROUTINE init_scalapack( MB )
    implicit none
    integer,intent(INOUT) :: MB
    integer :: NPCOL0,i,j,n,loop

    call write_border( 0, " init_scalapack(start)" )

    MBSIZE = 0
    NBSIZE = 0
    iblacs = .false.

    if ( NPROW<1 .or. NPCOL<1 ) then
       NPCOL0 = node_partition(1)*node_partition(2) &
               *node_partition(3)*node_partition(4)
       NPCOL  = NPCOL0
       NPROW  = 1
       do i=2,node_partition(1)*node_partition(2)*node_partition(3)
          j=NPCOL0/i
          if ( i*j==NPCOL0 .and. i<=j .and. j-i<NPCOL-NPROW ) then
             NPCOL=j
             NPROW=i
          end if
       end do
    else
       n=node_partition(1)*node_partition(2)*node_partition(3)*node_partition(4)
       if ( NPROW*NPCOL > n ) then
          write(*,*) "NPROW,NPCOL,np_band,np_grid=",NPROW,NPCOL &
               ,node_partition(4) &
               ,node_partition(1)*node_partition(2)*node_partition(3),myrank
          stop
       end if
    end if

    if ( MBSIZE<1 .or. NBSIZE<1 ) then
       i=(MB+NPROW-1)/NPROW
       j=(MB+NPCOL-1)/NPCOL
       MBSIZE=min(i,j)
       NBSIZE=MBSIZE
    else
       if ( MBSIZE/=NBSIZE ) then
          write(*,*) "MBSIZE,NBSIZE=",MBSIZE,NBSIZE,myrank
          write(*,*) "MBSIZE /= NBSIZE may not work well"
          stop "stop@init_scalapack"
       end if
    end if

    if ( disp_switch_parallel ) then
       write(*,*) "NPROW,NPCOL,MBSIZE,NBSIZE"
       write(*,*) NPROW,NPCOL,MBSIZE,NBSIZE
    end if

    if ( NBSIZE*NPCOL/=MB ) then
       write(*,*) "NBSIZE*NPCOL/=MB!"
       n=max( NBSIZE*NPCOL, MB )
       n=min( n, (NBSIZE+1)*NPCOL )
       write(*,*) "recommended value for MB =",n
       write(*,*) "replace MB"
       MB=n
       MBSIZE=0
       NBSIZE=0
    end if

    call write_border( 0, " init_scalapack(end)" )

  END SUBROUTINE init_scalapack


  SUBROUTINE read_scalapack
    implicit none
    integer :: itmp(2)
    itmp(:)=-1
    call IOTools_readIntegerKeywords( "SCL", itmp )
    if ( itmp(1) > -1 ) NPROW=itmp(1)
    if ( itmp(2) > -1 ) NPCOL=itmp(2)
  END SUBROUTINE read_scalapack


  SUBROUTINE prep_scalapack( MB )
    implicit none
    integer,intent(INOUT) :: MB
    integer :: ierr,NPCOL0,i,j,n,is,ik,m,ib,l,i1,i2,i3,i7
    integer :: MXLLD,MYROW,MYCOL,mm,mchk
    integer,save :: icount_visit=0, ICTXT=0, ICTXT0=0
    integer :: NUMROC

#ifndef _LAPACK_

    if ( iblacs ) return

    call write_border( 0, " prep_scalapack(start)" )

    if ( icount_visit>0 ) then
       call blacs_gridexit(ICTXT)
    end if
    icount_visit=icount_visit+1

    iblacs = .true.

! --- NPROW,NPCOL,MBSIZE,NBSIZE ---

    if ( NPROW<1 .or. NPCOL<1 ) then
       NPCOL0 = np_band*np_grid
       NPCOL  = np_band*np_grid
       NPROW  = 1
       do i=2,np_grid
          j=NPCOL0/i
          if ( i*j==NPCOL0 .and. i<=j .and. j-i<NPCOL-NPROW ) then
             NPCOL=j
             NPROW=i
          end if
       end do
    else
       if ( NPROW*NPCOL > np_band*np_grid ) then
          write(*,*) "NPROW,NPCOL,np_band,np_grid=",NPROW,NPCOL,np_band,np_grid,myrank
          stop
       end if
    end if

    if ( MBSIZE<1 .or. NBSIZE<1 ) then
       i=(MB+NPROW-1)/NPROW
       j=(MB+NPCOL-1)/NPCOL
       MBSIZE=min(i,j)
       NBSIZE=MBSIZE
    else
       if ( MBSIZE/=NBSIZE ) then
          write(*,*) "MBSIZE,NBSIZE=",MBSIZE,NBSIZE,myrank
          stop
       end if
    end if

    if ( disp_switch_parallel ) then
       write(*,*) "NPROW,NPCOL,MBSIZE,NBSIZE"
       write(*,*) NPROW,NPCOL,MBSIZE,NBSIZE
    end if

    if ( NBSIZE*NPCOL/=MB ) then
       n=max( NBSIZE*NPCOL, MB )
       n=min( n, (NBSIZE+1)*NPCOL )
       if ( disp_switch_parallel ) then
          write(*,*) "NBSIZE*NPCOL/=MB!"
          write(*,*) "recommended value for MB =",n
          write(*,*) "MB is replaced"
       end  if
       MB=n
    end if

! --- preparation for ScaLAPACK ---

    if ( .not.allocated(usermap) ) then

       allocate( usermap(0:NPROW-1,0:NPCOL-1,2) )
       usermap(:,:,:)=MPI_PROC_NULL

       n=-1
       do i7=0,node_partition(7)-1
       do is=0,node_partition(6)-1
       do ik=0,node_partition(5)-1
          m=-1 ; mchk=-NPROW*NPCOL
          do ib=0,node_partition(4)-1
             l=-1
             do i3=0,node_partition(3)-1
             do i2=0,node_partition(2)-1
             do i1=0,node_partition(1)-1
                n=n+1
                m=m+1 ; if ( mod(m,NPROW*NPCOL) == 0 ) mchk=mchk+NPROW*NPCOL
                l=l+1
                i=mod(m+NPROW,NPROW)
                j=m/NPROW
                mm=myrank_g+nprocs_g*myrank_b+1
                if ( id_class(myrank,5)==ik .and. &
                     id_class(myrank,6)==is .and. &
                     id_class(myrank,7)==i7 .and. mm > mchk ) then
                   usermap(i,j,1)=n
                   usermap(i,j,2)=l
                end if
             end do ! i1
             end do ! i2
             end do ! i3
          end do ! ib
       end do ! ik
       end do ! is
       end do ! i7

    end if

    call blacs_get(0,0,ICTXT0)
    ICTXT = ICTXT0

    call blacs_gridmap(ICTXT,usermap(0,0,1),NPROW,NPROW,NPCOL)
    call blacs_gridinfo(ICTXT,NPROW,NPCOL,MYROW,MYCOL)

    LLD_R = NUMROC(MB,MBSIZE,MYROW,0,NPROW)
    LLD_C = NUMROC(MB,NBSIZE,MYCOL,0,NPCOL)
    MXLLD = LLD_R

    call descinit(DESCA,MB,MB,MBSIZE,NBSIZE,0,0,ICTXT,MXLLD,ierr)
    call descinit(DESCB,MB,MB,MBSIZE,NBSIZE,0,0,ICTXT,MXLLD,ierr)
    call descinit(DESCZ,MB,MB,MBSIZE,NBSIZE,0,0,ICTXT,MXLLD,ierr)

    NP0=NUMROC(MB,MBSIZE,0,0,NPROW)
    NQ0=NUMROC(MB,NBSIZE,0,0,NPCOL)

    NPX=NUMROC(MB,MBSIZE,MYROW,0,NPROW)
    NQX=NUMROC(MB,NBSIZE,MYCOL,0,NPCOL)

    if ( disp_switch_parallel ) then
       write(*,*) "ICTXTX,ICTXT0    =",ICTXT,ICTXT0
       write(*,*) "NPROW,NPCOL      =",NPROW,NPCOL
       write(*,*) "MBSIZE,NBSIZE    =",MBSIZE,NBSIZE
       write(*,*) "MYROW,MYCOL      =",MYROW,MYCOL
       write(*,*) "LLD_R,LLD_C,MXLLD=",LLD_R,LLD_C,MXLLD
       write(*,*) "NP0,NQ0          =",NP0,NQ0
       write(*,*) "NPX,NQX          =",NPX,NQX
       write(*,*) "iblacs           =",iblacs
    end if

    call write_border( 0, " prep_scalapack(end)" )

#endif

  END SUBROUTINE prep_scalapack


END MODULE scalapack_module
