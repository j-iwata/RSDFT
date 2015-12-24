MODULE gram_schmidt_t_module

  use rgrid_module, only: dV,zdV
  use wf_module, only: unk,hunk,iflag_hunk
  use array_bound_module, only: ML_0,ML_1,MB,MB_0
  use parallel_module

  implicit none

  PRIVATE
  PUBLIC :: gram_schmidt_t, read_gram_schmidt_t, read_oldformat_gram_schmidt_t

#ifdef _DRSDFT_
  integer,parameter :: TYPE_MAIN = MPI_REAL8
  character(1),parameter :: TRANSA='T'
  character(1),parameter :: TRANSB='N'
  real(8),allocatable :: utmp2(:,:),utmp(:)
  real(8),parameter :: zero=0.d0,one=1.d0
#else
  integer,parameter :: TYPE_MAIN = MPI_COMPLEX16
  character(1),parameter :: TRANSA='C'
  character(1),parameter :: TRANSB='N'
  complex(8),allocatable :: utmp2(:,:),utmp(:)
  complex(8),parameter :: zero=(0.d0,0.d0),one=(1.d0,0.d0)
#endif

  integer :: NBLK,NBLK1

CONTAINS

!---------------------------------------------------------------------------------------
  SUBROUTINE read_gram_schmidt_t(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i
    character(5) :: cbuf,ckey
    NBLK=0
    NBLK1=0
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:5) == "NBLK " ) then
             backspace(unit)
             read(unit,*) cbuf,NBLK
          else if ( ckey(1:5) == "NBLK1" ) then
             backspace(unit)
             read(unit,*) cbuf,NBLK1
          end if
       end do
999    continue
       write(*,*) "NBLK =",NBLK
       write(*,*) "NBLK1=",NBLK1
    end if
    call send_gram_schmidt_t(0)
  END SUBROUTINE read_gram_schmidt_t

!-----------------------------------------------------------------------

  SUBROUTINE read_oldformat_gram_schmidt_t(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    if ( rank == 0 ) then
       read(unit,*) NBLK,NBLK1
       write(*,*) "NBLK,NBLK1=",NBLK,NBLK1
    end if
    call send_gram_schmidt_t(0)
  END SUBROUTINE read_oldformat_gram_schmidt_t

!-----------------------------------------------------------------------

  SUBROUTINE send_gram_schmidt_t(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(NBLK ,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(NBLK1,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_gram_schmidt_t

!-----------------------------------------------------------------------

  SUBROUTINE gram_schmidt_t(n0,n1,k,s)
    implicit none
    integer,intent(IN) :: n0,n1,k,s
    integer :: irank_b,ns,ne,ms,me,n,ML0,ierr
    integer :: nbss,k1,ib,NBAND_BLK,ncycle,mrnk
    integer,allocatable :: ir(:),id(:)

    ML0   = ML_1-ML_0+1
    mrnk  = id_class(myrank,4)

    if ( NBLK  == 0 ) NBLK =MB
    if ( NBLK1 == 0 ) NBLK1=4

    allocate( ir(0:np_band-1),id(0:np_band-1) )
    ir(0:np_band-1)=ir_band(0:np_band-1)*ML0
    id(0:np_band-1)=id_band(0:np_band-1)*ML0

    NBAND_BLK = NBLK
    ncycle    = (MB-1)/NBAND_BLK+1

    call mpi_allgatherv(unk(ML_0,MB_0,k,s),ir(mrnk),TYPE_MAIN &
          ,unk(ML_0,1,k,s),ir,id,TYPE_MAIN,comm_band,ierr)
    if ( iflag_hunk >= 1 ) then
       call mpi_allgatherv(hunk(ML_0,MB_0,k,s),ir(mrnk),TYPE_MAIN &
            ,hunk(ML_0,1,k,s),ir,id,TYPE_MAIN,comm_band,ierr)
    end if

    do k1=1,ncycle

       irank_b=mod(k1-1,np_band)

       ns=NBAND_BLK*(k1-1)+1
       ne=min(ns+NBAND_BLK-1,MB)

       if ( id_class(myrank,4) == irank_b ) then

          call gram_schmidt_sub(ns,ne,ns,ne,NBLK,k,s)

       end if

       n=ML0*(ne-ns+1)
       call mpi_bcast(unk(ML_0,ns,k,s),n,TYPE_MAIN,irank_b,comm_band,ierr)
       if ( iflag_hunk >= 1 ) then
          call mpi_bcast(hunk(ML_0,ns,k,s),n,TYPE_MAIN,irank_b,comm_band,ierr)
       end if

       if ( ns <= MB-NBAND_BLK ) then

          do ib=1,(ncycle-1)/np_band+1

             nbss=(ib-1)*np_band+myrank_b+1

             if ( nbss<=ncycle .and. nbss>= k1+1 ) then

                ms=NBAND_BLK*(nbss-1)+1
                me=min(ms+NBAND_BLK-1,MB)

                if ( ms<=me ) then
                   call gram_schmidt_sub(ms,me,ns,ne,NBLK,k,s)
                end if

             end if

          enddo

       end if

    end do ! k1

    deallocate( id,ir )

    return

  END SUBROUTINE gram_schmidt_t

!-----------------------------------------------------------------------

  RECURSIVE SUBROUTINE gram_schmidt_sub(mm1,mm2,nn1,nn2,MBLK,k,s)
    implicit none
    integer,intent(IN) :: mm1,mm2,nn1,nn2,MBLK,k,s
    integer :: n,ns,ne,m,ms,me,m1,m2,mm,nn,MBLKH,ierr,ML0,i
    real(8) :: c

    ML0 = ML_1-ML_0+1

    do ms=mm1,mm2,MBLK

       me=min(ms+MBLK-1,mm2)
       mm=me-ms+1

       do ns=nn1,nn2,MBLK

          ne=min(ns+MBLK-1,nn2)
          ne=min(ne,me-1)
          nn=ne-ns+1

          if ( nn<=0 ) cycle

          if ( ms>=ne+1 ) then

             allocate( utmp2(ns:ne,ms:me) )

#ifdef _DRSDFT_
             call dgemm(TRANSA,TRANSB,nn,mm,ML0, -dV,unk(ML_0,ns,k,s) &
                  ,ML0,unk(ML_0,ms,k,s),ML0,zero,utmp2,nn)
#else
             call zgemm(TRANSA,TRANSB,nn,mm,ML0,-zdV,unk(ML_0,ns,k,s) &
                  ,ML0,unk(ML_0,ms,k,s),ML0,zero,utmp2,nn)
#endif

             call mpi_allreduce(MPI_IN_PLACE,utmp2,nn*mm,TYPE_MAIN,mpi_sum &
                  ,comm_grid,ierr)

#ifdef _DRSDFT_
             call dgemm(TRANSB,TRANSB,ML0,mm,nn,one,unk(ML_0,ns,k,s) &
                  ,ML0,utmp2,nn,one,unk(ML_0,ms,k,s),ML0)
             if ( iflag_hunk >= 1 ) then
                call dgemm(TRANSB,TRANSB,ML0,mm,nn,one,hunk(ML_0,ns,k,s) &
                     ,ML0,utmp2,nn,one,hunk(ML_0,ms,k,s),ML0)
             end if
#else
             call zgemm(TRANSB,TRANSB,ML0,mm,nn,one,unk(ML_0,ns,k,s) &
                  ,ML0,utmp2,nn,one,unk(ML_0,ms,k,s),ML0)
             if ( iflag_hunk >= 1 ) then
                call zgemm(TRANSB,TRANSB,ML0,mm,nn,one,hunk(ML_0,ns,k,s) &
                     ,ML0,utmp2,nn,one,hunk(ML_0,ms,k,s),ML0)
             end if
#endif

             deallocate( utmp2 )

             if ( ms==ne+1 ) then

                c=0.d0
                do i=ML_0,ML_1
                   c=c+abs(unk(i,ms,k,s))**2
                end do

                call mpi_allreduce(MPI_IN_PLACE,c,1,mpi_real8,mpi_sum,comm_grid,ierr)

                c=1.d0/sqrt(c*dV)

                do i=ML_0,ML_1
                   unk(i,ms,k,s)=c*unk(i,ms,k,s)
                end do
                if ( iflag_hunk >= 1 ) then
                   do i=ML_0,ML_1
                      hunk(i,ms,k,s)=c*hunk(i,ms,k,s)
                   end do
                endif

             end if

          else if ( mm<=NBLK1 ) then

             allocate( utmp(NBLK1) )

             do m=ms,me

                n=min(m-1,ne)

                if ( n-ns+1>0 ) then

#ifdef _DRSDFT_
                   call dgemv(TRANSA,ML0,n-ns+1, -dV,unk(ML_0,ns,k,s) &
                        ,ML0,unk(ML_0,m,k,s),1,zero,utmp,1)
#else
                   call zgemv(TRANSA,ML0,n-ns+1,-zdV,unk(ML_0,ns,k,s) &
                        ,ML0,unk(ML_0,m,k,s),1,zero,utmp,1)
#endif

                   call mpi_allreduce(MPI_IN_PLACE,utmp,n-ns+1,TYPE_MAIN &
                        ,mpi_sum,comm_grid,ierr)

#ifdef _DRSDFT_
                   call dgemv(TRANSB,ML0,n-ns+1,one,unk(ML_0,ns,k,s) &
                        ,ML0,utmp,1,one,unk(ML_0,m,k,s),1)
                   if ( iflag_hunk >= 1 ) then
                      call dgemv(TRANSB,ML0,n-ns+1,one,hunk(ML_0,ns,k,s) &
                           ,ML0,utmp,1,one,hunk(ML_0,m,k,s),1)
                   end if
#else
                   call zgemv(TRANSB,ML0,n-ns+1,one,unk(ML_0,ns,k,s) &
                        ,ML0,utmp,1,one,unk(ML_0,m,k,s),1)
                   if ( iflag_hunk >= 1 ) then
                      call zgemv(TRANSB,ML0,n-ns+1,one,hunk(ML_0,ns,k,s) &
                           ,ML0,utmp,1,one,hunk(ML_0,m,k,s),1)
                   end if
#endif

                end if

                if ( m==1 .or. (n==m-1 .and. m/=ns) ) then

                   c=0.d0

                   do i=ML_0,ML_1
                      c=c+abs(unk(i,m,k,s))**2
                   end do

                   call mpi_allreduce(MPI_IN_PLACE,c,1,mpi_real8,mpi_sum,comm_grid,ierr)

                   c=1.d0/sqrt(c*dV)

                   do i=ML_0,ML_1
                      unk(i,m,k,s)=c*unk(i,m,k,s)
                   end do
                   if ( iflag_hunk >= 1 ) then
                      do i=ML_0,ML_1
                         hunk(i,m,k,s)=c*hunk(i,m,k,s)
                      end do
                   end if

                end if

             end do ! m

             deallocate( utmp )

          else

             MBLKH=max(MBLK/2,NBLK1)
             call gram_schmidt_sub(ms,me,ns,ne,MBLKH,k,s)

          end if

       end do ! ns

    end do ! ms

    return

  END SUBROUTINE gram_schmidt_sub

END MODULE gram_schmidt_t_module
