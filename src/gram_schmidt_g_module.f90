MODULE gram_schmidt_g_module

  use rgrid_module, only: dV,zdV
  use wf_module, only: unk,Sunk,ML_0_WF,ML_1_WF
  use array_bound_module, only: ML_0,ML_1,MB,MB_0
  use parallel_module, only: comm_grid,comm_band,ir_band,id_band &
                            ,id_class,myrank,np_band,myrank_b
  use real_complex_module, only: zero,one,TYPE_MAIN,TRANSA,TRANSB
  use inner_product_module
  use var_sys_parameter, only: pp_kind

  implicit none

  include 'mpif.h'

  PRIVATE
  PUBLIC :: init_gram_schmidt_g
  PUBLIC :: gram_schmidt_g

  integer :: NBLK, NBLK1

CONTAINS


  SUBROUTINE init_gram_schmidt_g( iswitch_algorithm )
    implicit none
    integer,intent(INOUT) :: iswitch_algorithm
    if ( pp_kind == "USPP" ) iswitch_algorithm=101
  END SUBROUTINE init_gram_schmidt_g


  ! Gram-Schmidt orthogonalization
  ! ( Takahashi, Block cyclic )
  SUBROUTINE gram_schmidt_g(ni,nb,k,s)
    implicit none
    integer,intent(IN) :: ni,nb,k,s
    integer :: irank_b,ns,ne,ms,me,n,ML0,ierr
    integer :: nbss,k1,ib,NBAND_BLK,ncycle,mrnk
    integer,allocatable :: ir(:),id(:)
    integer :: n0,n1

    ML0 = ML_1 - ML_0 + 1
    mrnk = id_class(myrank,4)
    n0=ML_0_WF
    n1=ML_1_WF

    if ( NBLK  == 0 ) NBLK =MB
    if ( NBLK1 == 0 ) NBLK1=4

    allocate( ir(0:np_band-1),id(0:np_band-1) ) ; ir=0 ; id=0

    ir(0:np_band-1)=ir_band(0:np_band-1)*ML0
    id(0:np_band-1)=id_band(0:np_band-1)*ML0

    NBAND_BLK=NBLK
    ncycle=(MB-1)/NBAND_BLK+1

    call MPI_ALLGATHERV( unk(ML_0,MB_0,k,s),ir(mrnk),TYPE_MAIN,unk(ML_0,1,k,s),ir,id,TYPE_MAIN,COMM_BAND,ierr )
    do k1=1,ncycle
      irank_b = mod(k1-1,np_band)

      ns = NBAND_BLK*(k1-1) + 1
      ne = min(ns+NBAND_BLK-1,MB)

      if ( id_class(myrank,4)==irank_b ) then
        allocate( Sunk(n0:n1,ns:ne) ) ; Sunk=zero
        call get_Sunk_Mat( n0,n1,ns,ne,k,s )
        call GramSchmidtGSub(ns,ne,ns,ne,NBLK,k,s,NBLK,NBLK1)
        deallocate( Sunk )
      end if
      n=ML0*(ne-ns+1)
      call MPI_BCAST( unk(ML_0,ns,k,s),n,TYPE_MAIN,irank_b,COMM_BAND,ierr )

      if ( ns <= MB-NBAND_BLK ) then
        do ib=1,(ncycle-1)/np_band+1
          nbss = (ib-1)*np_band + myrank_b + 1
          if ( nbss<=ncycle .and. nbss>= k1+1 ) then
            ms=NBAND_BLK*(nbss-1)+1
            me=min(ms+NBAND_BLK-1,MB)

            if ( ms<=me ) then
              allocate( Sunk(n0:n1,ms:me) ) ; Sunk=zero
              call get_Sunk_Mat( n0,n1,ms,me,k,s )
              call GramSchmidtGSub(ms,me,ns,ne,NBLK,k,s,NBLK,NBLK1)
              deallocate( Sunk )
            end if
          end if
        end do
      end if
    end do ! k1
    deallocate( id,ir )

    return
  END SUBROUTINE gram_schmidt_g

!---------------------------------------------------------------------------------------
  ! Gram-Schmidt orthogonalization
  RECURSIVE SUBROUTINE GramSchmidtGSub(mm1,mm2,nn1,nn2,MBLK,k,s,NBLK,NBLK1)
    implicit none
    integer,intent(IN) :: mm1,mm2,nn1,nn2,MBLK,k,s
    integer,intent(IN) :: NBLK,NBLK1
    integer :: n,ns,ne,m,ms,me,m1,m2,mm,nn,MBLKH,ierr
    integer :: ML0,i
    integer :: iii,jjj
    real(8) :: c

#ifdef _DRSDFT_
    real(8),allocatable :: utmp2(:,:),utmp(:)
    real(8) :: d
#else
    complex(8),allocatable :: utmp2(:,:),utmp(:)
    complex(8) :: d
#endif

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

          allocate( utmp2(ns:ne,ms:me) ) ; utmp2=zero

#ifdef _DRSDFT_
          call dgemm(TRANSA,TRANSB,nn,mm,ML0, -dV,unk(ML_0,ns,k,s),ML0,Sunk(ML_0,ms),ML0,zero,utmp2,nn)
#else
          call zgemm(TRANSA,TRANSB,nn,mm,ML0,-zdV,unk(ML_0,ns,k,s),ML0,Sunk(ML_0,ms),ML0,zero,utmp2,nn)
#endif

          call MPI_ALLREDUCE( MPI_IN_PLACE,utmp2,nn*mm,TYPE_MAIN,MPI_SUM,COMM_GRID,ierr )

#ifdef _DRSDFT_
          call dgemm(TRANSB,TRANSB,ML0,mm,nn,one,unk(ML_0,ns,k,s),ML0,utmp2,nn,one,unk(ML_0,ms,k,s),ML0)
#else
          call zgemm(TRANSB,TRANSB,ML0,mm,nn,one,unk(ML_0,ns,k,s),ML0,utmp2,nn,one,unk(ML_0,ms,k,s),ML0)
#endif

          deallocate( utmp2 )

          if ( ms==ne+1 ) then
            call get_gSf(unk(ML_0,ms,k,s),unk(ML_0,ms,k,s),ML_0,ML_1,k,d,0)
            c=real(d,kind=8)
            call MPI_ALLREDUCE( MPI_IN_PLACE,c,1,mpi_real8,MPI_SUM,COMM_GRID,ierr )
            c=1.d0/sqrt(c)
            unk(ML_0:ML_1,ms,k,s) = c*unk(ML_0:ML_1,ms,k,s)
          end if

        else if ( mm<=NBLK1 ) then
          allocate( utmp(NBLK1) ) ; utmp=zero
          do m=ms,me
            n = min(m-1,ne)
            if ( n-ns+1>0 ) then

#ifdef _DRSDFT_
              call dgemv(TRANSA,ML0,n-ns+1, -dV,unk(ML_0,ns,k,s),ML0,Sunk(ML_0,m),1,zero,utmp,1)
#else
              call zgemv(TRANSA,ML0,n-ns+1,-zdV,unk(ML_0,ns,k,s),ML0,Sunk(ML_0,m),1,zero,utmp,1)
#endif

              call mpi_allreduce(MPI_IN_PLACE,utmp,n-ns+1,TYPE_MAIN,mpi_sum,comm_grid,ierr)

#ifdef _DRSDFT_
              call dgemv(TRANSB,ML0,n-ns+1,one,unk(ML_0,ns,k,s),ML0,utmp,1,one,unk(ML_0,m,k,s),1)
#else
              call zgemv(TRANSB,ML0,n-ns+1,one,unk(ML_0,ns,k,s),ML0,utmp,1,one,unk(ML_0,m,k,s),1)
#endif

            end if
            if ( m==1 .or. (n==m-1 .and. m/=ns) ) then
              call get_gSf(unk(ML_0,m,k,s),unk(ML_0,m,k,s),ML_0,ML_1,k,d,0)
              c=real(d,kind=8)
              call MPI_ALLREDUCE( MPI_IN_PLACE,c,1,mpi_real8,MPI_SUM,COMM_GRID,ierr )
              c=1.d0/sqrt(c)
              unk(ML_0:ML_1,m,k,s)=c*unk(ML_0:ML_1,m,k,s)
            end if
          end do ! m
          deallocate( utmp )

        else
          MBLKH=max(MBLK/2,NBLK1)
          call GramSchmidtGSub(ms,me,ns,ne,MBLKH,k,s,NBLK,NBLK1)
        end if
      end do ! ns
    end do ! ms

  END SUBROUTINE GramSchmidtGSub


#ifdef _GS_SIMPLE_
  SUBROUTINE gram_schmidt_u(k,s)
    implicit none
    integer,intent(IN) :: k,s
    integer :: m,n,ierr
#ifdef _DRSDFT_
    real(8),allocatable :: uu(:)
#else
    complex(8),allocatable :: uu(:)
#endif
    real(8) :: c,d

    allocate( uu(MB-1) ) ; uu=zero

    allocate( Sunk(ML_0:ML_1,MB) ) ; Sunk=zero
    call get_Sunk_Mat( ML_0, ML_1, 1, MB, k, s )

    do n=1,MB

       if ( n > 1 ) then
          do m=1,n-1
#ifdef _DRSDFT_
             uu(m) = sum( unk(:,m,k,s)*Sunk(:,n) )*dV
#else
             uu(m) = sum( conjg(unk(:,m,k,s))*Sunk(:,n) )*dV
#endif
          end do
          call MPI_ALLREDUCE( MPI_IN_PLACE, uu, n-1, TYPE_MAIN &
                             ,MPI_SUM, comm_grid, ierr )
          do m=1,n-1
             unk(:,n,k,s) = unk(:,n,k,s) - unk(:,m,k,s)*uu(m)
          end do
       end if

       call get_gSf( unk(:,n,k,s), unk(:,n,k,s), ML_0, ML_1, k, uu(1), 0 )
       c=uu(1)
       call MPI_ALLREDUCE( c,d,1,MPI_REAL8,MPI_SUM,comm_grid,ierr )
       c=1.0d0/sqrt(d)
       unk(:,n,k,s) = c*unk(:,n,k,s)

    end do ! n

    deallocate( Sunk )
    deallocate( uu )

  END SUBROUTINE gram_schmidt_u
#endif

END MODULE gram_schmidt_g_module
