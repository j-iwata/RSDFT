MODULE rsdft_mpi_module

  implicit none

  PRIVATE
  PUBLIC :: rsdft_allgatherv,rsdft_allreduce

  include 'mpif.h'

#ifdef _DRSDFT_
  integer,parameter :: TYPE_MAIN=MPI_REAL8
#else
  integer,parameter :: TYPE_MAIN=MPI_COMPLEX16
#endif

CONTAINS

  SUBROUTINE rsdft_allgatherv(n0,m1,m2,mt,u,np,comm,mrnk,npa)
    implicit none
    integer,intent(IN) :: n0,m1,m2,mt,np,comm,mrnk,npa
    integer :: i,ierr,m,m0,i0,i1,ipart,j
    integer,allocatable :: ir(:),id(:),ir0(:),id0(:),id1(:)
#ifdef _DRSDFT_
    real(8),intent(INOUT) :: u(n0*mt)
    real(8),allocatable :: v(:)
#else
    complex(8),intent(INOUT) :: u(n0*mt)
    complex(8),allocatable :: v(:)
#endif

    allocate( ir(0:np-1),id(0:np-1) ) ; ir=0 ; id=0

    ir(mrnk)=n0*(m2-m1+1)
    call mpi_allgather(ir(mrnk),1,mpi_integer,ir,1,mpi_integer,comm,ierr)
    do i=0,np-1
       id(i)=sum(ir(0:i))-ir(i)
    end do

    if ( npa<=1 ) then
       call mpi_allgatherv(u(id(mrnk)+1),ir(mrnk),TYPE_MAIN &
            ,u,ir,id,TYPE_MAIN,comm,ierr)
    else
       m=(maxval(ir)*np+npa-1)/npa
       allocate( v(m) )
       allocate( ir0(0:np-1) ) ; ir0=0
       allocate( id0(0:np-1) ) ; id0=0
       allocate( id1(0:np-1) ) ; id1=0

       do j=0,np-1
          ir0(j)=(ir(j)+npa-1)/npa
       end do
       do j=0,np-1
          id0(j)=sum(ir0(0:j))-ir0(j)
       end do
       id1(:)=id(:)

       do ipart=1,npa
          if ( ipart>1 ) then
             id1(:)=id1(:)+ir0(:)
             if ( ipart==npa ) then
                do j=0,np-1
                   ir0(j)=id(j)+ir(j)-id1(j)
                end do
                do j=0,np-1
                   id0(j)=sum(ir0(0:j))-ir0(j)
                end do
             end if
          end if

          call mpi_allgatherv(u(id1(mrnk)+1),ir0(mrnk),TYPE_MAIN &
               ,v,ir0,id0,TYPE_MAIN,comm,ierr)

          do j=0,np-1
             if ( j==mrnk ) cycle
             u(id1(j)+1:id1(j)+ir0(j))=v(id0(j)+1:id0(j)+ir0(j))
          end do

       end do

       deallocate( id1,id0,ir0,v )

    end if

    deallocate( id,ir )

  END SUBROUTINE rsdft_allgatherv


  SUBROUTINE rsdft_allreduce(comm,reduce1,n_size,num)
    implicit none
    integer,intent(IN) :: comm
    integer :: i,reduce_id,reduce_ir,n_size,num,ierr
#ifdef _DRSDFT_
    real(8) :: reduce1(n_size)
    real(8),allocatable :: reduce2(:)
#else
    complex(8) :: reduce1(n_size)
    complex(8),allocatable :: reduce2(:)
#endif

    ierr=0
    if ( num<=1 ) then  
       allocate( reduce2(n_size) ) ; reduce2=(0.0d0,0.0d0)
       call mpi_allreduce(reduce1(1),reduce2(1),n_size,TYPE_MAIN &
            ,mpi_sum,comm,ierr)
    else
       allocate( reduce2(n_size/num) ) ; reduce2=(0.0d0,0.0d0)
       reduce_ir=n_size/num
       do i=0,num
          if ( i == num ) then
             if ( mod(n_size,num)==0 ) then
                exit
             else
                reduce_ir=mod(n_size,num)
             end if
          end if
          reduce_id=n_size/num*i+1
          call mpi_allreduce(reduce1(reduce_id),reduce2,reduce_ir &
               ,TYPE_MAIN,mpi_sum,comm,ierr)
          reduce1(reduce_id:reduce_id+reduce_ir-1)=reduce2(1:reduce_ir)
       end do
    end if

    deallocate( reduce2 )

    return

  END SUBROUTINE rsdft_allreduce


  SUBROUTINE rsdft_allreduce_old(comm1,comm2,reduce1,n_size,num)
    implicit none
    integer,intent(IN) :: comm1,comm2
    integer :: i,reduce_id,reduce_ir,n_size,num,ierr
#ifdef _DRSDFT_
    real(8) :: reduce1(n_size)
    real(8),allocatable :: reduce2(:)
#else
    complex(8) :: reduce1(n_size)
    complex(8),allocatable :: reduce2(:)
#endif


    ierr=0
    if ( num<=1 ) then  
       allocate( reduce2(n_size) ) ; reduce2=(0.0d0,0.0d0)
       call mpi_allreduce(reduce1(1),reduce2(1),n_size,TYPE_MAIN &
            ,mpi_sum,comm1,ierr)
       call mpi_allreduce(reduce2(1),reduce1(1),n_size,TYPE_MAIN &
            ,mpi_sum,comm2,ierr)
    else
       allocate( reduce2(n_size/num) ) ; reduce2=(0.0d0,0.0d0)
       reduce_ir=n_size/num
       do i=0,num
          if ( i==num ) then
             if ( mod(n_size,num)==0 ) then
                exit
             else
                reduce_ir=mod(n_size,num)
             end if
          end if
          reduce_id=n_size/num*i+1
          call mpi_allreduce(reduce1(reduce_id),reduce2,reduce_ir &
               ,TYPE_MAIN,mpi_sum,comm1,ierr)
          reduce1(reduce_id:reduce_id+reduce_ir-1)=reduce2(1:reduce_ir)
          call mpi_allreduce(reduce1(reduce_id),reduce2,reduce_ir &
               ,TYPE_MAIN,mpi_sum,comm2,ierr)
          reduce1(reduce_id:reduce_id+reduce_ir-1)=reduce2(1:reduce_ir)
       end do
    end if

    deallocate( reduce2 )

    return

  END SUBROUTINE rsdft_allreduce_old


END MODULE rsdft_mpi_module
