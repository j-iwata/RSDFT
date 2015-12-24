MODULE mixing_broyden_module

  implicit none

  PRIVATE
  PUBLIC :: broyden_mixing

CONTAINS


  SUBROUTINE broyden_mixing(ML0,MSP,comm,iter,mmix,alpha0,g,zXin,zXou)

    implicit none

    include 'mpif.h'

    integer,intent(IN)       :: ML0,MSP,comm,mmix
    integer,intent(INOUT)    :: iter
    real(8),intent(IN)       :: alpha0
    real(8),intent(INOUT)    :: g(ML0,MSP)
    complex(8),intent(INOUT) :: zXin(ML0,MSP,mmix)
    complex(8),intent(INOUT) :: zXou(ML0,MSP,mmix)
    integer :: i,j,k,l,mmix0,ierr,mmm
    integer,allocatable :: ipiv(:)
    real(8),allocatable :: F(:,:,:)
    real(8),allocatable :: dF(:,:,:)
    real(8),allocatable :: u(:,:,:)
    real(8),allocatable :: dX(:,:,:)
    real(8),allocatable :: a(:,:),beta(:,:),cm(:),gamma(:)
    real(8),allocatable :: w(:),Xnew(:,:)
    real(8),allocatable :: Xin(:,:,:),Xou(:,:,:)
    real(8),allocatable :: s0(:),s1(:),work(:)
    real(8) :: alpha,dFnorm,c1

    call write_border( 1, " broyden_mixing(start)" )

    zXou(:,:,mmix) = g(:,:)

    iter = iter + 1

    mmix0 = min( iter, mmix )

!    if ( mmix == 1 .or. mmix0 < mmix ) then
!       do k=2,mmix
!          zXin(:,:,k-1) = zXin(:,:,k)
!          zXou(:,:,k-1) = zXou(:,:,k)
!       end do
!       alpha=1.0d0-alpha0
!       zXin(:,:,mmix) = alpha*zXin(:,:,mmix) + alpha0*zXou(:,:,mmix)
!       g(:,:) = zXin(:,:,mmix)
!       return
!    end if
    if ( mmix == 1 .or. mmix0 == 1 ) then
       do k=2,mmix
          zXin(:,:,k-1) = zXin(:,:,mmix)
          zXou(:,:,k-1) = zXou(:,:,mmix)
       end do
       alpha=1.0d0-alpha0
       zXin(:,:,mmix) = alpha*zXin(:,:,mmix) + alpha0*zXou(:,:,mmix)
       g(:,:) = zXin(:,:,mmix)
       return
    else if ( mmix0 < mmix ) then
       do k=2,mmix0
          zXin(:,:,k-1) = zXin(:,:,k)
          zXou(:,:,k-1) = zXou(:,:,k)
       end do
       zXin(:,:,mmix0) = zXin(:,:,mmix)
       zXou(:,:,mmix0) = zXou(:,:,mmix)
    end if

    allocate( Xin(ML0,MSP,mmix)     )
    allocate( Xou(ML0,MSP,mmix)     )
    allocate( F(ML0,MSP,mmix0)      )
    allocate( dF(ML0,MSP,mmix0)     )
    allocate( s0(mmix0)             )
    allocate( s1(mmix0)             )
    allocate( w(0:mmix0)            )
    allocate( dX(ML0,MSP,mmix0)     )
    allocate( u(ML0,MSP,mmix0)      )
    allocate( cm(mmix0)             )
    allocate( a(mmix0-1,mmix0-1)    )
    allocate( beta(mmix0-1,mmix0-1) )
    allocate( gamma(mmix0)          )
    allocate( Xnew(ML0,MSP)         )
    allocate( ipiv(mmix0-1)         )
    allocate( work(mmix0-1)         )

    Xin(:,:,:) = real( zXin(:,:,:) )
    Xou(:,:,:) = real( zXou(:,:,:) )

    do k=1,mmix0
       F(:,:,k) = Xou(:,:,k) - Xin(:,:,k)
    end do

! Weights

    alpha      = alpha0
    w(0)       = alpha0
    w(1:mmix0) = 1.0d0
    do k=1,mmix0
       s0(k) = sum( F(:,:,k)*F(:,:,k) )
    end do
    call MPI_ALLREDUCE(s0,s1,mmix0,MPI_REAL8,MPI_SUM,comm,ierr)
    do k=1,mmix0
       if ( s1(k) <= 0.0d0 ) then
          w(k) = 1.0d0
       else
          w(k)=1.0d0/sqrt( s1(k) )
          if ( w(k) > 1.0d0 ) w(k) = 1.0d0
       end if
    end do

    do k=1,mmix0-1
       dF(:,:,k) = F(:,:,k+1) - F(:,:,k)
       dFnorm    = sum( dF(:,:,k)*dF(:,:,k) )
       call mpi_allreduce(dFnorm,c1,1,mpi_real8,mpi_sum,comm,ierr)
       dFnorm    = 1.0d0/sqrt(c1)
       dF(:,:,k) = dF(:,:,k)*dFnorm
       dX(:,:,k) = ( Xin(:,:,k+1)-Xin(:,:,k) )*dFnorm
    end do

    do k=1,mmix0-1
       u(:,:,k) = alpha*dF(:,:,k) + dX(:,:,k)
    end do

    do k=1,mmix0-1
       s0(k) = w(k)*sum( dF(:,:,k)*F(:,:,mmix0) )
    end do
    call mpi_allreduce(s0,cm,mmix0-1,MPI_REAL8,MPI_SUM,comm,ierr)

    do l=1,mmix0-1
       do k=1,l
          s0(k)=sum( dF(:,:,k)*dF(:,:,l) )
       end do
       call mpi_allreduce(s0,s1,l,MPI_REAL8,MPI_SUM,comm,ierr)
       do k=1,l
          a(k,l) = w(k)*w(l)*s1(k)
          a(l,k) = a(k,l)
       end do
       a(l,l) = a(l,l) + w(0)*w(0)
    end do

    beta(:,:)=a(:,:)
    call dgetrf( mmix0-1, mmix0-1, beta, mmix0-1, ipiv, ierr )
    call dgetri( mmix0-1, beta, mmix0-1, ipiv, work, mmix0-1, ierr )

    gamma(:)=0.0d0
    do l=1,mmix0-1
       do k=1,mmix0-1
          gamma(l) = gamma(l) + cm(k)*beta(k,l)
       end do
    end do

    Xnew(:,:) = Xin(:,:,mmix0) + alpha*F(:,:,mmix0)
    do k=1,mmix0-1
       Xnew(:,:) = Xnew(:,:) - w(k)*gamma(k)*u(:,:,k)
    end do

! ---

    if ( mmix0 == mmix ) then

       do k=1,mmix-1
          zXin(:,:,k) = Xin(:,:,k+1)
          zXou(:,:,k) = Xou(:,:,k+1)
       end do

    end if

    zXin(:,:,mmix) = Xnew(:,:)

    g(:,:) = Xnew(:,:)

    deallocate( work )
    deallocate( ipiv )
    deallocate( Xnew )
    deallocate( gamma )
    deallocate( beta )
    deallocate( a )
    deallocate( cm )
    deallocate( u )
    deallocate( dX )
    deallocate( w )
    deallocate( s1 )
    deallocate( s0 )
    deallocate( dF )
    deallocate( F )
    deallocate( Xou )
    deallocate( Xin )

    call write_border( 1, " broyden_mixing(end)" )

    return
  END SUBROUTINE broyden_mixing


END MODULE mixing_broyden_module
