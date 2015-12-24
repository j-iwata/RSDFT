MODULE mixing_pulay_module

  use fft_module
  use grid_module
  use bb_module
  use parallel_module, only: comm_grid

  implicit none

  PRIVATE
  PUBLIC :: pulay_mixing, pulay_g_mixing
  PUBLIC :: set_init_pulay

  include 'mpif.h' 

  complex(8),allocatable :: zwork0(:,:,:), zwork1(:,:,:)

CONTAINS

  SUBROUTINE set_init_pulay( fd_in, fz_out )
    implicit none
    real(8),intent(IN) :: fd_in(:,:)
    complex(8),intent(OUT) :: fz_out(:,:)
    integer :: s
    call pulay_sub1_fft( fd_in, fz_out )
  END SUBROUTINE set_init_pulay

  SUBROUTINE pulay_mixing( m,n,comm,mmix_count,mmix,beta,f,Xin,Xou )

    implicit none
    integer,intent(IN)       :: m,n,comm,mmix
    integer,intent(INOUT)    :: mmix_count
    real(8),intent(IN)       :: beta
    real(8),intent(INOUT)    :: f(m,n)
    complex(8),intent(INOUT) :: Xin(m,n,mmix)
    complex(8),intent(INOUT) :: Xou(m,n,mmix)
    integer :: s,mmix0,ierr,i,i0,j0,mm,j
    integer,allocatable :: ipiv(:)
    real(8),allocatable :: rwork(:,:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8) :: zc
    complex(8),allocatable :: A0(:,:),A1(:,:),b1(:),X(:),Y(:)

    Xou(:,:,mmix) = f(:,:)
!    call pulay_sub1_fft( f, Xou(:,:,mmix) )

    mmix_count = mmix_count + 1

    mmix0 = min( mmix_count, mmix )

    if ( mmix == 1 .or. mmix0 < mmix ) then
       do i=2,mmix
          Xin(:,:,i-1)=Xin(:,:,i)
          Xou(:,:,i-1)=Xou(:,:,i)
       end do
       Xin(:,:,mmix) = Xin(:,:,mmix) + beta*( Xou(:,:,mmix)-Xin(:,:,mmix) )
!       Xin(:,:,mmix) = Xou(:,:,mmix)
       f(:,:) = Xin(:,:,mmix)
!       call pulay_sub2_fft( Xin(:,:,mmix), f )
       return
    end if

    allocate( ipiv(mmix0) )
    allocate( X(m),Y(m) )
    allocate( b1(mmix0) )
    allocate( A1(mmix0,mmix0) )
    allocate( A0(mmix0,mmix0) )

    b1(:)   = zero
    A1(:,:) = zero
    A0(:,:) = zero
    mm      = size(A1)

    do j0=1 ,mmix0
    do i0=j0,mmix0
       i=mmix-mmix0+i0
       j=mmix-mmix0+j0
       A0(i0,j0)=sum( conjg(Xou(:,:,i)-Xin(:,:,i)) &
                          *(Xou(:,:,j)-Xin(:,:,j)) )
       A0(j0,i0)=conjg( A0(i0,j0) )
    end do
    end do

    call mpi_allreduce(A0,A1,mm,mpi_complex16,mpi_sum,comm,ierr)

    b1(1:mmix0) = (1.d0,0.d0)
    A0(:,:)     = A1(:,:)

    call zgesv(mmix0,1,A1,mmix0,ipiv,b1,mmix0,ierr)

    zc=1.d0/sum( b1(1:mmix0) )
    b1(1:mmix0)=zc*b1(1:mmix0)

    do s=1,n

       X(:)=zero
       Y(:)=zero

       do i0=1,mmix0
          i=mmix-mmix0+i0
          X(:)=X(:)+b1(i0)*Xin(:,s,i)
          Y(:)=Y(:)+b1(i0)*Xou(:,s,i)
       end do

       do i=max(1,mmix-mmix0),mmix-1
          Xin(:,s,i)=Xin(:,s,i+1)
          Xou(:,s,i)=Xou(:,s,i+1)
       end do

       Xin(:,s,mmix) = X(:) + beta*( Y(:)-X(:) )

    end do ! s

    f(:,:) = real( Xin(:,:,mmix) )
    !call pulay_sub2_fft( Xin(:,:,mmix), f )

    deallocate( A0,A1,b1,Y,X,ipiv )
    return

  END SUBROUTINE pulay_mixing


  SUBROUTINE pulay_g_mixing( m,n,comm,mmix_count,mmix,beta,f,Xin,Xou )

    implicit none
    integer,intent(IN)       :: m,n,comm,mmix
    integer,intent(INOUT)    :: mmix_count
    real(8),intent(IN)       :: beta
    real(8),intent(INOUT)    :: f(m,n)
    complex(8),intent(INOUT) :: Xin(m,n,mmix)
    complex(8),intent(INOUT) :: Xou(m,n,mmix)
    integer :: s,mmix0,ierr,i,i0,j0,mm,j
    integer,allocatable :: ipiv(:)
    real(8),allocatable :: rwork(:,:), wq(:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8) :: zc
    complex(8),allocatable :: A0(:,:),A1(:,:),b1(:),X(:),Y(:)

    call write_border( 1, " pulay_g_mixing(start)" )

!    Xou(:,:,mmix) = f(:,:)
    call pulay_sub1_fft( f, Xou(:,:,mmix) )

    allocate( wq(m) ) ; wq=0.0d0
    call get_Kerker_weight( wq )

    mmix_count = mmix_count + 1

    mmix0 = min( mmix_count, mmix )

    if ( mmix == 1 .or. mmix0 < mmix ) then
       do i=2,mmix
          Xin(:,:,i-1)=Xin(:,:,i)
          Xou(:,:,i-1)=Xou(:,:,i)
       end do
       Xin(:,:,mmix) = Xin(:,:,mmix) + beta*( Xou(:,:,mmix)-Xin(:,:,mmix) )
!       Xin(:,:,mmix) = Xou(:,:,mmix)
!       f(:,:) = Xin(:,:,mmix)
       call pulay_sub2_fft( Xin(:,:,mmix), f )
       return
    end if

    allocate( ipiv(mmix0) )
    allocate( X(m),Y(m) )
    allocate( b1(mmix0) )
    allocate( A1(mmix0,mmix0) )
    allocate( A0(mmix0,mmix0) )

    b1(:)   = zero
    A1(:,:) = zero
    A0(:,:) = zero
    mm      = size(A1)

    do j0=1 ,mmix0
    do i0=j0,mmix0
       i=mmix-mmix0+i0
       j=mmix-mmix0+j0
!       A0(i0,j0)=sum( conjg(Xou(:,:,i)-Xin(:,:,i)) &
!                          *(Xou(:,:,j)-Xin(:,:,j)) )
       do s=1,n
          A0(i0,j0)=A0(i0,j0)+sum( conjg(Xou(:,s,i)-Xin(:,s,i)) &
                                       *(Xou(:,s,j)-Xin(:,s,j))/wq(:) )
       end do
       A0(j0,i0)=conjg( A0(i0,j0) )
    end do
    end do

    call mpi_allreduce(A0,A1,mm,mpi_complex16,mpi_sum,comm,ierr)

    b1(1:mmix0) = (1.d0,0.d0)
    A0(:,:)     = A1(:,:)

    call zgesv(mmix0,1,A1,mmix0,ipiv,b1,mmix0,ierr)

    zc=1.d0/sum( b1(1:mmix0) )
    b1(1:mmix0)=zc*b1(1:mmix0)

    do s=1,n

       X(:)=zero
       Y(:)=zero

       do i0=1,mmix0
          i=mmix-mmix0+i0
          X(:)=X(:)+b1(i0)*Xin(:,s,i)
          Y(:)=Y(:)+b1(i0)*Xou(:,s,i)
       end do

       do i=max(1,mmix-mmix0),mmix-1
          Xin(:,s,i)=Xin(:,s,i+1)
          Xou(:,s,i)=Xou(:,s,i+1)
       end do

       Xin(:,s,mmix) = X(:) + beta*( Y(:)-X(:) )

    end do ! s

    !f(:,:) = real( Xin(:,:,mmix) )
    call pulay_sub2_fft( Xin(:,:,mmix), f )

    deallocate( A0,A1,b1,Y,X,ipiv )
    deallocate( wq )

    call write_border( 1, " pulay_g_mixing(end)" )

    return

  END SUBROUTINE pulay_g_mixing


  SUBROUTINE pulay_sub1_fft( fr, fz )

    implicit none
    real(8),intent(IN) :: fr(:,:)
    complex(8),intent(OUT) :: fz(:,:)
    integer :: s

    call init_fft

    do s=1,size(fr,2)
       call d1_to_z3_fft( fr(:,s), zwork0 )
       call forward_fft( zwork0, zwork1 )
       call z3_to_z1_fft( zwork0, fz(:,s) )
    end do

    call finalize_fft

    deallocate( zwork1, zwork0 )

  END SUBROUTINE pulay_sub1_fft

  SUBROUTINE pulay_sub2_fft( fz, fr )

    implicit none
    complex(8),intent(IN) :: fz(:,:)
    real(8),intent(OUT) :: fr(:,:)
    integer :: s

    call init_fft

    do s=1,size(fr,2)
       call z1_to_z3_fft( fz(:,s), zwork0 )
       call backward_fft( zwork0, zwork1 )
       call z3_to_d1_fft( zwork0, fr(:,s) )
    end do

    call finalize_fft

    deallocate( zwork0, zwork1 )

  END SUBROUTINE pulay_sub2_fft
       

  SUBROUTINE get_Kerker_weight( wq )
    implicit none
    real(8),intent(OUT) :: wq(:)
    integer :: i1,i2,i3,i,j1,j2,j3,k1,k2,k3,ML1,ML2,ML3
    real(8) :: gx,gy,gz,gg,gg_min
    type(grid) :: rgrid
    call get_range_rgrid( rgrid )
    ML1=rgrid%g3%x%size_global
    ML2=rgrid%g3%y%size_global
    ML3=rgrid%g3%z%size_global
    wq=0.0d0
    gg_min=1.d100
    i=0
    do i3=rgrid%g3%z%head,rgrid%g3%z%tail
    do i2=rgrid%g3%y%head,rgrid%g3%y%tail
    do i1=rgrid%g3%x%head,rgrid%g3%x%tail
       i=i+1
       j1 = i1-ML1 ; if ( j1 < -(ML1-1)/2 ) j1=i1
       j2 = i2-ML2 ; if ( j2 < -(ML2-1)/2 ) j2=i2
       j3 = i3-ML3 ; if ( j3 < -(ML3-1)/2 ) j3=i3
       k1 = mod(j1+ML1,ML1)
       k2 = mod(j2+ML2,ML2)
       k3 = mod(j3+ML3,ML3)
       if ( k1/=i1 .or. k2/=i2 .or. k3/=i3 ) then
          write(*,*) k1,k2,k3,i1,i2,i3
          call stop_program_f( "stop@get_wq_fft" )
       end if
       gx=j1*bb(1,1)+j2*bb(1,2)+j3*bb(1,3)
       gy=j1*bb(2,1)+j2*bb(2,2)+j3*bb(2,3)
       gz=j1*bb(3,1)+j2*bb(3,2)+j3*bb(3,3)
       gg=gx*gx+gy*gy+gz*gz
       if ( gg > 1.d-10 ) gg_min=min(gg,gg_min)
       wq(i)=gg
    end do
    end do
    end do
    call mpi_allreduce(gg_min,gg,1,MPI_REAL8,MPI_MIN,comm_grid,i)
    do i=1,size(wq)
       wq(i) = max( 0.5d0, wq(i)/(wq(i)+gg) )
    end do
  END SUBROUTINE get_Kerker_weight


END MODULE mixing_pulay_module
