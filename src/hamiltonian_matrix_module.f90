MODULE hamiltonian_matrix_module

  use kinetic_sol_module
  use construct_matrix_ps_nloc2_module
  use localpot_module
  use symmetry_module
  use rgrid_module

  implicit none

  PRIVATE
  PUBLIC :: construct_hamiltonian_matrix

#ifdef _DRSDFT_
  real(8),parameter :: zero=0.0d0
  real(8),allocatable :: Hmat(:,:)
  real(8),allocatable :: w1(:,:), w2(:,:)
#else
  complex(8),parameter :: zero=(0.0d0,0.0d0)
  complex(8),allocatable :: Hmat(:,:)
  complex(8),allocatable :: w1(:,:), w2(:,:)
#endif

  integer :: ML

  real(8),allocatable :: SymMat(:,:)

CONTAINS


  SUBROUTINE construct_hamiltonian_matrix( ML_in )
    implicit none
    integer,intent(IN) :: ML_in
    integer :: i,n,k,s,isym,i1,i2,i3,j,i0
    real(8) :: c,r0(3),r1(3),r

    write(*,'(a50," construct_hamiltonian_matrix")') repeat("-",50)

    ML = ML_in

    c = 1.0d0/nsym

    allocate( Hmat(ML,ML)   ) ; Hmat=zero
    allocate( SymMat(ML,ML) ) ; SymMat=0.0d0
    allocate( w1(ML,ML)     ) ; w1=zero
    allocate( w2(ML,ML)     ) ; w2=zero

    do s=1,1
    do k=1,1

       Hmat(:,:) = zero

       call construct_matrix_kinetic_sol( k, ML, Hmat )
!       call construct_matrix_localpot( s, ML, Hmat )
!       call construct_matrix_ps_nloc2( k, ML, Hmat )

!       n=count(abs(Hmat)>1.d-10)
!       write(*,*) n,dble(n)/dble(ML)
!       do i=1,ML
!          write(*,'(1x,2i5,2f20.15)') i,count(abs(Hmat(i,:))>1.d-10) &
!          ,Hmat(i,i)
!       end do
!       i1=Ngrid(1)/2
!       i2=Ngrid(2)/2
!       i3=Ngrid(3)/2
!       i0=1+i1+i2*Ngrid(1)+i3*Ngrid(1)*Ngrid(2)
!       write(*,*) i0,i1,i2,i3
!       j=0
!       do i3=0,Ngrid(3)-1
!       do i2=0,Ngrid(2)-1
!       do i1=0,Ngrid(1)-1
!          i=1+i1+i2*Ngrid(1)+i3*Ngrid(1)*Ngrid(2)
!          if ( abs(Hmat(i0,i)) > 1.d-10 ) then
!             j=j+1
!             w1(1,j)=i1
!             w1(2,j)=i2
!             w1(3,j)=i3
!             w1(4,j)=Hmat(i0,i)
!             if ( i == i0 ) then
!                r0(1)=i1
!                r0(2)=i2
!                r0(3)=i3
!             end if
!          end if
!       end do
!       end do
!       end do
!       do i=1,j
!          r1(1:3)=real(w1(1:3,i))
!          r=sqrt(sum((r1-r0)**2))
!          write(*,'(1x,4i5,2f20.15)') i,nint(real(w1(1:3,i))),real(w1(4,i)),r
!       end do

       if ( isymmetry /= 0 ) then
          w2(:,:)=zero
          do isym=1,nsym
             write(*,'(1x,"isym/nsym=",i3," /",i3)') isym,nsym
             call construct_matrix_symmetry( isym, ML, SymMat )
             w1(:,:) = matmul( Hmat, transpose(SymMat) )
             w2(:,:) = w2(:,:) + matmul( SymMat, w1 )
          end do
          Hmat(:,:) = c*w2(:,:)
       end if
!       n=count(abs(Hmat)>1.d-10)
!       write(*,*) n,dble(n)/dble(ML)
!       do i=1,ML
!          write(*,'(1x,2i5,2f20.15)') i,count(abs(Hmat(i,:))>1.d-10) &
!          ,Hmat(i,i)
!       end do
!       j=0
!       do i3=0,Ngrid(3)-1
!       do i2=0,Ngrid(2)-1
!       do i1=0,Ngrid(1)-1
!          i=1+i1+i2*Ngrid(1)+i3*Ngrid(1)*Ngrid(2)
!          if ( abs(Hmat(i0,i)) > 1.d-10 ) then
!             j=j+1
!             w1(1,j)=i1
!             w1(2,j)=i2
!             w1(3,j)=i3
!             w1(4,j)=Hmat(i0,i)
!          end if
!       end do
!       end do
!       end do
!       do i=1,j
!          r1(1:3)=real(w1(1:3,i))
!          r=sqrt(sum((r1-r0)**2))
!          write(*,'(1x,4i5,2f20.15)') i,nint(real(w1(1:3,i))),real(w1(4,i)),r
!       end do

!       call construct_matrix_kinetic_sol( k, ML, Hmat )
       call construct_matrix_localpot( s, ML, Hmat )
       call construct_matrix_ps_nloc2( k, ML, Hmat )

       call diag_Hmat

    end do ! k
    end do ! s

    deallocate( w2 )
    deallocate( w1 )
    deallocate( SymMat )
    deallocate( Hmat )

  END SUBROUTINE construct_hamiltonian_matrix


  SUBROUTINE diag_Hmat
    implicit none
    integer :: LWORK,LRWORK,ierr,i
    real(8),allocatable :: rwork(:),eigval(:)
    complex(8),allocatable :: work(:)
    LWORK=2*ML-1
    LRWORK=3*ML-2
    allocate( work(LWORK),rwork(LRWORK) )
    allocate( eigval(ML) )
    call zheev('N','L',ML,Hmat,ML,eigval,work,LWORK,rwork,ierr)
    do i=1,10
       write(*,'(1x,i5,f20.15)') i,eigval(i)
    end do
    do i=1,3
       write(*,'(1x,5x,".")')
    end do
    do i=ML-9,ML
       write(*,'(1x,i5,f20.15)') i,eigval(i)
    end do
    deallocate( eigval )
    deallocate( rwork,work )
  END SUBROUTINE diag_Hmat


END MODULE hamiltonian_matrix_module
