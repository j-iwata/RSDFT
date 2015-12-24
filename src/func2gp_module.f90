MODULE func2gp_module

 !use esm_rgrid_module
  use rgrid_module
  use rgrid_mol_module, LL_mol => LL
  use parallel_module
  use fd_module
  
  implicit none

  PRIVATE
  PUBLIC :: func2gp_c_esm,func2gp_r_esm,func2gp_r,func2gp_r_mol &
       ,func2gp_devr,func2gp_diff

CONTAINS


  SUBROUTINE func2gp_c_esm(unit,n1,n2,f)
    implicit none
    integer,intent(IN) :: unit,n1,n2
    complex(8),intent(IN) :: f(n1:n2)
    integer :: i,i1,i2,i3,a1,a2,a3,b1,b2,b3,ierr,m
    real(8),allocatable :: w(:,:,:,:)
    real(8),parameter :: zero=(0.d0,0.d0)

    a1 = -Ngrid(1)/2
    b1 =  Ngrid(1)/2+1
    a2 = -Ngrid(2)/2
    b2 =  Ngrid(2)/2+1
    a3 = -Ngrid(3)/2
    b3 =  Ngrid(3)/2+1

    allocate( w(a1:b1,a2:b2,a3:b3,2) )
    w=zero

    do i=n1,n2
       !i1=LL_ESM(1,i)
       !i2=LL_ESM(2,i)
       !i3=LL_ESM(3,i)
       w(i1,i2,i3,1) = abs( f(i) )**2
    end do

    m=size( w(:,:,:,1) )
    call mpi_allreduce(w(:,:,:,1),w(:,:,:,2),m,mpi_real8,mpi_sum,comm_grid,ierr)

    if ( myrank == 0 ) then
       rewind unit
       do i=a1,b1
          write(unit,'(1x,4g20.10)') i*Hgrid(1),w(i,0,0,2),w(0,i,0,2),w(0,0,i,2)
       end do
    end if

    deallocate( w )

  END SUBROUTINE func2gp_c_esm

  SUBROUTINE func2gp_r_esm(unit,n1,n2,f)
    implicit none
    integer,intent(IN) :: unit,n1,n2
    real(8),intent(IN) :: f(n1:n2)
    integer :: i,i1,i2,i3,a1,a2,a3,b1,b2,b3,ierr,m
    real(8),allocatable :: w(:,:,:,:)
    real(8),parameter :: zero=(0.d0,0.d0)

    a1 = -Ngrid(1)/2
    b1 =  Ngrid(1)/2+1
    a2 = -Ngrid(2)/2
    b2 =  Ngrid(2)/2+1
    a3 = -Ngrid(3)/2
    b3 =  Ngrid(3)/2+1

    allocate( w(a1:b1,a2:b2,a3:b3,2) )
    w=zero

    do i=n1,n2
!       i1=LL_ESM(1,i)
!       i2=LL_ESM(2,i)
!       i3=LL_ESM(3,i)
       w(i1,i2,i3,1) = f(i) 
    end do

    m=size( w(:,:,:,1) )
    call mpi_allreduce(w(:,:,:,1),w(:,:,:,2),m,mpi_real8,mpi_sum,comm_grid,ierr)

    if ( myrank == 0 ) then
       rewind unit
       do i=a1,b1
          write(unit,'(1x,4g20.10)') i*Hgrid(1),w(i,0,0,2),w(0,i,0,2),w(0,0,i,2)
       end do
    end if

    deallocate( w )

  END SUBROUTINE func2gp_r_esm


  SUBROUTINE func2gp_r(unit,n1,n2,f)
    implicit none
    integer,intent(IN) :: unit,n1,n2
    real(8),intent(IN) :: f(n1:n2)
    integer :: i,i1,i2,i3,a1,a2,a3,b1,b2,b3,ierr,m
    real(8),allocatable :: w(:,:,:,:)
    real(8),parameter :: zero=0.0d0

    a1 = 0
    b1 = Ngrid(1)-1
    a2 = 0
    b2 = Ngrid(2)-1
    a3 = 0
    b3 = Ngrid(3)-1

    allocate( w(a1:b1,a2:b2,a3:b3,2) )
    w=zero

    i=n1-1
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       i=i+1
       w(i1,i2,i3,1) = f(i) 
    end do
    end do
    end do

    m=size( w(:,:,:,1) )
    call mpi_allreduce(w(:,:,:,1),w(:,:,:,2),m,mpi_real8,mpi_sum,comm_grid,ierr)

    if ( myrank == 0 ) then
       rewind unit
       do i=a1,b1
          write(unit,'(1x,4g20.10)') i*Hgrid(1),w(i,0,0,2),sum(w(i,:,:,2))
       end do
       write(unit,*)
       write(unit,*)
       do i=a2,b2
          write(unit,'(1x,4g20.10)') i*Hgrid(2),w(0,i,0,2),sum(w(:,i,:,2))
       end do
       write(unit,*)
       write(unit,*)
       do i=a3,b3
          write(unit,'(1x,4g20.10)') i*Hgrid(3),w(0,0,i,2),sum(w(:,:,i,2))
       end do
    end if

    deallocate( w )

  END SUBROUTINE func2gp_r


  SUBROUTINE func2gp_r_mol(unit,n1,n2,f)
    implicit none
    integer,intent(IN) :: unit,n1,n2
    real(8),intent(IN) :: f(n1:n2)
    integer :: i,i1,i2,i3,a1,a2,a3,b1,b2,b3,ierr,m
    real(8),allocatable :: w(:,:,:,:)
    real(8),parameter :: zero=0.0d0

    a1 =-Ngrid(1)
    b1 = Ngrid(1)
    a2 =-Ngrid(2)
    b2 = Ngrid(2)
    a3 =-Ngrid(3)
    b3 = Ngrid(3)

    allocate( w(a1:b1,a2:b2,a3:b3,2) )
    w=zero

    do i=n1,n2
       w(LL_mol(1,i),LL_mol(2,i),LL_mol(3,i),1) = f(i) 
    end do

    m=size( w(:,:,:,1) )
    call mpi_allreduce(w(:,:,:,1),w(:,:,:,2),m,mpi_real8,mpi_sum,comm_grid,ierr)

    if ( myrank == 0 ) then
       rewind unit
       do i=a1,b1
          write(unit,'(1x,4g20.10)') i*Hgrid(1),w(i,0,0,2),sum(w(i,:,:,2))
       end do
       write(unit,*)
       write(unit,*)
       do i=a2,b2
          write(unit,'(1x,4g20.10)') i*Hgrid(2),w(0,i,0,2),sum(w(:,i,:,2))
       end do
       write(unit,*)
       write(unit,*)
       do i=a3,b3
          write(unit,'(1x,4g20.10)') i*Hgrid(3),w(0,0,i,2),sum(w(:,:,i,2))
       end do
    end if

    deallocate( w )

  END SUBROUTINE func2gp_r_mol


  SUBROUTINE func2gp_devr( unit, n1, n2, Md, f )
    implicit none
    integer,intent(IN) :: unit,n1,n2,Md
    real(8),intent(IN) :: f(n1:n2)
    integer :: a1,a2,a3,b1,b2,b3,ierr,m
    integer :: i1,i2,i3,j1,j2,j3,i
    real(8) :: c1,c2,c3
    real(8),allocatable :: w(:,:,:,:),u(:,:,:)
    real(8),parameter :: zero=0.0d0
    type(fd) :: nabla, lapla

    a1 = 0
    b1 = Ngrid(1)-1
    a2 = 0
    b2 = Ngrid(2)-1
    a3 = 0
    b3 = Ngrid(3)-1

    allocate( w(a1-Md:b1+Md,a2-Md:b2+Md,a3-Md:b3+Md,0:3) )
    w=zero
    allocate( u(a1-Md:b1+Md,a2-Md:b2+Md,a3-Md:b3+Md) )
    u=zero

    i=n1-1
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       i=i+1
       w(i1,i2,i3,1) = f(i) 
    end do
    end do
    end do

    m=size( w(:,:,:,1) )
    call mpi_allreduce(w(:,:,:,1),w(:,:,:,0),m,mpi_real8,mpi_sum,comm_grid,ierr)

    do i2=a2,b2
    do i1=a1,b1
       do i3=a3-Md,a3-1
          j3=mod(i3+Ngrid(3),Ngrid(3))
          w(i1,i2,i3,0) = w(i1,i2,j3,0)
       end do
       do i3=b3+1,b3+Md
          j3=mod(i3+Ngrid(3),Ngrid(3))
          w(i1,i2,i3,0) = w(i1,i2,j3,0)
       end do
    end do
    end do
    do i3=a3,b3
    do i1=a1,b1
       do i2=a2-Md,a2-1
          j2=mod(i2+Ngrid(2),Ngrid(2))
          w(i1,i2,i3,0) = w(i1,j2,i3,0)
       end do
       do i2=b2+1,b2+Md
          j2=mod(i2+Ngrid(2),Ngrid(2))
          w(i1,i2,i3,0) = w(i1,j2,i3,0)
       end do
    end do
    end do
    do i3=a3,b3
    do i2=a2,b2
       do i1=a1-Md,a1-1
          j1=mod(i1+Ngrid(1),Ngrid(1))
          w(i1,i2,i3,0) = w(j1,i2,i3,0)
       end do
       do i1=b1+1,b1+Md
          j1=mod(i1+Ngrid(1),Ngrid(1))
          w(i1,i2,i3,0) = w(j1,i2,i3,0)
       end do
    end do
    end do

    call construct_nabla_fd( nabla )

    u(:,:,:) = w(:,:,:,0)**(4.0d0/3.0d0)

    w(:,:,:,1:3)=zero
    do m=1,Md
       c1=nabla%coef(m) / Hgrid(1)
       c2=nabla%coef(m) / Hgrid(2)
       c3=nabla%coef(m) / Hgrid(3)
       do i3=a3,b3
       do i2=a2,b2
       do i1=a1,b1
          w(i1,i2,i3,1)=w(i1,i2,i3,1)+c1*(w(i1+m,i2,i3,0)-w(i1-m,i2,i3,0))
          w(i1,i2,i3,2)=w(i1,i2,i3,2)+c2*(w(i1,i2+m,i3,0)-w(i1,i2-m,i3,0))
          w(i1,i2,i3,3)=w(i1,i2,i3,3)+c3*(w(i1,i2,i3+m,0)-w(i1,i2,i3-m,0))
!          w(i1,i2,i3,1)=w(i1,i2,i3,1)+c1*( (w(i1+m,i2,i3,0)/u(i1,i2,i3)) &
!                                          -(w(i1-m,i2,i3,0)/u(i1,i2,i3)) )
!          w(i1,i2,i3,2)=w(i1,i2,i3,2)+c2*( (w(i1,i2+m,i3,0)/u(i1,i2,i3)) &
!                                          -(w(i1,i2-m,i3,0)/u(i1,i2,i3)) )
!          w(i1,i2,i3,3)=w(i1,i2,i3,3)+c3*( (w(i1,i2,i3+m,0)/u(i1,i2,i3)) &
!                                          -(w(i1,i2,i3-m,0)/u(i1,i2,i3)) )
       end do
       end do
       end do
    end do

    do i=1,3
       do i2=a2,b2
       do i1=a1,b1
          do i3=a3-Md,a3-1
             j3=mod(i3+Ngrid(3),Ngrid(3))
             w(i1,i2,i3,i) = w(i1,i2,j3,i)
          end do
          do i3=b3+1,b3+Md
             j3=mod(i3+Ngrid(3),Ngrid(3))
             w(i1,i2,i3,i) = w(i1,i2,j3,i)
          end do
       end do
       end do
       do i3=a3,b3
       do i1=a1,b1
          do i2=a2-Md,a2-1
             j2=mod(i2+Ngrid(2),Ngrid(2))
             w(i1,i2,i3,i) = w(i1,j2,i3,i)
          end do
          do i2=b2+1,b2+Md
             j2=mod(i2+Ngrid(2),Ngrid(2))
             w(i1,i2,i3,i) = w(i1,j2,i3,i)
          end do
       end do
       end do
       do i3=a3,b3
       do i2=a2,b2
          do i1=a1-Md,a1-1
             j1=mod(i1+Ngrid(1),Ngrid(1))
             w(i1,i2,i3,i) = w(j1,i2,i3,i)
          end do
          do i1=b1+1,b1+Md
             j1=mod(i1+Ngrid(1),Ngrid(1))
             w(i1,i2,i3,i) = w(j1,i2,i3,i)
          end do
       end do
       end do
    end do ! i


!    u(:,:,:) = w(:,:,:,1)**2 + w(:,:,:,2)**2 + w(:,:,:,3)**2

    u(:,:,:) = (w(:,:,:,1)**2+w(:,:,:,2)**2+w(:,:,:,3)**2 ) &
         /( w(:,:,:,0)**(8.0d0/3.0d0) )

    if ( myrank == 0 ) then
       rewind unit
       do i=a1-Md,b1+Md
          write(unit,'(1x,i4,4g20.10)') i,i*Hgrid(1),w(i,0,0,0),w(i,0,0,1) &
               ,u(i,0,0)
       end do
       write(unit,*)
       write(unit,*)
       do i=a2-Md,b2+Md
          write(unit,'(1x,i4,4g20.10)') i,i*Hgrid(2),w(0,i,0,0),w(0,i,0,2) &
               ,u(i,0,0)
       end do
       write(unit,*)
       write(unit,*)
       do i=a3-Md,b3+Md
          write(unit,'(1x,i4,4g20.10)') i,i*Hgrid(3),w(0,0,i,0),w(0,0,i,3) &
               ,u(i,0,0)
       end do
    end if

    call destruct_nabla_fd( nabla )

    deallocate( u )
    deallocate( w )

  END SUBROUTINE func2gp_devr


  SUBROUTINE func2gp_diff( unit, n1, n2, f, g )
    implicit none
    integer,intent(IN) :: unit,n1,n2
    real(8),intent(IN) :: f(n1:n2), g(n1:n2)
    real(8),allocatable :: w(:,:,:)
    integer :: a1,a2,a3,b1,b2,b3,i,i1,i2,i3
    real(8) :: x,y,z,r
    include 'mpif.h'
    a1 = 0
    b1 = Ngrid(1)-1
    a2 = 0
    b2 = Ngrid(2)-1
    a3 = 0
    b3 = Ngrid(3)-1

    allocate( w(a1:b1,a2:b2,a3:b3) ) ; w=0.0d0

    i=n1-1
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       i=i+1
       w(i1,i2,i3) = f(i) - g(i)
    end do
    end do
    end do

    call mpi_allreduce(MPI_IN_PLACE,w,size(w),MPI_REAL8,MPI_SUM,comm_grid,i)

    if ( myrank == 0 ) then
       rewind unit
       do i3=a3,b3
       do i2=a2,b2
       do i1=a1,b1
          x=i1*Hgrid(1)
          y=i2*Hgrid(2)
          z=i3*Hgrid(3)
          r=sqrt(x*x+y*y+z*z)
          write(unit,*) r,w(i1,i2,i3),abs(w(i1,i2,i3))
       end do
       end do
       end do
    end if

    deallocate( w )

  END SUBROUTINE func2gp_diff


END MODULE func2gp_module
