MODULE gradient_module

  use grid_module, only: grid
  use bc_module, only: www, bcset
  use fd_module, only: fd,construct_nabla_fd,destruct_nabla_fd
  use lattice_module, only: lattice,get_aa_lattice,get_reciprocal_lattice
  use basic_type_factory

  implicit none

  PRIVATE
  PUBLIC :: gradient, construct_gradient, destruct_gradient &
           ,gradient16, construct_gradient16, destruct_gradient16

  integer,parameter :: DP=kind(0.0d0)
  integer,parameter :: QP=kind(0.0q0)

  type gradient
     real(DP),allocatable :: gx(:),gy(:),gz(:)
     real(DP),allocatable :: gg(:)
  end type gradient

  type gradient16
     real(QP),allocatable :: gx(:),gy(:),gz(:)
     real(QP),allocatable :: gg(:)
  end type gradient16

CONTAINS


  SUBROUTINE construct_gradient( rgrid, rho, grad )

    implicit none
    type( grid ),intent(IN) :: rgrid
    type( GSArray ),intent(IN) :: rho
    type( gradient ),intent(OUT) :: grad
    type(fd) :: nabla
    type(lattice) :: aa,bb
    integer :: i,i1,i2,i3,s,m,Md,m0,m1
    real(DP) :: g1,g2,g3,pi2,b(3,3)

! ---

    call construct_nabla_fd( nabla )

    Md = nabla%md

! ---

    call get_aa_lattice( aa )
    call get_reciprocal_lattice( aa, bb )
    pi2      = 2.0d0*acos(-1.0d0)
    b(1:3,1) = aa%Length(1)*bb%LatticeVector(1:3,1)/( pi2*rgrid%spacing(1) )
    b(1:3,2) = aa%Length(2)*bb%LatticeVector(1:3,2)/( pi2*rgrid%spacing(2) )
    b(1:3,3) = aa%Length(3)*bb%LatticeVector(1:3,3)/( pi2*rgrid%spacing(3) )

! ---

    m0 = rgrid%g1%head
    m1 = rgrid%g1%tail

    allocate( grad%gx(m0:m1) ) ; grad%gx=0.0d0
    allocate( grad%gy(m0:m1) ) ; grad%gy=0.0d0
    allocate( grad%gz(m0:m1) ) ; grad%gz=0.0d0
    allocate( grad%gg(m0:m1) ) ; grad%gg=0.0d0

    www(:,:,:,:)=0.0d0
    do s=rho%s_range%head_global,rho%s_range%tail_global
       i=m0-1
       do i3=rgrid%g3%z%head,rgrid%g3%z%tail
       do i2=rgrid%g3%y%head,rgrid%g3%y%tail
       do i1=rgrid%g3%x%head,rgrid%g3%x%tail
          i=i+1
          www(i1,i2,i3,1) = www(i1,i2,i3,1) + rho%val(i,s)
       end do
       end do
       end do
    end do

    call bcset(1,1,Md,0)

    i=m0-1
    do i3=rgrid%g3%z%head,rgrid%g3%z%tail
    do i2=rgrid%g3%y%head,rgrid%g3%y%tail
    do i1=rgrid%g3%x%head,rgrid%g3%x%tail
       g1=0.0d0
       g2=0.0d0
       g3=0.0d0
       do m=1,Md
          g1 = g1 - nabla%coef(m)*( www(i1-m,i2,i3,1) - www(i1+m,i2,i3,1) )
          g2 = g2 - nabla%coef(m)*( www(i1,i2-m,i3,1) - www(i1,i2+m,i3,1) )
          g3 = g3 - nabla%coef(m)*( www(i1,i2,i3-m,1) - www(i1,i2,i3+m,1) )
       end do
       i=i+1
       grad%gx(i) = b(1,1)*g1 + b(1,2)*g2 + b(1,3)*g3
       grad%gy(i) = b(2,1)*g1 + b(2,2)*g2 + b(2,3)*g3
       grad%gz(i) = b(3,1)*g1 + b(3,2)*g2 + b(3,3)*g3
    end do
    end do
    end do

    do i=m0,m1
       grad%gg(i) = grad%gx(i)**2 + grad%gy(i)**2 + grad%gz(i)**2
    end do

    call destruct_nabla_fd( nabla )

  END SUBROUTINE construct_gradient


  SUBROUTINE destruct_gradient( grad )
    implicit none
    type(gradient) :: grad
    deallocate( grad%gg )
    deallocate( grad%gz )
    deallocate( grad%gy )
    deallocate( grad%gx )
  END SUBROUTINE destruct_gradient


  SUBROUTINE construct_gradient16( rgrid, rho, grad )

    implicit none

    type( grid ),intent(IN) :: rgrid
    type( GSArray ),intent(IN)   :: rho
    type( gradient16 ),intent(OUT) :: grad

    type(fd) :: nabla
    type(lattice) :: aa,bb
    integer :: i,i1,i2,i3,s,m,Md,m0,m1
    integer :: a1,a2,a3,b1,b2,b3
    real(QP) :: g1,g2,g3,pi2,b(3,3)
    real(QP),allocatable :: w(:,:,:)

! ---

    call construct_nabla_fd( nabla )

    Md = nabla%md

! ---

    call get_aa_lattice( aa )
    call get_reciprocal_lattice( aa, bb )
    pi2      = 2.0q0*acos(-1.0q0)
    b(1:3,1) = aa%Length(1)*bb%LatticeVector(1:3,1)/( pi2*rgrid%spacing(1) )
    b(1:3,2) = aa%Length(2)*bb%LatticeVector(1:3,2)/( pi2*rgrid%spacing(2) )
    b(1:3,3) = aa%Length(3)*bb%LatticeVector(1:3,3)/( pi2*rgrid%spacing(3) )

! ---

    m0 = rgrid%g1%head
    m1 = rgrid%g1%tail
    allocate( grad%gx(m0:m1) ) ; grad%gx=0.0q0
    allocate( grad%gy(m0:m1) ) ; grad%gy=0.0q0
    allocate( grad%gz(m0:m1) ) ; grad%gz=0.0q0
    allocate( grad%gg(m0:m1) ) ; grad%gg=0.0q0

    www(:,:,:,:)=0.0d0
    do s=rho%s_range%head_global,rho%s_range%tail_global
       i=m0-1
       do i3=rgrid%g3%z%head,rgrid%g3%z%tail
       do i2=rgrid%g3%y%head,rgrid%g3%y%tail
       do i1=rgrid%g3%x%head,rgrid%g3%x%tail
          i=i+1
          www(i1,i2,i3,1) = www(i1,i2,i3,1) + rho%val(i,s)
       end do
       end do
       end do
    end do

    call bcset(1,1,Md,0)

    a1 = rgrid%g3%x%head - Md
    b1 = rgrid%g3%x%tail + Md
    a2 = rgrid%g3%y%head - Md
    b2 = rgrid%g3%y%tail + Md
    a3 = rgrid%g3%z%head - Md
    b3 = rgrid%g3%z%tail + Md
    allocate( w(a1:b1,a2:b2,a3:b3) ) ; w=0.0q0

    w=www(:,:,:,1)

    i=m0-1
    do i3=rgrid%g3%z%head,rgrid%g3%z%tail
    do i2=rgrid%g3%y%head,rgrid%g3%y%tail
    do i1=rgrid%g3%x%head,rgrid%g3%x%tail
       g1=0.0q0
       g2=0.0q0
       g3=0.0q0
       do m=1,Md
          g1 = g1 - nabla%coef(m)*( w(i1-m,i2,i3) - w(i1+m,i2,i3) )
          g2 = g2 - nabla%coef(m)*( w(i1,i2-m,i3) - w(i1,i2+m,i3) )
          g3 = g3 - nabla%coef(m)*( w(i1,i2,i3-m) - w(i1,i2,i3+m) )
       end do
       i=i+1
       grad%gx(i) = b(1,1)*g1 + b(1,2)*g2 + b(1,3)*g3
       grad%gy(i) = b(2,1)*g1 + b(2,2)*g2 + b(2,3)*g3
       grad%gz(i) = b(3,1)*g1 + b(3,2)*g2 + b(3,3)*g3
    end do
    end do
    end do

    do i=m0,m1
       grad%gg(i) = grad%gx(i)**2 + grad%gy(i)**2 + grad%gz(i)**2
    end do

    call destruct_nabla_fd( nabla )

  END SUBROUTINE construct_gradient16


  SUBROUTINE destruct_gradient16( grad )
    implicit none
    type(gradient16) :: grad
    deallocate( grad%gg )
    deallocate( grad%gz )
    deallocate( grad%gy )
    deallocate( grad%gx )
  END SUBROUTINE destruct_gradient16


END MODULE gradient_module
