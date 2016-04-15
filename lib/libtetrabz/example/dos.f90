program dos
  !
  use libtetrabz, only : libtetrabz_fermieng, libtetrabz_dos
  implicit none
  !
  integer :: ltetra, nb, ng, nge(3), ngw(3), i1, i2, i3, ik, ne, ie, nke, nkw
  real(8) :: bvec(3,3), ef, nelec, kvec(3), pi
  real(8),allocatable :: eig(:,:), wght(:,:), e0(:), wght_dos(:,:,:)
  !
  write(*,'(a)', advance = "no") "Which tetrahedron method ?(1 = Linear, 2 = Optimized): "
  read(*,*) ltetra
  write(*,'(a)', advance = "no") "k-point mesh ?: "
  read(*,*) ng
  !
  pi = acos(-1d0)
  nb = 1
  nge(1:3) = ng
  ngw(1:3) = ng
  nelec = 0.5
  bvec(1:3,1) = (/1d0, 0d0, 0d0/)
  bvec(1:3,2) = (/0d0, 1d0, 0d0/)
  bvec(1:3,3) = (/0d0, 0d0, 1d0/)
  nke = product(nge(1:3))
  nkw = product(ngw(1:3))
  !
  allocate(eig(nb,nke), wght(nb,nkw))
  !
  ik = 0
  do i3 = 0, nge(3) - 1
     do i2 = 0, nge(2) - 1
        do i1 = 0, nge(1) - 1
           !
           ik = ik + 1
           kvec(1:3) = 2d0 * pi * dble((/i1, i2, i3/) - nge(1:3) / 2) / dble(nge(1:3))
           !
           eig(1,ik) = - sum(cos(kvec(1:3)))
           !
        end do
     end do
  end do
  !
  call libtetrabz_fermieng(ltetra,bvec,nb,nge,eig,ngw,wght,ef,nelec)
  !
  write(*,*) "  E_F = ", ef
  !
  ne = 100
  allocate(wght_dos(ne,nb,nkw), e0(ne))
  !
  do ie = 1, ne
     e0(ie) = 6d0 / dble(ne - 1) * dble(ie - 1) - 3d0
  end do
  !
  call libtetrabz_dos(ltetra,bvec,nb,nge,eig,ngw,wght_dos,ne,e0)
  !
  open(10, file = "dos.dat")
  do ie = 1, ne
     write(10,*) e0(ie), sum(wght_dos(ie,1:nb,1:nkw))
  end do
  !
  close(10)
  !
end program dos
