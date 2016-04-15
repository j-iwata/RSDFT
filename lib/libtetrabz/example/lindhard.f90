program dos
  !
  use libtetrabz, only : libtetrabz_polstat, libtetrabz_dos
  implicit none
  !
  integer :: ltetra, nb, ng, nge(3), ngw(3), i1, i2, i3, ik, iq, nq, &
  &          nke, nkw, ikvec(3), ne = 1
  real(8) :: bvec(3,3), ef, kvec(3), qvec(3), qmax, e0(1), VBZ, pi
  real(8),allocatable :: eig(:,:), eig1(:,:), eig2(:,:), wght(:,:,:), wght_dos(:,:,:)
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
  bvec(1:3,1) = (/3d0, 0d0, 0d0/)
  bvec(1:3,2) = (/0d0, 3d0, 0d0/)
  bvec(1:3,3) = (/0d0, 0d0, 3d0/)
  nke = product(nge(1:3))
  nkw = product(ngw(1:3))
  ef = 0.5d0
  nq = 30
  qmax = 4d0
  VBZ = bvec(1,1) * bvec(2,2) * bvec(3,3) + bvec(1,2) * bvec(2,3) * bvec(3,1) &
  &   + bvec(1,3) * bvec(2,1) * bvec(3,2) - bvec(1,3) * bvec(2,2) * bvec(3,1) &
  &   + bvec(1,2) * bvec(2,1) * bvec(3,3) - bvec(1,1) * bvec(2,3) * bvec(3,2)
  !
  allocate(eig1(nb,nke), eig2(nb,nke), wght(nb,nb,nkw), wght_dos(ne,nb,nkw))
  !
  open(10, file = "lindhard.dat")
  !
  do iq = 0, nq
     !
     qvec(1:3) = (/qmax * dble(iq) / dble(nq), 0d0, 0d0/)
     !
     ik = 0
     do i3 = 0, nge(3) - 1
        do i2 = 0, nge(2) - 1
           do i1 = 0, nge(1) - 1
              !
              ik = ik + 1
              !
              ikvec(1:3) = (/i1, i2, i3/)
              ikvec(1:3) = modulo(ikvec(1:3) + nge(1:3) / 2, nge(1:3)) - nge(1:3) / 2
              kvec(1:3) = dble(ikvec(1:3)) / dble(nge(1:3))
              kvec(1:3) = matmul(bvec(1:3,1:3), kvec(1:3))
              !
              eig1(1,ik) = 0.5d0 * dot_product(kvec(1:3), kvec(1:3)) - ef
              !
              kvec(1:3) = kvec(1:3) + qvec(1:3)
              eig2(1,ik) = 0.5d0 * dot_product(kvec(1:3), kvec(1:3)) - ef
              !
           end do
        end do
     end do
     !
     if(iq == 0) then
        e0(1) = 0d0
        call libtetrabz_dos(ltetra,bvec,nb,nge,eig1,ngw,wght_dos,ne,e0)
        write(10,*) 0d0, sum(wght_dos(1:ne,1:nb,1:nkw)) * VBZ / (4d0 * pi)
     else
        call libtetrabz_polstat(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght)
        write(10,*) qvec(1), 2d0 * sum(wght(1:nb,1:nb,1:nkw)) * VBZ / (4d0 * pi)
     end if
     !
  end do ! iq
  !
  close(10)
  !
  write(*,*) "To compare with the analytical result, type below in gnuplot"
  write(*,*) ""
  write(*,*) 'plot "lindhard.dat" u 1:2 w p, 0.5+0.5/x*(1-0.25*x**2)*log(abs((x+2)/(x-2)))'
  write(*,*) ""
  !
end program dos
