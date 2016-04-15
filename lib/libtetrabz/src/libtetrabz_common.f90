module libtetrabz_vals
  !
  integer,save :: &
  & fst,          & !
  & lst,          & !
  & nk,           & !
  & nk0,          & !
  & nt,           & !
  & nb,           & !
  & ne,           & !
  & ng(3),        & !
  & ltetra,       & !
  & ivvec(3,20,6)   !
  !
  real(8),save :: &
  & wlsm(4,20)      !
  !
  integer,save,allocatable :: &
  & indx1(:,:), &
  & indx2(:,:), &
  & indx3(:)
  !
end module libtetrabz_vals
!
!
!
module libtetrabz_routines
  !
  implicit none
  !
contains
!
! define shortest diagonal line & define type of tetragonal
!
subroutine libtetrabz_initialize(bvec)
  !
  use libtetrabz_vals, only : ltetra, ng, wlsm, ivvec, nt, nk, ng
  !
  real(8),intent(in) :: bvec(3,3)
  !
  integer :: itype, i1, i2, i3, it, ii, divvec(4,4), ivvec0(4)
  real(8) :: l(4), bvec2(3,3), bvec3(3,4)
  !
  nk = product(ng(1:3))
  nt = nk * 6
  !
  do i1 = 1, 3
     bvec2(1:3,i1) = bvec(1:3,i1) / dble(ng(i1))
  end do
  !
  bvec3(1:3,1) = -bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
  bvec3(1:3,2) =  bvec2(1:3,1) - bvec2(1:3,2) + bvec2(1:3,3)
  bvec3(1:3,3) =  bvec2(1:3,1) + bvec2(1:3,2) - bvec2(1:3,3)
  bvec3(1:3,4) =  bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
  !
  ! length of delta bvec
  !
  do i1 = 1, 4
     l(i1) = dot_product(bvec3(1:3,i1),bvec3(1:3,i1))
  end do
  !
  itype = minloc(l(1:4),1)
  !
  ! start & last
  !
  ivvec0(1:4) = (/ 0, 0, 0, 0 /)
  !
  divvec(1:4,1) = (/ 1, 0, 0, 0 /)
  divvec(1:4,2) = (/ 0, 1, 0, 0 /)
  divvec(1:4,3) = (/ 0, 0, 1, 0 /)
  divvec(1:4,4) = (/ 0, 0, 0, 1 /)
  !
  ivvec0(itype) = 1
  divvec(itype, itype) = - 1
  !
  ! Corners of tetrahedra
  !
  it = 0
  do i1 = 1, 3
     do i2 = 1, 3
        if(i2 == i1) cycle
        do i3 = 1, 3
           if(i3 == i1 .or. i3 == i2) cycle
           !
           it = it + 1
           !
           ivvec(1:3,1,it) = ivvec0(1:3)
           ivvec(1:3,2,it) = ivvec(1:3,1,it) + divvec(1:3,i1)
           ivvec(1:3,3,it) = ivvec(1:3,2,it) + divvec(1:3,i2)
           ivvec(1:3,4,it) = ivvec(1:3,3,it) + divvec(1:3,i3)
           !
        end do
     end do
  end do
  !
  ! Additional points
  !
  ivvec(1:3, 5,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,2,1:6)
  ivvec(1:3, 6,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,3,1:6)
  ivvec(1:3, 7,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,4,1:6)
  ivvec(1:3, 8,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,1,1:6)
  !
  ivvec(1:3, 9,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,3,1:6)
  ivvec(1:3,10,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,4,1:6)
  ivvec(1:3,11,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,1,1:6)
  ivvec(1:3,12,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,2,1:6)
  !
  ivvec(1:3,13,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,4,1:6)
  ivvec(1:3,14,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,1,1:6)
  ivvec(1:3,15,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,2,1:6)
  ivvec(1:3,16,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,3,1:6)
  !
  ivvec(1:3,17,1:6) =  ivvec(1:3,4,1:6) - ivvec(1:3,1,1:6) + ivvec(1:3,2,1:6)
  ivvec(1:3,18,1:6) =  ivvec(1:3,1,1:6) - ivvec(1:3,2,1:6) + ivvec(1:3,3,1:6)
  ivvec(1:3,19,1:6) =  ivvec(1:3,2,1:6) - ivvec(1:3,3,1:6) + ivvec(1:3,4,1:6)
  ivvec(1:3,20,1:6) =  ivvec(1:3,3,1:6) - ivvec(1:3,4,1:6) + ivvec(1:3,1,1:6)
  !
  if(ltetra == 1) then
     !
     !write(*,*) "[libtetrabz] Linear tetrahedron method is used."
     !
     wlsm(1:4,1:20) = 0.0d0
     wlsm(1,1) = 1.0d0
     wlsm(2,2) = 1.0d0
     wlsm(3,3) = 1.0d0
     wlsm(4,4) = 1.0d0
     !
  else if(ltetra == 2) then
     !
     !write(*,*) "[libtetrabz] Improved tetrahedron method is used."
     !
     wlsm(1, 1: 4) = dble((/1440,    0,   30,    0/))
     wlsm(2, 1: 4) = dble((/   0, 1440,    0,   30/))
     wlsm(3, 1: 4) = dble((/  30,    0, 1440,    0/))
     wlsm(4, 1: 4) = dble((/   0,   30,    0, 1440/))
     !
     wlsm(1, 5: 8) = dble((/ -38,    7,   17,  -28/))
     wlsm(2, 5: 8) = dble((/ -28,  -38,    7,   17/))
     wlsm(3, 5: 8) = dble((/  17,  -28,  -38,    7/))
     wlsm(4, 5: 8) = dble((/   7,   17,  -28,  -38/))
     !
     wlsm(1, 9:12) = dble((/ -56,    9,  -46,    9/))
     wlsm(2, 9:12) = dble((/   9,  -56,    9,  -46/))
     wlsm(3, 9:12) = dble((/ -46,    9,  -56,    9/))
     wlsm(4, 9:12) = dble((/   9,  -46,    9,  -56/))
     !
     wlsm(1,13:16) = dble((/ -38,  -28,   17,    7/))
     wlsm(2,13:16) = dble((/   7,  -38,  -28,   17/))
     wlsm(3,13:16) = dble((/  17,    7,  -38,  -28/))
     wlsm(4,13:16) = dble((/ -28,   17,    7,  -38/))
     !
     wlsm(1,17:20) = dble((/ -18,  -18,   12,  -18/))
     wlsm(2,17:20) = dble((/ -18,  -18,  -18,   12/))
     wlsm(3,17:20) = dble((/  12,  -18,  -18,  -18/))
     wlsm(4,17:20) = dble((/ -18,   12,  -18,  -18/))
     !
     wlsm(1:4,1:20) = wlsm(1:4,1:20) / 1260d0
     !
  else
     !
     write(*,*) "[libtetrabz] STOP! ltetrta is invalid."
     stop
     !
  end if
  !
end subroutine libtetrabz_initialize
!
! Main subroutine for occupation : Theta(EF - E1)
!
subroutine libtetrabz_occ1(ef,eig,occ)
  !
  use libtetrabz_vals, only : nb, nk, nk0, fst, lst, wlsm, indx1, indx2
  !
  real(8),intent(in) :: ef, eig(nb,nk)
  real(8),intent(out) :: occ(nb,nk0)
  !
  integer :: it, ib, ii
  real(8) :: e(4), a(4,4), V, tmp(5,4), &
  &          w0(4,4), w1(4), ei(4,nb)
  !
  w0(1:4,1:4) = 0d0
  do ii = 1, 4
     w0(ii,ii) = 1d0
  end do
  !
  occ(1:nb,1:nk0) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(fst,lst,nb,wlsm,indx1,indx2,eig,ef,w0,occ) &
  !$OMP & PRIVATE(it,ii,ib,ei,tmp,a,e,V,w1)
  !
  do it = fst, lst
     !
     ei(1:4,1:nb) = 0d0
     do ii = 1, 20
        do ib = 1, nb
           ei(1:4,ib) = ei(1:4,ib) + wlsm(1:4,ii) * eig(ib,indx1(ii,it))
        end do
     end do
     !
     !$OMP DO
     do ib = 1, nb
        !
        w1(1:4) = 0d0
        !
        tmp(  1,1:4) = ei(1:4,  ib)
        tmp(2:5,1:4) = w0(1:4, 1:4)
        call libtetrabz_sort(5,4,tmp)
        !
        e(1:4) = tmp(1,1:4)
        !
        do ii = 1, 4
           a(1:4,ii) = (ef - e(ii)) / (e(1:4) - e(ii))
        end do
        !
        if(e(1) <= ef .and. ef < e(2)) then
           !
           ! A - 1
           !
           V = 0.25d0 * a(2,1) * a(3,1) * a(4,1)
           !
           w1(1:4) = V * ( tmp(2:5,1)                                &
           &             + tmp(2:5,1) * a(1,2) + tmp(2:5,2) * a(2,1) &
           &             + tmp(2:5,1) * a(1,3) + tmp(2:5,3) * a(3,1) &
           &             + tmp(2:5,1) * a(1,4) + tmp(2:5,4) * a(4,1) )
           !
        else if(e(2) <= ef .and. ef < e(3)) then
           !
           ! B - 1
           !
           V = 0.25d0 * a(3,1) * a(4,1) * a(2,4)
           !
           w1(1:4) = V * ( tmp(2:5,1)                                &
           &             + tmp(2:5,1) * a(1,3) + tmp(2:5,3) * a(3,1) &
           &             + tmp(2:5,1) * a(1,4) + tmp(2:5,4) * a(4,1) &
           &             + tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2) )
           !
           ! B - 2
           !
           V = 0.25d0 * a(3,2) * a(4,2)
           !
           w1(1:4) = w1(1:4) &
           &       + V * ( tmp(2:5,1)                                & 
           &             + tmp(2:5,2)                                &
           &             + tmp(2:5,2) * a(2,3) + tmp(2:5,3) * a(3,2) &
           &             + tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2) )
           !
           ! B - 3
           !
           V = 0.25d0 * a(2,3) * a(3,1) * a(4,2)
           !
           w1(1:4) = w1(1:4)                                         &
           &       + V * ( tmp(2:5,1)                                &
           &             + tmp(2:5,1) * a(1,3) + tmp(2:5,3) * a(3,1) &
           &             + tmp(2:5,2) * a(2,3) + tmp(2:5,3) * a(3,2) &
           &             + tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2) )
           !
        else if(e(3) <= ef .and. ef < e(4)) then
           !
           ! C - 1
           !
           V = 0.25d0 * a(4,3)
           !
           w1(1:4) = V * ( tmp(2:5,1)                                &
           &             + tmp(2:5,2)                                &
           &             + tmp(2:5,3)                                &
           &             + tmp(2:5,3) * a(3,4) + tmp(2:5,4) * a(4,3) )
           !
           ! C - 2
           !
           V = 0.25d0 * a(3,4) * a(4,2)
           !
           w1(1:4) = w1(1:4)                                         &
           &       + V * ( tmp(2:5,1)                                &
           &             + tmp(2:5,2)                                &
           &             + tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2) &
           &             + tmp(2:5,3) * a(3,4) + tmp(2:5,4) * a(4,3) )
           !
           ! C - 3
           !
           V = 0.25d0 * a(3,4) * a(2,4) * a(4,1)
           !
           w1(1:4) = w1(1:4)                                         &
           &       + V * ( tmp(2:5,1)                                &
           &             + tmp(2:5,1) * a(1,4) + tmp(2:5,4) * a(4,1) &
           &             + tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2) &
           &             + tmp(2:5,3) * a(3,4) + tmp(2:5,4) * a(4,3) )
           !
        else if(e(4) <= ef) then
           !
           ! D - 1
           !
           V = 0.25d0
           !
           w1(1:4) = V * ( tmp(2:5,1) &
           &             + tmp(2:5,2) &
           &             + tmp(2:5,3) &
           &             + tmp(2:5,4) )
           !
        end if
        !
        do ii = 1, 20
           occ(ib,indx2(ii,it)) = occ(ib,indx2(ii,it)) + dot_product(wlsm(1:4,ii), w1(1:4))
        end do ! ii
        !
     end do ! ib
     !$OMP END DO NOWAIT
     !
  end do ! it
  !
  !$OMP END PARALLEL
  !
  occ(1:nb,1:nk0) = occ(1:nb,1:nk0) / dble(6 * nk)
  !
end subroutine libtetrabz_occ1
!
! Simple sort
!
subroutine libtetrabz_sort(n1,n2,a)
  !
  !
  integer,intent(in) :: n1, n2
  real(8),intent(inout) :: a(n1,n2) 
  !
  integer :: i, m
  real(8) :: am, atmp(n1)
  !
  do i = 1, n2 - 1
     am = minval(a(1,i+1:n2) )
     m  = minloc(a(1,i+1:n2),1) + i
     if(a(1,i) .gt. am) then
        atmp(1:n1) = a(1:n1, m)
        a(1:n1, m) = a(1:n1,i)
        a(1:n1, i) = atmp(1:n1)
     end if
  end do
  !
end subroutine libtetrabz_sort
!
! Main subroutine for Dos : Delta(E - E1)
!
subroutine libtetrabz_dos1(eig,e0,dos)
  !
  use libtetrabz_vals, only : nb, nk, nk0, ne, fst, lst, wlsm, indx1, indx2
  !
  real(8),intent(in) :: eig(nb,nk), e0(ne)
  real(8),intent(out) :: dos(ne,nb,nk0)
  !
  integer :: ib, it, ii, ie
  real(8) :: e(4), a(4,4), tmp(5,4), V, &
  &          w0(4,4), w1(4,ne), ei(4,nb)
  !
  w0(1:4,1:4) = 0d0
  do ii = 1, 4
     w0(ii,ii) = 1d0
  end do
  !
  dos(1:ne, 1:nb, 1:nk0) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(fst,lst,nb,ne,indx1,indx2,wlsm,eig,w0,e0,dos) &
  !$OMP & PRIVATE(ib,it,ii,ie,e,a,tmp,w1,ei,V)
  !
  do it = fst, lst
     !
     ei(1:4,1:nb) = 0d0
     do ii = 1, 20
        do ib = 1, nb
           ei(1:4,ib) = ei(1:4,ib) + wlsm(1:4,ii) * eig(ib,indx1(ii,it))
        end do
     end do
     !
     !$OMP DO
     do ib = 1, nb
        !
        tmp(  1,1:4) = ei(1:4,ib)
        tmp(2:5,1:4) = w0(1:4,1:4)
        call libtetrabz_sort(5, 4, tmp)
        e(1:4) = tmp(1,1:4)
        !
        w1(1:4,1:ne) = 0d0
        !
        do ie = 1, ne
           !
           do ii = 1, 4
              a(1:4,ii) = (e0(ie) - e(ii)) / (e(1:4) - e(ii))
           end do
           !
           if(e(1) < e0(ie) .and. e0(ie) <= e(2)) then
              !
              ! A
              !
              V = a(2,1) * a(3,1) * a(4,1) / (e0(ie) - e(1))
              !
              w1(1:4,ie) = V * ( tmp(2:5,1) * a(1,2) + tmp(2:5,2) * a(2,1) &
              &                + tmp(2:5,1) * a(1,3) + tmp(2:5,3) * a(3,1) &
              &                + tmp(2:5,1) * a(1,4) + tmp(2:5,4) * a(4,1) )
              !
           else if(e(2) < e0(ie) .and. e0(ie) <= e(3)) then
              !
              ! B - 1
              !
              V = a(3,1) * a(4,1) * a(2,4) / (e0(ie) - e(1))
              !
              w1(1:4,ie) = V * ( tmp(2:5,1) * a(1,3) + tmp(2:5,3) * a(3,1) &
              &                + tmp(2:5,1) * a(1,4) + tmp(2:5,4) * a(4,1) &
              &                + tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2) )
              !
              ! B - 2
              !
              V = a(2,3) * a(3,1) * a(4,2) / (e0(ie) - e(1))
              !
              w1(1:4,ie) = w1(1:4,ie)                                      &
              &          + V * ( tmp(2:5,1) * a(1,3) + tmp(2:5,3) * a(3,1) &
              &                + tmp(2:5,2) * a(2,3) + tmp(2:5,3) * a(3,2) &
              &                + tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2) )
              !
           else if(e(3) < e0(ie) .and. e0(ie) < e(4)) then
              !
              ! C
              !
              V = a(1,4) * a(2,4) * a(3,4) / (e(4) - e0(ie))
              !
              w1(1:4,ie) = V * ( tmp(2:5,1) * a(1,4) + tmp(2:5,4) * a(4,1) &
              &                + tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2) &
              &                + tmp(2:5,3) * a(3,4) + tmp(2:5,4) * a(4,3) )
              !
           end if
           !
        end do ! ie
        !
        do ii = 1, 20
           dos(1:ne,ib,indx2(ii,it)) = dos(1:ne,ib,indx2(ii,it)) + &
           &   matmul(wlsm(1:4,ii), w1(1:4,1:ne))
        end do ! ii
        !
     end do ! ib
     !$OMP END DO NOWAIT
     !
  end do ! it
  !
  !$OMP END PARALLEL
  !
  dos(1:ne,1:nb,1:nk0) = dos(1:ne,1:nb,1:nk0) / dble(6 * nk)
  !
end subroutine libtetrabz_dos1
!
! Main subroutine for Delta(E1) * Delta(E2)
!
subroutine libtetrabz_doubledelta1(eig1,eig2,ddel)
  !
  use libtetrabz_vals, only : nb, nk, nk0, fst, lst, indx1, indx2, wlsm
  !
  real(8),intent(in) :: eig1(nb,nk), eig2(nb,nk)
  real(8),intent(out) :: ddel(nb,nb,nk0)
  !
  integer :: it, ib, ii, nn
  real(8) :: e(4), a(4,4), V, &
  &          ei(4,nb), ej(nb,4), ej2(nb,4), &
  &          w0(4,nb,4), w1(4,nb), w2(4,nb,3), &
  &          tmp(1 + nb + 4 * nb, 4), tmp2(nb + 4 * nb, 3)
  !
  nn = 1 + nb + 4 * nb
  !
  w0(1:4,1:nb,1:4) = 0d0
  do ii = 1, 4
     w0(ii,1:nb,ii) = 1d0
  end do
  !
  ddel(1:nb,1:nb,1:nk0) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(fst,lst,nb,nn,indx1,indx2,wlsm,eig1,eig2,w0,ddel) &
  !$OMP & PRIVATE(ib,it,ii,e,a,tmp,tmp2,w1,w2,ei,ej,ej2,V)
  !
  do it = fst, lst
     !
     ei(1:4, 1:nb) = 0d0
     ej(1:nb,1:4) = 0d0
     do ii = 1, 20
        do ib = 1, nb
           ei(1:4,ib) = ei(1:4,ib) + wlsm(1:4,ii) * eig1(ib, indx1(ii,it))
           ej(ib,1:4) = ej(ib,1:4) + wlsm(1:4,ii) * eig2(ib, indx1(ii,it))
        end do
     end do
     !
     !$OMP DO
     do ib = 1, nb
        !
        w1(1:4,1:nb) = 0d0
        !
        tmp(1,                      1:4) = ei(1:4,ib)
        tmp(2:1 + nb,               1:4) = ej(1:nb,1:4)
        tmp(2 + nb:1 + nb + 4 * nb, 1:4) = reshape(w0(1:4,1:nb,1:4), (/4 * nb, 4/))
        !
        call libtetrabz_sort(nn, 4, tmp)
        !
        e(1:4) = tmp(1, 1:4)
        !
        do ii = 1, 4
           a(1:4,ii) = (0d0 - e(ii)) / (e(1:4) - e(ii))
        end do
        !
        if(e(1) < 0d0 .and. 0d0 <= e(2)) then
           !
           ! A
           !
           !V = 3d0 * a(2,1) * a(3,1) * a(4,1) / (0d0 - e(1))
           V = 3d0 * a(2,1) * a(3,1)           / (e(4) - e(1))
           !
           tmp2(1:nn - 1,1) = tmp(2:nn,1) * a(1,2) + tmp(2:nn,2) * a(2,1)
           tmp2(1:nn - 1,2) = tmp(2:nn,1) * a(1,3) + tmp(2:nn,3) * a(3,1) 
           tmp2(1:nn - 1,3) = tmp(2:nn,1) * a(1,4) + tmp(2:nn,4) * a(4,1)
           !
           ej2(1:nb,1:3) = tmp2(1:nb,1:3)
           w2(1:4,1:nb,1:3) = reshape(tmp2(nb + 1:nb + 4 * nb,1:3), (/4, nb, 3/))
           !
           call libtetrabz_doubledelta2(ej2,w2)
           !
           w1(1:4,1:nb) = w1(1:4,1:nb) + V * sum(w2(1:4,1:nb,1:3), 3)
           !
        else if( e(2) < 0d0 .and. 0d0 <= e(3)) then
           !
           ! B - 1
           !
           !V = 3d0 * a(3,1) * a(4,1) * a(2,4) / (0d0 - e(1))
           V = 3d0           * a(4,1) * a(2,4) / (e(3) - e(1))
           !
           tmp2(1:nn - 1,1) = tmp(2:nn,1) * a(1,3) + tmp(2:nn,3) * a(3,1) 
           tmp2(1:nn - 1,2) = tmp(2:nn,1) * a(1,4) + tmp(2:nn,4) * a(4,1) 
           tmp2(1:nn - 1,3) = tmp(2:nn,2) * a(2,4) + tmp(2:nn,4) * a(4,2) 
           !
           ej2(1:nb,1:3) = tmp2(1:nb,1:3)
           w2(1:4,1:nb,1:3) = reshape(tmp2(nb + 1:nb + 4 * nb,1:3), (/4, nb, 3/))
           !
           call libtetrabz_doubledelta2(ej2,w2)
           !
           w1(1:4,1:nb) = w1(1:4,1:nb) + V * sum(w2(1:4,1:nb,1:3), 3)
           !
           ! B - 2
           !
           !V = 3d0 * a(2,3) * a(3,1) * a(4,2) / (0d0 - e(1))
           V = 3d0 * a(2,3)           * a(4,2) / (e(3) - e(1))
           !
           tmp2(1:nn - 1,1) = tmp(2:nn,1) * a(1,3) + tmp(2:nn,3) * a(3,1) 
           tmp2(1:nn - 1,2) = tmp(2:nn,2) * a(2,3) + tmp(2:nn,3) * a(3,2) 
           tmp2(1:nn - 1,3) = tmp(2:nn,2) * a(2,4) + tmp(2:nn,4) * a(4,2) 
           !
           ej2(1:nb,1:3) = tmp2(1:nb,1:3)
           w2(1:4,1:nb,1:3) = reshape(tmp2(nb + 1:nb + 4 * nb,1:3), (/4, nb, 3/))
           !
           call libtetrabz_doubledelta2(ej2,w2)
           !
           w1(1:4,1:nb) = w1(1:4,1:nb) + V * sum(w2(1:4,1:nb,1:3), 3)
           !
        else if(e(3) < 0d0 .and. 0d0 < e(4)) then
           !
           ! C
           !
           !V = 3d0 * a(1,4) * a(2,4) * a(3,4) / (e(4) - 0d0)
           V = 3d0 * a(1,4) * a(2,4)           / (e(4) - e(3))
           !
           tmp2(1:nn - 1,1) = tmp(2:nn,1) * a(1,4) + tmp(2:nn,4) * a(4,1) 
           tmp2(1:nn - 1,2) = tmp(2:nn,2) * a(2,4) + tmp(2:nn,4) * a(4,2) 
           tmp2(1:nn - 1,3) = tmp(2:nn,3) * a(3,4) + tmp(2:nn,4) * a(4,3) 
           !
           ej2(1:nb,1:3) = tmp2(1:nb,1:3)
           w2(1:4,1:nb,1:3) = reshape(tmp2(nb + 1:nb + 4 * nb,1:3), (/4, nb, 3/))
           ! 
           call libtetrabz_doubledelta2(ej2,w2)
           !
           w1(1:4,1:nb) = w1(1:4,1:nb) + V * sum(w2(1:4,1:nb,1:3), 3)
           !
        end if
        !
        do ii = 1, 20
           ddel(1:nb,ib,indx2(ii,it)) = ddel(1:nb,ib,indx2(ii,it)) &
           &                         + matmul(wlsm(1:4,ii), w1(1:4,1:nb))
        end do ! ii
        !
     end do ! ib
     !$OMP END DO NOWAIT
     !
  end do ! it
  !
  !$OMP END PARALLEL
  !
  ddel(1:nb,1:nb,1:nk0) = ddel(1:nb,1:nb,1:nk0) / dble(6 * nk)
  !
end subroutine libtetrabz_doubledelta1
!
! 2nd step of tetrahedra method.
!
subroutine libtetrabz_doubledelta2(ej,w)
  !
  use libtetrabz_vals, only : nb
  !
  real(8),intent(in) :: ej(nb,3)
  real(8),intent(inout) :: w(4,nb,3)
  !
  integer :: ib, ii
  real(8) :: tmp(5,3), e(3), a(3,3), V
  !
  do ib = 1, nb
     !
     if(maxval(abs(ej(ib,1:3))) < 1d-10) stop "Nesting !!"
     !
     tmp(  1, 1:3) = ej(     ib, 1:3)
     tmp(2:5, 1:3) = w(1:4, ib, 1:3)
     !
     call libtetrabz_sort(5, 3, tmp)
     !
     e(1:3) = tmp(1,1:3)
     w(1:4, ib, 1:3) = 0d0
     !
     do ii = 1, 3
        a(1:3,ii) = (0d0 - e(ii)) / (e(1:3) - e(ii))
     end do
     !
     if((e(1) < 0d0 .and. 0d0 <= e(2)) .or. (e(1) <= 0d0 .and. 0d0 < e(2))) then
        !
        !V = a(2,1) * a(3,1) / (0d0 - e(1)) 
        V = a(2,1)           / (e(3) - e(1)) 
        !
        w(1:4,ib,1) = tmp(2:5,1) * V * (a(1,2) + a(1,3))
        w(1:4,ib,2) = tmp(2:5,2) * V * a(2,1)
        w(1:4,ib,3) = tmp(2:5,3) * V * a(3,1)
        !
     else if((e(2) <= 0d0 .and. 0d0 < e(3)) .or. (e(2) < 0d0 .and. 0d0 <= e(3))) then
        !
        !V = a(1,3) * a(2,3) / (e(3) - 0d0) 
        V = a(1,3)           / (e(3) - e(2)) 
        !
        w(1:4,ib,1) = tmp(2:5,1) * V * a(1,3)
        w(1:4,ib,2) = tmp(2:5,2) * V * a(2,3)
        w(1:4,ib,3) = tmp(2:5,3) * V * (a(3,1) + a(3,2))
        !
     end if
     !
  end do ! ib
  !
end subroutine libtetrabz_doubledelta2
!
! Main subroutine for Theta(- E1) * Theta(E1 - E2)
!
subroutine libtetrabz_occstep1(eig1,eig2,ocst)
  !
  use libtetrabz_vals, ONLY : nb, nk, nk0, fst, lst, indx1, indx2, wlsm
  !
  real(8),intent(in) :: eig1(nb,nk), eig2(nb,nk)
  real(8),intent(out) :: ocst(nb,nb,nk)
  !
  integer :: it, ib, ii, nn
  real(8) :: e(4), a(4,4), V, thr = 1d-10, &
  &          ei(4,nb), ej(nb,4), de(nb,4), &
  &          w0(4,nb,4), w1(4,nb), w2(4,nb,4), &
  &          tmp(1 + nb + 4 * nb,4), tmp2(nb + 4 * nb, 4)
  !
  nn = 1 + nb + 4 * nb
  !
  w0(1:4,1:nb,1:4) = 0d0
  do ii = 1, 4
     w0(ii,1:nb,ii) = 1d0
  end do
  !
  ocst(1:nb,1:nb,1:nk0) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(fst,lst,nb,nn,indx1,indx2,wlsm,eig1,eig2,w0,thr,ocst) &
  !$OMP & PRIVATE(ib,it,ii,e,a,tmp,tmp2,w1,w2,ei,ej,de,V)
  !
  do it = fst, lst
     !
     ei(1:4, 1:nb) = 0d0
     ej(1:nb, 1:4) = 0d0
     do ii = 1, 20
        do ib = 1, nb
           ei(1:4, ib) = ei(1:4, ib) + wlsm(1:4,ii) * eig1(ib, indx1(ii,it))
           ej(ib, 1:4) = ej(ib, 1:4) + wlsm(1:4,ii) * eig2(ib, indx1(ii,it))
        end do
     end do
     !
     !$OMP DO
     do ib = 1, nb
        !
        w1(1:4,1:nb) = 0d0
        !
        tmp(1,                      1:4) = ei(1:4,ib)
        do ii = 1, 4
           tmp(2:1 + nb, ii) = ej(1:nb,ii) - ei(ii,ib)
        end do
        tmp(2 + nb:1 + nb + 4 * nb, 1:4) = reshape(w0(1:4,1:nb,1:4), (/4 * nb, 4/))
        !
        call libtetrabz_sort(nn, 4, tmp)
        !
        e(1:4) = tmp(1, 1:4)
        !
        do ii = 1, 4
           a(1:4,ii) = (0d0 - e(ii) ) / (e(1:4) - e(ii) )
        end do
        !
        if( e(1) <= 0d0 .and. 0d0 < e(2) ) then
           !
           ! A - 1
           !
           V = a(2,1) * a(3,1) * a(4,1)
           !
           if(V > thr) then
              !
              tmp2(1:nn - 1,1) = tmp(2:nn,1)
              tmp2(1:nn - 1,2) = tmp(2:nn,1) * a(1,2) + tmp(2:nn,2) * a(2,1) 
              tmp2(1:nn - 1,3) = tmp(2:nn,1) * a(1,3) + tmp(2:nn,3) * a(3,1)
              tmp2(1:nn - 1,4) = tmp(2:nn,1) * a(1,4) + tmp(2:nn,4) * a(4,1)
              !
              de(1:nb,1:4) = tmp2(1:nb, 1:4)
              w2(1:4, 1:nb, 1:4) = reshape(tmp2(nb + 1:nb + 4 * nb, 1:4), (/4, nb, 4/))
              !
              call libtetrabz_occstep2(de,w2)
              !
              w1(1:4,1:nb) = w1(1:4,1:nb) + V * sum(w2(1:4,1:nb,1:4), 3)
              !
           end if
           !
        else if( e(2) <= 0d0 .and. 0d0 < e(3)) then
           !
           ! B - 1
           !
           V = a(3,1) * a(4,1) * a(2,4)
           !
           if(V > thr) then
              !
              tmp2(1:nn - 1,1) = tmp(2:nn,1)
              tmp2(1:nn - 1,2) = tmp(2:nn,1) * a(1,3) + tmp(2:nn,3) * a(3,1) 
              tmp2(1:nn - 1,3) = tmp(2:nn,1) * a(1,4) + tmp(2:nn,4) * a(4,1) 
              tmp2(1:nn - 1,4) = tmp(2:nn,2) * a(2,4) + tmp(2:nn,4) * a(4,2) 
              !
              de(1:nb,1:4) = tmp2(1:nb, 1:4)
              w2(1:4, 1:nb, 1:4) = reshape(tmp2(nb + 1:nb + 4 * nb, 1:4), (/4, nb, 4/))
              !
              call libtetrabz_occstep2(de,w2)
              !
              w1(1:4,1:nb) = w1(1:4,1:nb) + V * sum(w2(1:4,1:nb,1:4), 3)
              !
           end if
           !
           ! B - 2
           !
           V = a(3,2) * a(4,2)
           !
           if(V > thr) then
              !
              tmp2(1:nn - 1,1:2) = tmp(2:nn,1:2)
              tmp2(1:nn - 1,3)   = tmp(2:nn,2) * a(2,3) + tmp(2:nn,3) * a(3,2) 
              tmp2(1:nn - 1,4)   = tmp(2:nn,2) * a(2,4) + tmp(2:nn,4) * a(4,2) 
              !
              de(1:nb,1:4) = tmp2(1:nb, 1:4)
              w2(1:4, 1:nb, 1:4) = reshape(tmp2(nb + 1:nb + 4 * nb, 1:4), (/4, nb, 4/))
              !
              call libtetrabz_occstep2(de,w2)
              !
              w1(1:4,1:nb) = w1(1:4,1:nb) + V * sum(w2(1:4,1:nb,1:4), 3)
              !
           end if
           !
           ! B - 3
           !
           V = a(2,3) * a(3,1) * a(4,2)
           !
           if(V > thr) then
              !
              tmp2(1:nn - 1,1) = tmp(2:nn,1)
              tmp2(1:nn - 1,2) = tmp(2:nn,1) * a(1,3) + tmp(2:nn,3) * a(3,1) 
              tmp2(1:nn - 1,3) = tmp(2:nn,2) * a(2,3) + tmp(2:nn,3) * a(3,2) 
              tmp2(1:nn - 1,4) = tmp(2:nn,2) * a(2,4) + tmp(2:nn,4) * a(4,2) 
              !
              de(1:nb,1:4) = tmp2(1:nb, 1:4)
              w2(1:4, 1:nb, 1:4) = reshape(tmp2(nb + 1:nb + 4 * nb, 1:4), (/4, nb, 4/))
              !
              call libtetrabz_occstep2(de,w2)
              !
              w1(1:4,1:nb) = w1(1:4,1:nb) + V * sum(w2(1:4,1:nb,1:4), 3)
              !
           end if
           !
        else if( e(3) <= 0d0 .and. 0d0 < e(4)) then
           !
           ! C - 1
           !
           V = a(4,3)
           !
           if(V > thr) then
              !
              tmp2(1:nn - 1,1:3) = tmp(2:nn,1:3)
              tmp2(1:nn - 1,4)   = tmp(2:nn,3) * a(3,4) + tmp(2:nn,4) * a(4,3) 
              !
              de(1:nb,1:4) = tmp2(1:nb, 1:4)
              w2(1:4, 1:nb, 1:4) = reshape(tmp2(nb + 1:nb + 4 * nb, 1:4), (/4, nb, 4/))
              !
              call libtetrabz_occstep2(de,w2)
              !
              w1(1:4,1:nb) = w1(1:4,1:nb) + V * sum(w2(1:4,1:nb,1:4), 3)
              !
           end if
           !
           ! C - 2
           !
           V = a(3,4) * a(4,2)
           !
           if(V > thr) then
              !
              tmp2(1:nn - 1,1:2) = tmp(2:nn,1:2)
              tmp2(1:nn - 1,3)   = tmp(2:nn,2) * a(2,4) &
              &                    + tmp(2:nn,4) * a(4,2) 
              tmp2(1:nn - 1,4)   = tmp(2:nn,3) * a(3,4) &
              &                    + tmp(2:nn,4) * a(4,3) 
              !
              de(1:nb,1:4) = tmp2(1:nb, 1:4)
              w2(1:4, 1:nb, 1:4) = reshape(tmp2(nb + 1:nb + 4 * nb, 1:4), (/4, nb, 4/))
              !
              call libtetrabz_occstep2(de,w2)
              !
              w1(1:4,1:nb) = w1(1:4,1:nb) + V * sum(w2(1:4,1:nb,1:4), 3)
              !
           end if
           !
           ! C - 3
           !
           V = a(3,4) * a(2,4) * a(4,1)
           !
           if(V > thr) then
              !
              tmp2(1:nn - 1,1) = tmp(2:nn,1)
              tmp2(1:nn - 1,2) = tmp(2:nn,1) * a(1,4) + tmp(2:nn,4) * a(4,1) 
              tmp2(1:nn - 1,3) = tmp(2:nn,2) * a(2,4) + tmp(2:nn,4) * a(4,2) 
              tmp2(1:nn - 1,4) = tmp(2:nn,3) * a(3,4) + tmp(2:nn,4) * a(4,3) 
              !
              de(1:nb,1:4) = tmp2(1:nb, 1:4)
              w2(1:4, 1:nb, 1:4) = reshape(tmp2(nb + 1:nb + 4 * nb, 1:4), (/4, nb, 4/))
              !
              call libtetrabz_occstep2(de,w2)
              !
              w1(1:4,1:nb) = w1(1:4,1:nb) + V * sum(w2(1:4,1:nb,1:4), 3)
              !
           end if
           !
        else if( e(4) <= 0d0 ) then
           !
           ! D - 1
           !
           V = 1d0
           !
           tmp2(1:nn - 1,1:4) = tmp(2:nn,1:4)
           !
           de(1:nb,1:4) = tmp2(1:nb, 1:4)
           w2(1:4, 1:nb, 1:4) = reshape(tmp2(nb + 1:nb + 4 * nb, 1:4), (/4, nb, 4/))
           !
           call libtetrabz_occstep2(de,w2)
           !
           w1(1:4,1:nb) = w1(1:4,1:nb) + V * sum(w2(1:4,1:nb,1:4), 3)
           !
        end if
        !
        do ii = 1, 20
           ocst(1:nb,ib,indx2(ii,it)) = ocst(1:nb,ib,indx2(ii,it)) &
           &                         + matmul(wlsm(1:4,ii), w1(1:4,1:nb))
        end do ! ii
        !
     end do ! ib
     !$OMP END DO NOWAIT
     !
  end do ! it
  !
  !$OMP END PARALLEL
  !
  ocst(1:nb,1:nb,1:nk0) = ocst(1:nb,1:nb,1:nk0) / dble(6 * nk)
  !
end subroutine libtetrabz_occstep1
!
! Tetrahedra method for theta( - de)
!
subroutine libtetrabz_occstep2(de,w)
  !
  use libtetrabz_vals, ONLY : nb
  !
  real(8),intent(in) :: de(nb,4)
  real(8),intent(inout) :: w(4,nb,4)
  !
  integer :: ii, jj, ib
  real(8) :: V, w2(4,4), thr = 1d-8
  real(8) :: tmp(5,4), e(4), a(4,4)
  !
  do ib = 1, nb
     !
     tmp(1,1:4) = de(ib,1:4)
     tmp(2:5,1:4) = w(1:4,ib,1:4)
     call libtetrabz_sort(5, 4, tmp)
     e(1:4) = tmp(1,1:4)
     w(1:4,ib,1:4) = 0d0
     !
     do ii = 1, 4
        a(1:4,ii) = ( 0d0 - e(ii) ) / (e(1:4) - e(ii) )
     end do
     !
     if(abs(e(1)) < thr .and. abs(e(4)) < thr) then
        !
        ! Theta(0) = 0.5
        !
        V = 0.25d0 * 0.5d0
        !
        w2(1:4,1:4) = tmp(2:5, 1:4)
        !
        w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
        !
     else if((e(1) <= 0d0 .and. 0d0 < e(2)) .or. (e(1) < 0d0 .and. 0d0 <= e(2))) then
        !
        ! A - 1
        !
        V = 0.25d0 * a(2,1) * a(3,1) * a(4,1)
        !
        if(V > thr) then
           !
           w2(1:4,1) = tmp(2:5,1)
           w2(1:4,2) = tmp(2:5,1) * a(1,2) + tmp(2:5,2) * a(2,1) 
           w2(1:4,3) = tmp(2:5,1) * a(1,3) + tmp(2:5,3) * a(3,1)
           w2(1:4,4) = tmp(2:5,1) * a(1,4) + tmp(2:5,4) * a(4,1)
           !
           w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
           !
        end if
        !
     else if((e(2) <= 0d0 .and. 0d0 < e(3)) .or. (e(2) < 0d0 .and. 0d0 <= e(3))) then
        !
        ! B - 1
        !
        V = 0.25d0 * a(3,1) * a(4,1) * a(2,4)
        !
        if(V > thr) then
           !
           w2(1:4,1) = tmp(2:5,1)
           w2(1:4,2) = tmp(2:5,1) * a(1,3) + tmp(2:5,3) * a(3,1) 
           w2(1:4,3) = tmp(2:5,1) * a(1,4) + tmp(2:5,4) * a(4,1) 
           w2(1:4,4) = tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2) 
           !
           w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
           !
        end if
        !
        ! B - 2
        !
        V = 0.25d0 * a(3,2) * a(4,2)
        !
        if(V > thr) then
           !
           w2(1:4,1:2) = tmp(2:5,1:2)
           w2(1:4,3)   = tmp(2:5,2) * a(2,3) + tmp(2:5,3) * a(3,2) 
           w2(1:4,4)   = tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2) 
           !
           w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
           !
        end if
        !
        ! B - 3
        !
        V = 0.25d0 * a(2,3) * a(3,1) * a(4,2)
        !
        if(V > thr) then
           !
           w2(1:4,1) = tmp(2:5,1)
           w2(1:4,2) = tmp(2:5,1) * a(1,3) + tmp(2:5,3) * a(3,1) 
           w2(1:4,3) = tmp(2:5,2) * a(2,3) + tmp(2:5,3) * a(3,2) 
           w2(1:4,4) = tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2) 
           !
           w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
           !
        end if
        !
     else if((e(3) <= 0d0 .and. 0d0 < e(4)) .or. (e(3) < 0d0 .and. 0d0 <= e(4))) then
        !
        ! C - 1
        !
        V = 0.25d0 * a(4,3)
        !
        if(V > thr) then
           !
           w2(1:4,1:3) = tmp(2:5,1:3)
           w2(1:4,4)   = tmp(2:5,3) * a(3,4) + tmp(2:5,4) * a(4,3) 
           !
           w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
           !
        end if
        !
        ! C - 2
        !
        V = 0.25d0 * a(3,4) * a(4,2)
        !
        if(V > thr) then
           !
           w2(1:4,1:2) = tmp(2:5,1:2)
           w2(1:4,3)   = tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2) 
           w2(1:4,4)   = tmp(2:5,3) * a(3,4) + tmp(2:5,4) * a(4,3) 
           !
           w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
           !
        end if
        !
        ! C - 3
        !
        V = 0.25d0 * a(3,4) * a(2,4) * a(4,1)
        !
        if(V > thr) then
           !
           w2(1:4,1) = tmp(2:5,1)
           w2(1:4,2) = tmp(2:5,1) * a(1,4) + tmp(2:5,4) * a(4,1) 
           w2(1:4,3) = tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2) 
           w2(1:4,4) = tmp(2:5,3) * a(3,4) + tmp(2:5,4) * a(4,3) 
           !
           w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
           !
        end if
        !
     else if( e(4) <= 0d0 ) then
        !
        ! D - 1
        !
        V = 0.25d0
        !
        w2(1:4,1:4) = tmp(2:5,1:4)
        !
        w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
        !
     end if
     !
  end do
  !
end subroutine libtetrabz_occstep2
!
! Main subroutine for polalization function : Theta(- E1) * Theta(E2) / (E2 - E1)
!
subroutine libtetrabz_polstat1(eig1,eig2,pols)
  !
  use libtetrabz_vals, ONLY : nk, nk0, nb, fst, lst, indx1, indx2, wlsm
  !
  real(8),intent(in) :: eig1(nb,nk), eig2(nb,nk)
  real(8),intent(out) :: pols(nb,nb,nk0)
  !
  integer :: it, ib, ii, nn
  real(8) :: e(4), a(4,4), V, thr = 1d-10, &
  &          ei(4,nb), ei2(4), ej(nb,4), ej2(nb,4), &
  &          w0(4,nb,4), w1(4,nb), w2(4,nb,4), &
  &          tmp(1 + nb + 4 * nb,4), tmp2(1 + nb + 4 * nb,4)
  !
  nn = 1 + nb + 4 * nb
  !
  w0(1:4,1:nb,1:4) = 0d0
  do ii = 1, 4
     w0(ii,1:nb,ii) = 1d0
  end do
  !
  pols(1:nb,1:nb,1:nk0) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(fst,lst,nb,nn,indx1,indx2,wlsm,eig1,eig2,w0,pols,thr) &
  !$OMP & PRIVATE(ib,it,ii,e,a,tmp,tmp2,w1,w2,ei,ei2,ej,ej2,V)
  !
  do it = fst, lst
     !
     ei(1:4, 1:nb) = 0d0
     ej(1:nb, 1:4) = 0d0
     do ii = 1, 20
        do ib = 1, nb
           ei(1:4, ib) = ei(1:4, ib) + wlsm(1:4,ii) * eig1(ib, indx1(ii,it))
           ej(ib, 1:4) = ej(ib, 1:4) + wlsm(1:4,ii) * eig2(ib, indx1(ii,it))
        end do
     end do
     !
     !$OMP DO
     do ib = 1, nb
        !
        w1(1:4,1:nb) = 0d0
        !
        tmp(1,                      1:4) = ei(1:4,ib)
        tmp(2:1 + nb,               1:4) = ej(1:nb,1:4)
        tmp(2 + nb:1 + nb + 4 * nb, 1:4) = reshape(w0(1:4,1:nb,1:4), (/4 * nb, 4/))
        !
        call libtetrabz_sort(nn, 4, tmp)
        !
        e(1:4) = tmp(1, 1:4)
        !
        do ii = 1, 4
           a(1:4,ii) = (0d0 - e(ii) ) / (e(1:4) - e(ii) )
        end do
        !
        if( e(1) <= 0d0 .and. 0d0 < e(2) ) then
           !
           ! A - 1
           !
           V = a(2,1) * a(3,1) * a(4,1)
           !
           if(V > thr) then
              !
              tmp2(1:nn,1) = tmp(1:nn,1)
              tmp2(1:nn,2) = tmp(1:nn,1) * a(1,2) + tmp(1:nn,2) * a(2,1) 
              tmp2(1:nn,3) = tmp(1:nn,1) * a(1,3) + tmp(1:nn,3) * a(3,1)
              tmp2(1:nn,4) = tmp(1:nn,1) * a(1,4) + tmp(1:nn,4) * a(4,1)
              !
              ei2(          1:4) =         tmp2(       1,               1:4)
              ej2(    1:nb, 1:4) =         tmp2(     2:1 + nb,          1:4)
              w2(1:4, 1:nb, 1:4) = reshape(tmp2(nb + 2:1 + nb + 4 * nb, 1:4), (/4, nb, 4/))
              ! 
              call libtetrabz_polstat2(ei2,ej2,w2)
              !
              w1(1:4,1:nb) = w1(1:4,1:nb) + V * sum(w2(1:4,1:nb,1:4), 3)
              !
           end if
           !
        else if( e(2) <= 0d0 .and. 0d0 < e(3)) then
           !
           ! B - 1
           !
           V = a(3,1) * a(4,1) * a(2,4)
           !
           if(V > thr) then
              !
              tmp2(1:nn,1) = tmp(1:nn,1)
              tmp2(1:nn,2) = tmp(1:nn,1) * a(1,3) + tmp(1:nn,3) * a(3,1) 
              tmp2(1:nn,3) = tmp(1:nn,1) * a(1,4) + tmp(1:nn,4) * a(4,1) 
              tmp2(1:nn,4) = tmp(1:nn,2) * a(2,4) + tmp(1:nn,4) * a(4,2) 
              !
              ei2(          1:4) =         tmp2(1,                      1:4)
              ej2(    1:nb, 1:4) =         tmp2(2:1 + nb,               1:4)
              w2(1:4, 1:nb, 1:4) = reshape(tmp2(nb + 2:1 + nb + 4 * nb, 1:4), (/4, nb, 4/))
              ! 
              call libtetrabz_polstat2(ei2,ej2,w2)
              !
              w1(1:4,1:nb) = w1(1:4,1:nb) + V * sum(w2(1:4,1:nb,1:4), 3)
              !
           end if
           !
           ! B - 2
           !
           V = a(3,2) * a(4,2)
           !
           if(V > thr) then
              !
              tmp2(1:nn,1:2) = tmp(1:nn,1:2)
              tmp2(1:nn,3)   = tmp(1:nn,2) * a(2,3) + tmp(1:nn,3) * a(3,2) 
              tmp2(1:nn,4)   = tmp(1:nn,2) * a(2,4) + tmp(1:nn,4) * a(4,2) 
              !
              ei2(          1:4) =         tmp2(                     1, 1:4)
              ej2(    1:nb, 1:4) =         tmp2(2:1 + nb,               1:4)
              w2(1:4, 1:nb, 1:4) = reshape(tmp2(nb + 2:1 + nb + 4 * nb, 1:4), (/4, nb, 4/))
              ! 
              call libtetrabz_polstat2(ei2,ej2,w2)
              !
              w1(1:4,1:nb) = w1(1:4,1:nb) + V * sum(w2(1:4,1:nb,1:4), 3)
              !
           end if
           !
           ! B - 3
           !
           V = a(2,3) * a(3,1) * a(4,2)
           !
           if(V > thr) then
              !
              tmp2(1:nn,1) = tmp(1:nn,1)
              tmp2(1:nn,2) = tmp(1:nn,1) * a(1,3) + tmp(1:nn,3) * a(3,1) 
              tmp2(1:nn,3) = tmp(1:nn,2) * a(2,3) + tmp(1:nn,3) * a(3,2) 
              tmp2(1:nn,4) = tmp(1:nn,2) * a(2,4) + tmp(1:nn,4) * a(4,2) 
              !
              ei2(          1:4) =         tmp2(1,                      1:4)
              ej2(    1:nb, 1:4) =         tmp2(2:1 + nb,               1:4)
              w2(1:4, 1:nb, 1:4) = reshape(tmp2(nb + 2:1 + nb + 4 * nb, 1:4), (/4, nb, 4/))
              ! 
              call libtetrabz_polstat2(ei2,ej2,w2)
              !
              w1(1:4,1:nb) = w1(1:4,1:nb) + V * sum(w2(1:4,1:nb,1:4), 3)
              !
           end if
           !
        else if( e(3) <= 0d0 .and. 0d0 < e(4)) then
           !
           ! C - 1
           !
           V = a(4,3)
           !
           if(V > thr) then
              !
              tmp2(1:nn,1:3) = tmp(1:nn,1:3)
              tmp2(1:nn,4)   = tmp(1:nn,3) * a(3,4) + tmp(1:nn,4) * a(4,3) 
              !
              ei2(          1:4) =         tmp2(1,                      1:4)
              ej2(    1:nb, 1:4) =         tmp2(2:1 + nb,               1:4)
              w2(1:4, 1:nb, 1:4) = reshape(tmp2(nb + 2:1 + nb + 4 * nb, 1:4), (/4, nb, 4/))
              ! 
              call libtetrabz_polstat2(ei2,ej2,w2)
              !
              w1(1:4,1:nb) = w1(1:4,1:nb) + V * sum(w2(1:4,1:nb,1:4), 3)
              !
           end if
           !
           ! C - 2
           !
           V = a(3,4) * a(4,2)
           !
           if(V > thr) then
              !
              tmp2(1:nn,1:2) = tmp(1:nn,1:2)
              tmp2(1:nn,3)   = tmp(1:nn,2) * a(2,4) + tmp(1:nn,4) * a(4,2) 
              tmp2(1:nn,4)   = tmp(1:nn,3) * a(3,4) + tmp(1:nn,4) * a(4,3) 
              !
              ei2(          1:4) =         tmp2(1,                      1:4)
              ej2(    1:nb, 1:4) =         tmp2(2:1 + nb,               1:4)
              w2(1:4, 1:nb, 1:4) = reshape(tmp2(nb + 2:1 + nb + 4 * nb, 1:4), (/4, nb, 4/))
              ! 
              call libtetrabz_polstat2(ei2,ej2,w2)
              !
              w1(1:4,1:nb) = w1(1:4,1:nb) + V * sum(w2(1:4,1:nb,1:4), 3)
              !
           end if
           !
           ! C - 3
           !
           V = a(3,4) * a(2,4) * a(4,1)
           !
           if(V > thr) then
              !
              tmp2(1:nn,1) = tmp(1:nn,1)
              tmp2(1:nn,2) = tmp(1:nn,1) * a(1,4) + tmp(1:nn,4) * a(4,1) 
              tmp2(1:nn,3) = tmp(1:nn,2) * a(2,4) + tmp(1:nn,4) * a(4,2) 
              tmp2(1:nn,4) = tmp(1:nn,3) * a(3,4) + tmp(1:nn,4) * a(4,3) 
              !
              ei2(          1:4) =         tmp2(1,                      1:4)
              ej2(    1:nb, 1:4) =         tmp2(2:1 + nb,               1:4)
              w2(1:4, 1:nb, 1:4) = reshape(tmp2(nb + 2:1 + nb + 4 * nb, 1:4), (/4, nb, 4/))
              ! 
              call libtetrabz_polstat2(ei2,ej2,w2)
              !
              w1(1:4,1:nb) = w1(1:4,1:nb) + V * sum(w2(1:4,1:nb,1:4), 3)
              !
           end if
           !
        else if( e(4) <= 0d0 ) then
           !
           ! D - 1
           !
           V = 1d0
           !
           tmp2(1:nn,1:4) = tmp(1:nn,1:4)
           !
           ei2(          1:4) =         tmp2(1,                      1:4)
           ej2(    1:nb, 1:4) =         tmp2(2:1 + nb,               1:4)
           w2(1:4, 1:nb, 1:4) = reshape(tmp2(nb + 2:1 + nb + 4 * nb, 1:4), (/4, nb, 4/))
           !
           call libtetrabz_polstat2(ei2,ej2,w2)
           !
           w1(1:4,1:nb) = w1(1:4,1:nb) + V * sum(w2(1:4,1:nb,1:4), 3)
           !
        end if
        !
        do ii = 1, 20
           pols(1:nb,ib,indx2(ii,it)) = pols(1:nb,ib,indx2(ii,it)) &
           &                         + matmul(wlsm(1:4,ii), w1(1:4,1:nb))
        end do ! ii
        !
     end do ! ib
     !$OMP END DO NOWAIT
     !
  end do ! it
  !
  !$OMP END PARALLEL
  !
  pols(1:nb,1:nb,1:nk0) = pols(1:nb,1:nb,1:nk0) / dble(6 * nk)
  !
end subroutine libtetrabz_polstat1
!
! Tetrahedra method for theta( - E2)
!
subroutine libtetrabz_polstat2(ei,ej,w)
  !
  use libtetrabz_vals, ONLY : nb
  !
  real(8),intent(in) :: ei(4), ej(nb,4)
  real(8),intent(inout) :: w(4,nb,4)
  !
  integer :: ii, jj, ib
  real(8) :: V, de(4), w2(4,4), thr = 1d-8, &
  &          tmp(6,4), tmp2(5,4), e(4), a(4,4)
  !
  do ib = 1, nb
     !
     tmp(1,1:4) = - ej(ib,1:4)
     tmp(2,1:4) =   ej(ib,1:4) - ei(1:4)
     tmp(3:6,1:4) = w(1:4,ib,1:4)
     call libtetrabz_sort(6, 4, tmp)
     e(1:4) = tmp(1,1:4)
     w(1:4,ib,1:4) = 0d0
     !
     do ii = 1, 4
        a(1:4,ii) = (0d0 - e(ii) ) / (e(1:4) - e(ii))
     end do
     !
     if((e(1) <= 0d0 .and. 0d0 < e(2)) .or. (e(1) < 0d0 .and. 0d0 <= e(2))) then
        !
        ! A - 1
        !
        V = a(2,1) * a(3,1) * a(4,1)
        !
        if(V > thr) then
           !
           tmp2(1:5,1) = tmp(2:6,1)
           tmp2(1:5,2) = tmp(2:6,1) * a(1,2) + tmp(2:6,2) * a(2,1) 
           tmp2(1:5,3) = tmp(2:6,1) * a(1,3) + tmp(2:6,3) * a(3,1)
           tmp2(1:5,4) = tmp(2:6,1) * a(1,4) + tmp(2:6,4) * a(4,1)
           !
           de(    1:4) = tmp2(  1,1:4)
           w2(1:4,1:4) = tmp2(2:5,1:4)
           !
           call libtetrabz_polstat3(de,w2)
           !
           w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
           !
        end if
        !
     else if((e(2) <= 0d0 .and. 0d0 < e(3)) .or. (e(2) < 0d0 .and. 0d0 <= e(3))) then
        !
        ! B - 1
        !
        V = a(3,1) * a(4,1) * a(2,4)
        !
        if(V > thr) then
           !
           tmp2(1:5,1) = tmp(2:6,1)
           tmp2(1:5,2) = tmp(2:6,1) * a(1,3) + tmp(2:6,3) * a(3,1) 
           tmp2(1:5,3) = tmp(2:6,1) * a(1,4) + tmp(2:6,4) * a(4,1) 
           tmp2(1:5,4) = tmp(2:6,2) * a(2,4) + tmp(2:6,4) * a(4,2) 
           !
           de(    1:4) = tmp2(  1,1:4)
           w2(1:4,1:4) = tmp2(2:5,1:4)
           !
           call libtetrabz_polstat3(de,w2)
           !
           w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
           !
        end if
        !
        ! B - 2
        !
        V = a(3,2) * a(4,2)
        !
        if(V > thr) then
           !
           tmp2(1:5,1:2) = tmp(2:6,1:2)
           tmp2(1:5,3)   = tmp(2:6,2) * a(2,3) + tmp(2:6,3) * a(3,2) 
           tmp2(1:5,4)   = tmp(2:6,2) * a(2,4) + tmp(2:6,4) * a(4,2) 
           !
           de(    1:4) = tmp2(  1,1:4)
           w2(1:4,1:4) = tmp2(2:5,1:4)
           !
           call libtetrabz_polstat3(de,w2)
           !
           w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
           !
        end if
        !
        ! B - 3
        !
        V = a(2,3) * a(3,1) * a(4,2)
        !
        if(V > thr) then
           !
           tmp2(1:5,1) = tmp(2:6,1)
           tmp2(1:5,2) = tmp(2:6,1) * a(1,3) + tmp(2:6,3) * a(3,1) 
           tmp2(1:5,3) = tmp(2:6,2) * a(2,3) + tmp(2:6,3) * a(3,2) 
           tmp2(1:5,4) = tmp(2:6,2) * a(2,4) + tmp(2:6,4) * a(4,2) 
           !
           de(    1:4) = tmp2(  1,1:4)
           w2(1:4,1:4) = tmp2(2:5,1:4)
           !
           call libtetrabz_polstat3(de,w2)
           !
           w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
           !
        end if
        !
     else if((e(3) <= 0d0 .and. 0d0 < e(4)) .or. (e(3) < 0d0 .and. 0d0 <= e(4))) then
        !
        ! C - 1
        !
        V = a(4,3)
        !
        if(V > thr) then
           !
           tmp2(1:5,1:3) = tmp(2:6,1:3)
           tmp2(1:5,4)   = tmp(2:6,3) * a(3,4) + tmp(2:6,4) * a(4,3) 
           !
           de(    1:4) = tmp2(  1,1:4)
           w2(1:4,1:4) = tmp2(2:5,1:4)
           !
           call libtetrabz_polstat3(de,w2)
           !
           w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
           !
        end if
        !
        ! C - 2
        !
        V = a(3,4) * a(4,2)
        !
        if(V > thr) then
           !
           tmp2(1:5,1:2) = tmp(2:6,1:2)
           tmp2(1:5,3)   = tmp(2:6,2) * a(2,4) + tmp(2:6,4) * a(4,2) 
           tmp2(1:5,4)   = tmp(2:6,3) * a(3,4) + tmp(2:6,4) * a(4,3) 
           !
           de(    1:4) = tmp2(  1,1:4)
           w2(1:4,1:4) = tmp2(2:5,1:4)
           !
           call libtetrabz_polstat3(de,w2)
           !
           w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
           !
        end if
        !
        ! C - 3
        !
        V = a(3,4) * a(2,4) * a(4,1)
        !
        if(V > thr) then
           !
           tmp2(1:5,1) = tmp(2:6,1)
           tmp2(1:5,2) = tmp(2:6,1) * a(1,4) + tmp(2:6,4) * a(4,1) 
           tmp2(1:5,3) = tmp(2:6,2) * a(2,4) + tmp(2:6,4) * a(4,2) 
           tmp2(1:5,4) = tmp(2:6,3) * a(3,4) + tmp(2:6,4) * a(4,3) 
           !
           de(    1:4) = tmp2(  1,1:4)
           w2(1:4,1:4) = tmp2(2:5,1:4)
           !
           call libtetrabz_polstat3(de,w2)
           !
           w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
           !
        end if
        !
     else if( e(4) <= 0d0 ) then
        !
        ! D - 1
        !
        V = 1d0
        !
        tmp2(1:5,1:4) = tmp(2:6,1:4)
        !
        de(    1:4) = tmp2(  1,1:4)
        w2(1:4,1:4) = tmp2(2:5,1:4)
        !
        call libtetrabz_polstat3(de,w2)
        !
        w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
        !
     end if
     !
  end do
  !
end subroutine libtetrabz_polstat2
!
! Tetarahedra method for delta(om - ep + e)
!
subroutine libtetrabz_polstat3(de0,w)
  !
  real(8),intent(in) :: de0(4)
  real(8),intent(inout) :: w(4,4)
  !
  integer :: ii
  real(8) :: tmp(5,4), w2(4), de(4), lnd(4), thr, thr2
  !
  tmp(  1, 1:4) = de0(   1:4)
  tmp(2:5, 1:4) = w( 1:4,1:4)
  call libtetrabz_sort(5, 4, tmp)
  de(   1:4) = tmp(  1, 1:4)
  w(1:4,1:4) = tmp(2:5, 1:4)
  !
  thr = maxval(de(1:4)) * 1d-3
  thr2 = 1d-8
  !
  do ii = 1, 4
     if(de(ii) < thr2) then
        if(ii == 3) then
           stop "  Nesting ! "
        end if
        lnd(ii) = 0d0
        de(ii) = 0d0
     else
        lnd(ii) = log(de(ii))
     end if
  end do
  !
  if(abs(de(4) - de(3)) < thr ) then
     if(abs(de(4) - de(2)) < thr ) then
        if(abs(de(4) - de(1)) < thr ) then
           !
           ! de(4) = de(3) = de(2) = de(1)
           !
           w2(4) = 0.25d0 / de(4)
           w2(3) = w2(4)
           w2(2) = w2(4)
           w2(1) = w2(4)
           !
        else
           !
           ! de(4) = de(3) = de(2)
           !
           w2(4) = libtetrabz_polstat_1211(de(4),de(1),lnd(4),lnd(1))
           w2(3) = w2(4)
           w2(2) = w2(4)
           w2(1) = libtetrabz_polstat_1222(de(1),de(4),lnd(1),lnd(4))
           !
           if(any(w2(1:4) < 0d0)) then
              write(*,'(100e15.5)') de(1:4)
              write(*,'(100e15.5)') w2(1:4)
              stop "weighting 4=3=2"
           end if
           !
        end if
     else if(abs(de(2) - de(1)) < thr) then
        !
        ! de(4) = de(3), de(2) = de(1)
        !
        w2(4) = libtetrabz_polstat_1221(de(4),de(2), lnd(4),lnd(2))
        w2(3) = w2(4)
        w2(2) = libtetrabz_polstat_1221(de(2),de(4), lnd(2),lnd(4))
        w2(1) = w2(2)
        !
        if(any(w2(1:4) < 0d0)) then
           write(*,'(100e15.5)') de(1:4)
           write(*,'(100e15.5)') w2(1:4)
           stop "weighting 4=3 2=1"
        end if
        !
     else
        !
        ! de(4) = de(3)
        !
        w2(4) = libtetrabz_polstat_1231(de(4),de(1),de(2),lnd(4),lnd(1),lnd(2))
        w2(3) = w2(4)
        w2(2) = libtetrabz_polstat_1233(de(2),de(1),de(4),lnd(2),lnd(1),lnd(4))
        w2(1) = libtetrabz_polstat_1233(de(1),de(2),de(4),lnd(1),lnd(2),lnd(4))
        !
        if(any(w2(1:4) < 0d0)) then
           write(*,'(100e15.5)') de(1:4)
           write(*,'(100e15.5)') w2(1:4)
           stop "weighting 4=3"
        end if
        !
     end if
  else if(abs(de(3) - de(2)) < thr) then
     if(abs(de(3) - de(1)) < thr) then
        !
        ! de(3) = de(2) = de(1)
        !
        w2(4) = libtetrabz_polstat_1222(de(4),de(3), lnd(4),lnd(3))
        w2(3) = libtetrabz_polstat_1211(de(3),de(4), lnd(3),lnd(4))
        w2(2) = w2(3)
        w2(1) = w2(3)
        !
        if(any(w2(1:4) < 0d0)) then
           write(*,'(100e15.5)') de(1:4)
           write(*,'(100e15.5)') w2(1:4)
           stop "weighting 3=2=1"
        end if
        !
     else
        !
        ! de(3) = de(2)
        !
        w2(4) = libtetrabz_polstat_1233(de(4),de(1),de(3),lnd(4),lnd(1),lnd(3))
        w2(3) = libtetrabz_polstat_1231(de(3),de(1),de(4),lnd(3),lnd(1),lnd(4))
        w2(2) = w2(3)
        w2(1) = libtetrabz_polstat_1233(de(1),de(4),de(3),lnd(1),lnd(4),lnd(3))
        !
        if(any(w2(1:4) < 0d0)) then
           write(*,'(100e15.5)') de(1:4)
           write(*,'(100e15.5)') w2(1:4)
           stop "weighting 3=2"
        end if
        !
     end if
  else if(abs(de(2) - de(1)) < thr) then
     !
     ! de(2) = de(1)
     !
     w2(4) = libtetrabz_polstat_1233(de(4),de(3),de(2),lnd(4),lnd(3),lnd(2))
     w2(3) = libtetrabz_polstat_1233(de(3),de(4),de(2),lnd(3),lnd(4),lnd(2))
     w2(2) = libtetrabz_polstat_1231(de(2),de(3),de(4),lnd(2),lnd(3),lnd(4))
     w2(1) = w2(2)
     !
     if(any(w2(1:4) < 0d0)) then
        write(*,'(100e15.5)') de(1:4)
        write(*,'(100e15.5)') w2(1:4)
        stop "weighting 2=1"
     end if
     !
  else
     !
     ! Different each other.
     !
     w2(4) = libtetrabz_polstat_1234(de(4),de(1),de(2),de(3),lnd(4),lnd(1),lnd(2),lnd(3))
     w2(3) = libtetrabz_polstat_1234(de(3),de(1),de(2),de(4),lnd(3),lnd(1),lnd(2),lnd(4))
     w2(2) = libtetrabz_polstat_1234(de(2),de(1),de(3),de(4),lnd(2),lnd(1),lnd(3),lnd(4))
     w2(1) = libtetrabz_polstat_1234(de(1),de(2),de(3),de(4),lnd(1),lnd(2),lnd(3),lnd(4))
     !
     if(any(w2(1:4) < 0d0)) then
        write(*,'(100e15.5)') de(1:4)
        write(*,'(100e15.5)') w2(1:4)
        stop "weighting"
     end if
     !
  end if
  !
  do ii = 1, 4
     w(1:4,ii) = w2(ii) * w(1:4,ii)
  end do
  !
end subroutine libtetrabz_polstat3
!
! Results of Integration (1-x-y-z)/(g0+(g1-g0)x+(g2-g0)y+(g3-g0))
!  for 0<x<1, 0<y<1-x, 0<z<1-x-y
!
! 1, Different each other
!
function libtetrabz_polstat_1234(g1,g2,g3,g4,lng1,lng2,lng3,lng4) result(w)
  !
  real(8),intent(in) :: g1,g2,g3,g4,lng1,lng2,lng3,lng4
  real(8) :: w
  !
  real(8) :: w2, w3, w4
  !
  w2 = ((lng2 - lng1)/(g2 - g1)*g2 - 1d0)*g2/(g2 - g1)
  w3 = ((lng3 - lng1)/(g3 - g1)*g3 - 1d0)*g3/(g3 - g1)
  w4 = ((lng4 - lng1)/(g4 - g1)*g4 - 1d0)*g4/(g4 - g1)
  w2 = ((w2 - w3)*g2)/(g2 - g3)
  w4 = ((w4 - w3)*g4)/(g4 - g3)
  w = (w4 - w2)/(g4 - g2)
  !
end function libtetrabz_polstat_1234
!
! 2, g4 = g1
!
function libtetrabz_polstat_1231(g1,g2,g3,lng1,lng2,lng3) result(w)
  !
  real(8),intent(in) :: g1,g2,g3,lng1,lng2,lng3
  real(8) :: w
  !
  real(8) :: w2, w3
  !
  w2 = ((lng2 - lng1)/(g2 - g1)*g2 - 1d0)*g2**2/(g2 - g1) - g1/( &
  &   2d0)
  w2 = w2/(g2 - g1)
  w3 = ((lng3 - lng1)/(g3 - g1)*g3 - 1d0)*g3**2/(g3 - g1) - g1/( &
  &   2d0)
  w3 = w3/(g3 - g1)
  w = (w3 - w2)/(g3 - g2)
  !
end function libtetrabz_polstat_1231
!
! 3, g4 = g3
!
function libtetrabz_polstat_1233(g1,g2,g3,lng1,lng2,lng3) result(w)
  !
  real(8),intent(in) :: g1,g2,g3,lng1,lng2,lng3
  real(8) :: w
  !
  real(8) :: w2, w3
  !
  w2 = (lng2 - lng1)/(g2 - g1)*g2 - 1d0
  w2 = (g2*w2)/(g2 - g1)
  w3 = (lng3 - lng1)/(g3 - g1)*g3 - 1d0
  w3 = (g3*w3)/(g3 - g1)
  w2 = (w3 - w2)/(g3 - g2)
  w3 = (lng3 - lng1)/(g3 - g1)*g3 - 1d0
  w3 = 1d0 - (2d0*w3*g1)/(g3 - g1)
  w3 = w3/(g3 - g1)
  w = (g3*w3 - g2*w2)/(g3 - g2)
  !
end function libtetrabz_polstat_1233
!
! 4, g4 = g1 and g3 = g2
!
function libtetrabz_polstat_1221(g1,g2,lng1,lng2) result(w)
  !
  real(8),intent(in) :: g1, g2, lng1, lng2
  real(8) :: w
  !
  w = 1d0 - (lng2 - lng1)/(g2 - g1)*g1
  w = -1d0 + (2d0*g2*w)/(g2 - g1)
  w = -1d0 + (3d0*g2*w)/(g2 - g1)
  w = w/(2d0*(g2 - g1))
  !
end function libtetrabz_polstat_1221
!
! 5, g4 = g3 = g2
!
function libtetrabz_polstat_1222(g1,g2,lng1,lng2) result(w)
  !
  real(8),intent(in) :: g1, g2, lng1, lng2
  real(8) :: w
  !
  w = (lng2 - lng1)/(g2 - g1)*g2 - 1d0
  w = (2d0*g1*w)/(g2 - g1) - 1d0
  w = (3d0*g1*w)/(g2 - g1) + 1d0
  w = w/(2d0*(g2 - g1))
  !
end function libtetrabz_polstat_1222
!
! 6, g4 = g3 = g1
!
function libtetrabz_polstat_1211(g1,g2,lng1,lng2) result(w)
  !
  real(8),intent(in) :: g1,g2,lng1,lng2
  real(8) :: w
  !
  w = -1d0 + (lng2 - lng1)/(g2 - g1)*g2
  w = -1d0 + (2d0*g2*w)/(g2 - g1)
  w = -1d0 + (3d0*g2*w)/(2d0*(g2 - g1))
  w = w/(3d0*(g2 - g1))
  !
end function libtetrabz_polstat_1211
!
! Main subroutine for Fermi's Gorlden rule : Theta(- E1) * Theta(E2) * Delta(E2 - E1 - w)
!
subroutine libtetrabz_fermigr1(eig1,eig2,e0,fgr)
  !
  use libtetrabz_vals, ONLY : nb, nk, nk0, ne, fst, lst, indx1, indx2, wlsm
  !
  real(8),intent(in) :: eig1(nb,nk), eig2(nb,nk), e0(ne)
  real(8),intent(out) :: fgr(ne,nb,nb,nk0)
  !
  integer :: it, ib, ii, nn
  real(8) :: e(4), a(4,4), V, thr = 1d-10, ei2(4), &
  &          ei(4,nb), ej(nb,4), ej2(nb,4), &
  &          w0(4,ne,nb,4), w1(4,ne,nb), w2(4,ne,nb,4), &
  &          tmp(1 + nb + 4 * nb * ne,4), tmp2(1 + nb + 4 * nb * ne,4)
  !
  nn = 1 + nb + 4 * nb * ne
  !
  w0(1:4,1:ne,1:nb,1:4) = 0d0
  do ii = 1, 4
     w0(ii,1:ne,1:nb,ii) = 1d0
  end do
  !
  fgr(1:ne,1:nb,1:nb,1:nk0) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(fst,lst,nb,ne,nn,indx1,indx2,wlsm,eig1,eig2,w0,fgr,thr,e0) &
  !$OMP & PRIVATE(ib,it,ii,e,a,tmp,tmp2,w1,w2,ei,ej,ei2,ej2,V)
  !
  do it = fst, lst
     !
     ei(1:4, 1:nb) = 0d0
     ej(1:nb, 1:4) = 0d0
     do ii = 1, 20
        do ib = 1, nb
           ei(1:4, ib) = ei(1:4, ib) + wlsm(1:4,ii) * eig1(ib, indx1(ii,it))
           ej(ib, 1:4) = ej(ib, 1:4) + wlsm(1:4,ii) * eig2(ib, indx1(ii,it))
        end do
     end do
     !
     !$OMP DO
     do ib = 1, nb
        !
        w1(1:4,1:ne,1:nb) = 0d0
        !
        tmp(1,                           1:4) = ei(1:4,ib)
        tmp(2:1 + nb,                    1:4) = ej(1:nb,1:4)
        tmp(2 + nb:1 + nb + 4 * nb * ne, 1:4) = reshape(w0(1:4,1:ne,1:nb,1:4), (/4 * nb * ne, 4/))
        !
        call libtetrabz_sort(nn, 4, tmp)
        !
        e(1:4) = tmp(1, 1:4)
        !
        do ii = 1, 4
           a(1:4,ii) = (0d0 - e(ii) ) / (e(1:4) - e(ii) )
        end do
        !
        if( e(1) <= 0d0 .and. 0d0 < e(2) ) then
           !
           ! A - 1
           !
           V = a(2,1) * a(3,1) * a(4,1)
           !
           if(V > thr) then
              !
              tmp2(1:nn,1) = tmp(1:nn,1)
              tmp2(1:nn,2) = tmp(1:nn,1) * a(1,2) + tmp(1:nn,2) * a(2,1) 
              tmp2(1:nn,3) = tmp(1:nn,1) * a(1,3) + tmp(1:nn,3) * a(3,1)
              tmp2(1:nn,4) = tmp(1:nn,1) * a(1,4) + tmp(1:nn,4) * a(4,1)
              !
              ei2(1:4) = tmp2(1, 1:4)
              ej2(1:nb,1:4) = tmp2(2:1 + nb, 1:4)
              w2(1:4, 1:ne, 1:nb, 1:4) = reshape(tmp2(nb + 2:1 + nb + 4 * nb * ne, 1:4), &
              &                                  (/4, ne, nb, 4/))
              !
              call libtetrabz_fermigr2(e0,ei2,ej2,w2)
              !
              w1(1:4,1:ne,1:nb) = w1(1:4,1:ne,1:nb) + V * sum(w2(1:4,1:ne,1:nb,1:4), 4)
              !
           end if
           !
        else if( e(2) <= 0d0 .and. 0d0 < e(3)) then
           !
           ! B - 1
           !
           V = a(3,1) * a(4,1) * a(2,4)
           !
           if(V > thr) then
              !
              tmp2(1:nn,1) = tmp(1:nn,1)
              tmp2(1:nn,2) = tmp(1:nn,1) * a(1,3) + tmp(1:nn,3) * a(3,1) 
              tmp2(1:nn,3) = tmp(1:nn,1) * a(1,4) + tmp(1:nn,4) * a(4,1) 
              tmp2(1:nn,4) = tmp(1:nn,2) * a(2,4) + tmp(1:nn,4) * a(4,2) 
              !
              ei2(1:4) = tmp2(1, 1:4)
              ej2(1:nb,1:4) = tmp2(2:1 + nb, 1:4)
              w2(1:4, 1:ne, 1:nb, 1:4) = reshape(tmp2(nb + 2:1 + nb + 4 * nb * ne, 1:4), &
              &                                  (/4, ne, nb, 4/))
              ! 
              call libtetrabz_fermigr2(e0,ei2,ej2,w2)
              !
              w1(1:4,1:ne,1:nb) = w1(1:4,1:ne,1:nb) + V * sum(w2(1:4,1:ne,1:nb,1:4), 4)
              !
           end if
           !
           ! B - 2
           !
           V = a(3,2) * a(4,2)
           !
           if(V > thr) then
              !
              tmp2(1:nn,1:2) = tmp(1:nn,1:2)
              tmp2(1:nn,3)   = tmp(1:nn,2) * a(2,3) + tmp(1:nn,3) * a(3,2) 
              tmp2(1:nn,4)   = tmp(1:nn,2) * a(2,4) + tmp(1:nn,4) * a(4,2) 
              !
              ei2(1:4) = tmp2(1, 1:4)
              ej2(1:nb,1:4) = tmp2(2:1 + nb, 1:4)
              w2(1:4, 1:ne, 1:nb, 1:4) = reshape(tmp2(nb + 2:1 + nb + 4 * nb * ne, 1:4), &
              &                                  (/4, ne, nb, 4/))
              ! 
              call libtetrabz_fermigr2(e0,ei2,ej2,w2)
              !
              w1(1:4,1:ne,1:nb) = w1(1:4,1:ne,1:nb) + V * sum(w2(1:4,1:ne,1:nb,1:4), 4)
              !
           end if
           !
           ! B - 3
           !
           V = a(2,3) * a(3,1) * a(4,2)
           !
           if(V > thr) then
              !
              tmp2(1:nn,1) = tmp(1:nn,1)
              tmp2(1:nn,2) = tmp(1:nn,1) * a(1,3) + tmp(1:nn,3) * a(3,1) 
              tmp2(1:nn,3) = tmp(1:nn,2) * a(2,3) + tmp(1:nn,3) * a(3,2) 
              tmp2(1:nn,4) = tmp(1:nn,2) * a(2,4) + tmp(1:nn,4) * a(4,2) 
              !
              ei2(1:4) = tmp2(1, 1:4)
              ej2(1:nb,1:4) = tmp2(2:1 + nb, 1:4)
              w2(1:4, 1:ne, 1:nb, 1:4) = reshape(tmp2(nb + 2:1 + nb + 4 * nb * ne, 1:4), &
              &                                  (/4, ne, nb, 4/))
              ! 
              call libtetrabz_fermigr2(e0,ei2,ej2,w2)
              !
              w1(1:4,1:ne,1:nb) = w1(1:4,1:ne,1:nb) + V * sum(w2(1:4,1:ne,1:nb,1:4), 4)
              !
           end if
           !
        else if( e(3) <= 0d0 .and. 0d0 < e(4)) then
           !
           ! C - 1
           !
           V = a(4,3)
           !
           if(V > thr) then
              !
              tmp2(1:nn,1:3) = tmp(1:nn,1:3)
              tmp2(1:nn,4)   = tmp(1:nn,3) * a(3,4) + tmp(1:nn,4) * a(4,3) 
              !
              ei2(1:4) = tmp2(1, 1:4)
              ej2(1:nb,1:4) = tmp2(2:1 + nb, 1:4)
              w2(1:4, 1:ne, 1:nb, 1:4) = reshape(tmp2(nb + 2:1 + nb + 4 * nb * ne, 1:4), &
              &                                  (/4, ne, nb, 4/))
              ! 
              call libtetrabz_fermigr2(e0,ei2,ej2,w2)
              !
              w1(1:4,1:ne,1:nb) = w1(1:4,1:ne,1:nb) + V * sum(w2(1:4,1:ne,1:nb,1:4), 4)
              !
           end if
           !
           ! C - 2
           !
           V = a(3,4) * a(4,2)
           !
           if(V > thr) then
              !
              tmp2(1:nn,1:2) = tmp(1:nn,1:2)
              tmp2(1:nn,3)   = tmp(1:nn,2) * a(2,4) + tmp(1:nn,4) * a(4,2) 
              tmp2(1:nn,4)   = tmp(1:nn,3) * a(3,4) + tmp(1:nn,4) * a(4,3) 
              !
              ei2(1:4) = tmp2(1, 1:4)
              ej2(1:nb,1:4) = tmp2(2:1 + nb, 1:4)
              w2(1:4, 1:ne, 1:nb, 1:4) = reshape(tmp2(nb + 2:1 + nb + 4 * nb * ne, 1:4), &
              &                                  (/4, ne, nb, 4/))
              ! 
              call libtetrabz_fermigr2(e0,ei2,ej2,w2)
              !
              w1(1:4,1:ne,1:nb) = w1(1:4,1:ne,1:nb) + V * sum(w2(1:4,1:ne,1:nb,1:4), 4)
              !
           end if
           !
           ! C - 3
           !
           V = a(3,4) * a(2,4) * a(4,1)
           !
           if(V > thr) then
              !
              tmp2(1:nn,1) = tmp(1:nn,1)
              tmp2(1:nn,2) = tmp(1:nn,1) * a(1,4) + tmp(1:nn,4) * a(4,1) 
              tmp2(1:nn,3) = tmp(1:nn,2) * a(2,4) + tmp(1:nn,4) * a(4,2) 
              tmp2(1:nn,4) = tmp(1:nn,3) * a(3,4) + tmp(1:nn,4) * a(4,3) 
              !
              ei2(1:4) = tmp2(1, 1:4)
              ej2(1:nb,1:4) = tmp2(2:1 + nb, 1:4)
              w2(1:4, 1:ne, 1:nb, 1:4) = reshape(tmp2(nb + 2:1 + nb + 4 * nb * ne, 1:4), &
              &                                  (/4, ne, nb, 4/))
              ! 
              call libtetrabz_fermigr2(e0,ei2,ej2,w2)
              !
              w1(1:4,1:ne,1:nb) = w1(1:4,1:ne,1:nb) + V * sum(w2(1:4,1:ne,1:nb,1:4), 4)
              !
           end if
           !
        else if( e(4) <= 0d0 ) then
           !
           ! D - 1
           !
           V = 1d0
           !
           tmp2(1:nn,1:4) = tmp(1:nn,1:4)
           !
           ei2(1:4) = tmp2(1, 1:4)
           ej2(1:nb,1:4) = tmp2(2:1 + nb, 1:4)
           w2(1:4, 1:ne, 1:nb, 1:4) = reshape(tmp2(nb + 2:1 + nb + 4 * nb * ne, 1:4), &
           &                                  (/4, ne, nb, 4/))
           ! 
           call libtetrabz_fermigr2(e0,ei2,ej2,w2)
           !
           w1(1:4,1:ne,1:nb) = w1(1:4,1:ne,1:nb) + V * sum(w2(1:4,1:ne,1:nb,1:4), 4)
           !
        end if
        !
        do ii = 1, 20
           fgr(1:ne,1:nb,ib,indx2(ii,it)) = fgr(1:ne,1:nb,ib,indx2(ii,it)) + reshape( &
           &  matmul(wlsm(1:4,ii), reshape(w1(1:4,1:ne,1:nb), (/4, ne * nb/))), (/ne, nb/))
        end do ! ii
        !
     end do ! ib
     !$OMP END DO NOWAIT
     !
  end do ! it
  !
  !$OMP END PARALLEL
  !
  fgr(1:ne,1:nb,1:nb,1:nk0) = fgr(1:ne,1:nb,1:nb,1:nk0) / dble(6 * nk)
  !
end subroutine libtetrabz_fermigr1
!
! Tetrahedra method for theta( - E2)
!
subroutine libtetrabz_fermigr2(e0,ei,ej,w)
  !
  use libtetrabz_vals, ONLY : nb, ne
  !
  real(8),intent(in) :: e0(ne), ei(4), ej(nb,4)
  real(8),intent(inout) :: w(4,ne,nb,4)
  !
  integer :: ii, jj, ib, nn
  real(8) :: V, de(4), w2(4,ne,4), thr = 1d-8, &
  &          tmp(2 + ne * 4,4), tmp2(1 + ne * 4,4), e(4), a(4,4)
  !
  nn = 2 + ne * 4
  !
  do ib = 1, nb
     !
     tmp(1,1:4) = - ej(ib,1:4)
     tmp(2,1:4) = ej(ib,1:4) - ei(1:4)
     tmp(3:2 + 4 * ne,1:4) = reshape(w(1:4,1:ne,ib,1:4), (/4 * ne, 4/))
     call libtetrabz_sort(nn, 4, tmp)
     e(1:4) = tmp(1,1:4)
     w(1:4,1:ne,ib,1:4) = 0d0
     !
     do ii = 1, 4
        a(1:4,ii) = (0d0 - e(ii) ) / (e(1:4) - e(ii))
     end do
     !
     if((e(1) <= 0d0 .and. 0d0 < e(2)) .or. (e(1) < 0d0 .and. 0d0 <= e(2))) then
        !
        ! A - 1
        !
        V = a(2,1) * a(3,1) * a(4,1)
        !
        if(V > thr) then
           !
           tmp2(1:nn-1,1) = tmp(2:nn,1)
           tmp2(1:nn-1,2) = tmp(2:nn,1) * a(1,2) + tmp(2:nn,2) * a(2,1) 
           tmp2(1:nn-1,3) = tmp(2:nn,1) * a(1,3) + tmp(2:nn,3) * a(3,1)
           tmp2(1:nn-1,4) = tmp(2:nn,1) * a(1,4) + tmp(2:nn,4) * a(4,1)
           !
           de(         1:4) =         tmp2(1,           1:4)
           w2(1:4,1:ne,1:4) = reshape(tmp2(2:1 + 4 * ne,1:4), (/4, ne, 4/))
           !
           call libtetrabz_fermigr3(e0,de,w2)
           !
           w(1:4, 1:ne, ib,1:4) = w(1:4, 1:ne, ib, 1:4) + w2(1:4, 1:ne, 1:4) * V
           !
        end if
        !
     else if((e(2) <= 0d0 .and. 0d0 < e(3)) .or. (e(2) < 0d0 .and. 0d0 <= e(3))) then
        !
        ! B - 1
        !
        V = a(3,1) * a(4,1) * a(2,4)
        !
        if(V > thr) then
           !
           tmp2(1:nn-1,1) = tmp(2:nn,1)
           tmp2(1:nn-1,2) = tmp(2:nn,1) * a(1,3) + tmp(2:nn,3) * a(3,1) 
           tmp2(1:nn-1,3) = tmp(2:nn,1) * a(1,4) + tmp(2:nn,4) * a(4,1) 
           tmp2(1:nn-1,4) = tmp(2:nn,2) * a(2,4) + tmp(2:nn,4) * a(4,2) 
           !
           de(         1:4) =         tmp2(1,           1:4)
           w2(1:4,1:ne,1:4) = reshape(tmp2(2:1 + 4 * ne,1:4), (/4, ne, 4/))
           !
           call libtetrabz_fermigr3(e0,de,w2)
           !
           w(1:4, 1:ne, ib,1:4) = w(1:4, 1:ne, ib, 1:4) + w2(1:4, 1:ne, 1:4) * V
           !
        end if
        !
        ! B - 2
        !
        V = a(3,2) * a(4,2)
        !
        if(V > thr) then
           !
           tmp2(1:nn-1,1:2) = tmp(2:nn,1:2)
           tmp2(1:nn-1,3)   = tmp(2:nn,2) * a(2,3) + tmp(2:nn,3) * a(3,2) 
           tmp2(1:nn-1,4)   = tmp(2:nn,2) * a(2,4) + tmp(2:nn,4) * a(4,2) 
           !
           de(         1:4) =         tmp2(1,           1:4)
           w2(1:4,1:ne,1:4) = reshape(tmp2(2:1 + 4 * ne,1:4), (/4, ne, 4/))
           !
           call libtetrabz_fermigr3(e0,de,w2)
           !
           w(1:4, 1:ne, ib,1:4) = w(1:4, 1:ne, ib, 1:4) + w2(1:4, 1:ne, 1:4) * V
           !
        end if
        !
        ! B - 3
        !
        V = a(2,3) * a(3,1) * a(4,2)
        !
        if(V > thr) then
           !
           tmp2(1:nn-1,1) = tmp(2:nn,1)
           tmp2(1:nn-1,2) = tmp(2:nn,1) * a(1,3) + tmp(2:nn,3) * a(3,1) 
           tmp2(1:nn-1,3) = tmp(2:nn,2) * a(2,3) + tmp(2:nn,3) * a(3,2) 
           tmp2(1:nn-1,4) = tmp(2:nn,2) * a(2,4) + tmp(2:nn,4) * a(4,2) 
           !
           de(         1:4) =         tmp2(1,           1:4)
           w2(1:4,1:ne,1:4) = reshape(tmp2(2:1 + 4 * ne,1:4), (/4, ne, 4/))
           !
           call libtetrabz_fermigr3(e0,de,w2)
           !
           w(1:4, 1:ne, ib,1:4) = w(1:4, 1:ne, ib, 1:4) + w2(1:4, 1:ne, 1:4) * V
           !
        end if
        !
     else if((e(3) <= 0d0 .and. 0d0 < e(4)) .or. (e(3) < 0d0 .and. 0d0 <= e(4))) then
        !
        ! C - 1
        !
        V = a(4,3)
        !
        if(V > thr) then
           !
           tmp2(1:nn-1,1:3) = tmp(2:nn,1:3)
           tmp2(1:nn-1,4)   = tmp(2:nn,3) * a(3,4) + tmp(2:nn,4) * a(4,3) 
           !
           de(         1:4) =         tmp2(1,           1:4)
           w2(1:4,1:ne,1:4) = reshape(tmp2(2:1 + 4 * ne,1:4), (/4, ne, 4/))
           !
           call libtetrabz_fermigr3(e0,de,w2)
           !
           w(1:4, 1:ne, ib,1:4) = w(1:4, 1:ne, ib, 1:4) + w2(1:4, 1:ne, 1:4) * V
           !
        end if
        !
        ! C - 2
        !
        V = a(3,4) * a(4,2)
        !
        if(V > thr) then
           !
           tmp2(1:nn-1,1:2) = tmp(2:nn,1:2)
           tmp2(1:nn-1,3)   = tmp(2:nn,2) * a(2,4) + tmp(2:nn,4) * a(4,2) 
           tmp2(1:nn-1,4)   = tmp(2:nn,3) * a(3,4) + tmp(2:nn,4) * a(4,3) 
           !
           de(         1:4) =         tmp2(1,           1:4)
           w2(1:4,1:ne,1:4) = reshape(tmp2(2:1 + 4 * ne,1:4), (/4, ne, 4/))
           !
           call libtetrabz_fermigr3(e0,de,w2)
           !
           w(1:4, 1:ne, ib,1:4) = w(1:4, 1:ne, ib, 1:4) + w2(1:4, 1:ne, 1:4) * V
           !
        end if
        !
        ! C - 3
        !
        V = a(3,4) * a(2,4) * a(4,1)
        !
        if(V > thr) then
           !
           tmp2(1:nn-1,1) = tmp(2:nn,1)
           tmp2(1:nn-1,2) = tmp(2:nn,1) * a(1,4) + tmp(2:nn,4) * a(4,1) 
           tmp2(1:nn-1,3) = tmp(2:nn,2) * a(2,4) + tmp(2:nn,4) * a(4,2) 
           tmp2(1:nn-1,4) = tmp(2:nn,3) * a(3,4) + tmp(2:nn,4) * a(4,3) 
           !
           de(         1:4) =         tmp2(1,           1:4)
           w2(1:4,1:ne,1:4) = reshape(tmp2(2:1 + 4 * ne,1:4), (/4, ne, 4/))
           !
           call libtetrabz_fermigr3(e0,de,w2)
           !
           w(1:4, 1:ne, ib,1:4) = w(1:4, 1:ne, ib, 1:4) + w2(1:4, 1:ne, 1:4) * V
           !
        end if
        !
     else if( e(4) <= 0d0 ) then
        !
        ! D - 1
        !
        V = 1d0
        !
        tmp2(1:nn-1,1:4) = tmp(2:nn,1:4)
        !
        de(         1:4) =         tmp2(1,           1:4)
        w2(1:4,1:ne,1:4) = reshape(tmp2(2:1 + 4 * ne,1:4), (/4, ne, 4/))
        !
        call libtetrabz_fermigr3(e0,de,w2)
        !
        w(1:4, 1:ne, ib,1:4) = w(1:4, 1:ne, ib, 1:4) + w2(1:4, 1:ne, 1:4) * V
        !
     end if
     !
  end do
  !
end subroutine libtetrabz_fermigr2
!
!
!
subroutine libtetrabz_fermigr3(e0,de,w)
  !
  use libtetrabz_vals, only : ne
  !
  real(8),intent(in) :: e0(ne), de(4)
  real(8),intent(inout) :: w(4,ne,4)
  !
  integer :: ie, ii
  real(8) :: tmp(5,4), a(4,4), e(4), V
  !
  tmp(  1,1:4) = de(1:4)
  tmp(2:5,1:4) = w(1:4,1,1:4)
  call libtetrabz_sort(5, 4, tmp)
  e(1:4) = tmp(1,1:4)
  w(1:4,1:ne,1:4) = 0d0
  !
  do ie = 1, ne
     !
     do ii = 1, 4
        a(1:4,ii) = (e0(ie) - e(ii)) / (e(1:4) - e(ii))
     end do
     !
     if(e(1) < e0(ie) .and. e0(ie) <= e(2)) then
        !
        ! A
        !
        V = a(2,1) * a(3,1) * a(4,1) / (e0(ie) - e(1))
        !
        w(1:4,ie,1) = w(1:4,ie,1) + V * tmp(2:5,1) * (a(1,2) + a(1,3) + a(1,4))   
        w(1:4,ie,2) = w(1:4,ie,2) + V * tmp(2:5,2) * a(2,1)
        w(1:4,ie,3) = w(1:4,ie,3) + V * tmp(2:5,3) * a(3,1)
        w(1:4,ie,4) = w(1:4,ie,4) + V * tmp(2:5,4) * a(4,1)
        !
     else if(e(2) < e0(ie) .and. e0(ie) <= e(3)) then
        !
        ! B - 1
        !
        V = a(3,1) * a(4,1) * a(2,4) / (e0(ie) - e(1))
        !
        w(1:4,ie,1) = w(1:4,ie,1) + V * tmp(2:5,1) * (a(1,3) + a(1,4))   
        w(1:4,ie,2) = w(1:4,ie,2) + V * tmp(2:5,2) * a(2,4)
        w(1:4,ie,3) = w(1:4,ie,3) + V * tmp(2:5,3) * a(3,1)
        w(1:4,ie,4) = w(1:4,ie,4) + V * tmp(2:5,4) * (a(4,1) + a(4,2))
        !
        ! B - 2
        !
        V = a(2,3) * a(3,1) * a(4,2) / (e0(ie) - e(1))
        !
        w(1:4,ie,1) = w(1:4,ie,1) + V * tmp(2:5,1) * a(1,3)  
        w(1:4,ie,2) = w(1:4,ie,2) + V * tmp(2:5,2) * (a(2,3) + a(2,4))
        w(1:4,ie,3) = w(1:4,ie,3) + V * tmp(2:5,3) * (a(3,1) + a(3,2))
        w(1:4,ie,4) = w(1:4,ie,4) + V * tmp(2:5,4) * a(4,2)
        !
     else if(e(3) < e0(ie) .and. e0(ie) < e(4)) then
        !
        ! C
        !
        V = a(1,4) * a(2,4) * a(3,4) / (e(4) - e0(ie))
        !
        w(1:4,ie,1) = w(1:4,ie,1) + V * tmp(2:5,1) * a(1,4)
        w(1:4,ie,2) = w(1:4,ie,2) + V * tmp(2:5,2) * a(2,4)
        w(1:4,ie,3) = w(1:4,ie,3) + V * tmp(2:5,3) * a(3,4)
        w(1:4,ie,4) = w(1:4,ie,4) + V * tmp(2:5,4) * (a(4,1) + a(4,2) + a(4,3))
        !
     end if
     !
  end do ! ie
  !
end subroutine libtetrabz_fermigr3
!
! Main subroutine for Polaization (Imaginaly axis) : Theta(- E1) * Theta(E2) / (E2 - E1 - iw)
!
subroutine libtetrabz_polimg1(eig1,eig2,e0,poli)
  !
  use libtetrabz_vals, ONLY : nb, nk, nk0, ne, fst, lst, indx1, indx2, wlsm
  !
  real(8),intent(in) :: eig1(nb,nk), eig2(nb,nk), e0(ne)
  real(8),intent(out) :: poli(2,ne,nb,nb,nk0)
  !
  integer :: it, ib, ii, nn
  real(8) :: e(4), a(4,4), V, thr = 1d-8, &
  &          ei(4,nb), ei2(4), ej(nb,4), ej2(nb,4), &
  &          w0(4,2,ne,nb,4), w1(4,2,ne,nb), w2(4,2,ne,nb,4), &
  &          tmp(1 + nb + 8 * nb * ne,4), tmp2(1 + nb + 8 * nb * ne,4)
  !
  nn = 1 + nb + 8 * nb * ne
  !
  w0(1:4,1:2,1:ne,1:nb,1:4) = 0d0
  do ii = 1, 4
     w0(ii,1:2,1:ne,1:nb,ii) = 1d0
  end do
  !
  poli(1:2,1:ne,1:nb,1:nb,1:nk0) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(fst,lst,nb,ne,nn,indx1,indx2,wlsm,eig1,eig2,w0,poli,thr,e0) &
  !$OMP & PRIVATE(ib,it,ii,e,a,tmp,tmp2,w1,w2,ei,ei2,ej,ej2,V)
  !
  do it = fst, lst
     !
     ei(1:4, 1:nb) = 0d0
     ej(1:nb, 1:4) = 0d0
     do ii = 1, 20
        do ib = 1, nb
           ei(1:4, ib) = ei(1:4, ib) + wlsm(1:4,ii) * eig1(ib, indx1(ii,it))
           ej(ib, 1:4) = ej(ib, 1:4) + wlsm(1:4,ii) * eig2(ib, indx1(ii,it))
        end do
     end do
     !
     !$OMP DO
     do ib = 1, nb
        !
        w1(1:4,1:2,1:ne,1:nb) = 0d0
        !
        tmp(1,                           1:4) = ei(1:4,ib)
        tmp(2:1 + nb,                    1:4) = ej(1:nb,1:4)
        tmp(2 + nb:1 + nb + 8 * nb * ne, 1:4) = &
        &      reshape(w0(1:4,1:2,1:ne,1:nb,1:4), (/8 * nb * ne, 4/))
        !
        call libtetrabz_sort(nn, 4, tmp)
        !
        e(1:4) = tmp(1, 1:4)
        !
        do ii = 1, 4
           a(1:4,ii) = (0d0 - e(ii) ) / (e(1:4) - e(ii) )
        end do
        !
        if( e(1) <= 0d0 .and. 0d0 < e(2) ) then
           !
           ! A - 1
           !
           V = a(2,1) * a(3,1) * a(4,1)
           !
           if(V > thr) then
              !
              tmp2(1:nn,1) = tmp(1:nn,1)
              tmp2(1:nn,2) = tmp(1:nn,1) * a(1,2) + tmp(1:nn,2) * a(2,1) 
              tmp2(1:nn,3) = tmp(1:nn,1) * a(1,3) + tmp(1:nn,3) * a(3,1)
              tmp2(1:nn,4) = tmp(1:nn,1) * a(1,4) + tmp(1:nn,4) * a(4,1)
              !
              ei2(                     1:4) =         tmp2(1,                           1:4)
              ej2(               1:nb, 1:4) =         tmp2(2:1 + nb,                    1:4)
              w2(1:4, 1:2, 1:ne, 1:nb, 1:4) = reshape(tmp2(nb + 2:1 + nb + 8 * nb * ne, 1:4), &
              &                                  (/4, 2, ne, nb, 4/))
              !
              call libtetrabz_polimg2(e0,ei2,ej2,w2)
              !
              w1(1:4,1:2,1:ne,1:nb) = w1(1:4,1:2,1:ne,1:nb) &
              &             + V * sum(w2(1:4,1:2,1:ne,1:nb,1:4), 5)
              !
           end if
           !
        else if( e(2) <= 0d0 .and. 0d0 < e(3)) then
           !
           ! B - 1
           !
           V = a(3,1) * a(4,1) * a(2,4)
           !
           if(V > thr) then
              !
              tmp2(1:nn,1) = tmp(1:nn,1)
              tmp2(1:nn,2) = tmp(1:nn,1) * a(1,3) + tmp(1:nn,3) * a(3,1) 
              tmp2(1:nn,3) = tmp(1:nn,1) * a(1,4) + tmp(1:nn,4) * a(4,1) 
              tmp2(1:nn,4) = tmp(1:nn,2) * a(2,4) + tmp(1:nn,4) * a(4,2) 
              !
              ei2(                     1:4) =         tmp2(1,                           1:4)
              ej2(               1:nb, 1:4) =         tmp2(2:1 + nb,                    1:4)
              w2(1:4, 1:2, 1:ne, 1:nb, 1:4) = reshape(tmp2(nb + 2:1 + nb + 8 * nb * ne, 1:4), &
              &                                  (/4, 2, ne, nb, 4/))
              !
              call libtetrabz_polimg2(e0,ei2,ej2,w2)
              !
              w1(1:4,1:2,1:ne,1:nb) = w1(1:4,1:2,1:ne,1:nb) &
              &             + V * sum(w2(1:4,1:2,1:ne,1:nb,1:4), 5)
              !
           end if
           !
           ! B - 2
           !
           V = a(3,2) * a(4,2)
           !
           if(V > thr) then
              !
              tmp2(1:nn,1:2) = tmp(1:nn,1:2)
              tmp2(1:nn,3)   = tmp(1:nn,2) * a(2,3) + tmp(1:nn,3) * a(3,2) 
              tmp2(1:nn,4)   = tmp(1:nn,2) * a(2,4) + tmp(1:nn,4) * a(4,2) 
              !
              ei2(                     1:4) =         tmp2(1,                           1:4)
              ej2(               1:nb, 1:4) =         tmp2(2:1 + nb,                    1:4)
              w2(1:4, 1:2, 1:ne, 1:nb, 1:4) = reshape(tmp2(nb + 2:1 + nb + 8 * nb * ne, 1:4), &
              &                                  (/4, 2, ne, nb, 4/))
              !
              call libtetrabz_polimg2(e0,ei2,ej2,w2)
              !
              w1(1:4,1:2,1:ne,1:nb) = w1(1:4,1:2,1:ne,1:nb) &
              &             + V * sum(w2(1:4,1:2,1:ne,1:nb,1:4), 5)
              !
           end if
           !
           ! B - 3
           !
           V = a(2,3) * a(3,1) * a(4,2)
           !
           if(V > thr) then
              !
              tmp2(1:nn,1) = tmp(1:nn,1)
              tmp2(1:nn,2) = tmp(1:nn,1) * a(1,3) + tmp(1:nn,3) * a(3,1) 
              tmp2(1:nn,3) = tmp(1:nn,2) * a(2,3) + tmp(1:nn,3) * a(3,2) 
              tmp2(1:nn,4) = tmp(1:nn,2) * a(2,4) + tmp(1:nn,4) * a(4,2) 
              !
              ei2(                     1:4) =         tmp2(1,                           1:4)
              ej2(               1:nb, 1:4) =         tmp2(2:1 + nb,                    1:4)
              w2(1:4, 1:2, 1:ne, 1:nb, 1:4) = reshape(tmp2(nb + 2:1 + nb + 8 * nb * ne, 1:4), &
              &                                  (/4, 2, ne, nb, 4/))
              ! 
              call libtetrabz_polimg2(e0,ei2,ej2,w2)
              !
              w1(1:4,1:2,1:ne,1:nb) = w1(1:4,1:2,1:ne,1:nb) &
              &             + V * sum(w2(1:4,1:2,1:ne,1:nb,1:4), 5)
              !
           end if
           !
        else if( e(3) <= 0d0 .and. 0d0 < e(4)) then
           !
           ! C - 1
           !
           V = a(4,3)
           !
           if(V > thr) then
              !
              tmp2(1:nn,1:3) = tmp(1:nn,1:3)
              tmp2(1:nn,4)   = tmp(1:nn,3) * a(3,4) + tmp(1:nn,4) * a(4,3) 
              !
              ei2(                     1:4) =         tmp2(1,                           1:4)
              ej2(               1:nb, 1:4) =         tmp2(2:1 + nb,                    1:4)
              w2(1:4, 1:2, 1:ne, 1:nb, 1:4) = reshape(tmp2(nb + 2:1 + nb + 8 * nb * ne, 1:4), &
              &                                  (/4, 2, ne, nb, 4/))
              ! 
              call libtetrabz_polimg2(e0,ei2,ej2,w2)
              !
              w1(1:4,1:2,1:ne,1:nb) = w1(1:4,1:2,1:ne,1:nb) &
              &             + V * sum(w2(1:4,1:2,1:ne,1:nb,1:4), 5)
              !
           end if
           !
           ! C - 2
           !
           V = a(3,4) * a(4,2)
           !
           if(V > thr) then
              !
              tmp2(1:nn,1:2) = tmp(1:nn,1:2)
              tmp2(1:nn,3)   = tmp(1:nn,2) * a(2,4) + tmp(1:nn,4) * a(4,2) 
              tmp2(1:nn,4)   = tmp(1:nn,3) * a(3,4) + tmp(1:nn,4) * a(4,3) 
              !
              ei2(                     1:4) =         tmp2(1,                           1:4)
              ej2(               1:nb, 1:4) =         tmp2(2:1 + nb,                    1:4)
              w2(1:4, 1:2, 1:ne, 1:nb, 1:4) = reshape(tmp2(nb + 2:1 + nb + 8 * nb * ne, 1:4), &
              &                                  (/4, 2, ne, nb, 4/))
              ! 
              call libtetrabz_polimg2(e0,ei2,ej2,w2)
              !
              w1(1:4,1:2,1:ne,1:nb) = w1(1:4,1:2,1:ne,1:nb) &
              &             + V * sum(w2(1:4,1:2,1:ne,1:nb,1:4), 5)
              !
           end if
           !
           ! C - 3
           !
           V = a(3,4) * a(2,4) * a(4,1)
           !
           if(V > thr) then
              !
              tmp2(1:nn,1) = tmp(1:nn,1)
              tmp2(1:nn,2) = tmp(1:nn,1) * a(1,4) + tmp(1:nn,4) * a(4,1) 
              tmp2(1:nn,3) = tmp(1:nn,2) * a(2,4) + tmp(1:nn,4) * a(4,2) 
              tmp2(1:nn,4) = tmp(1:nn,3) * a(3,4) + tmp(1:nn,4) * a(4,3) 
              !
              ei2(                     1:4) =         tmp2(1,                           1:4)
              ej2(               1:nb, 1:4) =         tmp2(2:1 + nb,                    1:4)
              w2(1:4, 1:2, 1:ne, 1:nb, 1:4) = reshape(tmp2(nb + 2:1 + nb + 8 * nb * ne, 1:4), &
              &                                  (/4, 2, ne, nb, 4/))
              ! 
              call libtetrabz_polimg2(e0,ei2,ej2,w2)
              !
              w1(1:4,1:2,1:ne,1:nb) = w1(1:4,1:2,1:ne,1:nb) &
              &             + V * sum(w2(1:4,1:2,1:ne,1:nb,1:4), 5)
              !
           end if
           !
        else if( e(4) <= 0d0 ) then
           !
           ! D - 1
           !
           V = 1d0
           !
           tmp2(1:nn,1:4) = tmp(1:nn,1:4)
           !
           ei2(                     1:4) =         tmp2(1,                           1:4)
           ej2(               1:nb, 1:4) =         tmp2(2:1 + nb,                    1:4)
           w2(1:4, 1:2, 1:ne, 1:nb, 1:4) = reshape(tmp2(nb + 2:1 + nb + 8 * nb * ne, 1:4), &
           &                                  (/4, 2, ne, nb, 4/))
           !
           call libtetrabz_polimg2(e0,ei2,ej2,w2)
           !
           w1(1:4,1:2,1:ne,1:nb) = w1(1:4,1:2,1:ne,1:nb) &
           &             + V * sum(w2(1:4,1:2,1:ne,1:nb,1:4), 5)
           !
        end if
        !
        do ii = 1, 20
           poli(1:2,1:ne,1:nb,ib,indx2(ii,it)) = poli(1:2,1:ne,1:nb,ib,indx2(ii,it)) + reshape( &
           &  matmul(wlsm(1:4,ii), reshape(w1(1:4,1:2,1:ne,1:nb), (/4, 2 * ne * nb/))), &
           &  (/2, ne, nb/))
        end do ! ii
        !
     end do ! ib
     !$OMP END DO NOWAIT
     !
  end do ! it
  !
  !$OMP END PARALLEL
  !
  poli(1:2,1:ne,1:nb,1:nb,1:nk0) = poli(1:2,1:ne,1:nb,1:nb,1:nk0) / dble(6 * nk)
  !
end subroutine libtetrabz_polimg1
!
! Tetrahedra method for theta( - E2)
!
subroutine libtetrabz_polimg2(e0,ei,ej,w)
  !
  use libtetrabz_vals, ONLY : nb, ne
  !
  real(8),intent(in) :: e0(ne), ei(4), ej(nb,4)
  real(8),intent(inout) :: w(4,2,ne,nb,4)
  !
  integer :: ii, jj, ib, nn
  real(8) :: V, de(4), w2(4,2,ne,4), thr = 1d-8, &
  &          tmp(2 + 8 * ne,4), tmp2(1 + 8 * ne,4), e(4), a(4,4)
  !
  nn = 2 + 8 * ne
  !
  do ib = 1, nb
     !
     tmp(1,           1:4) = - ej(ib,1:4)
     tmp(2,           1:4) = ej(ib,1:4) - ei(1:4)
     tmp(3:2 + 8 * ne,1:4) = reshape(w(1:4,1:2,1:ne,ib,1:4), (/8 * ne, 4/))
     call libtetrabz_sort(nn, 4, tmp)
     e(                1:4) = tmp(1,1:4)
     w(1:4,1:2,1:ne,ib,1:4) = 0d0
     !
     do ii = 1, 4
        a(1:4,ii) = (0d0 - e(ii) ) / (e(1:4) - e(ii))
     end do
     !
     if((e(1) <= 0d0 .and. 0d0 < e(2)) .or. (e(1) < 0d0 .and. 0d0 <= e(2))) then
        !
        ! A - 1
        !
        V = a(2,1) * a(3,1) * a(4,1)
        !
        if(V > thr) then
           !
           tmp2(1:nn-1,1) = tmp(2:nn,1)
           tmp2(1:nn-1,2) = tmp(2:nn,1) * a(1,2) + tmp(2:nn,2) * a(2,1) 
           tmp2(1:nn-1,3) = tmp(2:nn,1) * a(1,3) + tmp(2:nn,3) * a(3,1)
           tmp2(1:nn-1,4) = tmp(2:nn,1) * a(1,4) + tmp(2:nn,4) * a(4,1)
           !
           de(             1:4) =         tmp2(1,           1:4)
           w2(1:4,1:2,1:ne,1:4) = reshape(tmp2(2:1 + 8 * ne,1:4), (/4, 2, ne, 4/))
           !
           call libtetrabz_polimg3(e0,de,w2)
           !
           w(1:4, 1:2, 1:ne, ib, 1:4) = w(1:4, 1:2, 1:ne, ib, 1:4) + w2(1:4, 1:2, 1:ne, 1:4) * V
           !
        end if
        !
     else if((e(2) <= 0d0 .and. 0d0 < e(3)) .or. (e(2) < 0d0 .and. 0d0 <= e(3))) then
        !
        ! B - 1
        !
        V = a(3,1) * a(4,1) * a(2,4)
        !
        if(V > thr) then
           !
           tmp2(1:nn-1,1) = tmp(2:nn,1)
           tmp2(1:nn-1,2) = tmp(2:nn,1) * a(1,3) + tmp(2:nn,3) * a(3,1) 
           tmp2(1:nn-1,3) = tmp(2:nn,1) * a(1,4) + tmp(2:nn,4) * a(4,1) 
           tmp2(1:nn-1,4) = tmp(2:nn,2) * a(2,4) + tmp(2:nn,4) * a(4,2) 
           !
           de(             1:4) =         tmp2(1,           1:4)
           w2(1:4,1:2,1:ne,1:4) = reshape(tmp2(2:1 + 8 * ne,1:4), (/4, 2, ne, 4/))
           !
           call libtetrabz_polimg3(e0,de,w2)
           !
           w(1:4, 1:2, 1:ne, ib, 1:4) = w(1:4, 1:2, 1:ne, ib, 1:4) + w2(1:4, 1:2, 1:ne, 1:4) * V
           !
        end if
        !
        ! B - 2
        !
        V = a(3,2) * a(4,2)
        !
        if(V > thr) then
           !
           tmp2(1:nn-1,1:2) = tmp(2:nn,1:2)
           tmp2(1:nn-1,3)   = tmp(2:nn,2) * a(2,3) + tmp(2:nn,3) * a(3,2) 
           tmp2(1:nn-1,4)   = tmp(2:nn,2) * a(2,4) + tmp(2:nn,4) * a(4,2) 
           !
           de(             1:4) =         tmp2(1,           1:4)
           w2(1:4,1:2,1:ne,1:4) = reshape(tmp2(2:1 + 8 * ne,1:4), (/4, 2, ne, 4/))
           !
           call libtetrabz_polimg3(e0,de,w2)
           !
           w(1:4, 1:2, 1:ne, ib, 1:4) = w(1:4, 1:2, 1:ne, ib, 1:4) + w2(1:4, 1:2, 1:ne, 1:4) * V
           !
        end if
        !
        ! B - 3
        !
        V = a(2,3) * a(3,1) * a(4,2)
        !
        if(V > thr) then
           !
           tmp2(1:nn-1,1) = tmp(2:nn,1)
           tmp2(1:nn-1,2) = tmp(2:nn,1) * a(1,3) + tmp(2:nn,3) * a(3,1) 
           tmp2(1:nn-1,3) = tmp(2:nn,2) * a(2,3) + tmp(2:nn,3) * a(3,2) 
           tmp2(1:nn-1,4) = tmp(2:nn,2) * a(2,4) + tmp(2:nn,4) * a(4,2) 
           !
           de(             1:4) =         tmp2(1,           1:4)
           w2(1:4,1:2,1:ne,1:4) = reshape(tmp2(2:1 + 8 * ne,1:4), (/4, 2, ne, 4/))
           !
           call libtetrabz_polimg3(e0,de,w2)
           !
           w(1:4, 1:2, 1:ne, ib, 1:4) = w(1:4, 1:2, 1:ne, ib, 1:4) + w2(1:4, 1:2, 1:ne, 1:4) * V
           !
        end if
        !
     else if((e(3) <= 0d0 .and. 0d0 < e(4)) .or. (e(3) < 0d0 .and. 0d0 <= e(4))) then
        !
        ! C - 1
        !
        V = a(4,3)
        !
        if(V > thr) then
           !
           tmp2(1:nn-1,1:3) = tmp(2:nn,1:3)
           tmp2(1:nn-1,4)   = tmp(2:nn,3) * a(3,4) + tmp(2:nn,4) * a(4,3) 
           !
           de(             1:4) =         tmp2(1,           1:4)
           w2(1:4,1:2,1:ne,1:4) = reshape(tmp2(2:1 + 8 * ne,1:4), (/4, 2, ne, 4/))
           !
           call libtetrabz_polimg3(e0,de,w2)
           !
           w(1:4, 1:2, 1:ne, ib, 1:4) = w(1:4, 1:2, 1:ne, ib, 1:4) + w2(1:4, 1:2, 1:ne, 1:4) * V
           !
        end if
        !
        ! C - 2
        !
        V = a(3,4) * a(4,2)
        !
        if(V > thr) then
           !
           tmp2(1:nn-1,1:2) = tmp(2:nn,1:2)
           tmp2(1:nn-1,3)   = tmp(2:nn,2) * a(2,4) + tmp(2:nn,4) * a(4,2) 
           tmp2(1:nn-1,4)   = tmp(2:nn,3) * a(3,4) + tmp(2:nn,4) * a(4,3) 
           !
           de(             1:4) =         tmp2(1,           1:4)
           w2(1:4,1:2,1:ne,1:4) = reshape(tmp2(2:1 + 8 * ne,1:4), (/4, 2, ne, 4/))
           !
           call libtetrabz_polimg3(e0,de,w2)
           !
           w(1:4, 1:2, 1:ne, ib, 1:4) = w(1:4, 1:2, 1:ne, ib, 1:4) + w2(1:4, 1:2, 1:ne, 1:4) * V
           !
        end if
        !
        ! C - 3
        !
        V = a(3,4) * a(2,4) * a(4,1)
        !
        if(V > thr) then
           !
           tmp2(1:nn-1,1) = tmp(2:nn,1)
           tmp2(1:nn-1,2) = tmp(2:nn,1) * a(1,4) + tmp(2:nn,4) * a(4,1) 
           tmp2(1:nn-1,3) = tmp(2:nn,2) * a(2,4) + tmp(2:nn,4) * a(4,2) 
           tmp2(1:nn-1,4) = tmp(2:nn,3) * a(3,4) + tmp(2:nn,4) * a(4,3) 
           !
           de(             1:4) =         tmp2(1,           1:4)
           w2(1:4,1:2,1:ne,1:4) = reshape(tmp2(2:1 + 8 * ne,1:4), (/4, 2, ne, 4/))
           !
           call libtetrabz_polimg3(e0,de,w2)
           !
           w(1:4, 1:2, 1:ne, ib, 1:4) = w(1:4, 1:2, 1:ne, ib, 1:4) + w2(1:4, 1:2, 1:ne, 1:4) * V
           !
        end if
        !
     else if( e(4) <= 0d0 ) then
        !
        ! D - 1
        !
        V = 1d0
        !
        tmp2(1:nn-1,1:4) = tmp(2:nn,1:4)
        !
        de(             1:4) =         tmp2(1,           1:4)
        w2(1:4,1:2,1:ne,1:4) = reshape(tmp2(2:1 + 8 * ne,1:4), (/4, 2, ne, 4/))
        !
        call libtetrabz_polimg3(e0,de,w2)
        !
        w(1:4, 1:2, 1:ne, ib, 1:4) = w(1:4, 1:2, 1:ne, ib, 1:4) + w2(1:4, 1:2, 1:ne, 1:4) * V
        !
     end if
     !
  end do
  !
end subroutine libtetrabz_polimg2
!
! Tetarahedra method for delta(om - ep + e)
!
subroutine libtetrabz_polimg3(e0,de,w)
  !
  use libtetrabz_vals, ONLY : ne
  !
  real(8),intent(in) :: e0(ne), de(4)
  real(8),intent(inout) :: w(4,2,ne,4)
  !
  integer :: ii, ie, ierr
  real(8) :: tmp(1 + 8 * ne,4), w2(2,4), e(4), lnd(4), thr
  !
  tmp(  1, 1:4) = de(1:4)
  tmp(2:1 + 8 * ne, 1:4) = reshape(w(1:4, 1:2, 1:ne, 1:4), (/8 * ne, 4/))
  call libtetrabz_sort(1 + 8 * ne, 4, tmp)
  w(1:4,1:2,1:ne,1:4) = reshape(tmp(2:1 + 8 * ne, 1:4), (/4, 2, ne, 4/))
  !
  do ie = 1, ne
     !
     e(1:4) = tmp(1, 1:4) / e0(ie)
     !thr = maxval(de(1:4)) * 1d-3
     thr = max(1d-3,  maxval(e(1:4)) * 1d-2)
     !
     if(abs(e(4) - e(3)) < thr ) then
        if(abs(e(4) - e(2)) < thr ) then
           if(abs(e(4) - e(1)) < thr ) then
              !
              ! e(4) = e(3) = e(2) = e(1)
              !
              w2(1,4) = 0.25d0 * e(4) / ((1d0 + e(4)**2))
              w2(2,4) = 0.25d0        / ((1d0 + e(4)**2))
              w2(1:2,3) = w2(1:2,4)
              w2(1:2,2) = w2(1:2,4)
              w2(1:2,1) = w2(1:2,4)
              !
           else
              !
              ! e(4) = e(3) = e(2)
              !
              w2(1:2,4) = libtetrabz_polimg_1211(e(4),e(1))
              w2(1:2,3) = w2(1:2,4)
              w2(1:2,2) = w2(1:2,4)
              w2(1:2,1) = libtetrabz_polimg_1222(e(1),e(4))
              !
              if(any(w2(1:2,1:4) < 0d0)) then
                 write(*,*) ie
                 write(*,'(100e15.5)') e(1:4)
                 write(*,'(2e15.5)') w2(1:2,1:4)
                 stop "weighting 4=3=2"
              end if
              !
           end if
        else if(abs(e(2) - e(1)) < thr ) then
           !
           ! e(4) = e(3), e(2) = e(1)
           !
           w2(1:2,4) = libtetrabz_polimg_1221(e(4),e(2))
           w2(1:2,3) = w2(1:2,4)
           w2(1:2,2) = libtetrabz_polimg_1221(e(2),e(4))
           w2(1:2,1) = w2(1:2,2)
           !
           if(any(w2(1:2,1:4) < 0d0)) then
              write(*,*) ie
              write(*,'(100e15.5)') e(1:4)
              write(*,'(2e15.5)') w2(1:2,1:4)
              stop "weighting 4=3 2=1"
           end if
           !
        else
           !
           ! e(4) = e(3)
           !
           w2(1:2,4) = libtetrabz_polimg_1231(e(4),e(1),e(2))
           w2(1:2,3) = w2(1:2,4)
           w2(1:2,2) = libtetrabz_polimg_1233(e(2),e(1),e(4))
           w2(1:2,1) = libtetrabz_polimg_1233(e(1),e(2),e(4))
           !
           if(any(w2(1:2,1:4) < 0d0)) then
              write(*,*) ie
              write(*,'(100e15.5)') e(1:4)
              write(*,'(2e15.5)') w2(1:2,1:4)
              stop "weighting 4=3"
           end if
           !
        end if
     else if(abs(e(3) - e(2)) < thr) then
        if(abs(e(3) - e(1)) < thr) then
           !
           ! e(3) = e(2) = e(1)
           !
           w2(1:2,4) = libtetrabz_polimg_1222(e(4),e(3))
           w2(1:2,3) = libtetrabz_polimg_1211(e(3),e(4))
           w2(1:2,2) = w2(1:2,3)
           w2(1:2,1) = w2(1:2,3)
           !
           if(any(w2(1:2,1:4) < 0d0)) then
              write(*,*) ie
              write(*,'(100e15.5)') e(1:4)
              write(*,'(2e15.5)') w2(1:2,1:4)
              stop "weighting 3=2=1"
           end if
           !
        else
           !
           ! e(3) = e(2)
           !
           w2(1:2,4) = libtetrabz_polimg_1233(e(4),e(1),e(3))
           w2(1:2,3) = libtetrabz_polimg_1231(e(3),e(1),e(4))
           w2(1:2,2) = w2(1:2,3)
           w2(1:2,1) = libtetrabz_polimg_1233(e(1),e(4),e(3))
           !
           if(any(w2(1:2,1:4) < 0d0)) then
              write(*,*) ie
              write(*,'(100e15.5)') e(1:4)
              write(*,'(2e15.5)') w2(1:2,1:4)
              stop "weighting 3=2"
           end if
           !
        end if
     else if(abs(e(2) - e(1)) < thr) then
        !
        ! e(2) = e(1)
        !
        w2(1:2,4) = libtetrabz_polimg_1233(e(4),e(3),e(2))
        w2(1:2,3) = libtetrabz_polimg_1233(e(3),e(4),e(2))
        w2(1:2,2) = libtetrabz_polimg_1231(e(2),e(3),e(4))
        w2(1:2,1) = w2(1:2,2)
        !
        if(any(w2(1:2,1:4) < 0d0)) then
           write(*,*) ie
           write(*,'(100e15.5)') e(1:4)
           write(*,'(2e15.5)') w2(1:2,1:4)
           stop "weighting 2=1"
        end if
        !
     else
        !
        ! Different each other.
        !
        w2(1:2,4) = libtetrabz_polimg_1234(e(4),e(1),e(2),e(3))
        w2(1:2,3) = libtetrabz_polimg_1234(e(3),e(1),e(2),e(4))
        w2(1:2,2) = libtetrabz_polimg_1234(e(2),e(1),e(3),e(4))
        w2(1:2,1) = libtetrabz_polimg_1234(e(1),e(2),e(3),e(4))
        !
        if(any(w2(1:2,1:4) < 0d0)) then
           write(*,*) ie
           write(*,'(100e15.5)') e(1:4)
           write(*,'(2e15.5)') w2(1:2,1:4)
           stop "weighting"
        end if
        !
     end if
     !
     do ii = 1, 4
        w(1:4,1,ie,ii) = w2(1,ii) * w(1:4,1,ie,ii) /    e0(ie)
        w(1:4,2,ie,ii) = w2(2,ii) * w(1:4,2,ie,ii) / (- e0(ie))
     end do ! ii
     !
  end do ! ie
  !
end subroutine libtetrabz_polimg3
!
! Results of Integration (1-x-y-z)/(g0+(g1-g0)x+(g2-g0)y+(g3-g0))
!  for 0<x<1, 0<y<1-x, 0<z<1-x-y
!
! 1, Different each other
!
function libtetrabz_polimg_1234(g1,g2,g3,g4) result(w)
  !
  real(8),intent(in) :: g1, g2, g3, g4
  real(8) :: w(2)
  !
  real(8) :: w2, w3, w4
  !
  ! Real
  !
  w2 = 2d0*(3d0*g2**2 - 1d0)*(atan(g2) - atan(g1)) + (g2**2 - &
  &      3d0)*g2*log((1d0 + g2**2)/( 1d0 + g1**2))
  w2 = -2d0*(g2**2 - 1d0) + w2/(g2 - g1 )
  w2 = w2/(g2 - g1 )
  w3 = 2d0*(3d0*g3**2 - 1d0)*(atan(g3) - atan(g1)) + (g3**2 -  &
  &      3d0)*g3*log((1d0 + g3**2)/( 1d0 + g1**2))
  w3 = -2d0*(g3**2 - 1d0) + w3/(g3 - g1 )
  w3 = w3/(g3 - g1 )
  w4 = 2d0*(3d0*g4**2 - 1d0)*(atan(g4) - atan(g1)) + (g4**2 -  &
  &      3d0)*g4*log((1d0 + g4**2)/( 1d0 + g1**2))
  w4 = -2d0*(g4**2 - 1d0) + w4/(g4 - g1 )
  w4 = w4/(g4 - g1 )
  w2 = (w2 - w3)/(g2 - g3)
  w4 = (w4 - w3)/(g4 - g3)
  w(1) = (w4 - w2)/(2d0*(g4 - g2))
  !
  ! Imaginal
  !
  w2 = 2d0*(3d0 - g2**2)* &
  &    g2*(atan(g2) - atan(g1)) + (3d0*g2**2 - 1d0)* &
  &    log((1d0 + g2**2)/(1d0 + g1**2))
  w2 = 4d0*g2 - w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2d0*(3d0 - g3**2)* &
  &    g3*(atan(g3) - atan(g1)) + (3d0*g3**2 - 1d0)* &
  &    log((1d0 + g3**2)/(1d0 + g1**2))
  w3 = 4d0*g3 - w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w4 = 2d0*(3d0 - g4**2)* &
  &    g4*(atan(g4) - atan(g1)) + (3d0*g4**2 - 1d0)* &
  &    log((1d0 + g4**2)/(1d0 + g1**2))
  w4 = 4d0*g4 - w4/(g4 - g1)
  w4 = w4/(g4 - g1)
  w2 = (w2 - w3)/(g2 - g3)
  w4 = (w4 - w3)/(g4 - g3)
  w(2) = (w4 - w2)/(2d0*(g4 - g2))
  !
end function libtetrabz_polimg_1234
!
! 2, g4 = g1
!
function libtetrabz_polimg_1231(g1,g2,g3) result(w)
  !
  real(8),intent(in) :: g1, g2, g3
  real(8) :: w(2)
  !
  real(8) :: w2, w3
  !
  ! Real
  !
  w2 = 2d0*(-1d0 + 3d0*g2**2)*(atan(g2) - atan(g1)) +  &
  &   g2*(-3d0 + g2**2)*log((1d0 + g2**2)/(1d0 + g1**2))
  w2 = 2d0*(1d0 - g2**2) + w2/(g2 - g1)
  w2 = -g1 + w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2d0*(-1d0 + 3d0*g3**2)*(atan(g3) - atan(g1)) +  &
  &   g3*(-3d0 + g3**2)*log((1d0 + g3**2)/(1d0 + g1**2))
  w3 = 2d0*(1 - g3**2) + w3/(g3 - g1)
  w3 = -g1 + w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w(1) = (w3 - w2)/(2d0*(g3 - g2))
  !
  ! Imaginal
  !
  w2 = 2d0* &
  &    g2*(3d0 - g2**2)*(atan(g2) - atan(g1)) + (-1d0 + 3d0*g2**2)* &
  &    log((1d0 + g2**2)/(1d0 + g1**2))
  w2 = 4d0*g2 - w2/(g2 - g1)
  w2 = 1 + w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2d0* &
  &    g3*(3d0 - g3**2)*(atan(g3) - atan(g1)) + (-1d0 + 3d0*g3**2)* &
  &    log((1d0 + g3**2)/(1d0 + g1**2))
  w3 = 4d0*g3 - w3/(g3 - g1)
  w3 = 1 + w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w(2) = (w3 - w2)/(2d0*(g3 - g2))
  !
end function libtetrabz_polimg_1231
!
! 3, g4 = g3
!
function libtetrabz_polimg_1233(g1, g2, g3) result(w)
  !
  real(8),intent(in) :: g1, g2, g3
  real(8) :: w(2)
  !
  real(8) :: w2, w3
  !
  ! Real
  !
  w2 = 2d0*(1d0 - 3d0*g2**2)*(atan(g2) - atan(g1)) +  &
  &   g2*(3d0 - g2**2)*log((1d0 + g2**2)/(1d0 + g1**2))
  w2 = 2d0*(1 - g2**2) - w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2d0*(1d0 - 3d0*g3**2)*(atan(g3) - atan(g1)) +  &
  &   g3*(3d0 - g3**2)*log((1d0 + g3**2)/(1d0 + g1**2))
  w3 = 2d0*(1 - g3**2) - w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w2 = (w3 - w2)/(g3 - g2)
  w3 = 4d0*(1d0 - 3d0*g1*g3)*(atan(g3) - atan(g1)) + (3d0*g1 +  &
  &      3d0*g3 - 3d0*g1*g3**2 + g3**3) * log((1d0 + g3**2)/( &
  &     1d0 + g1**2))
  w3 = -4d0*(1d0 - g1**2) + w3/(g3 - g1)
  w3 = 4d0*g1 + w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w(1) = (w3 - w2)/(2d0*(g3 - g2))
  !
  ! Imaginal
  !
  w2 = 2d0* &
  &    g2*(3d0 - g2**2)*(atan(g2) - atan(g1)) + (-1d0 + 3d0*g2**2)* &
  &    log((1d0 + g2**2)/(1d0 + g1**2))
  w2 = 4d0*g2 - w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2d0* &
  &    g3*(3d0 - g3**2)*(atan(g3) - atan(g1)) + (-1d0 + 3d0*g3**2)* &
  &    log((1d0 + g3**2)/(1d0 + g1**2))
  w3 = 4d0*g3 - w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w2 = (w3 - w2)/(g3 - g2)
  w3 = (3d0*g1 - 3d0*g1*g3**2 + 3d0*g3 + g3**3)*(atan(g3) -  &
  &      atan(g1)) + (3d0*g1*g3 - 1d0)* &
  &    log((1d0 + g3**2)/(1d0 + g1**2))
  w3 = w3/(g3 - g1) - 4d0*g1
  w3 = w3/(g3 - g1) - 2d0
  w3 = (2d0*w3)/(g3 - g1)
  w(2) = (w3 - w2)/(2d0*(g3 - g2))
  !
end function libtetrabz_polimg_1233
!
! 4, g4 = g1 and g3 = g2
!
function libtetrabz_polimg_1221(g1,g2) result(w)
  !
  real(8),intent(in) :: g1, g2
  real(8) :: w(2)
  !
  ! Real
  !
  w(1) = -2d0*(-1d0 + 2d0*g1*g2 + g2**2)*(atan(g2) -  &
  &      atan(g1)) + (g1 + 2d0*g2 - g1*g2**2)* &
  &    log((1d0 + g2**2)/(1d0 + g1**2))
  w(1) = 2d0*(-1d0 + g1**2) + w(1)/(g2 - g1)
  w(1) = 3d0*g1 + w(1)/(g2 - g1)
  w(1) = 2d0 + (3d0*w(1))/(g2 - g1)
  w(1) = w(1)/(2d0*(g2 - g1))
  !
  ! Imaginal
  !
  w(2) = 2d0*(g1 + 2d0*g2 - g1*g2**2)*(atan(g2) -  &
  &      atan(g1)) + (-1d0 + 2d0*g1*g2 + g2**2)* &
  &    log((1 + g2**2)/(1 + g1**2))
  w(2) = -4d0*g1 + w(2)/(g2 - g1)
  w(2) = -3d0 + w(2)/(g2 - g1)
  w(2) = (3d0*w(2))/(2d0*(g2 - g1)**2)
  !
end function libtetrabz_polimg_1221
!
! 5, g4 = g3 = g2
!
function libtetrabz_polimg_1222(g1,g2) result(w)
  !
  real(8),intent(in) :: g1, g2
  real(8) :: w(2)
  !
  ! Real
  !
  w(1) = 2d0*(-1d0 + g1**2 + 2d0*g1*g2)*(atan(g2) -  &
  &      atan(g1)) + (-2d0*g1 - g2 + g1**2*g2) * log((1d0 + g2**2)/( &
  &     1d0 + g1**2))
  w(1) = 2d0*(1d0 - g1**2) + w(1)/(g2 - g1)
  w(1) = g1 - w(1)/(g2 - g1)
  w(1) = 1d0 - (3d0*w(1))/(g2 - g1)
  w(1) = w(1)/(2d0*(g2 - g1))
  !
  ! Imaginal
  !
  w(2) = 2d0*(-2d0*g1 - g2 + g1**2*g2)*(atan(g2) - atan(g1)) + (1d0 - &
  &       g1**2 - 2d0*g1*g2) * log((1d0 + g2**2)/(1d0 + g1**2))
  w(2) = 4d0*g1 + w(2)/(g2 - g1)
  w(2) = 1d0 + w(2)/(g2 - g1)
  w(2) = (3d0*w(2))/(2d0*(g2 - g1)**2)
  !
end function libtetrabz_polimg_1222
!
! 6, g4 = g3 = g1
!
function libtetrabz_polimg_1211(g1,g2) result(w)
  !
  real(8),intent(in) :: g1, g2
  real(8) :: w(2)
  !
  ! Real
  !
  w(1) = 2d0*(3d0*g2**2 - 1d0)*(atan(g2) - atan(g1)) +  &
  &   g2*(g2**2 - 3d0)*log((1d0 + g2**2)/(1d0 + g1**2))
  w(1) = 2d0*(1d0 - g1**2) + w(1)/(g2 - g1)
  w(1) = -5d0*g1 + w(1)/(g2 - g1)
  w(1) = -11d0 + (3d0*w(1))/(g2 - g1)
  w(1) = w(1)/(6d0*(g2 - g1))
  !
  ! Imaginal
  !
  w(2) = 2d0*g2*(-3d0 + g2**2)*(atan(g2) - atan(g1)) + (1d0 -  &
  &      3d0*g2**2)*log((1d0 + g2**2)/(1d0 + g1**2))
  w(2) = 4d0*g2 + w(2)/(g2 - g1)
  w(2) = 1d0 + w(2)/(g2 - g1)
  w(2) = w(2)/(2d0*(g2 - g1)**2)
  !
end function libtetrabz_polimg_1211
!
! Interpolate integration weight
!
subroutine libtetrabz_interpol_weight(nb,ngc,ngd,wc,wd)
  !
  use libtetrabz_vals, only : nk0, indx3
  !
  integer,intent(in) :: nb, ngc(3), ngd(3)
  real(8),intent(in) :: wd(nb,nk0)
  real(8),intent(out) :: wc(nb,product(ngc(1:3)))
  !
  integer :: i1, i2, i3, ik, nkc, nkd
  real(8) :: kv(3, product(ngd(1:3)))
  !
  nkc = product(ngc(1:3))
  nkd = product(ngd(1:3))
  !
  ik = 0
  do i3 = 1, ngd(3)
     do i2 = 1, ngd(2)
        do i1 = 1, ngd(1)
           ik = ik + 1
           kv(1:3,ik) = dble((/i1, i2, i3/) - 1) / dble(ngd(1:3))
        end do
     end do
  end do
  !
  wc(1:nb,1:nkc) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nk0,nkc,nkd,nb,ngc,kv,wc,wd,indx3) &
  !$OMP PRIVATE(ik)
  !
  !$OMP DO REDUCTION(+: wc)
  do ik = 1, nk0
     call libtetrabz_interpol_weight2(nkc, nb, ngc, kv(1:3,indx3(ik)), wd(1:nb,ik), wc)
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  !
end subroutine libtetrabz_interpol_weight
!
! first or third order interpolation of weights
!
subroutine libtetrabz_interpol_weight2(nk,nb,ng,ko,wi,wo)
  !
  use libtetrabz_vals, only : ltetra, ivvec
  !
  integer,intent(in)  :: nk, nb, ng(3)
  real(8),intent(in)  :: ko(3)
  real(8),intent(in) :: wi(nb)
  real(8),intent(inout) :: wo(nb,nk)
  !
  integer :: ikv(3), ikv1(3), ik(20), ii, it, it0, ierr
  real(8) :: rot(3,3), res(3), prod(3), u, x, y, z, thr = 1d-10
  !
  rot(1:3,1) = (/  2d0, - 1d0,   0d0/)
  rot(1:3,2) = (/- 1d0,   2d0, - 1d0/)
  rot(1:3,3) = (/  0d0, - 1d0,   1d0/)
  !
  ! Search nearest neighbor grid points.
  !
  res(1:3) = ko(1:3) * dble(ng(1:3))
  ikv(1:3) = floor(res(1:3))
  res(1:3) = res(1:3) - dble(ikv(1:3))
  !
  do it = 1, 6
     !
     do ii = 1, 3
        prod(ii) = dot_product(dble(ivvec(1:3,1 + ii,it) - ivvec(1:3,1,it)), &
        &                                  res(1:3) - dble(ivvec(1:3,1,it))  )
     end do
     !
     prod(1:3) = matmul(rot(1:3,1:3), prod(1:3))
     !
     if(minval(prod(1:3)) > - thr .and. sum(prod(1:3)) < 1d0 + thr) then
        it0 = it
        goto 10
     end if
     !
  end do
  !
  stop "interpol"
  !
10 continue
  !
  x = prod(1)
  y = prod(2)
  z = prod(3)
  u = 1d0 - x - y - z
  !
  do ii = 1, 20
     !
     ikv1(1:3) = ikv(1:3) + ivvec(1:3,ii,it0)
     ikv1(1:3) = modulo(ikv1(1:3), ng(1:3))
     ik(ii) = 1 + ikv1(1) + ikv1(2) * ng(1) + ikv1(3) * ng(1) * ng(2)
     !
  end do
  !
  if(ltetra == 0 .or. ltetra == 1) then
     !
     wo(1:nb,ik(1)) = wo(1:nb,ik(1)) + wi(1:nb) * u
     wo(1:nb,ik(2)) = wo(1:nb,ik(2)) + wi(1:nb) * x
     wo(1:nb,ik(3)) = wo(1:nb,ik(3)) + wi(1:nb) * y
     wo(1:nb,ik(4)) = wo(1:nb,ik(4)) + wi(1:nb) * z
     !
  else if(ltetra == 2) then
     !
     wo(1:nb,ik( 1)) = wo(1:nb,ik( 1)) + wi(1:nb) * 0.5d0 * u * (2d0 + u * (1d0 - u) + 2d0 * y * (x + z))
     wo(1:nb,ik( 2)) = wo(1:nb,ik( 2)) + wi(1:nb) * 0.5d0 * x * (2d0 + x * (1d0 - x) + 2d0 * z * (u + y))
     wo(1:nb,ik( 3)) = wo(1:nb,ik( 3)) + wi(1:nb) * 0.5d0 * y * (2d0 + y * (1d0 - y) + 2d0 * u * (x + z))
     wo(1:nb,ik( 4)) = wo(1:nb,ik( 4)) + wi(1:nb) * 0.5d0 * z * (2d0 + z * (1d0 - z) + 2d0 * x * (u + y))
     wo(1:nb,ik( 5)) = wo(1:nb,ik( 5)) + wi(1:nb) * x * u * (2d0 * y - u - 1d0) / 6d0
     wo(1:nb,ik( 6)) = wo(1:nb,ik( 6)) + wi(1:nb) * x * y * (2d0 * z - x - 1d0) / 6d0
     wo(1:nb,ik( 7)) = wo(1:nb,ik( 7)) + wi(1:nb) * y * z * (2d0 * u - y - 1d0) / 6d0
     wo(1:nb,ik( 8)) = wo(1:nb,ik( 8)) + wi(1:nb) * z * u * (2d0 * x - z - 1d0) / 6d0
     wo(1:nb,ik( 9)) = wo(1:nb,ik( 9)) + wi(1:nb) * y * u * (2d0 * y + u - 3d0) / 6d0
     wo(1:nb,ik(10)) = wo(1:nb,ik(10)) + wi(1:nb) * x * z * (2d0 * z + x - 3d0) / 6d0
     wo(1:nb,ik(11)) = wo(1:nb,ik(11)) + wi(1:nb) * y * u * (2d0 * u + y - 3d0) / 6d0
     wo(1:nb,ik(12)) = wo(1:nb,ik(12)) + wi(1:nb) * x * z * (2d0 * x + z - 3d0) / 6d0
     wo(1:nb,ik(13)) = wo(1:nb,ik(13)) + wi(1:nb) * z * u * (2d0 * y - u - 1d0) / 6d0
     wo(1:nb,ik(14)) = wo(1:nb,ik(14)) + wi(1:nb) * x * u * (2d0 * z - x - 1d0) / 6d0
     wo(1:nb,ik(15)) = wo(1:nb,ik(15)) + wi(1:nb) * x * y * (2d0 * u - y - 1d0) / 6d0
     wo(1:nb,ik(16)) = wo(1:nb,ik(16)) + wi(1:nb) * y * z * (2d0 * x - z - 1d0) / 6d0
     wo(1:nb,ik(17)) = wo(1:nb,ik(17)) + wi(1:nb) * (- x * z * u)
     wo(1:nb,ik(18)) = wo(1:nb,ik(18)) + wi(1:nb) * (- x * y * u)
     wo(1:nb,ik(19)) = wo(1:nb,ik(19)) + wi(1:nb) * (- x * y * z)
     wo(1:nb,ik(20)) = wo(1:nb,ik(20)) + wi(1:nb) * (- y * z * u)
     !
  else
     !
     stop "interpol2"
     ! 
  end if
  !
end subroutine libtetrabz_interpol_weight2
!
!
!
end module libtetrabz_routines
