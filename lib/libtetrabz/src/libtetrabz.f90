module libtetrabz
  !
  implicit none
  !
contains
!
! Compute occupation
!
subroutine libtetrabz_occ(ltetra0,bvec,nb0,nge,eig,ngw,wght0)
  !
  use libtetrabz_vals, only : ltetra, ng, nb, nk0, indx1, indx2, indx3
  use libtetrabz_routines, only : libtetrabz_initialize, libtetrabz_interpol_weight, &
  &                               libtetrabz_occ1
  !
  integer,intent(in) :: ltetra0, nge(3), ngw(3), nb0
  real(8),intent(in) :: bvec(3,3)
  real(8),intent(in) :: eig(nb0,product(nge(1:3)))
  real(8),intent(out) :: wght0(nb0,product(ngw(1:3)))
  !
  integer :: ierr, nn
  real(8),allocatable :: wght1(:,:)
  !
  ltetra = ltetra0
  nb = nb0
  ng(1:3) = nge(1:3)
  nn = nb
  !
  call libtetrabz_initialize(bvec)
  call libtetrabz_kgrid()
  !
  if(any(nge(1:3) /= ngw(1:3))) then
     !
     allocate(wght1(nn, nk0))
     call libtetrabz_occ1(0d0,eig,wght1)
     !
     call libtetrabz_interpol_weight(nn,ngw,nge,wght0,wght1)
     deallocate(wght1)
     !
  else
     call libtetrabz_occ1(0d0,eig,wght0)
  end if
  !
  deallocate(indx1, indx2, indx3)
  !
end subroutine libtetrabz_occ
!
! Calculate Fermi energy
!
subroutine libtetrabz_fermieng(ltetra0,bvec,nb0,nge,eig,ngw,wght0,ef,nelec)
  !
  use libtetrabz_vals, only : ltetra, ng, nb, nk, nk0, indx1, indx2, indx3
  use libtetrabz_routines, only : libtetrabz_initialize, libtetrabz_interpol_weight, &
  &                               libtetrabz_occ1
  !
  integer,intent(in) :: ltetra0, nge(3), ngw(3), nb0
  real(8),intent(in) :: bvec(3,3), nelec
  real(8),intent(in) :: eig(nb0,product(nge(1:3)))
  real(8),intent(out) :: ef
  real(8),intent(out) :: wght0(nb0,product(ngw(1:3)))
  !
  integer :: ierr, nn
  real(8),allocatable :: wght1(:,:)
  !
  integer :: iter, maxiter = 300
  real(8) :: elw, eup, sumkmid, eps= 1.0d-10
  !
  ltetra = ltetra0
  nb = nb0
  ng(1:3) = nge(1:3)
  nn = nb
  !
  call libtetrabz_initialize(bvec)
  call libtetrabz_kgrid()
  !
  if(any(nge(1:3) /= ngw(1:3))) allocate(wght1(nn, nk0))
  !
  elw = minval(eig(1:nb,1:nk))
  eup = maxval(eig(1:nb,1:nk))
  !
  ! Bisection method
  !
  do iter = 1, maxiter
     !
     ef = (eup + elw) / 2.d0
     !
     ! Calc. # of electrons 
     !
     if(any(nge(1:3) /= ngw(1:3))) then
        call libtetrabz_occ1(ef, eig,wght1)
        sumkmid = sum(wght1(1:nb,1:nk0))
     else
        call libtetrabz_occ1(ef, eig,wght0)
        sumkmid = sum(wght0(1:nb,1:nk0))
     end if
     !
     !
     ! convergence check
     !
     if(abs(sumkmid - nelec) < eps) then
        exit
     elseif(sumkmid < nelec) then
        elw = ef
     else
        eup = ef
     endif
     !
  enddo ! iter
  !
  if(iter >= maxiter) stop "libtetrabz_omp_fermieng"
  !
  if(any(nge(1:3) /= ngw(1:3))) then
     call libtetrabz_interpol_weight(nn,ngw,nge,wght0,wght1)
     deallocate(wght1)
  end if
  !
  deallocate(indx1, indx2, indx3)
  !
end subroutine libtetrabz_fermieng
!
! Compute DOS
!
subroutine libtetrabz_dos(ltetra0,bvec,nb0,nge,eig,ngw,wght0,ne0,e0)
  !
  use libtetrabz_vals, only : ltetra, ng, nb, nk0, indx1, indx2, indx3, ne
  use libtetrabz_routines, only : libtetrabz_initialize, libtetrabz_interpol_weight, &
  &                               libtetrabz_dos1
  !
  integer,intent(in) :: ltetra0, nge(3), ngw(3), nb0, ne0
  real(8),intent(in) :: bvec(3,3), e0(ne0)
  real(8),intent(in) :: eig(nb0,product(nge(1:3)))
  real(8),intent(out) :: wght0(ne0,nb0,product(ngw(1:3)))
  !
  integer :: ierr, nn
  real(8),allocatable :: wght1(:,:)
  !
  ltetra = ltetra0
  nb = nb0
  ng(1:3) = nge(1:3)
  ne = ne0
  nn = ne * nb
  !
  call libtetrabz_initialize(bvec)
  call libtetrabz_kgrid()
  !
  if(any(nge(1:3) /= ngw(1:3))) then
     !
     allocate(wght1(nn, nk0))
     call libtetrabz_dos1(eig,e0,wght1)
     !
     call libtetrabz_interpol_weight(nn,ngw,nge,wght0,wght1)
     deallocate(wght1)
     !
  else
     call libtetrabz_dos1(eig,e0,wght0)
  end if
  !
  deallocate(indx1, indx2, indx3)
  !
end subroutine libtetrabz_dos
!
! Compute doubledelta
!
subroutine libtetrabz_doubledelta(ltetra0,bvec,nb0,nge,eig1,eig2,ngw,wght0)
  !
  use libtetrabz_vals, only : ltetra, ng, nb, nk0, indx1, indx2, indx3
  use libtetrabz_routines, only : libtetrabz_initialize, libtetrabz_interpol_weight, &
  &                               libtetrabz_doubledelta1
  !
  integer,intent(in) :: ltetra0, nge(3), ngw(3), nb0
  real(8),intent(in) :: bvec(3,3)
  real(8),intent(in) :: eig1(nb0,product(nge(1:3))), eig2(nb0,product(nge(1:3)))
  real(8),intent(out) :: wght0(nb0,nb0,product(ngw(1:3)))
  !
  integer :: ierr, nn
  real(8),allocatable :: wght1(:,:)
  !
  ltetra = ltetra0
  nb = nb0
  ng(1:3) = nge(1:3)
  nn = nb * nb
  !
  call libtetrabz_initialize(bvec)
  call libtetrabz_kgrid()
  !
  if(any(nge(1:3) /= ngw(1:3))) then
     !
     allocate(wght1(nn, nk0))
     call libtetrabz_doubledelta1(eig1,eig2,wght1)
     !
     call libtetrabz_interpol_weight(nn,ngw,nge,wght0,wght1)
     deallocate(wght1)
     !
  else
     call libtetrabz_doubledelta1(eig1,eig2,wght0)
  end if
  !
  deallocate(indx1, indx2, indx3)
  !
end subroutine libtetrabz_doubledelta
!
! Compute Occ * Step
!
subroutine libtetrabz_occstep(ltetra0,bvec,nb0,nge,eig1,eig2,ngw,wght0)
  !
  use libtetrabz_vals, only : ltetra, ng, nb, nk0, indx1, indx2, indx3
  use libtetrabz_routines, only : libtetrabz_initialize, libtetrabz_interpol_weight, &
  &                               libtetrabz_occstep1
  !
  integer,intent(in) :: ltetra0, nge(3), ngw(3), nb0
  real(8),intent(in) :: bvec(3,3)
  real(8),intent(in) :: eig1(nb0,product(nge(1:3))), eig2(nb0,product(nge(1:3)))
  real(8),intent(out) :: wght0(nb0,nb0,product(ngw(1:3)))
  !
  integer :: ierr, nn
  real(8),allocatable :: wght1(:,:)
  !
  ltetra = ltetra0
  nb = nb0
  ng(1:3) = nge(1:3)
  nn = nb * nb
  !
  call libtetrabz_initialize(bvec)
  call libtetrabz_kgrid()
  !
  if(any(nge(1:3) /= ngw(1:3))) then
     !
     allocate(wght1(nn, nk0))
     call libtetrabz_occstep1(eig1,eig2,wght1)
     !
     call libtetrabz_interpol_weight(nn,ngw,nge,wght0,wght1)
     deallocate(wght1)
     !
  else
     call libtetrabz_occstep1(eig1,eig2,wght0)
  end if
  !
  deallocate(indx1, indx2, indx3)
  !
end subroutine libtetrabz_occstep
!
! Compute Static polalization function
!
subroutine libtetrabz_polstat(ltetra0,bvec,nb0,nge,eig1,eig2,ngw,wght0)
  !
  use libtetrabz_vals, only : ltetra, ng, nb, nk0, indx1, indx2, indx3
  use libtetrabz_routines, only : libtetrabz_initialize, libtetrabz_interpol_weight, &
  &                               libtetrabz_polstat1
  !
  integer,intent(in) :: ltetra0, nge(3), ngw(3), nb0
  real(8),intent(in) :: bvec(3,3)
  real(8),intent(in) :: eig1(nb0,product(nge(1:3))), eig2(nb0,product(nge(1:3)))
  real(8),intent(out) :: wght0(nb0,nb0,product(ngw(1:3)))
  !
  integer :: ierr, nn
  real(8),allocatable :: wght1(:,:)
  !
  ltetra = ltetra0
  nb = nb0
  ng(1:3) = nge(1:3)
  nn = nb * nb
  !
  call libtetrabz_initialize(bvec)
  call libtetrabz_kgrid()
  !
  if(any(nge(1:3) /= ngw(1:3))) then
     !
     allocate(wght1(nn, nk0))
     call libtetrabz_polstat1(eig1,eig2,wght1)
     !
     call libtetrabz_interpol_weight(nn,ngw,nge,wght0,wght1)
     deallocate(wght1)
     !
  else
     call libtetrabz_polstat1(eig1,eig2,wght0)
  end if
  !
  deallocate(indx1, indx2, indx3)
  !
end subroutine libtetrabz_polstat
!
! Compute Fermi's goldn rule
!
subroutine libtetrabz_fermigr(ltetra0,bvec,nb0,nge,eig1,eig2,ngw,wght0,ne0,e0)
  !
  use libtetrabz_vals, only : ltetra, ng, nb, nk0, indx1, indx2, indx3, ne
  use libtetrabz_routines, only : libtetrabz_initialize, libtetrabz_interpol_weight, &
  &                               libtetrabz_fermigr1
  !
  integer,intent(in) :: ltetra0, nge(3), ngw(3), nb0, ne0
  real(8),intent(in) :: bvec(3,3), e0(ne0)
  real(8),intent(in) :: eig1(nb0,product(nge(1:3))), eig2(nb0,product(nge(1:3)))
  real(8),intent(out) :: wght0(ne0,nb0,nb0,product(ngw(1:3)))
  !
  integer :: ierr, nn
  real(8),allocatable :: wght1(:,:)
  !
  ltetra = ltetra0
  nb = nb0
  ng(1:3) = nge(1:3)
  ne = ne0
  nn = ne * nb * nb
  !
  call libtetrabz_initialize(bvec)
  call libtetrabz_kgrid()
  !
  if(any(nge(1:3) /= ngw(1:3))) then
     !
     allocate(wght1(nn, nk0))
     call libtetrabz_fermigr1(eig1,eig2,e0,wght1)
     !
     call libtetrabz_interpol_weight(nn,ngw,nge,wght0,wght1)
     deallocate(wght1)
     !
  else
     call libtetrabz_fermigr1(eig1,eig2,e0,wght0)
  end if
  !
  deallocate(indx1, indx2, indx3)
  !
end subroutine libtetrabz_fermigr
!
! Compute Polarization of imaginary frequency
!
subroutine libtetrabz_polimg(ltetra0,bvec,nb0,nge,eig1,eig2,ngw,wght0,ne0,e0)
  !
  use libtetrabz_vals, only : ltetra, ng, nb, nk0, indx1, indx2, indx3, ne
  use libtetrabz_routines, only : libtetrabz_initialize, libtetrabz_interpol_weight, &
  &                               libtetrabz_polimg1
  !
  integer,intent(in) :: ltetra0, nge(3), ngw(3), nb0, ne0
  real(8),intent(in) :: bvec(3,3), e0(ne0)
  real(8),intent(in) :: eig1(nb0,product(nge(1:3))), eig2(nb0,product(nge(1:3)))
  real(8),intent(out) :: wght0(2,ne0,nb0,nb0,product(ngw(1:3)))
  !
  integer :: ierr, nn
  real(8),allocatable :: wght1(:,:)
  !
  ltetra = ltetra0
  nb = nb0
  ng(1:3) = nge(1:3)
  ne = ne0
  nn = 2 * ne * nb * nb
  !
  call libtetrabz_initialize(bvec)
  call libtetrabz_kgrid()
  !
  if(any(nge(1:3) /= ngw(1:3))) then
     !
     allocate(wght1(nn, nk0))
     call libtetrabz_polimg1(eig1,eig2,e0,wght1)
     !
     call libtetrabz_interpol_weight(nn,ngw,nge,wght0,wght1)
     deallocate(wght1)
     !
  else
     call libtetrabz_polimg1(eig1,eig2,e0,wght0)
  end if
  !
  deallocate(indx1, indx2, indx3)
  !
end subroutine libtetrabz_polimg
!
! Initialize grid
!
subroutine libtetrabz_kgrid()
  !
  use libtetrabz_vals, only : nk, nk0, indx1, indx2, indx3, ng, ivvec, fst, lst
  !
  implicit none
  !
  integer :: it, i1, i2, i3, ii, ikv(3), nt, ik
  !
  allocate(indx1(20, 6 * nk), indx2(20, 6 * nk), indx3(20 * 6 * nk))
  !
  nt = 0
  do i3 = 1, ng(3)
     do i2  = 1, ng(2)
        do i1 = 1, ng(1)
           !
           do it = 1, 6
              !
              nt = nt + 1
              !
              do ii = 1, 20
                 !
                 ikv(1:3) = (/i1, i2, i3/) + ivvec(1:3,ii,it) - 1
                 ikv(1:3) = modulo(ikv(1:3), ng(1:3))
                 !
                 indx1(ii,nt) = 1 + ikv(1) + ng(1) * ikv(2) + ng(1) * ng(2) * ikv(3)
                 !
              end do
              !
           end do
           !
        end do
     end do
  end do
  !
  indx2(1:20,1:6 * nk) = indx1(1:20,1:6 * nk)
  indx3(1:20 * 6 * nk) = 0
  !
  nk0 = nk
  fst = 1
  lst = 6 * nk
  !
  do ik = 1, nk
     indx3(ik) = ik
  end do
  !
end subroutine libtetrabz_kgrid
!
!
!
end module libtetrabz
