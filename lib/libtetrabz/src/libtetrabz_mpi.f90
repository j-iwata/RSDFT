module libtetrabz_mpi_routines
  !
  implicit none
  !
  integer,save :: comm
  !
contains
!
! Initialize grid
!
subroutine libtetrabz_mpi_kgrid()
  !
  use libtetrabz_vals, only : nk, nk0, indx1, indx2, indx3, ng, ivvec, fst, lst
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
  indx2(1:20,1:6 * nk) = 0
  indx3(1:20 * 6 * nk) = 0
  !
  call libtetrabz_fst_and_lst()
  !
  nk0 = 0
  do it = fst, lst
     !
     do ii = 1, 20
        !
        do ik = 1, nk0
           !
           if(indx1(ii,it) == indx3(ik)) then
              !
              indx2(ii,it) = ik
              goto 10
              !
           end if
           !
        end do
        !
        nk0 = nk0 + 1
        indx2(ii,it) = nk0
        indx3(nk0) = indx1(ii,it)
        !
10      continue
        !
     end do
     !
  end do
  !
end subroutine libtetrabz_mpi_kgrid
!
! Compute cnt and dsp
!
subroutine libtetrabz_fst_and_lst()
  !
  use mpi, only : mpi_comm_size, mpi_comm_rank
  use libtetrabz_vals, only : fst, lst, nk
  !
  integer :: ii, petot, my_rank, ierr
  integer,allocatable :: cnt(:), dsp(:)
  !
  call MPI_COMM_SIZE(comm, petot, ierr)
  call MPI_COMM_RANK(comm, my_rank, ierr)
  !
  allocate(cnt(0:petot-1), dsp(0:petot-1))
  !
  cnt(0:petot-1)        = 6 * nk / petot
  cnt(0:mod(6 * nk,petot)-1) = 6 * nk / petot + 1
  dsp(0) = 0
  do ii = 1, petot - 1
     dsp(ii) = dsp(ii - 1) + cnt(ii - 1)
  end do
  !
  fst = dsp(my_rank) + 1
  lst = dsp(my_rank) + cnt(my_rank)
  !
  deallocate(cnt,dsp)
  !
end subroutine libtetrabz_fst_and_lst
!
end module libtetrabz_mpi_routines
!
!
!
module libtetrabz_mpi
  !
  implicit none
  !
contains
!
! Compute occupation
!
subroutine libtetrabz_mpi_occ(ltetra0,comm0,bvec,nb0,nge,eig,ngw,wght0)
  !
  use mpi, only : MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
  use libtetrabz_vals, only : ltetra, ng, nb, nk0, indx1, indx2, indx3
  use libtetrabz_routines, only : libtetrabz_initialize, libtetrabz_interpol_weight, &
  &                               libtetrabz_occ1
  !
  use libtetrabz_mpi_routines, only : comm, libtetrabz_mpi_kgrid
  !
  integer,intent(in) :: ltetra0, comm0, nge(3), ngw(3), nb0
  real(8),intent(in) :: bvec(3,3)
  real(8),intent(in) :: eig(nb0,product(nge(1:3)))
  real(8),intent(out) :: wght0(nb0,product(ngw(1:3)))
  !
  integer :: ierr, nn
  real(8),allocatable :: wght1(:,:)
  !
  ltetra = ltetra0
  comm = comm0
  nb = nb0
  ng(1:3) = nge(1:3)
  nn = nb
  !
  call libtetrabz_initialize(bvec)
  call libtetrabz_mpi_kgrid()
  !
  allocate(wght1(nn, nk0))
  call libtetrabz_occ1(0d0,eig,wght1)
  !
  call libtetrabz_interpol_weight(nn,ngw,nge,wght0,wght1)
  !
  deallocate(wght1, indx1, indx2, indx3)
  !
  call MPI_allREDUCE(MPI_IN_PLACE, wght0, nn * product(ngw(1:3)), &
  &                  MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
  !
end subroutine libtetrabz_mpi_occ
!
! Calculate Fermi energy
!
subroutine libtetrabz_mpi_fermieng(ltetra0,comm0,bvec,nb0,nge,eig,ngw,wght0,ef,nelec)
  !
  use mpi, only : MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
  use libtetrabz_vals, only : ltetra, ng, nb, nk, nk0, indx1, indx2, indx3
  use libtetrabz_routines, only : libtetrabz_initialize, libtetrabz_interpol_weight, &
  &                           libtetrabz_occ1
  !
  use libtetrabz_mpi_routines, only : comm, libtetrabz_mpi_kgrid
  !
  integer,intent(in) :: ltetra0, comm0, nge(3), ngw(3), nb0
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
  comm = comm0
  nb = nb0
  ng(1:3) = nge(1:3)
  nn = nb
  !
  call libtetrabz_initialize(bvec)
  call libtetrabz_mpi_kgrid()
  !
  allocate(wght1(nn, nk0))
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
     call libtetrabz_occ1(ef, eig,wght1)
     !
     sumkmid = sum(wght1(1:nb,1:nk0))
     call MPI_allREDUCE(MPI_IN_PLACE, sumkmid, 1, &
     &                  MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
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
  call libtetrabz_interpol_weight(nn,ngw,nge,wght0,wght1)
  !
  deallocate(wght1, indx1, indx2, indx3)
  !
  call MPI_allREDUCE(MPI_IN_PLACE, wght0, nn * product(ngw(1:3)), &
  &                  MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
  !
end subroutine libtetrabz_mpi_fermieng
!
! Compute DOS
!
subroutine libtetrabz_mpi_dos(ltetra0,comm0,bvec,nb0,nge,eig,ngw,wght0,ne0,e0)
  !
  use mpi, only : MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
  use libtetrabz_vals, only : ltetra, ng, nb, nk0, indx1, indx2, indx3, ne
  use libtetrabz_routines, only : libtetrabz_initialize, libtetrabz_interpol_weight, &
  &                           libtetrabz_dos1
  !
  use libtetrabz_mpi_routines, only : comm, libtetrabz_mpi_kgrid
  !
  integer,intent(in) :: ltetra0, comm0, nge(3), ngw(3), nb0, ne0
  real(8),intent(in) :: bvec(3,3), e0(ne0)
  real(8),intent(in) :: eig(nb0,product(nge(1:3)))
  real(8),intent(out) :: wght0(ne0,nb0,product(ngw(1:3)))
  !
  integer :: ierr, nn
  real(8),allocatable :: wght1(:,:)
  !
  ltetra = ltetra0
  comm = comm0
  nb = nb0
  ng(1:3) = nge(1:3)
  ne = ne0
  nn = ne * nb
  !
  call libtetrabz_initialize(bvec)
  call libtetrabz_mpi_kgrid()
  !
  allocate(wght1(nn, nk0))
  call libtetrabz_dos1(eig,e0,wght1)
  !
  call libtetrabz_interpol_weight(nn,ngw,nge,wght0,wght1)
  !
  deallocate(wght1, indx1, indx2, indx3)
  !
  call MPI_allREDUCE(MPI_IN_PLACE, wght0, nn * product(ngw(1:3)), &
  &                  MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
  !
end subroutine libtetrabz_mpi_dos
!
! Compute doubledelta
!
subroutine libtetrabz_mpi_doubledelta(ltetra0,comm0,bvec,nb0,nge,eig1,eig2,ngw,wght0)
  !
  use mpi, only : MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
  use libtetrabz_vals, only : ltetra, ng, nb, nk0, indx1, indx2, indx3
  use libtetrabz_routines, only : libtetrabz_initialize, libtetrabz_interpol_weight, &
  &                           libtetrabz_doubledelta1
  !
  use libtetrabz_mpi_routines, only : comm, libtetrabz_mpi_kgrid
  !
  integer,intent(in) :: ltetra0, comm0, nge(3), ngw(3), nb0
  real(8),intent(in) :: bvec(3,3)
  real(8),intent(in) :: eig1(nb0,product(nge(1:3))), eig2(nb0,product(nge(1:3)))
  real(8),intent(out) :: wght0(nb0,nb0,product(ngw(1:3)))
  !
  integer :: ierr, nn
  real(8),allocatable :: wght1(:,:)
  !
  ltetra = ltetra0
  comm = comm0
  nb = nb0
  ng(1:3) = nge(1:3)
  nn = nb * nb
  !
  call libtetrabz_initialize(bvec)
  call libtetrabz_mpi_kgrid()
  !
  allocate(wght1(nn, nk0))
  call libtetrabz_doubledelta1(eig1,eig2,wght1)
  !
  call libtetrabz_interpol_weight(nn,ngw,nge,wght0,wght1)
  !
  deallocate(wght1, indx1, indx2, indx3)
  !
  call MPI_allREDUCE(MPI_IN_PLACE, wght0, nn * product(ngw(1:3)), &
  &                  MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
  !
end subroutine libtetrabz_mpi_doubledelta
!
! Compute Occ * Step
!
subroutine libtetrabz_mpi_occstep(ltetra0,comm0,bvec,nb0,nge,eig1,eig2,ngw,wght0)
  !
  use mpi, only : MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
  use libtetrabz_vals, only : ltetra, ng, nb, nk0, indx1, indx2, indx3
  use libtetrabz_routines, only : libtetrabz_initialize, libtetrabz_interpol_weight, &
  &                           libtetrabz_occstep1
  !
  use libtetrabz_mpi_routines, only : comm, libtetrabz_mpi_kgrid
  !
  integer,intent(in) :: ltetra0, comm0, nge(3), ngw(3), nb0
  real(8),intent(in) :: bvec(3,3)
  real(8),intent(in) :: eig1(nb0,product(nge(1:3))), eig2(nb0,product(nge(1:3)))
  real(8),intent(out) :: wght0(nb0,nb0,product(ngw(1:3)))
  !
  integer :: ierr, nn
  real(8),allocatable :: wght1(:,:)
  !
  ltetra = ltetra0
  comm = comm0
  nb = nb0
  ng(1:3) = nge(1:3)
  !
  call libtetrabz_initialize(bvec)
  call libtetrabz_mpi_kgrid()
  nn = nb * nb
  !
  allocate(wght1(nn, nk0))
  call libtetrabz_occstep1(eig1,eig2,wght1)
  !
  call libtetrabz_interpol_weight(nn,ngw,nge,wght0,wght1)
  !
  deallocate(wght1, indx1, indx2, indx3)
  !
  call MPI_allREDUCE(MPI_IN_PLACE, wght0, nn * product(ngw(1:3)), &
  &                  MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
  !
end subroutine libtetrabz_mpi_occstep
!
! Compute Static polalization function
!
subroutine libtetrabz_mpi_polstat(ltetra0,comm0,bvec,nb0,nge,eig1,eig2,ngw,wght0)
  !
  use mpi, only : MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
  use libtetrabz_vals, only : ltetra, ng, nb, nk0, indx1, indx2, indx3
  use libtetrabz_routines, only : libtetrabz_initialize, libtetrabz_interpol_weight, &
  &                           libtetrabz_polstat1
  !
  use libtetrabz_mpi_routines, only : comm, libtetrabz_mpi_kgrid
  !
  integer,intent(in) :: ltetra0, comm0, nge(3), ngw(3), nb0
  real(8),intent(in) :: bvec(3,3)
  real(8),intent(in) :: eig1(nb0,product(nge(1:3))), eig2(nb0,product(nge(1:3)))
  real(8),intent(out) :: wght0(nb0,nb0,product(ngw(1:3)))
  !
  integer :: ierr, nn
  real(8),allocatable :: wght1(:,:)
  !
  ltetra = ltetra0
  comm = comm0
  nb = nb0
  ng(1:3) = nge(1:3)
  nn = nb * nb
  !
  call libtetrabz_initialize(bvec)
  call libtetrabz_mpi_kgrid()
  !
  allocate(wght1(nn, nk0))
  call libtetrabz_polstat1(eig1,eig2,wght1)
  !
  call libtetrabz_interpol_weight(nn,ngw,nge,wght0,wght1)
  !
  deallocate(wght1, indx1, indx2, indx3)
  !
  call MPI_allREDUCE(MPI_IN_PLACE, wght0, nn * product(ngw(1:3)), &
  &                  MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
  !
end subroutine libtetrabz_mpi_polstat
!
! Compute Fermi's goldn rule
!
subroutine libtetrabz_mpi_fermigr(ltetra0,comm0,bvec,nb0,nge,eig1,eig2,ngw,wght0,ne0,e0)
  !
  use mpi, only : MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
  use libtetrabz_vals, only : ltetra, ng, nb, nk0, ne, indx1, indx2, indx3
  use libtetrabz_routines, only : libtetrabz_initialize, libtetrabz_interpol_weight, &
  &                           libtetrabz_fermigr1
  !
  use libtetrabz_mpi_routines, only : comm, libtetrabz_mpi_kgrid
  !
  integer,intent(in) :: ltetra0, comm0, nge(3), ngw(3), nb0, ne0
  real(8),intent(in) :: bvec(3,3), e0(ne0)
  real(8),intent(in) :: eig1(nb0,product(nge(1:3))), eig2(nb0,product(nge(1:3)))
  real(8),intent(out) :: wght0(ne0,nb0,nb0,product(ngw(1:3)))
  !
  integer :: ierr, nn
  real(8),allocatable :: wght1(:,:)
  !
  ltetra = ltetra0
  comm = comm0
  nb = nb0
  ng(1:3) = nge(1:3)
  ne = ne0
  nn = ne * nb * nb
  !
  call libtetrabz_initialize(bvec)
  call libtetrabz_mpi_kgrid()
  !
  allocate(wght1(nn, nk0))
  call libtetrabz_fermigr1(eig1,eig2,e0,wght1)
  !
  call libtetrabz_interpol_weight(nn,ngw,nge,wght0,wght1)
  !
  deallocate(wght1, indx1, indx2, indx3)
  !
  call MPI_allREDUCE(MPI_IN_PLACE, wght0, nn * product(ngw(1:3)), &
  &                  MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
  !
end subroutine libtetrabz_mpi_fermigr
!
! Compute Polarization of imaginary frequency
!
subroutine libtetrabz_mpi_polimg(ltetra0,comm0,bvec,nb0,nge,eig1,eig2,ngw,wght0,ne0,e0)
  !
  use mpi, only : MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
  use libtetrabz_vals, only : ltetra, ng, nb, nk0, indx1, indx2, indx3, ne
  use libtetrabz_routines, only : libtetrabz_initialize, libtetrabz_interpol_weight, &
  &                               libtetrabz_polimg1
  !
  use libtetrabz_mpi_routines, only : comm, libtetrabz_mpi_kgrid
  !
  integer,intent(in) :: ltetra0, comm0, nge(3), ngw(3), nb0, ne0
  real(8),intent(in) :: bvec(3,3), e0(ne0)
  real(8),intent(in) :: eig1(nb0,product(nge(1:3))), eig2(nb0,product(nge(1:3)))
  real(8),intent(out) :: wght0(2,ne0,nb0,nb0,product(ngw(1:3)))
  !
  integer :: ierr, nn
  real(8),allocatable :: wght1(:,:)
  !
  ltetra = ltetra0
  comm = comm0
  nb = nb0
  ng(1:3) = nge(1:3)
  ne = ne0
  nn = 2 * ne * nb * nb
  !
  call libtetrabz_initialize(bvec)
  call libtetrabz_mpi_kgrid()
  !
  allocate(wght1(nn, nk0))
  call libtetrabz_polimg1(eig1,eig2,e0,wght1)
  !
  call libtetrabz_interpol_weight(nn,ngw,nge,wght0,wght1)
  !
  deallocate(wght1, indx1, indx2, indx3)
  !
  call MPI_allREDUCE(MPI_IN_PLACE, wght0, nn * product(ngw(1:3)), &
  &                  MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
  !
end subroutine libtetrabz_mpi_polimg
!
!
!
end module libtetrabz_mpi
