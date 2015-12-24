MODULE ps_prepNzqr_g_module

  use parallel_module, only: myrank
  use var_ps_member_g, only: nlop_pair, nzqr_pair,qij_f, ddi, k1map, atommap &
                         ,N_nzqr, N_nlop, Dij00, Dij, qqr, qij &
                         ,N_k1, k1max, k1_to_l, k1_to_m, k1_to_iorb
  use var_ps_member, only: no
  use atom_module, only: Natom,ki_atom
  use ps_nloc2_variables, only: nzlma, amap, lmap, mmap, iorbmap
  use electron_module, only: nspin

  implicit none
  
  PRIVATE
  PUBLIC :: prepNzqr

CONTAINS

  SUBROUTINE prepNzqr

    implicit none
    integer :: kk1,k1,iorb,lma
    integer :: i,j,ik
    integer :: lma1,lma2,a1,a2,l1,l2,m1,m2,i1,i2,a,l,m
    integer,allocatable :: k1a(:)

    call write_border( 0," prepNzqr(start)" )

!----- get N_nzqr -----

    kk1=0
    do lma1=1,nzlma
       if ( amap(lma1) == 0 .or. iorbmap(lma1) == 0 ) cycle
       a1=amap(lma1)
       l1=lmap(lma1)
       m1=mmap(lma1)
       i1=no(iorbmap(lma1),ki_atom(a1))
       do lma2=1,nzlma
          if ( amap(lma2) == 0 .or. iorbmap(lma2) == 0 ) cycle
          a2=amap(lma2)
          l2=lmap(lma2)
          m2=mmap(lma2)
          i2=no(iorbmap(lma2),ki_atom(a2))
          if ( a2 /= a1 ) cycle
          if ( l2 > l1 ) cycle
          if ( l1 == l2 .and. i2 > i1 ) cycle
          if ( l1 == l2 .and. i1 == i2 .and. m2 > m1 ) cycle
          kk1=kk1+1
       end do
    end do

    N_nzqr = kk1

!    if ( myrank == 0 ) write(*,*) "N_nzqr= ",N_nzqr

!===== get N_nzqr =====

!----- get N_nlop -----

    kk1=0
    do lma1=1,nzlma
       if ( amap(lma1) == 0 .or. iorbmap(lma1) == 0 ) cycle
       do lma2=1,nzlma
          if ( amap(lma2) == 0 .or. iorbmap(lma2) == 0 ) cycle
          if ( amap(lma1) /= amap(lma2) .or. &
               lmap(lma1) /= lmap(lma2) .or. &
               mmap(lma1) /= mmap(lma2) ) cycle
          kk1=kk1+1
       end do
    end do

    N_nlop = kk1

!    if ( myrank == 0 ) write(*,*) "N_nlop= ",N_nlop

!===== get N_nlop =====

    if ( allocated(nzqr_pair) ) deallocate( nzqr_pair )
    if ( allocated(k1map)     ) deallocate( k1map     )
    if ( allocated(nlop_pair) ) deallocate( nlop_pair )
    if ( allocated(Dij)       ) deallocate( Dij       )
    if ( allocated(Dij00)     ) deallocate( Dij00     )
    if ( allocated(qij)       ) deallocate( qij       )
    if ( allocated(qij_f)     ) deallocate( qij_f     )
    allocate( nzqr_pair(N_nzqr,2) ) ; nzqr_pair=0
    allocate( k1map(N_nzqr)       ) ; k1map=0
    allocate( nlop_pair(2,N_nlop) ) ; nlop_pair=0
    allocate( Dij(N_nzqr,nspin)   ) ; Dij=0.0d0
    allocate( Dij00(N_nzqr)       ) ; Dij00=0.0d0
    allocate( qij(N_nzqr)         ) ; qij=0.0d0
    allocate( qij_f(N_nzqr)       ) ; qij_f=0.0d0

!----- get nzqr_pair, atommap, k1map, kk1map -----

    kk1=0
    do lma1=1,nzlma
       if ( amap(lma1) == 0 .or. iorbmap(lma1) == 0 ) cycle
       a1=amap(lma1)
       l1=lmap(lma1)
       m1=mmap(lma1)
       i1=no(iorbmap(lma1),ki_atom(a1))
    do lma2=1,nzlma
       if ( amap(lma2) == 0 .or. iorbmap(lma2) == 0 ) cycle
       a2=amap(lma2)
       l2=lmap(lma2)
       m2=mmap(lma2)
       i2=no(iorbmap(lma2),ki_atom(a2))
       if ( a2 /= a1 ) cycle
       if ( l2 > l1 ) cycle
       if ( l1 == l2 .and. i2 > i1 ) cycle
       if ( l1 == l2 .and. i1 == i2 .and. m2 > m1 ) cycle
       kk1=kk1+1
       nzqr_pair(kk1,1) = lma1
       nzqr_pair(kk1,2) = lma2
       ik=ki_atom(a1)
       do k1=1,N_k1(ik)
          if ( k1_to_iorb(1,k1,ik) == iorbmap(lma1) .and. &
               k1_to_iorb(2,k1,ik) == iorbmap(lma2) .and. &
               k1_to_l(1,k1,ik) == l1 .and. &
               k1_to_l(2,k1,ik) == l2 .and. &
               k1_to_m(1,k1,ik) == m1 .and. &
               k1_to_m(2,k1,ik) == m2 ) then
             k1map(kk1) = k1
             exit
          end if
       end do ! k1
    end do ! lma2
    end do ! lma1

!    if ( myrank == 0 ) write(*,*) "N_nzqr,kk1 (check) =",N_nzqr,kk1

!===== get nzqr_pair, atommap, k1map, kk1map =====

!    if ( myrank == 0 ) write(*,*) "--- Dij00(1:N_nzqr) ---"
!    if ( myrank == 0 ) write(*,*) "  [ qij_f(1:N_nzqr) ]  "
    do kk1=1,N_nzqr
       i=nzqr_pair(kk1,1)
       j=nzqr_pair(kk1,2)
       a1=amap(i)
       a2=amap(j)
       l1=lmap(i)
       l2=lmap(j)
       m1=mmap(i)
       m2=mmap(j)
       if ( .not.( a1 == a2 .and. l1 == l2 .and. m1 == m2 ) ) cycle
       ik=ki_atom(a1)
       i1=no(iorbmap(i),ik)
       i2=no(iorbmap(j),ik)
       Dij00(kk1)=ddi(i1,i2,l1+1,ik)
       qij_f(kk1)=qqr(i1,i2,l1+1,ik)
    end do

!----- Nlop_type Matrix -----

!    if ( myrank == 0 ) write(*,*) "--- Dij0 (1:N_nlop) ---"
!    if ( myrank == 0 ) write(*,*) "   [ qij(1:N_nlop) ]   "

    kk1=0
    do lma1=1,nzlma
       if ( amap(lma1) == 0 .or. iorbmap(lma1) == 0 ) cycle
       a1=amap(lma1)
       l1=lmap(lma1)
       m1=mmap(lma1)
       ik=ki_atom(a1)
       do lma2=1,nzlma
          a2=amap(lma2)
          l2=lmap(lma2)
          m2=mmap(lma2)
          if ( a2 == 0 .or. iorbmap(lma2) == 0 ) cycle
          if ( a1 /= a2 .or. l1 /= l2 .or. m1 /= m2 ) cycle
          kk1=kk1+1
          nlop_pair(1,kk1)=lma1
          nlop_pair(2,kk1)=lma2
          i1=no(iorbmap(lma1),ik)
          i2=no(iorbmap(lma2),ik)
          qij(kk1)=qqr(i1,i2,l1+1,ik)
       end do
    end do

!    if ( myrank == 0 ) write(*,*) "N_nlop= ",N_nlop,kk1

!===== Nlop_type Matrix =====

    call write_border( 0," prepNzqr(end)" )

  END SUBROUTINE prepNzqr  

END MODULE ps_prepNzqr_g_module
