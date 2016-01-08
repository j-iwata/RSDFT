MODULE fock_module

  use xc_hybrid_module, only: occ_hf, unk_hf, iflag_hybrid, alpha_hf &
                             ,FKMB_0,FKMB_1,FKBZ_0,FKBZ_1,FOCK_0,FOCK_1 &
                             ,occ_factor,gamma_hf
  use array_bound_module, only: ML_0,ML_1,MB,MB_0,MB_1,MBZ,MBZ_0,MBZ_1 &
                               ,MSP_0,MSP_1
  use wf_module, only: unk, occ, hunk
  use fock_fft_module
  use fock_cg_module
  use parallel_module
  use watch_module
  use rsdft_mpi_module
  use hartree_mol_module, only: timer_reset_hartree_mol, timer_result_hartree_mol

  implicit none

  PRIVATE
  PUBLIC :: Fock, op_fock, UpdateWF_fock

#ifdef _DRSDFT_
  real(8),parameter :: zero = 0.0d0
  integer,parameter :: TYPE_MAIN = MPI_REAL8
#else
  complex(8),parameter :: zero = (0.0d0,0.0d0)
  integer,parameter :: TYPE_MAIN = MPI_COMPLEX16
#endif

  integer :: SYStype=0

CONTAINS


  SUBROUTINE Fock( n,k,s, n1,n2, psi, tpsi )
    implicit none
    integer,intent(IN) :: n,k,s,n1,n2
#ifdef _DRSDFT_
    real(8),intent(IN)  :: psi(n1:n2)
    real(8),intent(INOUT) :: tpsi(n1:n2)
    real(8),allocatable :: trho(:),tvht(:),tphi(:)
    integer,parameter :: TYPE_MAIN = MPI_REAL8
#else
    complex(8),intent(IN)  :: psi(n1:n2)
    complex(8),intent(INOUT) :: tpsi(n1:n2)
    complex(8),allocatable :: trho(:),tvht(:),tphi(:)
    integer,parameter :: TYPE_MAIN = MPI_COMPLEX16
#endif
    real(8) :: c
    integer :: ML0,m,q,i,t,ierr

#ifdef _DRSDFT_
    if ( SYStype == 0 ) then
       call Fock_Double( s, n1, n2, psi, tpsi )
       return
    end if
#endif

    ML0 = n2 - n1 + 1

    allocate( trho(n1:n2) ) ; trho=zero
    allocate( tvht(n1:n2) ) ; tvht=zero
    allocate( tphi(n1:n2) ) ; tphi=zero

! ---

    if ( gamma_hf == 0 ) then

#ifdef _DRSDFT_
       stop "stop@Fock: This routine is only for COMPLEX16 calculations"
#else
       do t=FOCK_0,FOCK_1

          if ( t == 1 ) then

             do q=FKBZ_0,FKBZ_1
             do m=FKMB_0,FKMB_1

                if ( abs(occ_hf(m,q,s)) < 1.d-10 ) cycle

                do i=n1,n2
                   trho(i) = conjg(unk_hf(i,m,q,s))*psi(i)
                end do ! i

                call Fock_fft(n1,n2,k,q,trho,tvht,t)

                c = occ_factor*occ_hf(m,q,s)*alpha_hf

                do i=n1,n2
                   tphi(i)=tphi(i)-c*tvht(i)*unk_hf(i,m,q,s)
                end do

             end do ! m
             end do ! q

          else ! [ t /= 1 ]

             do q=FKBZ_0,FKBZ_1
             do m=FKMB_0,FKMB_1

                if ( abs(occ_hf(m,q,s)) < 1.d-10 ) cycle

                do i=n1,n2
                   trho(i)=unk_hf(i,m,q,s)*psi(i)
                end do

                call Fock_fft(n1,n2,k,q,trho,tvht,t)

                c = occ_factor*occ_hf(m,q,s)*alpha_hf

                do i=n1,n2
                   tphi(i)=tphi(i)-c*tvht(i)*conjg(unk_hf(i,m,q,s))
                end do ! i

             end do ! m
             end do ! q

          end if ! [ t ]

       end do ! t
#endif

    else ! [ gamma_hf /= 0 ]

       q = k

       do m=FKMB_0,FKMB_1

          if ( abs(occ_hf(m,q,s)) < 1.d-10 ) cycle

          do i=n1,n2
#ifdef _DRSDFT_
             trho(i) = unk_hf(i,m,q,s)*psi(i)
#else   
             trho(i) = conjg(unk_hf(i,m,q,s))*psi(i)
#endif
          end do

          if ( SYStype == 1 ) then
             call Fock_cg( n1,n2,k,q,trho,tvht,1 )
          else
             call Fock_fft(n1,n2,k,q,trho,tvht,1)
          end if

          c = alpha_hf*(2.0d0*occ_factor)*occ_hf(m,q,s)

          do i=n1,n2
             tphi(i)=tphi(i)-c*tvht(i)*unk_hf(i,m,q,s)
          end do

       end do ! m

    end if ! [ gamma_hf ]

    call mpi_allreduce(tphi,trho,ML0,TYPE_MAIN,mpi_sum,comm_fkmb,ierr)
    do i=n1,n2
       tpsi(i) = tpsi(i) + trho(i)
    end do

    deallocate( tphi )
    deallocate( tvht ) 
    deallocate( trho )

    return

  END SUBROUTINE Fock


  SUBROUTINE op_fock(k,s,n1,n2,ib1,ib2,tpsi,htpsi)
    implicit none
    integer,intent(IN) :: k,s,n1,n2,ib1,ib2
#ifdef _DRSDFT_
    real(8),intent(IN) :: tpsi(n1:n2,ib1:ib2)
    real(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
#else
    complex(8),intent(IN) :: tpsi(n1:n2,ib1:ib2)
    complex(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
#endif
    integer :: ib,i

    if ( iflag_hybrid == 0 ) return

    if ( iflag_hybrid == 2 ) then

       do ib=ib1,ib2
          do i=n1,n2
             htpsi(i,ib)=htpsi(i,ib)+hunk(i,ib,k,s)
          end do
       end do

    else if ( iflag_hybrid > 0 ) then

       do ib=ib1,ib2
          call Fock( ib, k, s, n1, n2, tpsi(n1,ib), htpsi(n1,ib) )
       end do

    end if

  END SUBROUTINE op_fock


  SUBROUTINE UpdateWF_fock( SYStype_in )
    implicit none
    integer,optional,intent(IN) :: SYStype_in
    integer :: s,k,n,m,i_occ,i_orb,ierr

    if ( present(SYStype_in) ) SYStype = SYStype_in

!    if ( disp_switch_parallel ) write(*,*) "UpdateWF_fock"

! ---

    occ_hf(:,:,:)   = 0.0d0
    unk_hf(:,:,:,:) = zero

!#ifdef _DRSDFT_
!    do s=MSP_0,MSP_1
!    do k=MBZ_0,MBZ_1
!       if ( k < FKBZ_0 .or. FKBZ_1 < k ) cycle
!       i_orb=0
!       i_occ=0
!       do n=1,MB
!          if ( abs(occ(n,k,s)) <= 1.d-10 ) cycle
!          i_occ = i_occ + 1
!          if ( mod(i_occ-1,np_fkmb) == myrank_f ) then 
!             i_orb = i_orb + 1
!             if ( MB_0 <= n .and. n <= MB_1 ) unk_hf(:,i_orb,k,s) = unk(:,n,k,s)
!             occ_hf(i_orb,k,s) = occ(n,k,s)
!          end if
!       end do ! n
!    end do ! k
!    end do ! s
!#else
    do s=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1
    do n=MB_0,MB_1
       unk_hf(:,n,k,s) = unk(:,n,k,s)
       occ_hf(n,k,s) = occ(n,k,s)
    end do ! n
    end do ! k
    end do ! s
!#endif

    m = size(unk_hf,1)*size(unk_hf,2)
    do s=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1
       call MPI_ALLREDUCE( MPI_IN_PLACE, unk_hf(:,:,k,s), m, TYPE_MAIN &
            ,MPI_SUM, comm_band, ierr )
    end do ! k
    end do ! s

    m = size(unk_hf,1)*size(unk_hf,2)*size(unk_hf,3)
    do s=MSP_0,MSP_1
       call MPI_ALLREDUCE( MPI_IN_PLACE, unk_hf(:,:,:,s), m, TYPE_MAIN &
            ,MPI_SUM, comm_bzsm, ierr )
    end do ! s

    m = size( occ_hf,1 )
    do s=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1
       call MPI_ALLREDUCE( MPI_IN_PLACE, occ_hf(:,k,s), m, MPI_REAL8 &
            ,MPI_SUM, comm_band, ierr )
    end do ! k
    end do ! s

    m = size( occ_hf,1 )*size( occ_hf,2 )
    do s=MSP_0,MSP_1
       call MPI_ALLREDUCE( MPI_IN_PLACE, occ_hf(:,:,s), m, MPI_REAL8 &
            ,MPI_SUM, comm_bzsm, ierr )
    end do ! s

! ---

    hunk(:,:,:,:) = zero
    ct_fock_fft(:)=0.0d0
    et_fock_fft(:)=0.0d0

    do s=MSP_0,MSP_1
       if ( gamma_hf == 1 ) then
          call Fock_4( s,ML_0,ML_1 )
       else
          call Fock_5( s,ML_0,ML_1 )
       end if
    end do ! s

    iflag_hybrid = 2

!    if ( disp_switch_parallel ) then
!       write(*,'(1x,"iflag_hybrid=",i2)') iflag_hybrid
!       write(*,*) "time(fock_fft1)=",ct_fock_fft(1),et_fock_fft(1)
!       write(*,*) "time(fock_fft2)=",ct_fock_fft(2),et_fock_fft(2)
!       write(*,*) "time(fock_fft3)=",ct_fock_fft(3),et_fock_fft(3)
!       write(*,*) "time(fock_fft4)=",ct_fock_fft(4),et_fock_fft(4)
!       write(*,*) "time(fock_fft5)=",ct_fock_fft(5),et_fock_fft(5)
!       write(*,*) "time(fock_fft6)=",ct_fock_fft(6),et_fock_fft(6)
!       write(*,*) "time(fock_fft7)=",ct_fock_fft(7),et_fock_fft(7)
!       write(*,*) "time(fock_fft8)=",ct_fock_fft(8),et_fock_fft(8)
!    end if

  END SUBROUTINE UpdateWF_fock


  SUBROUTINE Fock_4( s,n1,n2 )
    implicit none
    integer,intent(IN) :: s,n1,n2
#ifdef _DRSDFT_
    real(8),allocatable :: trho(:),tvht(:)
#else
    complex(8),allocatable :: trho(:),tvht(:)
#endif
    real(8) :: c
    integer :: m,n,i,j,k,nwork,iwork,ierr,nwork_0,nwork_1
    integer,allocatable :: mapwork(:,:)

#ifdef _DRSDFT_
    if ( SYStype == 0 ) then
       call Fock_4_Double( s,n1,n2 )
       return
    end if
#endif

    k = 1

! ---

    nwork=0
    i=0
    do n=1,MB
       do m=1,n
          if ( abs(occ(n,k,s))<1.d-10 .and. abs(occ(m,k,s))<1.d-10 ) cycle
          i=i+1
          j=mod(i-1,np_band)
          if ( j == myrank_b ) nwork=nwork+1
       end do
    end do

!    if ( disp_switch_parallel ) then
!       write(*,*) "total # of me    =",i
!       write(*,*) "# of me of rank0 =",nwork
!    end if

    allocate( mapwork(2,nwork) ) ; mapwork=0

    nwork=0
    i=0
    do n=1,MB
       do m=1,n
          if ( abs(occ(n,k,s))<1.d-10 .and. abs(occ(m,k,s))<1.d-10 ) cycle
          i=i+1
          j=mod(i-1,np_band)
          if ( j == myrank_b ) then
             nwork=nwork+1
             mapwork(1,nwork)=m
             mapwork(2,nwork)=n
          end if
       end do
    end do

    ir_fkmb(:)=0
    id_fkmb(:)=0
    do i=1,nwork
       j=mod(i-1,np_fkmb)
       ir_fkmb(j)=ir_fkmb(j)+1
    end do
    do j=0,np_fkmb-1
       id_fkmb(j)=sum(ir_fkmb(0:j))-ir_fkmb(j)
    end do

    nwork_0 = id_fkmb(myrank_f) + 1
    nwork_1 = id_fkmb(myrank_f) + ir_fkmb(myrank_f)

! ---

    allocate( trho(n1:n2) ) ; trho=zero
    allocate( tvht(n1:n2) ) ; tvht=zero

! ---

    do iwork=nwork_0,nwork_1

       m = mapwork(1,iwork)
       n = mapwork(2,iwork)

       do i=n1,n2
#ifdef _DRSDFT_
          trho(i) = unk(i,m,k,s)*unk(i,n,k,s)
#else   
          trho(i) = conjg(unk(i,m,k,s))*unk(i,n,k,s)
#endif
       end do

       if ( SYStype == 1 ) then
          call Fock_cg( n1,n2,k,k,trho,tvht,1 )
       else
          call Fock_fft( n1,n2,k,k,trho,tvht,1 )
       end if

       if ( abs(occ(m,k,s)) >= 1.d-10 ) then
          c = alpha_hf*(2.0d0*occ_factor)*occ(m,k,s)
          do i=n1,n2
             hunk(i,n,k,s)=hunk(i,n,k,s)-c*tvht(i)*unk(i,m,k,s)
          end do
       end if

       if ( m == n ) cycle

       if ( abs(occ(n,k,s)) >= 1.d-10 ) then
          c = alpha_hf*(2.0d0*occ_factor)*occ(n,k,s)
          do i=n1,n2
#ifdef _DRSDFT_
             hunk(i,m,k,s)=hunk(i,m,k,s)-c*tvht(i)*unk(i,n,k,s)
#else
             hunk(i,m,k,s)=hunk(i,m,k,s)-c*conjg(tvht(i))*unk(i,n,k,s)
#endif
          end do
       end if

    end do ! iwork

    deallocate( tvht ) 
    deallocate( trho )

    deallocate( mapwork )

! ---

    call mpi_allreduce( MPI_IN_PLACE,hunk(n1,1,k,s),MB*(n2-n1+1) &
                       ,TYPE_MAIN,MPI_SUM,comm_band,ierr )
    call mpi_allreduce( MPI_IN_PLACE,hunk(n1,1,k,s),MB*(n2-n1+1) &
                       ,TYPE_MAIN,MPI_SUM,comm_fkmb,ierr )
!    call rsdft_allreduce( comm_band, hunk(n1,1,k,s), MB*(n2-n1+1), 512 )

    return

  END SUBROUTINE Fock_4

#ifdef _DRSDFT_
  SUBROUTINE Fock_4_Double( s, ml0, ml1 )
    implicit none
    integer,intent(IN) :: s,ml0,ml1
    complex(8),allocatable :: trho(:),tvht(:)
    real(8),parameter :: tol=1.d-10
    real(8) :: c
    integer :: m1,n1,m2,n2,m,n,i,j,k,nwork,iwork1,iwork2,ierr
    integer :: nwork_0,nwork_1
    integer,allocatable :: mapwork(:,:)

    k = 1

! ---

    nwork=0
    i=0
    do n=1,MB
       do m=1,n
          if ( abs(occ(n,k,s)) < tol .and. abs(occ(m,k,s)) < tol ) cycle
          i=i+1
          j=mod(i-1,np_band)
          if ( j == myrank_b ) nwork=nwork+1
       end do
    end do

!    if ( disp_switch_parallel ) then
!       write(*,*) "total # of me    =",i
!       write(*,*) "# of me of rank0 =",nwork
!    end if

    allocate( mapwork(2,nwork) ) ; mapwork=0

    nwork=0
    i=0
    do n=1,MB
       do m=1,n
          if ( abs(occ(n,k,s)) < tol .and. abs(occ(m,k,s)) < tol ) cycle
          i=i+1
          j=mod(i-1,np_band)
          if ( j == myrank_b ) then
             nwork=nwork+1
             mapwork(1,nwork)=m
             mapwork(2,nwork)=n
          end if
       end do
    end do

    ir_fkmb(:)=0
    id_fkmb(:)=0
    do i=1,nwork
       j=mod(i-1,np_fkmb)
       ir_fkmb(j)=ir_fkmb(j)+1
    end do
    do j=0,np_fkmb-1
       id_fkmb(j)=sum(ir_fkmb(0:j))-ir_fkmb(j)
    end do

    nwork_0 = id_fkmb(myrank_f) + 1
    nwork_1 = id_fkmb(myrank_f) + ir_fkmb(myrank_f)

! ---

    allocate( trho(ml0:ml1) ) ; trho=(0.0d0,0.0d0)
    allocate( tvht(ml0:ml1) ) ; tvht=(0.0d0,0.0d0)

! ---

    do iwork1=nwork_0,nwork_1,2

       iwork2=min(iwork1+1,nwork_1)

       m1 = mapwork(1,iwork1)
       n1 = mapwork(2,iwork1)
       m2 = mapwork(1,iwork2)
       n2 = mapwork(2,iwork2)

       if ( iwork1 < iwork2 ) then
          do i=ml0,ml1
             trho(i) = dcmplx( unk(i,m1,k,s)*unk(i,n1,k,s) &
                              ,unk(i,m2,k,s)*unk(i,n2,k,s) )
          end do
       else
          do i=ml0,ml1
             trho(i) = unk(i,m1,k,s)*unk(i,n1,k,s)
          end do
       end if

       call Fock_FFT_Double( ml0, ml1, trho, tvht )

       if ( abs(occ(m1,k,s)) >= tol ) then
          c = alpha_hf*(2.0d0*occ_factor)*occ(m1,k,s)
          do i=ml0,ml1
             hunk(i,n1,k,s)=hunk(i,n1,k,s)-c*real( tvht(i) )*unk(i,m1,k,s)
          end do
       end if

       if ( m1 /= n1 .and. abs(occ(n1,k,s)) >= tol ) then
          c = alpha_hf*(2.0d0*occ_factor)*occ(n1,k,s)
          do i=ml0,ml1
             hunk(i,m1,k,s)=hunk(i,m1,k,s)-c*real( tvht(i) )*unk(i,n1,k,s)
          end do
       end if

       if ( iwork1 < iwork2 ) then

          if ( abs(occ(m2,k,s)) >= tol ) then
             c = alpha_hf*(2.0d0*occ_factor)*occ(m2,k,s)
             do i=ml0,ml1
                hunk(i,n2,k,s)=hunk(i,n2,k,s)-c*aimag( tvht(i) )*unk(i,m2,k,s)
             end do
          end if

          if ( m2 /= n2 .and. abs(occ(n2,k,s)) >= tol ) then
             c = alpha_hf*(2.0d0*occ_factor)*occ(n2,k,s)
             do i=ml0,ml1
                hunk(i,m2,k,s)=hunk(i,m2,k,s)-c*aimag( tvht(i) )*unk(i,n2,k,s)
             end do
          end if

       end if

    end do ! iwork

    deallocate( tvht ) 
    deallocate( trho )

    deallocate( mapwork )

! ---

    call mpi_allreduce( MPI_IN_PLACE,hunk(ml0,1,k,s),MB*(ml1-ml0+1) &
                       ,TYPE_MAIN,MPI_SUM,comm_band,ierr )
    call mpi_allreduce( MPI_IN_PLACE,hunk(ml0,1,k,s),MB*(ml1-ml0+1) &
                       ,TYPE_MAIN,MPI_SUM,comm_fkmb,ierr )
!    call rsdft_allreduce( comm_band, hunk(ml0,1,k,s), MB*(ml1-ml0+1), 512 )

    return

  END SUBROUTINE Fock_4_Double
#endif

  SUBROUTINE Fock_5( s,n1,n2 )
    implicit none
    integer,intent(IN) :: s,n1,n2
#ifdef _DRSDFT_
    stop "stop@Fock_5: This routine is only for COMPLEX16 calculations"
#else
    complex(8),allocatable :: trho(:),tvht(:)
    real(8) :: c,ctt(0:3),ett(0:3)
    integer :: m,n,q,k,i,j,a,b,nwork,iwork,ierr,nwork_0,nwork_1
    integer,allocatable :: mapnk(:,:),mapwork(:,:)

    call watch(ctt(0),ett(0))

! ---

    allocate( mapnk(2,MB*MBZ) ) ; mapnk=0

    i=0
    do k=1,MBZ
       do n=1,MB
          i=i+1
          mapnk(1,i)=n
          mapnk(2,i)=k
       end do
    end do

! ---

    nwork=0
    a=0
    do j=1,MB*MBZ
       do i=1,j
          m=mapnk(1,i)
          q=mapnk(2,i)
          n=mapnk(1,j)
          k=mapnk(2,j)
          if ( abs(occ_hf(n,k,s))<1.d-10 .and. abs(occ_hf(m,q,s))<1.d-10 ) cycle
          a=a+1
          b=mod(a-1,np_band*np_bzsm)
          if ( b == myrank_k+myrank_b*np_bzsm ) nwork=nwork+1
       end do
    end do

!    if ( disp_switch_parallel ) then
!       write(*,*) "total # of me    =",a
!       write(*,*) "# of me of rank0 =",nwork
!    end if

    allocate( mapwork(2,nwork) ) ; mapwork=0

    nwork=0
    a=0
    do j=1,MB*MBZ
       do i=1,j
          m=mapnk(1,i)
          q=mapnk(2,i)
          n=mapnk(1,j)
          k=mapnk(2,j)
          if ( abs(occ_hf(n,k,s))<1.d-10 .and. abs(occ_hf(m,q,s))<1.d-10 ) cycle
          a=a+1
          b=mod(a-1,np_band*np_bzsm)
          if ( b == myrank_k+myrank_b*np_bzsm ) then
             nwork=nwork+1
             mapwork(1,nwork)=i
             mapwork(2,nwork)=j
          end if
       end do ! i
    end do ! j

    ir_fkmb(:)=0
    id_fkmb(:)=0
    do i=1,nwork
       j=mod(i-1,np_fkmb)
       ir_fkmb(j)=ir_fkmb(j)+1
    end do
    do j=0,np_fkmb-1
       id_fkmb(j)=sum(ir_fkmb(0:j))-ir_fkmb(j)
    end do

    nwork_0 = id_fkmb(myrank_f) + 1
    nwork_1 = id_fkmb(myrank_f) + ir_fkmb(myrank_f)

! ---

    allocate( trho(n1:n2) ) ; trho=zero
    allocate( tvht(n1:n2) ) ; tvht=zero

    call watch(ctt(1),ett(1))

! ---

    do iwork=nwork_0,nwork_1

       a = mapwork(1,iwork)
       b = mapwork(2,iwork)
       m = mapnk(1,a)
       q = mapnk(2,a)
       n = mapnk(1,b)
       k = mapnk(2,b)

! - normal part -

       do i=n1,n2
          trho(i) = conjg(unk_hf(i,m,q,s))*unk_hf(i,n,k,s)
       end do

       call Fock_fft( n1,n2,k,q,trho,tvht,1 )

       if ( abs(occ_hf(m,q,s)) >= 1.d-10 ) then
          c = alpha_hf*occ_factor*occ_hf(m,q,s)
          do i=n1,n2
             hunk(i,n,k,s)=hunk(i,n,k,s)-c*tvht(i)*unk_hf(i,m,q,s)
          end do
       end if

       if ( a /= b .and. abs(occ_hf(n,k,s)) >= 1.d-10 ) then
          c = alpha_hf*occ_factor*occ_hf(n,k,s)
          do i=n1,n2
             hunk(i,m,q,s)=hunk(i,m,q,s)-c*conjg(tvht(i))*unk_hf(i,n,k,s)
          end do
       end if

! - time-reversal part -

       do i=n1,n2
          trho(i) = unk_hf(i,m,q,s)*unk_hf(i,n,k,s)
       end do

       call Fock_fft( n1,n2,k,q,trho,tvht,2 )

       if ( abs(occ_hf(m,q,s)) >= 1.d-10 ) then
          c = alpha_hf*occ_factor*occ_hf(m,q,s)
          do i=n1,n2
             hunk(i,n,k,s)=hunk(i,n,k,s)-c*tvht(i)*conjg(unk_hf(i,m,q,s))
          end do
       end if

       if ( a /= b .and. abs(occ_hf(n,k,s)) >= 1.d-10 ) then
          c = alpha_hf*occ_factor*occ_hf(n,k,s)
          do i=n1,n2
             hunk(i,m,q,s)=hunk(i,m,q,s)-c*tvht(i)*conjg(unk_hf(i,n,k,s))
          end do
       end if

    end do ! iwork

    call watch(ctt(2),ett(2))

! ---

    m=size( hunk,1 )*size( hunk,2 )
    do k=1,MBZ
       call mpi_allreduce( MPI_IN_PLACE, hunk(n1,1,k,s), m, TYPE_MAIN, &
                           MPI_SUM, comm_band, ierr )
    end do

    m=size( hunk,1 )*size( hunk,2 )*size( hunk,3 )
    call mpi_allreduce( MPI_IN_PLACE, hunk(n1,1,1,s), m, TYPE_MAIN, &
                        MPI_SUM, comm_bzsm, ierr )

    m=size( hunk,1 )*size( hunk,2 )*size( hunk,3 )
    call mpi_allreduce( MPI_IN_PLACE, hunk(n1,1,1,s), m, TYPE_MAIN, &
                        MPI_SUM, comm_fkmb, ierr )

    call watch(ctt(3),ett(3))

    ct_focK_fft(6) = ctt(1) - ctt(0)
    et_focK_fft(6) = ett(1) - ett(0)
    ct_focK_fft(7) = ctt(2) - ctt(1)
    et_focK_fft(7) = ett(2) - ett(1)
    ct_focK_fft(8) = ctt(3) - ctt(2)
    et_focK_fft(8) = ett(3) - ett(2)

! ---

    deallocate( tvht ) 
    deallocate( trho )

    deallocate( mapwork )
    deallocate( mapnk )

! ---

    return
#endif
  END SUBROUTINE Fock_5

#ifdef _DRSDFT_
  SUBROUTINE Fock_Double( s, n1, n2, psi, tpsi )
    implicit none
    integer,intent(IN) :: s,n1,n2
    real(8),intent(IN)  :: psi(n1:n2)
    real(8),intent(INOUT) :: tpsi(n1:n2)
    complex(8),allocatable :: trho(:),tvht(:)
    real(8) :: c1,c2
    integer :: m1,m2,i,k

    k = 1

    allocate( trho(n1:n2) ) ; trho=(0.0d0,0.0d0)
    allocate( tvht(n1:n2) ) ; tvht=(0.0d0,0.0d0)

    do m1=FKMB_0,FKMB_1,2

       m2=min( m1+1, FKMB_1 )

       if ( abs(occ_hf(m1,k,s)) < 1.d-10 ) cycle

       if ( m1 < m2 ) then
          do i=n1,n2
             trho(i) = dcmplx( unk_hf(i,m1,k,s),unk_hf(i,m2,k,s) )*psi(i)
          end do
       else
          do i=n1,n2
             trho(i) = unk_hf(i,m1,k,s)*psi(i)
          end do
       end if

       call Fock_FFT_Double( n1, n2, trho, tvht )

       c1 = alpha_hf*(2.0d0*occ_factor)*occ_hf(m1,k,s)
       c2 = alpha_hf*(2.0d0*occ_factor)*occ_hf(m2,k,s)

       if ( m1 < m2 ) then
          do i=n1,n2
             tpsi(i) = tpsi(i) - c1*real( tvht(i) )*unk_hf(i,m1,k,s) &
                               - c2*aimag( tvht(i) )*unk_hf(i,m2,k,s)
          end do
       else
          do i=n1,n2
             tpsi(i) = tpsi(i) - c1*real( tvht(i) )*unk_hf(i,m1,k,s)
          end do
       end if

    end do ! m1

    deallocate( tvht ) 
    deallocate( trho )

    return

  END SUBROUTINE Fock_Double
#endif

END MODULE fock_module
