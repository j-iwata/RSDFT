MODULE cg_module

  use rgrid_module, only: zdV,dV
  use hamiltonian_module
  use cgpc_module
  use parallel_module
  use array_bound_module, only: ML_0,ML_1,MB_0,MB_1
  use cg_lobpcg_module, only: init_lobpcg, lobpcg
 !use cg_u_module, only: init_cg_u, cg_u
  use cggs_module
  use wf_module, only: hunk, iflag_hunk
  use kinetic_module, only: SYStype
  use watch_module
  use conjugate_gradient_g_module, only: conjugate_gradient_g, pp_kind
  use io_tools_module

  implicit none

  PRIVATE
  PUBLIC :: conjugate_gradient
  PUBLIC :: read_cg

  integer :: Ncg = 3
  integer :: iswitch_gs = 0
  integer :: iswitch_cg = 1
  logical :: flag_init_read = .true.
  logical :: flag_init_cg = .true.

CONTAINS


  SUBROUTINE read_cg
    implicit none
    call IOTools_readIntegerKeyword( "NCG" , Ncg )
    call IOTools_readIntegerKeyword( "ICG" , iswitch_cg )
    call IOTools_readIntegerKeyword( "SWGS", iswitch_gs )
    flag_init_read = .false.
  END SUBROUTINE read_cg


  SUBROUTINE init_cg
    implicit none
    logical :: disp
    if ( .not.flag_init_cg ) return
    call check_disp_switch( disp, 0 )
    if ( disp ) then
       write(*,*) "NCG =",Ncg
       write(*,*) "ICG =",iswitch_cg
       write(*,*) "SWGS=",iswitch_gs
    end if
    flag_init_cg=.false.
  END SUBROUTINE init_cg

  SUBROUTINE conjugate_gradient( n1,n2, MB, k,s, unk, esp, res )
    implicit none
    integer,intent(IN) :: n1,n2,MB,k,s
    real(8),intent(INOUT) :: esp(MB),res(MB)
#ifdef _DRSDFT_
    real(8),intent(INOUT) :: unk(n1:n2,MB)
#else
    complex(8),intent(INOUT) :: unk(n1:n2,MB)
#endif
    integer :: ipc

    call write_border( 1, " conjugate_gradient(start)" )

    call init_cg

    call init_cgpc( n1, n2, k, s, dV, SYStype, ipc )

    if ( pp_kind == "USPP" ) then

       call conjugate_gradient_g( n1,n2,MB,k,s,Ncg,unk,esp,res,iswitch_gs )

    else

       select case( iswitch_cg )
       case default
       case( 1 )

          call conjugate_gradient_1(n1,n2,MB,k,s,Ncg,unk,esp,res)

       case( 2 )

          call init_lobpcg( n1,n2,MB_0,MB_1,dV,MB_d,comm_grid )
          call lobpcg( k,s,Ncg,iswitch_gs,unk,esp,res )

       end select

    end if

    call write_border( 1, " conjugate_gradient(end)" )

  END SUBROUTINE conjugate_gradient

#ifdef _DRSDFT_
  SUBROUTINE conjugate_gradient_1(n1,n2,MB,k,s,Mcg,unk,esp,res)
    implicit none
    integer,intent(IN) :: n1,n2,MB,k,s,Mcg
    real(8),intent(INOUT) :: unk(n1:n2,MB)
    real(8),intent(INOUT) :: esp(MB),res(MB)
    integer :: ns,ne,nn,n,m,icg,ML0,Nhpsi,Npc,Ncgtot,ierr
    integer :: mm,icmp,i,TYPE_MAIN,timer_counter
    real(8),parameter :: ep0=0.d0
    real(8),parameter :: ep1=1.d-15
    real(8) :: rwork(9),W(2),c,d,r,c1,ct0,ct1,et0,et1,ctt(4),ett(4)
    real(8),allocatable :: sb(:),rb(:),E(:),E1(:),gkgk(:),bk(:)
    complex(8) :: work(9),zphase,ztmp
    real(8),parameter :: zero=0.d0
    real(8),allocatable :: hxk(:,:),hpk(:,:),gk(:,:),Pgk(:,:)
    real(8),allocatable :: pk(:,:),pko(:,:)
    real(8),allocatable :: vtmp2(:,:),wtmp2(:,:)
    real(8),allocatable :: utmp2(:,:),btmp2(:,:)
    real(8) :: timecg(2,16), ttmp(2), ttmp_cg(2)
    real(8) :: timecg_min(2,16),timecg_max(2,16)
    real(8) :: time_cgpc_min(2,16),time_cgpc_max(2,16)
    real(8) :: time_kine_min(2,16),time_kine_max(2,16)
    character(5) :: timecg_indx(7)
    character(5) :: time_cgpc_indx(13)

    call watchb( ttmp_cg ) ; ttmp(:)=ttmp_cg(:)

    TYPE_MAIN = MPI_REAL8

    timecg(:,:)=0.0d0
    time_hmlt(:,:)=0.0d0
    time_kine(:,:)=0.0d0
    time_nlpp(:,:)=0.0d0
    time_cgpc(:,:)=0.0d0

    timecg_indx(1:7) = (/"hamil","dotp","allr","prec","init","deall","tot"/)

    ML0 = ML_1-ML_0+1

    mm  = ML0  ; if (TYPE_MAIN==mpi_complex16) mm=2*ML0
    c1  = 2.d0 ; if (TYPE_MAIN==mpi_complex16) c1=1.d0
    icmp= 1    ; if (TYPE_MAIN==mpi_complex16) icmp=2

    Ncgtot = 0
    Nhpsi  = 0
    Npc    = 0

    allocate( hxk(n1:n2,MB_d), hpk(n1:n2,MB_d) )
    allocate( gk(n1:n2,MB_d) , Pgk(n1:n2,MB_d) )
    allocate( pk(n1:n2,MB_d) , pko(n1:n2,MB_d) )
    allocate( sb(MB_d),rb(MB_d) )
    allocate( E(MB_d),E1(MB_d),gkgk(MB_d),bk(MB_d) )
    allocate( vtmp2(6,MB_d),wtmp2(6,MB_d) )
    allocate( utmp2(2,2),btmp2(2,2) )

!$OMP parallel workshare
    res(:) = 0.0d0
    esp(:) = 0.0d0
!$OMP end parallel workshare

    call watchb( ttmp, timecg(:,5) )

    do ns=MB_0,MB_1
       ne=ns
       nn=ne-ns+1

       E1(:)=1.d10

       call watchb( ttmp )

       if ( iflag_hunk >= 1 ) then
!$OMP parallel workshare
          hxk(:,1:nn)=hunk(:,ns:ne,k,s)
!$OMP end parallel workshare
       else
          call hamiltonian(k,s,unk(n1,ns),hxk,n1,n2,ns,ne) ; Nhpsi=Nhpsi+1
       end if

       call watchb( ttmp, timecg(:,1) )

       do n=1,nn
          call dot_product(unk(n1,n+ns-1),hxk(n1,n),sb(n),dV,mm,1)
       end do

       call watchb( ttmp, timecg(:,2) )

       call mpi_allreduce(sb,E,nn,mpi_real8,mpi_sum,comm_grid,ierr)

       call watchb( ttmp, timecg(:,3) )

       do n=1,nn
!$OMP parallel do
          do i=n1,n2
             gk(i,n)=-c1*(hxk(i,n)-E(n)*unk(i,n+ns-1))
          end do
!$OMP end parallel do
          call dot_product(gk(n1,n),gk(n1,n),sb(n),dV,mm,1)
       end do

       call watchb( ttmp, timecg(:,2) )

       call mpi_allreduce(sb,rb,nn,mpi_real8,mpi_sum,comm_grid,ierr)

       call watchb( ttmp, timecg(:,3) )

       do icg=1,Mcg+1

          call watchb( ttmp )

          Ncgtot=Ncgtot+1

          do n=1,nn
!$OMP parallel workshare
             Pgk(n1:n2,n)=gk(n1:n2,n)
!$OMP end parallel workshare
          end do

          res(ns:ne)=rb(1:nn)/c1**2

! --- Convergence check ---

          if ( all(rb(1:nn)<ep0) ) exit
          if ( all(abs(E(1:nn)-E1(1:nn))<ep1) ) exit
          if ( icg==Mcg+1 ) exit

          call watchb( ttmp, timecg(:,2) )

! --- Preconditioning ---

          call preconditioning(E,k,s,nn,ML0,unk(n1,ns),gk,Pgk)

          call watchb( ttmp, timecg(:,4) )

! --- orthogonalization

          do n=ns,ne
             call cggs( iswitch_gs, ML0, MB, n, dV, unk(n1,1), Pgk(n1,n-ns+1) )
          end do

! ---

          do n=1,nn
             call dot_product(Pgk(n1,n),gk(n1,n),sb(n),dV,mm,1)
          end do

          call watchb( ttmp, timecg(:,2) )

          call mpi_allreduce(sb,rb,nn,mpi_real8,mpi_sum,comm_grid,ierr)

          call watchb( ttmp, timecg(:,3) )

          if ( icg==1 ) then
!$OMP parallel workshare
             pk(n1:n2,1:nn) = Pgk(n1:n2,1:nn)
!$OMP end parallel workshare
          else
             do n=1,nn
                bk(n)=rb(n)/gkgk(n)
!$OMP parallel do
                do i=n1,n2
                   pk(i,n)=Pgk(i,n)+bk(n)*pk(i,n)
                end do
!$OMP end parallel do
             end do
          end if
          gkgk(1:nn)=rb(1:nn)

          call watchb( ttmp, timecg(:,2) )

          call hamiltonian(k,s,pk,hpk,n1,n2,ns,ne) ; Nhpsi=Nhpsi+1

          call watchb( ttmp, timecg(:,1) )

          do n=1,nn
             vtmp2(1:6,n)=zero
             m=n+ns-1
             call dot_product(unk(n1,m),unk(n1,m),vtmp2(1,n),dV,mm,1)
             call dot_product(pk(n1,n),unk(n1,m),vtmp2(2,n),dV,mm,icmp)
             call dot_product(pk(n1,n),pk(n1,n),vtmp2(3,n),dV,mm,1)
             call dot_product(unk(n1,m),hxk(n1,n),vtmp2(4,n),dV,mm,1)
             call dot_product(pk(n1,n),hxk(n1,n),vtmp2(5,n),dV,mm,icmp)
             call dot_product(pk(n1,n),hpk(n1,n),vtmp2(6,n),dV,mm,1)
          end do

          call watchb( ttmp, timecg(:,2) )

          call mpi_allreduce(vtmp2,wtmp2,6*nn,TYPE_MAIN,mpi_sum,comm_grid,ierr)

          call watchb( ttmp, timecg(:,3) )

          do n=1,nn
             m=n+ns-1
             btmp2(1,1)=wtmp2(1,n)
             btmp2(2,1)=wtmp2(2,n)
             btmp2(1,2)=wtmp2(2,n)
             btmp2(2,2)=wtmp2(3,n)
             utmp2(1,1)=wtmp2(4,n)
             utmp2(2,1)=wtmp2(5,n)
             utmp2(1,2)=wtmp2(5,n)
             utmp2(2,2)=wtmp2(6,n)
             if (TYPE_MAIN==mpi_complex16) then
                ztmp=btmp2(1,2)
                ztmp=conjg(ztmp)
                btmp2(1,2)=ztmp
                ztmp=utmp2(1,2)
                ztmp=conjg(ztmp)
                utmp2(1,2)=ztmp
                call zhegv(1,'V','U',2,utmp2,2,btmp2,2,W,work,9,rwork,ierr)
                if ( abs(W(1)-E(n))>1.d-1 .and. abs(W(2)-E(n))<=1.d-1 ) then
                   utmp2(1,1)=utmp2(1,2)
                   utmp2(2,1)=utmp2(2,2)
                   W(1)=W(2)
                end if
!- Fix the phase -
                ztmp=utmp2(1,1)
                r=abs(ztmp)
                c=real(ztmp)/r
                d=aimag(ztmp)/r
                zphase=dcmplx(c,-d)
                utmp2(1,1)=utmp2(1,1)*zphase
                utmp2(2,1)=utmp2(2,1)*zphase
             else
                call dsygv(1,'V','U',2,utmp2,2,btmp2,2,W,rwork,9,ierr)
                if ( abs(W(1)-E(n))>1.d-1 .and. abs(W(2)-E(n))<=1.d-1 ) then
                   utmp2(1,1)=utmp2(1,2)
                   utmp2(2,1)=utmp2(2,2)
                   W(1)=W(2)
                end if
!- Fix the phase -
                c=utmp2(1,1)
                if( c<0.d0 ) then
                   utmp2(1,1)=-utmp2(1,1)
                   utmp2(2,1)=-utmp2(2,1)
                end if
             end if

             E1(n)=E(n)
             E(n) =W(1)

!$OMP parallel
!$OMP do
             do i=n1,n2
                hxk(i,n)=utmp2(1,1)*hxk(i,n)+utmp2(2,1)*hpk(i,n)
             end do
!$OMP end do
!$OMP do
             do i=n1,n2
                gk(i,n) = -c1*( hxk(i,n) &
                     -W(1)*(utmp2(1,1)*unk(i,m)+utmp2(2,1)*pk(i,n)) )
             end do
!$OMP end do
!$OMP end parallel

             call dot_product(gk(n1,n),gk(n1,n),sb(n),dV,mm,1)

          end do ! n

          call watchb( ttmp, timecg(:,2) )

          call mpi_allreduce(sb,rb,nn,mpi_real8,mpi_sum,comm_grid,ierr)

          call watchb( ttmp, timecg(:,3) )

          do n=1,nn
             m=n+ns-1
             if ( rb(n)/res(m)>1.d8 ) then
                E(n)=E1(n)
                cycle
             end if
!$OMP parallel do
             do i=n1,n2
                unk(i,m)=utmp2(1,1)*unk(i,m)+utmp2(2,1)*pk(i,n)
             end do
!$OMP end parallel do
             if ( iflag_hunk >= 1 ) then
!$OMP parallel do
                do i=n1,n2
                   hunk(i,m,k,s)=hxk(i,n)
                end do
!$OMP end parallel do
             end if
          end do

          call watchb( ttmp, timecg(:,2) )

       end do ! icg

       esp(ns:ne)=E(1:nn)

    end do  ! band-loop

    call watchb( ttmp )

    deallocate( btmp2,utmp2  )
    deallocate( wtmp2,vtmp2  )
    deallocate( bk,gkgk,E1,E )
    deallocate( rb,sb )
    deallocate( pko,pk  )
    deallocate( Pgk,gk  )
    deallocate( hpk,hxk )

    call watchb( ttmp, timecg(:,6) )
    call watchb( ttmp_cg, timecg(:,7) )

!    call get_time_min( 7, timecg, timecg_min )
!    call get_time_max( 7, timecg, timecg_max )
!    call get_time_min( 11, time_kine, time_kine_min )
!    call get_time_max( 11, time_kine, time_kine_max )
!    call get_time_min( 13, time_cgpc, time_cgpc_min )
!    call get_time_max( 13, time_cgpc, time_cgpc_max )

!    if ( disp_switch_parallel ) then
!       write(*,*) "time(hmlt_kin)",( time_hmlt(i,1), i=1,2 )
!       write(*,*) "time(hmlt_loc)",( time_hmlt(i,2), i=1,2 )
!       write(*,*) "time(hmlt_nlc)",( time_hmlt(i,3), i=1,2 )
!       write(*,*) "time(hmlt_exx)",( time_hmlt(i,4), i=1,2 )
!       write(*,'(a20," kine")') repeat("-",20)
!       call write_watchb( time_kine(1,6),6, time_kine_indx(6) ) 
!       write(*,*) "(min)"
!       call write_watchb( time_kine_min(1,6),6, time_kine_indx(6) ) 
!       write(*,*) "(max)"
!       call write_watchb( time_kine_max(1,6),6, time_kine_indx(6) ) 
!       call write_watchb( time_nlpp(1,1), 7, time_nlpp_indx ) 
!       write(*,'(a20," cgpc")') repeat("-",20)
!       call write_watchb( time_cgpc(1,8),6, time_cgpc_indx(8) ) 
!       write(*,*) "(min)"
!       call write_watchb( time_cgpc_min(1,8),6, time_cgpc_indx(8) ) 
!       write(*,*) "(max)"
!       call write_watchb( time_cgpc_max(1,8),6, time_cgpc_indx(8) ) 
!       write(*,'(a20," cg_1")') repeat("-",20)
!       call write_watchb( timecg(1,1), 7, timecg_indx ) 
!       write(*,*) "iswitch_gs=",iswitch_gs
!    end if

  END SUBROUTINE conjugate_gradient_1

#else

  SUBROUTINE conjugate_gradient_1(n1,n2,MB,k,s,Mcg,unk,esp,res)
    implicit none
    integer,intent(IN) :: n1,n2,MB,k,s,Mcg
    complex(8),intent(INOUT) :: unk(n1:n2,MB)
    real(8),intent(INOUT) :: esp(MB),res(MB)
    integer :: ns,ne,nn,n,m,icg,ML0,Nhpsi,Npc,Ncgtot,ierr
    integer :: mm,icmp,i,TYPE_MAIN
    real(8),parameter :: ep0=0.d0
    real(8),parameter :: ep1=1.d-15
    real(8) :: rwork(9),W(2),c,d,r,c1,ct0,ct1,et0,et1,ctt(4),ett(4)
    real(8),allocatable :: sb(:),rb(:),E(:),E1(:),gkgk(:),bk(:)
    complex(8) :: work(9),zphase,ztmp

    complex(8),parameter :: zero=(0.d0,0.d0)
    complex(8),allocatable :: hxk(:,:),hpk(:,:),gk(:,:),Pgk(:,:)
    complex(8),allocatable :: pk(:,:),pko(:,:)
    complex(8),allocatable :: vtmp2(:,:),wtmp2(:,:)
    complex(8),allocatable :: utmp2(:,:),btmp2(:,:)

    TYPE_MAIN = MPI_COMPLEX16

    ctt(:)=0.d0
    ett(:)=0.d0
    ctt_hamil(:)=0.d0
    ett_hamil(:)=0.d0

    ML0 = ML_1-ML_0+1

    mm  = ML0  ; if (TYPE_MAIN==mpi_complex16) mm=2*ML0
    c1  = 2.d0 ; if (TYPE_MAIN==mpi_complex16) c1=1.d0
    icmp= 1    ; if (TYPE_MAIN==mpi_complex16) icmp=2

    Ncgtot = 0
    Nhpsi  = 0
    Npc    = 0

    allocate( hxk(n1:n2,MB_d), hpk(n1:n2,MB_d) )
    allocate( gk(n1:n2,MB_d) , Pgk(n1:n2,MB_d) )
    allocate( pk(n1:n2,MB_d) , pko(n1:n2,MB_d) )
    allocate( sb(MB_d),rb(MB_d) )
    allocate( E(MB_d),E1(MB_d),gkgk(MB_d),bk(MB_d) )
    allocate( vtmp2(6,MB_d),wtmp2(6,MB_d) )
    allocate( utmp2(2,2),btmp2(2,2) )

!$OMP parallel workshare
    res(:)  = 0.d0
    esp(:)  = 0.d0
!$OMP end parallel workshare

    do ns=MB_0,MB_1
       ne=ns
       nn=ne-ns+1

       E1(:)=1.d10

       call watch(ct0,et0)

       if ( iflag_hunk >= 1 ) then
!$OMP parallel workshare
          hxk(:,1:nn)=hunk(:,ns:ne,k,s)
!$OMP end parallel workshare
       else
          call hamiltonian(k,s,unk(n1,ns),hxk,n1,n2,ns,ne) ; Nhpsi=Nhpsi+1
       end if

       call watch(ct1,et1) ; ctt(1)=ctt(1)+ct1-ct0 ; ett(1)=ett(1)+et1-et0

       do n=1,nn
          call dot_product(unk(n1,n+ns-1),hxk(n1,n),sb(n),dV,mm,1)
       end do

       call watch(ct0,et0) ; ctt(2)=ctt(2)+ct0-ct1 ; ett(2)=ett(2)+et0-et1

       call mpi_allreduce(sb,E,nn,mpi_real8,mpi_sum,comm_grid,ierr)

       call watch(ct1,et1) ; ctt(3)=ctt(3)+ct1-ct0 ; ett(3)=ett(3)+et1-et0

       do n=1,nn
!$OMP parallel do
          do i=n1,n2
             gk(i,n)=-c1*(hxk(i,n)-E(n)*unk(i,n+ns-1))
          end do
!$OMP end parallel do
          call dot_product(gk(n1,n),gk(n1,n),sb(n),dV,mm,1)
       end do

       call watch(ct0,et0) ; ctt(2)=ctt(2)+ct0-ct1 ; ett(2)=ett(2)+et0-et1

       call mpi_allreduce(sb,rb,nn,mpi_real8,mpi_sum,comm_grid,ierr)

       call watch(ct1,et1) ; ctt(3)=ctt(3)+ct1-ct0 ; ett(3)=ett(3)+et1-et0

       do icg=1,Mcg+1

          Ncgtot=Ncgtot+1

          do n=1,nn
!$OMP parallel workshare
             Pgk(n1:n2,n)=gk(n1:n2,n)
!$OMP end parallel workshare
          end do

          res(ns:ne)=rb(1:nn)/c1**2

! --- Convergence check ---

          if ( all(rb(1:nn)<ep0) ) exit
          if ( all(abs(E(1:nn)-E1(1:nn))<ep1) ) exit
          if ( icg==Mcg+1 ) exit

          call watch(ct1,et1) ; ctt(2)=ctt(2)+ct1-ct0 ; ett(2)=ett(2)+et1-et0

! --- Preconditioning ---

          call preconditioning(E,k,s,nn,ML0,unk(n1,ns),gk,Pgk)

          call watch(ct0,et0) ; ctt(4)=ctt(4)+ct0-ct1 ; ett(4)=ett(4)+et0-et1

! --- orthogonalization

          do n=ns,ne
             call cggs( iswitch_gs, ML0, MB, n, dV, unk(n1,1), Pgk(n1,n-ns+1) )
          end do

! ---

          do n=1,nn
             call dot_product(Pgk(n1,n),gk(n1,n),sb(n),dV,mm,1)
          end do

          call watch(ct1,et1) ; ctt(2)=ctt(2)+ct1-ct0 ; ett(2)=ett(2)+et1-et0

          call mpi_allreduce(sb,rb,nn,mpi_real8,mpi_sum,comm_grid,ierr)

          call watch(ct0,et0) ; ctt(3)=ctt(3)+ct0-ct1 ; ett(3)=ett(3)+et0-et1

          if ( icg==1 ) then
!$OMP parallel workshare
             pk(n1:n2,1:nn) = Pgk(n1:n2,1:nn)
!$OMP end parallel workshare
          else
             do n=1,nn
                bk(n)=rb(n)/gkgk(n)
!$OMP parallel do
                do i=n1,n2
                   pk(i,n)=Pgk(i,n)+bk(n)*pk(i,n)
                end do
!$OMP end parallel do
             end do
          end if
          gkgk(1:nn)=rb(1:nn)

          call watch(ct1,et1) ; ctt(2)=ctt(2)+ct1-ct0 ; ett(2)=ett(2)+et1-et0

          call hamiltonian(k,s,pk,hpk,n1,n2,ns,ne) ; Nhpsi=Nhpsi+1

          call watch(ct0,et0) ; ctt(1)=ctt(1)+ct0-ct1 ; ett(1)=ett(1)+et0-et1

          do n=1,nn
             vtmp2(1:6,n)=zero
             m=n+ns-1
             call dot_product(unk(n1,m),unk(n1,m),vtmp2(1,n),dV,mm,1)
             call dot_product(pk(n1,n),unk(n1,m),vtmp2(2,n),dV,mm,icmp)
             call dot_product(pk(n1,n),pk(n1,n),vtmp2(3,n),dV,mm,1)
             call dot_product(unk(n1,m),hxk(n1,n),vtmp2(4,n),dV,mm,1)
             call dot_product(pk(n1,n),hxk(n1,n),vtmp2(5,n),dV,mm,icmp)
             call dot_product(pk(n1,n),hpk(n1,n),vtmp2(6,n),dV,mm,1)
          end do

          call watch(ct1,et1) ; ctt(2)=ctt(2)+ct1-ct0 ; ett(2)=ett(2)+et1-et0

          call mpi_allreduce(vtmp2,wtmp2,6*nn,TYPE_MAIN,mpi_sum,comm_grid,ierr)

          call watch(ct0,et0) ; ctt(3)=ctt(3)+ct0-ct1 ; ett(3)=ett(3)+et0-et1

          do n=1,nn
             m=n+ns-1
             btmp2(1,1)=wtmp2(1,n)
             btmp2(2,1)=wtmp2(2,n)
             btmp2(1,2)=wtmp2(2,n)
             btmp2(2,2)=wtmp2(3,n)
             utmp2(1,1)=wtmp2(4,n)
             utmp2(2,1)=wtmp2(5,n)
             utmp2(1,2)=wtmp2(5,n)
             utmp2(2,2)=wtmp2(6,n)
             if (TYPE_MAIN==mpi_complex16) then
                ztmp=btmp2(1,2)
                ztmp=conjg(ztmp)
                btmp2(1,2)=ztmp
                ztmp=utmp2(1,2)
                ztmp=conjg(ztmp)
                utmp2(1,2)=ztmp
                call zhegv(1,'V','U',2,utmp2,2,btmp2,2,W,work,9,rwork,ierr)
                if ( abs(W(1)-E(n))>1.d-1 .and. abs(W(2)-E(n))<=1.d-1 ) then
                   utmp2(1,1)=utmp2(1,2)
                   utmp2(2,1)=utmp2(2,2)
                   W(1)=W(2)
                end if
!- Fix the phase -
                ztmp=utmp2(1,1)
                r=abs(ztmp)
                c=real(ztmp)/r
                d=aimag(ztmp)/r
                zphase=dcmplx(c,-d)
                utmp2(1,1)=utmp2(1,1)*zphase
                utmp2(2,1)=utmp2(2,1)*zphase
             else
                call dsygv(1,'V','U',2,utmp2,2,btmp2,2,W,rwork,9,ierr)
                if ( abs(W(1)-E(n))>1.d-1 .and. abs(W(2)-E(n))<=1.d-1 ) then
                   utmp2(1,1)=utmp2(1,2)
                   utmp2(2,1)=utmp2(2,2)
                   W(1)=W(2)
                end if
!- Fix the phase -
                c=utmp2(1,1)
                if( c<0.d0 ) then
                   utmp2(1,1)=-utmp2(1,1)
                   utmp2(2,1)=-utmp2(2,1)
                end if
             end if

             E1(n)=E(n)
             E(n) =W(1)

!$OMP parallel
!$OMP do
             do i=n1,n2
                hxk(i,n)=utmp2(1,1)*hxk(i,n)+utmp2(2,1)*hpk(i,n)
             end do
!$OMP end do
!$OMP do
             do i=n1,n2
                gk(i,n) = -c1*( hxk(i,n) &
                     -W(1)*(utmp2(1,1)*unk(i,m)+utmp2(2,1)*pk(i,n)) )
             end do
!$OMP end do
!$OMP end parallel

             call dot_product(gk(n1,n),gk(n1,n),sb(n),dV,mm,1)

          end do

          call watch(ct1,et1) ; ctt(2)=ctt(2)+ct1-ct0 ; ett(2)=ett(2)+et1-et0

          call mpi_allreduce(sb,rb,nn,mpi_real8,mpi_sum,comm_grid,ierr)

          call watch(ct0,et0) ; ctt(3)=ctt(3)+ct0-ct1 ; ett(3)=ett(3)+et0-et1

          do n=1,nn
             m=n+ns-1
             if ( rb(n)/res(m)>1.d8 ) then
                E(n)=E1(n)
                cycle
             end if
!$OMP parallel do
             do i=n1,n2
                unk(i,m)=utmp2(1,1)*unk(i,m)+utmp2(2,1)*pk(i,n)
             end do
!$OMP end parallel do
             if ( iflag_hunk >= 1 ) then
!$OMP parallel do
                do i=n1,n2
                   hunk(i,m,k,s)=hxk(i,n)
                end do
!$OMP end parallel do
             end if
          end do

          call watch(ct1,et1) ; ctt(2)=ctt(2)+ct1-ct0 ; ett(2)=ett(2)+et1-et0

       end do ! icg

       esp(ns:ne)=E(1:nn)

    end do  ! band-loop

    deallocate( btmp2,utmp2  )
    deallocate( wtmp2,vtmp2  )
    deallocate( bk,gkgk,E1,E )
    deallocate( rb,sb )
    deallocate( pko,pk  )
    deallocate( Pgk,gk  )
    deallocate( hpk,hxk )

!    if ( disp_switch_parallel ) then
!       write(*,*) "time(hamil_kin)",ctt_hamil(1),ett_hamil(1)
!       write(*,*) "time(hamil_loc)",ctt_hamil(2),ett_hamil(2)
!       write(*,*) "time(hamil_nlc)",ctt_hamil(3),ett_hamil(3)
!       write(*,*) "time(hamil_exx)",ctt_hamil(4),ett_hamil(4)
!       write(*,*) "time(hamil_cg)",ctt(1),ett(1)
!       write(*,*) "time(op_cg   )",ctt(2),ett(2)
!       write(*,*) "time(com_cg  )",ctt(3),ett(3)
!       write(*,*) "time(pc_cg   )",ctt(4),ett(4)
!       write(*,*) "iswitch_gs=",iswitch_gs
!    end if

  END SUBROUTINE conjugate_gradient_1
#endif

  SUBROUTINE get_time_min( n, t_in, t_min )
    implicit none
    integer,intent(IN)  :: n
    real(8),intent(IN)  :: t_in(2,n)
    real(8),intent(OUT) :: t_min(2,n)
    integer :: i
    real(8),allocatable :: t_tmp(:,:)
    allocate( t_tmp(2,n) ) ; t_tmp=0.0d0
    call mpi_allreduce( t_in, t_tmp, 2*n, mpi_real8, mpi_min, comm_grid, i )
    call mpi_allreduce( t_tmp, t_min, 2*n, mpi_real8, mpi_min, comm_band, i )
    deallocate( t_tmp )
  END SUBROUTINE get_time_min

  SUBROUTINE get_time_max( n, t_in, t_max )
    implicit none
    integer,intent(IN)  :: n
    real(8),intent(IN)  :: t_in(2,n)
    real(8),intent(OUT) :: t_max(2,n)
    integer :: i
    real(8),allocatable :: t_tmp(:,:)
    allocate( t_tmp(2,n) ) ; t_tmp=0.0d0
    call mpi_allreduce( t_in, t_tmp, 2*n, mpi_real8, mpi_max, comm_grid, i )
    call mpi_allreduce( t_tmp, t_max, 2*n, mpi_real8, mpi_max, comm_band, i )
    deallocate( t_tmp )
  END SUBROUTINE get_time_max


END MODULE cg_module
