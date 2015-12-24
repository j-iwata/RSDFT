MODULE conjugate_gradient_g_module
!----------------------------------------------------------------------------
! this module is conjugate gradient calculation for ultrasoft PS
!----------------------------------------------------------------------------

  use rgrid_module, only: zdV,dV
  use hamiltonian_module
  use cgpc_module
  use parallel_module
  use array_bound_module, only: ML_0,ML_1,MB_0,MB_1
  use watch_module
  use inner_product_module, only: get_Sf,get_gSf
  use real_complex_module, only: TYPE_MAIN,zero
  use var_sys_parameter, only: pp_kind

  implicit none

  PRIVATE
  PUBLIC :: conjugate_gradient_g
  PUBLIC :: pp_kind

CONTAINS


  SUBROUTINE conjugate_gradient_g( n1,n2,MB,k,s,Mcg,unk,esp,res,iswitch_gs )
    implicit none
    integer,intent(IN) :: n1,n2,MB,k,s,Mcg,iswitch_gs
#ifdef _DRSDFT_
    real(8),intent(INOUT) :: unk(n1:n2,MB)
#else
    complex(8),intent(INOUT) :: unk(n1:n2,MB)
#endif
    real(8),intent(INOUT) :: esp(MB),res(MB)
    integer :: ns,ne,nn,n,m,icg,ML0,Nhpsi,Npc,Ncgtot,ierr
    integer :: mm,icmp,i
    real(8),parameter :: ep0=0.d0
    real(8),parameter :: ep1=1.d-15
    real(8) :: rwork(9),W(2),c,d,r,c1,ct0,ct1,et0,et1,ctt(4),ett(4)
    real(8),allocatable :: sb(:),rb(:),E(:),E1(:),gkgk(:),bk(:)


#ifdef _DRSDFT_
    real(8) :: tmp
    real(8) :: zphase,ztmp

    real(8),allocatable :: hxk(:,:),hpk(:,:),gk(:,:),Pgk(:,:)
    real(8),allocatable :: pk(:,:),pko(:,:)
    real(8),allocatable :: vtmp2(:,:),wtmp2(:,:)
    real(8),allocatable :: utmp2(:,:),btmp2(:,:)
!----- USPP -----
    real(8) :: gSf
    real(8),allocatable :: unk_tmp(:),Sf(:)
!===== USPP =====
#else
    complex(8) :: tmp
    complex(8) :: work(9),zphase,ztmp

    complex(8),allocatable :: hxk(:,:),hpk(:,:),gk(:,:),Pgk(:,:)
    complex(8),allocatable :: pk(:,:),pko(:,:)
    complex(8),allocatable :: vtmp2(:,:),wtmp2(:,:)
    complex(8),allocatable :: utmp2(:,:),btmp2(:,:)
!----- USPP -----
    complex(8) :: gSf
    complex(8),allocatable :: unk_tmp(:),Sf(:)
!===== USPP =====
#endif
    

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

    allocate( hxk(n1:n2,MB_d), hpk(n1:n2,MB_d) ) ; hxk=zero ; hpk=zero
    allocate( gk(n1:n2,MB_d) , Pgk(n1:n2,MB_d) ) ; gk=zero ; Pgk=zero
    allocate( pk(n1:n2,MB_d) , pko(n1:n2,MB_d) ) ; pk=zero ; pko=zero
    allocate( sb(MB_d),rb(MB_d) ) ; sb=0.d0 ; rb=0.d0
    allocate( E(MB_d),E1(MB_d),gkgk(MB_d),bk(MB_d) )
    E=0.d0 ; E1=0.d0 ; gkgk=0.d0 ; bk=0.d0
    allocate( vtmp2(6,MB_d),wtmp2(6,MB_d) ) ; vtmp2=zero ; wtmp2=zero
    allocate( utmp2(2,2),btmp2(2,2) ) ; utmp2=zero ; btmp2=zero
    allocate( Sf(n1:n2) ) ; Sf=zero
    allocate( unk_tmp(n1:n2) ) ; unk_tmp=zero

!$OMP parallel workshare
    res(:)  = 0.d0
    esp(:)  = 0.d0
!$OMP end parallel workshare

    do ns=MB_0,MB_1
       ne=ns
       nn=ne-ns+1

       E1(:)=1.d10

       call watch(ct0,et0)
       call hamiltonian(k,s,unk(n1,ns),hxk,n1,n2,ns,ne) ; Nhpsi=Nhpsi+1
       call watch(ct1,et1) ; ctt(1)=ctt(1)+ct1-ct0 ; ett(1)=ett(1)+et1-et0

       do n=1,nn
          call dot_product(unk(n1,n+ns-1),hxk(n1,n),sb(n),dV,mm,1)
       end do

       call watch(ct0,et0) ; ctt(2)=ctt(2)+ct0-ct1 ; ett(2)=ett(2)+et0-et1
       call mpi_allreduce(sb,E,nn,mpi_real8,mpi_sum,comm_grid,ierr)
       call watch(ct1,et1) ; ctt(3)=ctt(3)+ct1-ct0 ; ett(3)=ett(3)+et1-et0

       do n=1,nn
!----- USPP -----
          call get_Sf( unk(n1,n+ns-1),n1,n2,k,Sf )
!$OMP parallel do
          do i=n1,n2
             gk(i,n)=-c1*(hxk(i,n)-E(n)*Sf(i))
          end do
!$OMP end parallel do
          call get_gSf( gk,gk,n1,n2,k,gSf,0 )
          sb(n) = real( gSf,kind=8 )
!===== USPP =====
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

! ---

          do n=1,nn
!----- USPP -----
!--- <RSDFT-original>
!             call get_gSf( Pgk(n1,n),gk(n1,n),n1,n2,k,gSf,0 )
!--- <TAPP-style>
             call get_gSf( Pgk(n1,n),Pgk(n1,n),n1,n2,k,gSf,0 )
! <Pgk|S|Pgk> is a real number
! get_gSf gets a complex number. This is a waste of time.
             sb(n) = real( gSf,kind=8 )
!===== USPP =====
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
!----- USPP -----
! need omp_parallel
!             call get_gSf( unk(n1,m,k,s),unk(n1,m,k,s),n1,n2,k,vtmp2(1:n),0 )
!             call get_gSf( pk(n1,n),unk(n1,m,k,s),n1,n2,k,vtmp2(2:n),0 )
             call get_Sf( unk(n1,m),n1,n2,k,Sf )
             gSf=zero
             do i=n1,n2
#ifdef _DRSDFT_
                gSf = gSf + unk(i,m)*Sf(i)
#else
                gSf = gSf + conjg(unk(i,m))*Sf(i)
#endif
             end do
! WARNING zdV or dV
             vtmp2(1,n) = dV*gSf
             gSf=zero
             do i=n1,n2
#ifdef _DRSDFT_
                gSf = gSf + pk(i,n)*Sf(i)
#else
                gSf = gSf + conjg(pk(i,n))*Sf(i)
#endif
             end do
             vtmp2(2,n) = dV*gSf
             call get_gSf( pk(n1,n),pk(n1,n),n1,n2,k,vtmp2(3,n),0 )
!===== USPP =====
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
#ifdef _DRSDFT_
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
#else
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
#endif
             E1(n)=E(n)
             E(n) =W(1)

!$OMP parallel
!$OMP do
             do i=n1,n2
                hxk(i,n)=utmp2(1,1)*hxk(i,n)+utmp2(2,1)*hpk(i,n)
             end do
!$OMP end do
!----- USPP -----
!$OMP do
             do i=n1,n2
                unk_tmp(i) = utmp2(1,1)*unk(i,m) + utmp2(2,1)*pk(i,n)
             end do
!$OMP end do
!$OMP end parallel
             call get_Sf( unk_tmp,n1,n2,k,Sf )
!$OMP parallel
!$OMP do
             do i=n1,n2
                gk(i,n) = -c1*( hxk(i,n) - W(1)*Sf(i) )
             end do
!$OMP end do
!$OMP end parallel

             call get_gSf( gk(n1,n),gk(n1,n),n1,n2,k,gSf,0 )
             sb(n) = real( gSf,kind=8 )
!===== USPP =====

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
    deallocate( Sf )
    deallocate( unk_tmp )

!    if ( disp_switch_parallel ) then
!       write(*,*) "time(hamil_kin)",ctt_hamil(1),ett_hamil(1)
!       write(*,*) "time(hamil_loc)",ctt_hamil(2),ett_hamil(2)
!       write(*,*) "time(hamil_nlc)",ctt_hamil(3),ett_hamil(3)
!       write(*,*) "time(hamil_exx)",ctt_hamil(4),ett_hamil(4)
!       write(*,*) "time(hamil_cg)",ctt(1),ett(1)
!       write(*,*) "time(op_cg   )",ctt(2),ett(2)
!       write(*,*) "time(com_cg  )",ctt(3),ett(3)
!       write(*,*) "time(pc_cg   )",ctt(4),ett(4)
!    end if

  END SUBROUTINE conjugate_gradient_g


END MODULE conjugate_gradient_g_module
