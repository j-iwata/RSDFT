MODULE atomopt_module

  use parallel_module
  use atom_module
  use total_energy_module
  use bb_module
  use scf_module
  use eion_module, only: calc_eion
  use strfac_module
  use ps_local_module
  use ps_pcc_module
  use pseudopot_module
  use ps_nloc2_module
  use ps_nloc3_module
  use ps_nloc_mr_module
  use force_module
  use kinetic_module, only: SYStype
  use ps_local_mol_module, only: construct_ps_local_mol
  use ps_nloc2_mol_module
  use ps_pcc_mol_module
  use eion_mol_module
  use ps_qrij_prep_module
  use ps_prepNzqr_g_module, only: prepNzqr
  use vdw_grimme_module

  implicit none

  PRIVATE
  PUBLIC :: ncycl,most,nrfr,okatom,eeps,feps,decr &
           ,read_atomopt,atomopt,read_oldformat_atomopt

  integer :: ncycl,most,nrfr
  real(8) :: okatom,eeps,feps,decr

  logical :: disp_switch_loc
  integer :: diter_opt

  integer :: strlog = 0
  integer,parameter :: unit_strlog = 297
  integer,parameter :: unit197 = 197
  integer,parameter :: unit97 = 97

CONTAINS


  SUBROUTINE read_atomopt(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i
    character(8) :: cbuf,ckey
    ncycl     = 0
    most      = 6
    nrfr      = 5
    diter_opt = 50
    okatom    = 0.5d0
    eeps      = 1.d-10
    feps      = 5.d-4
    decr      = 1.d-1

    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:8) == "ATOMOPT1" ) then
             backspace(unit)
             read(unit,*) cbuf,ncycl,most,nrfr
          else if ( ckey(1:8) == "ATOMOPT2" ) then
             backspace(unit)
             read(unit,*) cbuf,okatom,eeps,feps,decr
          else if ( ckey(1:8) == "ATOMOPT3" ) then
             backspace(unit)
             read(unit,*) cbuf,diter_opt
          else if ( ckey(1:8) == "STRLOG" ) then
             backspace(unit)
             read(unit,*) cbuf,strlog
          end if
       end do
999    continue
       write(*,*) "ncycl, most, nrfr =",ncycl,most,nrfr
       write(*,*) "okatom, eeps      =",okatom,eeps
       write(*,*) "feps, decr        =",feps,decr
       write(*,*) "diter_opt         =",diter_opt
       if ( diter_opt <= 0 ) then
          diter_opt=50
          write(*,*) "diter_opt         =",diter_opt
       end if
       write(*,*) "strlog            =",strlog
    end if
    call send_atomopt
  END SUBROUTINE read_atomopt


  SUBROUTINE read_oldformat_atomopt(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    if ( rank == 0 ) then
       read(unit,*) ncycl,most,nrfr,diter_opt
       read(unit,*) okatom,eeps,feps,decr
       write(*,*) "ncycl, most, nrfr =",ncycl,most,nrfr
       write(*,*) "diter_opt         =",diter_opt
       write(*,*) "okatom, eeps      =",okatom,eeps
       write(*,*) "feps, decr        =",feps,decr
       if ( diter_opt <= 0 ) then
          diter_opt=50
          write(*,*) "diter_opt         =",diter_opt
       end if
    end if
    call send_atomopt
  END SUBROUTINE read_oldformat_atomopt


  SUBROUTINE send_atomopt
    implicit none
    integer :: ierr
    call mpi_bcast(ncycl,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(most ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nrfr ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(okatom,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(eeps  ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(feps  ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(decr  ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(diter_opt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(strlog,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_atomopt


  SUBROUTINE atomopt(iswitch_opt,disp_switch)
    implicit none
    integer,intent(IN) :: iswitch_opt
    logical,intent(INOUT) :: disp_switch
    integer,parameter :: max_nhist=100000
    integer :: SCF_hist(max_nhist),ICY_hist(max_nhist)
    integer :: LIN_hist(max_nhist)
    integer :: a,nhist,most0,ncycl0,ifar,ierr,icy,nhist0
    integer :: i,itlin,amax,isafe,icflag,iter_final
    real(8) :: Fmax_hist(max_nhist),Etot_hist(max_nhist)
    real(8) :: dmax_hist(max_nhist)
    real(8) :: grad_hist(max_nhist),alpha_hist(max_nhist)
    real(8) :: Fmax,ss,Etot_save,dif,alpha1_0,okstep0
    real(8) :: al(3),all(3),ar(3),arr(3),alp,wd,xmin,emin
    real(8) :: Etot,Etot0,Fmax0_2,alpha0,Fmax0,Fmax00,Etot00,Etsave
    real(8) :: gh,alpha,tmp,ddmax,almax,okstep,alpha0_2,Etot0_2
    real(8) :: alpha2,alpha1,grad0,signh,dif0,gamma,gigi,hh
    real(8) :: alsave,safe,c0,c1,c2,c3,pi,ddmin,grad,safety
    real(8),allocatable :: Force(:,:),aa_atom_0(:,:)
    real(8),allocatable :: gi(:,:),hi(:,:)
    character(22) :: loop_info

    if ( disp_switch ) write(*,'(a60," atomopt")') repeat("-",60)

    call check_disp_switch( disp_switch_loc, 0 )
!    call check_disp_switch( .false., 1 )

    ddmin  = 1.d-8
    safe   = 0.01d0
    safety = 0.01d0
    pi     = acos(-1.d0)

!- allocate ---------------------------------------
    allocate( Force(3,Natom)    ) ; Force=0.d0
    allocate( aa_atom_0(3,Natom) ) ; aa_atom_0=0.d0
    allocate( gi(3,Natom)       ) ; gi=0.d0
    allocate( hi(3,Natom)       ) ; hi=0.d0
!--------------------------------------------------

    if ( iswitch_opt < 2  ) then

       if ( iswitch_opt == 1 ) then
          call calc_total_energy( .false., Etot )
          if ( disp_switch_loc ) write(*,*) "Etot(har)=",Etot
       end if

       call calc_force( Natom, Force )

       Fmax=0.d0
       do a=1,Natom
          ss=Force(1,a)**2+Force(2,a)**2+Force(3,a)**2
          Fmax=max(Fmax,ss)
       end do
       Fmax=sqrt(Fmax)

       Etot_save = 0.d0
       dif       = 0.0d0

       nhist                = 1
       Etot_hist(nhist)     = Etot
       alpha_hist(nhist)    = 0.d0
       Fmax_hist(nhist)     = Fmax
       SCF_hist(nhist)      = 0
       ICY_hist(nhist)      = 0
       LIN_hist(nhist)      = 0
       dmax_hist(nhist)     = 0

       ncycl0   = 1
       most0    = 1

! --- Convergence check 1 ---

       if ( Fmax <= feps ) then
          if ( disp_switch_loc ) then
             write(*,*) 'Fmax,feps =',Fmax,feps
          end if
          goto 999
       end if

! --- Read the previous optimization information ---

    else if ( iswitch_opt >= 2 ) then

       if ( myrank == 0 ) then
          open(1,file="wopt.dat",form='unformatted',status='old')
          read(1) ncycl0,most0,ifar
          read(1) Etot,dif,alpha1_0,okstep0
          read(1) al(1:3),all(1:3),ar(1:3),arr(1:3)
          read(1) Force(1:3,1:Natom)
          read(1) hi(1:3,1:Natom)
          read(1) aa_atom(1:3,1:Natom)
          close(1)
       end if
       call mpi_bcast(ncycl0,1,mpi_integer,0,mpi_comm_world,ierr)
       call mpi_bcast(most0,1,mpi_integer,0,mpi_comm_world,ierr)
       call mpi_bcast(ifar,1,mpi_integer,0,mpi_comm_world,ierr)
       call mpi_bcast(Etot,1,mpi_real8,0,mpi_comm_world,ierr)
       call mpi_bcast(dif,1,mpi_real8,0,mpi_comm_world,ierr)
       call mpi_bcast(alpha1_0,1,mpi_real8,0,mpi_comm_world,ierr)
       call mpi_bcast(okstep0,1,mpi_real8,0,mpi_comm_world,ierr)
       call mpi_bcast(al,3,mpi_real8,0,mpi_comm_world,ierr)
       call mpi_bcast(all,3,mpi_real8,0,mpi_comm_world,ierr)
       call mpi_bcast(ar,3,mpi_real8,0,mpi_comm_world,ierr)
       call mpi_bcast(arr,3,mpi_real8,0,mpi_comm_world,ierr)
       call mpi_bcast(Force,3*Natom,mpi_real8,0,mpi_comm_world,ierr)
       call mpi_bcast(hi,3*Natom,mpi_real8,0,mpi_comm_world,ierr)
       call mpi_bcast(aa_atom,3*Natom,mpi_real8,0,mpi_comm_world,ierr)

       Fmax=0.d0
       do a=1,Natom
          ss=Force(1,a)**2+Force(2,a)**2+Force(3,a)**2
          Fmax=max(Fmax,ss)
       end do
       Fmax=sqrt(Fmax)

       nhist = 0

    end if

    call write_atomic_coordinates_log(197,0,0,strlog,iswitch_opt)

!    disp_switch = .false.
!    disp_switch_parallel = .false.

    dif0 = dif

    if ( disp_switch_loc ) then
       write(*,*) "ncycl0,ncycl ",ncycl0,ncycl
       write(*,*) "most  ",most
       write(*,*) "okatom",okatom
    end if

!
! -------------------- CG-loop start --------------------
!

    opt_ion : do icy=ncycl0,ncycl0+ncycl-1

!       if ( disp_switch_loc ) then
!          write(*,'(a57," ICY (",i5,")")') repeat("-",57),icy
!       end if
       write(loop_info,'(" ICY    (",i5,")")') icy
       call write_border( 0, loop_info(1:len_trim(loop_info)) )

!
! --- Ion configuration ---
!

       aa_atom_0(1:3,1:Natom) = aa_atom(1:3,1:Natom)

!
! gi ---> gradient ( =-grad(Etot)=Force )
! hi ---> search direction
!
       if ( mod(icy-1,nrfr) == 0 ) then

          gamma=0.d0
          if ( icy>1 .and. disp_switch_loc ) then
             write(*,*) 'CG-direction is refreshed !!!'
          else
             if ( disp_switch_loc ) write(*,*) 'The first CG step !'
          end if

       else

          gigi=sum(gi(:,:)*gi(:,:))
          if ( gigi>0.d0 ) then
             gamma=sum((Force(:,:)-gi(:,:))*Force(:,:))/gigi
          end if

       end if

       gi(1:3,1:Natom)=Force(1:3,1:Natom)

       if ( iswitch_opt >= 2 .and. icy == ncycl0 ) then
       else
          hi(1:3,1:Natom)=gi(1:3,1:Natom)+gamma*hi(1:3,1:Natom)
       end if

       Etot_save = Etot
       hh        = sqrt( sum(hi(:,:)*hi(:,:)) )
       gh        = sum( gi(:,:)*hi(:,:) )
       alpha     = 2.d0*abs(dif/gh)
       if ( dif == 0.0d0 ) alpha = 0.5d0

!
! --- Check alpha ---
!

       okstep = okatom
!chstep
       tmp    = 0.d0
       amax   = 1
       do a=1,Natom
          ss=hi(1,a)*hi(1,a)+hi(2,a)*hi(2,a)+hi(3,a)*hi(3,a)
          if ( ss > tmp ) then
             tmp =ss
             amax=a
          end if
       end do
       ddmax=sqrt(tmp)*abs(alpha)
       if ( ddmax < 1.d-14 ) then
          if ( disp_switch_loc ) then
             write(*,*) "ddmax is too small : ddmax =",ddmax
             write(*,*) "alpha,amax =",alpha,amax
          end if
          ddmax=1.d-12
       end if
       almax=okstep/ddmax*abs(alpha)
!chstep
       if ( disp_switch_loc ) then
          write(*,'(1x,"Maximum displacement size =" &
               ,f16.7,3x,"( atom =",i5," )")') ddmax,amax
          write(*,'(1x,"okstep, alpha, almax =",3g16.6)') okstep,alpha,almax
       end if
       if ( ddmax > okstep ) then
          alpha=almax
          if ( disp_switch_loc ) then
             write(*,*) "alpha is too large and replaced by almax."
          end if
       end if

!
! --- hi-direction component of gi, and its sign ---
!

       grad0 = gh
       if ( disp_switch_loc ) write(*,*) "grad0 =",grad0
       if ( grad0 >= 0.d0 ) then
          signh=1.d0
       else
          if ( disp_switch_loc ) write(*,*) 'sign of hi is changed.'
          signh=-1.d0
          hi(1:3,1:Natom)=-hi(1:3,1:Natom)
          grad0=-grad0
       end if

       if ( iswitch_opt>=2 .and. icy==ncycl0 ) then

          alpha1   = alpha1_0
          alpha1_0 = 0.d0
          alpha2   = 0.d0
          okstep   = okstep0

       else

          all(1:3) =-1.d15
          ar(1:3)  = 1.d15
          arr(1:3) = 1.d15

          alpha1 = 0.d0
          alpha2 = 0.d0

          ifar = 0

          al(1) = 0.d0
          al(2) = Etot
          al(3) = grad0

       end if

       Etot00 = Etot
       Fmax00 = Fmax

       if ( disp_switch_loc ) then
          write(*,*) "initial configuration"
          if ( Natom <= 11 ) then
             do a=1,Natom
                write(*,'(1x,i4,3f15.5)') a,aa_atom(:,a)
             end do
          else
             do a=1,min(5,Natom)
                write(*,'(1x,i4,3f15.5)') a,aa_atom(:,a)
             end do
             write(*,'(1x,10x,".")')
             write(*,'(1x,10x,".")')
             write(*,'(1x,10x,".")')
             do a=Natom-5,Natom
                write(*,'(1x,i4,3f15.5)') a,aa_atom(:,a)
             end do
          end if
          write(*,*) "Etot =",Etot
          write(*,*) "grad =",grad0
          write(*,*) "Fmax =",Fmax
       end if

       nhist = nhist + 1

       Etot_hist(nhist)     = Etot
       alpha_hist(nhist)    = 0.d0
       Fmax_hist(nhist)     = Fmax
       grad_hist(nhist)     = grad0
       SCF_hist(nhist)      = 0
       ICY_hist(nhist)      = icy
       LIN_hist(nhist)      = 0
       dmax_hist(nhist)     = sqrt(maxval(hi(1,:)**2+hi(2,:)**2+hi(3,:)**2))

       nhist0 = nhist

       Fmax0  = Fmax
       alpha0 = 0.d0
       Etot0  = Etot

       Fmax0_2  = Fmax
       alpha0_2 = 0.d0
       Etot0_2  = Etot

!
! ---------- Line minimization start (along hi) ----------
!

       linmin : do itlin=most0,most

!          if ( disp_switch_loc ) then
!             write(*,'(a57," ICY    (",i5,")")') repeat("-",57),icy
!             write(*,'(a57," LINMIN (",i5,")")') repeat("-",57),itlin
!          end if
          write(loop_info,'(" ICY    (",i5,")")') icy
          call write_border( 0, loop_info(1:len_trim(loop_info)) )
          write(loop_info,'(" LINMIN (",i5,")")') itlin
          call write_border( 0, loop_info(1:len_trim(loop_info)) )

          Etsave = Etot
          emin   = 0.d0
          ddmax  = 1000.d0

          if ( disp_switch_loc ) write(*,*) "ifar     =",ifar

          alpha1_0 = alpha1
          okstep0  = okstep

          if ( itlin==1 ) then

             alpha1=alpha

          else

             if ( ifar==0 ) then

                okstep=okstep*2.d0
                alpha2=alpha1
                call get_min_parabola(all,al,xmin,emin)
                alpha1=xmin
                alp=alpha1-alpha2
                if ( disp_switch_loc ) then
                   write(*,*) "alpha1,alpha2=",alpha1,alpha2
                end if
!chstep
                tmp=0.d0
                amax=1
                do a=1,Natom
                   ss=hi(1,a)*hi(1,a)+hi(2,a)*hi(2,a)+hi(3,a)*hi(3,a)
                   if ( ss > tmp ) then
                      tmp=ss
                      amax=a
                   end if
                end do
                ddmax=sqrt(tmp)*abs(alp)
                if ( ddmax < 1.d-14 ) then
                   if ( disp_switch_loc ) then
                      write(*,*) "ddmax is too small : ddmax =",ddmax
                      write(*,*) "alp1,amax =",alp,amax
                   end if
                   ddmax=1.d-12
                end if
                almax=okstep/ddmax*abs(alp)
!chstep
                if ( disp_switch_loc ) then
                   write(*,'(1x,"Maximum displacement =",f16.7 &
                        ,3x,"( atom =",i5," )")') ddmax,amax
                   write(*,'(1x,"okstep, alp, almax =",3g16.6)') &
                        okstep,alp,almax
                end if
                if ( ddmax > okstep ) then
                   if ( disp_switch_loc ) then
                      write(*,*) "alpha1-alpha2 is " &
                           ,"too large and replaced by almax."
                   end if
                   alpha1=alpha2+almax
                   ddmax=okstep
                else if ( alpha1 < alpha2 ) then
                   if ( disp_switch_loc ) then
                      write(*,*) "alpha1 is irregular"
                   end if
                   tmp=abs( al(1)-all(1) )*10.d0
                   alpha1=alpha2+tmp
                   ddmax=tmp/almax*okstep
                end if

             else !( ifar/=0 )

                alpha2=alpha1
                call get_min_parabola(al,ar,xmin,emin)
                alpha1=xmin
                if ( .not.( (al(1)<alpha1).and.(alpha1<ar(1)) ) ) then
                   alpha1=0.5d0*(al(1)+ar(1))
                end if
                alp=alpha1-alpha2
!chstep
                tmp=0.d0
                amax=1
                do a=1,Natom
                   ss=hi(1,a)*hi(1,a)+hi(2,a)*hi(2,a)+hi(3,a)*hi(3,a)
                   if ( ss>tmp ) then
                      tmp=ss
                      amax=a
                   end if
                end do
                ddmax=sqrt(tmp)*abs(alp)
                if ( ddmax < 1.d-14 ) then
                   if ( disp_switch_loc ) then
                      write(*,*) "ddmax is too small : ddmax =",ddmax
                      write(*,*) "alp1,amax =",alp,amax
                   end if
                   ddmax=1.d-12
                end if
                almax=okstep/ddmax*abs(alp)
!chstep
                if ( disp_switch_loc ) then
                   write(*,'(1x,"Maximum displacement =" &
                        ,f16.7,3x,"( atom =",i5," )")') ddmax,amax
                   write(*,'(1x,"okstep, alp, almax =",3g16.6)') &
                        okstep,alp,almax
                end if

!
! --- safety mechanism for stable convergence ---
!

                wd=abs( ar(1)-al(1) )
                isafe=0
                alsave=alpha1
                if ( abs(alpha1-al(1)) < safety*wd ) then
                   isafe=1
                   alpha1=al(1)+safe*(ar(1)-al(1))
                else if ( abs(ar(1)-alpha1) < safety*wd ) then
                   isafe=1
                   alpha1=ar(1)-safe*(ar(1)-al(1))
                end if
                if ( isafe == 1 ) then
                   ddmax=ddmax*abs( (alpha1-alpha2)/alp )
                   if ( disp_switch_loc ) then
                      write(*,*) "Safety mechanism is applied."
                      write(*,*) "safety & safe =",safety,safe
                      write(*,*) "alpha1 is replaced by ",alpha1
                   end if
                end if

             end if

          end if

!
! ---- save the optimization information ---
!

          if ( myrank == 0 ) then
             open(1,file="wopt.dat",form='unformatted')
             write(1) icy,itlin,ifar
             write(1) Etot_save,dif0,alpha1_0,okstep0
             write(1) al,all,ar,arr
             write(1) gi(1:3,1:Natom)
             write(1) hi(1:3,1:Natom)
             write(1) aa_atom_0(1:3,1:Natom)
             close(1)
          end if

!
! --- Trial Configuration ---
!                  
          if ( SYStype == 0 ) then

             c0=alpha1/(2.d0*pi)
             do a=1,Natom
                c1=c0*hi(1,a)
                c2=c0*hi(2,a)
                c3=c0*hi(3,a)
                aa_atom(1,a)=aa_atom_0(1,a)+bb(1,1)*c1+bb(2,1)*c2+bb(3,1)*c3
                aa_atom(2,a)=aa_atom_0(2,a)+bb(1,2)*c1+bb(2,2)*c2+bb(3,2)*c3
                aa_atom(3,a)=aa_atom_0(3,a)+bb(1,3)*c1+bb(2,3)*c2+bb(3,3)*c3
             end do

          else if ( SYStype == 1 ) then

             do a=1,Natom
                aa_atom(1,a) = aa_atom_0(1,a) + alpha1*hi(1,a)
                aa_atom(2,a) = aa_atom_0(2,a) + alpha1*hi(2,a)
                aa_atom(3,a) = aa_atom_0(3,a) + alpha1*hi(3,a)
             end do

          end if

          if ( disp_switch_loc ) then
             write(*,*) 'Trial configuration (see fort.97)'
             if ( Natom <= 11 ) then
                do a=1,Natom
                   write(* ,'(1x,i5,3f20.12,i4)') &
                        ki_atom(a),aa_atom(:,a),md_atom(a)
                end do
             else
                do a=1,min(5,Natom)
                   write(* ,'(1x,i5,3f20.12,i4)') &
                        ki_atom(a),aa_atom(:,a),md_atom(a)
                end do
                write(*,'(1x,10x,".")')
                write(*,'(1x,10x,".")')
                write(*,'(1x,10x,".")')
                do a=Natom-5,Natom
                   write(* ,'(1x,i5,3f20.12,i4)') &
                        ki_atom(a),aa_atom(:,a),md_atom(a)
                end do
             end if
          end if

          call write_atomic_coordinates_log(97,icy,itlin,0,iswitch_opt)

!
! --- SCF ---
!

          select case(SYStype)
          case default

             call calc_eion

             call construct_strfac
             call construct_ps_local
             call construct_ps_pcc
             call destruct_strfac

             select case(pselect)
             case(2)
                call prep_ps_nloc2
             case(3)
                call prep_ps_nloc3
             case(5)
                call prep_ps_nloc_mr
             case(102)
                call prep_ps_nloc2
                call prepNzqr
                call prepQRijp102
             end select

          case(1)

             call calc_eion

             call construct_ps_local_mol
             call construct_ps_pcc_mol
             call prep_ps_nloc2_mol

          end select

          call calc_E_vdw_grimme( aa_atom )

          write(loop_info,'("( linmin:",i3,", cg:",i3," )")') itlin,icy
          call calc_scf( disp_switch, ierr, diter_opt, feps, loop_info, Etot )

          if ( ierr == -1 ) then
             if ( myrank == 0 ) write(*,*) "time limit !!!"
             exit opt_ion
          end if
          if ( ierr == -2 ) then
             if ( myrank == 0 ) write(*,*) "SCF is not converged"
          end if
          iter_final=ierr

          call calc_force( Natom, Force )

          if ( disp_switch_loc ) then
             write(*,'(1x,"# Force (total)")')
             if ( Natom <= 11 ) then
                do a=1,Natom
                   write(*,'(1x,i4,i3,3g21.12)') a,ki_atom(a),force(:,a)
                end do
             else
                do a=1,min(5,Natom)
                   write(*,'(1x,i4,i3,3g21.12)') a,ki_atom(a),force(:,a)
                end do
                write(*,'(1x,10x,".")')
                write(*,'(1x,10x,".")')
                write(*,'(1x,10x,".")')
                do a=Natom-5,Natom
                   write(*,'(1x,i4,i3,3g21.12)') a,ki_atom(a),force(:,a)
                end do
             end if
          end if

          grad=sum( Force(1:3,1:Natom)*hi(1:3,1:Natom) )

          Fmax=0.d0
          do a=1,Natom
             ss=Force(1,a)**2+Force(2,a)**2+Force(3,a)**2
             Fmax=max(Fmax,ss)
          end do
          Fmax=sqrt( Fmax )

          if ( (al(2)<Etot) .or. (grad<=0.d0) ) then
             ifar  = 1
             arr(1:3) = ar(1:3)
             ar(1)    = alpha1
             ar(2)    = Etot
             ar(3)    = grad
          else
             all(1:3) = al(1:3)
             al(1)    = alpha1
             al(2)    = Etot
             al(3)    = grad
          end if

!
! --- history ---
!

          nhist = nhist + 1

          Etot_hist(nhist)     = Etot
          alpha_hist(nhist)    = alpha1
          Fmax_hist(nhist)     = Fmax
          grad_hist(nhist)     = grad
          SCF_hist(nhist)      = iter_final
          ICY_hist(nhist)      = icy
          LIN_hist(nhist)      = itlin
          dmax_hist(nhist)     = &
               abs(alpha1)*sqrt(maxval(hi(1,:)**2+hi(2,:)**2+hi(3,:)**2))

          if ( disp_switch_loc ) then
             write(*,'(1x,3x,1x,a12,1x,a20,1x,3a13)') &
                  'alpha    ','Etot    ','grad    ','Fmax    ','dmax    '
             do i=nhist0,nhist
                write(*,'(1x,i3,1x,g13.6,1x,g20.10,1x,3g13.5,i5)') &
                     i-nhist0,alpha_hist(i),Etot_hist(i) &
                     ,grad_hist(i),Fmax_hist(i),dmax_hist(i),SCF_hist(i)
             end do
          end if

          call write_atomic_coordinates_log(97,icy,itlin,2,iswitch_opt)

!
! --- Convergence check 2 ---
!
          dif = Etot - Etsave

          if ( abs(grad) <= abs(grad0*decr) ) then
             icflag=0
             if ( disp_switch_loc ) then
                write(*,*) "--- Convergence achieved in LINMIN ---"
                write(*,*) "decr =",decr
                write(*,*) 'grad, grad0*decr =',grad,grad0*decr
             end if
          else if ( Fmax <= feps ) then
             icflag=0
             if ( disp_switch_loc ) then
                write(*,*) "--- Convergence achieved in LINMIN ---"
                write(*,*) "Fmax =",Fmax
             end if
          else if ( abs(dif) <= eeps ) then
             icflag=0
             if ( disp_switch_loc ) then
                write(*,*) "--- Convergence achieved in LINMIN ---"
                write(*,*) 'dif(Etot) =',Etot
             end if
          else if ( ddmax <= ddmin ) then
             icflag=1
             if ( disp_switch_loc ) then
                write(*,*) "--- ddmax is too small ---"
             end if
          else
             icflag=-1
          end if

          if ( Fmax < Fmax0 ) then
             Fmax0  = Fmax
             alpha0 = alpha1
             Etot0  = Etot
          end if

          if ( Etot < Etot0_2 ) then
             Fmax0_2  = Fmax
             alpha0_2 = alpha1
             Etot0_2  = Etot
          end if

          if ( icflag >= 0 ) then
             if ( Etot > Etot00 ) then
                if ( disp_switch_loc ) then
! why this kind of thing happens?
                   write(*,*) "can not exit linmin !!!"
                   write(*,*) "Etot, Etot00 =",Etot,Etot00
                end if
             else
                exit linmin
             end if
          end if

       end do linmin

       hi(1:3,1:Natom) = signh*hi(1:3,1:Natom)

       most0=1

!  Best structure on a line.

       if ( disp_switch_loc ) then
          write(*,*) 'Best structure on a line (see fort.197)'
       end if  
       call write_atomic_coordinates_log(197,icy,itlin,1,iswitch_opt)

!
! --- Convergence check 3 ---
!
       dif  = Etot-Etot_save
       dif0 = dif

       if ( Fmax <= feps ) then

          if ( disp_switch_loc ) then
             write(*,*) "Fmax,feps =",Fmax,feps
          end if
          exit opt_ion

       else if ( abs(dif) < eeps ) then

          if ( disp_switch_loc ) then
             write(*,*) "Etot,Etot_save =",Etot,Etot_save
          end if
          exit opt_ion

       end if

    end do opt_ion

999 continue

    if ( disp_switch_loc ) then
       write(*,*) "histry (Etot,Fmax,SCF)"
       do i=1,nhist
          write(*,'(1x,2i3,f16.8,f14.7,i6,f14.7)') ICY_hist(i) &
               ,LIN_hist(i),Etot_hist(i),Fmax_hist(i),SCF_hist(i),dmax_hist(i)
       end do
    end if

    call check_disp_switch( disp_switch_loc, 1 )
    call check_disp_switch( disp_switch, 0 )

    deallocate( Force )
    deallocate( aa_atom_0, gi, hi )

    return

  END SUBROUTINE atomopt


  SUBROUTINE get_min_parabola(g1,g2,xmin,emin)
    implicit none
    real(8),intent(IN)  :: g1(3),g2(3)
    real(8),intent(OUT) :: xmin,emin
    real(8) :: f1(3),f2(3)
    real(8) :: a,b,c,dl,d1,d2,x

    f1(1:3) = (/ g1(1), g1(2), -g1(3) /)
    f2(1:3) = (/ g2(1), g2(2), -g2(3) /)

    a=(f2(3)-f1(3))/(f2(1)-f1(1))
    dl=abs(f2(1)-f1(1))
    if ( abs(a) <= 0.0d0 ) then
       b=(f1(3)+f2(3))*0.5d0
       if ( b >= 0.0d0 ) then
          xmin=f2(1)-1.d6*dl
       else
          xmin=f1(1)+1.d6*dl
       end if
       emin=0.0d0
    else
       b=f1(3)-a*f1(1)
       xmin=-b/a
       d1=abs(xmin-f1(1))
       d2=abs(xmin-f2(1))
       if ( d1 < d2 ) then
          x=f1(1)
          c=f1(2)-0.5d0*a*x*x-b*x
       else
          x=f2(1)
          c=f2(2)-0.5d0*a*x*x-b*x
       end if
       x=xmin
       emin=0.5d0*a*x*x+b*x+c
    end if
    if ( disp_switch_loc ) then
       write(*,*) 'a   =',a
       write(*,*) 'xmin=',xmin
       write(*,*) 'emin=',emin
    end if
    return
  END SUBROUTINE get_min_parabola
  

  SUBROUTINE write_atomic_coordinates_log(unit,icy,itlin,flag,iswitch_opt)
    implicit none
    integer, intent(IN) :: unit,icy,itlin,flag,iswitch_opt

    if ( myrank == 0 ) call write_coordinates_atom( unit, 3 )

    if ( strlog /= 0 .and. flag == strlog ) then
       if ( icy == 0 ) then
          if ( myrank == 0 ) open(unit_strlog,file="strlog.dat")
          if ( iswitch_opt >= 2 ) return
          if ( myrank == 0 ) call write_coordinates_atom( unit_strlog, 1 )
       end if
       if ( myrank == 0 ) then
          write(unit_strlog,'("#_STRLOG_",a63," icy, itlin =",2(X,I3))') &
               repeat("-",63),icy,itlin
          call write_coordinates_atom( unit_strlog, 2 )
       end if
    end if

  END SUBROUTINE write_atomic_coordinates_log  

   
END MODULE atomopt_module
