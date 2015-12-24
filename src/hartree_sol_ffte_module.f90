MODULE hartree_sol_ffte_module

  use hartree_variables, only: E_hartree, Vh
  use bb_module, only: bb
  use ffte_sub_module, only: zwork1_ffte, zwork2_ffte, npuz, npuy, npux &
                            ,comm_fftx, comm_ffty, comm_fftz
  use rgrid_module, only: Ngrid, Igrid, dV
  use ggrid_module, only: NGgrid, LLG, construct_ggrid, destruct_ggrid
  use parallel_module
  use watch_module

  implicit none

  PRIVATE
  PUBLIC :: calc_hartree_sol_ffte

  logical :: first_time=.true.
  integer :: NGHT
  integer,allocatable :: LGHT(:,:)
  integer,allocatable :: IGHT(:,:)
  real(8),allocatable :: GGHT(:)

CONTAINS


  SUBROUTINE calc_hartree_sol_ffte(n1,n2,n3,rho)
    implicit none
    integer,intent(IN) :: n1,n2,n3
    real(8),intent(IN) :: rho(n1:n2,n3)
    integer :: i,i1,i2,i3,ierr,ispin,n
    real(8) :: Eh0,pi4,g2,ctt(0:5),ett(0:5)
    complex(8),parameter :: z0=(0.d0,0.d0)
    integer :: ML1,ML2,ML3,ML
    integer :: MG,ML_0,ML_1,a1b,b1b,a2b,b2b,a3b,b3b,ab1,ab12
    integer :: MG1,MG2,MG3,NG1,NG2,NG3,np1,np2,np3

    call write_border( 1, " calc_hartree_sol_ffte(start)" )

    pi4 = 4.d0*acos(-1.d0)

    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)
    MG  = NGgrid(0)
    ML_0= Igrid(1,0)
    ML_1= Igrid(2,0)
    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)
    ab1 = (b1b-a1b+1)
    ab12= (b1b-a1b+1)*(b2b-a2b+1)

    if ( first_time ) then
       call construct_Ggrid(0)
       n=0
       do i=1,NGgrid(0)
          i1=mod( Ngrid(1)+LLG(1,i), Ngrid(1) )
          i2=mod( Ngrid(2)+LLG(2,i), Ngrid(2) )
          i3=mod( Ngrid(3)+LLG(3,i), Ngrid(3) )
          if ( all(LLG(1:3,i)==0) ) cycle
          if ( a2b <= i2 .and. i2 <= b2b .and. a3b <= i3 .and. i3 <= b3b ) then
             n=n+1
          end if
       end do
       allocate( LGHT(3,n) ) ; LGHT=0
       allocate( IGHT(3,n) ) ; IGHT=0
       allocate( GGHT(n)   ) ; GGHT=0.0d0
       n=0
       do i=1,NGgrid(0)
          i1=mod( Ngrid(1)+LLG(1,i), Ngrid(1) )
          i2=mod( Ngrid(2)+LLG(2,i), Ngrid(2) )
          i3=mod( Ngrid(3)+LLG(3,i), Ngrid(3) )
          if ( all(LLG(1:3,i)==0) ) cycle
          if ( a2b <= i2 .and. i2 <= b2b .and. a3b <= i3 .and. i3 <= b3b ) then
             n=n+1
             LGHT(1,n)=i1
             LGHT(2,n)=i2
             LGHT(3,n)=i3
             g2=( bb(1,1)*LLG(1,i)+bb(1,2)*LLG(2,i)+bb(1,3)*LLG(3,i) )**2 &
               +( bb(2,1)*LLG(1,i)+bb(2,2)*LLG(2,i)+bb(2,3)*LLG(3,i) )**2 &
               +( bb(3,1)*LLG(1,i)+bb(3,2)*LLG(2,i)+bb(3,3)*LLG(3,i) )**2
             GGHT(n)=pi4/g2
          end if
       end do
       NGHT=n
       call destruct_Ggrid
       first_time=.false.
    end if
    
    ctt(:)=0.d0
    ett(:)=0.d0

    call watch(ctt(0),ett(0))

    zwork1_ffte(:,:,:)=z0
    do ispin=1,n3
!$OMP parallel do collapse(3) private(i)
       do i3=a3b,b3b
       do i2=a2b,b2b
       do i1=a1b,b1b
          i=ML_0+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
          zwork1_ffte(i1,i2,i3)=zwork1_ffte(i1,i2,i3)+rho(i,ispin)
       end do
       end do
       end do
!$OMP end parallel do
    end do

    call mpi_allreduce(zwork1_ffte,zwork2_ffte,ML1*(b2b-a2b+1)*(b3b-a3b+1) &
         ,mpi_complex16,mpi_sum,comm_fftx,ierr)

    call watch(ctt(1),ett(1))

    call pzfft3dv(zwork2_ffte,zwork1_ffte,ML1,ML2,ML3 &
         ,comm_ffty,comm_fftz,npuy,npuz,-1)

    call watch(ctt(2),ett(2))

    zwork2_ffte(:,:,:)=z0
    do i=1,NGHT
       i1=LGHT(1,i)
       i2=LGHT(2,i)
       i3=LGHT(3,i)
       zwork2_ffte(i1,i2,i3)=zwork1_ffte(i1,i2,i3)*GGHT(i)
    end do

    call watch(ctt(3),ett(3))

    call pzfft3dv(zwork2_ffte,zwork1_ffte,ML1,ML2,ML3 &
         ,comm_ffty,comm_fftz,npuy,npuz,1)

    call watch(ctt(4),ett(4))

!$OMP parallel do collapse(3) private(i)
    do i3=a3b,b3b
    do i2=a2b,b2b
    do i1=a1b,b1b
       i=ML_0+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
       Vh(i)=real( zwork1_ffte(i1,i2,i3) )
    end do
    end do
    end do
!$OMP end parallel do

    Eh0=0.d0
    do ispin=1,n3
!$OMP parallel do collapse(3) private(i) reduction(+:Eh0)
       do i3=a3b,b3b
       do i2=a2b,b2b
       do i1=a1b,b1b
          i=ML_0+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
          Eh0 = Eh0 + real( zwork1_ffte(i1,i2,i3) )*rho(i,ispin)
       end do
       end do
       end do
!$OMP end parallel do
    end do
    Eh0=0.5d0*Eh0*dV
    call mpi_allreduce(Eh0,E_hartree,1,mpi_real8,mpi_sum,comm_grid,ierr)

    call watch(ctt(5),ett(5))

!    if ( disp_switch_parallel ) then
!       write(*,*) "time(hatree1_ffte)=",ctt(1)-ctt(0),ett(1)-ett(0)
!       write(*,*) "time(hatree2_ffte)=",ctt(2)-ctt(1),ett(2)-ett(1)
!       write(*,*) "time(hatree3_ffte)=",ctt(3)-ctt(2),ett(3)-ett(2)
!       write(*,*) "time(hatree4_ffte)=",ctt(4)-ctt(3),ett(4)-ett(3)
!       write(*,*) "time(hatree5_ffte)=",ctt(5)-ctt(4),ett(5)-ett(4)
!    end if

    call write_border( 1, " calc_hartree_sol_ffte(end)" )

  END SUBROUTINE calc_hartree_sol_ffte


END MODULE hartree_sol_ffte_module
