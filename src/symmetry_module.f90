MODULE symmetry_module

  use parallel_module, np_2d => node_partition

  implicit none

  PRIVATE
  PUBLIC :: init_symmetry, read_symmetry, prep_symmetry, sym_rho, sym_force &
           ,nsym,isymmetry,rgb,construct_matrix_symmetry

  integer :: isymmetry
  character(30) :: file_symdat

  integer :: nsym, nnp
  integer,allocatable :: rga(:,:,:), pga(:,:)
  real(8),allocatable :: rgb(:,:,:), pgb(:,:)

  integer :: a1b,b1b,a2b,b2b,a3b,b3b
  integer :: ML,ML1,ML2,ML3
  real(8) :: dV
  real(8) :: aa(3,3), bb(3,3)

CONTAINS


  SUBROUTINE read_symmetry( rank, unit )
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i,ierr
    character(4) :: cbuf,ckey

    isymmetry   = 0
    file_symdat = ""

    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey == "ISYM" ) then
             backspace(unit)
             read(unit,*) cbuf, isymmetry, file_symdat
             exit
          end if
       end do
999    continue
       write(*,*) "isymmetry=",isymmetry
       write(*,*) "file_symdat=",file_symdat
    end if

    call MPI_BCAST(isymmetry,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(file_symdat,30,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

  END SUBROUTINE read_symmetry


  SUBROUTINE init_symmetry(Ngrid,dV_in,aa_in,bb_in,MI,Kion,asi)
    implicit none
    integer,intent(IN) :: Ngrid(0:3)
    real(8),intent(IN) :: dV_in,aa_in(3,3),bb_in(3,3)
    integer,intent(IN) :: MI,Kion(MI)
    real(8),intent(IN) :: asi(3,MI)

    if ( isymmetry == 0 ) return

    call write_border( 0," init_symmetry(start)")

    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)
    dV  = dV_in
    aa  = aa_in
    bb  = bb_in

    call input_symdat( MI, Kion, asi )

    call input_symdat_2( MI, Kion, asi )
!    call input_symdat_3( MI, Kion, asi )
!    call input_symdat_4( MI, Kion, asi )

    call write_border( 0," init_symmetry(end)")

  END SUBROUTINE init_symmetry

!--------1---------2---------3---------4---------5---------6---------7--
!
! "isymmetry>0" assume the atomic coordinates are given within [0,1]
! "isymmetry<0" assume the atomic coordinates are given within [-0.5,0.5]

  SUBROUTINE input_symdat( MI, Kion, asi )
    implicit none
    integer,intent(IN) :: MI
    integer,intent(IN) :: Kion(MI)
    real(8),intent(IN) :: asi(3,MI)
    integer,parameter :: u=2
    integer :: i,i1,i2,j,k,n,ierr,flag
    integer :: A11,A12,A13,A21,A22,A23,A31,A32,A33,gcounter
    integer :: AAA(9,19683)
    real(8) :: RR0(3),RR1(3),RR2(3),RRR
    real(8) :: R1q,R2q,R3q,R1p,R2p,R3p,R1s,R2s,R3s
    real(8) :: tmp2(3,3,2)
    real(8) :: c1
    real(8),allocatable :: XX(:,:), pga_tmp(:,:)
    character(12) :: label_sym

    if ( .not.(abs(isymmetry)==1 .or. abs(isymmetry)==2 .or. &
               abs(isymmetry)==3) )  return

    call write_border( 0, " input_symdat(start)" )

    if ( disp_switch_parallel ) then
       write(*,*) "isymmetry=",isymmetry
    end if

    nnp  = 1
    nsym = 0

!
! --- Read symmetry operation matrix ---
!
    if ( myrank == 0 ) then

       if ( abs(isymmetry) <= 2 ) open(u,file=file_symdat,status='old')

       select case( isymmetry )
!(1)
       case( 1, -1 )

          read(u,*) nsym, nnp
          allocate( rga(3,3,nsym),pga(3,nsym) ) ; rga=0 ; pga=0
          allocate( rgb(3,3,nsym),pgb(3,nsym) ) ; rgb=0 ; pgb=0
          do n=1,nsym
             read(u,*) ( (rga(i,j,n),j=1,3),i=1,3 ),pga(:,n)
          end do
!(2)
       case( 2, -2 )

          read(u,*) label_sym, nsym
          write(*,*) "label_sym= ",label_sym
          write(*,*) "nsym     = ",nsym
          allocate( rga(3,3,nsym),pga(3,nsym) ) ; rga=0 ; pga=0
          allocate( rgb(3,3,nsym),pgb(3,nsym) ) ; rgb=0 ; pgb=0
          do i=1,nsym,2
             j=min(i+1,nsym)
             n=j-i+1
             read(u,*) (((tmp2(i2,i1,j),i1=1,3),i2=1,3),j=1,n)
             do j=1,n
                rga(1:3,1:3,i+j-1)=nint( tmp2(1:3,1:3,j) )
             end do
          end do

          allocate( pga_tmp(3,nsym) ) ; pga_tmp=0.d0
          do i=1,nsym,3
             j=min(i+2,nsym)
             read(u,*) pga_tmp(1:3,i:j)
          end do
          c1=1.0d10
          do i=1,nsym
             do j=1,3
                if ( abs(pga_tmp(j,i)) > 1.0d-10 ) then
                   c1=min( c1, abs(pga_tmp(j,i)) )
                end if
             end do
          end do
          pga(1:3,1:nsym)=nint( pga_tmp(1:3,1:nsym)/c1 )
          nnp=nint( 1.0d0/c1 )
          deallocate( pga_tmp )
!(3)
       case( 3, -3 )

          allocate( XX(3,MI) ) ; XX=0.0d0

          do i=1,MI
             XX(1:3,i)=asi(1:3,i)
             do j=1,3
                do  while ( XX(j,i) >  0.5d0 )
                   XX(j,i)=XX(j,i)-1.0d0
                end do
                do  while ( XX(j,i) <= -0.5d0 )
                   XX(j,i)=XX(j,i)+1.0d0
                end do
             end do
          end do ! i
 
          do A11=-1,1
          do A12=-1,1
          do A13=-1,1
          do A21=-1,1
          do A22=-1,1
          do A23=-1,1
          do A31=-1,1
          do A32=-1,1
          do A33=-1,1

             gcounter=0

             do i=1,MI

                RR0(1:3) = XX(1:3,i)

                flag = 0

                do j=1,MI

                   if ( Kion(i) /= Kion(j) ) cycle

                   RR1(1:3) = XX(1:3,j)

                   RR2(1) = A11*RR1(1) + A12*RR1(2) + A13*RR1(3)
                   RR2(2) = A21*RR1(1) + A22*RR1(2) + A23*RR1(3)
                   RR2(3) = A31*RR1(1) + A32*RR1(2) + A33*RR1(3)

                   do k=1,3
                      do  while ( RR2(k) >  0.5d0 ) 
                         RR2(k)=RR2(k)-1.0d0
                      end do
                      do while ( RR2(k) <= -0.5d0 ) 
                         RR2(k)=RR2(k)+1.0d0
                      end do
                   end do ! k

                   RRR = sqrt( (RR2(1)-RR0(1))**2 &
                             + (RR2(2)-RR0(2))**2 + (RR2(3)-RR0(3))**2 )
                   if ( RRR < 1.d-6 ) then
                      flag=1
                   end if

                end do ! j

                gcounter = gcounter + flag

             end do ! i
           

             if ( gcounter == MI ) then

                R1q = sqrt( AA(1,1)**2 + AA(2,1)**2 + AA(3,1)**2 )
                R2q = sqrt( AA(1,2)**2 + AA(2,2)**2 + AA(3,2)**2 )
                R3q = sqrt( AA(1,3)**2 + AA(2,3)**2 + AA(3,3)**2 )

                R1s = real(A11)*1 + real(A12)*0 + real(A13)*0
                R2s = real(A21)*1 + real(A22)*0 + real(A23)*0
                R3s = real(A31)*1 + real(A32)*0 + real(A33)*0
                R1p = sqrt( (AA(1,1)*R1s + AA(1,2)*R2s + AA(1,3)*R3s)**2 &
                          + (AA(2,1)*R1s + AA(2,2)*R2s + AA(2,3)*R3s)**2 &
                          + (AA(3,1)*R1s + AA(3,2)*R2s + AA(3,3)*R3s)**2   )

                R1s = real(A11)*0 + real(A12)*1 + real(A13)*0
                R2s = real(A21)*0 + real(A22)*1 + real(A23)*0
                R3s = real(A31)*0 + real(A32)*1 + real(A33)*0
                R2p = sqrt( (AA(1,1)*R1s + AA(1,2)*R2s + AA(1,3)*R3s)**2 &
                          + (AA(2,1)*R1s + AA(2,2)*R2s + AA(2,3)*R3s)**2 &
                          + (AA(3,1)*R1s + AA(3,2)*R2s + AA(3,3)*R3s)**2   )

                R1s = real(A11)*0 + real(A12)*0 + real(A13)*1
                R2s = real(A21)*0 + real(A22)*0 + real(A23)*1
                R3s = real(A31)*0 + real(A32)*0 + real(A33)*1
                R3p = sqrt( (AA(1,1)*R1s + AA(1,2)*R2s + AA(1,3)*R3s)**2 &
                          + (AA(2,1)*R1s + AA(2,2)*R2s + AA(2,3)*R3s)**2 &
                          + (AA(3,1)*R1s + AA(3,2)*R2s + AA(3,3)*R3s)**2   )

                if ( abs(R1q-R1p) <= 1.d-6 .and. &
                     abs(R2q-R2p) <= 1.d-6 .and. &
                     abs(R3q-R3p) <= 1.d-6        ) then
                   nsym=nsym+1
                   AAA(1,nsym)=A11
                   AAA(2,nsym)=A12
                   AAA(3,nsym)=A13
                   AAA(4,nsym)=A21
                   AAA(5,nsym)=A22
                   AAA(6,nsym)=A23
                   AAA(7,nsym)=A31
                   AAA(8,nsym)=A32
                   AAA(9,nsym)=A33
                end if

             end if

          end do
          end do
          end do
          end do
          end do
          end do
          end do
          end do
          end do

          deallocate( XX )

          write(*,*) "------------------------------------------------------"
          write(*,'(I0,42X,I0)')  nsym, 1
          do i=1,nsym
             write(*,'(3(3I3,3X),6X,3I3)')  AAA(1:9,i),0,0,0
          end do !i
          write(*,*) "------------------------------------------------------"
          write(*,*) "BE CAREFUL!! Bubun heishin & Rasen are NOT considered."
          write(*,*) 'Intentional STOP, isymmetry=', isymmetry
          stop "stop@input_symdat"
!
!---
!
       end select

       close(u)

    end if

    call mpi_bcast(nsym,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(nnp ,1,mpi_integer,0,mpi_comm_world,ierr)

    if ( myrank /= 0 ) then
       allocate( rga(3,3,nsym),pga(3,nsym) ) ; rga=0 ; pga=0
       allocate( rgb(3,3,nsym),pgb(3,nsym) ) ; rgb=0 ; pgb=0
    end if

    call mpi_bcast(rga,9*nsym,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(pga,3*nsym,mpi_integer,0,mpi_comm_world,ierr)

    if ( disp_switch_parallel ) then
       write(*,*) "nsym,nnp=",nsym,nnp
    end if

    call write_border( 0, " input_symdat(end)" )

  END SUBROUTINE input_symdat


  SUBROUTINE chk_grp( rga, aa, bb )
    implicit none
    integer,intent(IN) :: rga(:,:,:)
    real(8),intent(IN) :: aa(3,3), bb(3,3)
    real(8) :: aa_inv(3,3),tmp0(3,3),tmp1(3,3)
    real(8),allocatable :: rtmp(:,:,:)
    integer :: n,ns
    ns=size(rga,3)
    allocate( rtmp(3,3,ns) ) ; rtmp=0.0d0
    aa_inv(:,:) = transpose( bb )/(2.0d0*acos(-1.0d0))
    do n=1,ns
       tmp0(:,:) = rga(:,:,n)
       tmp1(:,:) = matmul( tmp0, aa_inv )
       rtmp(:,:,n) = matmul( aa, tmp1 )
    end do
    do n=1,ns
       tmp0(:,:) = transpose( rtmp(:,:,n) )
       tmp1(:,:) = matmul( rtmp(:,:,n), tmp0 )
       write(*,*) n,sum(abs(tmp1))
    end do
    deallocate( rtmp )
    stop "stop@chk_grp"
  END SUBROUTINE chk_grp


  SUBROUTINE input_symdat_2( MI, Kion, asi )
    implicit none
    integer,intent(IN) :: MI, Kion(MI)
    real(8),intent(IN) :: asi(3,MI)
    real(8) :: pi2,c1,r0
    real(8) :: tmp(3,3),tmp1(3,3),tmp2(3,3,2),tmp0(3,3),tmp4(3)
    real(8) :: bb_inv(3,3),aa_inv(3,3)
    real(8),allocatable :: Rra(:,:)
    integer :: mtmp1(3,3),mtmp2(3,3),itmp1(3),itmp2(3)
    integer :: i,i1,i2,i3,isym,j,loop,m,m1,m2,n
    integer,allocatable :: itmp(:)

    call write_border( 0, " input_symdat_2(start)" )

    pi2 = 2.0d0*acos(-1.0d0)

! --- make rgb (symmetry operation matrix in bb-representation) ---
!
    bb_inv(1:3,1:3) = transpose( aa(1:3,1:3) )/pi2
    aa_inv(1:3,1:3) = transpose( bb(1:3,1:3) )/pi2
    tmp2(1:3,1:3,1) = matmul( bb_inv(1:3,1:3), aa(1:3,1:3) )
    tmp2(1:3,1:3,2) = matmul( aa_inv(1:3,1:3), bb(1:3,1:3) )
    do n=1,nsym
       tmp0(1:3,1:3)  = matmul( tmp2(1:3,1:3,1), rga(1:3,1:3,n) )
       tmp1(1:3,1:3)  = matmul( tmp0(1:3,1:3), tmp2(1:3,1:3,2) )
       rgb(1:3,1:3,n) = tmp1(1:3,1:3)
    end do

!
! --- Check the symmetry operations ---
!

! tr(rgb)=(rga)^-1
!
    do n=1,nsym
       tmp0(1:3,1:3) = transpose( rgb(1:3,1:3,n) )
       tmp1(1:3,1:3) = matmul( tmp0(1:3,1:3), rga(1:3,1:3,n) )
       itmp1(1)=nint( tmp1(1,1) )
       itmp1(2)=nint( tmp1(2,2) )
       itmp1(3)=nint( tmp1(3,3) )
       i1=sum(abs(nint(tmp1)))
       if ( .not.( all(itmp1==1) .and. i1==3 ) ) then
          write(*,*) "tr(rgb) is not the inverse of rga!"
          write(*,*) n,i1
          write(*,*) tmp1
          goto 900
       end if
    end do

! check duplication
!
    allocate( itmp(nsym) )
    do n=1,nsym
       j=0
       itmp=0
       do m=1,nsym
          mtmp1(1:3,1:3)=rga(1:3,1:3,n)-rga(1:3,1:3,m)
          itmp1(1:3)=pga(1:3,n)-pga(1:3,m)
          if ( all(mtmp1==0) .and. all(itmp1==0) ) then
             j=j+1
             itmp(j)=m
          end if
       end do
       if ( j /= 1 ) then
          write(*,*) "symmetry error (duplication)",j,n
          write(*,'(1x,i5,3x,9i3,3x,3i3)') n,rga(1:3,1:3,n),pga(1:3,n)
          do i=1,j
             write(*,'(1x,i5,3x,9i3,3x,3i3)') &
                  itmp(i),rga(1:3,1:3,itmp(i)),pga(1:3,itmp(i))
          end do
          goto 900
       end if
    end do
    deallocate( itmp )

! existence of the unit operator
!
    m1=0
    do n=1,nsym
       itmp1(1)=rga(1,1,n)
       itmp1(2)=rga(2,2,n)
       itmp1(3)=rga(3,3,n)
       i=sum(abs(rga(1:3,1:3,n)))
       if ( i==3 .and. all(itmp1==1) .and. all(pga(1:3,n)==0) ) then
          m1=m1+1
       end if
    end do
    if ( m1 /= 1 ) then
       write(*,*) "symmetry error (no unit operator)",m1
       goto 900
    end if

! (rga_1)(rga_2)=(rga)
!
    do n=1,nsym
       do m=1,nsym
          mtmp1(1:3,1:3)=matmul( rga(1:3,1:3,m),rga(1:3,1:3,n) )
          itmp1(1:3)=matmul( rga(1:3,1:3,m),pga(1:3,n) )+pga(1:3,m)
          do j=1,3
             do loop=1,10000
                if ( itmp1(j)<0 ) then
                   itmp1(j)=itmp1(j)+nnp
                else if ( itmp1(j)>=nnp ) then
                   itmp1(j)=itmp1(j)-nnp
                else
                   exit
                end if
             end do
          end do
          m2=0
          do i=1,nsym
             mtmp2(1:3,1:3)=mtmp1(1:3,1:3)-rga(1:3,1:3,i)
             itmp2(1:3)=itmp1(1:3)-pga(1:3,i)
             if ( all(mtmp2==0) .and. all(itmp2==0) ) then
                m2=m2+1
             end if
          end do ! i
          if ( m2 /= 1 ) then
             write(*,*) "symmetry error (this is not group)",m2,m,n
             write(*,'(1x,3(3i3,1x),2x,3i3)') mtmp1,itmp1
             goto 900
          end if
       end do ! m
    end do ! n

!
! --- Check the consistency of symmetry operations with atomic coordinates ---
!
    m1 = 1
    m2 = MI

    allocate( Rra(3,m1:m2) ) ; Rra=0.0d0

    c1=1.d0/dble(nnp)

    do isym=1,nsym

       do i=m1,m2

          Rra(1:3,i)=matmul( rga(1:3,1:3,isym),asi(1:3,i) )+c1*pga(1:3,isym)

          do j=1,3
             if ( abs(Rra(j,i)) >= 1.0d0 ) then
                n=int(Rra(j,i))
                Rra(j,i)=Rra(j,i)-n
             end if
          end do

          loop_j : do j=1,MI
             if( Kion(i) == Kion(j) )then
                do i3=-1,1
                do i2=-1,1
                do i1=-1,1
                   tmp4(1)=Rra(1,i)+i1
                   tmp4(2)=Rra(2,i)+i2
                   tmp4(3)=Rra(3,i)+i3
                   r0=sum( (tmp4(1:3)-asi(1:3,j))**2 )
                   if ( r0 < 1.d-10 ) then
                      Rra(1:3,i)=tmp4(1:3)
                      exit loop_j
                   end if
                end do
                end do
                end do
             end if
          end do loop_j
          if ( j > MI ) then
             write(*,*) "symmetry error (atomic coordinates)",isym,i
             goto 900
          end if

       end do ! i

    end do ! isym

    deallocate( Rra )

    call write_border( 0, " input_symdat_2(end)" )

    return

900 stop "stop@input_symdat_2"

  END SUBROUTINE input_symdat_2


  SUBROUTINE check_duplication(rga,pga)
    implicit none
    integer :: rga(:,:,:), pga(:,:)
    integer,allocatable :: itmp(:),ichk(:)
    integer :: n,j,m,i,mtmp1(3,3),itmp1(3)
    allocate( itmp(nsym) ) ; itmp=0
    allocate( ichk(nsym) ) ; ichk=0
    if ( myrank == 0 ) rewind 11
    do n=1,nsym
       ichk(n)=ichk(n)+1
       j=0
       itmp(:)=0
       do m=1,nsym
          if ( m == n ) cycle
          mtmp1(1:3,1:3)=rga(1:3,1:3,n)-rga(1:3,1:3,m)
          itmp1(1:3)=pga(1:3,n)-pga(1:3,m)
          if ( all(mtmp1==0) .and. all(itmp1==0) ) then
             j=j+1
             itmp(j)=m
             ichk(m)=ichk(m)+1
          end if
       end do ! m
       if ( j /= 0 ) then
          write(*,*) "symmetry error (duplication)",j,n
          write(*,'(1x,i5,3x,9i3,3x,3i3)') n,rga(1:3,1:3,n),pga(1:3,n)
          do i=1,j
             write(*,'(1x,i5,3x,9i3,3x,3i3)') &
                  itmp(i),rga(1:3,1:3,itmp(i)),pga(1:3,itmp(i))
          end do
       end if
       if ( ichk(n) == 1 ) then
          if ( myrank == 0 ) then
             write(11,'(1x,3(3i3,1x),2x,3i3)') &
                  ((rga(i,j,n),j=1,3),i=1,3),pga(1:3,n)
          end if
       end if
    end do ! n
    deallocate( ichk )
    deallocate( itmp )
  END SUBROUTINE check_duplication


  SUBROUTINE input_symdat_4( MI, Kion, asi )
    use aa_module, only: ax
    use bb_module, only: calc_bb
    implicit none
    integer,intent(IN) :: MI, Kion(MI)
    real(8),intent(IN) :: asi(3,MI)
    real(8) :: aa0(3,3),bb0(3,3),aa_inv(3,3),aa0_inv(3,3)
    real(8) :: pi2, tmp1(3,3),tmp2(3,3),tmp3(3,3),tmp4(3,3)
    integer :: isym, u, i, j

    pi2 = 2.0d0*acos(-1.0d0)

!    aa0(1:3,1) = (/ 0.0d0, 0.5d0, 0.5d0 /)
!    aa0(1:3,2) = (/ 0.5d0, 0.0d0, 0.5d0 /)
!    aa0(1:3,3) = (/ 0.5d0, 0.5d0, 0.0d0 /)
    aa0(1:3,1) = (/ 1.0d0, 0.0d0, 0.0d0 /)
    aa0(1:3,2) = (/ 0.0d0, 1.0d0, 0.0d0 /)
    aa0(1:3,3) = (/ 0.0d0, 0.0d0, 1.0d0 /)
    aa0=aa0*ax
    call calc_bb( aa0, bb0 )
    aa0_inv = transpose( bb0 )/pi2

    aa_inv = transpose( bb )/pi2

    u=10
    rewind u

    write(u,*) nsym,nnp

    do isym=1,nsym

       tmp1 = matmul( rga(:,:,isym), aa0_inv ) 
       tmp2 = matmul( aa0, tmp1 )

       tmp3 = matmul( aa_inv, tmp2 )
       tmp4 = matmul( tmp3, aa )

       write(*,*) "------------ isym=",isym
       write(*,*) tmp4,pga(:,isym)

       rga(:,:,isym) = nint( tmp4 )

       write(u,'(1x,3(3i3,1x),2x,3i3)') &
            ((rga(i,j,isym),j=1,3),i=1,3),pga(1:3,isym)

    end do ! isym
stop
  END SUBROUTINE input_symdat_4


  SUBROUTINE translation1(Rra)
    implicit none
    real(8),intent(INOUT) :: Rra(3)
    integer :: j,k
    do j=1,3
       if ( Rra(j) >= 1.0d0 ) then
          do k=1,1000
             Rra(j) = Rra(j) - k
             if ( 0.0d0 <= Rra(j) .and. Rra(j) < 1.0d0 ) exit
          end do
       else if ( Rra(j) < 0.0d0 ) then
          do k=1,1000
             Rra(j) = Rra(j) + k
             if ( 0.0d0 <= Rra(j) .and. Rra(j) < 1.0d0 ) exit
          end do
       end if
    end do
  END SUBROUTINE translation1

  SUBROUTINE translation2(Rra)
    implicit none
    real(8),intent(INOUT) :: Rra(3)
    integer :: j
    do j=1,3
       if ( Rra(j) >= 0.5d0 ) then
          Rra(j)=Rra(j)-1.0d0
       else if ( Rra(j) < -0.5d0 ) then
          Rra(j)=Rra(j)+1.0d0
       end if
    end do
  END SUBROUTINE translation2

  SUBROUTINE input_symdat_3( MI, Kion, asi )
    implicit none
    integer,intent(IN) :: MI, Kion(MI)
    real(8),intent(IN) :: asi(3,MI)
    real(8) :: tmp1(3), shift(3)
    real(8) :: fac,c1,r0
    real(8),allocatable :: Rra(:,:)
    integer :: i,isym,j,k,n,u,ierr,ierr2
    integer,allocatable :: ichk(:)

!
! --- make symmetry operation matrix ---
!

    allocate( ichk(MI) ) ; ichk=0

    fac = 4.0d0

    u = 20
    rewind u

    if ( myrank == 0 ) then
       write(*,*) "make symmetry opration matrix"
       rewind 10
       write(10,*) nsym,nint(fac)
    end if

    c1=1.d0/nnp
    allocate( Rra(3,MI) ) ; Rra=0.d0

    do isym=1,nsym

       ichk(:)=0
       ierr   =0
       ierr2  =0
       do i=1,MI

          Rra(1:3,i)=matmul( rga(1:3,1:3,isym), asi(1:3,i) )
          if ( isymmetry > 0 ) then
             call translation1( Rra(1,i) )
          else if ( isymmetry < 0 ) then
             call translation2( Rra(1,i) )
          end if

          if ( myrank == 0 ) then
             write(*,'(1x,2i6,2(3f10.5,2x))') isym,i, Rra(1:3,i), asi(1:3,i)
          end if

          do j=1,MI
             if ( Kion(i) /= Kion(j) ) cycle
             r0=sum( (Rra(:,i)-asi(:,j))**2 )
             if ( r0 < 1.d-10 ) then
                ichk(j)=ichk(j)+1
                exit
             end if
          end do ! j
          if ( j > MI ) then
             ierr=ierr+1
          end if
          if ( any( ichk >= 2 ) ) then
             ierr2=ierr2+1
          end if

       end do ! i

       if ( ierr2 /= 0 ) then
          do i=1,MI
             write(*,*) i,ichk(i)
          end do
          stop "stop@input_symda_3(1)"
       end if

       if ( ierr == 0 ) then
          if ( myrank == 0 ) then
             write(u,*) i,"okok"
             write(*,*) i,"okok"
             write(10,'(1x,3(3i3,1x),2x,3i3)') &
                  ((rga(i,j,isym),j=1,3),i=1,3),nint(fac*tmp1(1:3))
             write(*,*) "isym/nsym=",isym,nsym
          end if
       else
          write(*,*) "symmetry error (atomic coordinates)"
          write(*,*) "isym=",isym
          write(*,*) "ierr=",ierr
          stop "input_symdat_3"
          tmp1(1:3)=-0.25d0
          do i=1,MI
             Rra(1:3,i)=Rra(1:3,i)+tmp1(1:3)
             do j=1,3
                do k=1,1000
                   if ( 0.0d0 <= Rra(j,i) .and. Rra(j,i) < 1.0d0 ) exit
                   Rra(j,i) = Rra(j,i) - k
                end do
             end do
             write(*,'(1x,i6,2(3f10.5,2x))') i, Rra(1:3,i), asi(1:3,i)
             ierr=0
             do j=1,MI
                if( Kion(i) == Kion(j) )then
                   r0=sum( (Rra(:,i)-asi(:,j))**2 )
                   if( r0 < 1.d-10 )exit
                end if
             end do
             if ( j > MI ) then
                ierr=ierr+1
             end if
          end do ! i
          if ( ierr == 0 ) then
             write(u,*) i,"ok"
             write(*,*) i,"ok"
             if ( myrank == 0 ) then
                write(10,'(1x,3(3i3,1x),2x,3i3)') &
                     ((rga(i,j,isym),j=1,3),i=1,3),nint(fac*tmp1(1:3))
                write(*,*) "isym/nsym=",isym,nsym
             end if
          else
             write(u,*) i,"xx"
             write(*,*) i,"xx"
          end if
       end if

       pga(1:3,isym) = nint( fac*tmp1(1:3) )

    end do ! isym

    deallocate( Rra  )
    deallocate( ichk )

    call check_duplication(rga,pga)

    if ( myrank == 0 ) write(*,*) "see fort.10"
    stop "stop@input_symdat_3"

  END SUBROUTINE input_symdat_3

!--------1---------2---------3---------4---------5---------6---------7--

  SUBROUTINE prep_symmetry( Igrid )

    implicit none

    integer,intent(IN) :: Igrid(2,0:3)
    integer :: i,j,n,m,i1,i2,i3,isym,ierr,n1,n2,ML0
    integer :: j1,j2,j3,ii1,ii2,ii3,jj1,jj2,jj3,k1,k2,k3
    integer :: m1,m2,m3,ktmp(3,2),k
    integer,allocatable :: itmp(:),la(:,:),Rla(:,:),mtmp(:,:,:)
    integer,allocatable :: LLL2(:,:,:)
    real(8),allocatable :: Ra(:,:),Rra(:,:)
    real(8) :: r0,r1,r2,r3,tmp(3,3),c1

    if ( isymmetry == 0 ) return

    call write_border( 0, " prep_symmetry(start)" )

    n1  = idisp(myrank)+1
    n2  = idisp(myrank)+ircnt(myrank)
    ML0 = n2-n1+1
    c1  = 1.d0/dble(nnp)

    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)

    allocate( itmp(nsym) ) ; itmp=0
    allocate( LLL2(0:ML1-1,0:ML2-1,0:ML3-1) ) ; LLL2=0

    allocate( la(3,n1:n2)  ) ; la=0
    allocate( Rla(3,n1:n2) ) ; Rla=0

    n=-1
    i=0
    do j3=0,np_2d(3)-1
    do j2=0,np_2d(2)-1
    do j1=0,np_2d(1)-1
       n=n+1
       ii1=pinfo_grid(1,n) ; jj1=ii1+pinfo_grid(2,n)-1
       ii2=pinfo_grid(3,n) ; jj2=ii2+pinfo_grid(4,n)-1
       ii3=pinfo_grid(5,n) ; jj3=ii3+pinfo_grid(6,n)-1
       do i3=ii3,jj3
       do i2=ii2,jj2
       do i1=ii1,jj1
          i=i+1
          LLL2(i1,i2,i3)=i
       end do
       end do
       end do
    end do
    end do
    end do

    i=n1-1
    do i3=a3b,b3b
    do i2=a2b,b2b
    do i1=a1b,b1b
       i=i+1
       la(1,i)=i1
       la(2,i)=i2
       la(3,i)=i3
    end do
    end do
    end do

    if ( mod(ML1,nnp)/=0 .or. mod(ML2,nnp)/=0 .or. mod(ML2,nnp)/=0 ) then
       write(*,*) "ML1,ML2,ML3,nnp=",ML1,ML2,ML3,nnp
       write(*,*) "The grid is inconsistent with the symmetry"
       goto 900
    end if

    m1 = ML1/nnp
    m2 = ML2/nnp
    m3 = ML3/nnp

    do i=n1,n2

       itmp(1:nsym) = 0

       do isym=1,nsym

          Rla(1,i)=sum( rga(1,1:3,isym)*la(1:3,i) ) + m1*pga(1,isym)
          Rla(2,i)=sum( rga(2,1:3,isym)*la(1:3,i) ) + m2*pga(2,isym)
          Rla(3,i)=sum( rga(3,1:3,isym)*la(1:3,i) ) + m3*pga(3,isym)

          k1=Rla(1,i)/ML1 ; if ( Rla(1,i)<0 ) k1=(Rla(1,i)+1)/ML1-1
          k2=Rla(2,i)/ML2 ; if ( Rla(2,i)<0 ) k2=(Rla(2,i)+1)/ML2-1
          k3=Rla(3,i)/ML3 ; if ( Rla(3,i)<0 ) k3=(Rla(3,i)+1)/ML3-1
          Rla(1,i)=Rla(1,i)-k1*ML1
          Rla(2,i)=Rla(2,i)-k2*ML2
          Rla(3,i)=Rla(3,i)-k3*ML3
          i1=Rla(1,i)
          i2=Rla(2,i)
          i3=Rla(3,i)
          if ( i1<0 .or. i2<0 .or. i3<0 .or. &
               i1>=ML1 .or. i2>=ML2 .or. i3>=ML3 ) then
             write(*,*) i1,i2,i3
             goto 900
          end if

          itmp(isym)=LLL2(i1,i2,i3)

       end do ! isym

    end do ! i

    deallocate( Rla )
    deallocate( la )
    deallocate( LLL2 )
    deallocate( itmp )

    call write_border( 0, " prep_symmetry(end)" )

    return

900 stop "stop@prep_symmetry"

  END SUBROUTINE prep_symmetry

!--------1---------2---------3---------4---------5---------6---------7--

  SUBROUTINE sym_rho( n1,n2,n3,MSP_0,MSP_1,rho )
    implicit none
    integer,intent(IN) :: n1,n2,n3,MSP_0,MSP_1
    real(8),intent(INOUT) :: rho(n1:n2,n3)
    integer :: i,j,m,n,ierr,ispin,m1,m2,m3
    integer :: k1,k2,k3,i1,i2,i3,j1,j2,j3,l1,l2,l3,isym
    real(8) :: rho_i,c,c0,fac
    real(8),allocatable :: rho_tmp(:),rho3(:,:,:),work3(:,:,:)

    if ( isymmetry == 0 ) return

    call write_border( 1, " sym_rho(start)" )

    allocate( rho3(0:ML1-1,0:ML2-1,0:ML3-1)  ) ; rho3=0.0d0
    allocate( work3(0:ML1-1,0:ML2-1,0:ML3-1) ) ; work3=0.0d0

    m1  = ML1/nnp
    m2  = ML2/nnp
    m3  = ML3/nnp
    fac = 1.0d0/nsym

    do ispin=MSP_0,MSP_1

       work3(:,:,:)=0.0d0

       i=n1-1
       do i3=a3b,b3b
       do i2=a2b,b2b
       do i1=a1b,b1b
          i=i+1
          work3(i1,i2,i3)=rho(i,ispin)
       end do
       end do
       end do

       call mpi_allreduce(work3,rho3,ML,mpi_real8,mpi_sum,comm_grid,ierr)

       rho(n1:n2,ispin)=0.0d0

       i=n1-1
       do i3=a3b,b3b
       do i2=a2b,b2b
       do i1=a1b,b1b
          i=i+1
          do isym=1,nsym
             j1=i1-m1*pga(1,isym)
             j2=i2-m2*pga(2,isym)
             j3=i3-m3*pga(3,isym)
             l1=nint(rgb(1,1,isym)*j1+rgb(2,1,isym)*j2+rgb(3,1,isym)*j3)
             l2=nint(rgb(1,2,isym)*j1+rgb(2,2,isym)*j2+rgb(3,2,isym)*j3)
             l3=nint(rgb(1,3,isym)*j1+rgb(2,3,isym)*j2+rgb(3,3,isym)*j3)
             k1=l1/ML1 ; if ( l1<0 ) k1=(l1+1)/ML1-1
             k2=l2/ML2 ; if ( l2<0 ) k2=(l2+1)/ML2-1
             k3=l3/ML3 ; if ( l3<0 ) k3=(l3+1)/ML3-1
             l1=l1-k1*ML1
             l2=l2-k2*ML2
             l3=l3-k3*ML3
             rho(i,ispin) = rho(i,ispin) + fac*rho3(l1,l2,l3)
          end do
       end do
       end do
       end do

       c0=sum(rho(n1:n2,ispin))*dV
       call mpi_allreduce(c0,c,1,mpi_real8,mpi_sum,comm_grid,ierr)
!       if ( disp_switch_parallel ) write(*,*) "check sum(rho)=",c

    end do ! ispin

    deallocate( work3 )
    deallocate( rho3  )

    call write_border( 1, " sym_rho(end)" )

    return

  END SUBROUTINE sym_rho

!--------1---------2---------3---------4---------5---------6---------7--

  SUBROUTINE sym_force( MI, Kion, asi, force )
    implicit none
    integer,intent(IN) :: MI, Kion(MI)
    real(8),intent(IN) :: asi(3,MI)
    real(8),intent(INOUT) :: force(3,MI)
    integer :: a,b,isym,j,m,i1,i2,i3,n
    real(8) :: atmp(3),c1,c2,c3,r,r0,faa(3),tmp1(3)
    real(8),allocatable :: ftmp(:,:)

    if ( isymmetry /= 1 ) return

    call write_border( 1, " sym_force(start)" )

    c1 = 1.0d0/dble(nnp)
    c2 = 1.0d0/dble(nsym)
    c3 = 1.0d0/( 2.0d0*acos(-1.0d0) )

    allocate( ftmp(3,MI) )

    ftmp(:,:) = force(:,:)

    do a=1,MI
       force(1:3,a)=matmul( ftmp(1:3,a),bb(1:3,1:3) )*c3
    end do

    ftmp(:,:) = 0.0d0

    do a=1,MI

       do isym=1,nsym

          atmp(:) = matmul( rga(:,:,isym),asi(:,a) )+c1*pga(:,isym)

          do j=1,3
             if ( abs(atmp(j)) >= 1.0d0 ) then
                n=int(atmp(j))
                atmp(j)=atmp(j)-n
             end if
          end do

          loop_b : do b=1,MI

             if ( Kion(a) == Kion(b) ) then

                do i3=-1,1
                do i2=-1,1
                do i1=-1,1
                   tmp1(1)=atmp(1)+i1
                   tmp1(2)=atmp(2)+i2
                   tmp1(3)=atmp(3)+i3
                   r0=sum( (tmp1(1:3)-asi(1:3,b))**2 )
                   if ( r0 < 1.d-10 ) then
                      ftmp(1:3,b) = ftmp(1:3,b) &
                           + matmul( rga(:,:,isym),force(1:3,a) )
                      exit loop_b
                   end if
                end do
                end do
                end do
             end if

          end do loop_b

          if ( b > MI ) then
             write(*,*) "symmetry error (sym_force)",isym,a
             stop "stop@sym_force"
          end if

       end do ! isym

    end do ! a

    do a=1,MI
       force(1:3,a) = matmul( aa(1:3,1:3), ftmp(1:3,a) )*c2
    end do

    deallocate( ftmp )

    call write_border( 1, " sym_force(end)" )

    return
  END SUBROUTINE sym_force

!--------1---------2---------3---------4---------5---------6---------7--

  SUBROUTINE chk_sym_mat(ng,np,mat,vec)

    implicit none

    integer,intent(IN) :: ng,np
    integer,intent(IN) :: mat(3,3,ng),vec(3,ng)
    integer :: n,j,m,loop,i,m1,m2
    integer :: mtmp1(3,3),itmp1(3),mtmp2(3,3),itmp2(3)
    integer,allocatable :: itmp(:)

! check duplication
!
    allocate( itmp(ng) )

    do n=1,ng
       j=0
       itmp=0
       do m=1,ng
          mtmp1(1:3,1:3)=mat(1:3,1:3,n)-mat(1:3,1:3,m)
          itmp1(1:3)=vec(1:3,n)-vec(1:3,m)
          if ( all(mtmp1==0) .and. all(itmp1==0) ) then
             j=j+1
             itmp(j)=m
          end if
       end do
       if ( j/=1 ) then
          write(*,*) "symmetry error (duplication)",j,n
          write(*,'(1x,i5,3x,9i3,3x,3i3)') n,mat(1:3,1:3,n),vec(1:3,n)
          do i=1,j
             write(*,'(1x,i5,3x,9i3,3x,3i3)') &
                  itmp(i),mat(1:3,1:3,itmp(i)),vec(1:3,itmp(i))
          end do
          goto 900
       end if
    end do
    deallocate( itmp )

! existence of the unit operator
!
    m1=0
    do n=1,ng
       itmp1(1)=mat(1,1,n)
       itmp1(2)=mat(2,2,n)
       itmp1(3)=mat(3,3,n)
       i=sum(abs(mat(1:3,1:3,n)))
       if ( i==3 .and. all(itmp1==1) .and. all(vec(1:3,n)==0) ) then
          m1=m1+1
       end if
    end do
    if ( m1/=1 ) then
       write(*,*) "symmetry error (no unit operator)",m1
       goto 900
    end if

! (mat_1)(mat_2)=(mat)
!
    do n=1,ng
       do m=1,ng
          mtmp1(1:3,1:3)=matmul( mat(1:3,1:3,m),mat(1:3,1:3,n) )
          itmp1(1:3)=matmul( mat(1:3,1:3,m),vec(1:3,n) )+vec(1:3,m)
          do j=1,3
             do loop=1,10000
                if ( itmp1(j)<0 ) then
                   itmp1(j)=itmp1(j)+np
                else if ( itmp1(j)>=np ) then
                   itmp1(j)=itmp1(j)-np
                else
                   exit
                end if
             end do
          end do
          m2=0
          do i=1,ng
             mtmp2(1:3,1:3)=mtmp1(1:3,1:3)-mat(1:3,1:3,i)
             itmp2(1:3)=itmp1(1:3)-vec(1:3,i)
             if ( all(mtmp2==0) .and. all(itmp2==0) ) then
                m2=m2+1
             end if
          end do
          if ( m2/=1 ) then
             write(*,*) "symmetry error (this is not group)",m2,m,n
             goto 900
          end if
       end do ! m
    end do ! n

    return

900 stop "stop@chk_sym_mat"

  END SUBROUTINE chk_sym_mat


  SUBROUTINE construct_matrix_symmetry( isym, ML, SymMat )
    implicit none
    integer,intent(IN)  :: isym, ML
    real(8),intent(OUT) :: SymMat(ML,ML)
    integer :: i,j,n,m,i1,i2,i3,ierr,n1,n2,ML0
    integer :: j1,j2,j3,ii1,ii2,ii3,jj1,jj2,jj3,k1,k2,k3
    integer :: m1,m2,m3,ktmp(3,2),k
    integer,allocatable :: LLL2(:,:,:)
    real(8) :: la(3),Rla(3)
    real(8) :: r0,r1,r2,r3,tmp(3,3),c1

    if ( isymmetry == 0 ) return

    SymMat(:,:) = 0.0d0

    n1  = idisp(myrank)+1
    n2  = idisp(myrank)+ircnt(myrank)
    ML0 = n2-n1+1
    c1  = 1.d0/dble(nnp)

    allocate( LLL2(0:ML1-1,0:ML2-1,0:ML3-1) ) ; LLL2=0

    n=-1
    i=0
    do j3=0,np_2d(3)-1
    do j2=0,np_2d(2)-1
    do j1=0,np_2d(1)-1
       n=n+1
       ii1=pinfo_grid(1,n) ; jj1=ii1+pinfo_grid(2,n)-1
       ii2=pinfo_grid(3,n) ; jj2=ii2+pinfo_grid(4,n)-1
       ii3=pinfo_grid(5,n) ; jj3=ii3+pinfo_grid(6,n)-1
       do i3=ii3,jj3
       do i2=ii2,jj2
       do i1=ii1,jj1
          i=i+1
          LLL2(i1,i2,i3)=i
       end do
       end do
       end do
    end do
    end do
    end do

    m1 = ML1/nnp
    m2 = ML2/nnp
    m3 = ML3/nnp

    do i3=0,ML3-1
    do i2=0,ML2-1
    do i1=0,ML1-1

       i = 1 + i1 + i2*ML1 + i3*ML1*ML2

       la(1) = i1
       la(2) = i2
       la(3) = i3

       Rla(1)=sum( rga(1,1:3,isym)*la(1:3) ) + m1*pga(1,isym)
       Rla(2)=sum( rga(2,1:3,isym)*la(1:3) ) + m2*pga(2,isym)
       Rla(3)=sum( rga(3,1:3,isym)*la(1:3) ) + m3*pga(3,isym)

       k1=Rla(1)/ML1 ; if ( Rla(1)<0 ) k1=(Rla(1)+1)/ML1-1
       k2=Rla(2)/ML2 ; if ( Rla(2)<0 ) k2=(Rla(2)+1)/ML2-1
       k3=Rla(3)/ML3 ; if ( Rla(3)<0 ) k3=(Rla(3)+1)/ML3-1
       Rla(1)=Rla(1)-k1*ML1
       Rla(2)=Rla(2)-k2*ML2
       Rla(3)=Rla(3)-k3*ML3
       j1=Rla(1)
       j2=Rla(2)
       j3=Rla(3)

       j = LLL2(j1,j2,j3)

       SymMat(j,i) = 1.0d0

    end do ! i1
    end do ! i2
    end do ! i3

    deallocate( LLL2 )

    return

  END SUBROUTINE construct_matrix_symmetry


END MODULE symmetry_module
