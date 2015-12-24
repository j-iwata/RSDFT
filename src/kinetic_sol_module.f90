MODULE kinetic_sol_module

!$  use omp_lib
  use rgrid_module
  use omp_variables
  use bc_module, only: www, bcset_1, bcset_3
  use kinetic_variables, only: coef_lap0, coef_lap, zcoef_kin, coef_nab &
       ,flag_nab, flag_n12, flag_n23, flag_n31, const_k2, ggg, Md, wk
  use watch_module, only: watchb_omp, time_kine, time_bcfd

  implicit none

  PRIVATE
  PUBLIC :: op_kinetic_sol
  PUBLIC :: construct_matrix_kinetic_sol

CONTAINS


  SUBROUTINE op_kinetic_sol(k,tpsi,htpsi,n1,n2,ib1,ib2)
    implicit none
    integer,intent(IN) :: k,n1,n2,ib1,ib2
#ifdef _DRSDFT_
    real(8),intent(IN)    ::  tpsi(n1:n2,ib1:ib2)
    real(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
    real(8),parameter :: zero=0.d0
#else
    complex(8),intent(IN)    ::  tpsi(n1:n2,ib1:ib2)
    complex(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
    complex(8),parameter :: zero=(0.d0,0.d0)
#endif
    integer :: i,ib,i1,i2,i3,nb,m,n,j
    integer :: a1,a2,a3,b1,b2,b3,p,mm,nn
    integer :: a1b,b1b,a2b,b2b,a3b,b3b,ab1,ab12
    real(8) :: c,d,et0,et1, ttmp(2)
    integer,allocatable :: ic(:)
    integer :: a1b_omp,b1b_omp,a2b_omp,b2b_omp,a3b_omp,b3b_omp,n1_omp,n2_omp
    integer :: ib1_omp,ib2_omp,nb_omp

    !call watchb_omp( ttmp )
!$OMP master
    time_bcfd(:,:)=0.0d0
!$OMP end master

    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)
    ab1 = (b1b-a1b+1)
    ab12= (b1b-a1b+1)*(b2b-a2b+1)

    nb = ib2-ib1+1

!!$OMP parallel private(a3b_omp,b3b_omp,a2b_omp,b2b_omp,a1b_omp,b1b_omp &
!!$OMP                 ,n1_omp,n2_omp,j,p,d,mm,c)

    mm=0
!$  mm=omp_get_thread_num()

    n1_omp = Igrid_omp(1,0,mm)
    n2_omp = Igrid_omp(2,0,mm)

    a1b_omp = Igrid_omp(1,1,mm)
    b1b_omp = Igrid_omp(2,1,mm)
    a2b_omp = Igrid_omp(1,2,mm)
    b2b_omp = Igrid_omp(2,2,mm)
    a3b_omp = Igrid_omp(1,3,mm)
    b3b_omp = Igrid_omp(2,3,mm)

    c=coef_lap0+const_k2(k)
    do ib=ib1,ib2
       do i=n1_omp,n2_omp
          htpsi(i,ib) = c*tpsi(i,ib)
       end do
    end do

    !call watchb_omp( ttmp, time_kine(1,1) )

    do ib=ib1,ib2
       do i3=a3b_omp,b3b_omp
       do i2=a2b_omp,b2b_omp
       do i1=a1b_omp,b1b_omp
          j=n1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
          www(i1,i2,i3,ib-ib1+1) = tpsi(j,ib)
       end do
       end do
       end do
    end do

!$OMP barrier
    !call watchb_omp( ttmp, time_kine(1,2) )

    call bcset_3(1,nb,Md,0)

!$OMP barrier
    !call watchb_omp( ttmp, time_kine(1,3) )

    if ( flag_nab ) then

       do ib=ib1,ib2

          p = ib-ib1+1

          do m=1,Md

             do i3=a3b_omp,b3b_omp
             do i2=a2b_omp,b2b_omp
             do i1=a1b_omp,b1b_omp
                j=n1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
                htpsi(j,ib)=htpsi(j,ib) &
                     +zcoef_kin(1,m,k) *www(i1+m,i2,i3,p) &
               +conjg(zcoef_kin(1,m,k))*www(i1-m,i2,i3,p) &
                     +zcoef_kin(2,m,k) *www(i1,i2+m,i3,p) &
               +conjg(zcoef_kin(2,m,k))*www(i1,i2-m,i3,p) &
                     +zcoef_kin(3,m,k) *www(i1,i2,i3+m,p) &
               +conjg(zcoef_kin(3,m,k))*www(i1,i2,i3-m,p)   
             end do
             end do
             end do

          end do ! m

       end do ! ib

    else

       do ib=ib1,ib2

          p=ib-ib1+1

          do m=1,Md

             do i3=a3b_omp,b3b_omp
             do i2=a2b_omp,b2b_omp
             do i1=a1b_omp,b1b_omp
                j=n1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
                htpsi(j,ib)=htpsi(j,ib) &
                  +coef_lap(1,m)*( www(i1-m,i2,i3,p)+www(i1+m,i2,i3,p) ) &
                  +coef_lap(2,m)*( www(i1,i2-m,i3,p)+www(i1,i2+m,i3,p) ) &
                  +coef_lap(3,m)*( www(i1,i2,i3-m,p)+www(i1,i2,i3+m,p) )
             end do
             end do
             end do

          end do ! m

       end do ! ib

    end if

!$OMP barrier
    !call watchb_omp( ttmp, time_kine(1,4) )

    if ( flag_n12 .or. flag_n23 .or. flag_n31 ) then

!$OMP workshare
       wk=www
!$OMP end workshare

       if ( flag_n12 ) then

          do n=1,nb
             do i3=a3b_omp,b3b_omp
             do i2=a2b_omp,b2b_omp
             do i1=a1b_omp,b1b_omp
                www(i1,i2,i3,n)=zero
             end do
             end do
             end do
          end do
          do n=1,nb
             do m=1,Md
                d=coef_nab(1,m)
                do i3=a3b_omp,b3b_omp
                do i2=a2b_omp,b2b_omp
                do i1=a1b_omp,b1b_omp
                   www(i1,i2,i3,n)=www(i1,i2,i3,n) &
                        -d*( wk(i1-m,i2,i3,n)-wk(i1+m,i2,i3,n) )
                end do
                end do
                end do
             end do
          end do

!$OMP barrier
          call bcset_1(1,nb,Md,3)
!$OMP barrier

          do ib=ib1,ib2
             p=ib-ib1+1
             do m=1,Md
                d=-ggg(4)*coef_nab(2,m)
                do i3=a3b_omp,b3b_omp
                do i2=a2b_omp,b2b_omp
                do i1=a1b_omp,b1b_omp
                   j=n1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
                   htpsi(j,ib)=htpsi(j,ib) &
                        -d*(www(i1,i2-m,i3,p)-www(i1,i2+m,i3,p))
                end do
                end do
                end do
             end do
          end do

!$OMP barrier

       end if

       if ( flag_n23 ) then

          do n=1,nb
             do i3=a3b_omp,b3b_omp
             do i2=a2b_omp,b2b_omp
             do i1=a1b_omp,b1b_omp
                www(i1,i2,i3,n)=zero
             end do
             end do
             end do
          end do

          do n=1,nb
             do m=1,Md
                d=coef_nab(2,m)
                do i3=a3b_omp,b3b_omp
                do i2=a2b_omp,b2b_omp
                do i1=a1b_omp,b1b_omp
                   www(i1,i2,i3,n)=www(i1,i2,i3,n) &
                        -d*( wk(i1,i2-m,i3,n)-wk(i1,i2+m,i3,n) )
                end do
                end do
                end do
             end do
          end do

!$OMP barrier
          call bcset_1(1,nb,Md,5)
!$OMP barrier

          do ib=ib1,ib2
             p=ib-ib1+1
             do m=1,Md
                d=-ggg(5)*coef_nab(3,m)
                do i3=a3b_omp,b3b_omp
                do i2=a2b_omp,b2b_omp
                do i1=a1b_omp,b1b_omp
                   j=n1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
                   htpsi(j,ib)=htpsi(j,ib) &
                        -d*(www(i1,i2,i3-m,p)-www(i1,i2,i3+m,p))
                end do
                end do
                end do
             end do
          end do

!$OMP barrier

       end if

       if ( flag_n31 ) then

          do n=1,nb
             do i3=a3b_omp,b3b_omp
             do i2=a2b_omp,b2b_omp
             do i1=a1b_omp,b1b_omp
                www(i1,i2,i3,n)=zero
             end do
             end do
             end do
          end do

          do n=1,nb
             do m=1,Md
                d=coef_nab(3,m)
                do i3=a3b_omp,b3b_omp
                do i2=a2b_omp,b2b_omp
                do i1=a1b_omp,b1b_omp
                   www(i1,i2,i3,n)=www(i1,i2,i3,n) &
                        -d*( wk(i1,i2,i3-m,n)-wk(i1,i2,i3+m,n) )
                end do
                end do
                end do
             end do
          end do

!$OMP barrier
          call bcset_1(1,nb,Md,1)
!$OMP barrier

          do ib=ib1,ib2
             p=ib-ib1+1
             do m=1,Md
                d=-ggg(6)*coef_nab(1,m)
                do i3=a3b_omp,b3b_omp
                do i2=a2b_omp,b2b_omp
                do i1=a1b_omp,b1b_omp
                   j=n1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
                   htpsi(j,ib)=htpsi(j,ib) &
                        -d*(www(i1-m,i2,i3,p)-www(i1+m,i2,i3,p))
                end do
                end do
                end do
             end do
          end do

       end if

    end if

    !call watchb_omp( ttmp, time_kine(1,5) )

!$OMP master
    time_kine(1:2,6:11) = time_kine(1:2,6:11) + time_bcfd(1:2,1:6)
!$OMP end master

  END SUBROUTINE op_kinetic_sol


  SUBROUTINE construct_matrix_kinetic_sol( k, ML, Hmat )
    implicit none
    integer,intent(IN) :: k,ML
#ifdef _DRSDFT_
    real(8),intent(INOUT) :: Hmat(ML,ML)
#else
    complex(8),intent(INOUT) :: Hmat(ML,ML)
#endif
    integer :: i,i1,i2,i3,m,n,j,j1,j2,j3
    integer :: ML1,ML2,ML3
    real(8) :: d

    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    d = coef_lap0 + const_k2(k)

    do i=1,ML
       Hmat(i,i) = Hmat(i,i) + d
    end do

    if ( flag_nab ) then

       do i3=0,ML3-1
       do i2=0,ML2-1
       do i1=0,ML1-1

          i = 1 + i1 + i2*ML1 + i3*ML1*ML2

          do m=1,Md

             j1 = mod(i1+m+ML1,ML1)
             j  = 1 + j1 + i2*ML1 + i3*ML1*ML2
             Hmat(i,j) = Hmat(i,j) + zcoef_kin(1,m,k)

             j1 = mod(i1-m+ML1,ML1)
             j  = 1 + j1 + i2*ML1 + i3*ML1*ML2
             Hmat(i,j) = Hmat(i,j) + conjg( zcoef_kin(1,m,k) )

             j2 = mod(i2+m+ML2,ML2)
             j  = 1 + i1 + j2*ML1 + i3*ML1*ML2
             Hmat(i,j) = Hmat(i,j) + zcoef_kin(2,m,k)

             j2 = mod(i2-m+ML2,ML2)
             j  = 1 + i1 + j2*ML1 + i3*ML1*ML2
             Hmat(i,j) = Hmat(i,j) + conjg( zcoef_kin(2,m,k) )

             j3 = mod(i3+m+ML3,ML3)
             j  = 1 + i1 + i2*ML1 + j3*ML1*ML2
             Hmat(i,j) = Hmat(i,j) + zcoef_kin(3,m,k)

             j3 = mod(i3-m+ML3,ML3)
             j  = 1 + i1 + i2*ML1 + j3*ML1*ML2
             Hmat(i,j) = Hmat(i,j) + conjg( zcoef_kin(3,m,k) )

          end do ! m

       end do ! i1
       end do ! i2
       end do ! i3

    else

       do i3=0,ML3-1
       do i2=0,ML2-1
       do i1=0,ML1-1

          i = 1 + i1 + i2*ML1 + i3*ML1*ML2

          do m=1,Md

             j1 = mod(i1+m+ML1,ML1)
             j  = 1 + j1 + i2*ML1 + i3*ML1*ML2
             Hmat(i,j) = Hmat(i,j) + coef_lap(1,m)

             j1 = mod(i1-m+ML1,ML1)
             j  = 1 + j1 + i2*ML1 + i3*ML1*ML2
             Hmat(i,j) = Hmat(i,j) + coef_lap(1,m)

             j2 = mod(i2+m+ML2,ML2)
             j  = 1 + i1 + j2*ML1 + i3*ML1*ML2
             Hmat(i,j) = Hmat(i,j) + coef_lap(2,m)

             j2 = mod(i2-m+ML2,ML2)
             j  = 1 + i1 + j2*ML1 + i3*ML1*ML2
             Hmat(i,j) = Hmat(i,j) + coef_lap(2,m)

             j3 = mod(i3+m+ML3,ML3)
             j  = 1 + i1 + i2*ML1 + j3*ML1*ML2
             Hmat(i,j) = Hmat(i,j) + coef_lap(3,m)

             j3 = mod(i3-m+ML3,ML3)
             j  = 1 + i1 + i2*ML1 + j3*ML1*ML2
             Hmat(i,j) = Hmat(i,j) + coef_lap(3,m)

          end do ! m

       end do ! i1
       end do ! i2
       end do ! i3

    end if

    if ( flag_n12 .or. flag_n23 .or. flag_n31 ) then

       if ( flag_n12 ) then

          do i3=0,ML3-1
          do i2=0,ML2-1
          do i1=0,ML1-1

             i = 1 + i1 + i2*ML1 + i3*ML1*ML2

             do n=1,Md
             do m=1,Md

                d = -ggg(4)*coef_nab(1,m)*coef_nab(2,n)

                j1 = mod(i1+m+ML1,ML1)
                j2 = mod(i2+n+ML2,ML2)
                j  = 1 + j1 + j2*ML1 + i3*ML1*ML2
                Hmat(i,j) = Hmat(i,j) + d

                j1 = mod(i1-m+ML1,ML1)
                j2 = mod(i2+n+ML2,ML2)
                j  = 1 + j1 + j2*ML1 + i3*ML1*ML2
                Hmat(i,j) = Hmat(i,j) - d

                j1 = mod(i1+m+ML1,ML1)
                j2 = mod(i2-n+ML2,ML2)
                j  = 1 + j1 + j2*ML1 + i3*ML1*ML2
                Hmat(i,j) = Hmat(i,j) - d

                j1 = mod(i1-m+ML1,ML1)
                j2 = mod(i2-n+ML2,ML2)
                j  = 1 + j1 + j2*ML1 + i3*ML1*ML2
                Hmat(i,j) = Hmat(i,j) + d

             end do ! m
             end do ! n

          end do ! i1
          end do ! i2
          end do ! i3

       end if

       if ( flag_n23 ) then

          do i3=0,ML3-1
          do i2=0,ML2-1
          do i1=0,ML1-1

             i = 1 + i1 + i2*ML1 + i3*ML1*ML2

             do n=1,Md
             do m=1,Md

                d = -ggg(5)*coef_nab(2,m)*coef_nab(3,n)

                j2 = mod(i2+m+ML2,ML2)
                j3 = mod(i3+n+ML3,ML3)
                j  = 1 + i1 + j2*ML1 + j3*ML1*ML2
                Hmat(i,j) = Hmat(i,j) + d

                j2 = mod(i2-m+ML2,ML2)
                j3 = mod(i3+n+ML3,ML3)
                j  = 1 + i1 + j2*ML1 + j3*ML1*ML2
                Hmat(i,j) = Hmat(i,j) - d

                j2 = mod(i2+m+ML2,ML2)
                j3 = mod(i3-n+ML3,ML3)
                j  = 1 + i1 + j2*ML1 + j3*ML1*ML2
                Hmat(i,j) = Hmat(i,j) - d

                j2 = mod(i2-m+ML2,ML2)
                j3 = mod(i3-n+ML3,ML3)
                j  = 1 + i1 + j2*ML1 + j3*ML1*ML2
                Hmat(i,j) = Hmat(i,j) + d

             end do ! m
             end do ! n

          end do ! i1
          end do ! i2
          end do ! i3

       end if

       if ( flag_n31 ) then

          do i3=0,ML3-1
          do i2=0,ML2-1
          do i1=0,ML1-1

             i = 1 + i1 + i2*ML1 + i3*ML1*ML2

             do n=1,Md
             do m=1,Md

                d = -ggg(6)*coef_nab(1,m)*coef_nab(3,n)

                j1 = mod(i1+m+ML1,ML1)
                j3 = mod(i3+n+ML3,ML3)
                j  = 1 + j1 + i2*ML1 + j3*ML1*ML2
                Hmat(i,j) = Hmat(i,j) + d

                j1 = mod(i1-m+ML1,ML1)
                j3 = mod(i3+n+ML3,ML3)
                j  = 1 + j1 + i2*ML1 + j3*ML1*ML2
                Hmat(i,j) = Hmat(i,j) - d

                j1 = mod(i1+m+ML1,ML1)
                j3 = mod(i3-n+ML3,ML3)
                j  = 1 + j1 + i2*ML1 + j3*ML1*ML2
                Hmat(i,j) = Hmat(i,j) - d

                j1 = mod(i1-m+ML1,ML1)
                j3 = mod(i3-n+ML3,ML3)
                j  = 1 + j1 + i2*ML1 + j3*ML1*ML2
                Hmat(i,j) = Hmat(i,j) + d

             end do ! m
             end do ! n

          end do ! i1
          end do ! i2
          end do ! i3

       end if

    end if
 
  END SUBROUTINE construct_matrix_kinetic_sol


END MODULE kinetic_sol_module
