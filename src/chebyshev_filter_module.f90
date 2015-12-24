MODULE ChebyshevFilter_module

  use rgrid_module
  use parallel_module
  use hamiltonian_module
  use wf_module

  implicit none

  PRIVATE
  PUBLIC :: ChebyshevFilter, Init_ChebyshevFilter

  integer :: ML_0, ML_1
  integer :: kmax = 5
  integer :: m_poly = 10

#ifdef _DRSDFT_
  real(8),allocatable :: Y(:,:), Z(:,:)
  real(8),parameter :: zero=0.0d0
#else
  complex(8),allocatable :: Y(:,:), Z(:,:)
  complex(8),parameter :: zero=(0.0d0,0.0d0)
#endif

CONTAINS


  SUBROUTINE Init_ChebyshevFilter( rank, unit )
    implicit none
    integer,intent(IN) :: rank,unit
    call Read_ChebyshevFilter( rank, unit )
  END SUBROUTINE Init_ChebyshevFilter
    
  SUBROUTINE Read_ChebyshevFilter( rank, unit )
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i,p(2),ierr
    character(6) :: cbuf,ckey
    p(:)=0
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey == "CHEBFI" ) then
             backspace(unit)
             read(unit,*) cbuf,p(1),p(2)
          end if
       end do
999    continue
       if ( p(1) > 0 ) m_poly=p(1) 
       if ( p(2) > 0 ) kmax=p(2) 
       write(*,*) "m_poly=",m_poly
       write(*,*) "kmax  =",kmax
    end if
    call MPI_BCAST(m_poly,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(kmax  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  END SUBROUTINE Read_ChebyshevFilter


  SUBROUTINE ChebyshevFilter( k, s, MB_0, MB_1 )
    implicit none
    integer,intent(IN) :: k, s
    integer,intent(IN) :: MB_0, MB_1

    integer :: ns,ne,nn,m
    real(8) :: e_lower, e_upper, e_lowest
    real(8) :: e,c,sigma,sigma1,sigma2,c1,c2,cc

    if ( disp_switch_parallel ) then
       write(*,'(a20," ChebyshevFilter(m=",i2,")")') repeat("-",20),m_poly
    end if

    ML_0 = Igrid(1,0)
    ML_1 = Igrid(2,0)

! ---

    call get_upper_bound( k, s, e_upper )

! ---

    e_lower = maxval( esp(:,k,s) )
    if ( e_lower > e_upper ) e_lower = e_upper - abs(e_upper)*0.5d0

    e_lowest = minval( esp(:,k,s) )
    if ( e_lowest > e_lower ) e_lowest = e_lower - abs(e_lower)*0.5d0

!    if ( disp_switch_parallel ) then
       write(*,'(1x,"e(lowest,lower,upper)=",3f15.10)') e_lowest,e_lower,e_upper
!    end if

! ---

    allocate( Y(ML_0:ML_1,MB_d) ) ; Y=zero
    allocate( Z(ML_0:ML_1,MB_d) ) ; Y=zero

    e = ( e_upper - e_lower )*0.5d0
    c = ( e_upper + e_lower )*0.5d0
!   sigma = e/( e_lower - c )
    sigma = e/( e_lower - e_lowest )
    sigma1= sigma

    c1 = sigma1/e

    do ns=MB_0,MB_1,MB_d
       ne=min( ns+MB_d-1, MB_1 )
       nn=ne-ns+1

       call hamiltonian(k,s,unk(ML_0,ns,k,s),Y,ML_0,ML_1,ns,ne)
       Y(ML_0:ML_1,1:nn)=c1*( Y(ML_0:ML_1,1:nn)-c*unk(ML_0:ML_1,ns:ne,k,s) )

       do m=2,m_poly

          sigma2 = 1.0d0/(2.0d0/sigma1-sigma)

          c2 = 2.0d0*sigma2/e
          cc = sigma*sigma2
          call hamiltonian(k,s,Y,Z,ML_0,ML_1,ns,ne)
          Z(ML_0:ML_1,1:nn) = c2*( Z(ML_0:ML_1,1:nn) - c*Y(ML_0:ML_1,1:nn) ) &
               - cc*unk(ML_0:ML_1,ns:ne,k,s)

          if ( m == m_poly ) then
             unk(ML_0:ML_1,ns:ne,k,s) = Z(ML_0:ML_1,1:nn)
             exit
          else
             unk(ML_0:ML_1,ns:ne,k,s) = Y(ML_0:ML_1,1:nn)
             Y(ML_0:ML_1,1:nn) = Z(ML_0:ML_1,1:nn)
          end if

          sigma = sigma2

       end do ! m

    end do ! ns

    deallocate( Z )
    deallocate( Y )

! ---

    if ( disp_switch_parallel ) then
       write(*,'(a20," ChebyshevFilter(end)")') repeat("-",20)
    end if

  END SUBROUTINE ChebyshevFilter


  SUBROUTINE get_upper_bound( k, s, eu_guess )
    implicit none
    integer,intent(IN) :: k, s
    real(8),intent(OUT) :: eu_guess
    integer :: i,j,ierr
    real(8) :: c(2),d(3),alpha,beta
    real(8),allocatable :: T(:,:)
    real(8),allocatable :: Td(:),Te(:),work(:),w(:)
    real(8) :: vl,vu
    integer,allocatable :: iwork(:),iblock(:),isplit(:)
    integer :: m, nsplit, info, il, iu
#ifdef _DRSDFT_
    real(8),parameter :: zero=0.0d0
    real(8),allocatable :: v(:),v0(:),f(:)
#else
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8),allocatable :: v(:),v0(:),f(:)
#endif

    ML_0 = Igrid(1,0)
    ML_1 = Igrid(2,0)

    allocate( T(kmax,kmax)  ) ; T=0.0d0
    allocate( v(ML_0:ML_1)  ) ; v=zero
    allocate( v0(ML_0:ML_1) ) ; v0=zero
    allocate( f(ML_0:ML_1)  ) ; f=zero

    do i=1,ML_0-1
       call random_number(c)
    end do
    do i=ML_0,ML_1
       call random_number(c)
       v(i) = dcmplx( c(1), c(2) )
    end do
    c(1) = sum( abs(v)**2 )*dV
    call MPI_ALLREDUCE(c(1),c(2),1,MPI_REAL8,MPI_SUM,comm_grid,ierr)
    c(2) = 1.0d0/sqrt(c(2))
    v(:) = c(2)*v(:)

    call hamiltonian( k, s, v, f, ML_0, ML_1, 1, 1 )

#ifdef _DRSDFT_
    c(1) = sum( f*v )*dV
#else
    c(1) = sum( conjg(f)*v )*dV
#endif
    call MPI_ALLREDUCE(c,alpha,1,MPI_REAL8,MPI_SUM,comm_grid,ierr)

    f(:) = f(:) - alpha*v(:)

    T(1,1) = alpha

    do j=2,kmax

       c(1) = sum( abs(f)**2 )*dV
       call MPI_ALLREDUCE(c,beta,1,MPI_REAL8,MPI_SUM,comm_grid,ierr)
       beta = sqrt(beta)

       v0(:) = v(:)
       c(1)  = 1.0d0/beta
       v(:)  = c(1)*f(:)

       call hamiltonian( k, s, v, f, ML_0, ML_1, 1, 1 )
       f(:) = f(:) - beta*v0(:)

#ifdef _DRSDFT_
       c(1) = sum( f*v )*dV
#else
       c(1) = sum( conjg(f)*v )*dV
#endif
       call MPI_ALLREDUCE(c,alpha,1,MPI_REAL8,MPI_SUM,comm_grid,ierr)

       f(:) = f(:) - alpha*v(:)

       T(j,j-1) = beta
       T(j-1,j) = beta
       T(j,j)   = alpha

       c(1) = sum( T(1:j,1:j)**2 )
       c(2) = sum( abs(f)**2 )*dV
       call MPI_ALLREDUCE(c,d,2,MPI_REAL8,MPI_SUM,comm_grid,ierr)
       d(1:2)=sqrt( d(1:2) )

       allocate( Td(j), Te(j-1), work(4*j), iwork(3*j) )
       allocate( w(j), iblock(j), isplit(j) )
       do i=1,j
          Td(i) = T(i,i)
          if ( i < j ) Te(i)=T(i+1,i)
       end do
       call dstebz( 'A', 'E', j, vl, vu, il, iu, 0.0d0, Td, Te, m, nsplit &
            ,w, iblock, isplit, work, iwork, info )
       d(1)=maxval(w)
       d(3)=minval(w)
       deallocate( isplit, iblock, w )
       deallocate( iwork, work, Te, Td )

       eu_guess = d(1) + d(2) 

!       if ( disp_switch_parallel ) then
!          write(*,'(1x,i3,4f20.15)') j,eu_guess,d(1:3)
!       end if

    end do ! j

    deallocate( f  )
    deallocate( v  )
    deallocate( v0 )
    deallocate( T  )

  END SUBROUTINE get_upper_bound


END MODULE ChebyshevFilter_module
