MODULE wf_sub_module

  use parallel_module
  use rgrid_module
  use rgrid_mol_module, only: LL
  use fft_module

  implicit none

  PRIVATE
  PUBLIC :: fft_initial_wf_sub, random_initial_wf_sub

!  complex(8),parameter :: half=(0.5d0,0.5d0)
  complex(8),parameter :: half=(0.0d0,0.0d0)

CONTAINS

  SUBROUTINE fft_initial_wf_sub &
       ( ML,MB,MK,MS,ML_0,ML_1,MB_0,MB_1,MK_0,MK_1,MS_0,MS_1,unk )
    implicit none
    integer,intent(IN) :: ML,MB,MK,MS,ML_0,ML_1,MB_0,MB_1,MK_0,MK_1,MS_0,MS_1
#ifdef _DRSDFT_
    real(8),intent(OUT) :: unk(ML_0:ML_1,MB,MK_0:MK_1,MS_0:MS_1)
#else
    complex(8),intent(OUT) :: unk(ML_0:ML_1,MB,MK_0:MK_1,MS_0:MS_1)
#endif
    integer :: s,k,n,i1,i2,i3,i,m0,m1,m12,ML1,ML2,ML3
    complex(8),allocatable :: z0(:,:,:),z1(:,:,:)
    real(8) :: c(2)

    call init_random_number

    m0  = Igrid(1,0)
    m1  = Igrid(2,1)-Igrid(1,1)+1
    m12 = (Igrid(2,1)-Igrid(1,1)+1)*(Igrid(2,2)-Igrid(1,2)+1)

    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    call init_fft

    allocate( z0(0:ML1-1,0:ML2-1,0:ML3-1) ) ; z0=(0.0d0,0.0d0)
    allocate( z1(0:ML1-1,0:ML2-1,0:ML3-1) ) ; z1=(0.0d0,0.0d0)

    do s=1,MS
       do k=1,MK
          do n=1,MB
             do i3=0,ML3-1
             do i2=0,ML2-1
             do i1=0,ML1-1
                call random_number(c)
                z0(i1,i2,i3) = dcmplx(c(1),c(2)) - half
             end do
             end do
             end do
             if ( MB_0 <= n .and. n <= MB_1 .and. &
                  MK_0 <= k .and. k <= MK_1 .and. &
                  MS_0 <= s .and. s <= MS_1 ) then
                call backward_fft( z0, z1 )
                i=m0-1
                do i3=Igrid(1,3),Igrid(2,3)
                do i2=Igrid(1,2),Igrid(2,2)
                do i1=Igrid(1,1),Igrid(2,1)
                   i=i+1
                   unk(i,n,k,s) = z0(i1,i2,i3)
                end do ! i1
                end do ! i2
                end do ! i3
             end if
          end do ! n
       end do ! k
    end do ! s

    deallocate( z1 )
    deallocate( z0 )

    call finalize_fft

  END SUBROUTINE fft_initial_wf_sub


  SUBROUTINE random_initial_wf_sub &
       ( ML,MB,MK,MS,ML_0,ML_1,MB_0,MB_1,MK_0,MK_1,MS_0,MS_1,unk,SYStype )
    implicit none
    integer,intent(IN) :: ML,MB,MK,MS,ML_0,ML_1,MB_0,MB_1,MK_0,MK_1,MS_0,MS_1
    integer,intent(IN) :: SYStype
#ifdef _DRSDFT_
    real(8),intent(INOUT) :: unk(ML_0:ML_1,MB,MK_0:MK_1,MS_0:MS_1)
    real(8),parameter :: zero=0.0d0
    real(8),allocatable :: w(:,:,:)
#else
    complex(8),intent(INOUT) :: unk(ML_0:ML_1,MB,MK_0:MK_1,MS_0:MS_1)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8),allocatable :: w(:,:,:)
#endif
    integer :: s,k,n,i1,i2,i3,i,m0,m1,m12,ML1,ML2,ML3
    real(8) :: c(2)

    call init_random_number

    ML1=Ngrid(1)
    ML2=Ngrid(2)
    ML3=Ngrid(3)

    allocate( w(0:ML1-1,0:ML2-1,0:ML3-1) ) ; w=zero

    do s=1,MS
       do k=1,MK
          do n=1,MB
             if ( SYStype == 0 ) then
                do i3=0,ML3-1
                do i2=0,ML2-1
                do i1=0,ML1-1
                   call random_number(c)
                   w(i1,i2,i3)=dcmplx(c(1),c(2)) - half
                end do
                end do
                end do
                if ( MB_0 <= n .and. n <= MB_1 .and. &
                     MK_0 <= k .and. k <= MK_1 .and. &
                     MS_0 <= s .and. s <= MS_1 ) then
                   i=Igrid(1,0)-1
                   do i3=Igrid(1,3),Igrid(2,3)
                   do i2=Igrid(1,2),Igrid(2,2)
                   do i1=Igrid(1,1),Igrid(2,1)
                      i=i+1
                      unk(i,n,k,s)=w(i1,i2,i3)
                   end do
                   end do
                   end do
                end if
             else
                do i3=0,ML3-1
                do i2=0,ML2-1
                do i1=0,ML1-1
                   call random_number(c)
                   w(i1,i2,i3)=dcmplx(c(1),c(2)) - half
                end do
                end do
                end do
                if ( MB_0 <= n .and. n <= MB_1 .and. &
                     MK_0 <= k .and. k <= MK_1 .and. &
                     MS_0 <= s .and. s <= MS_1 ) then
                   do i=Igrid(1,0),Igrid(2,0)
                      i1=LL(1,i)+(ML1-1)/2
                      i2=LL(2,i)+(ML2-1)/2
                      i3=LL(3,i)+(ML3-1)/2
                      unk(i,n,k,s)=w(i1,i2,i3)
                   end do
                end if
             end if
          end do ! n
       end do ! k
    end do ! s

    deallocate( w )

  END SUBROUTINE random_initial_wf_sub

  SUBROUTINE init_random_number
    implicit none
    integer :: n
    integer,allocatable :: ir(:)
    call random_seed( size=n )
    allocate( ir(n) ) ; ir=0
    call random_seed( put=ir )
    deallocate( ir )
  END SUBROUTINE init_random_number

END MODULE wf_sub_module
