MODULE basic_type_methods

  USE basic_type_factory

  implicit none

  PRIVATE
  PUBLIC :: allocateGS
  PUBLIC :: allocateGSArray, allocateGSArray_v2
  PUBLIC :: allocateGBKS
  PUBLIC :: allocaterGBKS
  PUBLIC :: allocatecGBKS

#ifdef REAL_VER
  double precision,parameter :: zero = 0.d0
#else
  complex(kind(0d0)),parameter :: zero = (0.d0,0.d0)
#endif
  complex(kind(0d0)),parameter :: z0 = (0.d0,0.d0)

CONTAINS

  SUBROUTINE allocateGSArray( gs )
    implicit none
    type( GSArray ) :: gs
    integer :: m1,m2,n1,n2
    m1=gs%g_range%head
    m2=gs%g_range%tail
    n1=gs%s_range%head
    n2=gs%s_range%tail
    allocate( gs%val(m1:m2,n1:n2) ) ; gs%val=0.0d0
  END SUBROUTINE allocateGSArray

  SUBROUTINE allocateGSArray_v2( gs )
    implicit none
    type( GSArray_v2 ) :: gs
    integer :: m1,m2,n1,n2
    m1=gs%g_range%alloc%head
    m2=gs%g_range%alloc%tail
    n1=gs%s_range%alloc%head
    n2=gs%s_range%alloc%tail
    allocate( gs%val(m1:m2,n1:n2) ) ; gs%val=0.0d0
  END SUBROUTINE allocateGSArray_v2

  SUBROUTINE allocateGS( gs )
    implicit none
    type( GSArray    ) ::  gs
    allocate( gs%val(gs%g_prange%head:gs%g_prange%tail, gs%s_srange%head:gs%s_srange%tail) )
    gs%val = 0.d0
  END SUBROUTINE allocateGS

  SUBROUTINE allocateGBKS( gbks )
    implicit none
    type( GBKSArray  ) ::  gbks
    allocate( gbks%val(gbks%g_prange%head:gbks%g_prange%tail, gbks%b_prange%head:gbks%b_prange%tail, gbks%k_prange%head:gbks%k_prange%tail, gbks%s_prange%head:gbks%s_prange%tail) )
    gbks%val = zero
  END SUBROUTINE allocateGBKS

  SUBROUTINE allocaterGBKS( gbks )
    implicit none
    type( rGBKSArray ) :: gbks
    allocate( gbks%val(gbks%g_prange%head:gbks%g_prange%tail, gbks%b_prange%head:gbks%b_prange%tail, gbks%k_prange%head:gbks%k_prange%tail, gbks%s_prange%head:gbks%s_prange%tail) )
    gbks%val = 0.d0
  END SUBROUTINE allocaterGBKS

  SUBROUTINE allocatecGBKS( gbks )
    implicit none
    type( cGBKSArray ) :: gbks
    allocate( gbks%val(gbks%g_prange%head:gbks%g_prange%tail, gbks%b_prange%head:gbks%b_prange%tail, gbks%k_prange%head:gbks%k_prange%tail, gbks%s_prange%head:gbks%s_prange%tail) )
    gbks%val = z0
  END SUBROUTINE allocatecGBKS

END MODULE basic_type_methods
