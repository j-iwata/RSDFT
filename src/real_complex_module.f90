MODULE real_complex_module

    implicit none
    include 'mpif.h'
    
    
#ifdef _DRSDFT_
    integer,parameter :: TYPE_MAIN=MPI_REAL8
    character(1),parameter :: TRANSA='T'
    character(1),parameter :: TRANSB='N'
    real(8),parameter :: zero=0.d0
    real(8),parameter :: one=1.d0
#else
    integer,parameter :: TYPE_MAIN=MPI_COMPLEX16
    character(1),parameter :: TRANSA='C'
    character(1),parameter :: TRANSB='N'
    complex(8),parameter :: zero=(0.d0,0.d0)
    complex(8),parameter :: one=(1.d0,0.d0)
#endif

CONTAINS


  !---------------------------------------------------------------------------------------

  SUBROUTINE RCProduct( INA,INB,OUTC )
    implicit none
#ifdef _DRSDFT_
    real(8),intent(IN) :: INA,INB
    real(8),intent(OUT) :: OUTC

    OUTC = INA*INB
#else
    complex(8),intent(IN) :: INA,INB
    complex(8),intent(OUT) :: OUTC

    OUTC = conjg(INA)*INB
#endif
  END SUBROUTINE RCProduct

END MODULE real_complex_module
