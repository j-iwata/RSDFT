include ../make.sys

all:libtetrabz.a libtetrabz_mpi.a

clean:
	rm -f *.a *.o *.mod

libtetrabz.a:libtetrabz_common.o libtetrabz.o
	ar cr $@ libtetrabz_common.o libtetrabz.o

libtetrabz_mpi.a:libtetrabz_common.o libtetrabz_mpi.o
	ar cr $@ libtetrabz_common.o libtetrabz_mpi.o

libtetrabz.o:libtetrabz_common.o libtetrabz.f90
	$(F90) libtetrabz.f90 $(FFLAGS) -c

libtetrabz_mpi.o:libtetrabz_common.o libtetrabz_mpi.f90
	$(MPIF90) libtetrabz_mpi.f90 $(FFLAGS) -c

libtetrabz_common.o:libtetrabz_common.f90
	$(F90) libtetrabz_common.f90 $(FFLAGS) -c
