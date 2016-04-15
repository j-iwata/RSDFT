include ../make.sys

dos2gp: dos2gp.f90 ../src/libtetrabz_mpi.a
	$(MPIF90) $< -L $(TOPDIR)/src/ -I $(TOPDIR)/src/ -ltetrabz_mpi $(FFLAGS) -o $@ 

clean:
	rm -f *.x *.a *.o *.mod dos2gp

