include ../make.sys

CC=icc -lifcore

all:dos.x lindhard.x

clean:
	rm -f *.x *.a *.o *.mod

dos.x:dos.f90 ../src/libtetrabz.a
	$(F90) $< -L $(TOPDIR)/src/ -I $(TOPDIR)/src/ -ltetrabz $(FFLAGS) -o $@ 
	rm -f *.o *.mod

lindhard.x:lindhard.f90 ../src/libtetrabz.a
	$(F90) $< -L $(TOPDIR)/src/ -I $(TOPDIR)/src/ -ltetrabz $(FFLAGS) -o $@ 
	rm -f *.o *.mod

dos_c.x:dos.c ../src/libtetrabz.a
	$(CC) $< -L $(TOPDIR)/src/ -I $(TOPDIR)/src/ -ltetrabz $(FFLAGS) -o $@ 
	rm -f *.o

