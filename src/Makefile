include make.inc

########################################################################
########################################################################

include makefile.common

DIR0 = ext0
EXTOBJ0 = $(DIR0)/dotp.o \
          $(DIR0)/convert_capital.o \
          $(DIR0)/write_info.o \
          $(DIR0)/stop_program.o \

DIR1 = ext1
EXTOBJ1 = $(DIR1)/ext_sub_minpac.o

DIR2 = ext2
EXTOBJ2 = $(DIR2)/pzfft3dv.o \
          $(DIR2)/factor.o \
          $(DIR2)/fft235.o \
          $(DIR2)/kernel.o \
          $(DIR2)/zfft3d.o \
          $(DIR2)/pzfft3d.o \

########################################################################
########################################################################

.PHONY: all clean re test runtest cleartest

all :
	@$(MAKE) lda0
	cd $(DIR0) ; $(MAKE)
	cd $(DIR1) ; $(MAKE)
	cd $(DIR2) ; $(MAKE)
	@$(MAKE) main.o
	$(FC) $(LFLAGS) $(EXTOBJ0) $(EXTOBJ1) $(EXTOBJ2) $(MODS1) main.o $(FFTW_L) $(LAPACK_L) -o rsdft.x

lda0 : $(MODS1)

clean :
	rm -f *.o *.mod mpif.h *.lst *.optlog *.i90
	cd $(DIR0) ; $(MAKE) clean
	cd $(DIR1) ; $(MAKE) clean
	cd $(DIR2) ; $(MAKE) clean

