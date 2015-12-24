include ../make.inc

FFLAGS_FFTE = $(FFLAGS)
#FFLAGS_FFTE = -fast -xHOST -openmp

FFTE_DIR = ffte_dummy
#FFTE_DIR = ffte-6.0

all : kernel.o factor.o pzfft3dv.o fft235.o zfft3d.o pzfft3d.o

kernel.o:
	cp $(FFTE_DIR)/param.h . ; cp $(FFTE_DIR)/kernel.f . ; $(FC) $(FFLAGS_FFTE) -c kernel.f
factor.o:
	cp $(FFTE_DIR)/factor.f . ; $(FC) $(FFLAGS_FFTE) -c factor.f
pzfft3dv.o:
	cp $(FFTE_DIR)/mpi/pzfft3dv.f . ; $(FC) $(FFLAGS_FFTE) -c pzfft3dv.f
fft235.o:
	cp $(FFTE_DIR)/fft235.f . ; $(FC) $(FFLAGS_FFTE) -c fft235.f
zfft3d.o:
	cp $(FFTE_DIR)/zfft3d.f . ; $(FC) $(FFLAGS_FFTE) -c zfft3d.f
pzfft3d.o:
	cp $(FFTE_DIR)/mpi/pzfft3d.f . ; $(FC) $(FFLAGS_FFTE) -c pzfft3d.f

clean :
	rm -r -f *.o *.h *.lst *.f
