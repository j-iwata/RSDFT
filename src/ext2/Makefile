include ../make.inc

FFLAGS_FFTE = $(FFLAGS)
#FFLAGS_FFTE = -fast -xHOST -openmp

all : kernel.o factor.o pzfft3dv.o fft235.o zfft3d.o pzfft3d.o

clean :
	rm -r -f *.o *.lst
