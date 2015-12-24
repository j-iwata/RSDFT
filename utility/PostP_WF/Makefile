include make.inc

MODS0 = parallel_module.o watch_module.o atom_module.o
MODS1 = wf_data_module.o
OBJS0 = wf_data_main.o convert_capital.o

all : $(MODS0) $(MODS1) $(OBJS0)
	$(FC) $(LFLAGS) $(MODS0) $(MODS1) $(OBJS0)



clean :
	rm -r -f *.o *~ a.out *.mod *.lst
