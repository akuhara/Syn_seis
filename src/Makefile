#
BINDIR  = ../bin
FC      = gfortran
FCFLAGS = -g -Wall -fbounds-check -O -Wuninitialized \
          -ffpe-trap=invalid,zero,overflow -fbacktrace
SAC     = -I/usr/local/sac/include -L/usr/local/sac/lib -lsacio
FFTW    = -I/usr/local/include -lfftw3
LIB_INC = $(SAC) $(FFTW)
OBJS    =  syn_seis.o read_input.o param.o fwd_seis.o output_sac.o
TARGET  = $(BINDIR)/syn_seis


$(TARGET): $(OBJS)
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(FC) $^ $(LIB_INC) -o $@

syn_seis.o: param.o
read_input.o: param.o
fwd_seis.o: param.o
output_sac.o: param.o
%.o: %.f90
	$(FC) $(FCFLAGS) -c $< $(LIB_INC)
%.mod: %.f90 %.o
	@:

.PHONY: clean
clean: 
	rm *.o *.mod $(BINDIR)/syn_seis
