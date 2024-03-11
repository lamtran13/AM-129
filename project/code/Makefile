FC = gfortran
FFLAGS = -Wall -Wextra -Wno-surprising -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace -J ./obj
LFLAGS = -llapack

ODIR=obj
_OBJ = utility.o setup_module.o leapfrog_module.o output_module.o fput.o   
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

fput.ex: $(OBJ)
	$(FC) -o $@ $^ $(LFLAGS) $(FFLAGS)

$(ODIR)/%.o: %.f90
	@mkdir -p $(ODIR)
	$(FC) -c -o $@ $< $(FFLAGS)

.PHONY: clean

clean:
	rm -f fput.ex *.o *.mod *~
