#------Makefile for ffind_peaks  ----------

# modules
SRC = ffind_peaks.f90
OBJ = ffind_peaks.o

INCLUDE    = ./

F2PY     = f2py

F2PY_F1   = --include-paths $(INCLUDE) --overwrite-signature -m 
F2PY_F2   = -c --fcompiler=gfortran 

PROGRAM = ffind_peaks

LIBS = 

#----------------------------------------------------------
all: $(OBJ)

%.o: %.f90
	$(F2PY) $(F2PY_F1) $(PROGRAM) -h sgn_$(PROGRAM).pyf $<
	$(F2PY) $(F2PY_F2) sgn_$(PROGRAM).pyf $(SRC) $(LIBS)

#----------------------------------------------------------
.PHONY : clean

clean:
	rm -f *.so *.pyf fluxpy.o
	rm -rf *.dSYM




