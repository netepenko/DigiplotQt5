#------Windows (mingw32) Makefile for ffind_peaks  ----------
# use: mingw32-make -f ffind_peaks_win.mak to compile the code
# 
# for a fresh start:
#       mingw32-make -f ffind_peaks_win.mak clean
# 
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
	copy .\ffind_peaks\.libs\* . 

%.o: %.f90
	$(F2PY) $(F2PY_F1) $(PROGRAM) -h sgn_$(PROGRAM).pyf $<
	$(F2PY) $(F2PY_F2) sgn_$(PROGRAM).pyf $(SRC) $(LIBS)

#----------------------------------------------------------
.PHONY : clean

clean:
	del /Q *.pyd *.pyf *.dll
	del /Q ffind_peaks\.libs\*



