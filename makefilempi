#-----------------------------------
#          MPI MAKEFILE
# ----------------------------------
#FORTRAN_COMPILER=mpif90
objects= global.o randomgen.o emax.o opt.o act.o

ifeq ($(ifortran),1)
	FORTRAN_COMPILER=ifort
	ifeq ($(d),1)
		switch= -llapack -O2 
	else
		switch= -llapack -debug -heap-arrays -CB 
	endif
else
	FORTRAN_COMPILER=mpif90
	ifeq ($(f),1)
		switch= -llapack -O2 
	else
		switch= -llapack -debug -heap-arrays -CB 
	endif
endif
exec: $(objects)
	$(FORTRAN_COMPILER) $(switch) -o exec $(objects)
global.mod: global.o global.f90
	$(FORTRAN_COMPILER) $(switch) -c global.f90
global.o: global.f90
	$(FORTRAN_COMPILER) $(switch)  -c global.f90
randomgen.mod: randomgen.o randomgen.f90
	$(FORTRAN_COMPILER) $(switch)  -c randomgen.f90
randomgen.o: randomgen.f90
	$(FORTRAN_COMPILER) $(switch) -c randomgen.f90
emax.mod: emax.o emax.f90
	$(FORTRAN_COMPILER) $(switch)  -c emax.f90
emax.o: emax.f90
	$(FORTRAN_COMPILER) $(switch) -c emax.f90
opt.mod: opt.o opt.f90
	$(FORTRAN_COMPILER) $(switch)  -c opt.f90
opt.o: opt.f90
	$(FORTRAN_COMPILER) $(switch) -c opt.f90
act.o: global.mod randomgen.mod emax.mod opt.mod act.f90
	$(FORTRAN_COMPILER) $(switch) -c act.f90
clean:
	rm *.mod *.o


