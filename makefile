#-----------------------------------
#          MPI MAKEFILE
# ----------------------------------
#FORTRAN_COMPILER=mpif90
objects= global.o randomgen.o emax.o optu.o act.o
FORTRAN_COMPILER=mpif90

ifeq ($(f),1)
	switch= -O2 -DOPENBLAS_ARCH_X86_64=1 -DOPENBLAS___64BIT__=1 -DOPENBLAS_ARCH_X86_64=1 -DOPENBLAS___64BIT__=1 -l/usr/include/atlas -l/usr/include/openblas -lptlapack -lopenblas_openmp
else
	switch=  -O2 -DOPENBLAS_ARCH_X86_64=1 -DOPENBLAS___64BIT__=1 -DOPENBLAS_ARCH_X86_64=1 -DOPENBLAS___64BIT__=1 -l/usr/include/atlas -l/usr/include/openblas -lptlapack -lopenblas_openmp -debug -heap-arrays -CB 
endif

exec: $(objects)
	$(FORTRAN_COMPILER) $(switch) -o exec $(objects) -I/home/utkusu/nlopt/include -L/home/utkusu/nlopt/lib -lnlopt -lm 
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
optu.mod: optu.o optu.f90
	$(FORTRAN_COMPILER) $(switch)  -c optu.f90
optu.o: optu.f90
	$(FORTRAN_COMPILER) $(switch) -c optu.f90
act.o: global.mod randomgen.mod emax.mod optu.mod act.f90
	$(FORTRAN_COMPILER) $(switch) act.f90 -I/home/utkusu/nlopt/include -L/home/utkusu/nlopt/lib -lnlopt -lm -c 
clean:
	rm *.mod *.o


