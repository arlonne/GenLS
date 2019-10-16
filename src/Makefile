.SUFFIXES: .o .f90

FC = ifort
OFLAG = -FR
FORT = 

#SOURCE = modpermu_part.o rmreppermu.o splitblock.o modaddmove.o modinterface.o polymer_main-test.o
SOURCE = modpermu_part.o modpermu.o rmreppermu.o modrdpermu.o splitblock.o modaddmove.o modconstruct.o modinterface.o gendefect.o genstack.o GenLS-main.o

default: GenLS
GenLS: $(SOURCE)
	/bin/rm -f ../GenLS.ne
	$(FC) $(FORT) $(OFLAG) $(SOURCE) -o ../GenLS.ne

clean:
	/bin/rm -f *.mod *.o 

.f90.o:
	$(FC) $(FORT) $(OFLAG) -c $<
