CH4_directInt: dcuhre.o pot_CH4_modify.o CH4_directInt.o 
	gfortran -o CH4_directInt dcuhre.o pot_CH4_modify.o CH4_directInt.o 
CH4_directInt.o: dcuhre.o pot_CH4_modify.o CH4_directInt.f03
	gfortran -c CH4_directInt.f03
pot_CH4_modify.o: pot_CH4_modify.f90
	gfortran -c pot_CH4_modify.f90
dcuhre.o: dcuhre.f
	gfortran -c dcuhre.f
clean: 
	rm CH4_directInt *.o
