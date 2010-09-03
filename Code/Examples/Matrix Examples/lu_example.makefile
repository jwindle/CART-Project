# lu_example.makefile
# To run: make -f lu_example.makefile

lu_example.out : lu_example.o
	g++ -static lu_example.o -lgsl -lgslcblas -lm \
	-o lu_example.out

lu_example.o : lu_example.c
	g++ -Wall -c lu_example.c

clean:
	rm lu_example.o lu_example.out