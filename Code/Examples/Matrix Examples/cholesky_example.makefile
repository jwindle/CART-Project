# cholesky_example.makefile
# To run: make -f cholesky_example.makefile

cholesky_example.out : cholesky_example.o
	g++ -static cholesky_example.o -lgsl -lgslcblas -lm \
	-o cholesky_example.out

cholesky_example.o : cholesky_example.c
	g++ -Wall -c cholesky_example.c

clean:
	rm cholesky_example.o cholesky_example.out