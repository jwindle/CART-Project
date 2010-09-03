# svd_example.makefile
# To run: make -f svd_example.makefile

svd_example.out : svd_example.o
	g++ -static svd_example.o -lgsl -lgslcblas -lm \
	-o svd_example.out

svd_example.o : svd_example.c
	g++ -Wall -c svd_example.c

clean:
	rm svd_example.o svd_example.out