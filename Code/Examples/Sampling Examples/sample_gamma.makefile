# sample_gamma.makefile
# To run: make -f sample_gamma.makefile

sample_gamma.out : sample_gamma.o
	g++ -static sample_gamma.o -lgsl -lgslcblas -lm \
	-o sample_gamma.out

sample_gamma.o : sample_gamma.c
	g++ -Wall -c sample_gamma.c

clean:
	rm sample_gamma.o sample_gamma.out