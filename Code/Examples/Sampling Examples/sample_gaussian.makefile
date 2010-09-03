# makefile.sample_gaussian
# To run: make -f makefile.sample_gaussian

sample_gaussian.out : sample_gaussian.o
	g++ -static sample_gaussian.o -lgsl -lgslcblas -lm \
	-o sample_gaussian.out

sample_gaussian.o : sample_gaussian.c
	g++ -Wall -c sample_gaussian.c

clean:
	rm sample_gaussian.o sample_gaussian.out