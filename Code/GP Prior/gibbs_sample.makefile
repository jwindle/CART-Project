# gibbs_sample.makefile
# To run: make -f gibbs_sample.makefile

gibbs_sample.out : gibbs_sample.o
	g++ -static gibbs_sample.o -lgsl -lgslcblas -lm \
	-o gibbs_sample.out

gibbs_sample.o : gibbs_sample.c
	g++ -Wall -c gibbs_sample.c

clean:
	rm gibbs_sample.o gibbs_sample.out