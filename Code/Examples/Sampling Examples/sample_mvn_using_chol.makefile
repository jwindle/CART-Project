# makefile.sample_mvn_using_chol
# To run: make -f makefile.sample_mvn_using_chol

sample_mvn_using_chol.out : sample_mvn_using_chol.o
	g++ -static sample_mvn_using_chol.o -lgsl -lgslcblas -lm \
	-o sample_mvn_using_chol.out

sample_mvn_using_chol.o : sample_mvn_using_chol.c
	g++ -Wall -c sample_mvn_using_chol.c

clean:
	rm sample_mvn_using_chol.o sample_mvn_using_chol.out