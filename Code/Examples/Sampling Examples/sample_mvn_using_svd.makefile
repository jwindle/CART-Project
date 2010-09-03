# makefile.sample_mvn_using_svd
# To run: make -f makefile.sample_mvn_using_svd

sample_mvn_using_svd.out : sample_mvn_using_svd.o
	g++ -static sample_mvn_using_svd.o -lgsl -lgslcblas -lm \
	-o sample_mvn_using_svd.out

sample_mvn_using_svd.o : sample_mvn_using_svd.c
	g++ -Wall -c sample_mvn_using_svd.c

clean:
	rm sample_mvn_using_svd.o sample_mvn_using_svd.out