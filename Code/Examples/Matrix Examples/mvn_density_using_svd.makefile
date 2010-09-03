# mvn_density_using_svd.makefile
# To run: make -f mvn_density_using_svd.makefile

mvn_density_using_svd.out : mvn_density_using_svd.o
	g++ -static mvn_density_using_svd.o -lgsl -lgslcblas -lm \
	-o mvn_density_using_svd.out

mvn_density_using_svd.o : mvn_density_using_svd.c
	g++ -Wall -c mvn_density_using_svd.c

clean:
	rm mvn_density_using_svd.o mvn_density_using_svd.out