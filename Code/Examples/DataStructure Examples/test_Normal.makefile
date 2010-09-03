# test_Normal.makefile
# To run: make -f test_Normal.makefile

test_Normal.out : test_Normal.o Normal.o 
		g++ -static Normal.o test_Normal.o -lgsl -lgslcblas -lm \
		-o test_Normal.out

test_Normal.o : test_Normal.c Normal.h
	      g++ -Wall -c test_Normal.c
	   
Normal.o : Normal.cpp Normal.h
	 g++ -Wall -c Normal.cpp

clean:
	rm Normal.o test_Normal.o