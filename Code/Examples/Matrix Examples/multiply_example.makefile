# multiply_example.makefile
# To run: make -f multiply_example.makefile

multiply_example.out : multiply_example.o
	g++ -static multiply_example.o -lgsl -lgslcblas -lm \
	-o multiply_example.out

multiply_example.o : multiply_example.c
	g++ -Wall -c multiply_example.c

clean:
	rm multiply_example.o multiply_example.out