# example.makefile
# http://www.cs.umd.edu/class/spring2002/cmsc214/Tutorial/makefile.html
# To run: make -f example.makefile

# A makefile maps a set of dependencies.  For instance,
# <example.out> depends on <example.o>.  When make
# determines that <example.out> is not up to date, given
# these dependencies, it runs the command g++ etc.
example.out : example.o
	g++ -static example.o -lgsl -lm -o example.out

example.o : example.c
	g++ -Wall -c example.c

# To clean up files
clean :
	rm example.o example.out