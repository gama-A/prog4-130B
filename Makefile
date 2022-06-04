# Makefile
CXX_FLAG = --std=c++11 -g

prog4: prog4.o
	g++ $(CXX_FLAG) -o prog4 prog4.o -fopenmp

prog4.o: prog4.cpp
	g++ $(CXX_FLAG) prog4.cpp -c

clean:
	rm -f prog4 *.o