main:main.o Multigrid.o
	g++ -o main main.o Multigrid.o -llapacke -ljsoncpp
main.o:main.cpp
	g++ -c main.cpp -std=c++11
Multigrid.o:Multigrid.cpp
	g++ -c Multigrid.cpp -std=c++11
clean:
	rm *.o main
