main:main.o BVP.o
	g++ -o main main.o BVP.o -llapacke -ljsoncpp
main.o:main.cpp
	g++ -c main.cpp
BVP.o:BVP.h BVP.cpp 
	g++ -c BVP.cpp 
clean:
	rm *.o main
