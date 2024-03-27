main:main.o AdamsBashforth.o ClassicalRK.o FehlbergRK.o DormandPrinceRK.o AdamsMoulton.o BDFs.o ESDIRK.o GaussLegendre.o
	g++ -o main main.o AdamsBashforth.o AdamsMoulton.o ClassicalRK.o FehlbergRK.o DormandPrinceRK.o BDFs.o ESDIRK.o GaussLegendre.o -ljsoncpp
main.o:main.cpp 
	g++ -c main.cpp -std=c++11
AdamsBashforth.o:AdamsBashforth.cpp
	g++ -c AdamsBashforth.cpp -std=c++11
AdamsMoulton.o:AdamsMoulton.cpp
	g++ -c AdamsMoulton.cpp -std=c++11
BDFs.o:BDFs.cpp
	g++ -c BDFs.cpp -std=c++11
ESDIRK.o:ESDIRK.cpp
	g++ -c ESDIRK.cpp -std=c++11
GaussLegendre.o:GaussLegendre.cpp
	g++ -c GaussLegendre.cpp -std=c++11
ClassicalRK.o:ClassicalRK.cpp
	g++ -c ClassicalRK.cpp -std=c++11
FehlbergRK.o:FehlbergRK.cpp
	g++ -c FehlbergRK.cpp -std=c++11
DormandPrinceRK.o:DormandPrinceRK.cpp
	g++ -c DormandPrinceRK.cpp -std=c++11
clean:
	rm *.o main 
