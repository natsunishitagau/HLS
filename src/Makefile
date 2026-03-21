all: CouplingGraph.o QASMReader.o HLS clean

CouplingGraph.o: CouplingGraph.h CouplingGraph.cpp
	g++ -O3 -static -std=c++17 -c CouplingGraph.cpp -o CouplingGraph.o

QASMReader.o: QASMReader.h QASMReader.cpp
	g++ -O3 -static -std=c++17 -c QASMReader.cpp -o QASMReader.o

HLS: CouplingGraph.o QASMReader.o gate.h WeightedGate.h main.cpp
	g++ -O3 -static -std=c++17 CouplingGraph.o QASMReader.o gate.h WeightedGate.h main.cpp -o HLS

clean:
	rm -f *~
	rm -f *.o

cleanup: clean
	rm -f HLS