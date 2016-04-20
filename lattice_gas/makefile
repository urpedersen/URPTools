all: src/lattice_gas.cc
	mkdir -p ./bin
	g++ -O3 ./src/lattice.h ./src/lattice.cc ./src/lattice_gas.cc -lboost_iostreams -o ./bin/lattice_gas

clean: 
	rm ./bin/make_lattice
	rmdir bin
