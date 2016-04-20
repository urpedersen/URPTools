all: src/lattice_gas.cc
	mkdir -p ./bin
	g++ -O3 ./src/lattice.h ./src/lattice.cc ./src/lattice_gas.cc -lboost_program_options -o ./bin/lattice_gas

clean: 
	rm ./bin/lattice_gas
	rmdir bin
