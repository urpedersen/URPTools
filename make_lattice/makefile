all: src/make_lattice.cc
	mkdir -p bin
	g++ -O3 ./src/split.h ./src/split.cc ./src/lattice.h ./src/lattice.cc ./src/make_lattice.cc -lboost_iostreams -o ./bin/make_lattice

clean: 
	rm bin/make_lattice
	rmdir bin
