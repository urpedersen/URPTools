PRG=xyz2sofq
CC=g++
CFLAGS=-O3

all: 
	mkdir -p ./bin
	$(CC) $(CFLAGS) ./src/helpers_io.h ./src/helpers_io.cc ./src/$(PRG).cc -o ./bin/$(PRG)

clean: 
	rm ./bin/$(PRG)
	rmdir bin
