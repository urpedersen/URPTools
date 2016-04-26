PRG=dat2wav
all: src/$(PRG).cc
	mkdir -p ./bin
	g++ -O3 ./src/$(PRG).cc -o ./bin/$(PRG)

clean: 
	rm ./bin/$(PRG)
	rmdir bin
