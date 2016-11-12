CC = g++
CFLAGS  = -O3 -Wall
INCLUDES = 
LFLAGS = 
LIBS = -lboost_iostreams

SRCS = ./src/make_lattice.cc ./src/lattice.cc ./src/split.cc 
OBJS = $(SRCS:.cc=.o)
MAIN = ./bin/make_lattice

.PHONY: clean

all:    $(MAIN)
	@echo   The program is compiled. The executable is located in the ./bin directory.

$(MAIN): $(OBJS) 
	mkdir -p bin
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

.cc.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) ./src/*.o src./*~ $(MAIN)

