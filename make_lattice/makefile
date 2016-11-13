CC = g++
CFLAGS  = -O3 -Wall
INCLUDES = 
LFLAGS = 
LIBS = -lboost_iostreams

SRCDIR = src
SRCS = $(SRCDIR)/make_lattice.cc $(SRCDIR)/lattice.cc $(SRCDIR)/split.cc 
OBJS = $(SRCS:.cc=.o)
MAINDIR = bin
MAIN = $(MAINDIR)/make_lattice

MKDIR = mkdir -p

.PHONY: clean subdirs $(SUBDIRS)

all:    $(MAIN)
	@echo   The program is compiled. The executable is located in the ./bin directory.

$(MAIN): $(OBJS)
	$(MKDIR) $(MAINDIR)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

.cc.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) ./src/*.o src./*~ $(MAIN)

