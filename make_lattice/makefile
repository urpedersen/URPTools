CXX = g++
CXXFLAGS  = -O3 -Wall
LDFLAGS =
CPPFLAGS = 
LIBS = -lboost_iostreams

SRCDIR = src
SRCS = $(SRCDIR)/make_lattice.cc $(SRCDIR)/lattice.cc $(SRCDIR)/split.cc 
OBJS = $(SRCS:.cc=.o)
MAINDIR = bin
MAIN = $(MAINDIR)/make_lattice

MKDIR = mkdir -p

.PHONY: clean

all:    $(MAIN)
	@echo   The program is compiled. The executable is located in the ./bin directory.

$(MAIN): $(OBJS)
	$(MKDIR) $(MAINDIR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $(MAIN) $(OBJS) $(LDFLAGS) $(LIBS)

.cc.o:
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $<  -o $@

clean:
	$(RM) -r $(MAINDIR) $(SRCDIR)/*.o $(SRCDIR)/*~ $(MAIN)
