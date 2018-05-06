CXX = g++
CXXFLAGS  = -O3 -Wall
LDFLAGS =
CPPFLAGS = 
LIBS = -lboost_iostreams

SRCDIR = src
SRCS = $(SRCDIR)/hard_spheres.cc $(SRCDIR)/HSSim.cc $(SRCDIR)/Cell_list.cc $(SRCDIR)/split.cc 
OBJS = $(SRCS:.cc=.o)
MAINDIR = bin
MAIN = $(MAINDIR)/hard_spheres

MKDIR = mkdir -p

.PHONY: clean

all:    $(MAIN)
	@echo   The program $(MAIN) is successfully compiled.

$(MAIN): $(OBJS)
	$(MKDIR) $(MAINDIR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $(MAIN) $(OBJS) $(LDFLAGS) $(LIBS)

.cc.o:
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $<  -o $@

clean:
	$(RM) -r $(MAINDIR) $(SRCDIR)/*.o $(SRCDIR)/*~ $(MAIN)
