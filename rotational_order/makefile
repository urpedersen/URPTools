CXX = g++
CXXFLAGS  = -O3 -Wall
LDFLAGS =
CPPFLAGS = -I../cluster_analysis/src
LIBS = -lboost_iostreams

SRCDIR = src
SRCS = $(SRCDIR)/Q6.cc $(SRCDIR)/Rotational_order.cc $(SRCDIR)/Cell_list.cc $(SRCDIR)/split.cc 
OBJS = $(SRCS:.cc=.o)
MAINDIR = bin
MAIN = $(MAINDIR)/Q6

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
