CXX = g++
MPICXX = mpic++

STD = -std=c++14
OPTIM = -O2
WALL = -Wall
PROFILE = -g -pg
GDB = -g

SIM5INC = ./sim5/src
SIM5OBJ = ./sim5/src/sim5lib.o
DEBUG = -UDEBUG

# compilation flags
CFLAGS = $(STD) $(OPTIM) $(WALL) -I$(SIM5INC)
LDFLAGS = $(SIM5OBJ) -lstdc++fs #-lCCfits -lcfitsio

# directories
BINDIR = ./bin
OBJDIR = ./obj

.PHONY: clean sim5

# default compilation targets
default: objs bins

# compile sim5 library
sim5:
	$(MAKE) -C sim5

# make OBJ and BIN foldres
$(OBJDIR):
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR)
$(BINDIR):
	@[ -d $(BINDIR) ] || mkdir -p $(BINDIR)

# make object files
OBJS = detector diskobj geoinf gridgeod kerr quadroots scatter tridgeo electron_population electron_population_utils utils photon_dist_utils
OBJFILES = $(patsubst %, $(OBJDIR)/%.o, $(OBJS))
$(OBJS): %: %.cpp | $(OBJDIR)
	$(CXX) $(CFLAGS) $< -c -o $(OBJDIR)/$@.o
objs: sim5 $(OBJS)

# make binaries
bins: calspec 3dcorona 3dcorona_mpi

# make calspec
calspec: $(BINDIR) objs
	$(CXX) $(CFLAGS) calspec.cpp $(OBJFILES) $(LDFLAGS) -o $(BINDIR)/$@

# make 3dcorona
3dcorona: $(BINDIR) objs
	$(CXX) $(CFLAGS) 3dcorona.cpp $(OBJFILES) $(LDFLAGS) -o $(BINDIR)/$@

# make 3dcorona_mpi
3dcorona_mpi: $(BINDIR) objs
	$(MPICXX) $(CFLAGS) 3dcorona_mpi.cpp $(OBJFILES) $(LDFLAGS) -o $(BINDIR)/3dcorona_mpi

clean:
	rm -f $(OBJDIR)/* $(SIM5OBJ)
	rm -f $(BINDIR)/*
