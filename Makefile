CXX = g++
MPICXX = mpic++

STD = -std=c++14
OPTIM = -O2
WALL = -Wall
PROFILE = -g -pg
GDB = -g

MONKDIR = $(HOME)/data/sim5corona/monk
SIM5INC = $(MONKDIR)/src/sim5/src
SIM5OBJ = $(MONKDIR)/src/sim5/src/sim5lib.o
DEBUG = -UDEBUG

# compilation flags
CFLAGS = $(STD) $(OPTIM) $(WALL) -I$(SIM5INC)
LDFLAGS = $(SIM5OBJ) -lstdc++fs #-lCCfits -lcfitsio

# directories
BINDIR = ../bin
OBJDIR = ../obj
TRASHDIR = ~/.trash
INSTALLDIR = ~/.local/bin

.PHONY: bins $(BINS) objs $(OBJS) clean clinstall
# recipes for object files
OBJS = detector diskobj geoinf gridgeod kerr quadroots scatter tridgeo electron_population electron_population_utils utils photon_dist_utils
OBJSRCS = $(patsubst %, %.cpp, $(OBJS))
OBJFILES = $(patsubst %, $(OBJDIR)/%.o, $(OBJS))
objs: $(OBJS)
$(OBJS): %: %.cpp | $(OBJDIR)
	$(CXX) $(CFLAGS) $< -c -o $(OBJDIR)/$@.o

# recipes for binaries
BINS = calspec 3dcorona
BINFILES = $(patsubst %, $(BINDIR)/%, $(BINS))
INSTALLEDBINS = $(patsubst %, $(INSTALLDIR)/%, $(BINS))
$(BINS): %: %.cpp $(OBJSRCS)
	$(CXX) $(CFLAGS) $< $(OBJFILES) $(LDFLAGS) -o $(BINDIR)/$@
	cp $(BINDIR)/$@ $(INSTALLDIR)

# mpi version of 3dcorona
3dcorona_mpi:
	$(MPICXX) $(CFLAGS) 3dcorona_mpi.cpp $(OBJFILES) $(LDFLAGS) -o $(BINDIR)/3dcorona_mpi
	cp $(BINDIR)/3dcorona_mpi $(INSTALLDIR)

clean:
	mv $(OBJDIR)/*.o $(TRASHDIR)
	mv $(BINDIR)/* $(TRASHDIR)
