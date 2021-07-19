CXX = g++
MPICXX = mpic++

STD = -std=c++17
OPTIM = -O2
WALL = -Wall
PROFILE = -g -pg
GDB = -g

SIM5INC = $(HOME)/data/sim5corona/monk/src/sim5/src
SIM5OBJ = $(HOME)/data/sim5corona/monk/src/sim5/src/sim5lib.o
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
BINS = calspec 3dcorona ns_surface_write ns_surface nscorona
BINFILES = $(patsubst %, $(BINDIR)/%, $(BINS))
INSTALLEDBINS = $(patsubst %, $(INSTALLDIR)/%, $(BINS))
$(BINS): %: %.cpp $(OBJSRCS)
	$(CXX) $(CFLAGS) $< $(OBJFILES) $(LDFLAGS) -o $(BINDIR)/$@
	cp $(BINDIR)/$@ $(INSTALLDIR)

# mpi version of 3dcorona
3dcorona_mpi:
	$(MPICXX) $(CFLAGS) 3dcorona_mpi.cpp $(OBJFILES) $(LDFLAGS) -o $(BINDIR)/3dcorona_mpi
	cp $(BINDIR)/3dcorona_mpi $(INSTALLDIR)

nscorona_mpi:
	$(MPICXX) $(CFLAGS) nscorona_mpi.cpp $(OBJFILES) $(LDFLAGS) -o $(BINDIR)/nscorona_mpi
	cp $(BINDIR)/nscorona_mpi $(INSTALLDIR)

clean:
	mv $(OBJDIR)/*.o $(TRASHDIR)
	mv $(BINDIR)/* $(TRASHDIR)
