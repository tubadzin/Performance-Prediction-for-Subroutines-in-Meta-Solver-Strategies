# ──────────────────────────────────────────────────────────────────────
# Makefile for MIPfeature on Apple Silicon macOS with CPLEX Studio 2211
# on mac: brew install boost needed.
# changes: replaced: CPXFILEptr, CPXfopen, CPXgetclqcnt, CPXgetcovcnt
# ──────────────────────────────────────────────────────────────────────
BOOSTINCDIR   = /opt/homebrew/include

# target system and lib type
SYSTEM      = arm64_osx
LIBFORMAT   = static_pic

# point these at your CPLEX-Studio install
CPLEXDIR      = /Applications/CPLEX_Studio2211/cplex
CONCERTDIR    = /Applications/CPLEX_Studio2211/concert

# include dirs
CPLEXINCDIR   = $(CPLEXDIR)/include
CONCERTINCDIR = $(CONCERTDIR)/include

# lib dirs (static .a libraries under arm64_osx/static_pic)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

# linker flags
CCLNFLAGS = \
  -L$(CPLEXLIBDIR) -lilocplex -lcplex \
  -L$(CONCERTLIBDIR) -lconcert \
  -lm -lpthread

# objects and executable
OBJ      = comp_feat.o stopwatch.o
EXECNAME = MIPfeature

# optimization flags
RELEASEFLAGS = -O3 -DNDEBUG

# compile rule
.cpp.o:
	clang++ -g -c $< \
	  -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I$(BOOSTINCDIR)


# link rule
all: $(OBJ)
	clang++ -o $(EXECNAME) $(OBJ) $(CCLNFLAGS) $(RELEASEFLAGS)

clean:
	rm -f $(OBJ) $(EXECNAME)
