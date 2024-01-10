ifeq ($(OS),Windows_NT)
    $(error Operation system not support)
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        CC = icpc
        ifeq ($(CGAL),)
            $(CGAL not loaded)
        endif
        COMPFLAGS = -I$(INCLUDE)
        CGALLIBFLAG = -L$(CGAL)/lib64
    endif
    ifeq ($(UNAME_S),Darwin)
        CC = clang++
        CGALLIBFLAG = 
        COMPFLAGS = 
    endif
endif


CPPFLAGS = -std=c++11 -O3 -Wall -Wno-unused-but-set-variable -Wno-sign-compare
LINKERFLAGS = $(CGALLIBFLAG) -lboost_system -lboost_thread-mt -lCGAL -lgmp
FILE = large_point_set.cpp

SRCDIR = RLG_Delaunay_3D
CMNSRCDIR = shared

BUILDDIR = build/$(SRCDIR)
CMNBUILDDIR = build/$(CMNSRCDIR)
BINDIR = ../bin
BINNAME = $(BINDIR)/RLG_Delaunay_3D

SOURCES = read_input.cpp cell_disjset.cpp facet_decorate.cpp graph_construct_3d.cpp main.cpp parameters.cpp result_record.cpp traverse_facet.cpp utilities.cpp
CMNSOURCES = 

OBJECTS := $(addprefix $(BUILDDIR)/,$(SOURCES:.cpp=.o))
CMNOBJECTS := $(addprefix $(CMNBUILDDIR)/,$(CMNSOURCES:.cpp=.o))

.PHONY: clean

all: RLG_Delaunay_3D

RLG_Delaunay_3D: $(OBJECTS) $(CMNOBJECTS)
	$(CC) $(CPPFLAGS) -o $(BINNAME) $(OBJECTS) $(CMNOBJECTS) $(LINKERFLAGS)

$(OBJECTS): $(BUILDDIR)/%.o: $(SRCDIR)/%.cpp | $(BUILDDIR)
	$(CC) $(CPPFLAGS) $(COMPFLAGS) -I$(SRCDIR) -I$(CMNSRCDIR) -c $< -o $@

$(CMNOBJECTS): $(CMNBUILDDIR)/%.o: $(CMNSRCDIR)/%.cpp | $(CMNBUILDDIR)
	$(CC) $(CPPFLAGS) -I$(CMNSRCDIR) -c $< -o $@ 

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(CMNBUILDDIR):
	mkdir -p $(CMNBUILDDIR)

clean:
	rm -rf $(BUILDDIR) $(CMNBUILDDIR) $(BINNAME)
