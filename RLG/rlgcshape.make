ifeq ($(OS),Windows_NT)
    $(error Operation system not support)
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux) # cluster in school
        CC = g++
    endif
    ifeq ($(UNAME_S),Darwin)
        CC = clang++
    endif

    ifneq ($(INCLUDE),)
        COMPFLAGS = -I$(INCLUDE) -D CGAL_EIGEN3_ENABLED=1
    else
        COMPFLAGS = -D CGAL_EIGEN3_ENABLED=1
    endif

    ifneq ($(CGAL),)
        CGALLIBFLAG = -L$(CGAL)/lib64
    endif
endif


CPPFLAGS = -std=c++11 -O3 -Wall -Wno-unused-but-set-variable -Wno-sign-compare
LINKERFLAGS = $(CGALLIBFLAG) -lboost_system -lboost_thread-mt -lCGAL -lgmp

SRCDIR = RLG_CageShape
CMNSRCDIR = shared

BUILDDIR = build/$(SRCDIR)/$(DIM)
CMNBUILDDIR = build/$(CMNSRCDIR)
BINDIR = ../bin
BINNAME = $(BINDIR)/RLG_CegeShape_$(DIM)

SOURCES = cagestat.cpp cellsample.cpp facet_decorate.cpp graph_construct_cage.cpp main.cpp myutility.cpp read_cageshape.cpp
CMNSOURCES = 

OBJECTS := $(addprefix $(BUILDDIR)/,$(SOURCES:.cpp=.o))
CMNOBJECTS := $(addprefix $(CMNBUILDDIR)/,$(CMNSOURCES:.cpp=.o))

DEPENDS = $(OBJECTS:%.o=%.d)
CMNDEPENDS = $(CMNOBJECTS:%.o=%.d)

.PHONY: clean

all: RLG_CageShape

-include $(DEPENDS)
-include $(CMNDEPENDS)

RLG_CageShape: $(OBJECTS) $(CMNOBJECTS)
	$(CC) $(CPPFLAGS) -I../DIM/$(DIM) -o $(BINNAME) $(OBJECTS) $(CMNOBJECTS) $(LINKERFLAGS)

$(OBJECTS): $(BUILDDIR)/%.o: $(SRCDIR)/%.cpp | $(BUILDDIR)
	$(CC) $(CPPFLAGS) $(COMPFLAGS) -I../DIM/$(DIM) -I$(SRCDIR) -I$(CMNSRCDIR) -c $< -o $@

$(CMNOBJECTS): $(CMNBUILDDIR)/%.o: $(CMNSRCDIR)/%.cpp | $(CMNBUILDDIR)
	$(CC) $(CPPFLAGS) -I../DIM/$(DIM) -I$(CMNSRCDIR) -c $< -o $@ 

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(CMNBUILDDIR):
	mkdir -p $(CMNBUILDDIR)

clean:
	rm -rf $(BUILDDIR) $(CMNBUILDDIR) $(BINNAME)
