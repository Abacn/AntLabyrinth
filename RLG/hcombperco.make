ifeq ($(OS),Windows_NT)
    $(error Operation system not support)
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux) # cluster in school
        ifeq (, $(shell which icpc 2>/dev/null))
            CC = g++
            # condor
            ifneq (, $(shell which condor_submit 2>/dev/null))
                LINKERFLAGS = -static
            else
                LINKERFLAGS = 
            endif
        else
            CC = icpc
            LINKERFLAGS = 
        endif
        CPPFLAGSB = -std=c++11
        ifeq (, $(INCLUDE))
            COMPFLAGB = 
        else
            COMPFLAGB = -I$(INCLUDE)
        endif
    endif
    ifeq ($(UNAME_S),Darwin)
        CC = clang++
        CPPFLAGSB = -std=c++11
        COMPFLAGB = 
        LINKERFLAGS = 
    endif
endif


CPPFLAGS = -O3 -Wall -Wno-unused-but-set-variable -Wno-sign-compare $(CPPFLAGSB)

SRCDIR = hCombPercolation
DIMNAME = $(DIM:D=)

CMNDIR = qhshared

BUILDDIR = build/$(SRCDIR)/$(DIM)
BINDIR = ../bin
BINNAME = $(BINDIR)/RLG_D$(DIMNAME)
SOURCES = cellset.cpp graph_common.cpp graph_construct.cpp main.cpp read_input.cpp execqhull.cpp commondef.cpp \
 parameters.cpp myutility.cpp
OBJECTS := $(addprefix $(BUILDDIR)/,$(SOURCES:.cpp=.o))
DEPENDS = $(OBJECTS:%.o=%.d)

BUILDDIR_MT = build/$(SRCDIR)_mt/$(DIM)
BINNAME_MT = $(BINDIR)/RLG_mt_D$(DIMNAME)
SOURCES_MT = $(SOURCES) execmt.cpp
CPPFLAGS_MT = $(CPPFLAGS) -D_MULTITHREAD  -pthread

OBJECTS_MT := $(addprefix $(BUILDDIR_MT)/,$(SOURCES_MT:.cpp=.o))
DEPENDS_MT = $(OBJECTS_MT:%.o=%.d)

.PHONY: clean

all: $(BINNAME)

all_mt: $(BINNAME_MT)

-include $(DEPENDS)
-include $(DEPENDS_MT)

# target RLG_qh: calculate threshold
$(BINNAME): $(OBJECTS) | $(BINDIR)
	$(CC) $(CPPFLAGS) $(LINKERFLAGS) -o $(BINNAME) $(OBJECTS)

$(BINNAME_MT): $(OBJECTS_MT) | $(BINDIR)
	$(CC) $(CPPFLAGS) $(CPPFLAGS_MT) $(LINKERFLAGS) -o $(BINNAME_MT) $(OBJECTS_MT)

$(OBJECTS): $(BUILDDIR)/%.o: $(SRCDIR)/%.cpp | $(BUILDDIR)
	$(CC) $(CPPFLAGS) $(COMPFLAGB) -I../DIM/$(DIM) -I$(SRCDIR) -I$(CMNDIR) -MMD -c $< -o $@

$(OBJECTS_MT): $(BUILDDIR_MT)/%.o: $(SRCDIR)/%.cpp | $(BUILDDIR_MT)
	$(CC) $(CPPFLAGS) $(COMPFLAGB) $(CPPFLAGS_MT) -I../DIM/$(DIM) -I$(SRCDIR) -I$(CMNDIR) -MMD -c $< -o $@

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(BUILDDIR_MT):
	mkdir -p $(BUILDDIR_MT)

$(BINDIR):
	mkdir -p $(BINDIR)

clean:
	rm -rf $(BUILDDIR) $(BINNAME)
	rm -rf $(BUILDDIR_MT) $(BINNAME_MT)
