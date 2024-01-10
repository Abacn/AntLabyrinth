ifeq ($(OS),Windows_NT)
    $(error Operation system not support)
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux) # cluster in school
        ifeq (, $(shell which icpc 2>/dev/null))
            CC = g++
            C0 = gcc
            # condor
            ifneq (, $(shell which condor_submit 2>/dev/null))
                LINKERFLAGS = -static
            else
                LINKERFLAGS = 
            endif
        else
            CC = icpc
            C0 = icc
            LINKERFLAGS = 
        endif
        CPPFLAGSB = -std=c++11 -Wno-unused-but-set-variable
        ifeq (, $(INCLUDE))
            COMPFLAGB = 
        else
            COMPFLAGB = -I$(INCLUDE)
        endif
    endif
    ifeq ($(UNAME_S),Darwin)
        CC = clang++
        C0 = clang
        LINKERFLAGS = 
        CPPFLAGSB = -std=c++11
        COMPFLAGB = 
    endif
endif

CPPFLAGS = -O3 -march=haswell -Wall -Wno-sign-compare
C0FLAGS = -O3 -march=haswell -Wall -std=c99
SRCDIR = RLGDynamics
CMNSRCDIR = dyshared
SQSRCDIR = SqWellDynamics

BUILDDIR = build/$(SRCDIR)/$(DIM)
CMNBUILDDIR = build/$(CMNSRCDIR)/$(DIM)
SQBUILDDIR = build/$(SQSRCDIR)/$(DIM)
BINDIR = ../bin

SOURCES = box.cpp main.cpp msdrecorder.cpp gridlist.cpp read_input.cpp hcomb.cpp hcombpbc.cpp sphere.cpp
CMNSOURCES = commondef.cpp myutility.cpp brownian.cpp

OBJECTS := $(addprefix $(BUILDDIR)/,$(SOURCES:.cpp=.o))
CMNOBJECTS := $(addprefix $(CMNBUILDDIR)/,$(CMNSOURCES:.cpp=.o))

SQSOURCES = sqwellbox.cpp main.cpp
SQSHAREDSOURCES = msdrecorder.cpp gridlist.cpp

SQOBJECTS := $(addprefix $(SQBUILDDIR)/,$(SQSOURCES:.cpp=.o))
SQSHAREDOBJECTS := $(addprefix $(BUILDDIR)/,$(SQSHAREDSOURCES:.cpp=.o))

DEPENDS = $(OBJECTS:%.o=%.d) $(CMNOBJECTS:%.o=%.d) $(SQOBJECTS:%.o=%.d)

ifeq (24D,$(DIM))
    LEECHOBJECT=$(BUILDDIR)/leech.o
endif

.PHONY: clean clean_sq

all: RLG_dy

-include $(DEPENDS)

clean_all: clean
# target RLG_dy: calculate threshold
RLG_dy: $(CMNOBJECTS) $(OBJECTS) $(LEECHOBJECT) | $(BINDIR)
	$(CC) $(CPPFLAGS) $(CPPFLAGSB) -o $(BINDIR)/RLG_dy_$(DIM) $(CMNOBJECTS) $(OBJECTS) $(LEECHOBJECT) $(LINKERFLAGS)

SQRLG_dy: $(CMNOBJECTS) $(SQSHAREDOBJECTS) $(SQOBJECTS) | $(BINDIR)
	$(CC) $(CPPFLAGS) $(CPPFLAGSB) -o $(BINDIR)/SQRLG_dy_$(DIM) $(CMNOBJECTS) $(SQSHAREDOBJECTS) $(SQOBJECTS) $(LINKERFLAGS)

$(OBJECTS): $(BUILDDIR)/%.o: $(SRCDIR)/%.cpp | $(BUILDDIR)
	$(CC) $(CPPFLAGS) $(CPPFLAGSB) $(COMPFLAGB) -I$(SRCDIR) -I$(CMNSRCDIR) -I../DIM/$(DIM) -MMD -c $< -o $@

$(CMNOBJECTS): $(CMNBUILDDIR)/%.o: $(CMNSRCDIR)/%.cpp | $(CMNBUILDDIR)
	$(CC) $(CPPFLAGS) $(CPPFLAGSB) $(COMPFLAGB) -I../DIM/$(DIM) -MMD -c $< -o $@ 

$(SQOBJECTS): $(SQBUILDDIR)/%.o: $(SQSRCDIR)/%.cpp | $(SQBUILDDIR)
	$(CC) $(CPPFLAGS) $(CPPFLAGSB) $(COMPFLAGB) -I$(SRCDIR) -I$(CMNSRCDIR) -I$(SQSRCDIR) -I../DIM/$(DIM) -MMD -c $< -o $@

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(CMNBUILDDIR):
	mkdir -p $(CMNBUILDDIR)

$(SQBUILDDIR):
	mkdir -p $(SQBUILDDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

build/RLGDynamics/24D/leech.o: RLGDynamics/leech.c | build/RLGDynamics/24D
	$(C0) $(C0FLAGS) -o $@ -c $<

clean:
	rm -rf $(BUILDDIR) $(CMNBUILDDIR) $(DEPENDS) $(BINDIR)/RLG_dy_$(DIM)

clean_sq:
	rm -rf $(BUILDDIR) $(CMNBUILDDIR) $(SQBUILDDIR) $(DEPENDS) $(BINDIR)/SQRLG_dy_$(DIM)
