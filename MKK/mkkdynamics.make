
ifeq ($(OS),Windows_NT)
    $(error Operation system not support)
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux) # cluster in school
        ifeq (, $(shell which icpc 2>/dev/null))
            CC = g++
            ifeq (, $(INCLUDE))
                COMPFLG =
            else
                COMPFLG = -I$(INCLUDE)
            endif
            # condor
            ifneq (, $(shell which condor_submit 2>/dev/null))
                LINKERFLAGS = -static
            else
                LINKERFLAGS = 
            endif
        else
            CC = icpc
            ifeq (, $(INCLUDE))
                COMPFLG =
            else
                COMPFLG = -I$(INCLUDE)
            endif
            LINKERFLAGS = 
        endif
    endif
    ifeq ($(UNAME_S),Darwin)
        CC = clang++
        COMPFLG = 
        LINKERFLAGS = 
    endif
endif

CPPFLAGS = -std=c++11 -march=haswell -O2 -Wall -Wno-unused-variable -Wno-overloaded-virtual $(COMPFLG)
SOURCES = neighbor.cpp spheres.cpp box.cpp sphere.cpp event.cpp heap.cpp nlist.cpp read_input.cpp utility.cpp displacements.cpp pressuretensor.cpp

SRCDIR = MKKDynamics
BUILDDIR = build/$(DIM)
BINDIR = ../bin
BINNAME = mkkpoly_$(DIM)
BINFULLNAME = $(BINDIR)/$(BINNAME)
OBJECTS := $(addprefix $(BUILDDIR)/,$(SOURCES:.cpp=.o))
DEPENDS = $(OBJECTS:%.o=%.d)

-include $(DEPENDS)
.PHONY: clean
.DEFAULT_GOAL := $(BINFULLNAME)

$(BINFULLNAME) : $(OBJECTS)
	$(CC) $(CPPFLAGS) $(OBJECTS) -o $(BINFULLNAME)

$(OBJECTS) : $(BUILDDIR)/%.o: $(SRCDIR)/%.cpp | $(BUILDDIR)
	$(CC) $(CPPFLAGS) -I$(SRCDIR) -I../DIM/$(DIM) -MMD -c $< -o $@

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

clean :
	rm -f $(BINFULLNAME) build/$(DIM)/*.o build/$(DIM)/*.d  *.test
