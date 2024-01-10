ifeq ($(OS),Windows_NT)
    $(error Operation system not support)
else
    UNAME_S := $(shell uname -s)
	ifeq (, $(shell which icpc 2>/dev/null))
	    ifeq ($(UNAME_S),Linux)
	      CC = g++
		endif
		ifeq ($(UNAME_S),Darwin)
		  CC = clang++
		endif
	else
	    CC = icpc
		MKLROOT = $(INTEL_ROOT)/mkl
		ifeq ($(UNAME_S),Linux)
			MKLFLAGS =  -mkl=sequential -lpthread -lm -ldl
		endif
		ifeq ($(UNAME_S),Darwin)
			MKLFLAGS = -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
		endif
	endif
endif

BINDIR = ../bin
INC = AntLatticeWalk
HBDINC = HybridAnt
IVDINC = InvadePercolation
FCOMMON = $(INC)
SCOMMON = $(INC)/read_input.cpp 
SANT_A = $(INC)/ant_ord.cpp $(INC)/ant_main.cpp $(INC)/ant_global.cpp
SANT_B = $(INC)/ant_und.cpp $(INC)/ant_main.cpp $(INC)/ant_global.cpp
SANT_H = $(INC)/ant_und.cpp $(HBDINC)/ant_hbd.cpp $(INC)/ant_global.cpp $(HBDINC)/ant_hbd_main.cpp
SMSD = $(INC)/cluster_counter.cpp $(INC)/msdlimit.cpp
SLEATH = $(INC)/leath_counter.cpp $(INC)/leath.cpp
SINVADE = $(IVDINC)/invade.cpp $(IVDINC)/basevec.cpp $(IVDINC)/baseleech.cpp

FLAGS = -std=c++14 -Wall -O3 -march=native

all: ant msd leath

ant: ant_ord ant_und

ant_ord: $(BINDIR)
	$(CC) $(FLAGS) -I../DIM/$(DIM) -I$(INC) -o $(BINDIR)/ant_ord_$(DIM) $(SCOMMON) $(SANT_A)

ant_und: $(BINDIR)
	$(CC) $(FLAGS) -I../DIM/$(DIM) -I$(INC) -o $(BINDIR)/ant_und_$(DIM) $(SCOMMON) $(SANT_B)

ant_hbd: $(BINDIR)
	$(CC) $(FLAGS) -I../DIM/$(DIM) -I$(INC) -I$(HBDINC) -I$(MKLROOT)/include -c -o $(BINDIR)/ant_mat_$(DIM).o  $(HBDINC)/ant_mat.cpp
	$(CC) $(FLAGS) -I../DIM/$(DIM) -I$(INC) -I$(HBDINC) -I$(MKLROOT)/include -o $(BINDIR)/ant_hbd_$(DIM) $(BINDIR)/ant_mat_$(DIM).o $(SCOMMON) $(SANT_H) $(MKLFLAGS)

msd: $(BINDIR)
	$(CC) $(FLAGS) -I../DIM/$(DIM) -I$(INC) -o $(BINDIR)/msdlim_$(DIM) $(SCOMMON) $(SMSD)

leath: $(BINDIR)
	$(CC) $(FLAGS) -I../DIM/$(DIM) -I$(INC) -o $(BINDIR)/leath_$(DIM) $(SCOMMON) $(SLEATH)

invade: $(BINDIR)
	$(CC) $(FLAGS) -I../DIM/$(DIM) -I$(INC) -I$(IVDINC) -o $(BINDIR)/invade_$(DIM) $(SCOMMON) $(SINVADE)

usage:
	@echo "Usage make [ant/msd] DIM=3D"

$(BINDIR):
	mkdir -p $(BINDIR)

clean:
	rm -f *.o

.PHONY: ant_ord ant_und ant_hbd msd leath usage clean

