# Dn lattice animal enumeration

CC=g++
FLAGS=--std=c++11 -O3 -IDnLatticeAnimal

all: perimeter sitecount bondcount

perimeter:
	$(CC) $(FLAGS) -o ../bin/lperi DnLatticeAnimal/animal.cpp DnLatticeAnimal/main_perimeter.cpp

sitecount:
	$(CC) $(FLAGS) -o ../bin/lanimal DnLatticeAnimal/animal.cpp DnLatticeAnimal/main_count.cpp

bondcount:
	$(CC) $(FLAGS) -o ../bin/lbnimal DnLatticeAnimal/animal.cpp DnLatticeAnimal/main_bondcount.cpp
