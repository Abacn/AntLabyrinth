# Dn lattice animal

all: perimeter sitecount bondcount

perimeter:
	clang++ --std=c++11 -O3 -IDnLatticeAnimal -o lperi DnLatticeAnimal/animal.cpp DnLatticeAnimal/main_perimeter.cpp 

sitecount:
	clang++ --std=c++11 -O3 -IDnLatticeAnimal -o lanimal DnLatticeAnimal/animal.cpp DnLatticeAnimal/main_count.cpp 

bondcount:
	clang++ --std=c++11 -O3 -IDnLatticeAnimal -o lbnimal DnLatticeAnimal/animal.cpp DnLatticeAnimal/main_bondcount.cpp 