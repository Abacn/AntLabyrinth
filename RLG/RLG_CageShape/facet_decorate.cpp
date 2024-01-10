//
//  facet_decorate.c
//  VoidPercolation
//
//  Created by Yi Hu on 8/7/18.
//

// decoration class of facet

#include "graph_construct_cage.hpp"
#include "myutility.hpp"

Facet_decorate::Facet_decorate() // uninitialized constructor
{}

Facet_decorate::Facet_decorate(double weight, Full_cell_handle cell): cell(cell), weight(weight)
{
}

Facet_decorate::~Facet_decorate(){}

bool Facet_decorate::operator<(const Facet_decorate &facetB) const
{
  return weight<facetB.weight;
}
