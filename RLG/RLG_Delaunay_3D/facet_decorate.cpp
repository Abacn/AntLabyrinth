//
//  facet_decorate.c
//  VoidPercolation
//
//  Created by Yi Hu on 8/7/18.
//

// decoration class of facet

#include "graph_construct_3d.hpp"

Facet_decorate::Facet_decorate(double weight, const Cell_handle cell1, const Cell_handle cell2):
weight(weight), cell1(cell1), cell2(cell2)
{
}

Facet_decorate::~Facet_decorate(){}

bool Facet_decorate::operator<(const Facet_decorate &fb) const
{
  return weight < fb.weight;
}
