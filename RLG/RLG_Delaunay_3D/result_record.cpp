//
//  result_record.cpp
//  VoidPercolation
//
//  Created by Yi Hu on 10/6/18.
//

#include "graph_construct_3d.hpp"

Result_record::Result_record()
{
  n_allfacet = n_checkedfacet = n_allcell = n_checkedcell = 0;
  perco_phi = 0.0;
}

void Result_record::set_value(double s_perco_phi, long s_n_allfacet, long s_n_checkedfacet, long s_n_allcell, long s_n_checkedcell)
{
  perco_phi = s_perco_phi;
  n_allfacet = s_n_allfacet;
  n_checkedfacet = s_n_checkedfacet;
  n_allcell = s_n_allcell;
  n_checkedcell = s_n_checkedcell;
}

void Result_record::print_value(std::ofstream &fout)
{
  fout << std::setprecision(10) << perco_phi << "\t";
  fout << n_allfacet << "\t" << n_checkedfacet << "\t" << n_allcell << "\t" << n_checkedcell << std::endl;
}
