//
//  ant_hbd.hpp
//  AntLatticeWalk
//
//  Created by Yi Hu on 9/30/18.
//

#ifndef ant_hbd_hpp
#define ant_hbd_hpp

#include <string>
#include <unordered_map>
#include <vector>
#include "ant.hpp"

typedef unsigned int mat_limit_type;
typedef std::unordered_map<LeathSiteNode, mat_limit_type> MapType;

class HybridAntWalker: public AntWalker{
public:
  HybridAntWalker(const read_input_ant &inp);
  ~HybridAntWalker();
  int antrun();
  int antdump();
  long getCount(int idx=0);
  std::string getCountStr();
protected:
  mat_limit_type mat_limit;  // limit of matrix size
  LeathSiteNode unitvc[DIM];
  mat_limit_type run_leath(MapType&, std::vector<LeathSiteNode> &);
  int run_static(const MapType&, std::vector<LeathSiteNode> &);
  int run_matrix(const MapType&, std::vector<LeathSiteNode> &);
  int simu_len; // length of simulation
  int run_dyn(MapType&, long);
  bool optnotice_flag = false;
  long count_stat, count_dyn , count_partdyn; // counter on static dynamic partial dynamic runs
};
#endif /* ant_hbd_h */
