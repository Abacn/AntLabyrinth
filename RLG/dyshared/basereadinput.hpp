//
//  basereadinput.hpp
//  RLGDynamics
//

#ifndef basereadinput_hpp
#define basereadinput_hpp

#include <cstdint>
#include <string>

#include <boost/algorithm/string.hpp>

#define NAME_LEN 256

enum System_type{
  SYSTEM_TYPE_BOX,
  SYSTEM_TYPE_HCOMB,
  SYSTEM_TYPE_SPHERE,
  SYSTEM_TYPE_ERROR
};

enum Output_mode{
  OUTPUT_MODE_SUMMARY,
  OUTPUT_MODE_PRODUCT,
  OUTPUT_MODE_SAVEGRAPH,
  OUTPUT_MODE_DEBUG,
  OUTPUT_MODE_ALL
};

class ReadInput {
public:
  virtual int read(const char *inputf) = 0;
protected:
  static System_type getSystemType(const std::string stype)
  {
    if(boost::iequals(stype, "box"))
    {
      return SYSTEM_TYPE_BOX;
    }
    else if(boost::iequals(stype, "hcomb"))
    {
      return SYSTEM_TYPE_HCOMB;
    }
    else if(boost::iequals(stype, "sphere"))
    {
      return SYSTEM_TYPE_SPHERE;
    }
    else
    {
      return SYSTEM_TYPE_ERROR;
    }
  }

  static Output_mode getOutputMode(const std::string omode)
  {
    if(boost::iequals(omode, "debug"))
    {
      return OUTPUT_MODE_DEBUG;
    }
    else if(boost::iequals(omode, "savegraph"))
    {
      return OUTPUT_MODE_SAVEGRAPH;
    }
    else if(boost::iequals(omode, "product"))
    {
      return OUTPUT_MODE_PRODUCT;
    }
    else if(boost::iequals(omode, "all"))
    {
      return OUTPUT_MODE_ALL;
    }
    else
    {
      return OUTPUT_MODE_PRODUCT;
    }
  }
};

#endif /* basereadinput_hpp */
