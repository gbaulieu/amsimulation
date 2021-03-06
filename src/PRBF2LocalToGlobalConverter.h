#ifndef _PRBF2LOCALTOGLOBALCONVERTER_H_
#define _PRBF2LOCALTOGLOBALCONVERTER_H_

#include "LocalToGlobalConverter.h"
#include "Sector.h"
#include "CommonTools.h"

using namespace std;

/**
   \brief Converts local coordinates (layer, ladder, module, segment, hit local to the sector) to global cartesian coordinates
 **/
class PRBF2LocalToGlobalConverter:public LocalToGlobalConverter{

 private:
  /**
     module_pos[prbf2_layer][prbf2_ladder][prbf2_module] contains 
     [
     module_center_X, 
     module_center_Y, 
     module_center_Y, 
     strip_pitch_X, 
     strip_pitch_Y, 
     strip_pitch_Z,
     seg_pitch_X, 
     seg_pitch_Y, 
     seg_pitch_Z
     ]
  **/
  vector<float> ***module_pos;
  const Sector* sector;

 public:
  /**
     \brief Constructor
     \param sectorDefinition The sector concerned by the convertion
     \param geometryFile Text file containing the global coordinates of all modules
  **/
  PRBF2LocalToGlobalConverter(const Sector* sectorDefinition, string geometryFile);

  /**
     \brief Destructor
  **/
  ~PRBF2LocalToGlobalConverter();

  /**
     \brief Convertion of local coordinates (from PRBF2 format) to cartesian coordinates global to the tracker
     \param layer The PRBF2 value of the layer (0 to 15)
     \param ladder The PRBF2 value of the ladder (position of the ladder in the trigger tower)
     \param module The PRBF2 value of the module (position of the module in the ladder for the current trigger tower)
     \param segment ThE PRBF2 CIC bit for 2S modules, the PRBF2 CIC bit concatenated with the PRBF2 Z bits value for PS modules
     \param strip The PRBF2 strip value
     \return A vector of 3 float containing the X/Y/Z cartesian coordinates
   **/
  vector<float> toGlobal(int layer, int ladder, int module, int segment, float strip) const throw (std::runtime_error);

  /**
     \brief Convertion of local coordinates contained in a Hit object to cartesian coordinates global to the tracker
     \param h A Hit object
     \return A vector of 3 float containing the X/Y/Z cartesian coordinates
   **/
  vector<float> toGlobal(const Hit* h) const throw (std::runtime_error);

};
#endif
