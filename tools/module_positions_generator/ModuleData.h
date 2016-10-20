#ifndef _MODULEDATA_H_
#define _MODULEDATA_H_

#include <iostream>
#include <ostream>

using namespace std;

class ModuleData{

 private:
  int layer;
  int ladder;
  int module;

  int firstStrip;
  int firstSegment;
  float firstX;
  float firstY;
  float firstZ;
  
  int lastStrip;
  int lastSegment;
  float lastX;
  float lastY;
  float lastZ;

 public:
  ModuleData(int lay, int lad, int mod);
  ~ModuleData();
  void update(int seg, int str, float x, float y, float z);
  friend ostream& operator<<(ostream& out, const ModuleData& m);
};
#endif
