#include "ModuleData.h"

ModuleData::ModuleData(int lay, int lad, int mod){
  layer = lay;
  ladder = lad;
  module = mod;

  firstStrip=10000;
  firstSegment=10000;
  firstX=10000;
  firstY=10000;
  firstZ=10000;
  
  lastStrip=-1;
  lastSegment=-1;
  lastX=-10000;
  lastY=-10000;
  lastZ=-10000;
}

ModuleData::~ModuleData(){}

void ModuleData::update(int seg, int str, float x, float y, float z){
  
  if((seg+str)<=(firstSegment+firstStrip)){
    firstX=x;
    firstY=y;
    firstZ=z;
  }

  if((seg+str)>=(lastSegment+lastStrip)){
    lastX=x;
    lastY=y;
    lastZ=z;
  }

  if(seg<firstSegment)
    firstSegment=seg;
  if(seg>lastSegment)
    lastSegment=seg;

  if(str<firstStrip)
    firstStrip=str;
  if(str>lastStrip)
    lastStrip=str;

  /*
  if(x<firstX)
    firstX=x;
  if(x>lastX)
    lastX=x;

  if(y<firstY)
    firstY=y;
  if(y>lastY)
    lastY=y;

  if(z<firstZ)
    firstZ=z;
  if(z>lastZ)
    lastZ=z;
  */
}

ostream& operator<<(ostream& out, const ModuleData& m){
  out<<m.layer << " / "<< m.ladder <<" / "<<m.module<<" / "
     << (float(m.firstStrip) + float(m.lastStrip))/2.0 << " / "
     << (float(m.firstSegment) + float(m.lastSegment))/2.0   << " / "
     << m.firstX - m.lastX << " / "
     << m.firstY - m.lastY << " / "
     << m.firstZ - m.lastZ << " / "
    //The coordinate of the module middle is the mean of the first and the last coordinates of the module (because inputs coordinates are sorted by strips and segments)
     << (m.firstX + m.lastX)/2.0     << " / "
     << (m.firstY + m.lastY)/2.0     << " / "
     << (m.firstZ + m.lastZ)/2.0     ;
  //out<<"tmp:"<<endl;
  //out<<m.firstStrip<<" / "<<m.lastStrip<<" / "<<m.firstSegment<<" / "<<m.lastSegment<<" / "<<m.firstX<<" / "<<m.lastX<<" / "<<m.firstY<<" / "<<m.lastY<<" / "<<m.firstZ<<" / "<<m.lastZ<<endl;
  return out;
}
