#include "GradedPattern.h"

GradedPattern::GradedPattern():Pattern(){
  grade=0;
  averagePt=0;
  minPT=200;
  maxPT=0;
  charge = -2;
}

GradedPattern::GradedPattern(const Pattern& p):Pattern(p){
  grade=0;
  averagePt=0;
  minPT=200;
  maxPT=0;
  charge = -2;
}

int GradedPattern::getGrade() const{
  return grade;
}

float GradedPattern::getAveragePt() const{
  return averagePt;
}

float GradedPattern::getMinPt() const{
  return minPT;
}

float GradedPattern::getMaxPt() const{
  return maxPT;
}

int GradedPattern::getCharge() const{
  return charge;
}

void GradedPattern::increment(){
  grade++;
}

void GradedPattern::increment(float pt, int pdg){
  averagePt=(averagePt*grade+pt)/(grade+1);
  if(pt>maxPT)
    maxPT=pt;
  if(pt<minPT)
    minPT=pt;
  int partCharge = (pdg>0)?1:-1;
  if(charge==-2)//this is a new pattern
    charge=partCharge;
  else if(charge!=partCharge)
    charge = 0;
  grade++;
}

int GradedPattern::operator<(const GradedPattern& gp){
  return grade<gp.getGrade();
}
