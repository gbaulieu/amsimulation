#include "GradedPattern.h"

GradedPattern::GradedPattern():Pattern(){
  grade=0;
  averagePt=0;
  sign = -2;
}

GradedPattern::GradedPattern(const Pattern& p):Pattern(p){
  grade=0;
  averagePt=0;
  sign = -2;
}

int GradedPattern::getGrade() const{
  return grade;
}

float GradedPattern::getAveragePt() const{
  return averagePt;
}

int GradedPattern::getSign() const{
  return sign;
}

void GradedPattern::increment(){
  grade++;
}

void GradedPattern::increment(float pt, int pdg){
  averagePt=(averagePt*grade+pt)/(grade+1);
  int partSign = (pdg>0)?1:-1;
  if(sign==-2)//this is a new pattern
    sign=partSign;
  else if(sign!=partSign)
    sign = 0;
  grade++;
}

int GradedPattern::operator<(const GradedPattern& gp){
  return grade<gp.getGrade();
}
