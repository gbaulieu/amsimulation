#ifndef _GRADEDPATTERN_H_
#define _GRADEDPATTERN_H_

#include "Pattern.h"

#include <boost/serialization/base_object.hpp>

using namespace std;

/**
   \brief A Pattern wih a grade (the number of tracks corresponding to this pattern).
   Can also contain (optionnaly) the average Pt of the tracks
**/

class GradedPattern : public Pattern{
 public:
  /**
     \brief Constructor
  **/
  GradedPattern();
  /**
     \brief Copy Constructor
  **/
  GradedPattern(const Pattern& p);
  /**
     \brief Get the grade of the Pattern
     \return The number of tracks having generated the pattern
  **/
  int getGrade() const;
  /**
     \brief Get the average Pt of the tracks having generated the pattern (if used)
     \return The average Pt
  **/
  float getAveragePt() const;
  /**
     \brief Get the sign of the generating particles
     \return -1 for PDG<0, 1 for PDG>0, 0 if different sign were used
  **/
  int getSign() const;
  /**
     Increment the grade (tracks occurences + 1)
  **/
  void increment();
  /**
     Increment the grade (tracks occurences + 1), add a Pt value to the average Pt, add the sign according to the PDG
     @param pt The Pt value of the last track
  **/
  void increment(float pt, int pdg);
  /**
     \brief Allows to compare 2 patterns on their grade
     \param gp The second pattern
     \return -1 if the pattern has a lower grade
  **/
  int operator<(const GradedPattern& gp);

 private:
  int grade;
  float averagePt;
  int sign;

  friend class boost::serialization::access;
  
  template<class Archive> void save(Archive & ar, const unsigned int version) const{
    ar << boost::serialization::base_object<Pattern>(*this);
    ar << grade;
    ar << averagePt;
    ar << sign;
  }
  
  template<class Archive> void load(Archive & ar, const unsigned int version){
    ar >> boost::serialization::base_object<Pattern>(*this);
    ar >> grade;
    ar >> averagePt;
    if(version>0){
      ar >> sign;
    }
  }
  
  BOOST_SERIALIZATION_SPLIT_MEMBER()

};
BOOST_CLASS_VERSION(GradedPattern, 1)
#endif
