#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE ALL_TESTS tests
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "PatternTree.h"
#include "PatternGenerator.h"
#include "PatternFinder.h"
#include "SectorTree.h"
#include "Detector.h"


#ifndef __APPLE__
BOOST_CLASS_EXPORT_IMPLEMENT(CMSPatternLayer) 
#endif

using namespace std;

