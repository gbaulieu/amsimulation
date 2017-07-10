#include "UnitTest.h"

BOOST_AUTO_TEST_CASE( CMSPatternLayer_constructor_test )
{
  const int NB_LAYERS = 6;
  const int SEGMENT = 1;
  const int MODULE = 13;
  const int LADDER = 14;
  const int SSTRIP = 126;

  GradedPattern p(NB_LAYERS);

  for(int i=0;i<NB_LAYERS;i++){
    CMSPatternLayer test;
    test.setValues(MODULE,LADDER,SSTRIP,SEGMENT);
    p.setLayerStrip(i, &test);
  }

  BOOST_CHECK_EQUAL( p.getNbLayers() , NB_LAYERS );

  for(int i=0;i<NB_LAYERS;i++){
    CMSPatternLayer* pl = dynamic_cast<CMSPatternLayer*>(p.getLayerStrip(i));
    BOOST_CHECK_EQUAL(pl->getModule(),MODULE);
    BOOST_CHECK_EQUAL(pl->getPhi(),LADDER);
    BOOST_CHECK_EQUAL(pl->getStrip(),SSTRIP);
    BOOST_CHECK_EQUAL(pl->getSegment(),SEGMENT);
  }
}


BOOST_AUTO_TEST_CASE( bank_compatibility_test ){

  const int BANK_SIZE = 131072;
  int* pattern100[6];
  int patternLayer0[4] = {0,0,0,1};
  int patternLayer1[4] = {0,0,1,6};
  int patternLayer2[4] = {0,0,2,25};
  int patternLayer3[4] = {0,1,1,6};
  int patternLayer4[4] = {0,1,3,15};
  int patternLayer5[4] = {0,1,3,1};
  
  pattern100[0] = patternLayer0;
  pattern100[1] = patternLayer1;
  pattern100[2] = patternLayer2;
  pattern100[3] = patternLayer3;
  pattern100[4] = patternLayer4;
  pattern100[5] = patternLayer5;

  SectorTree st;
  SectorTree::loadBank(st,"./test_data/test_bank.pbk");
  vector<Sector*> sectors = st.getAllSectors();

  BOOST_CHECK_EQUAL((int)sectors.size(),1);

  for(unsigned int i=0;i<sectors.size();i++){
    Sector* mySector = sectors[i];
    vector<GradedPattern*> patterns = mySector->getPatternTree()->getLDPatterns();

    BOOST_CHECK_EQUAL((int)patterns.size(),BANK_SIZE);

    Pattern* p = patterns[100];
    for(int k=0;k<p->getNbLayers();k++){
      CMSPatternLayer* mp = dynamic_cast<CMSPatternLayer*>(p->getLayerStrip(k));
      BOOST_CHECK_EQUAL(mp->getModule(),pattern100[k][0]);
      BOOST_CHECK_EQUAL(mp->getSegment(),pattern100[k][1]);
      BOOST_CHECK_EQUAL(mp->getPhi(),pattern100[k][2]);
      BOOST_CHECK_EQUAL(mp->getStrip(),pattern100[k][3]);
    }
  }
  
}

BOOST_AUTO_TEST_CASE( pattern_finding_test ){
  SectorTree st;
  SectorTree::loadBank(st,"./test_data/test_bank.pbk");

  ///////////////////////////////////////////////////////////////
  // If we don't have a fitter -> create a TCBuilder default one
  vector<Sector*> sectors = st.getAllSectors();
  for(unsigned int i=0;i<sectors.size();i++){
    if(sectors[i]->getFitter()==NULL){
      //TrackFitter* fitter = new KarimakiTrackFitter(sectors[i]->getNbLayers());
      TrackFitter* fitter = new TCBuilder(sectors[i]->getNbLayers());
      (dynamic_cast<TCBuilder*>(fitter))->setHardwareEmulation(false);
      sectors[i]->setFitter(fitter);
      sectors[i]->updateFitterPhiRotation();
    }
  }
  ///////////////////////////////////////////////////////////////

  PatternFinder pf(0, &st,  "test_data/PU4T_sample.root",  "test_data/output.root");
  {
    pf.useMissingHitThreshold(1);
    int stop = 10;
    pf.find(0, stop);
  }
}
