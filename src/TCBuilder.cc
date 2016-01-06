/**
C++ implementation of the TC builder

S.Viret, G.Galbit : 09/12/15
**/

#include "TCBuilder.h"

TCBuilder::TCBuilder():TrackFitter(0){

}

TCBuilder::TCBuilder(int nb):TrackFitter(nb)
{

  //Layer infos
  m_nLayer = 6;

  //
  // Initialisation of the matching windows
  //

  for (int i=0;i<3;++i)
  {
    for (int j=0;j<6;++j)
    {
      fenetre_b_z[i][j]=-1;
      fenetre_b_phi[i][j]=-1;
    }    

    for (int j=0;j<15;++j)
    {
      fenetre_e_r[i][j]=-1;
      fenetre_e_phi[i][j]=-1;
    } 
  }

  // Windows for the barrel towers
  // 1 z window per layer (in cm)
  // 1 phi window per layer (in rad)

  fenetre_b_z[2][0]  = 0.28;  
  fenetre_b_z[2][1]  = 0.14;  
  fenetre_b_z[2][2]  = 0.34;  
  fenetre_b_z[2][3]  = 2.8;  
  fenetre_b_z[2][4]  = 3.0;  
  fenetre_b_z[2][5]  = 3.2;  
  fenetre_b_phi[2][0]= 0.0026;  
  fenetre_b_phi[2][1]= 0.0013;  
  fenetre_b_phi[2][2]= 0.003;  
  fenetre_b_phi[2][3]= 0.0065;  
  fenetre_b_phi[2][4]= 0.0115;  
  fenetre_b_phi[2][5]= 0.016;  
   
  // Windows for the endcap towers
  // For barrel layers
  // 1 z window per layer (in cm)
  // 1 phi window per layer (in rad)

  fenetre_b_z[0][0]  = 1.6;  
  fenetre_b_z[0][1]  = 0.9;  
  fenetre_b_phi[0][0]= 0.0038;  
  fenetre_b_phi[0][1]= 0.0016;  
 
  // For disks
  // 1 r window per ring (in cm)
  // 1 phi window per ring (in rad)
  fenetre_e_r[0][0]  = 0.1;  
  fenetre_e_r[0][1]  = 0.25;  
  fenetre_e_r[0][2]  = 0.33;  
  fenetre_e_r[0][3]  = 0.38;  
  fenetre_e_r[0][4]  = 0.46;  
  fenetre_e_r[0][5]  = 0.5;  
  fenetre_e_r[0][6]  = 0.65;  
  fenetre_e_r[0][7]  = 0.7;  
  fenetre_e_r[0][8]  = 0.75;
  fenetre_e_r[0][9]  = 2.9;
  fenetre_e_r[0][10] = 2.8;
  fenetre_e_r[0][11] = 3.0;
  fenetre_e_r[0][12] = 3.0;
  fenetre_e_r[0][13] = 3.1;
  fenetre_e_r[0][14] = 3.2;
  fenetre_e_phi[0][0]= 0.0006;  
  fenetre_e_phi[0][1]= 0.0018;  
  fenetre_e_phi[0][2]= 0.0025;  
  fenetre_e_phi[0][3]= 0.0038;  
  fenetre_e_phi[0][4]= 0.005;  
  fenetre_e_phi[0][5]= 0.0055;  
  fenetre_e_phi[0][6]= 0.007;  
  fenetre_e_phi[0][7]= 0.008;  
  fenetre_e_phi[0][8]= 0.008; 
  fenetre_e_phi[0][9]= 0.010; 
  fenetre_e_phi[0][10]= 0.010; 
  fenetre_e_phi[0][11]= 0.012; 
  fenetre_e_phi[0][12]= 0.012; 
  fenetre_e_phi[0][13]= 0.015; 
  fenetre_e_phi[0][14]= 0.017; 

 // Windows for the hybrid towers

  fenetre_b_z[1][0]  = 0.28;  
  fenetre_b_z[1][1]  = 0.14;  
  fenetre_b_z[1][2]  = 0.34;  
  fenetre_b_z[1][3]  = 2.8;  
  fenetre_b_z[1][4]  = 3.0;  
  fenetre_b_z[1][5]  = 3.2;  
  fenetre_b_phi[1][0]= 0.0026;  
  fenetre_b_phi[1][1]= 0.0013;  
  fenetre_b_phi[1][2]= 0.003;  
  fenetre_b_phi[1][3]= 0.0065;  
  fenetre_b_phi[1][4]= 0.0115;  
  fenetre_b_phi[1][5]= 0.016; 
  
  fenetre_e_r[1][3]  = 0.08;  
  fenetre_e_r[1][4]  = 0.09;  
  fenetre_e_r[1][5]  = 0.11;  
  fenetre_e_r[1][6]  = 0.12;  
  fenetre_e_r[1][7]  = 0.14;  
  fenetre_e_r[1][8]  = 0.16;
  fenetre_e_r[1][9]  = 2.45;
  fenetre_e_r[1][10] = 2.6;
  fenetre_e_r[1][11] = 2.65;
  fenetre_e_r[1][12] = 2.7;
  fenetre_e_r[1][13] = 2.8;
  fenetre_e_r[1][14] = 3.0;
  fenetre_e_phi[1][3]= 0.00065;  
  fenetre_e_phi[1][4]= 0.00095;  
  fenetre_e_phi[1][5]= 0.0018;  
  fenetre_e_phi[1][6]= 0.0022;  
  fenetre_e_phi[1][7]= 0.0033;  
  fenetre_e_phi[1][8]= 0.0036; 
  fenetre_e_phi[1][9]= 0.007; 
  fenetre_e_phi[1][10]= 0.0095; 
  fenetre_e_phi[1][11]= 0.012; 
  fenetre_e_phi[1][12]= 0.012; 
  fenetre_e_phi[1][13]= 0.014; 
  fenetre_e_phi[1][14]= 0.016; 

  // cout<<"TC builder init done"<<endl;
}

TCBuilder::~TCBuilder()
{
}

void TCBuilder::initialize(){

}

void TCBuilder::mergePatterns(){
  //cout<<"Merging of patterns not implemented"<<endl;
}

// Tentative duplicate removal
// Not used for the moment

void TCBuilder::mergeTracks(){
  cout<<"We have "<<tracks.size()<<" tracks before merging..."<<endl;
  unsigned int index = 0;
  vector<Track*>::iterator it = tracks.begin();
  while(it!=tracks.end()){
    Track* newTrack = *it;
    bool found = false;
    for(unsigned int i=0;i<index;i++){
      Track* ref = tracks[i];
      float dpt,dphi,dz,deta;
      dpt = fabs(newTrack->getCurve()-ref->getCurve());
      dphi = fabs(newTrack->getPhi0()-ref->getPhi0());
      dz = fabs(newTrack->getZ0()-ref->getZ0());
      deta = fabs(newTrack->getEta0()-ref->getEta0());
      found = (deta<0.02) &&
	(dphi<0.005) &&
	(dpt<0.1) &&
	(dz<0.3);
      if(found)
	break;
    }
    if(found)
      tracks.erase(it);
    else{
      index++;
      it++;
    }
  }
  cout<<"We have "<<tracks.size()<<" tracks after merging..."<<endl;
}

Track* TCBuilder::createFittedTrack(vector <Hit*> &bestTC)
{
  int lay;
  bool barrelmod;

  int size = bestTC.size();
    
  float r;
  float rmin,rmax;
  
  float phi_est,z_est,eta_est,pt_est;
  
  double parZX[2][2];
  double resZX[2];
  double invZX[2][2];
  double detZX = 0;
     
  double x,y,z;
  
  int npt=0;
  
  float x1,x2,y1,y2;
  float x0,y0;
  
  for (int i=0;i<2;++i)
  {
    resZX[i] = 0.;
    for (int jj=0;jj<2;++jj) parZX[i][jj] = 0.;
    for (int jj=0;jj<2;++jj) invZX[i][jj] = 0.;
  }
    
  rmin = 1000;
  rmax = 0;
  int imax=-1;
  int imin=-1;
  
  x0=0;
  y0=0;
  
  x2=0;
  y2=0;

  int n2Sb=0;
  int nPSb=0;
  int n2Se=0;
  int nPSe=0;

  // Loop over stubs in the TC
  // In order to determine the lowest/largest radius of the 
  // TC, and the TC composition
  //
  // The stub with the largest radius is our reference in 
  // the conformal space


  for (int kk=0;kk<size;++kk) 
  {
    lay  = bestTC.at(kk)->getLayer();
    
    if (lay<=7)
    {
      (lay<=2)
	    ? nPSb++
	    : nPSe++;
    }
    else
    {
      (lay<11)
	    ? n2Sb++
	    : n2Se++;
    }
            
    if (lay>10) continue; // 2S disks stubs not used
    ++npt;
        
    x = bestTC.at(kk)->getX();
    y = bestTC.at(kk)->getY();
    r = sqrt(x*x+y*y);
    
    if (r<rmin)
    {
      rmin = r;
      imin = kk;
      x2   = x;
      y2   = y;
    }
        
    if (r>rmax)
    {
      rmax = r;
      imax = kk;
      x0   = x;
      y0   = y;
    }
  }

  float rmax2 = 0;
    
  x1 = 0;
  y1 = 0;
  
  int nc=0;
    
  float xtemp1,ytemp1;
  float xtemp2,ytemp2;
  
  xtemp1=0.;
  ytemp1=0.;

  // Loop 2 over stubs in the TC
  // In order to determine the point with the second largest radius 

  for (int kk=0;kk<size;++kk) // Loop over stubs in the TC
  {
    if (kk==imax || kk==imin) continue; // Already used

    lay  = bestTC.at(kk)->getLayer();

    barrelmod=0;
    if (lay<=2 || (lay>=8 && lay<=10)) barrelmod=1;
    if (!barrelmod && (nPSb+n2Sb)>=3) continue; // Barrel modules have a priority
    if (lay>10 && (nPSb+nPSe)>=3) continue;     // Don't use 2S modules in the disks if possible
    
    x = bestTC.at(kk)->getX();
    y = bestTC.at(kk)->getY();
    r = sqrt(x*x+y*y);
    
    if (r>rmax2)
    {
      rmax2  = r;
      x1     = x;
      y1     = y;
      nc++;
    }
  }

  // Now get the coordinates in the conformal space.

  xtemp1 = (x1-x0)/((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
  ytemp1 = (y1-y0)/((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
  
  xtemp2 = (x2-x0)/((x2-x0)*(x2-x0)+(y2-y0)*(y2-y0));
  ytemp2 = (y2-y0)/((x2-x0)*(x2-x0)+(y2-y0)*(y2-y0));
  
  x1=xtemp1;
  y1=ytemp1;
  x2=xtemp2;
  y2=ytemp2;
  

  // Now we got everything for the r/phi plane   

  double a = x0-0.5*(y1-y2)/(y2*x1-y1*x2);
  double b = y0-0.5*(x1-x2)/(y1*x2-y2*x1);
    
  int charge =-b/fabs(b);
   
  phi_est = atan2(charge*a,-charge*b)+sec_phi;
  pt_est  = 0.003*3.833*sqrt(a*a+b*b);

  // Then we do the RZ fit (LS)

  float wght;
  int cnt=0;
  for (int kk=0;kk<size;++kk) // Loop over stubs in the TC
  {
    lay  = bestTC.at(kk)->getLayer();

    if (lay>7) continue; // Don't use 2S modules
    if (lay>2 && nPSb>=2) continue; // Don't use PS modules of the disks if number of good point in the barrel is sufficient        

    ++cnt;
    x = bestTC.at(kk)->getX();
    y = bestTC.at(kk)->getY();
    z = bestTC.at(kk)->getZ();        
    r = sqrt(x*x+y*y);
            
    wght=1;
    if (lay>2) wght=1/7.;

    parZX[0][0] += wght*r*r;
    parZX[1][1] += wght*1;
    parZX[1][0] += wght*r;
    
    resZX[0] += wght*r*z;
    resZX[1] += wght*z;       
  } // End of stub loop
    

  detZX = parZX[0][0]*parZX[1][1]-parZX[1][0]*parZX[1][0];

  invZX[0][0] =  parZX[1][1]/detZX;
  invZX[1][0] = -parZX[1][0]/detZX;
  invZX[1][1] =  parZX[0][0]/detZX;

  // Finally estimate track parameters in the R/Z plane 

  eta_est = std::asinh((invZX[0][0]*resZX[0] + invZX[1][0]*resZX[1]));
  z_est   = invZX[1][0]*resZX[0] + invZX[1][1]*resZX[1];

  Track* fit_track = new Track();
  fit_track->setCurve(pt_est);
  fit_track->setPhi0(phi_est);
  fit_track->setEta0(eta_est);
  fit_track->setZ0(z_est);
      
  for(unsigned int hitIndex=0;hitIndex < bestTC.size();hitIndex++)
  {
    fit_track->addStubIndex(bestTC[hitIndex]->getID());
  }

  return fit_track;
}

float TCBuilder::binning(float fNumber, int nFractionnalPartWidth)
{
  if (nFractionnalPartWidth == 0)
		//If this parameter is set to zero, no binning is applied
		return fNumber;
		
	//The number is divided by the power of 2 corresponding to the fractionnal part width
	float divRes = fNumber / pow(2 , -nFractionnalPartWidth);

	//The result is rounded to the nearest integer and then multiplied by the power of 2
	float fBinnedNumber = round(divRes) * pow(2 , -nFractionnalPartWidth);
	
	return fBinnedNumber;
}

void TCBuilder::alignScore(Hit& hSeed1, Hit& hSeed2, Hit& hTestStub, float tScores[], int nFractionnalPartWidth)
{
  float fRPHI_Score , fRZ_Score;

  float X1, Y1, Z1, R1, PHI1;
  float X2, Y2, Z2, R2, PHI2;
  float X3, Y3, Z3, R3, PHI3;

  float RPHI_S1, RPHI_S2, RZ_S1, RZ_S2;

  X1 = binning(hSeed1.getX(), nFractionnalPartWidth);
	Y1 = binning(hSeed1.getY(), nFractionnalPartWidth);
  Z1 = binning(hSeed1.getZ(), nFractionnalPartWidth);
	
	X2 = binning(hSeed2.getX(), nFractionnalPartWidth);
	Y2 = binning(hSeed2.getY(), nFractionnalPartWidth);
	Z2 = binning(hSeed2.getZ(), nFractionnalPartWidth);

	X3 = binning(hTestStub.getX(), nFractionnalPartWidth);
	Y3 = binning(hTestStub.getY(), nFractionnalPartWidth);
  Z3 = binning(hTestStub.getZ(), nFractionnalPartWidth);
	
	R1 = binning(sqrt(X1*X1 + Y1*Y1), nFractionnalPartWidth);
	R2 = binning(sqrt(X2*X2 + Y2*Y2), nFractionnalPartWidth);
	R3 = binning(sqrt(X3*X3 + Y3*Y3), nFractionnalPartWidth);

	//RPHI plan
  PHI1 = binning(atan(Y1/X1), nFractionnalPartWidth);
  PHI2 = binning(atan(Y2/X2), nFractionnalPartWidth);
  PHI3 = binning(atan(Y3/X3), nFractionnalPartWidth);

  RPHI_S1 = binning((PHI2 - PHI1) * (R3 - R2), nFractionnalPartWidth);
  RPHI_S2 = binning((PHI2 - PHI3) * (R2 - R1), nFractionnalPartWidth);

  fRPHI_Score = binning(RPHI_S1 + RPHI_S2, nFractionnalPartWidth);

  //RZ plan

  RZ_S1 = binning((Z2 - Z1) * (R3 - R2), nFractionnalPartWidth);
  RZ_S2 = binning((Z2 - Z3) * (R2 - R1), nFractionnalPartWidth);

  fRZ_Score = binning(RZ_S1 + RZ_S2, nFractionnalPartWidth);

  tScores[0] = fRPHI_Score;
  tScores[1] = fRZ_Score;
}

void TCBuilder::getThresholds(int nLaySeed1, int nLaySeed2, int nLayTestStub, float tabThresh[])
{
  tabThresh[0] = 25.0;
  tabThresh[1] = 25.0;
}

char TCBuilder::transcodeLayer(Hit * pHit)
{
  int nOrigLayer = pHit->getLayer();

  //layer transcoding of the disks is based on the radius, it can be optimized by using the ladder_ID
  int X = pHit->getX();
  int Y = pHit->getY();

  int nTransLayer;
		
	if (nOrigLayer <= 10)
	{
		//If the stub is in the barrel

		if (nOrigLayer <= 7)
    {    
      //If layer 5, 6, 7
      nTransLayer = nOrigLayer - 5;     //5->0, 6->1, ...
    }    
    else
    {
      //If layer 8, 9, 10
      nTransLayer = nOrigLayer;         //no change
    }
  }
	else if (sqrt(X*X + Y*Y) >= 62)
	{
		//If the stub is on an outer ring of a disk !!! (2S modules)

		if (nOrigLayer <= 15)
    {
      //If layer 11, 12, 13, 14, 15
      nTransLayer = nOrigLayer;         //no change
    }
    else
    {
      //If layer 18, 19, 20, 21, 22
      nTransLayer = nOrigLayer - 7;     //18->11, 19->12, ...
    }
	}
	else
  {
		//If the stub is on an inner ring of a disk !!! (PS modules)

		if (nOrigLayer <= 15)
    {
      //If layer 11, 12, 13, 14, 15
      nTransLayer = nOrigLayer - 8;     //11->3, 12->4, ...
    }
    else
    {
      //If layer 18, 19, 20, 21, 22
      nTransLayer = nOrigLayer - 15;    //18->3, 19->4, ...
    }
	}

  return nTransLayer;	
}

// TC builder module
//
// Take as input the list of stubs contained in a matched road

void TCBuilder::fit(vector<Hit*> originalHits)
{
  int nFractionnalPartWidth = 0;
  int nMissingHits = 1;

  cout<<"trying to fit "<<originalHits.size()<<" points"<<endl;

  int tow = sector_id; // The tower ID, necessary to get the phi shift
  
  //Process the starting phi of the tower
  float sec_phi = (tow%8) * M_PI / 4 - 0.4;

  //cos and sin values for a rotation of an angle -sec_phi
  float ci = cos(-sec_phi);
  float si = sin(-sec_phi);

  float rotatedX, rotatedY;

  //Create a new vector to store the custom hits parameters
  vector <Hit> hits;
  
  Hit* pOrigHit;
  
  //For each hit of the lists
  for (unsigned int origHitIndex = 0; origHitIndex<originalHits.size(); origHitIndex++)
  {
    pOrigHit = originalHits[origHitIndex];

    //Process the rotated coordinnates
    rotatedX = pOrigHit->getX() * ci - pOrigHit->getY() * si;
    rotatedY = pOrigHit->getX() * si + pOrigHit->getY() * ci;
    
    //Add the modified hit to the hits vector
    hits.push_back( Hit(transcodeLayer(pOrigHit),
                        pOrigHit->getLadder(),
                        pOrigHit->getModule(),
                        pOrigHit->getSegment(),
                        pOrigHit->getStripNumber(),
                        pOrigHit->getID(),
                        pOrigHit->getParticuleID(),
                        pOrigHit->getParticulePT(),
                        pOrigHit->getParticuleIP(),
                        pOrigHit->getParticuleETA(),
                        pOrigHit->getParticulePHI0(),
                        binning(rotatedX, nFractionnalPartWidth),
                        binning(rotatedY, nFractionnalPartWidth),
                        binning(pOrigHit->getZ(), nFractionnalPartWidth),
                        pOrigHit->getX0(),
                        pOrigHit->getY0(),
                        pOrigHit->getZ0(),
                        pOrigHit->getBend())
                        );

  }

  //Sort the hits by ascending order of layerID
  //(using a lambda definition of the sorting criteria which return a boolean)
  sort(hits.begin(), hits.end(), [ ]( const Hit& lhs, const Hit& rhs ) { return lhs.getLayer() < rhs.getLayer(); });

  
  int nLayersCurrentPattern = 0;
  int lastAddedLayer = -1;
  //Count the number of layers present in the pattern
  for (unsigned int hitIndex=0; hitIndex < hits.size(); hitIndex++)
  {
    if (lastAddedLayer != hits[hitIndex].getLayer())
    {
      nLayersCurrentPattern++;
      lastAddedLayer = hits[hitIndex].getLayer();
    }
  }


  vector <Hit*> vecCurrentCandidateHits;
  vector <float> vecCurrentCandidateScore;

  vector <Hit*> vecBestCandidateHits;
  float fBestCandidateScore = 0.0;

  for (unsigned int seed1Index=0; seed1Index<hits.size(); seed1Index++)
  {
    Hit& hSeed1   = hits[seed1Index];
    int nLaySeed1 = hSeed1.getLayer();

    if (nLaySeed1 == 2) continue; //layer 2 can't be the innermost seed stub
    if (nLaySeed1 > 3) break;     //no more possible combinations for this pattern

    //We have a correct Seed1


    for (unsigned int seed2Index = seed1Index+1; seed2Index<hits.size(); seed2Index++)
    {
      Hit& hSeed2   = hits[seed2Index];
      int nLaySeed2 = hSeed2.getLayer();

      if (nLaySeed1 == nLaySeed2) continue; //The seed layers have to be differents
      if (nLaySeed2 > 4) break;             //no more possible combinations for the current seed1


      //We have a correct Seed1/Seed2 combination !!!

      //Current candidate initialization (the 2 seeds)
      vecCurrentCandidateHits.clear();
      vecCurrentCandidateHits.push_back(&hSeed1);
      vecCurrentCandidateHits.push_back(&hSeed2);

      vecCurrentCandidateScore.clear();


      for (unsigned int testStubIndex = seed2Index+1; testStubIndex<hits.size(); testStubIndex++)
      {

        Hit& hTestStub    = hits[testStubIndex];
        int nLayTestStub  = hTestStub.getLayer();

        if (nLayTestStub == nLaySeed2) continue; //The layers have to be differents

        
        //Score processing of the Seed1/Seed2/testStub combination
        float tabScore[2];
        alignScore(hSeed1, hSeed2, hTestStub, tabScore, nFractionnalPartWidth);
        
        //cout<< "Score RPHI = "<<tabScore[0]<<"   Score RZ = "<<tabScore[1]<<endl;

        //Get the thresholds corresponding to the current layer combination
        float tabThresh[2];
        getThresholds(nLaySeed1, nLaySeed2, nLayTestStub, tabThresh);

        if (tabScore[0] <= tabThresh[0] && tabScore[1] <= tabThresh[1])
        {
          //The stub is in the window defined by the seed projection (correct stub candidate !)
          
          if (nLayTestStub != vecCurrentCandidateHits.back()->getLayer())
          {
            //The current testStub layer is not yet in the TC
            vecCurrentCandidateHits.push_back(&hTestStub);
            vecCurrentCandidateScore.push_back(tabScore[0]);
          }
          else if (tabScore[0] < vecCurrentCandidateScore.back())
          {
            //The layer is already in the TC but the Phi score of the current stub is better than the previous one
            vecCurrentCandidateHits.back()   = &hTestStub;            
            vecCurrentCandidateScore.back()  = tabScore[0];
          }
        }
      }

      //If the current candidate own more than 6 stubs, the lasts (outtermost) are removed
      while (vecCurrentCandidateHits.size() > 6)
      {
        vecCurrentCandidateHits.pop_back();
        vecCurrentCandidateScore.pop_back();
      }

      //All the stubs have been tested for the current Seeds combination

      if (int(vecCurrentCandidateHits.size()) >= nLayersCurrentPattern - nMissingHits)
      {
        //The current candidate has enough stubs to be a candidate
 
        //Process the score of the track candidate
        float fCurrentCandidateScore = 0.0;
        while (vecCurrentCandidateScore.empty() == false)
        {
          fCurrentCandidateScore += abs(vecCurrentCandidateScore.back());  //TODO norm?
          vecCurrentCandidateScore.pop_back();
        }

        if (vecCurrentCandidateHits.size() > vecBestCandidateHits.size() || (vecCurrentCandidateHits.size() == vecBestCandidateHits.size() && fCurrentCandidateScore < fBestCandidateScore))
        {
          //The current candidate is better than the previous best one
          vecBestCandidateHits = vecCurrentCandidateHits;
          fBestCandidateScore = fCurrentCandidateScore;
        }
      }
    }
  }

  //All the Seeds combinations have been tested

  if (vecBestCandidateHits.empty() == false)
  {
    //If there is a recorded best candidate

    //Fit the parameters and create the corresponding track object
    Track * fit_track;
    fit_track = createFittedTrack(vecBestCandidateHits);
    
    cout<<"adding one track..."<<endl;
    tracks.push_back(fit_track);
  }

}

void TCBuilder::fit(){
  for(unsigned int i=0;i<patterns.size();i++){
    vector<Hit*> allHits = patterns[i]->getHits();
    fit(allHits);
  }
}

TrackFitter* TCBuilder::clone(){
  TCBuilder* fit = new TCBuilder(nb_layers);
  fit->setPhiRotation(sec_phi);
  fit->setSectorID(sector_id);
  return fit;
}

