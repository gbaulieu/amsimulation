// Class for the data filtering
// For more info, look at the header file

#include "filter.h"

// If you don't have the PR output in the same file than the rest


filter::filter(std::string filename, std::string secfilename, 
	       std::string outfile, int secid, int hit_lim, int format)
{  

  m_tilted=true;

  int n_tilted_rings[6];
  int n_flat_rings[6];

  for (int i=0; i < 6; ++i) n_tilted_rings[i]=0;
  for (int i=0; i < 6; ++i) n_flat_rings[i]=0;

  if (m_tilted)
  {
    n_tilted_rings[0]=11;
    n_tilted_rings[1]=12;
    n_tilted_rings[2]=12;
    n_flat_rings[0]=7;
    n_flat_rings[1]=11;
    n_flat_rings[2]=15;
  }

  for (int i=0; i < 6; ++i)
  {
    for (int j=0; j < 3; ++j)
    {
      limits[i][j]=0;

      if (n_tilted_rings[i]==0) continue;

      limits[i][j]=(j%2)*n_flat_rings[i]+(j>0)*n_tilted_rings[i];
    }
  }


  filter::initTuple(filename,outfile,format);

  if (!filter::convert(secfilename)) return; // Don't go further if there is no sector file

  filter::do_filter(secid,hit_lim); // Launch the filter loop

  m_outfile->Write();
  m_outfile->Close();      
}

/////////////////////////////////////////////////////////////////////////////////
//
// ==> filter::do_filter(int secid,int hit_lim)
//
// Main method, where the filtering is made
//
/////////////////////////////////////////////////////////////////////////////////

void filter::do_filter(int secid,int hit_lim)
{
  const int m_nsec = m_sec_mult; // How many sectors are in the file

  int ndat = m_L1TT->GetEntries(); // How many events will we test

  cout << "Starting a filtering loop over " << ndat << " events..." << endl;
  cout << "... using " << m_nsec << " trigger sectors..." << endl;
  cout << "Looking events in sector " << secid << " with nhits >= " << hit_lim << "..." << endl;

  int is_sec_there[m_nsec][20];
  int modid;
  int layer;

  int n_hits[m_nsec];

  std::vector<int> list_sec;

  int new_nb_stub;

  int nentr=0;

  bool insec;

  // Loop over the events
 
  for (int i=0;i<ndat;++i)
  {    

    m_L1TT->GetEntry(i);
    
    if (i%100000==0) 
      cout << "Processed " << i << "/" << ndat << endl;

    // We have a muon and an anti-muon in each event
    // We treat each particle on after the other
    for(int pdg=-13;pdg<14;pdg+=26)
    {
      new_nb_stub = 0;
      filter::reset();

      if (m_stub < 4) continue; // Not enough stubs anyway, don't go further

      for (int j=0;j<m_nsec;++j)
      {
	for (int k=0;k<20;++k) is_sec_there[j][k] = 0;
	n_hits[j] = 0;
      }
      
      list_sec.clear();
      
      for (int j=0;j<m_stub;++j)
      {  
	insec=false;
	
	//Takes only stubs with the correct PDG ID
	if(m_stub_pdg[j]!=pdg)
	  continue;
	
	modid  = int(m_stub_modid[j]/100); 
	layer  = int(m_stub_modid[j]/1000000); 
	
	//	cout << modid << " / " << m_modules.at(modid).size() << endl;

	if (m_modules.at(modid).size()==0) continue;

	for (unsigned int kk=0;kk<m_modules.at(modid).size();++kk) // In which sector the module is
	{
	  if (m_modules.at(modid).at(kk)==secid) insec=true;
	  //if (insec) cout << m_modules.at(modid).at(kk) << " / " << layer-5 << " / " << m_stub_detid[j] << " / " << insec << endl;
	  ++is_sec_there[m_modules.at(modid).at(kk)][layer-5]; 
	}

	if (!insec) continue; // Just write the stubs which are in the tower

	new_nb_stub++;
	mf_stub++;
	
	mf_stub_ptGEN->push_back(m_stub_ptGEN[j]);
	mf_stub_etaGEN->push_back(m_stub_etaGEN[j]);
	mf_stub_strip->push_back(m_stub_strip[j]);
	mf_stub_modid->push_back(m_stub_modid[j]);
	mf_stub_detid->push_back(m_stub_detid[j]);
	mf_stub_x->push_back(m_stub_x[j]); 
	mf_stub_y->push_back(m_stub_y[j]); 
	mf_stub_z->push_back(m_stub_z[j]); 
	mf_stub_bend->push_back(m_stub_bend[j]); 
	mf_stub_X0->push_back(m_stub_X0[j]); 
	mf_stub_Y0->push_back(m_stub_Y0[j]); 
	mf_stub_Z0->push_back(m_stub_Z0[j]); 
	mf_stub_PHI0->push_back(m_stub_PHI0[j]);
	mf_stub_pdg->push_back(m_stub_pdg[j]);

      } // End of loop over stubs

      if (new_nb_stub < 4) continue; // Not enough stubs anyway, don't go further

      mf_stub = new_nb_stub;

      // Check if the sector we are interested in contains the track
      // in a sufficiently large number of layers
      
      int n_hits_max=hit_lim;
      std::vector<int> sec_max;
      sec_max.clear();
      
      for (int j=0;j<m_nsec;++j)
      {
	n_hits[j]=0;
	  
	for (int k=0;k<20;++k)
	{
	  if (is_sec_there[j][k]>0) ++n_hits[j]; 
	}
	  
	if (n_hits[j]>=n_hits_max)
	{
	  if (n_hits[j]==n_hits_max) sec_max.push_back(j);
	  
	  if (n_hits[j]>n_hits_max)
	  {
	    n_hits_max=n_hits[j];
	    sec_max.clear();
	    sec_max.push_back(j);
	  }
	}
      } // End of loop on towers
      
      // sec_max contains the tower(s) containing most of the hits
      // this size is given by n_hits_max
      
      // If n_hits_max is equal to 5, we just have to be careful not to give the track to a barrel tower
      
      if (sec_max.size()==0) continue;

      bool keepit  = true;

      if(std::find(sec_max.begin(), sec_max.end(), secid) == sec_max.end()) //This track is not in our sector
	keepit=false;

      if(keepit)
      {
	if (secid<16) // -Z side Hybrid>Endcap
	{
	  for (int j=sec_max.size()-1;j>=0;--j)
	  {
	    if (secid==sec_max.at(j) && keepit) break;
	  
	    // There is a tower before the one (higher id) we are looking for which has 6 layers, barrel has priority
	    if (secid!=sec_max.at(j) && n_hits_max>=6) keepit=false;
	
	    // Max is 5 stubs and a tower with larger ID is before, we give the priority (hybrid goes first)
	    if (secid!=sec_max.at(j) && n_hits_max==5 && (sec_max.at(j)<16 || sec_max.at(j)>=32)) keepit=false;
	  }
	}
	else if (secid>=32) // +Z side Hybrid>Endcap
	{
	  for (unsigned int j=0;j<sec_max.size();++j)
	  {
	    if (secid==sec_max.at(j) && keepit) break;
	  
	    // There is a tower before the one we are looking for which has 6 layers, barrel has priority
	    if (secid!=sec_max.at(j) && n_hits_max>=6) keepit=false;
	
	    //  Max is 5 stubs and a tower with lower ID is before, we give the priority (hybrid goes first)
	    if (secid!=sec_max.at(j) && n_hits_max==5 && (sec_max.at(j)<16 || sec_max.at(j)>=32)) keepit=false;
	  }
	}
	else // Barrel
	{
	  for (unsigned int j=0;j<sec_max.size();++j)
	  {
	    if (secid==sec_max.at(j) && keepit) break;
	  
	    // There is a barrel tower before the one we are looking for which has 6 layers, barrel+ has priority
	    if (secid!=sec_max.at(j) && n_hits_max>=6  && (sec_max.at(j)>=16 && sec_max.at(j)<32)) keepit=false;
	  }
	}

      }

      //      cout << keepit << endl;

      if (!keepit) continue; // This track is better in another tower
      ++nentr;
      m_efftree->Fill(); // If yes fill the skimmed tree  
    }
  } // End of loop on tracks 

  cout << ">>> " << nentr << endl;
}


/////////////////////////////////////////////////////////////////////////////////
//
// ==> filter::initTuple(std::string test,std::string out)
//
// This method opens and creates the differe rootuples involved
//
/////////////////////////////////////////////////////////////////////////////////


void filter::initTuple(std::string test,std::string out, int format)
{
  m_L1TT   = new TChain("BankStubs"); 

  // Input data file 

  std::size_t found = test.find(".root");

  // Case 1, it's a root file
  if (found!=std::string::npos)
  {
    m_L1TT->Add(test.c_str());
  }
  else // This is a list provided into a text file
  {
    std::string STRING;
    std::ifstream in(test.c_str());
    if (!in) 
    {
      std::cout << "Please provide a valid data filename list" << std::endl; 
      return;
    }    
  
    while (!in.eof()) 
    {
      getline(in,STRING);

      found = STRING.find(".root");
      if (found!=std::string::npos) m_L1TT->Add(STRING.c_str());   
    }

    in.close();
  }

  pm_stub_modid=&m_stub_modid;
  pm_stub_strip=&m_stub_strip;
  pm_stub_ptGEN=&m_stub_ptGEN;
  pm_stub_etaGEN=&m_stub_etaGEN;
  pm_stub_pdg=&m_stub_pdg;
  pm_stub_detid=&m_stub_detid;  
  pm_stub_x=&m_stub_x;      
  pm_stub_y=&m_stub_y;     
  pm_stub_z=&m_stub_z;    
  pm_stub_bend=&m_stub_bend;  
  pm_stub_X0=&m_stub_X0;
  pm_stub_Y0=&m_stub_Y0;
  pm_stub_Z0=&m_stub_Z0;
  pm_stub_PHI0=&m_stub_PHI0;


  m_L1TT->SetBranchAddress("STUB_n",         &m_stub);
  m_L1TT->SetBranchAddress("STUB_modid",     &pm_stub_modid);
  m_L1TT->SetBranchAddress("STUB_detid",     &pm_stub_detid);
  m_L1TT->SetBranchAddress("STUB_strip",     &pm_stub_strip);
  m_L1TT->SetBranchAddress("STUB_ptGEN",     &pm_stub_ptGEN);
  m_L1TT->SetBranchAddress("STUB_etaGEN",    &pm_stub_etaGEN);
  m_L1TT->SetBranchAddress("STUB_pdg",       &pm_stub_pdg);
  m_L1TT->SetBranchAddress("STUB_x",           &pm_stub_x);
  m_L1TT->SetBranchAddress("STUB_y",           &pm_stub_y);
  m_L1TT->SetBranchAddress("STUB_z",           &pm_stub_z);
  m_L1TT->SetBranchAddress("STUB_bend",        &pm_stub_bend);
  m_L1TT->SetBranchAddress("STUB_X0",          &pm_stub_X0);
  m_L1TT->SetBranchAddress("STUB_Y0",          &pm_stub_Y0);
  m_L1TT->SetBranchAddress("STUB_Z0",          &pm_stub_Z0);
  m_L1TT->SetBranchAddress("STUB_PHI0",        &pm_stub_PHI0);

  // Output file definition (see the header)

  m_outfile = new TFile(out.c_str(),"recreate");

  m_efftree = new TTree("BankStubs","Stubs for bank");

  mf_stub_etaGEN  = new  std::vector<float>; 
  mf_stub_strip   = new  std::vector<float>; 
  mf_stub_ptGEN   = new  std::vector<float>;  
  mf_stub_x       = new  std::vector<float>;  
  mf_stub_y       = new  std::vector<float>;  
  mf_stub_z       = new  std::vector<float>;
  mf_stub_bend    = new  std::vector<float>;    
  mf_stub_modid   = new  std::vector<int>;  
  mf_stub_detid   = new  std::vector<int>;
  mf_stub_X0      = new  std::vector<float>;  
  mf_stub_Y0      = new  std::vector<float>;  
  mf_stub_Z0      = new  std::vector<float>;
  mf_stub_PHI0    = new  std::vector<float>;
  mf_stub_pdg  = new std::vector<int>; 

  filter::reset();

  m_efftree->Branch("STUB_n",           &mf_stub);
  m_efftree->Branch("STUB_ptGEN",       &mf_stub_ptGEN);
  m_efftree->Branch("STUB_etaGEN",      &mf_stub_etaGEN);
  m_efftree->Branch("STUB_modid",       &mf_stub_modid);
  m_efftree->Branch("STUB_detid",       &mf_stub_detid);
  m_efftree->Branch("STUB_strip",       &mf_stub_strip);
  m_efftree->Branch("STUB_pdg",         &mf_stub_pdg);
 
  if (format==1)
  {
    m_efftree->Branch("STUB_x",           &mf_stub_x);
    m_efftree->Branch("STUB_y",           &mf_stub_y);
    m_efftree->Branch("STUB_z",           &mf_stub_z);
    m_efftree->Branch("STUB_bend",        &mf_stub_bend);
    m_efftree->Branch("STUB_X0",          &mf_stub_X0);
    m_efftree->Branch("STUB_Y0",          &mf_stub_Y0);
    m_efftree->Branch("STUB_Z0",          &mf_stub_Z0);
    m_efftree->Branch("STUB_PHI0",        &mf_stub_PHI0);
    m_efftree->Branch("STUB_pdg",         &mf_stub_pdg);
  }
}



/////////////////////////////////////////////////////////////////////////////////
//
// ==> filter::convert(std::string sectorfilename) 
//
// Here we retrieve info from the TKLayout CSV file containing the sector definition
//
// This file contains, for each sector, the ids of the modules contained in the sector
//
// The role of this method is to create the opposite, ie a vector containing, for every module the list of sectors belonging to it
//
/////////////////////////////////////////////////////////////////////////////////

bool filter::convert(std::string sectorfilename) 
{
  int modid,lay,lad,mod,disk,type;

  //  std::cout << "Starting the conversion" << std::endl;

  m_sec_mult = 0;

  std::vector<int> module;

  m_modules.clear();

  for (unsigned int i=0;i<230000;++i)
  {
    module.clear();
    m_modules.push_back(module);
  }
 
  std::string STRING;
  std::ifstream in(sectorfilename.c_str());
  if (!in) 
  {
    std::cout << "Please provide a valid csv sector filename" << std::endl; 
    return false;
  }    
  
  int npar = 0;

  while (!in.eof()) 
  {
    ++m_sec_mult;
    getline(in,STRING);

    if (m_sec_mult<2) continue;

    std::istringstream ss(STRING);
    npar = 0;

    while (ss)
    {
      std::string s;
      if (!getline( ss, s, ',' )) break;

      ++npar;
      if (npar<=2) continue;

      modid = atoi(s.c_str());

      std::bitset<32> detid = modid; // Le detid

      int rmodid; // Le modid que l'on utilise

      if (detid[25]) // barrel;
      {
	lay  = 8*detid[23]+4*detid[22]+2*detid[21]+detid[20]+4;
	type = 2*detid[19]+detid[18];

	if (type==3) // Pas tilté
	{
	  lad  = 128*detid[17]+64*detid[16]+32*detid[15]+16*detid[14]+
	    8*detid[13]+4*detid[12]+2*detid[11]+detid[10]-1;
	  mod  = 128*detid[9]+64*detid[8]+32*detid[7]+16*detid[6]+
	    8*detid[5]+4*detid[4]+2*detid[3]+detid[2]-1+limits[lay-5][type-1];
	}
	else // tilté
	{
	  mod  = 128*detid[17]+64*detid[16]+32*detid[15]+16*detid[14]+
	    8*detid[13]+4*detid[12]+2*detid[11]+detid[10]-1+limits[lay-5][type-1];
	  lad  = 128*detid[9]+64*detid[8]+32*detid[7]+16*detid[6]+
	    8*detid[5]+4*detid[4]+2*detid[3]+detid[2]-1;
	}
      }
      else // endcap
      {
	disk  = 8*detid[21]+4*detid[20]+2*detid[19]+detid[18];
	lay   = 10+disk+abs(2-(2*detid[24]+detid[23]))*7;
	lad   = 32*detid[17]+16*detid[16]+8*detid[15]+4*detid[14]+2*detid[13]+detid[12]-1;

	if (disk>=3 && m_tilted) lad += 2;

	mod  = 128*detid[9]+64*detid[8]+32*detid[7]+16*detid[6]+
	  8*detid[5]+4*detid[4]+2*detid[3]+detid[2]-1;
      }


      rmodid = 10000*lay+100*lad+mod;

      //      if (m_sec_mult-2==0)
      //	std::cout << modid << " / " << rmodid << std::endl; 

      module = m_modules.at(rmodid);
      module.push_back(m_sec_mult-2);

      m_modules.at(rmodid) = module;
    }
  }

  //  std::cout << "Found " << m_modules.size() << " modules" << endl;

  in.close();

  m_sec_mult -= 2;

  return true;
}


void filter::reset() 
{
  mf_stub = 0;

  mf_stub_etaGEN->clear();  
  mf_stub_strip->clear();   
  mf_stub_ptGEN->clear();   
  mf_stub_x->clear();       
  mf_stub_y->clear();       
  mf_stub_z->clear();       
  mf_stub_bend->clear();    
  mf_stub_modid->clear();   
  mf_stub_detid->clear();   
  mf_stub_X0->clear();      
  mf_stub_Y0->clear();      
  mf_stub_Z0->clear();      
  mf_stub_PHI0->clear();    
  mf_stub_pdg->clear();  

}
