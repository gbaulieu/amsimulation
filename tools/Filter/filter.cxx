// Class for the data filtering
// For more info, look at the header file

#include "filter.h"

// If you don't have the PR output in the same file than the rest


filter::filter(std::string filename, std::string secfilename, 
	       std::string outfile, int secid, int hit_lim)
{  
  filter::initTuple(filename,outfile);

  if (!filter::convert(secfilename)) return; // Don't go further if there is no sector file

  filter::do_filter(secid,hit_lim); // Launch the filter loop
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

  int ladder;
  int module;

  int n_hits[m_nsec];

  std::vector<int> list_sec;

  int new_nb_stub;

  // Loop over the events
 
  for (int i=0;i<ndat;++i)
  {    

    m_L1TT->GetEntry(i);
    
    if (i%100000==0) 
      cout << "Processed " << i << "/" << ndat << endl;

    // We have a muon and an anti-muon in each event
    // We treat each particle on after the other
    for(int pdg=-13;pdg<14;pdg+=26){
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
	  //Takes only stubs with the correct PDG ID
	  if(m_stub_pdg[j]!=pdg)
	    continue;

	  new_nb_stub++;
	  
	  mf_stub_ptGEN->push_back(m_stub_ptGEN[j]);
	  mf_stub_etaGEN->push_back(m_stub_etaGEN[j]);
	  mf_stub_strip->push_back(m_stub_strip[j]);
	  mf_stub_modid->push_back(m_stub_modid[j]);
	  
	  modid  = m_stub_modid[j]; 
	  layer  = int(modid/1000000); 
	  ladder  = int(modid%1000000/10000); 
	  module  = int(modid%10000/100); 
	  
	  if(layer<=7)
	    module=module/2;
	  
	  //modid  = int(modid/100); 
	  modid  = layer*10000+ladder*100+module; 
	  
	  if (m_modules.at(modid).size()>1)
	    {
	      for (unsigned int kk=1;kk<m_modules.at(modid).size();++kk) // In which sector the module is
		{
		  ++is_sec_there[m_modules.at(modid).at(kk)][layer-5]; 
		}
	    }
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

      if (!keepit) continue; // This track is better in another tower
      
      m_efftree->Fill(); // If yes fill the skimmed tree  
    }
  } // End of loop on tracks 

  m_outfile->Write();
  m_outfile->Close();
}


/////////////////////////////////////////////////////////////////////////////////
//
// ==> filter::initTuple(std::string test,std::string out)
//
// This method opens and creates the differe rootuples involved
//
/////////////////////////////////////////////////////////////////////////////////


void filter::initTuple(std::string test,std::string out)
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

  m_L1TT->SetBranchAddress("STUB_n",         &m_stub);
  m_L1TT->SetBranchAddress("STUB_modid",     &pm_stub_modid);
  m_L1TT->SetBranchAddress("STUB_strip",     &pm_stub_strip);
  m_L1TT->SetBranchAddress("STUB_ptGEN",     &pm_stub_ptGEN);
  m_L1TT->SetBranchAddress("STUB_etaGEN",    &pm_stub_etaGEN);
  m_L1TT->SetBranchAddress("STUB_pdg",       &pm_stub_pdg);

  // Output file definition (see the header)

  m_outfile = new TFile(out.c_str(),"recreate");

  m_efftree = new TTree("BankStubs","Stubs for bank");

  mf_stub_etaGEN = new std::vector<float>;
  mf_stub_strip  = new std::vector<float>;
  mf_stub_ptGEN  = new std::vector<float>;
  mf_stub_modid  = new std::vector<int>; 

  filter::reset();

  m_efftree->Branch("STUB_n",      &mf_stub);
  m_efftree->Branch("STUB_ptGEN",  &mf_stub_ptGEN);
  m_efftree->Branch("STUB_etaGEN", &mf_stub_etaGEN);
  m_efftree->Branch("STUB_modid",  &mf_stub_modid);
  m_efftree->Branch("STUB_strip",  &mf_stub_strip);
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
  int n_rods[6] = {16,24,34,48,62,76};

  int modid,lay,lad,mod;

  m_sec_mult = 0;

  std::vector<int> module;

  m_modules.clear();

  for (int i=0;i<230000;++i)
  {
    module.clear();
    module.push_back(-1);
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
      
      lay   = int(modid/10000); 
      modid-= 10000*lay;
      lad   = int(modid/100); 
      modid-= 100*lad;
      mod   = modid; 

      ///////
      // This hack is temporary and is due to a numbering problem in the TkLayout tool
      if (lay<=10) lad = (lad+n_rods[lay-5]/4)%(n_rods[lay-5]);
      if (lay<=7)  mod = mod/2;
      ///////

      modid = 10000*lay+100*lad+mod;

      m_modules.at(modid).push_back(m_sec_mult-2);
    }
  }

  in.close();

  m_sec_mult -= 1;

  return true;
}


void filter::reset() 
{
  mf_stub = 0;

  mf_stub_ptGEN->clear();  
  mf_stub_etaGEN->clear();  
  mf_stub_strip->clear();  
  mf_stub_modid->clear();  


}
