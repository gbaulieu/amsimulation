// Class for the data filtering
// For more info, look at the header file

#include "filter_stub.h"

// If you don't have the PR output in the same file than the rest


filter_stub::filter_stub(std::string filename, std::string secfilename, 
	       std::string outfile, int secid)
{  
  filter_stub::initTuple(filename,outfile);

  if (!filter_stub::convert(secfilename,secid)) return; // Don't go further if there is no sector file

  filter_stub::do_filter(secid); // Launch the filter loop

  m_outfile->Write();
  m_outfile->Close();      
}

/////////////////////////////////////////////////////////////////////////////////
//
// ==> filter_stub::do_filter(int secid,int hit_lim)
//
// Main method, where the filtering is made
//
/////////////////////////////////////////////////////////////////////////////////

void filter_stub::do_filter(int secid)
{
  const int m_nsec = m_sec_mult; // How many sectors are in the file

  int ndat = m_L1TT->GetEntries(); // How many events will we test

  cout << "Starting a filtering loop over " << ndat << " events..." << endl;
  cout << "... using " << m_nsec << " trigger sectors..." << endl;
  cout << "Looking events in sector " << secid << "..." << endl;

  int modid;
 
  float eta,phi;

  std::vector<int> list_sec;

  bool insec;

  // Loop over the events
 
  for (int i=0;i<ndat;++i)
  {    

    m_L1TT->GetEntry(i);
    
    if (i%100000==0) 
      cout << "Processed " << i << "/" << ndat << endl;

    filter_stub::reset();

    if (m_stub < 1) continue; // Not enough stubs anyway, don't go further
      
    for (int j=0;j<m_stub;++j)
    {  	
      modid = m_stub_detid[j]; 
      insec = false;
      
      for (unsigned int k=0;k<m_modules.size();++k)
      { 
	if (m_modules.at(k).at(0)!=modid) continue;
	
	for (unsigned int kk=1;kk<m_modules.at(k).size();++kk) // In which sector the module is
	  if (m_modules.at(k).at(kk)==secid) insec=true;;
	  
	if (insec) break;
	  
      }

      if (!insec) continue;
      
      mf_stub++;
      
      eta = -log(tan(atan2(sqrt(m_stub_x[j]*m_stub_x[j]+m_stub_y[j]*m_stub_y[j]),m_stub_z[j])/2));
      phi = atan2(m_stub_y[j],m_stub_x[j]);

      mf_stub_ptGEN->push_back(m_stub_ptGEN[j]);
      mf_stub_etaGEN->push_back(m_stub_etaGEN[j]);
      mf_stub_strip->push_back(m_stub_strip[j]);
      mf_stub_modid->push_back(m_stub_modid[j]);
      mf_stub_detid->push_back(m_stub_detid[j]);
      mf_stub_x->push_back(m_stub_x[j]); 
      mf_stub_y->push_back(m_stub_y[j]); 
      mf_stub_z->push_back(m_stub_z[j]); 
      mf_stub_eta->push_back(eta); 
      mf_stub_phi->push_back(phi);
      mf_stub_bend->push_back(m_stub_bend[j]); 
      mf_stub_X0->push_back(m_stub_X0[j]); 
      mf_stub_Y0->push_back(m_stub_Y0[j]); 
      mf_stub_Z0->push_back(m_stub_Z0[j]); 
      mf_stub_PHI0->push_back(m_stub_PHI0[j]);
      
    } // End of loop over stubs
 
    m_efftree->Fill(); // If yes fill the skimmed tree  
    
  } // End of loop on events; 
}


/////////////////////////////////////////////////////////////////////////////////
//
// ==> filter_stub::initTuple(std::string test,std::string out)
//
// This method opens and creates the differe rootuples involved
//
/////////////////////////////////////////////////////////////////////////////////


void filter_stub::initTuple(std::string test,std::string out)
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
  mf_stub_eta     = new  std::vector<float>;  
  mf_stub_phi     = new  std::vector<float>;
  mf_stub_bend    = new  std::vector<float>;    
  mf_stub_modid   = new  std::vector<int>;  
  mf_stub_detid   = new  std::vector<int>;
  mf_stub_X0      = new  std::vector<float>;  
  mf_stub_Y0      = new  std::vector<float>;  
  mf_stub_Z0      = new  std::vector<float>;
  mf_stub_PHI0    = new  std::vector<float>;

  filter_stub::reset();


  m_efftree->Branch("STUB_n",           &mf_stub);
  m_efftree->Branch("STUB_ptGEN",       &mf_stub_ptGEN);
  m_efftree->Branch("STUB_etaGEN",      &mf_stub_etaGEN);
  m_efftree->Branch("STUB_modid",       &mf_stub_modid);
  m_efftree->Branch("STUB_detid",       &mf_stub_detid);
  m_efftree->Branch("STUB_strip",       &mf_stub_strip);
  m_efftree->Branch("STUB_x",           &mf_stub_x);
  m_efftree->Branch("STUB_y",           &mf_stub_y);
  m_efftree->Branch("STUB_z",           &mf_stub_z);
  m_efftree->Branch("STUB_eta",         &mf_stub_eta);
  m_efftree->Branch("STUB_phi",         &mf_stub_phi);
  m_efftree->Branch("STUB_bend",        &mf_stub_bend);
  m_efftree->Branch("STUB_X0",          &mf_stub_X0);
  m_efftree->Branch("STUB_Y0",          &mf_stub_Y0);
  m_efftree->Branch("STUB_Z0",          &mf_stub_Z0);
  m_efftree->Branch("STUB_PHI0",        &mf_stub_PHI0);
}



/////////////////////////////////////////////////////////////////////////////////
//
// ==> filter_stub::convert(std::string sectorfilename) 
//
// Here we retrieve info from the TKLayout CSV file containing the sector definition
//
// This file contains, for each sector, the ids of the modules contained in the sector
//
// The role of this method is to create the opposite, ie a vector containing, for every module the list of sectors belonging to it
//
/////////////////////////////////////////////////////////////////////////////////

bool filter_stub::convert(std::string sectorfilename,int secid) 
{
  int modid;

  m_sec_mult = 0;

  std::vector<int> module;

  m_modules.clear();

  std::string STRING;
  std::ifstream in(sectorfilename.c_str());
  if (!in) 
  {
    std::cout << "Please provide a valid csv sector filename" << std::endl; 
    return false;
  }    
  
  int npar = 0;
  bool found;

  while (!in.eof()) 
  {
    ++m_sec_mult;
    getline(in,STRING);

    if (m_sec_mult!=secid+2) continue;

    std::istringstream ss(STRING);
    npar = 0;

    while (ss)
    {
      std::string s;
      if (!getline( ss, s, ',' )) break;

      ++npar;
      if (npar<=2) continue;

      modid = atoi(s.c_str());
      found=false;

      for (unsigned int i=0;i< m_modules.size();++i)
      {
	if (m_modules.at(i).at(0) == modid)
	{
	  found = true;
	  m_modules.at(i).push_back(m_sec_mult-2);
	  break;
	}
      }
      
      if (!found)
      {
	module.clear();
	module.push_back(modid);
	module.push_back(m_sec_mult-2);
	m_modules.push_back(module);
      }


      /*
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
      */
      //      m_modules.at(modid).push_back(m_sec_mult-2);
    }
  }

  std::cout << "Found " << m_modules.size() << " modules" << endl;

  in.close();

  m_sec_mult -= 3;

  return true;
}


void filter_stub::reset() 
{
  mf_stub = 0;

  mf_stub_etaGEN->clear();  
  mf_stub_strip->clear();   
  mf_stub_ptGEN->clear();   
  mf_stub_x->clear();       
  mf_stub_y->clear();       
  mf_stub_z->clear();     
  mf_stub_eta->clear();       
  mf_stub_phi->clear();    
  mf_stub_bend->clear();    
  mf_stub_modid->clear();   
  mf_stub_detid->clear();   
  mf_stub_X0->clear();      
  mf_stub_Y0->clear();      
  mf_stub_Z0->clear();      
  mf_stub_PHI0->clear();    

}
