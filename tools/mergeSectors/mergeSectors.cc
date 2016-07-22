#include <iostream>
#include <sstream>
#include <cmath>
#include <TChain.h>
#include <TTree.h>
#include <TPolyMarker.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TAxis.h>
#include <TFile.h>
#include <TKey.h>
#include <TROOT.h>

using namespace std;

/**
   After AMPR analysis using the AMSimulation standalone package or the AMTestBench.py script, the ROOT output file contains one TTree per sector.
   This script allows to merge all the L1tracks_sec* TTree to one single L1tracks TTree.
 **/

class mergeSectors{

private:

  class SectorTree{
  private :

  public:
    int event_id;
    int m_patt;
    TTree* ttree;
    std::vector< std::vector<int> > m_patt_links;
    std::vector< std::vector<int> > m_patt_ss_size;
    std::vector<int> m_patt_secid;
    std::vector<int> m_patt_miss;
    std::vector<int> m_patt_id;
    std::vector<int> m_patt_nbstubs;
    std::vector<int> m_patt_area;
    
    int nb_tc;
    std::vector<float> m_tc_pt;
    std::vector<float> m_tc_eta;
    std::vector<float> m_tc_phi;
    std::vector<float> m_tc_z;
    std::vector< std::vector<int> > m_tc_links;
    std::vector<int> m_tc_secid;
    std::vector<int> m_tc_pattid;
    
    std::vector< std::vector<int> > *p_m_patt_links;
    std::vector< std::vector<int> > *p_m_patt_ss_size;
    std::vector<int> *p_m_patt_secid;
    std::vector<int> *p_m_patt_miss;
    std::vector<int> *p_m_patt_id;
    std::vector<int> *p_m_patt_nbstubs;
    std::vector<int> *p_m_patt_area;
    
    std::vector<float> *p_m_tc_pt;
    std::vector<float> *p_m_tc_eta;
    std::vector<float> *p_m_tc_phi;
    std::vector<float> *p_m_tc_z;
    std::vector< std::vector<int> > *p_m_tc_links;
    std::vector<int> *p_m_tc_secid; 
    std::vector<int> *p_m_tc_pattid; 

    int nb_tracks;
    std::vector<float> m_trk_pt;
    std::vector<float> m_trk_eta;
    std::vector<float> m_trk_chi2;
    std::vector<float> m_trk_phi;
    std::vector<float> m_trk_z;
    std::vector< std::vector<int> > m_trk_links;
    std::vector<int> m_trk_secid;
    std::vector<int> m_trk_pattid;
    
    std::vector<float> *p_m_trk_pt;
    std::vector<float> *p_m_trk_eta;
    std::vector<float> *p_m_trk_chi2;
    std::vector<float> *p_m_trk_phi;
    std::vector<float> *p_m_trk_z;
    std::vector< std::vector<int> > *p_m_trk_links;
    std::vector<int> *p_m_trk_secid; 
    std::vector<int> *p_m_trk_pattid; 

    SectorTree(TTree* t){
      ttree = t;
      m_patt = 0;
      nb_tc = 0;
      nb_tracks = 0;
      p_m_patt_links = &m_patt_links;
      p_m_patt_ss_size = &m_patt_ss_size;
      p_m_patt_secid = &m_patt_secid;
      p_m_patt_miss = &m_patt_miss;
      p_m_patt_id = &m_patt_id;
      p_m_patt_nbstubs = &m_patt_nbstubs;
      p_m_patt_area = &m_patt_area;
      p_m_tc_pt = &m_tc_pt;
      p_m_tc_eta = &m_tc_eta;
      p_m_tc_phi = &m_tc_phi;
      p_m_tc_z = &m_tc_z;
      p_m_tc_links = &m_tc_links;
      p_m_tc_secid = &m_tc_secid;
      p_m_tc_pattid = &m_tc_pattid;
      p_m_trk_pt = &m_trk_pt;
      p_m_trk_eta = &m_trk_eta;
      p_m_trk_chi2 = &m_trk_chi2;
      p_m_trk_phi = &m_trk_phi;
      p_m_trk_z = &m_trk_z;
      p_m_trk_links = &m_trk_links;
      p_m_trk_secid = &m_trk_secid;
      p_m_trk_pattid = &m_trk_pattid;

      t->SetBranchAddress("L1evt",              &event_id); // Simple evt number or event ID
      t->SetBranchAddress("L1PATT_n",           &m_patt);
      t->SetBranchAddress("L1PATT_links",       &p_m_patt_links);
      t->SetBranchAddress("L1PATT_ss_size",      &p_m_patt_ss_size);
      t->SetBranchAddress("L1PATT_secid",       &p_m_patt_secid);
      t->SetBranchAddress("L1PATT_nmiss",       &p_m_patt_miss);
      t->SetBranchAddress("L1PATT_pattid",          &p_m_patt_id);
      t->SetBranchAddress("L1PATT_nbstubs",          &p_m_patt_nbstubs);
      t->SetBranchAddress("L1PATT_area",          &p_m_patt_area);
      t->SetBranchAddress("L1TC_n",            &nb_tc);
      t->SetBranchAddress("L1TC_links",        &p_m_tc_links);
      t->SetBranchAddress("L1TC_secid",        &p_m_tc_secid);
      t->SetBranchAddress("L1TC_pattid",        &p_m_tc_pattid);
      t->SetBranchAddress("L1TC_pt",           &p_m_tc_pt);
      t->SetBranchAddress("L1TC_phi",          &p_m_tc_phi);
      t->SetBranchAddress("L1TC_z",            &p_m_tc_z);
      t->SetBranchAddress("L1TC_eta",          &p_m_tc_eta);
      t->SetBranchAddress("L1TRK_n",            &nb_tracks);
      t->SetBranchAddress("L1TRK_links",        &p_m_trk_links);
      t->SetBranchAddress("L1TRK_secid",        &p_m_trk_secid);
      t->SetBranchAddress("L1TRK_pattid",        &p_m_trk_pattid);
      t->SetBranchAddress("L1TRK_pt",           &p_m_trk_pt);
      t->SetBranchAddress("L1TRK_phi",          &p_m_trk_phi);
      t->SetBranchAddress("L1TRK_z",            &p_m_trk_z);
      t->SetBranchAddress("L1TRK_eta",          &p_m_trk_eta);
      t->SetBranchAddress("L1TRK_chi2",          &p_m_trk_chi2);
    }

    ~SectorTree(){
      delete ttree;
    }

    int getNbEntries(){
      return ttree->GetEntries();
    }

    void getEntry(int index){
      m_patt = 0;
      nb_tc = 0;
      ttree->GetEntry(index);
    }

  };
  
  
public:

  mergeSectors(){
    //used to generate the root dictionnary to support vector<vector<int>> and vector<float>
    gROOT->ProcessLine(".L Loader.C+");
  }

  void do_ana(std::string filename)
  {

    TFile* f = new TFile(filename.c_str(),"update");  
    TIter next(f->GetListOfKeys());
    TKey *key;
    vector<SectorTree*> sector_trees;
    int nb_records = -1;
    while ((key=(TKey*)next())) {
      string name(key->GetName());
      if(name.find("L1tracks_sec")!=std::string::npos){
	SectorTree* s = new SectorTree((TTree*)(f->Get(key->GetName())));
	int nb = s->getNbEntries();
	if(nb_records==-1)
	  nb_records = nb;
	if(nb_records!=nb){
	  cout<<"All sector TTree must have the same number of records! ("<<name<<" differs)"<<endl;
	  return;
	}
	sector_trees.push_back(s);
      }
    }

    f->Delete("L1tracks;*");  

    TTree* PATT = new TTree("L1tracks","Official L1-AM tracks info");
  
    int event_id;
    int m_patt=0;
    std::vector< std::vector<int> > m_patt_links;
    std::vector< std::vector<int> > m_patt_ss_size;
    std::vector<int> m_patt_secid;
    std::vector<int> m_patt_miss;
    std::vector<int> m_patt_id;
    std::vector<int> m_patt_nbstubs;
    std::vector<int> m_patt_area;
    
    int nb_tc=0;
    std::vector<float> m_tc_pt;
    std::vector<float> m_tc_eta;
    std::vector<float> m_tc_phi;
    std::vector<float> m_tc_z;
    std::vector< std::vector<int> > m_tc_links;
    std::vector<int> m_tc_secid;
    std::vector<int> m_tc_pattid;

    int nb_tracks=0;
    std::vector<float> m_trk_pt;
    std::vector<float> m_trk_eta;
    std::vector<float> m_trk_chi2;
    std::vector<float> m_trk_phi;
    std::vector<float> m_trk_z;
    std::vector< std::vector<int> > m_trk_links;
    std::vector<int> m_trk_secid;
    std::vector<int> m_trk_pattid;

    std::vector< std::vector<int> > *p_m_patt_links = &m_patt_links;
    std::vector< std::vector<int> > *p_m_patt_ss_size = &m_patt_ss_size;
    std::vector<int> *p_m_patt_secid = &m_patt_secid;
    std::vector<int> *p_m_patt_miss = &m_patt_miss;
    std::vector<int> *p_m_patt_id = &m_patt_id;
    std::vector<int> *p_m_patt_nbstubs = &m_patt_nbstubs;
    std::vector<int> *p_m_patt_area = &m_patt_area;
    
    std::vector<float> *p_m_tc_pt = &m_tc_pt;
    std::vector<float> *p_m_tc_eta = &m_tc_eta;
    std::vector<float> *p_m_tc_phi = &m_tc_phi;
    std::vector<float> *p_m_tc_z = &m_tc_z;
    std::vector< std::vector<int> > *p_m_tc_links = &m_tc_links;
    std::vector<int> *p_m_tc_secid = &m_tc_secid;
    std::vector<int> *p_m_tc_pattid = &m_tc_pattid;

    std::vector<float> *p_m_trk_pt = &m_trk_pt;
    std::vector<float> *p_m_trk_eta = &m_trk_eta;
    std::vector<float> *p_m_trk_chi2 = &m_trk_chi2;
    std::vector<float> *p_m_trk_phi = &m_trk_phi;
    std::vector<float> *p_m_trk_z = &m_trk_z;
    std::vector< std::vector<int> > *p_m_trk_links = &m_trk_links;
    std::vector<int> *p_m_trk_secid = &m_trk_secid;
    std::vector<int> *p_m_trk_pattid = &m_trk_pattid;

    /////////////////////////////////////////
    
    // Branches definition
    
    PATT->Branch("L1evt", &event_id); // Simple evt number or event ID
    
    PATT->Branch("L1PATT_n",           &m_patt);
    PATT->Branch("L1PATT_links",       &p_m_patt_links);
    PATT->Branch("L1PATT_ss_size",     &p_m_patt_ss_size);
    PATT->Branch("L1PATT_secid",       &p_m_patt_secid);
    PATT->Branch("L1PATT_nmiss",       &p_m_patt_miss);
    PATT->Branch("L1PATT_pattid",       &p_m_patt_id);
    PATT->Branch("L1PATT_nbstubs",       &p_m_patt_nbstubs);
    PATT->Branch("L1PATT_area",       &p_m_patt_area);
    
    PATT->Branch("L1TC_n",            &nb_tc);
    PATT->Branch("L1TC_links",        &p_m_tc_links);
    PATT->Branch("L1TC_secid",        &p_m_tc_secid);
    PATT->Branch("L1TC_pattid",        &p_m_tc_pattid);
    PATT->Branch("L1TC_pt",           &p_m_tc_pt);
    PATT->Branch("L1TC_phi",          &p_m_tc_phi);
    PATT->Branch("L1TC_z",            &p_m_tc_z);
    PATT->Branch("L1TC_eta",          &p_m_tc_eta);

    PATT->Branch("L1TRK_n",            &nb_tracks);
    PATT->Branch("L1TRK_links",        &p_m_trk_links);
    PATT->Branch("L1TRK_secid",        &p_m_trk_secid);
    PATT->Branch("L1TRK_pattid",        &p_m_trk_pattid);
    PATT->Branch("L1TRK_pt",           &p_m_trk_pt);
    PATT->Branch("L1TRK_phi",          &p_m_trk_phi);
    PATT->Branch("L1TRK_z",            &p_m_trk_z);
    PATT->Branch("L1TRK_eta",          &p_m_trk_eta);
    PATT->Branch("L1TRK_chi2",          &p_m_trk_chi2);

    for(int i=0;i<nb_records;i++){
      for(unsigned int j=0;j<sector_trees.size();j++){
	SectorTree* st = sector_trees[j];
	st->getEntry(i);
	event_id = st->event_id;
	m_patt = m_patt + st->m_patt;
	nb_tc = nb_tc + st->nb_tc;
	nb_tracks = nb_tracks + st->nb_tracks;
	
	for(int k=0;k<st->m_patt;k++){
	  vector<int> vec;
	  vector<int> vec_size;
	  for(unsigned int l=0;l<st->m_patt_links[k].size();l++){
	    vec.push_back(st->m_patt_links[k][l]);
	  }
	  for(unsigned int l=0;l<st->m_patt_ss_size[k].size();l++){
	    vec_size.push_back(st->m_patt_ss_size[k][l]);
	  }
	  m_patt_links.push_back(vec);
	  m_patt_ss_size.push_back(vec_size);
	  m_patt_secid.push_back(st->m_patt_secid[k]);
	  m_patt_miss.push_back(st->m_patt_miss[k]);
	  m_patt_id.push_back(st->m_patt_id[k]);
	  m_patt_nbstubs.push_back(st->m_patt_nbstubs[k]);
	  m_patt_area.push_back(st->m_patt_area[k]);
	}
	for(int k=0;k<st->nb_tc;k++){
	  vector<int> vec;
	  for(unsigned int l=0;l<st->m_tc_links[k].size();l++){
	    vec.push_back(st->m_tc_links[k][l]);
	  }
	  m_tc_links.push_back(vec);
	}
	for(int k=0;k<st->nb_tc;k++){
	  m_tc_secid.push_back(st->m_tc_secid[k]);
	}
	for(int k=0;k<st->nb_tc;k++){
	  m_tc_pattid.push_back(st->m_tc_pattid[k]);
	}
	for(int k=0;k<st->nb_tc;k++){
	  m_tc_pt.push_back(st->m_tc_pt[k]);
	  m_tc_phi.push_back(st->m_tc_phi[k]);
	  m_tc_z.push_back(st->m_tc_z[k]);
	  m_tc_eta.push_back(st->m_tc_eta[k]);
	}
	for(int k=0;k<st->nb_tracks;k++){
	  vector<int> vec;
	  for(unsigned int l=0;l<st->m_trk_links[k].size();l++){
	    vec.push_back(st->m_trk_links[k][l]);
	  }
	  m_trk_links.push_back(vec);
	  m_trk_secid.push_back(st->m_trk_secid[k]);
	  m_trk_pattid.push_back(st->m_trk_pattid[k]);
	  m_trk_pt.push_back(st->m_trk_pt[k]);
	  m_trk_phi.push_back(st->m_trk_phi[k]);
	  m_trk_z.push_back(st->m_trk_z[k]);
	  m_trk_eta.push_back(st->m_trk_eta[k]);
	  m_trk_chi2.push_back(st->m_trk_chi2[k]);
	}
      }
      
      PATT->Fill();

      m_patt = 0;
      nb_tc = 0;
      nb_tracks = 0;
      m_patt_links.clear();
      m_patt_ss_size.clear();
      m_patt_secid.clear();
      m_patt_miss.clear();
      m_patt_id.clear();
      m_patt_nbstubs.clear();
      m_patt_area.clear();
      m_tc_links.clear();
      m_tc_secid.clear();
      m_tc_pattid.clear();
      m_tc_pt.clear();
      m_tc_phi.clear();
      m_tc_z.clear();
      m_tc_eta.clear();
      m_trk_links.clear();
      m_trk_secid.clear();
      m_trk_pattid.clear();
      m_trk_pt.clear();
      m_trk_phi.clear();
      m_trk_z.clear();
      m_trk_eta.clear();
      m_trk_chi2.clear();
    }
    
    PATT->Write();

    for(unsigned int i=0;i<sector_trees.size();i++){
      SectorTree* st = sector_trees[i];
      delete st;
    }
    sector_trees.clear();
    
    delete PATT;
    f->Close();
    delete f;
    
  }

  void removeSectorTrees(std::string filename)
  {
    TFile* f = new TFile(filename.c_str(),"update");  
    TIter next(f->GetListOfKeys());
    TKey *key;
    while ((key=(TKey*)next())) {
      string name(key->GetName());
      if(name.find("L1tracks_sec")!=std::string::npos){
	name = name+";*";
	cout<<"removing "<<name<<endl;
	f->Delete(name.c_str());  
      }
    }
    f->Write();
    f->Close();
  }

};

int main(int argv, char* args[]){
  if(argv>1){
    string f(args[1]);
    mergeSectors a;
    if(argv>2 && strcmp(args[2],"--clean")==0){
      a.removeSectorTrees(f);      
    }
    else{
      a.do_ana(f);
    }
  }
  else{
    cout<<"mergeSectors <ROOT File> : Merge all L1tracks_sec* TTrees into one single L1tracks TTree."<<endl;
    cout<<"mergeSectors <ROOT File> --clean : Remove all L1tracks_sec* TTrees from the file (NO MERGING DONE)."<<endl;
  }
  return 0;
}
