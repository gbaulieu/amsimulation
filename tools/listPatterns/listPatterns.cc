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

using namespace std;

class listPatterns{

public:
  void do_ana(std::string filename)
  {

    /***************** INPUT FILE ****************/
    TChain* TT = new TChain("TkStubs");
    TT->Add(filename.c_str());

    TChain* PATT = new TChain("L1tracks");
    PATT->Add(filename.c_str());
  
    int               n_evt;

    int m_stub;

    vector<int>           m_stub_layer;  // Layer du stub (5 a 10 pour les 6 layers qui nous interessent)
    vector<int>           m_stub_module; // Position en Z du module contenant le stub
    vector<int>           m_stub_ladder; // Position en PHI du module contenant le stub
    vector<int>           m_stub_seg;    // Segment du module contenant le stub
    vector<float>         m_stub_strip;  // Strip du cluster interne du stub
    vector<int>           m_stub_tp;     // particule du stub
    vector<float>         m_stub_px_gen; // pt initial de la particule ayant genere le stub
    vector<float>         m_stub_py_gen; // pt initial de la particule ayant genere le stub
    vector<float>         m_stub_deltas; // bend of the stub
    vector<float>         m_stub_x0;     // utilise pour calculer la distance au point d'interaction
    vector<float>         m_stub_y0;     // utilise pour calculer la distance au point d'interaction
    vector<float>         m_stub_z0;
    vector<float>         m_stub_phi0;
    vector<float>         m_stub_eta_gen;
    vector<float>         m_stub_x;      // x coordinate of the hit
    vector<float>         m_stub_y;      // y coordinate of the hit
    vector<float>         m_stub_z;      // z coordinate of the hit

    vector<int>           *p_m_stub_layer =  &m_stub_layer;
    vector<int>           *p_m_stub_module = &m_stub_module;
    vector<int>           *p_m_stub_ladder = &m_stub_ladder;
    vector<int>           *p_m_stub_seg =    &m_stub_seg;
    vector<float>         *p_m_stub_strip =  &m_stub_strip;
    vector<int>           *p_m_stub_tp =     &m_stub_tp;
    vector<float>         *p_m_stub_pxGEN = &m_stub_px_gen;  
    vector<float>         *p_m_stub_pyGEN = &m_stub_py_gen;  
    vector<float>         *p_m_stub_deltas = &m_stub_deltas;  
    vector<float>         *p_m_stub_x0 =     &m_stub_x0;
    vector<float>         *p_m_stub_y0 =     &m_stub_y0;
    vector<float>         *p_m_stub_z0 =     &m_stub_z0;
    vector<float>         *p_m_stub_phi0 =   &m_stub_phi0;
    vector<float>         *p_m_stub_etaGEN = &m_stub_eta_gen;
    vector<float>         *p_m_stub_x =      &m_stub_x;
    vector<float>         *p_m_stub_y =      &m_stub_y;
    vector<float>         *p_m_stub_z =      &m_stub_z;


    TT->SetBranchAddress("L1Tkevt",            &n_evt);
    TT->SetBranchAddress("L1TkSTUB_n",         &m_stub);
    TT->SetBranchAddress("L1TkSTUB_layer",     &p_m_stub_layer);
    TT->SetBranchAddress("L1TkSTUB_module",    &p_m_stub_module);
    TT->SetBranchAddress("L1TkSTUB_ladder",    &p_m_stub_ladder);
    TT->SetBranchAddress("L1TkSTUB_seg",       &p_m_stub_seg);
    TT->SetBranchAddress("L1TkSTUB_strip",     &p_m_stub_strip);
    TT->SetBranchAddress("L1TkSTUB_tp",        &p_m_stub_tp);
    TT->SetBranchAddress("L1TkSTUB_X0",        &p_m_stub_x0);
    TT->SetBranchAddress("L1TkSTUB_Y0",        &p_m_stub_y0);
    TT->SetBranchAddress("L1TkSTUB_Z0",        &p_m_stub_z0);
    TT->SetBranchAddress("L1TkSTUB_PHI0",      &p_m_stub_phi0);
    TT->SetBranchAddress("L1TkSTUB_etaGEN",    &p_m_stub_etaGEN);
    TT->SetBranchAddress("L1TkSTUB_pxGEN",     &p_m_stub_pxGEN);
    TT->SetBranchAddress("L1TkSTUB_pyGEN",     &p_m_stub_pyGEN);
    TT->SetBranchAddress("L1TkSTUB_deltas",    &p_m_stub_deltas);
    TT->SetBranchAddress("L1TkSTUB_x",         &p_m_stub_x);
    TT->SetBranchAddress("L1TkSTUB_y",         &p_m_stub_y);
    TT->SetBranchAddress("L1TkSTUB_z",         &p_m_stub_z);

    int event_id;
    int m_patt=0;
    std::vector< std::vector<int> > m_patt_links;
    std::vector<int> m_patt_secid;
    std::vector<int> m_patt_id;
    std::vector<int> m_patt_miss;
    
    int nb_tc=0;
    std::vector<float> m_tc_pt;
    std::vector<float> m_tc_eta;
    std::vector<float> m_tc_phi;
    std::vector<float> m_tc_z;
    std::vector< std::vector<int> > m_tc_links;
    std::vector<int> m_tc_secid;

    int nb_tracks=0;
    std::vector<float> m_trk_pt;
    std::vector<float> m_trk_eta;
    std::vector<float> m_trk_chi2;
    std::vector<float> m_trk_phi;
    std::vector<float> m_trk_z;
    std::vector< std::vector<int> > m_trk_links;
    std::vector<int> m_trk_secid;

    std::vector< std::vector<int> > *p_m_patt_links = &m_patt_links;
    std::vector<int> *p_m_patt_secid = &m_patt_secid;
    std::vector<int> *p_m_patt_id = &m_patt_id;
    std::vector<int> *p_m_patt_miss = &m_patt_miss;
    
    std::vector<float> *p_m_tc_pt = &m_tc_pt;
    std::vector<float> *p_m_tc_eta = &m_tc_eta;
    std::vector<float> *p_m_tc_phi = &m_tc_phi;
    std::vector<float> *p_m_tc_z = &m_tc_z;
    std::vector< std::vector<int> > *p_m_tc_links = &m_tc_links;
    std::vector<int> *p_m_tc_secid = &m_tc_secid;

    std::vector<float> *p_m_trk_pt = &m_trk_pt;
    std::vector<float> *p_m_trk_eta = &m_trk_eta;
    std::vector<float> *p_m_trk_chi2 = &m_trk_chi2;
    std::vector<float> *p_m_trk_phi = &m_trk_phi;
    std::vector<float> *p_m_trk_z = &m_trk_z;
    std::vector< std::vector<int> > *p_m_trk_links = &m_trk_links;
    std::vector<int> *p_m_trk_secid = &m_trk_secid;

    /////////////////////////////////////////
    
    // Branches definition
    
    PATT->SetBranchAddress("L1evt", &event_id); // Simple evt number or event ID
    
    PATT->SetBranchAddress("L1PATT_n",           &m_patt);
    PATT->SetBranchAddress("L1PATT_links",       &p_m_patt_links);
    PATT->SetBranchAddress("L1PATT_secid",       &p_m_patt_secid);
    PATT->SetBranchAddress("L1PATT_id",       &p_m_patt_id);
    PATT->SetBranchAddress("L1PATT_nmiss",       &p_m_patt_miss);
    
    PATT->SetBranchAddress("L1TC_n",            &nb_tc);
    PATT->SetBranchAddress("L1TC_links",        &p_m_tc_links);
    PATT->SetBranchAddress("L1TC_secid",        &p_m_tc_secid);
    PATT->SetBranchAddress("L1TC_pt",           &p_m_tc_pt);
    PATT->SetBranchAddress("L1TC_phi",          &p_m_tc_phi);
    PATT->SetBranchAddress("L1TC_z",            &p_m_tc_z);
    PATT->SetBranchAddress("L1TC_eta",          &p_m_tc_eta);

    PATT->SetBranchAddress("L1TRK_n",            &nb_tracks);
    PATT->SetBranchAddress("L1TRK_links",        &p_m_trk_links);
    PATT->SetBranchAddress("L1TRK_secid",        &p_m_trk_secid);
    PATT->SetBranchAddress("L1TRK_pt",           &p_m_trk_pt);
    PATT->SetBranchAddress("L1TRK_phi",          &p_m_trk_phi);
    PATT->SetBranchAddress("L1TRK_z",            &p_m_trk_z);
    PATT->SetBranchAddress("L1TRK_eta",          &p_m_trk_eta);
    PATT->SetBranchAddress("L1TRK_chi2",         &p_m_trk_chi2);

    int n_entries_MC = TT->GetEntries();
    cout<<n_entries_MC<<" events found"<<endl<<endl;

    // Loop on events
    for (int j=0;j<n_entries_MC;++j){
      TT->GetEntry(j); // Load entries
      PATT->GetEntry(j); // Load entries
      cout<<"event "<<n_evt<<" (index "<<j<<")"<<endl;
      if(m_patt>0){
	cout<<m_patt<<" pattern(s) found : "<<endl;
	cout<<endl;
	for(unsigned int k=0;k<m_patt_links.size();k++){
	  cout<<"Pattern "<<m_patt_id[k]<<endl;
	  for(unsigned int l=0;l<m_patt_links[k].size();l++){
	    int index = m_patt_links[k][l];
	    float pt = sqrt(m_stub_px_gen[index]*m_stub_px_gen[index]+m_stub_py_gen[index]*m_stub_py_gen[index]);
	    cout<<"Layer "<<m_stub_layer[index]<<" Ladder "<<m_stub_ladder[index]<<" Module "<<m_stub_module[index]<<" Seg "<<m_stub_seg[index]<<" Strip "<<int(m_stub_strip[index])<<" (TP="<<m_stub_tp[index]<<" PT="<<pt<<" GeV ETA="<<m_stub_eta_gen[index]<<" PHI="<<m_stub_phi0[index]<<")"<<endl;
	  }
	  cout<<endl;
	}
      }
      if(nb_tc>0){
	cout<<m_tc_links.size()<<" TC(s) found : "<<endl;
	cout<<endl;
	for(unsigned int k=0;k<m_tc_links.size();k++){
	  cout<<m_tc_links[k].size()<<" selected stubs :"<<endl;
	  for(unsigned int l=0;l<m_tc_links[k].size();l++){
	    int index = m_tc_links[k][l];
	    float pt = sqrt(m_stub_px_gen[index]*m_stub_px_gen[index]+m_stub_py_gen[index]*m_stub_py_gen[index]);
	    cout<<"Layer "<<m_stub_layer[index]<<" Ladder "<<m_stub_ladder[index]<<" Module "<<m_stub_module[index]<<" Seg "<<m_stub_seg[index]<<" Strip "<<int(m_stub_strip[index])<<" (TP="<<m_stub_tp[index]<<" PT="<<pt<<" GeV ETA="<<m_stub_eta_gen[index]<<" PHI="<<m_stub_phi0[index]<<")"<<endl;
	  }
	  cout<<"Parameters estimation : PT="<<m_tc_pt[k]<<" GeV - ETA="<<m_tc_eta[k]<<" - PHI="<<m_tc_phi[k]<<" - Z0="<<m_tc_z[k]<<endl;
	  cout<<endl;
	}
      }
      if(nb_tracks>0){
	cout<<m_trk_links.size()<<" tracks(s) found : "<<endl;
	cout<<endl;
	for(unsigned int k=0;k<m_trk_links.size();k++){
	  cout<<m_trk_links[k].size()<<" selected stubs :"<<endl;
	  for(unsigned int l=0;l<m_trk_links[k].size();l++){
	    int index = m_trk_links[k][l];
	    float pt = sqrt(m_stub_px_gen[index]*m_stub_px_gen[index]+m_stub_py_gen[index]*m_stub_py_gen[index]);
	    cout<<"Layer "<<m_stub_layer[index]<<" Ladder "<<m_stub_ladder[index]<<" Module "<<m_stub_module[index]<<" Seg "<<m_stub_seg[index]<<" Strip "<<int(m_stub_strip[index])<<" (TP="<<m_stub_tp[index]<<" PT="<<pt<<" GeV ETA="<<m_stub_eta_gen[index]<<" PHI="<<m_stub_phi0[index]<<")"<<endl;
	  }
	  cout<<"Parameters : CHI2="<<m_trk_chi2[k]<<" PT="<<m_trk_pt[k]<<" GeV - ETA="<<m_trk_eta[k]<<" - PHI="<<m_trk_phi[k]<<" - Z0="<<m_trk_z[k]<<endl;
	  cout<<endl;
	}
      }
      m_patt=0;
      nb_tc=0;
      nb_tracks=0;
    }

  delete TT;
  delete PATT;
  }

};

int main(int argv, char* args[]){
  if(argv>1){
    cout<<args[1]<<endl;
    string f(args[1]);
    listPatterns a;
    a.do_ana(f);
  }
  return 0;
}
