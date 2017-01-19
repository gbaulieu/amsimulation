#include "Event.h"

Event::Event(string filename){
  TT = new TChain("TkStubs");
  TT->Add(filename.c_str());

  PATT = new TChain("L1tracks");
  PATT->Add(filename.c_str());

  p_m_stub_layer =  &m_stub_layer;
  p_m_stub_module = &m_stub_module;
  p_m_stub_ladder = &m_stub_ladder;
  p_m_stub_seg =    &m_stub_seg;
  p_m_stub_strip =  &m_stub_strip;
  p_m_stub_tp =     &m_stub_tp;
  p_m_stub_pxGEN = &m_stub_px_gen;  
  p_m_stub_pyGEN = &m_stub_py_gen;  
  p_m_stub_deltas = &m_stub_deltas;  
  p_m_stub_x0 =     &m_stub_x0;
  p_m_stub_y0 =     &m_stub_y0;
  p_m_stub_z0 =     &m_stub_z0;
  p_m_stub_phi0 =   &m_stub_phi0;
  p_m_stub_etaGEN = &m_stub_eta_gen;
  p_m_stub_x =      &m_stub_x;
  p_m_stub_y =      &m_stub_y;
  p_m_stub_z =      &m_stub_z;
  

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

  p_m_patt_links = &m_patt_links;
  p_m_patt_secid = &m_patt_secid;
  p_m_patt_id = &m_patt_id;
  p_m_patt_miss = &m_patt_miss;
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
  
  /////////////////////////////////////////
  
  // Branches definition
  
  PATT->SetBranchAddress("L1evt", &event_id); // Simple evt number or event ID
  PATT->SetBranchAddress("L1PATT_n",           &m_patt);
  PATT->SetBranchAddress("L1PATT_links",       &p_m_patt_links);
  PATT->SetBranchAddress("L1PATT_secid",       &p_m_patt_secid);
  PATT->SetBranchAddress("L1PATT_pattid",          &p_m_patt_id);
  PATT->SetBranchAddress("L1PATT_nmiss",       &p_m_patt_miss);
  PATT->SetBranchAddress("L1TC_n",            &nb_tc);
  PATT->SetBranchAddress("L1TC_links",        &p_m_tc_links);
  PATT->SetBranchAddress("L1TC_secid",        &p_m_tc_secid);
  PATT->SetBranchAddress("L1TC_pattid",       &p_m_tc_pattid);
  PATT->SetBranchAddress("L1TC_pt",           &p_m_tc_pt);
  PATT->SetBranchAddress("L1TC_phi",          &p_m_tc_phi);
  PATT->SetBranchAddress("L1TC_z",            &p_m_tc_z);
  PATT->SetBranchAddress("L1TC_eta",          &p_m_tc_eta);
  PATT->SetBranchAddress("L1TRK_n",            &nb_tracks);
  PATT->SetBranchAddress("L1TRK_links",        &p_m_trk_links);
  PATT->SetBranchAddress("L1TRK_secid",        &p_m_trk_secid);
  PATT->SetBranchAddress("L1TRK_pattid",       &p_m_trk_pattid);
  PATT->SetBranchAddress("L1TRK_pt",           &p_m_trk_pt);
  PATT->SetBranchAddress("L1TRK_phi",          &p_m_trk_phi);
  PATT->SetBranchAddress("L1TRK_z",            &p_m_trk_z);
  PATT->SetBranchAddress("L1TRK_eta",          &p_m_trk_eta);
  PATT->SetBranchAddress("L1TRK_chi2",         &p_m_trk_chi2);
}

Event::~Event(){
  delete TT;
  delete PATT;
}

int Event::getEntries(){
  return TT->GetEntries();
}

void Event::getEntry(int i){

  m_stub=0;
  n_evt=0;
  event_id=0;
  m_patt=0;
  nb_tc=0;
  nb_tracks=0;

  TT->GetEntry(i);
  PATT->GetEntry(i);
}
