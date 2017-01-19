#include <iostream>
#include <TChain.h>
#include <vector>

using namespace std;

class Event{

 public:
  TChain* TT;
  TChain* PATT;

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
  
  vector<int>           *p_m_stub_layer;
  vector<int>           *p_m_stub_module;
  vector<int>           *p_m_stub_ladder;
  vector<int>           *p_m_stub_seg;
  vector<float>         *p_m_stub_strip;
  vector<int>           *p_m_stub_tp;
  vector<float>         *p_m_stub_pxGEN;
  vector<float>         *p_m_stub_pyGEN;
  vector<float>         *p_m_stub_deltas;
  vector<float>         *p_m_stub_x0;
  vector<float>         *p_m_stub_y0;
  vector<float>         *p_m_stub_z0;
  vector<float>         *p_m_stub_phi0;
  vector<float>         *p_m_stub_etaGEN;
  vector<float>         *p_m_stub_x;
  vector<float>         *p_m_stub_y;
  vector<float>         *p_m_stub_z;

  int event_id;
  int m_patt;
  std::vector< std::vector<int> > m_patt_links;
  std::vector<int> m_patt_secid;
  std::vector<int> m_patt_id;
  std::vector<int> m_patt_miss;
  int nb_tc;
  std::vector<float> m_tc_pt;
  std::vector<float> m_tc_eta;
  std::vector<float> m_tc_phi;
  std::vector<float> m_tc_z;
  std::vector< std::vector<int> > m_tc_links;
  std::vector<int> m_tc_secid;
  std::vector<int> m_tc_pattid;
  int nb_tracks;
  std::vector<float> m_trk_pt;
  std::vector<float> m_trk_eta;
  std::vector<float> m_trk_chi2;
  std::vector<float> m_trk_phi;
  std::vector<float> m_trk_z;
  std::vector< std::vector<int> > m_trk_links;
  std::vector<int> m_trk_secid;
  std::vector<int> m_trk_pattid;

  std::vector< std::vector<int> > *p_m_patt_links;
  std::vector<int> *p_m_patt_secid;
  std::vector<int> *p_m_patt_id;
  std::vector<int> *p_m_patt_miss;
  std::vector<float> *p_m_tc_pt;
  std::vector<float> *p_m_tc_eta;
  std::vector<float> *p_m_tc_phi;
  std::vector<float> *p_m_tc_z;
  std::vector< std::vector<int> > *p_m_tc_links;
  std::vector<int> *p_m_tc_secid;
  std::vector<int> *p_m_tc_pattid;
  std::vector<float> *p_m_trk_pt;
  std::vector<float> *p_m_trk_eta;
  std::vector<float> *p_m_trk_chi2;
  std::vector<float> *p_m_trk_phi;
  std::vector<float> *p_m_trk_z;
  std::vector< std::vector<int> > *p_m_trk_links;
  std::vector<int> *p_m_trk_secid;
  std::vector<int> *p_m_trk_pattid;

  Event(string filename);
  ~Event();
  void getEntry(int i);
  int getEntries();

};
