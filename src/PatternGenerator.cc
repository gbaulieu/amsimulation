#include "PatternGenerator.h"

PatternGenerator::PatternGenerator(){
  variableRes = 0;
  ptMin=2;
  ptMax=100;
  etaMin=0.0f;
  etaMax=1.0f;
}

void PatternGenerator::setMinPT(float minp){
  if(minp>0)
    ptMin = minp;
}

void PatternGenerator::setMaxPT(float maxp){
  if(maxp>0)
    ptMax = maxp;
}

void PatternGenerator::setMinEta(float mine){
  etaMin = mine;
}

void PatternGenerator::setMaxEta(float maxe){
  etaMax = maxe;
}

void PatternGenerator::setMaxFakeSuperstrips(int mf){
  if(mf>-1)
    nbMaxFakeSuperstrips = mf;
  else
    nbMaxFakeSuperstrips = 0;
}

void PatternGenerator::setLayers(vector<int> l){
  tracker_layers = l;
  sort(tracker_layers.begin(),tracker_layers.end());
}

void PatternGenerator::setInactiveLayers(vector<int> l){
  inactive_layers = l;
}

void PatternGenerator::setParticuleDirName(string f){
  particuleDirName = f;
}

void PatternGenerator::setVariableResolution(int nb){
  if(nb<0 || nb>3)
    nb = 0;
  variableRes = nb;
}

int PatternGenerator::getVariableResolutionState(){
  return variableRes;
}

TChain* PatternGenerator::createTChain(string directoryName, string tchainName){

  cout<<"Loading files from "<<directoryName<<" (this may take several minutes)..."<<endl;

  TChain* TT = new TChain("BankStubs");

  TSystemDirectory dir(directoryName.c_str(), directoryName.c_str());
  TList *files = dir.GetListOfFiles();

  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);

    while((file=(TSystemFile*)next())){
      fname = file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(".root")) {
	TT->Add((directoryName + fname.Data()).c_str());
      }
    }
  }
  else{//TSystemDirectory not supported for xrootd
    /*
     If using xrootd, directoryName should be a file (not a directory) containing one input file name per line 
    */
    ifstream in(directoryName.c_str());
    string fname;
    
    while (std::getline(in,fname)){
      TT->Add(fname.c_str());      
    }
    
    in.close();
  }

  p_m_stub_modid = &m_stub_modid; 
  p_m_stub_strip = &m_stub_strip;
  p_m_stub_ptGEN = &m_stub_ptGEN;  
  p_m_stub_etaGEN = &m_stub_etaGEN;  
  
  TT->SetBranchAddress("STUB_n",         &m_stub);
  TT->SetBranchAddress("STUB_modid",     &p_m_stub_modid);
  TT->SetBranchAddress("STUB_strip",     &p_m_stub_strip);
  TT->SetBranchAddress("STUB_ptGEN",     &p_m_stub_ptGEN);
  TT->SetBranchAddress("STUB_etaGEN",    &p_m_stub_etaGEN);
  TT->SetBranchStatus("*",0);
  TT->SetBranchStatus("STUB_n",1);
  TT->SetBranchStatus("STUB_modid",1);
  TT->SetBranchStatus("STUB_strip",1); 
  TT->SetBranchStatus("STUB_ptGEN",1); 
  TT->SetBranchStatus("STUB_etaGEN",1);

  int nb_entries = TT->GetEntries();
  cout<<nb_entries<<" events found."<<endl;

  return TT;
}

int PatternGenerator::generate(TChain* TT, int* evtIndex, int evtNumber, int* nbTrackUsed, SectorTree* sectors, map<int,pair<float,float> > eta_limits, int* coverageEstimation){

  if(tracker_layers.size()==0){
    cout<<"No layer defined!"<<endl;
    return -1;
  }
  vector<Pattern*> patterns;

  //--> Signification (et dimension) des variables

  int n_entries_TT = TT->GetEntries();

  int nbInLayer=0;
  int nbInSector = 0;
  int nbModuleOk = 0;

  int ld_fd_factor = (int)pow(2.0,(double)variableRes);

  int layers[tracker_layers.size()];
  /*
  int layers_backup[tracker_layers.size()];
  TH1I* sstrips_error = new TH1I("PHI errors","Missmatch on phi (in superstrips)", 201, -100, 100);
  TH1I* sstrips_error2 = new TH1I("Z Errors","Missmatch on Z (in outer segments)", 101, -50, 50);
  */
  vector<int> ladder_per_layer(tracker_layers.size());
  vector<int> module_per_layer(tracker_layers.size());

  while(nbModuleOk<evtNumber && (*evtIndex)<n_entries_TT){
    TT->GetEntry((*evtIndex));
    (*evtIndex)++;

    //cout<<"index "<<*evtIndex<<endl;

    //initialize arrays
    for(unsigned int j=0;j<tracker_layers.size();j++){
      layers[j]=-1;
    }
    /*
    for(unsigned int j=0;j<tracker_layers.size();j++){
      layers_backup[j]=-1;
    }
    */
    for(unsigned int j=0;j<tracker_layers.size();j++){
      ladder_per_layer[j]=-1;
    }
    for(unsigned int j=0;j<tracker_layers.size();j++){
      module_per_layer[j]=-1;
    }

    float current_eta = -10;

    //check the layers of the stubs
    for(int j=0;j<m_stub;j++){
    
      if(m_stub_etaGEN[j]<etaMin){// eta of the generating particule is bellow the threshold -> we do not use it for pattern generation
	continue;
      }
      if(m_stub_etaGEN[j]>etaMax){// eta of the generating particule is above the threshold -> we do not use it for pattern generation
	continue;
      }
      if(m_stub_ptGEN[j]<ptMin){// The PT of the generating particule is below the minimum required -> we do not use it for pattern generation
        continue;
      }
      if(m_stub_ptGEN[j]>ptMax){// The PT of the generating particule is above the maximum accepted -> we do not use it for pattern generation
        continue;
      }

      int value = m_stub_modid[j];
      //cout<<value<<endl;
      int layer = value/1000000;
      value = value-layer*1000000;
      int ladder = value/10000;
      value = value-ladder*10000;
      int module = value/100;
      value = value-module*100;
      //cout<<"layer : "<<layer<<" ladder : "<<ladder<<" module : "<<module<<" segment : "<<value<<endl;

      vector<int>::iterator iter;
      iter=find(inactive_layers.begin(),inactive_layers.end(),layer);
      if(iter!=inactive_layers.end()){
	continue;
      }

      int layer_position=-1;
      for(unsigned int cpt=0;cpt<tracker_layers.size();cpt++){
	if(layer==tracker_layers[cpt]){
	  layer_position=(int)cpt;
	  break;
	}
      }

      if(layer_position!=-1){ // is this layer in the layer list?
	/*
	if(layers[layer_position]!=-1)
	  layers_backup[layer_position]=layers[layer_position];
	*/
	layers[layer_position]=j;
	ladder_per_layer[layer_position]=CMSPatternLayer::getLadderCode(layer, ladder);
	module = CMSPatternLayer::getModuleCode(layer, module);
	module_per_layer[layer_position]=module;
      }

      current_eta = m_stub_etaGEN[j];

    }


    /**************************************
    Selection on the stubs/layer
    We need at least one stub per layer
    **************************************/
    bool missing_stub = false;
    int nbFakeSuperstrip = 0;
    for(unsigned int j=0;j<tracker_layers.size();j++){
      if(layers[j]==-1){
	missing_stub=true;
	if(sectors->getNbSectors()==1 && current_eta!=-10){ // we can use fake superstrips if we know the sector in which to add the tracks and if we have at least one stub (current_eta!=10)
	  if(eta_limits.find(tracker_layers[j])!=eta_limits.end()){//we have eta boundaries for this layer
	    pair<float,float> limits = eta_limits[tracker_layers[j]];
	    if(current_eta<limits.first || current_eta>limits.second){ // we are outside the eta limits for this layer -> we will add a fake superstrip for this layer
	      if(nbFakeSuperstrip<nbMaxFakeSuperstrips){//we don't want to have more than nbMaxFakeSuperstrips fake superstrips in the pattern
		//cout<<"missing hit on layer "<<tracker_layers[j]<<" for track with eta="<<current_eta<<endl;
		layers[j]=-2;
		//we put a ladder and a module just to be inside the sector
		ladder_per_layer[j]=sectors->getAllSectors()[0]->getLadders(j)[0];
		module_per_layer[j]=sectors->getAllSectors()[0]->getModules(j,ladder_per_layer[j])[0];
		//cout<<"Add stub for sector : "<<ladder_per_layer[j]<<" / "<<module_per_layer[j]<<endl;
		//debug=true;
		missing_stub=false;//we will create a fake superstrip, so the stub is not missing
		nbFakeSuperstrip++;
	      }
	    }
	  }
	}
      }
      if(missing_stub)
	break;
    }

    if(missing_stub){
      /*
      cout<<"stubs manquants ";
      for(unsigned int j=0;j<tracker_layers.size();j++){
	cout<<layers[j]<<",";
      }
      cout<<endl;
      */
      continue;//no stub on each layer -> drop the event    
    }

    nbInLayer++;

    //    cout<<"trace ok"<<endl;

    /****************************************
    Check that the track is part of a sector
    ****************************************/

    Sector* sector = sectors->getSector(ladder_per_layer, module_per_layer);
    
    /*
    for(unsigned int j=0;j<tracker_layers.size();j++){
      cout<<"stubs "<<ladder_per_layer[j]<<"/"<<module_per_layer[j]<<"-";
    }
    cout<<endl;
    cout<<"Secteur : "<<sector<<endl;
    */

    if(sector==NULL){
      //cout<<"No sector found"<<endl;
      continue;
    }
    else{
      nbInSector++;
    }

    float last_pt = 0;
    Pattern* p = new Pattern(tracker_layers.size());
    Pattern* lowDef_p=NULL;
    
    if(variableRes){ // we use variable resolution patterns so we create 2 patterns with different resolution
      lowDef_p = new Pattern(tracker_layers.size());
    }
    for(unsigned int j=0;j<tracker_layers.size();j++){
      int stub_number = layers[j];
      short module = -1;
      short ladder = -1;
      short strip = -1;
      short seg = -1;
      /*
      int stub_number_backup = layers_backup[j];
      short module_backup = -1;
      short ladder_backup = -1;
      short strip_backup = -1;
      short seg_backup = -1;
      */
      if(stub_number==-2){//creation of a fake superstrip
	module=0;
	strip=0;
	seg=0;
	ladder=1;//should never happen so this superstrip will never be activated
      }
      else{
	int value = m_stub_modid[stub_number];
	//	cout<<value<<endl;
	value = value-(value/1000000)*1000000;
	ladder = value/10000;
	value = value-ladder*10000;
	module = value/100;
	value = value-module*100;
	seg = value;
	/*
	if(stub_number_backup!=-1){
	  int value_backup = m_stub_modid[stub_number_backup];
	  //	cout<<value<<endl;
	  value_backup = value_backup-(value_backup/1000000)*1000000;
	  ladder_backup = value_backup/10000;
	  value_backup = value_backup-ladder_backup*10000;
	  module_backup = value_backup/100;
	  value_backup = value_backup-module_backup*100;
	  seg_backup = value_backup;
	  strip_backup = m_stub_strip[stub_number_backup];
	}
	*/
	seg = CMSPatternLayer::getSegmentCode(tracker_layers[j], CMSPatternLayer::getLadderCode(tracker_layers[j],ladder), seg);	
	//convertion to sector relative positions
	//module = sector->getModuleCode(tracker_layers[j], CMSPatternLayer::getLadderCode(tracker_layers[j],ladder), CMSPatternLayer::getModuleCode(tracker_layers[j],module));
	//ladder=sector->getLadderCode(tracker_layers[j],CMSPatternLayer::getLadderCode(tracker_layers[j],ladder));
	
	strip = m_stub_strip[stub_number];

	//if(variableRes){
	//  stripLD = strip/ld_fd_factor;
	//}
      }
      
      
      //cout<<"Event "<<*evtIndex<<endl;
      //cout<<"    Layer "<<tracker_layers[j]<<" segment "<<seg<<" module "<<module<<" ladder "<<ladder<<" strip "<<strip<<endl;
      //cout<<endl;
      

      if(stub_number!=-2)//this is not a fake stub
	last_pt = m_stub_ptGEN[stub_number];
      CMSPatternLayer pat;
      CMSPatternLayer lowDef_layer;
      try{
	pat.computeSuperstrip(tracker_layers[j], module, ladder, strip, seg, sectors->getSuperStripSize(tracker_layers[j]));
      }
      catch (const std::out_of_range& oor) {
	std::cerr << "One of the stubs can not be linked to a superstrip"<<endl;
	continue;
      }
      /*
      if(stub_number_backup!=-1 && ladder!=ladder_backup && module==module_backup){
	try{
	  CMSPatternLayer pat_backup;
	  pat_backup.computeSuperstrip(tracker_layers[j], module_backup, ladder_backup, strip_backup, seg_backup, sectors->getSuperStripSize(tracker_layers[j]));
	  cout<<"we have an overlap :"<<endl;
	  cout<<"    Layer "<<tracker_layers[j]<<" segment "<<seg<<" module "<<module<<" ladder "<<ladder<<" strip "<<strip<<endl;
	  cout<<pat.toString()<<endl;
	  cout<<"    Layer "<<tracker_layers[j]<<" segment "<<seg_backup<<" module "<<module_backup<<" ladder "<<ladder_backup<<" strip "<<strip_backup<<endl;
	  cout<<pat_backup.toString()<<endl;
	  int diff_phi = pat.getStrip()-pat_backup.getStrip();
	  int diff_eta = pat.getModule()-pat_backup.getModule();
	  cout<<"diff PHI is "<<diff_phi<<endl;
	  cout<<"diff ETA is "<<diff_eta<<endl;
	  cout<<endl;
	  sstrips_error->Fill(diff_phi);
	  sstrips_error2->Fill(diff_eta);
	}
	catch (const std::out_of_range& oor) {

	}
      }
      */
      p->setLayerStrip(j, &pat);

      if(variableRes){
	lowDef_layer.computeSuperstrip(tracker_layers[j], module, ladder, strip, seg, sectors->getSuperStripSize(tracker_layers[j])*ld_fd_factor);
	lowDef_p->setLayerStrip(j, &lowDef_layer);
      }

    }
    if(p==NULL){
      continue;
    }
    nbModuleOk++;

    //cout<<"creation pattern : "<<endl;
    //cout<<*lowDef_p<<endl;
    
    if(coverageEstimation==NULL){
      if(variableRes){
	sector->getPatternTree()->addPattern(lowDef_p,p, last_pt);
      }
      else{
	sector->getPatternTree()->addPattern(p,NULL, last_pt);
      }
    }
    else{
      if(variableRes){
	if(sector->getPatternTree()->checkPattern(lowDef_p, p))//does the bank contains the pattern?
	  (*coverageEstimation)++;
      }
    }
    
    delete p;
    if(lowDef_p!=NULL)
        delete lowDef_p;

  }

  *nbTrackUsed=nbModuleOk;
  /*
  sstrips_error->SetFillColor(41);
  sstrips_error->Write();
  delete sstrips_error;
  sstrips_error2->SetFillColor(41);
  sstrips_error2->Write();
  delete sstrips_error2;
  */
  if(coverageEstimation==NULL){
    cout<<"Event index : "<<*evtIndex<<endl;
    cout<<"Nb Events with stubs on all layers : "<<nbInLayer<<endl;
    cout<<"Nb Events in sectors : "<<nbInSector<<endl;

    if(variableRes)
      return sectors->getFDPatternNumber();
    else
      return sectors->getLDPatternNumber();
  }
  else
    return 0;

}



void PatternGenerator::generate(SectorTree* sectors, int step, float threshold, map<int,pair<float,float> > eta_limits){
  int nbPatterns = 1;
  int newCount = 0;
  int indexPart = 0;
  float dif=1;
 
  int loop=0;

  int nbTracks = 0;
  int trackUsed = 0;

  float tracks[10000];
  float patterns[10000];
  int iterationNbTracks=1;

  //creates the TChain from the directory name
  TChain* tc = createTChain(particuleDirName, particuleDirName);

  cout<<"Starting patterns generation using iterations of "<<step<<" events"<<endl;
  
  while((1-dif)<threshold && iterationNbTracks>0){
    threshold=0;
    loop++;
    iterationNbTracks=0;
    newCount = generate(tc, &indexPart, step, &nbTracks, sectors, eta_limits);
    trackUsed+=nbTracks;
    iterationNbTracks+=nbTracks;
    dif=(newCount-nbPatterns)/(float)iterationNbTracks;//% of coverage for this iteration
    nbPatterns=newCount;
    if(iterationNbTracks!=0){
      tracks[loop-1]=nbPatterns;
      patterns[loop-1]=(1-dif)*100;
    }
    if(loop>1){
      if(iterationNbTracks==0)
	cout<<"No more tracks to use : final patterns bank size : "<<nbPatterns<<endl;
      else
	cout<<"Current patterns bank size : "<<nbPatterns<<" (Coverage : "<<(1-dif)*100<<"%)"<<endl;
    }
  }
  
  TGraph* nbPatt = new TGraph(loop-1,tracks,patterns);
  nbPatt->SetTitle("Pattern Bank Generation");
  nbPatt->GetXaxis()->SetTitle("Patterns bank size");
  nbPatt->GetYaxis()->SetTitle("Coverage (%)");
  cout<<trackUsed<<" tracks used to generate patterns"<<endl;
  cout<<"Last event index : "<<indexPart<<endl;
  nbPatt->Write();
  delete nbPatt;
  
  if(variableRes){
    cout<<"Creating variable resolution bank..."<<endl;
    sectors->computeAdaptativePatterns(variableRes);
    if(iterationNbTracks==step){
      cout<<"Estimating coverage..."<<endl;
      int recognizedTracks=0;
      generate(tc, &indexPart, step, &nbTracks, sectors, eta_limits,&recognizedTracks);
      float cov = recognizedTracks*100/nbTracks;
      cout<<"Estimated coverage with variable resolution : "<<cov<<"% (estimated on "<<nbTracks<<" tracks)"<<endl;
    }
    else{
      cout<<"Not enough tracks to estimate the coverage."<<endl;
    }
  }

  //free the TChain
  delete tc;

  vector<Sector*> v_sector = sectors->getAllSectors();
  for(unsigned int k=0;k<v_sector.size();k++){
    vector<int> PT = v_sector[k]->getPatternTree()->getPTHisto();
    TH1I* pt_histo = new TH1I("PT sector "+k,"PT of pattern generating tracks", 110, 0, 110);
    for(int i=0;i<101;i++){
      for(int j=0;j<PT[i];j++){
	pt_histo->Fill(i);
      }
    }
    pt_histo->SetFillColor(41);
    pt_histo->Write();
    delete pt_histo;
  }
}
