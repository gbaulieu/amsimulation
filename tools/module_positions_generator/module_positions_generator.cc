#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <TChain.h>
#include <TTree.h>
#include <map>
#include <memory>
#include "ModuleData.h"

using namespace std;

class module_positions_generator{

public:
  void do_ana(std::string filename,std::string output_filename)
  {

    /***************** INPUT FILE ****************/
    TChain* TT = new TChain("StripCoords");
    TT->Add(filename.c_str());

    ofstream output_file;
    output_file.open (output_filename);

    int nLayer;
    int nModule;
    int nLadder;
    int nStrip;
    int nSegment;
    float fX;
    float fY;
    float fZ;

    TT->SetBranchAddress("layer",   &nLayer);
    TT->SetBranchAddress("ladder",  &nLadder);
    TT->SetBranchAddress("module",  &nModule);
    TT->SetBranchAddress("strip",   &nStrip);
    TT->SetBranchAddress("segment", &nSegment);
    TT->SetBranchAddress("x",       &fX);
    TT->SetBranchAddress("y",       &fY);
    TT->SetBranchAddress("z",       &fZ);

    unsigned int n_entries_MC = TT->GetEntries();
    float factor = 100.0/n_entries_MC;
    
    int current_layer=-1;
    int current_ladder=-1;
    int current_module=-1;

    map<string,shared_ptr<ModuleData> > modules;
    map<string,shared_ptr<ModuleData> >::iterator it;

    cout<<"Generating modules coordinates..."<<endl;

    // Loop on events
    for (unsigned int j=0;j<n_entries_MC;++j){

      if(j%1000000==0){
	cout<<fixed;
	cout<<"\r"<<std::setprecision(1)<<j*factor<<" % "<<std::flush;
	cout<<std::setprecision(6);
	cout.unsetf(std::ios_base::floatfield);
      }

      TT->GetEntry(j); // Load entries
      if (current_layer!=nLayer || current_ladder!=nLadder || current_module!=nModule){//are we still on the same module?
	//The current module hase changed, check if we already know the new one
	//The key of the map is 2 digits for layer, 2 digits for ladder and 2 digits for module
	ostringstream key;
	key<<std::setfill('0');
	key<<setw(2)<<nLayer;
	key<<setw(2)<<nLadder;
	key<<setw(2)<<nModule;
	
	it = modules.find(key.str());
	if (it == modules.end()){//the module is not yet in the map, we add it
	  shared_ptr<ModuleData> mod = make_shared<ModuleData>(nLayer, nLadder, nModule);
	  modules[key.str()]=mod;
	  it = modules.find(key.str());
	}
      }
      
      it->second->update(nSegment, nStrip, fX, fY, fZ);

      current_layer=nLayer;
      current_ladder=nLadder;
      current_module=nModule;
    }

    output_file<<"Detector geometry infos"<<endl;
    output_file<<"Layer / Ladder / Module / Mid_Strip / Mid_Segment / Xrange / Yrange / Zrange / Xmid / Ymid / Zmid"<<endl;
    for(it=modules.begin();it!=modules.end();it++){
      output_file<<*it->second<<endl;
    }
    output_file.close();

    cout<<"\rData for "<<modules.size()<<" modules written in "<<output_filename<<endl;

    modules.clear();
    delete TT;
  }

};

int main(int argv, char* args[]){
  if(argv>2){
    string f(args[1]);
    string f2(args[2]);
    module_positions_generator a;
    a.do_ana(f,f2);
  }
  else{
    cout<<args[0]<<" <input ROOT file> <output TXT file>"<<endl;
  }
  return 0;
}
