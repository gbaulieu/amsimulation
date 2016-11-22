#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <TChain.h>
#include <TTree.h>
#include <TPolyMarker.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TAxis.h>
#include <TFile.h>

#include "Event.h"

using namespace std;

class comparePatterns{

public:
  int do_ana(std::string filename, std::string filename2)
  {

    Event file1(filename);
    Event file2(filename2);
    int n_entries_MC = file1.getEntries();
    int n_entries_MC2 = file2.getEntries();

    if(n_entries_MC != n_entries_MC2){
      cout<<"Both files do not have the same number of events!"<<endl;
      return -1;
    }

    cout<<n_entries_MC<<" events found"<<endl<<endl;

    unsigned int nb_patterns=0;
    unsigned int added_patterns=0;
    unsigned int missing_patterns=0;
    unsigned int common_patterns=0;

    cout<<"Comparing lists of patterns..."<<endl;
    for (int j=0;j<n_entries_MC;++j){ 
      file1.getEntry(j); // Load entries
      file2.getEntry(j); // Load entries

      map<int,int> file1_pattern_index;
      map<int,int> file2_pattern_index;
      for(unsigned int k=0;k<file1.m_patt_id.size();k++){
	file1_pattern_index[file1.m_patt_id[k]]=k;
      }
      for(unsigned int k=0;k<file2.m_patt_id.size();k++){
	file2_pattern_index[file2.m_patt_id[k]]=k;
      }

      sort(file1.m_patt_id.begin(),file1.m_patt_id.end());
      sort(file2.m_patt_id.begin(),file2.m_patt_id.end());

      unsigned int idx1 = 0;
      unsigned int idx2 = 0;

      bool first = true;

      nb_patterns+=file1.m_patt_id.size();

      // We check that the 2 lists of patterns IDs are the same
      while(idx1<file1.m_patt_id.size() && idx2<file2.m_patt_id.size()){
	if(file1.m_patt_id[idx1]<file2.m_patt_id[idx2]){
	  // The value in list 1 is smaller so it's missing in list 2
	  if(first){
	    cout<<endl;
	    cout<<"Event index "<<j<<endl;
	    first=false;
	  }
	  missing_patterns++;
	  cout<<file1.m_patt_id[idx1]<<" is missing!"<<endl;
	  idx1++;
	}
	else if(file1.m_patt_id[idx1]>file2.m_patt_id[idx2]){
	  // The value in list 2 is smaller so it's missing in list 1
	  if(first){
	    cout<<endl;
	    cout<<"Event index "<<j<<endl;
	    first=false;
	  }
	  added_patterns++;
	  cout<<file2.m_patt_id[idx2]<<" is added!"<<endl;
	  idx2++;
	}
	else{
	  // The pattern is in both lists (good)
	  // we now check the stubs
	  int links_index_1 = file1_pattern_index[file1.m_patt_id[idx1]];
	  int links_index_2 = file2_pattern_index[file2.m_patt_id[idx2]];
	  sort(file1.m_patt_links[links_index_1].begin(),file1.m_patt_links[links_index_1].end());
	  sort(file2.m_patt_links[links_index_2].begin(),file2.m_patt_links[links_index_2].end());

	  bool error = false;
	  if(file1.m_patt_links[links_index_1].size()==file2.m_patt_links[links_index_2].size()){
	    for(unsigned int k=0;k<file1.m_patt_links[links_index_1].size();k++){
	      if(file1.m_patt_links[links_index_1][k]!=file2.m_patt_links[links_index_2][k]){
		error=true;
		break;
	      }
	    }
	  }
	  else{
	    error=true;
	  }
	  if(error){
	    cout<<"Event index "<<j<<endl;
	    cout<<"Pattern "<<file1.m_patt_id[links_index_1]<<" : we do not have the same lists of stubs!"<<endl;
	    for(unsigned int k=0;k<file1.m_patt_links[links_index_1].size();k++){
	      cout<<file1.m_patt_links[links_index_1][k]<<" ";
	    }
	    cout<<endl;
	    for(unsigned int k=0;k<file2.m_patt_links[links_index_2].size();k++){
	      cout<<file2.m_patt_links[links_index_2][k]<<" ";
	    }
	    cout<<endl;
	    cout<<"layer : ";
	    for(unsigned int k=0;k<file1.m_patt_links[links_index_1].size();k++){
	      cout<<file1.m_stub_layer[file1.m_patt_links[links_index_1][k]]<<" ";
	    }
	    cout<<"\t->\t";
	    for(unsigned int k=0;k<file2.m_patt_links[links_index_2].size();k++){
	      cout<<file2.m_stub_layer[file2.m_patt_links[links_index_2][k]]<<" ";
	    }
	    cout<<endl;
	  }
	  idx1++;
	  idx2++;
	  common_patterns++;
	}
      }
      if(idx1==file1.m_patt_id.size()){
	// If there are remaining values in list 1 they are missing in list 2
	while(idx2<file2.m_patt_id.size()){
	  if(first){
	    cout<<endl;
	    cout<<"Event index "<<j<<endl;
	    first=false;
	  }	  
	  added_patterns++;
	  cout<<file2.m_patt_id[idx2]<<" is added!"<<endl;
	  idx2++;
	}
      }
      if(idx2==file2.m_patt_id.size()){
	// If there are remaining values in list 2 they are missing in list 1
	while(idx1<file1.m_patt_id.size()){
	  if(first){
	    cout<<endl;
	    cout<<"Event index "<<j<<endl;
	    first=false;
	  }
	  missing_patterns++;
	  cout<<file1.m_patt_id[idx1]<<" is missing!"<<endl;
	  idx1++;
	}
      }
    }

    unsigned int nb_tc=0;
    unsigned int added_tc=0;
    unsigned int missing_tc=0;
    unsigned int common_tc=0;
    cout<<endl;
    cout<<"Comparing lists of TCs..."<<endl;
    cout<<endl;
    for (int j=0;j<n_entries_MC;++j){ 
      file1.getEntry(j); // Load entries
      file2.getEntry(j); // Load entries
      
      map<int,int> file1_tc_index;
      map<int,int> file2_tc_index;
      for(unsigned int k=0;k<file1.m_tc_pattid.size();k++){
	file1_tc_index[file1.m_tc_pattid[k]]=k;
      }
      for(unsigned int k=0;k<file2.m_tc_pattid.size();k++){
	file2_tc_index[file2.m_tc_pattid[k]]=k;
      }

      sort(file1.m_tc_pattid.begin(),file1.m_tc_pattid.end());
      sort(file2.m_tc_pattid.begin(),file2.m_tc_pattid.end());

      unsigned int idx1 = 0;
      unsigned int idx2 = 0;

      bool first = true;

      nb_tc += file1.m_tc_pattid.size();

      // We check that the 2 lists of patterns leading to TCs are the same
      while(idx1<file1.m_tc_pattid.size() && idx2<file2.m_tc_pattid.size()){
	if(file1.m_tc_pattid[idx1]<file2.m_tc_pattid[idx2]){
	  // The value in list 1 is smaller so it's missing in list 2
	  if(first){
	    cout<<endl;
	    cout<<"Event index "<<j<<endl;
	    first=false;
	  }
	  missing_tc++;
	  cout<<file1.m_tc_pattid[idx1]<<" is missing!"<<endl;
	  idx1++;
	}
	else if(file1.m_tc_pattid[idx1]>file2.m_tc_pattid[idx2]){
	  // The value in list 2 is smaller so it's missing in list 1
	  if(first){
	    cout<<endl;
	    cout<<"Event index "<<j<<endl;
	    first=false;
	  }
	  added_tc++;
	  cout<<file2.m_tc_pattid[idx2]<<" is added!"<<endl;
	  idx2++;
	}
	else{
	  // The TC is in both lists (good)
	  // we now check the stubs
	  int links_index_1 = file1_tc_index[file1.m_tc_pattid[idx1]];
	  int links_index_2 = file2_tc_index[file2.m_tc_pattid[idx2]];
	  sort(file1.m_tc_links[links_index_1].begin(),file1.m_tc_links[links_index_1].end());
	  sort(file2.m_tc_links[links_index_2].begin(),file2.m_tc_links[links_index_2].end());

	  bool error = false;
	  if(file1.m_tc_links[links_index_1].size()==file2.m_tc_links[links_index_2].size()){
	    for(unsigned int k=0;k<file1.m_tc_links[links_index_1].size();k++){
	      if(file1.m_tc_links[links_index_1][k]!=file2.m_tc_links[links_index_2][k]){
		error=true;
		break;
	      }
	    }
	  }
	  else{
	    error=true;
	  }
	  if(error){
	    cout<<"Event index "<<j<<endl;
	    cout<<"Pattern "<<file1.m_tc_pattid[idx1]<<" : we do not have the same lists of stubs!"<<endl;
	    for(unsigned int k=0;k<file1.m_tc_links[links_index_1].size();k++){
	      cout<<file1.m_tc_links[links_index_1][k]<<" ";
	    }
	    cout<<endl;
	    for(unsigned int k=0;k<file2.m_tc_links[links_index_2].size();k++){
	      cout<<file2.m_tc_links[links_index_2][k]<<" ";
	    }
	    cout<<endl;
	    cout<<"layer : ";
	    for(unsigned int k=0;k<file1.m_tc_links[links_index_1].size();k++){
	      cout<<file1.m_stub_layer[file1.m_tc_links[links_index_1][k]]<<" ";
	    }
	    cout<<"\t->\t";
	    for(unsigned int k=0;k<file2.m_tc_links[links_index_2].size();k++){
	      cout<<file2.m_stub_layer[file2.m_tc_links[links_index_2][k]]<<" ";
	    }
	    cout<<endl;
	  }
	  common_tc++;
	  idx1++;
	  idx2++;
	}
      }
      if(idx1==file1.m_tc_pattid.size()){
	// If there are remaining values in list 1 they are missing in list 2
	while(idx2<file2.m_tc_pattid.size()){
	  if(first){
	    cout<<endl;
	    cout<<"Event index "<<j<<endl;
	    first=false;
	  }
	  added_tc++;
	  cout<<file2.m_tc_pattid[idx2]<<" is added!"<<endl;
	  idx2++;
	}
      }
      if(idx2==file2.m_tc_pattid.size()){
	// If there are remaining values in list 2 they are missing in list 1
	while(idx1<file1.m_tc_pattid.size()){
	  if(first){
	    cout<<endl;
	    cout<<"Event index "<<j<<endl;
	    first=false;
	  }
	  missing_tc++;
	  cout<<file1.m_tc_pattid[idx1]<<" is missing!"<<endl;
	  idx1++;
	}
      }

    }
    
    cout<<endl;
    cout<<"***** PATTERNS *****"<<endl;
    cout<<"Common patterns : "<<(common_patterns*100)/(float)nb_patterns<<"% of first file"<<endl;
    cout<<"Missing patterns : "<<(missing_patterns*100)/(float)nb_patterns<<"% of first file"<<endl;
    cout<<"Additional patterns : "<<(added_patterns*100)/(float)nb_patterns<<"% of first file"<<endl;
    cout<<endl;
    cout<<"***** TRACK CANDIDATES *****"<<endl;
    cout<<"Common TCs : "<<(common_tc*100)/(float)nb_tc<<"% of first file"<<endl;
    cout<<"Missing TCs : "<<(missing_tc*100)/(float)nb_tc<<"% of first file"<<endl;
    cout<<"Additional TCs : "<<(added_tc*100)/(float)nb_tc<<"% of first file"<<endl;
    cout<<endl;

    if(common_patterns==nb_patterns && common_tc==nb_tc && missing_patterns==0 && added_patterns==0 && missing_tc==0 && added_tc==0)
      return 0;
    else
      return -1;
  }

};

int main(int argv, char* args[]){
  int res = -1;
  if(argv>2){
    cout<<args[1]<<endl;
    cout<<args[2]<<endl;
    string f(args[1]);
    string f2(args[2]);
    comparePatterns a;
    res = a.do_ana(f,f2);
  }
  return res;
}
