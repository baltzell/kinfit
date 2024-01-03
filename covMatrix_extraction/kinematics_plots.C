#include <vector>
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "reader.h"
#include <dirent.h>
#include <sys/types.h>
#include <cstring>
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TF1.h"
#include "TStyle.h"
#include "TMath.h"
#include <iostream>
#include "TPaveStats.h"
#include "TLatex.h"
#include <TFile.h>
#include "TTree.h"
#include "TCut.h"
#include <TMatrixDSymEigen.h>
#include "TMatrixDSym.h"
#include <numeric>
#include <random>
#include <unordered_map>
#include <TLine.h>
#include "TRobustEstimator.h"
#include "TError.h"



// Struct for storing info saved from Hipo files
struct partInfoStruct
{
  std::vector<float> mc_p_vec, mc_theta_deg_vec, mc_phi_deg_vec;
  std::vector<float> rec_p_vec, rec_theta_deg_vec, rec_phi_deg_vec;
  std::vector<float> chi2_vec, ndf_vec; 
};



void read_Hipo(char in_data[256], int part_charge, float part_mass, partInfoStruct &partInfo);
void read_Part_Bank(hipo::bank PartBank, std::vector<int>* pid, std::vector<float>* px, std::vector<float>* py, std::vector<float>* pz);
void read_Part_Bank(hipo::bank PartBank, hipo::bank MCMatchBank, std::vector<int>* pid, std::vector<float>* px, std::vector<float>* py, std::vector<float>* pz, std::vector<int> *rec_status, std::vector<float> *chi2, std::vector<int> *ndf);
void get_track_parameters(int index, const hipo::bank TrackBank, int& track_sec, float& track_chi2, int& track_ndf);


// Utility function to generate reverse index map
std::map<int, std::vector<int>> loadMapByIndex(const hipo::bank fromBank, const char* idxVarName)
{
  std::map<int, std::vector<int>> map;
  if (fromBank.getRows() > 0) {
    for (int iFrom = 0; iFrom < fromBank.getRows(); ++iFrom)
      {
	const int iTo = fromBank.getInt(idxVarName, iFrom);
	map[iTo].push_back(iFrom);
      }
  }
  return map;
}

// Select particle type
const int particle_id = 211;
//const int particle_id = -211


//Initialize the particle name
//TString part_name;
//if (particle_id == 211) {
//  part_name = "pip";
//}
//else if (particle_id == -211) {
//  part_name = "pim";
// }

//These variables are used for the plots that combine all P, theta, and phi bins
const std::vector<TString> KINES = {"P", "#theta", "#phi"};
const std::vector<TString> UNITS = {"[GeV]", "[deg]", "[deg]"};

//Use these values to set a maximum on the number of files that the reader loop will read and to set a number of elements for the vectors to be initilaized to
const int max_files = 10;
const int N_events_per_file = 10000;
const int N_events = max_files*N_events_per_file;
//const TString out_name = Form("%s_test_%i_in_files", part_name, max_files);

int covMatrix_extraction()
{
  // Input Data
  //std::vector<std::string> input_directories = {"/volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/pim_sec2_test/cooked/"};
  std::vector<std::string> input_directories = {"/volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/pip_sec2_truth/cooked/"};

  partInfoStruct partInfo;
  int event_count = 0;

  // Initialize particle charge and mass
  //if ((strcmp(part_type, "#pi^{+}") == 0) || (strcmp(part_type, "#pi^{-}") == 0)) {
  //  part_mass = 0.13957;
  //  if (strcmp(part_type, "#pi^{+}") == 0) {
  //    part_charge = 1;
  //  }
  //  else if (strcmp(part_type, "#pi^{-}") == 0) {
  //    part_charge = -1;
  //  }
  //}

  TString part_name;
  TString part_type;
  float part_mass;
  int part_charge;
  if (particle_id == 211) {
    part_mass = 0.13957;
    part_type = "#pi^{+}";
    part_name = "pip";
    part_charge = 1;
  }
  else if (particle_id == -211) {
    part_mass = 0.13957;
    part_type = "#pi^{-}";
    part_name = "pim";
    part_charge = -1;
  }

  TString out_name = Form("%s_test_%i_in_files", part_name.Data(), max_files);

  partInfo.mc_p_vec.reserve(N_events);
  partInfo.rec_p_vec.reserve(N_events);
  partInfo.mc_theta_deg_vec.reserve(N_events);
  partInfo.rec_theta_deg_vec.reserve(N_events);
  partInfo.mc_phi_deg_vec.reserve(N_events);
  partInfo.rec_phi_deg_vec.reserve(N_events);
  partInfo.chi2_vec.reserve(N_events);
  partInfo.ndf_vec.reserve(N_events);

  TFile* kinInfo_file = new TFile(Form("kinInfo_root_files/kinInfo_%s.root", out_name.Data()), "RECREATE");

  TTree* kinInfo_tree = new TTree("kinInfoTree", "Kinematics Tree");
  float mc_p, mc_theta, mc_phi, rec_p, rec_theta, rec_phi, chi2, ndf;
  kinInfo_tree->Branch("mc_p", &partInfo.mc_p_vec);
  kinInfo_tree->Branch("mc_theta", &partInfo.mc_theta_deg_vec);
  kinInfo_tree->Branch("mc_phi", &partInfo.mc_phi_deg_vec);
  kinInfo_tree->Branch("rec_p", &partInfo.rec_p_vec);
  kinInfo_tree->Branch("rec_theta", &partInfo.rec_theta_deg_vec);
  kinInfo_tree->Branch("rec_phi", &partInfo.rec_phi_deg_vec);
  kinInfo_tree->Branch("chi2", &partInfo.chi2_vec);
  kinInfo_tree->Branch("ndf", &partInfo.ndf_vec);

  int in_file_count = 0;
  for (const std::string& in_data : input_directories) {
    DIR *dr;
    struct dirent *en;
    dr = opendir(in_data.c_str());
    //If in_data is a directory, loop through all files within it  
    if (dr) {
      while (((en = readdir(dr)) != NULL) && (in_file_count < max_files)) {
	if ((strcmp(en->d_name, ".") != 0) && (strcmp(en->d_name, "..") != 0)) {
	  in_file_count++;
	  std::string dir_path_str(in_data);
	  dir_path_str.append(en->d_name);
	  //Convert file path string back to char to be used by read_Hipo    
	  int string_l = dir_path_str.length();
	  char dir_path_char[256];
	  strcpy(dir_path_char, dir_path_str.c_str());
	  std::cout << "input file = " << dir_path_char << std::endl;

	  try {
	    read_Hipo(dir_path_char, part_charge, part_mass, partInfo);
	  }
	  catch (const std::exception& e) {
	    std::cerr << "Error reading HIPO file: " << dir_path_char << std::endl;
	    continue; // Skip to the next file
	  }
	  //Print list of all read hipo files to input_files.txt   
	  //InFileList << in_data << en->d_name << std::endl;
	}
      }
      closedir(dr); //close all directory
      //std::cout << in_file_count << " files read" << std::endl;
      //std::cout << partInfo.p_diff_vec.size() << " events saved out of " << event_count << " events read." << std::endl;
    }
    //If in_data is not a directory, then get the single input file    
    //else {
    //  std::cout << "Not a directory. Only single input file: " << in_data << std::endl;
    //  read_Hipo(in_data, part_charge, part_mass, &mc_p_vec, &rec_p_vec, &mc_theta_deg_vec, &rec_theta_deg_vec, &mc_phi_deg_vec, &rec_phi_deg_vec, &p_diff_vec, &theta_diff_rad_vec, &phi_diff_rad_vec, &chi2_vec, &ndf_vec, event_count);
    //std::cout << p_diff_vec.size() << " events saved out of " << event_count << " events read." << std::endl;
    //}
  }


  //  std::vector<double> rec_p_this_bin;
  //  rec_p_this_bin.reserve(tot_N_bins);


  //std::cout << partInfo.mc_p_vec[1] << std::endl;
  //kinInfo_tree->Branch("mc_p", &partInfo.mc_p_vec);
  kinInfo_tree->Fill();
  // Save the tree to the ROOT file
  kinInfo_tree->Write();

  // Close the ROOT file
  kinInfo_file->Close();
  
  return 0;
}


void read_Hipo(char inputFile[256], int part_charge, float part_mass, partInfoStruct &partInfo)
{
  TLorentzVector target(0.0, 0.0, 0.0, 0.0);
  TLorentzVector beam(0.0, 0.0, 10.45, 11.004);
  TLorentzVector W = beam + target;

  hipo::reader  reader;
  reader.open(inputFile);
  hipo::dictionary  factory;
  reader.readDictionary(factory);
  hipo::event      event;
  hipo::bank RECPART(factory.getSchema("REC::Particle"));
  hipo::bank MCPART(factory.getSchema("MC::Particle"));
  hipo::bank MC_MATCH(factory.getSchema("MC::GenMatch"));
      
  int event_num;
  int good_events = 0;
  int total_events = 0;
  int rejected_events = 0;
  //while (nevents < max_events)
  while(reader.next()==true)
  //while(reader.next()==true && (total_events < 1000))
  {
    total_events++;
    reader.read(event);
    event.getStructure(RECPART);
    event.getStructure(MCPART);
    event.getStructure(MC_MATCH);
    
    //Get reconstructed particle info
    std::vector<int> rec_pid;
    std::vector<float> rec_px;
    std::vector<float> rec_py;
    std::vector<float> rec_pz;
    std::vector<int> rec_status;
    std::vector<int> rec_sector;
    std::vector<float> chi2;
    std::vector<int> ndf;
    read_Part_Bank(RECPART, MC_MATCH, &rec_pid, &rec_px, &rec_py, &rec_pz, &rec_sector, &chi2, &ndf);

    //Get generated particle info
    std::vector<int> mc_pid;
    std::vector<float> mc_px;
    std::vector<float> mc_py;
    std::vector<float> mc_pz;
    read_Part_Bank(MCPART, &mc_pid, &mc_px, &mc_py, &mc_pz);
    
    //Get the kinematics
    TLorentzVector mc;
    TLorentzVector rec;
    if (rec_pid.size() == mc_pid.size()) {
	float rec_E = sqrt(rec_px[0]*rec_px[0] + rec_py[0]*rec_py[0] + rec_pz[0]*rec_pz[0] + part_mass*part_mass);
	float mc_E= sqrt(mc_px[0]*mc_px[0] + mc_py[0]*mc_py[0] + mc_pz[0]*mc_pz[0] + part_mass*part_mass);
	rec.SetPxPyPzE(rec_px[0], rec_py[0], rec_pz[0], rec_E);
	mc.SetPxPyPzE(mc_px[0], mc_py[0], mc_pz[0], mc_E);
	float rec_p = rec.P();
	float mc_p = mc.P();
	float p_diff = rec_p - mc_p;
	float mc_MM2 = (W - mc).M2();
	float rec_MM2 = (W - rec).M2();
	float deltaP_over_p = p_diff / rec_p;
	
	//Theta and Phi
	float rec_theta_rad = rec.Theta();
	float rec_theta_deg = rec_theta_rad*180/TMath::Pi();
	float rec_phi_rad = rec.Phi();
	float rec_phi_deg = rec_phi_rad*180/TMath::Pi();

	float mc_theta_rad = mc.Theta();
	float mc_theta_deg = mc_theta_rad*180/TMath::Pi();
	float mc_phi_rad = mc.Phi();
	float mc_phi_deg = mc_phi_rad*180/TMath::Pi();
	float theta_diff_rad = rec_theta_rad - mc_theta_rad;
	float theta_diff_deg = rec_theta_deg - mc_theta_deg;
	float phi_diff_rad = rec_phi_rad - mc_phi_rad;
	float phi_diff_deg = rec_phi_deg - mc_phi_deg;

	//Save the kinematic info
	partInfo.rec_p_vec.push_back(rec_p);
	partInfo.mc_p_vec.push_back(mc_p);
	partInfo.chi2_vec.push_back(chi2[0]);
	//partInfo.ndf_vec.push_back(ndf[0]);
	partInfo.rec_theta_deg_vec.push_back(rec_theta_deg);
	partInfo.mc_theta_deg_vec.push_back(mc_theta_deg);
	partInfo.rec_phi_deg_vec.push_back(rec_phi_deg);
	partInfo.mc_phi_deg_vec.push_back(mc_phi_deg);
    }      
  }
}
//This is intended to be used for reading the MC info
void read_Part_Bank(hipo::bank PartBank, std::vector<int>* pid, std::vector<float>* px, std::vector<float>* py, std::vector<float>* pz){  
  int nrows = PartBank.getRows();
  for(int i = 0; i < nrows; i++){
    int   current_pid = PartBank.getInt("pid",i);
    float  current_px = PartBank.getFloat("px",i);
    float  current_py = PartBank.getFloat("py",i);
    float  current_pz = PartBank.getFloat("pz",i);
    pid->push_back(current_pid);
    px->push_back(current_px);
    py->push_back(current_py);
    pz->push_back(current_pz);
  }  
}
//This is intended to be used for reading Rec particle info 
void read_Part_Bank(hipo::bank PartBank, hipo::bank MCMatchBank, std::vector<int>* pid, std::vector<float>* px, std::vector<float>* py, std::vector<float>* pz, std::vector<int> *rec_sector, std::vector<float> *chi2, std::vector<int> *ndf){

  // fail safe for truth matching 
  if (MCMatchBank.getRows() < 1)
    return;

  // This assumes only 1 particle in MC::GenMatch bank
  int index_truth = MCMatchBank.getInt("pindex", 0);
  int nrows = PartBank.getRows();

  //event_count->push_back(nrows);

  //Temporary vectors
  std::vector<float> px_vec;
  std::vector<float> py_vec;
  std::vector<float> pz_vec;
  std::vector<int> pid_vec;
  std::vector<int> sector_vec;
  std::vector<float> chi2_vec;
  std::vector<int> ndf_vec;
  std::vector<int> status_vec;
  std::vector<int> charge_vec;
  //int clusters_recorded = 0;

  //int count = 0;
  float current_chi2;
  int current_ndf;
  float current_chi2_ndf;
  int current_sec;
  int   current_pid = PartBank.getInt("pid",index_truth);
  float  current_px = PartBank.getFloat("px",index_truth);
  float  current_py = PartBank.getFloat("py",index_truth);
  float  current_pz = PartBank.getFloat("pz",index_truth);
  float  current_chi2pid = PartBank.getFloat("chi2pid",index_truth);
  int current_status = PartBank.getInt("status",index_truth);
  int charge = PartBank.getInt("charge",index_truth);
  //std::cout << current_pid << std::endl;
    
  px->push_back(current_px);      
  py->push_back(current_py);   
  pz->push_back(current_pz);   
  pid->push_back(current_pid);
  chi2->push_back(current_chi2pid);
  //ndf->push_back(current_ndf);
}

void get_track_parameters(int index, const hipo::bank TrackBank, int& track_sec, float& track_chi2, int& track_ndf)
{
  std::map<int, std::vector<int>> trackBankMap = loadMapByIndex(TrackBank, "id");
  auto it = trackBankMap.find(index);
  if (it != trackBankMap.end())
    {
      const std::vector<int> &indices = it->second;
      if (!indices.empty())
	{
	  // This will always take the first column whose "pindex" value matches the given "index" 
	  track_sec = TrackBank.getInt("sector", indices[0]);
	  track_chi2 = TrackBank.getFloat("chi2", indices[0]);
	  track_ndf = TrackBank.getInt("ndf", indices[0]);
	}
    }
  return; // Index not found in the map                                                                             
}

int main() {
    gErrorIgnoreLevel = kBreak;      // Ignore errors and warnings but display abort messages
    return covMatrix_extraction();
}

