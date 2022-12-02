#include <vector>
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "KinFitter.h"
#include "KinParticle.h"
#include "KinFitTest.h"
#include "reader.h"
#include <dirent.h>
#include <sys/types.h>
#include <cstring>


void read_Hipo(char input_hipo[256], std::vector<int> required_pids, std::vector<double> masses, TLorentzVector target, TLorentzVector beam, KinFitTest test, std::vector<TString> parts);

int test_4C_hipo()
{
  TLorentzVector target(0.0, 0.0, 0.0, 0.938272);                                                                                                                                                    
  TLorentzVector beam(0.0, 0.0, 10.6, 10.6);
  TLorentzVector W = beam + target;
  const std::vector<TString> parts = {"e", "p", "#pi^{+}", "#pi^{-}"};
  const std::vector<double> masses = {0.938, 0.139, 0.139};                                                                                                                                        
  const std::vector<int> required_pids = {11, -211, 211, 2212};

  KinFitTest test("4C", parts, W);

  //The variable dir_path can be specified to either be the path to the output directory containing                                                                                                  
  //multiple osg outputs or can be full path to a specific hipo file.
  std::ofstream InFileList("input_files.txt");

  //---------------Specify input file or directory path---------------//
  char dir_path[256] = "/volatile/clas12/osg2/reedtg/job_5429/output/simu_0/dst.hipo";
  //char dir_path[256] = "/volatile/clas12/osg2/reedtg/job_5429/output/";
  //------------------------------------------------------------------//

  DIR *dr;
  struct dirent *en;
  dr = opendir(dir_path); //open all or present directory
  
//If dir_path is a directory, loop through all files within it
  int in_file_count = 0;
  if (dr) {
    while ((en = readdir(dr)) != NULL) {
      if ( (strcmp(en->d_name, ".") != 0) && (strcmp(en->d_name, "..") != 0)) {
	in_file_count++;
	std::string dir_path_str(dir_path);
        dir_path_str.append(en->d_name);
        dir_path_str.append("/dst.hipo");
	//Convert file path string back to char to be used by read_Hipo
	int string_l = dir_path_str.length();
	char dir_path_char[256];
	strcpy(dir_path_char, dir_path_str.c_str());
	read_Hipo(dir_path_char, required_pids, masses, target, beam, test, parts);
	//Print list of all read hipo files to input_files.txt
	InFileList << dir_path << en->d_name << "/dst.hipo" << std::endl;
      }
    }
      closedir(dr); //close all directory
  }
  else {
    std::cout << "Not a directory. Only single file " << dir_path << std::endl;
    read_Hipo(dir_path, required_pids, masses, target, beam, test, parts);
  }
  test.plot();

  return 0;
}
    
std::ofstream InfoFile("event_info.txt");
    
void read_Hipo(char inputFile[256], std::vector<int> required_pids, std::vector<double> masses, TLorentzVector target, TLorentzVector beam, KinFitTest test, std::vector<TString> parts)
{
  TLorentzVector W = beam + target;

  std::vector<double> resolutions;
  std::vector<int> constraint_idx;
  for (int i = 0; i < parts.size(); i++)
    {
      constraint_idx.push_back(i);
      for (int j = 0; j < KINES.size(); j++)
        {
	  resolutions.push_back(RESO[j]);
        }
    }
  //std::cout << "Opening input file " << std::endl;
  hipo::reader  reader;
  reader.open(inputFile);
  hipo::dictionary  factory;
  reader.readDictionary(factory);
  hipo::event      event;
  hipo::bank PART(factory.getSchema("REC::Particle"));
  hipo::bank COV(factory.getSchema("REC::CovMat"));
  hipo::bank MCPART(factory.getSchema("MC::Particle"));
  hipo::bank RUN(factory.getSchema("RUN::config"));
      
  int event_num;
  int good_events = 0;
  int total_events = 0;
  //while (nevents < max_events)
  while(reader.next()==true)
  //while(reader.next()==true && (total_events < 10))
  {
    total_events++;
    reader.read(event);
    event.getStructure(PART);
    event.getStructure(COV);
    event.getStructure(MCPART);
    event.getStructure(RUN);
	  
    //Get reconstructed particle info  -- I should make a method to read either Particle bank
    std::vector<int> pid_list;
    std::vector<int> input_part_index;
    std::vector<float> px_list;
    std::vector<float> py_list;
    std::vector<float> pz_list;
    std::vector<TLorentzVector> parts_rec;
    int nrows_rec = PART.getRows();
    for(int i = 0; i < nrows_rec; i++){
      int   pid = PART.getInt("pid",i);
      float  px = PART.getFloat("px",i);
      float  py = PART.getFloat("py",i);
      float  pz = PART.getFloat("pz",i);
      for (int j = 0; j < required_pids.size(); j++){
	if (pid == required_pids.at(j)){
	  pid_list.push_back(pid);
	  input_part_index.push_back(j);
	  px_list.push_back(px);
	  py_list.push_back(py);
	  pz_list.push_back(pz);
	}
      }
    }
	  
    //Get generated particle info
    std::vector<int> mc_pid_list;
    std::vector<int> mc_input_part_index;
    std::vector<float> mc_px_list;
    std::vector<float> mc_py_list;
    std::vector<float> mc_pz_list;
    std::vector<TLorentzVector> parts_mc;
    int nrows_mc = MCPART.getRows();
    for(int i = 0; i < nrows_mc; i++){
      int mc_pid = MCPART.getInt("pid",i);
      float mc_px = MCPART.getFloat("px",i);
      float mc_py = MCPART.getFloat("py",i);
      float mc_pz = MCPART.getFloat("pz",i);
      for (int j = 0; j < required_pids.size(); j++){
        if (mc_pid == required_pids.at(j)){
	  mc_pid_list.push_back(mc_pid);
	  mc_input_part_index.push_back(j);
	  mc_px_list.push_back(mc_px);
	  mc_py_list.push_back(mc_py);
	  mc_pz_list.push_back(mc_pz);
	}
      }
    }
	  
    //Check if event (reconstructed) contains exactly one of each of required particles
    bool contains_required_parts = false;
    for (int k = 0; k < required_pids.size(); k++){
      int pid_count = count(pid_list.begin(), pid_list.end(), required_pids.at(k));
      //InfoFile << "pid count " << pid_count << std::endl;
      if (pid_count == 1){
	contains_required_parts = true;
	//InfoFile << "PID from 'required pid' list: " << required_pids.at(k) << std::endl;
      }
      else {
	contains_required_parts = false;
	//InfoFile << "Break from event for-loop" << std::endl;
	break;
      }
    }

    //Proceed to kinematic fitting if reconstructed event meets particle requirements
    if (contains_required_parts){
      good_events++;
      //std::cout << "Good event" << std::endl;
      event_num = RUN.getInt("event", 0);
      //InfoFile << "Event Number: " << event_num << std::endl;
      std::vector<KinParticle> kin_parts;
      TLorentzVector current_part_vec;
      TLorentzVector current_part_mc_vec;
      for (int ipart = 0; ipart < pid_list.size(); ++ipart)
      {
	double px_part = px_list[ipart], py_part = py_list[ipart], pz_part = pz_list[ipart];
	double current_part_P_mag = sqrt(px_part*px_part + py_part*py_part + pz_part*pz_part);
	double current_part_E = sqrt(current_part_P_mag*current_part_P_mag + masses[input_part_index[ipart]]*masses[input_part_index[ipart]]);
	current_part_vec.SetPxPyPzE(px_part, py_part, pz_part, current_part_E);
	parts_rec.push_back(current_part_vec);
	kin_parts.push_back(KinParticle(current_part_vec, masses[input_part_index[ipart]], RESO));
	//InfoFile << "Reconstructed particle (PID " << pid_list[ipart] << ") 3-momentum: " << px_part << " " <<  py_part << " " <<  pz_part << std::endl;
      }
      for (int ipart = 0; ipart < mc_pid_list.size(); ++ipart)
      {
	double mc_px_part = mc_px_list[ipart], mc_py_part = mc_py_list[ipart], mc_pz_part = mc_pz_list[ipart];
	double current_mc_part_P_mag= sqrt(mc_px_part*mc_px_part + mc_py_part*mc_py_part + mc_pz_part*mc_pz_part);
	double current_mc_part_E = sqrt(current_mc_part_P_mag*current_mc_part_P_mag + masses[input_part_index[ipart]]*masses[input_part_index[ipart]]);
	current_part_mc_vec.SetPxPyPzE(mc_px_part, mc_py_part, mc_pz_part, current_mc_part_E);
	parts_mc.push_back(current_part_mc_vec);
      }
	    
      //InfoFile << "Reconstructed particle 4-vector: " << current_part_vec[0] << " " << current_part_vec[1] << " " << current_part_vec[2] << " " << current_part_vec[3] << std::endl;
      //InfoFile << "Reconstructed PID list: " << pid_list[0] << " " << pid_list[1] << " " << pid_list[2] << " " << pid_list[3] << std::endl;
      auto kin = new KinFitter({KinParticle(target), KinParticle(beam)}, kin_parts);
      kin->Add_EnergyMomentum_Constraint(constraint_idx);
      kin->DoFitting(100);
      bool is_background = false;
      double weight = 1.0;
      test.fill_MissingMass(kin, parts_mc, parts_rec, weight, is_background);
    }
  }
}

int main() {
    return test_4C_hipo();
}

