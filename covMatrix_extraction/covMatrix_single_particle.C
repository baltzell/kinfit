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
  std::vector<float> p_diff_vec, theta_diff_deg_vec, phi_diff_deg_vec, chi2_vec, ndf_vec; 
};

const int P_bins = 25;
const int theta_bins = 25;
const int phi_bins = 25;

const double P_bin_min = 0.2;
const double P_bin_max = 12;
const double P_bin_size = (P_bin_max - P_bin_min)/P_bins;
const double theta_bin_min = 5;    // Same for both pi- and pi+. Theta max gets initialized at the beginning of covMatrix_extraction() depeneding on particle type
//const double theta_bin_max = 35;
//const double theta_bin_size = (theta_bin_max - theta_bin_min) / theta_bins;
const double phi_bin_min = 30;
const double phi_bin_max = 90;
const double phi_bin_size = (phi_bin_max - phi_bin_min) / phi_bins;


// Struct for storing binned kinematic resolutions
struct binnedKinInfo
{
  //4-dimensional vectors for saving kinematic info. From innermost outward: (Rec - MC) for given kinematic, phi bin, theta bin, and then P bin
  std::vector<std::vector<std::vector<std::vector<double>>>> p_diff_binned;
  std::vector<std::vector<std::vector<std::vector<double>>>> phi_diff_binned;
  std::vector<std::vector<std::vector<std::vector<double>>>> theta_diff_binned;
  std::vector<std::vector<std::vector<std::vector<double>>>> rec_p_binned;
  std::vector<std::vector<std::vector<std::vector<double>>>> deltaP_over_p_binned;

  // Constructor to initialize vectors
  binnedKinInfo() :
    p_diff_binned(P_bins, std::vector<std::vector<std::vector<double>>>(theta_bins, std::vector<std::vector<double>>(phi_bins, std::vector<double>()))),
    phi_diff_binned(P_bins, std::vector<std::vector<std::vector<double>>>(theta_bins, std::vector<std::vector<double>>(phi_bins, std::vector<double>()))),
    theta_diff_binned(P_bins, std::vector<std::vector<std::vector<double>>>(theta_bins, std::vector<std::vector<double>>(phi_bins, std::vector<double>()))),
    rec_p_binned(P_bins, std::vector<std::vector<std::vector<double>>>(theta_bins, std::vector<std::vector<double>>(phi_bins, std::vector<double>()))),
    deltaP_over_p_binned(P_bins, std::vector<std::vector<std::vector<double>>>(theta_bins, std::vector<std::vector<double>>(phi_bins, std::vector<double>())))
  {}
};



void read_Hipo(char in_data[256], int part_charge, float part_mass, partInfoStruct &partInfo, int& event_count, binnedKinInfo &binnedKins, TH3* P_theta_and_phi);
void read_MC_Part_Bank(hipo::bank PartBank, std::vector<int>* pid, std::vector<float>* px, std::vector<float>* py, std::vector<float>* pz);
void read_Rec_Part_Bank(hipo::bank PartBank, hipo::bank MCMatchBank, hipo::bank TBTrack, std::vector<int>* pid, std::vector<float>* px, std::vector<float>* py, std::vector<float>* pz, std::vector<int> *rec_status, std::vector<float> *chi2, std::vector<int> *ndf, int& event_count);
void get_track_parameters(int index, const hipo::bank TrackBank, int& track_sec, float& track_chi2, int& track_ndf);


//Output file for saving any information for trouble shooting 
std::ofstream InfoFile("event_info.txt");

// Function to get index of an element in vector
int getIndex(std::vector<int> v, int K) {
  int index = -1;
  for (int i = 0; i < v.size(); i++) {
    if (v[i] == K) {
      index = i;
    }
  }
  return index;
}


//========================== Fitting Functions ==============================//
// Quadratic background function
Double_t p2Background(Double_t *x, Double_t *par) {
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}
//Linear background function  
Double_t p1Background(Double_t *x, Double_t *par) {
  return par[0] + par[1]*x[0];
}
//Gaussian Fit
Double_t GaussianPeak(Double_t *x, Double_t *par)
{
  Double_t arg = 0;
  arg = (x[0] - par[1])/par[2];
  Double_t sig = par[0]*TMath::Exp(-0.5*arg*arg);
  return sig;
}
//Sum of p2 background and peak function
Double_t fitFunction_p2_bg(Double_t *x, Double_t *par) {
  return p2Background(x,par) + GaussianPeak(x,&par[3]);
}
//Sum of p1 background and peak function
Double_t fitFunction_p1_bg(Double_t *x, Double_t *par) {
  return p1Background(x,par) + GaussianPeak(x,&par[2]);
}
// Get parameters from a gaussian fit for initial params in total fit
Double_t* gausfit_init(TH1 *hist, float min, float max) {
  TF1 *gausfit = new TF1("gausfit", "gaus", min, max);
  hist->Fit("gausfit", "0", "", min, max);
  //gausfit->GetParameters(gausfit_params);
  Double_t* gausfit_params = gausfit->GetParameters();
  return gausfit_params;
}
//===========================================================================//


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

void zoomOutliers(TH1* h, float n_rms) {
  const float mpv = h->GetBinCenter(h->GetMaximumBin());
  const float rms = h->GetRMS();
  const int min_bin = h->FindBin(mpv - n_rms*rms);
  const int max_bin = h->FindBin(mpv + n_rms*rms);
  h->GetXaxis()->SetRange(min_bin,max_bin);
}

void iqr_zoom(TH1* h, float iqr, float q25, float q75, float N) {
  const int min_bin = h->FindBin(q25 - N * iqr);
  const int max_bin = h->FindBin(q75 + N * iqr);

  // Get the current range
  h->GetXaxis()->SetRange(min_bin,max_bin);
  const int min_bin_current = h->GetXaxis()->GetFirst();
  const int max_bin_current = h->GetXaxis()->GetLast();
  
  // Only perform change of histogram range if it's being decreased
  if (min_bin > min_bin_current & max_bin < max_bin_current) {
    h->GetXaxis()->SetRange(min_bin,max_bin);
  }
  
}

struct CovarianceResult {
  const TMatrixDSym* covariance;
  const TVectorD* means;
  int status;
  int good_counts;
};

CovarianceResult calculateCovarianceRobust(const std::vector<double>& p_diff, const std::vector<double>& theta_diff, const std::vector<double>& phi_diff, double p_low, double p_upp, double theta_low, double theta_upp, double phi_low, double phi_upp) {
  std::vector<double> p_good;
  std::vector<double> th_good;
  std::vector<double> phi_good;
  std::vector<Double_t> row;
  TRobustEstimator estimator = TRobustEstimator(3,3);
  Int_t numValues = p_diff.size(); // Get the number of values in the arrays (any of the kinematic variables can be used for this)
  int good_counts = 0;  //good_counts are those that are within the kinematic (rec - MC) bounds, to remove outliers 
  for (int j = 0; j < numValues; j++) {
    if ((p_diff[j] > p_low && p_diff[j] < p_upp) && (theta_diff[j] > theta_low && theta_diff[j] < theta_upp) && (phi_diff[j] > phi_low && phi_diff[j] < phi_upp)) {
      good_counts++;
      p_good.push_back(p_diff[j]);
      th_good.push_back(theta_diff[j]);
      phi_good.push_back(phi_diff[j]);
      std::vector<Double_t> row = {p_diff[j], theta_diff[j], phi_diff[j]};
      estimator.AddRow(row.data());
    }
  }

  estimator.Evaluate();
  CovarianceResult result;
  //const TMatrixDSym* covmat = estimator.GetCovariance();
  //result.covariance = covmat;
  TMatrixDSym covmat = *(estimator.GetCovariance());
  result.covariance = new TMatrixDSym(covmat);
  if (result.covariance) {
    result.status = 0;
  }
  else {
    result.status = 2;
  }
  TVectorD means = *(estimator.GetMean());
  result.means = new TVectorD(means);
  Int_t N_obs = estimator.GetNumberObservations();
  Int_t N_out = estimator.GetNOut();
  TArrayI outlier_array = *(estimator.GetOuliers());
  //std::cout << "N_obs: " << N_obs << std::endl;
  //std::cout << "N_out: " << N_out << std::endl;
  result.good_counts = good_counts;
  //for (int i = 0; i < outlier_array.GetSize(); i++) {
  //  std::cout << outlier_array[i] << std::endl;
  //}

//result.covariance->Print();
  return result;
}

// Define a struct to hold input parameters for calculateCovariance
struct CovInputParams {
  int kin_id1;
  int kin_id2;
  double kinvar1_mean;
  double kinvar2_mean;
  double p_bin_center;
  std::vector<double> p_diff;
  std::vector<double> theta_diff;
  std::vector<double> phi_diff;
  double p_low;
  double p_upp;
  double theta_low;
  double theta_upp;
  double phi_low;
  double phi_upp;

  // Constructor for initialization of input params
  CovInputParams(int id1, int id2, double mean1, double mean2, double bin_center,
		 const std::vector<double>& pdiff, const std::vector<double>& thetadiff, const std::vector<double>& phidiff,
		 double plow, double pupp, double thetalow, double thetaupp, double philow, double phiupp)
    : kin_id1(id1), kin_id2(id2), kinvar1_mean(mean1), kinvar2_mean(mean2), p_bin_center(bin_center),
      p_diff(pdiff), theta_diff(thetadiff), phi_diff(phidiff),
      p_low(plow), p_upp(pupp), theta_low(thetalow), theta_upp(thetaupp), phi_low(philow), phi_upp(phiupp) {}
};

void calculateCovariance(double& covariance, int& good_counts, const CovInputParams& params) {
  std::vector<double> x;
  std::vector<double> y;
  int calculate_p_mean = 0;
  double p_sum = 0;
  double p_mean = 0;

  if (params.kin_id1 == 1) {
    x = params.p_diff;
  }
  else if (params.kin_id1 == 2) {
    x = params.theta_diff;
  }
  else if (params.kin_id1 == 3) {
    x = params.phi_diff;
  }
  if (params.kin_id2 == 1) {
    y = params.p_diff;
  }
  else if (params.kin_id2 == 2){
    y = params.theta_diff;
  }
  else if (params.kin_id2 == 3) {
    y = params.phi_diff;
  }

  Int_t numValues = x.size(); // Get the number of values in the arrays (any of the kinematic variables can be used for this)   
  covariance = 0;
  good_counts = 0;  //good_counts are those that are within the kinematic (rec - MC) bounds, to remove outlier
  for (int i = 0; i < numValues; i++) {
    if ((params.p_diff[i] > params.p_low && params.p_diff[i] < params.p_upp) && 
	(params.theta_diff[i] > params.theta_low && params.theta_diff[i] < params.theta_upp) && 
	(params.phi_diff[i] > params.phi_low && params.phi_diff[i] < params.phi_upp)) {
      covariance += (x[i] - params.kinvar1_mean) * (y[i] - params.kinvar2_mean);
      good_counts++;
    }
  }
  if (good_counts > 0) {
    covariance /= good_counts;
  } 
  else {
    covariance = 0;
  }
}



//Initialize the particle type and binning info
const TString part_type = "#pi^{+}";
//const TString part_type = "#pi^{-}";
float part_mass;
int part_charge;

//Set boundaries for the kinematics. This is used to set the initial bounds for the (rec - mc) distributions of the kinematics
const double p_low_init = -1;
const double p_upp_init = 1;
const double theta_low_init = -0.1 * 180 / 3.14159;
const double theta_upp_init = 0.1 * 180 / 3.14159;
const double phi_low_init = -0.2 * 180 / 3.14159;
const double phi_upp_init = 0.2 * 180 / 3.14159;

//These variables are used for the plots that combine all P, theta, and phi bins
const std::vector<TString> KINES = {"P", "#theta", "#phi"};
const std::vector<TString> UNITS = {"[GeV]", "[deg]", "[deg]"};
const std::vector<int> kin_plots_low = {0, 0, 25};
const std::vector<int> kin_plots_high = {12, 40, 95};
const std::vector<double> kin_delta_plots_low = {p_low_init, theta_low_init, phi_low_init};
const std::vector<double> kin_delta_plots_high = {p_upp_init, theta_upp_init, phi_upp_init};

//Use these values to set a maximum on the number of files that the reader loop will read and to set a number of elements for the vectors to be initilaized to
const int max_files = 10000;
const int N_events_per_file = 10000;
const int N_events = max_files*N_events_per_file;
const TString out_name = Form("pip_%i_in_files", max_files);

int covMatrix_extraction()
{
  //---------------Specify input file or directory path---------------//      
  
  // Pi-
  //std::vector<std::string> input_directories = {"/volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/pim-sec2-100mil_events_6-15-23/cooked/"};
  //std::vector<std::string> input_directories = {"volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/pim_sec2_100mil_events/testdir/"};
  //std::vector<std::string> input_directories = {"volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/pim_sec2_100mil_events/cooked/",
  //						"/volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/pim_sec2_10mil_events/cooked/",
  //						"/volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/pim_sec2_50mil_events/cooked/"};
  
  // Pi+
  //std::vector<std::string> input_directories = {"/volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/pip-sec2-100mil_events_6-13-23/cooked/"};
  //std::vector<std::string> input_directories = {"/volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/pip-sec2_small_region/cooked/",
  //						"/volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/pip-sec2_small_region_2/cooked/"};
  //std::vector<std::string> input_directories = {"/volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/pip-sec2_medium_region/cooked/"};
  //std::vector<std::string> input_directories = {
  //  "/volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/pip-sec2-100mil_events/cooked/",
  //  "/volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/pip-sec2-100mil_events_6-13-23/cooked/"};
  //std::vector<std::string> input_directories = {"/volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/pip-sec2-100mil_events_6-13-23/cooked/"};
  std::vector<std::string> input_directories = {"/volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/pip_sec2_truth/cooked/"};
  
  //-----------------------------------------------------------------//

  partInfoStruct partInfo;
  binnedKinInfo binnedKins;
  int event_count = 0;

  // Initialize particle charge and mass
  if ((strcmp(part_type, "#pi^{+}") == 0) || (strcmp(part_type, "#pi^{-}") == 0)) {
    part_mass = 0.13957;
    if (strcmp(part_type, "#pi^{+}") == 0) {
      part_charge = 1;
      theta_bin_max = 37;
    }
    else if (strcmp(part_type, "#pi^{-}") == 0) {
      part_charge = -1;
      theta_bin_max = 40;
    }
  }
  const double theta_bin_size = (theta_bin_max - theta_bin_min) / theta_bins; 

  partInfo.p_diff_vec.reserve(N_events);
  partInfo.theta_diff_deg_vec.reserve(N_events);
  partInfo.phi_diff_deg_vec.reserve(N_events);
  partInfo.mc_p_vec.reserve(N_events);
  partInfo.rec_p_vec.reserve(N_events);
  partInfo.mc_theta_deg_vec.reserve(N_events);
  partInfo.rec_theta_deg_vec.reserve(N_events);
  partInfo.mc_phi_deg_vec.reserve(N_events);
  partInfo.rec_phi_deg_vec.reserve(N_events);
  partInfo.chi2_vec.reserve(N_events);
  partInfo.ndf_vec.reserve(N_events);


  //Set up TH3 to store the binning info 
  TH3* P_theta_and_phi = new TH3D("P_theta_and_phi", "P_theta_and_phi", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);


  int in_file_count = 0;
  for (const std::string& in_data : input_directories) {
    DIR *dr;
    struct dirent *en;
    //dr = opendir(in_data);                                
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
	    read_Hipo(dir_path_char, part_charge, part_mass, partInfo, event_count, binnedKins, P_theta_and_phi);
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
      std::cout << in_file_count << " files read" << std::endl;
      std::cout << partInfo.p_diff_vec.size() << " events saved out of " << event_count << " events read." << std::endl;
    }
    //If in_data is not a directory, then get the single input file    
    //else {
    //  std::cout << "Not a directory. Only single input file: " << in_data << std::endl;
    //  read_Hipo(in_data, part_charge, part_mass, &mc_p_vec, &rec_p_vec, &mc_theta_deg_vec, &rec_theta_deg_vec, &mc_phi_deg_vec, &rec_phi_deg_vec, &p_diff_vec, &theta_diff_rad_vec, &phi_diff_rad_vec, &chi2_vec, &ndf_vec, event_count);
    //std::cout << p_diff_vec.size() << " events saved out of " << event_count << " events read." << std::endl;
    //}
  }


  // Create a ROOT file to store the data
  TFile* C_file = new TFile(Form("cov_root_files/covariances_%s.root", out_name.Data()), "RECREATE");

  // Create a TTree
  TTree* C_tree = new TTree("covarianceTree", "Covariance Tree");
  double P, theta, phi, C_P, C_P_err, C_theta, C_theta_err, C_phi, C_phi_err, C_P_phi, C_P_phi_err, C_P_theta, C_P_theta_err, C_theta_phi, C_theta_phi_err, entries, P_offset, theta_offset, phi_offset;
  C_tree->Branch("P", &P, "P/D");
  C_tree->Branch("theta", &theta, "theta/D");
  C_tree->Branch("phi", &phi, "phi/D");
  C_tree->Branch("C_P", &C_P, "C_P/D");
  C_tree->Branch("C_theta", &C_theta, "C_theta/D");
  C_tree->Branch("C_phi", &C_phi, "C_phi/D");
  C_tree->Branch("C_P_phi", &C_P_phi, "C_P_phi/D");
  C_tree->Branch("C_P_theta", &C_P_theta, "C_P_theta/D");
  C_tree->Branch("C_theta_phi", &C_theta_phi, "C_theta_phi/D");
  C_tree->Branch("C_P_err", &C_P_err, "C_P_err/D");
  C_tree->Branch("C_theta_err", &C_theta_err, "C_theta_err/D");
  C_tree->Branch("C_phi_err", &C_phi_err, "C_phi_err/D");
  C_tree->Branch("C_P_phi_err", &C_P_phi_err, "C_P_phi_err/D");
  C_tree->Branch("C_P_theta_err", &C_P_theta_err, "C_P_theta_err/D");
  C_tree->Branch("C_theta_phi_err", &C_theta_phi_err, "C_theta_phi_err/D");
  C_tree->Branch("event_count", &event_count, "event_count/I");
  C_tree->Branch("P_offset", &P_offset, "P_offset/D");
  C_tree->Branch("theta_offset", &theta_offset,"theta_offset/D");
  C_tree->Branch("phi_offset", &phi_offset,"phi_offset/D");

  //Create a TH3 for each covariance parameter and each covariance error parameter
  TH3D* C_P_hist = new TH3D("C_P_hist", "C_P Histogram", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);
  TH3D* C_phi_hist = new TH3D("C_phi_hist", "C_phi Histogram", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);
  TH3D* C_theta_hist = new TH3D("C_theta_hist", "C theta_Hist", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);
  TH3D* C_P_phi_hist = new TH3D("C_P_phi_hist", "C_P_phi Hist", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);
  TH3D* C_P_theta_hist = new TH3D("C_P_theta_hist", "C_P_theta Hist", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);
  TH3D* C_theta_phi_hist = new TH3D("C_theta_phi_hist", "C_theta_phi Hist", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);
  TH3D* C_P_err_hist = new TH3D("C_P_err_hist", "C_P_err Histogram", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);
  TH3D* C_phi_err_hist = new TH3D("C_phi_err_hist", "C_phi_err Histogram", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);
  TH3D* C_theta_err_hist = new TH3D("C_theta_err_hist", "C theta_err_Hist", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);
  TH3D* C_P_phi_err_hist = new TH3D("C_P_phi_err_hist", "C_P_phi_err Hist", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);
  TH3D* C_P_theta_err_hist = new TH3D("C_P_theta_err_hist", "C_P_theta_err Hist", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);
  TH3D* C_theta_phi_err_hist = new TH3D("C_theta_phi_err_hist", "C_theta_phi_err Hist", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);
  TH3D* entries_hist = new TH3D("entries_hist", "Entries Histogram", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);
  TH3D* P_offset_hist = new TH3D("P_offset_hist", "P Offset Hist", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);
  TH3D* theta_offset_hist = new TH3D("theta_offset_hist", "Theta Offset Hist", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);
  TH3D* phi_offset_hist = new TH3D("phi_offset_hist", "Phi Offset Hist", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);


  double P_bin_low;
  double P_bin_high;
  double theta_bin_low;
  double theta_bin_high;
  double phi_bin_low;
  double phi_bin_high;
  std::vector<double> p_diff_this_bin;
  std::vector<double> theta_diff_this_bin;
  std::vector<double> phi_diff_this_bin;
  std::vector<double> rec_p_this_bin;
  std::vector<double> deltaP_over_p_this_bin;
  int tot_N_bins = P_bins * theta_bins * phi_bins;
  p_diff_this_bin.reserve(tot_N_bins);
  theta_diff_this_bin.reserve(tot_N_bins);
  phi_diff_this_bin.reserve(tot_N_bins);
  rec_p_this_bin.reserve(tot_N_bins);
  deltaP_over_p_this_bin.reserve(tot_N_bins);
  double variance_P;
  double variance_theta;
  double variance_phi;
  double covariance_P_phi;
  double covariance_P_theta;
  double covariance_theta_phi;
  double variance_P_err;
  double variance_theta_err;
  double variance_phi_err;
  double covariance_P_phi_err;
  double covariance_P_theta_err;
  double covariance_theta_phi_err;
  int N_entries;    //This number could include outliers
  int N_good_events;
  int N_events_C_P;
  int N_events_C_theta;
  int N_events_C_phi;
  int N_events_C_P_theta;
  int N_events_C_P_phi;
  int N_events_C_theta_phi;
  double current_P_bin_center;
  double current_theta_bin_center;
  double current_phi_bin_center;

  // These are the amount by which the mean of the kinematic distributions (Rec - MC) is offset from 0
  double p_diff_offset;
  double theta_diff_offset;
  double phi_diff_offset;
  double holding_var;
  
  //TString delta_kin_pdf_name = Form("plots/rec_mc_diff_plots_%s.pdf", out_name.Data());
  auto delta_kin_canvas = new TCanvas("delta_kin_canvas", "delta_kin_canvas", 800, 1000);
  //delta_kin_canvas->Print(delta_kin_pdf_name + "[");
  TFile* delta_kin_root_file = TFile::Open(Form("delta_kin_plots_root_files/delta_kin_plots_%s.root", out_name.Data()), "RECREATE");

  // Create a directory in the ROOT file to organize your histograms
  TDirectory* histDir = delta_kin_root_file->mkdir("Histograms");
  histDir->cd();
  TDirectory* deltaP_histDir = histDir->mkdir("delta p Histograms");
  TDirectory* deltaTheta_histDir = histDir->mkdir("delta theta Histograms");
  TDirectory* deltaPhi_histDir = histDir->mkdir("delta phi Histograms");
  TDirectory* deltaP_deltaTheta_histDir = histDir->mkdir("delta p delta theta Histograms");
  TDirectory* deltaP_deltaPhi_histDir = histDir->mkdir("delta p delta phi Histograms");
  TDirectory* deltaTheta_deltaPhi_histDir = histDir->mkdir("delta theta delta phi Histograms");

  delta_kin_canvas->Divide(2, 3);
  std::ofstream covMatrix_root_method_txt_file;
  //const int N_bootstrapSamples = 1000;   // This is only used in caluclating the error of the covariance values

  //Loop through the bins and calculate the covariances
  for (int pBin = 0; pBin < P_bins; pBin++) {
    P_bin_low = P_bin_min + pBin*P_bin_size;  
    P_bin_high = P_bin_low + P_bin_size;
    for (int thetaBin = 0; thetaBin < theta_bins; thetaBin++) {
      theta_bin_low = theta_bin_min + thetaBin*theta_bin_size; 
      theta_bin_high = theta_bin_low + theta_bin_size;
      for (int phiBin = 0; phiBin < phi_bins; phiBin++) {
	phi_bin_low = phi_bin_min + phiBin*phi_bin_size;
	phi_bin_high = phi_bin_low + phi_bin_size;

	p_diff_this_bin = binnedKins.p_diff_binned[pBin][thetaBin][phiBin];
	theta_diff_this_bin = binnedKins.theta_diff_binned[pBin][thetaBin][phiBin];
	phi_diff_this_bin = binnedKins.phi_diff_binned[pBin][thetaBin][phiBin];
	rec_p_this_bin = binnedKins.rec_p_binned[pBin][thetaBin][phiBin];
	deltaP_over_p_this_bin = binnedKins.deltaP_over_p_binned[pBin][thetaBin][phiBin];
	//InfoFile << "p_diff_this_bin = " << p_diff_this_bin.size() << std::endl; 

	// ===================================== Make 1D and 2D histograms of rec-mc for each kinematic ========================================= //
	std::vector<TH1*> mc_kin_plots;
	std::vector<TH1*> rec_kin_plots;
	std::vector<TH1*> kin_delta_plots;
	std::vector<TH2*> kin_delta_2D_plots;
	int delta_p_plots_bins = 2000;
	int delta_theta_plots_bins = 1000;
	int delta_phi_plots_bins = 1000;

	//TH1* delta_p_hist = new TH1F(Form("delta_p_plot"), Form(";#frac{#DeltaP}{P};Counts"), delta_p_plots_bins, p_low_init, p_upp_init);
	TH1* delta_p_hist = new TH1F(Form("delta_p_plot"), Form(";#DeltaP;Counts"), delta_p_plots_bins, p_low_init, p_upp_init);
	delta_p_hist->GetXaxis()->SetTitleSize(0.05);
	delta_p_hist->GetYaxis()->SetTitleSize(0.05);
	delta_p_hist->GetXaxis()->SetTitleOffset(0.8);
	delta_p_hist->GetYaxis()->SetTitleOffset(0.9);
	kin_delta_plots.push_back(delta_p_hist);

	TH1* delta_theta_hist = new TH1F(Form("delta_theta_plot"), Form(";#Delta #theta (Rec - Gen) [deg];Counts"), delta_theta_plots_bins, theta_low_init, theta_upp_init);
        delta_theta_hist->GetXaxis()->SetTitleSize(0.05);
        delta_theta_hist->GetYaxis()->SetTitleSize(0.05);
        delta_theta_hist->GetXaxis()->SetTitleOffset(0.8);
        delta_theta_hist->GetYaxis()->SetTitleOffset(0.9);
        kin_delta_plots.push_back(delta_theta_hist);

	TH1* delta_phi_hist = new TH1F(Form("delta_phi_plot"), Form(";#Delta #phi (Rec - Gen) [deg];Counts"), delta_phi_plots_bins, phi_low_init, phi_upp_init);
	delta_phi_hist->GetXaxis()->SetTitleSize(0.05);
        delta_phi_hist->GetYaxis()->SetTitleSize(0.05);
        delta_phi_hist->GetXaxis()->SetTitleOffset(0.8);
        delta_phi_hist->GetYaxis()->SetTitleOffset(0.9);
        kin_delta_plots.push_back(delta_phi_hist);

	//TH2* delta_p_delta_theta_hist = new TH2F(Form("delta_P_delta_theta_plot"), Form(";#Delta #theta (Rec - Gen) [deg];#frac{#DeltaP}{P}"), delta_theta_plots_bins, theta_low_init, theta_upp_init, delta_p_plots_bins, p_low_init, p_upp_init);
	TH2* delta_p_delta_theta_hist = new TH2F(Form("delta_P_delta_theta_plot"), Form(";#Delta #theta (Rec - Gen) [deg];#DeltaP"), delta_theta_plots_bins, theta_low_init, theta_upp_init, delta_p_plots_bins, p_low_init, p_upp_init);
	delta_p_delta_theta_hist->GetXaxis()->SetTitleSize(0.05);
	delta_p_delta_theta_hist->GetYaxis()->SetTitleSize(0.05);
	delta_p_delta_theta_hist->GetXaxis()->SetTitleOffset(0.8);
	delta_p_delta_theta_hist->GetYaxis()->SetTitleOffset(0.9);
	kin_delta_2D_plots.push_back(delta_p_delta_theta_hist);

	TH2* delta_theta_delta_phi_hist = new TH2F(Form("delta_theta_delta_phi_plot"), Form(";#Delta #phi (Rec - Gen) [deg];#Delta #theta (Rec - Gen) GeV"), delta_phi_plots_bins, phi_low_init, phi_upp_init, delta_theta_plots_bins, theta_low_init, theta_upp_init);
        delta_theta_delta_phi_hist->GetXaxis()->SetTitleSize(0.05);
        delta_theta_delta_phi_hist->GetYaxis()->SetTitleSize(0.05);
        delta_theta_delta_phi_hist->GetXaxis()->SetTitleOffset(0.8);
        delta_theta_delta_phi_hist->GetYaxis()->SetTitleOffset(0.9);
        kin_delta_2D_plots.push_back(delta_theta_delta_phi_hist);

	//TH2* delta_p_delta_phi_hist = new TH2F(Form("delta_P_delta_phi_plot"), Form(";#Delta #phi (Rec - Gen) [deg];#frac{#DeltaP}{P}"), delta_phi_plots_bins, phi_low_init, phi_upp_init, delta_p_plots_bins, p_low_init, p_upp_init);
	TH2* delta_p_delta_phi_hist = new TH2F(Form("delta_P_delta_phi_plot"), Form(";#Delta #phi (Rec - Gen) [deg];#DeltaP"), delta_phi_plots_bins, phi_low_init, phi_upp_init, delta_p_plots_bins, p_low_init, p_upp_init);
        delta_p_delta_phi_hist->GetXaxis()->SetTitleSize(0.05);
        delta_p_delta_phi_hist->GetYaxis()->SetTitleSize(0.05);
        delta_p_delta_phi_hist->GetXaxis()->SetTitleOffset(0.8);
        delta_p_delta_phi_hist->GetYaxis()->SetTitleOffset(0.9);
        kin_delta_2D_plots.push_back(delta_p_delta_phi_hist);

        for (int i = 0; i < p_diff_this_bin.size(); i++) {
          kin_delta_plots[0]->Fill(p_diff_this_bin[i]);
	  //kin_delta_plots[0]->Fill(p_diff_this_bin[i] / rec_p_this_bin[i]);
	  //kin_delta_plots[0]->Fill(deltaP_over_p_this_bin[i]);

          kin_delta_plots[1]->Fill(theta_diff_this_bin[i]);
          kin_delta_plots[2]->Fill(phi_diff_this_bin[i]);
          //kin_delta_2D_plots[2]->Fill(phi_diff_this_bin[i], p_diff_this_bin[i] / rec_p_this_bin[i]);
          //kin_delta_2D_plots[0]->Fill(theta_diff_this_bin[i], p_diff_this_bin[i] / rec_p_this_bin[i]);
	  kin_delta_2D_plots[2]->Fill(phi_diff_this_bin[i], p_diff_this_bin[i]);
          kin_delta_2D_plots[0]->Fill(theta_diff_this_bin[i], p_diff_this_bin[i]);
	  //kin_delta_2D_plots[2]->Fill(phi_diff_this_bin[i], deltaP_over_p_this_bin[i]);
	  //kin_delta_2D_plots[0]->Fill(theta_diff_this_bin[i], deltaP_over_p_this_bin[i]);
          kin_delta_2D_plots[1]->Fill(phi_diff_this_bin[i], theta_diff_this_bin[i]);
        }
	// ========================================================================================================================================= //

	//Get bin centers
	current_P_bin_center = (P_bin_high + P_bin_low) / 2;
        current_theta_bin_center = (theta_bin_high + theta_bin_low) / 2;
        current_phi_bin_center = (phi_bin_high + phi_bin_low) / 2;
	/*
	// Use mean (or position of maximum bin) and standard deviation from initial 1D plots to set the outlier cuts
	//Double_t p_mean = kin_delta_plots[0]->GetMean();
	Double_t p_rms = kin_delta_plots[0]->GetRMS();
	Int_t p_binmax = kin_delta_plots[0]->GetMaximumBin();
	Double_t p_xmax = kin_delta_plots[0]->GetXaxis()->GetBinCenter(p_binmax);
	//Double_t theta_mean = kin_delta_plots[1]->GetMean();
        Double_t theta_rms = kin_delta_plots[1]->GetRMS();
	Int_t theta_binmax = kin_delta_plots[1]->GetMaximumBin();
	Double_t theta_xmax = kin_delta_plots[1]->GetXaxis()->GetBinCenter(theta_binmax);
	//Double_t phi_mean = kin_delta_plots[2]->GetMean();
        Double_t phi_rms = kin_delta_plots[2]->GetRMS();
	Int_t phi_binmax = kin_delta_plots[2]->GetMaximumBin();
	Double_t phi_xmax = kin_delta_plots[2]->GetXaxis()->GetBinCenter(phi_binmax);
	*/

	double p_binWidth = kin_delta_plots[0]->GetXaxis()->GetBinWidth(1);
	double theta_binWidth = kin_delta_plots[1]->GetXaxis()->GetBinWidth(1);
	double phi_binWidth = kin_delta_plots[2]->GetXaxis()->GetBinWidth(1);
	
	// Specify the quantiles (0.25 and 0.75 means 50% of the data lies in between these two points)
	Double_t delta_p_quantiles[2];
	Double_t delta_theta_quantiles[2];
	Double_t delta_phi_quantiles[2];
	Double_t delta_p_peak_quantiles[2];
        Double_t delta_theta_peak_quantiles[2];
        Double_t delta_phi_peak_quantiles[2];
	Double_t probs[2] = {0.25, 0.75};

	// Calculate the quantiles
	kin_delta_plots[0]->GetQuantiles(2, delta_p_quantiles, probs);
	kin_delta_plots[1]->GetQuantiles(2, delta_theta_quantiles, probs);
	kin_delta_plots[2]->GetQuantiles(2, delta_phi_quantiles, probs);

	// Calculate the interquartile range (IQR)
	Double_t delta_p_iqr = delta_p_quantiles[1] - delta_p_quantiles[0];
	Double_t delta_theta_iqr = delta_theta_quantiles[1] - delta_theta_quantiles[0];
	Double_t delta_phi_iqr = delta_phi_quantiles[1] - delta_phi_quantiles[0];

	// Print the results
	//std::cout << "Q1: " << quantiles[0] << std::endl;
	//std::cout << "Q3: " << quantiles[1] << std::endl;
	//std::cout << "IQR: " << iqr << std::endl;

	int bins_in_p_iqr = (int) (delta_p_iqr / p_binWidth + 0.5);
	int p_rebin_num = bins_in_p_iqr / 7;  // The numerator is a factor that dictates the desired number of bins over peak
	int bins_in_theta_iqr = (int) (delta_theta_iqr / theta_binWidth + 0.5);
	int theta_rebin_num = bins_in_theta_iqr / 7;
	int bins_in_phi_iqr = (int) (delta_phi_iqr / phi_binWidth + 0.5);
	int phi_rebin_num = bins_in_phi_iqr / 7;

	if (p_rebin_num > 1) {
	  kin_delta_plots[0]->Rebin(p_rebin_num);
	}
	if (theta_rebin_num > 1) {
          kin_delta_plots[1]->Rebin(theta_rebin_num);
        }
	if (phi_rebin_num > 1) {
          kin_delta_plots[2]->Rebin(theta_rebin_num);
	}

	delta_kin_canvas->cd(1);
	iqr_zoom(kin_delta_plots[0], delta_p_iqr, delta_p_quantiles[0], delta_p_quantiles[1], 5.0);
	kin_delta_plots[0]->Draw("SAME");   
	double p_mean = kin_delta_plots[0]->GetMean();     
	double p_rms = kin_delta_plots[0]->GetRMS();
	
	delta_kin_canvas->cd(2);
	iqr_zoom(kin_delta_plots[1], delta_theta_iqr, delta_theta_quantiles[0], delta_theta_quantiles[1], 5.0);
	kin_delta_plots[1]->Draw("SAME");
        double theta_mean = kin_delta_plots[1]->GetMean();
        double theta_rms = kin_delta_plots[1]->GetRMS();

	delta_kin_canvas->cd(3);
	iqr_zoom(kin_delta_plots[2], delta_phi_iqr, delta_phi_quantiles[0], delta_phi_quantiles[1], 5.0);
	kin_delta_plots[2]->Draw("SAME");
        double phi_mean = kin_delta_plots[2]->GetMean();
        double phi_rms = kin_delta_plots[2]->GetRMS();


	//Get bin centers 
	current_P_bin_center = (P_bin_high + P_bin_low) / 2;
        current_theta_bin_center = (theta_bin_high + theta_bin_low) / 2;
        current_phi_bin_center = (phi_bin_high + phi_bin_low) / 2;

	TPaveText *titleText = new TPaveText(0.1, 0.97, 0.9, 1.0, "NDC");  //x1,y1,x2,y2   
	titleText->SetFillColor(0); // Set background color (0 = white) 
	titleText->SetTextFont(42);
        titleText->SetTextSize(0.03);
        titleText->SetBorderSize(0);
        TString title;
        delta_kin_canvas->cd();
        title.Form("%s: P=%.2f (GeV), #theta=%.2f (deg), #phi=%.2f (deg)", part_type.Data(), current_P_bin_center, current_theta_bin_center, current_phi_bin_center);
        titleText->AddText(title);
        titleText->Draw();
	
	TString deltaKin_histTitle = TString::Format("%s: P=%.2f (GeV), #theta=%.2f (deg), #phi=%.2f (deg)", part_type.Data(), current_P_bin_center, current_theta_bin_center, current_phi_bin_center);
        kin_delta_plots[0]->SetTitle(deltaKin_histTitle.Data());
        kin_delta_plots[1]->SetTitle(deltaKin_histTitle.Data());
        kin_delta_plots[2]->SetTitle(deltaKin_histTitle.Data());
        kin_delta_2D_plots[0]->SetTitle(deltaKin_histTitle.Data());
        kin_delta_2D_plots[1]->SetTitle(deltaKin_histTitle.Data());
        kin_delta_2D_plots[2]->SetTitle(deltaKin_histTitle.Data());

	deltaP_histDir->cd();
        kin_delta_plots[0]->Write();
        deltaTheta_histDir->cd();
        kin_delta_plots[1]->Write();
	deltaPhi_histDir->cd();
	kin_delta_plots[2]->Write();
	
	// P and theta: (rec - mc) 2D histogram
	delta_kin_canvas->cd(2);
        kin_delta_2D_plots[0]->Draw("Colz");
        deltaP_deltaTheta_histDir->cd();
        kin_delta_2D_plots[0]->Write();

	// theta and phi: (rec - mc) 2D histogram
	delta_kin_canvas->cd(4);
        kin_delta_2D_plots[1]->Draw("Colz");
        deltaTheta_deltaPhi_histDir->cd();
        kin_delta_2D_plots[1]->Write();

	// P and phi: (rec - mc) 2D histogram
	delta_kin_canvas->cd(6);
        kin_delta_2D_plots[2]->Draw("Colz");
        deltaP_deltaPhi_histDir->cd();
        kin_delta_2D_plots[2]->Write();

	double rms_multiple = 2;
        double p_low = p_mean - rms_multiple * p_rms;
        double p_upp = p_mean + rms_multiple * p_rms;
        double theta_low = theta_mean - rms_multiple * theta_rms;
        double theta_upp = theta_mean + rms_multiple * theta_rms;
        double phi_low = phi_mean - rms_multiple * phi_rms;
        double phi_upp = phi_mean + rms_multiple * phi_rms;
	
	// ------------------ Calculate covariance with reset ranges--------------------//
	// Create an instance of CovInputParams and initialize its members
	CovInputParams C_P_params(1, 1, p_mean, p_mean, current_P_bin_center, p_diff_this_bin, theta_diff_this_bin, phi_diff_this_bin, p_low, p_upp, theta_low, theta_upp, phi_low, phi_upp);;
        CovInputParams C_theta_params(2, 2, theta_mean, theta_mean, current_P_bin_center, p_diff_this_bin, theta_diff_this_bin, phi_diff_this_bin, p_low, p_upp, theta_low, theta_upp, phi_low, phi_upp);
        CovInputParams C_phi_params(3, 3, phi_mean, phi_mean, current_P_bin_center, p_diff_this_bin, theta_diff_this_bin, phi_diff_this_bin, p_low, p_upp, theta_low, theta_upp, phi_low, phi_upp);
        CovInputParams C_P_phi_params(1, 3, p_mean, phi_mean, current_P_bin_center, p_diff_this_bin, theta_diff_this_bin, phi_diff_this_bin, p_low, p_upp, theta_low, theta_upp, phi_low, phi_upp);
        CovInputParams C_P_theta_params(1, 2, p_mean, theta_mean, current_P_bin_center, p_diff_this_bin, theta_diff_this_bin, phi_diff_this_bin, p_low, p_upp, theta_low, theta_upp, phi_low, phi_upp);
        CovInputParams C_theta_phi_params(2, 3, theta_mean, phi_mean, current_P_bin_center, p_diff_this_bin, theta_diff_this_bin, phi_diff_this_bin, p_low, p_upp, theta_low, theta_upp, phi_low, phi_upp);

	calculateCovariance(variance_P, N_events_C_P, C_P_params);
        calculateCovariance(variance_theta, N_events_C_theta, C_theta_params);
        calculateCovariance(variance_phi, N_events_C_phi, C_phi_params);
        calculateCovariance(covariance_P_phi, N_events_C_P_phi, C_P_phi_params);
        calculateCovariance(covariance_P_theta, N_events_C_P_theta, C_P_theta_params);
        calculateCovariance(covariance_theta_phi, N_events_C_theta_phi, C_theta_phi_params);

	N_good_events = N_events_C_P;  // This should be the same for all covariance elements

      	// Calculate errors on covariance values      
	variance_P_err = sqrt(2.0 / (N_good_events - 1)) * variance_P; 
	variance_theta_err = sqrt(2.0 / (N_good_events - 1)) * variance_theta;      
	variance_phi_err = sqrt(2.0 / (N_good_events - 1)) * variance_phi; 
	covariance_P_theta_err = sqrt((covariance_P_theta * covariance_P_theta + variance_P * variance_theta) / (N_good_events - 1));  
	covariance_P_phi_err = sqrt((covariance_P_phi * covariance_P_phi + variance_P * variance_phi) / (N_good_events - 1)); 
	covariance_theta_phi_err = sqrt((covariance_theta_phi * covariance_theta_phi + variance_phi * variance_theta) / (N_good_events - 1));

	// Get means (offsets)
	p_diff_offset = kin_delta_plots[0]->GetMean();
	theta_diff_offset = kin_delta_plots[1]->GetMean();
	phi_diff_offset = kin_delta_plots[2]->GetMean();
	// ------------------------------------------------------------------------------------- // 
	

	/*
	// ============ Calculate covariance using TRobustEstimator ===============//
	CovarianceResult result = calculateCovarianceRobust(p_diff_this_bin, theta_diff_this_bin, phi_diff_this_bin, p_low, p_upp, theta_low, theta_upp, phi_low, phi_upp);
	try {
	  const TMatrixDSym* covMatrix = result.covariance;
	  const TVectorD* means = result.means;
	  int N_good_events = result.good_counts;
	  //std::cout << "counts: " << N_good_events << std::endl;	  
	  variance_P = (*covMatrix)(0, 0);
	  variance_theta = (*covMatrix)(1, 1);
	  variance_phi = (*covMatrix)(2, 2);
	  covariance_P_theta = (*covMatrix)(0, 1);
	  covariance_P_phi = (*covMatrix)(0, 2);
	  covariance_theta_phi = (*covMatrix)(1, 2);

	  // Calculate errors on covariance values
	  variance_P_err = sqrt(2.0 / (N_good_events - 1)) * variance_P;
	  variance_theta_err = sqrt(2.0 / (N_good_events - 1)) * variance_theta;
	  variance_phi_err = sqrt(2.0 / (N_good_events - 1)) * variance_phi;
	  covariance_P_theta_err = sqrt((covariance_P_theta * covariance_P_theta + variance_P * variance_theta) / (N_good_events - 1));
	  covariance_P_phi_err = sqrt((covariance_P_phi * covariance_P_phi + variance_P * variance_phi) / (N_good_events - 1));
	  covariance_theta_phi_err = sqrt((covariance_theta_phi * covariance_theta_phi + variance_phi * variance_theta) / (N_good_events - 1));

	  p_diff_offset = (*means)[0];
	  theta_diff_offset= (*means)[1];
	  phi_diff_offset= (*means)[2];

	} catch (const std::exception& e) {
	  std::cerr << "Covariance matrix computation failed: " << e.what() << std::endl;
	  std::cout << "counts: " << result.good_counts << std::endl;
	}
	// ======================================================================//
	*/


	//Get the number of entries in each bin
	//N_entries = p_diff_binned[pBin][thetaBin][phiBin].size();
	//InfoFile << "N_good_events = " << N_good_events << std::endl;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Save covariance matrix elements in C_tree and TH3's~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	//Get bin centers
        current_P_bin_center = (P_bin_high + P_bin_low) / 2;  
	current_theta_bin_center = (theta_bin_high + theta_bin_low) / 2;
	current_phi_bin_center = (phi_bin_high + phi_bin_low) / 2;   

	InfoFile << "current_P_bin_center = " << current_P_bin_center << ", current_theta_bin_center = " << current_theta_bin_center << ", current_phi_bin_center = " << current_phi_bin_center << std::endl;
	
	//Fill the tree with this bin's covariance values
	P = current_P_bin_center;
	theta = current_theta_bin_center;
	phi = current_phi_bin_center;
	C_P = variance_P;
	C_theta = variance_theta;
	C_phi = variance_phi;
	C_P_phi = covariance_P_phi;
	C_P_theta = covariance_P_theta;
	C_theta_phi = covariance_theta_phi;
	C_P_err = variance_P_err;
        C_theta_err = variance_theta_err;
        C_phi_err = variance_phi_err;
        C_P_phi_err = covariance_P_phi_err;
        C_P_theta_err = covariance_P_theta_err;
        C_theta_phi_err = covariance_theta_phi_err;
	event_count = N_good_events;
	P_offset = p_diff_offset;
        theta_offset = theta_diff_offset;
        phi_offset = phi_diff_offset;
	C_tree->Fill();
	
	//Fill TH3's by setting the bin content of each bin to be equal to the covariance value
	C_P_hist->Fill(current_P_bin_center, current_theta_bin_center, current_phi_bin_center, variance_P);
	C_theta_hist->Fill(current_P_bin_center, current_theta_bin_center, current_phi_bin_center, variance_theta);
	C_phi_hist->Fill(current_P_bin_center, current_theta_bin_center, current_phi_bin_center, variance_phi);
	C_P_phi_hist->Fill(current_P_bin_center, current_theta_bin_center, current_phi_bin_center, covariance_P_phi);
	C_P_theta_hist->Fill(current_P_bin_center, current_theta_bin_center, current_phi_bin_center, covariance_P_theta);
	C_theta_phi_hist->Fill(current_P_bin_center, current_theta_bin_center, current_phi_bin_center, covariance_theta_phi);
	C_P_err_hist->Fill(current_P_bin_center, current_theta_bin_center, current_phi_bin_center, variance_P_err);
        C_theta_err_hist->Fill(current_P_bin_center, current_theta_bin_center, current_phi_bin_center, variance_theta_err);
        C_phi_err_hist->Fill(current_P_bin_center, current_theta_bin_center, current_phi_bin_center, variance_phi_err);
        C_P_phi_err_hist->Fill(current_P_bin_center, current_theta_bin_center, current_phi_bin_center, covariance_P_phi_err);
        C_P_theta_err_hist->Fill(current_P_bin_center, current_theta_bin_center, current_phi_bin_center, covariance_P_theta_err);
        C_theta_phi_err_hist->Fill(current_P_bin_center, current_theta_bin_center, current_phi_bin_center, covariance_theta_phi_err);
	entries_hist->Fill(current_P_bin_center, current_theta_bin_center, current_phi_bin_center, N_good_events);
	P_offset_hist->Fill(current_P_bin_center, current_theta_bin_center, current_phi_bin_center, p_diff_offset);
        theta_offset_hist->Fill(current_P_bin_center, current_theta_bin_center, current_phi_bin_center, theta_diff_offset);
        phi_offset_hist->Fill(current_P_bin_center, current_theta_bin_center, current_phi_bin_center, phi_diff_offset);
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	
	
	for (TH1* hist : kin_delta_plots) {
	  delete hist;
	}

	for (TH2* hist2D : kin_delta_2D_plots) {
	  delete hist2D;
	}
	kin_delta_plots.clear();
	kin_delta_2D_plots.clear();
      }
    }
  } 
  
  //delta_kin_canvas->Print(delta_kin_pdf_name + "]");
  //std::cout << "Closing rec_mc_diff_plots.pdf file." << std::endl;
  delete delta_kin_canvas;

  // Close the (rec-mc) histogram root file
  delta_kin_root_file->Close();

  // Change to the covariance root file
  C_file->cd();
  // Save the TH3's to the ROOT file
  C_P_hist->Write();
  C_theta_hist->Write();
  C_phi_hist->Write();
  C_P_phi_hist->Write();
  C_P_theta_hist->Write();
  C_theta_phi_hist->Write();
  C_P_err_hist->Write();
  C_theta_err_hist->Write();
  C_phi_err_hist->Write();
  C_P_phi_err_hist->Write();
  C_P_theta_err_hist->Write();
  C_theta_phi_err_hist->Write();
  entries_hist->Write();
  P_offset_hist->Write();
  theta_offset_hist->Write();
  phi_offset_hist->Write();

  // Save the tree to the ROOT file
  C_tree->Write();

  // Close the ROOT file
  C_file->Close();
  
  return 0;
}


void read_Hipo(char inputFile[256], int part_charge, float part_mass, partInfoStruct &partInfo, int &event_count, binnedKinInfo &binnedKins, TH3* P_theta_and_phi)
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
  //hipo::bank RECCAL(factory.getSchema("REC::Calorimeter"));
  //hipo::bank RECTRACK(factory.getSchema("REC::Track"));
  hipo::bank TBTRACK(factory.getSchema("TimeBasedTrkg::TBTracks"));
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
    //event.getStructure(RECCAL);
    //event.getStructure(RECTRACK);
    event.getStructure(TBTRACK);
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
    read_Rec_Part_Bank(RECPART, MC_MATCH, TBTRACK, &rec_pid, &rec_px, &rec_py, &rec_pz, &rec_sector, &chi2, &ndf, event_count);

    //Get generated particle info
    std::vector<int> mc_pid;
    std::vector<float> mc_px;
    std::vector<float> mc_py;
    std::vector<float> mc_pz;
    read_MC_Part_Bank(MCPART, &mc_pid, &mc_px, &mc_py, &mc_pz);
    
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

	Int_t P_bin_num;
	Int_t theta_bin_num;
	Int_t phi_bin_num;
	P_bin_num =  P_theta_and_phi->GetXaxis()->FindBin(rec_p);
	theta_bin_num =  P_theta_and_phi->GetYaxis()->FindBin(rec_theta_deg);
	phi_bin_num =  P_theta_and_phi->GetZaxis()->FindBin(rec_phi_deg);

	//If the conditions in the if-statement below are not met, it means that at least one of the reconstructed kinematics is outside the allowed range 
	if (((P_bin_num > 0) && (P_bin_num <= P_bins)) && ((theta_bin_num > 0) && (theta_bin_num <= theta_bins)) && ((phi_bin_num > 0) && (phi_bin_num <= phi_bins))) {

	  //Save the kinematic variable info for this bin to be used in the calculateCovariance method defined at the beginning of this file 
	  binnedKins.p_diff_binned[P_bin_num-1][theta_bin_num-1][phi_bin_num-1].push_back(p_diff);
	  binnedKins.theta_diff_binned[P_bin_num-1][theta_bin_num-1][phi_bin_num-1].push_back(theta_diff_deg);
	  binnedKins.phi_diff_binned[P_bin_num-1][theta_bin_num-1][phi_bin_num-1].push_back(phi_diff_deg);
	  binnedKins.rec_p_binned[P_bin_num-1][theta_bin_num-1][phi_bin_num-1].push_back(rec_p);
	  binnedKins.deltaP_over_p_binned[P_bin_num-1][theta_bin_num-1][phi_bin_num-1].push_back(deltaP_over_p);

	  partInfo.rec_p_vec.push_back(rec_p);
	  partInfo.mc_p_vec.push_back(mc_p);
	  partInfo.p_diff_vec.push_back(p_diff);
	  partInfo.chi2_vec.push_back(chi2[0]);
	  //partInfo.ndf_vec.push_back(ndf[0]);
	  partInfo.theta_diff_deg_vec.push_back(theta_diff_deg);
	  partInfo.phi_diff_deg_vec.push_back(phi_diff_deg);
	  partInfo.rec_theta_deg_vec.push_back(rec_theta_deg);
	  partInfo.mc_theta_deg_vec.push_back(mc_theta_deg);
	  partInfo.rec_phi_deg_vec.push_back(rec_phi_deg);
	  partInfo.mc_phi_deg_vec.push_back(mc_phi_deg);
	}
	else {
	  rejected_events++;
	}
    }      
  }
  //std::cout << rejected_events << " rejected events" << std::endl;
}
//This is intended to be used for reading the MC info
void read_MC_Part_Bank(hipo::bank PartBank, std::vector<int>* pid, std::vector<float>* px, std::vector<float>* py, std::vector<float>* pz){  
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
void read_Rec_Part_Bank(hipo::bank PartBank, hipo::bank MCMatchBank, hipo::bank TBTrack, std::vector<int>* pid, std::vector<float>* px, std::vector<float>* py, std::vector<float>* pz, std::vector<int> *rec_sector, std::vector<float> *chi2, std::vector<int> *ndf, int &event_count){

  // Fail safe for truth matching
  if (MCMatchBank.getRows() < 1)
    return;

  // This assumes only 1 particle in MC::GenMatch bank 
  int index_truth = MCMatchBank.getInt("pindex", 0);

  if (index_truth < 0)
    return;


  event_count++; 
  int nrows = PartBank.getRows();
  //int calBank_rows = CalBank.getRows();
  //int rectrack_rows = RecTrack.getRows();
  int tbtrack_rows = TBTrack.getRows();
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
  int clusters_recorded = 0;


  int   current_pid = PartBank.getInt("pid",index_truth);
  float  current_px = PartBank.getFloat("px",index_truth);
  float  current_py = PartBank.getFloat("py",index_truth);
  float  current_pz = PartBank.getFloat("pz",index_truth);
  float  current_chi2pid = PartBank.getFloat("chi2pid",index_truth);
  int current_status = PartBank.getInt("status",index_truth);
  int charge = PartBank.getInt("charge",index_truth);
  int current_pindex = -2;


  //Save particle info
  px->push_back(px_vec[0]);      
  py->push_back(py_vec[0]);   
  pz->push_back(pz_vec[0]);   
  pid->push_back(pid_vec[0]);
  chi2->push_back(chi2_vec[0]);
  ndf->push_back(ndf_vec[0]);
  
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

