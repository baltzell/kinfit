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
#include <numeric>
#include <random>


void read_Hipo(char in_data[256], int part_charge, float part_mass, std::vector<float>* mc_p_vec, std::vector<float>* rec_p_vec, std::vector<float>* mc_theta_deg_vec, std::vector<float>* rec_theta_deg_vec, std::vector<float>* mc_phi_deg_vec, std::vector<float>* rec_phi_deg_vec, std::vector<float>* p_diff_vec, std::vector<float>* theta_diff_rad_vec, std::vector<float>* phi_diff_rad_vec, std::vector<float>* chi2_vec, std::vector<int> *ndf_vec, int& event_count);
void read_Part_Bank(hipo::bank PartBank, std::vector<int>* pid, std::vector<float>* px, std::vector<float>* py, std::vector<float>* pz);
void read_Part_Bank(hipo::bank PartBank, hipo::bank CalBank, hipo::bank RecTrack, std::vector<int>* pid, std::vector<float>* px, std::vector<float>* py, std::vector<float>* pz, std::vector<int> *rec_status, std::vector<float> *chi2, std::vector<int> *ndf, int& event_count);

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
//Use Sturge's Rule to determine the number of bins to be used in histogram based on number of entries
int Sturges_rule(Int_t entries) {
  int N_bins = 1 + 3.322*TMath::Log10(entries);
  return N_bins;
}

double calculateCovariance(const std::vector<double>& x, const std::vector<double>& y, double x_low, double x_upp, double y_low, double y_upp) {
  Int_t numValues = x.size(); // Get the number of values in the arrays                                                                                                                                      
  double covariance = 0;
  int good_counts = 0;
  double x_tot = 0;
  double y_tot = 0;
  std::vector<double> x_good;
  std::vector<double> y_good;
  for (int j = 0; j < numValues; j++) {
    if ((x[j] > x_low && x[j] < x_upp) && (y[j] > y_low && y[j] < y_upp)) {
      x_tot += x[j];
      y_tot += y[j];
      good_counts++;
      x_good.push_back(x[j]);
      y_good.push_back(y[j]);
    }
  }
  double x_mean = x_tot / good_counts;
  double y_mean = y_tot / good_counts;
  for (Int_t i = 0; i < good_counts; ++i) {
     covariance += (x_good[i] - x_mean) * (y_good[i] - y_mean);
  }
  covariance /= (good_counts - 1);
  return covariance;
}

double calculateCovarianceError(const std::vector<double>& x, const std::vector<double>& y, int N_bootstrapSamples, double x_low, double x_upp, double y_low, double y_upp) {
  double observedCovariance = calculateCovariance(x, y, x_low, x_upp, y_low, y_upp);

  // Perform bootstrapping   
  int N_events = x.size();
  std::vector<double> bootstrapCovariances;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> dist(0, N_events - 1);
  for (int b = 0; b < N_bootstrapSamples; ++b) {
    std::vector<double> bootstrapX;
    std::vector<double> bootstrapY;
    // Generate a bootstrap sample by resampling with replacement         
    for (int i = 0; i < N_events; ++i) {
      int randomIndex = dist(gen);
      bootstrapX.push_back(x[randomIndex]);
      bootstrapY.push_back(y[randomIndex]);
    }
    // Calculate the covariance of the bootstrap sample  
    double bootstrapCovariance = calculateCovariance(bootstrapX, bootstrapY, x_low, x_upp, y_low, y_upp);
    bootstrapCovariances.push_back(bootstrapCovariance);
  }
  // Calculate the variance of the bootstrapped covariances 
  double bootstrapVariance = 0.0;
  for (double covariance : bootstrapCovariances) {
    double diff = covariance - observedCovariance;
    bootstrapVariance += (diff * diff);
  }
  bootstrapVariance /= (N_bootstrapSamples - 1);
  // Estimate the error (uncertainty) of the covariance       
  double covarianceError = std::sqrt(bootstrapVariance);
  return covarianceError;
}

//Output file for saving any information for trouble shooting 
std::ofstream InfoFile("event_info.txt");

//Initialize the particle type and binning info
const TString part_type = "#pi^{+}";
const float part_mass = 0.13957;
const int part_charge = 1;
const int P_bins = 2;
const int theta_bins = 2;
const int phi_bins = 2;

//Set boundaries for the kinematics. This is necessaery to exclude outliers
const double p_low = -1;
const double p_upp = 1;
const double theta_low = -0.05;
const double theta_upp = 0.05;
const double phi_low = -0.1;
const double phi_upp = 0.1;

//These variables are used for the plots that combine all P, theta, and phi bins
const std::vector<TString> KINES = {"P", "#theta", "#phi"};
const std::vector<TString> UNITS = {"[GeV]", "[rad]", "[rad]"};
const std::vector<int> kin_plots_low = {0, 0, 25};
const std::vector<int> kin_plots_high = {12, 40, 95};
const std::vector<double> kin_delta_plots_low = {p_low, theta_low, phi_low};
const std::vector<double> kin_delta_plots_high = {p_upp, theta_upp, phi_upp};

//Use these values to set a maximum on the number of files that the reader loop will read and to set a number of elements for the vectors to be initilaized to
const int max_files = 10;
const int N_events_per_file = 10000;
const int N_events = max_files*N_events_per_file;

int covMatrix_extraction()
{
  //---------------Specify input file or directory path---------------//                                                                                                                                  
  //char in_data[256] = "/volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/pip-sec2/cooked/out_pip-sec22.rec.hipo";
  //char in_data[256] = "/volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/pip-sec2-/cooked/"; 
  //char in_data[256] = "/volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/test_dir/";
  char in_data[256] = "/volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/pip-sec2-100mil_events/cooked/";
  //------------------------------------------------------------------//                                                                                                                                   


  //Declare vectors for storing kinematic info                                                                                                                                                             
  std::vector<float> p_diff_vec;
  std::vector<float> theta_diff_rad_vec;
  std::vector<float> phi_diff_rad_vec;
  std::vector<float> mc_p_vec;
  std::vector<float> rec_p_vec;
  std::vector<float> mc_theta_deg_vec;
  std::vector<float> rec_theta_deg_vec;
  std::vector<float> mc_phi_deg_vec;
  std::vector<float> rec_phi_deg_vec;
  std::vector<float> chi2_vec;
  std::vector<int> ndf_vec;
  int event_count = 0;

  p_diff_vec.reserve(N_events);
  phi_diff_rad_vec.reserve(N_events);
  mc_p_vec.reserve(N_events);
  rec_p_vec.reserve(N_events);
  mc_theta_deg_vec.reserve(N_events);
  rec_theta_deg_vec.reserve(N_events);
  mc_phi_deg_vec.reserve(N_events);
  rec_phi_deg_vec.reserve(N_events);
  chi2_vec.reserve(N_events);
  ndf_vec.reserve(N_events);

  DIR *dr;
  struct dirent *en;
  dr = opendir(in_data);                                                                                                                                                
  //If in_data is a directory, loop through all files within it  
  int in_file_count = 0;
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
        read_Hipo(dir_path_char, part_charge, part_mass, &mc_p_vec, &rec_p_vec, &mc_theta_deg_vec, &rec_theta_deg_vec, &mc_phi_deg_vec, &rec_phi_deg_vec, &p_diff_vec, &theta_diff_rad_vec, &phi_diff_rad_vec, &chi2_vec, &ndf_vec, event_count);
        //Print list of all read hipo files to input_files.txt   
        //InFileList << in_data << en->d_name << std::endl;
      }
    }
    closedir(dr); //close all directory
    std::cout << in_file_count << " files read" << std::endl;
    std::cout << p_diff_vec.size() << " events saved out of " << event_count << " events read." << std::endl;
  }
  //If in_data is not a directory, then get the single input file    
  else {
    std::cout << "Not a directory. Only single input file: " << in_data << std::endl;
    read_Hipo(in_data, part_charge, part_mass, &mc_p_vec, &rec_p_vec, &mc_theta_deg_vec, &rec_theta_deg_vec, &mc_phi_deg_vec, &rec_phi_deg_vec, &p_diff_vec, &theta_diff_rad_vec, &phi_diff_rad_vec, &chi2_vec, &ndf_vec, event_count);
  }



  //---------------------------------------------Plots Over All Bins-------------------------------------------------//
  auto chi2_ndf_can = new TCanvas("chi2_ndf_can", "chi2_ndf_can", 800, 800);
  chi2_ndf_can->Divide(2,3);
  auto chi2_can = new TCanvas("chi2_can", "chi2_can", 800, 800);
  chi2_can->Divide(2,3);
  auto ndf_can = new TCanvas("ndf_can", "ndf_can", 800, 800);
  ndf_can->Divide(2,3);
  //Create Histograms
  std::vector<TH2*> chi2_ndf_vs_rec_kin_plots;
  std::vector<TH2*> chi2_ndf_vs_delta_kin_plots;
  std::vector<TH2*> chi2_vs_rec_kin_plots;
  std::vector<TH2*> chi2_vs_delta_kin_plots;
  std::vector<TH2*> ndf_vs_rec_kin_plots;
  std::vector<TH2*> ndf_vs_delta_kin_plots;
  for (int j = 0; j < KINES.size(); j++) {
    TH2* chi2_ndf_vs_rec_kin_temp = new TH2D(Form("chi2_ndf_%s_plot", KINES[j].Data()), Form("%s;Rec %s %s;#Chi^{2}/NDF",  part_type.Data(), KINES[j].Data(), UNITS[j].Data()), 200, kin_plots_low[j], kin_plots_high[j], 200, 0, 10);
    chi2_ndf_vs_rec_kin_temp->GetXaxis()->SetTitleSize(0.07);
    chi2_ndf_vs_rec_kin_temp->GetYaxis()->SetTitleSize(0.07);
    chi2_ndf_vs_rec_kin_plots.push_back(chi2_ndf_vs_rec_kin_temp);

    TH2* chi2_ndf_vs_delta_kin_temp = new TH2D(Form("chi2_ndf_delta_%s_plot", KINES[j].Data()), Form("%s;#Delta %s (Rec - Gen) %s;#Chi^{2}/NDF", part_type.Data(), KINES[j].Data(), UNITS[j].Data()), 200, kin_delta_plots_low[j], kin_delta_plots_high[j], 200, 0, 10);
    chi2_ndf_vs_delta_kin_temp->GetXaxis()->SetTitleSize(0.07);
    chi2_ndf_vs_delta_kin_temp->GetYaxis()->SetTitleSize(0.07);
    chi2_ndf_vs_delta_kin_plots.push_back(chi2_ndf_vs_delta_kin_temp);

    TH2* chi2_vs_rec_kin_temp = new TH2D(Form("chi2_%s_plot", KINES[j].Data()), Form("%s;Rec %s %s;#Chi^{2}",  part_type.Data(), KINES[j].Data(), UNITS[j].Data()), 200, kin_plots_low[j], kin_plots_high[j], 200, 0, 2000);
    chi2_vs_rec_kin_temp->GetXaxis()->SetTitleSize(0.07);
    chi2_vs_rec_kin_temp->GetYaxis()->SetTitleSize(0.07);
    chi2_vs_rec_kin_plots.push_back(chi2_vs_rec_kin_temp);

    TH2* chi2_vs_delta_kin_temp = new TH2D(Form("chi2_delta_%s_plot", KINES[j].Data()), Form("%s;#Delta %s (Rec - Gen) %s;#Chi^{2}", part_type.Data(), KINES[j].Data(), UNITS[j].Data()), 200, kin_delta_plots_low[j], kin_delta_plots_high[j], 200, 0, 2000);
    chi2_vs_delta_kin_temp->GetXaxis()->SetTitleSize(0.07);
    chi2_vs_delta_kin_temp->GetYaxis()->SetTitleSize(0.07);
    chi2_vs_delta_kin_plots.push_back(chi2_vs_delta_kin_temp);

    TH2* ndf_vs_rec_kin_temp = new TH2D(Form("ndf_%s_plot", KINES[j].Data()), Form("%s;Rec %s %s;NDF",  part_type.Data(), KINES[j].Data(), UNITS[j].Data()), 200, kin_plots_low[j], kin_plots_high[j], 50, 0, 50);
    ndf_vs_rec_kin_temp->GetXaxis()->SetTitleSize(0.07);
    ndf_vs_rec_kin_temp->GetYaxis()->SetTitleSize(0.07);
    ndf_vs_rec_kin_plots.push_back(ndf_vs_rec_kin_temp);

    TH2* ndf_vs_delta_kin_temp = new TH2D(Form("ndf_delta_%s_plot", KINES[j].Data()), Form("%s;#Delta %s (Rec - Gen) %s;NDF", part_type.Data(), KINES[j].Data(), UNITS[j].Data()), 200, kin_delta_plots_low[j], kin_delta_plots_high[j], 50, 0, 50);
    ndf_vs_delta_kin_temp->GetXaxis()->SetTitleSize(0.07);
    ndf_vs_delta_kin_temp->GetYaxis()->SetTitleSize(0.07);
    ndf_vs_delta_kin_plots.push_back(ndf_vs_delta_kin_temp);
  }
  //Fill Histograms
  for (int i = 0; i < chi2_vec.size(); i++) {
    chi2_ndf_vs_rec_kin_plots[0]->Fill(rec_p_vec[i], chi2_vec[i]/ndf_vec[i]);
    chi2_ndf_vs_rec_kin_plots[1]->Fill(rec_theta_deg_vec[i], chi2_vec[i]/ndf_vec[i]);
    chi2_ndf_vs_rec_kin_plots[2]->Fill(rec_phi_deg_vec[i], chi2_vec[i]/ndf_vec[i]);
    chi2_ndf_vs_delta_kin_plots[0]->Fill(p_diff_vec[i], chi2_vec[i]/ndf_vec[i]);
    chi2_ndf_vs_delta_kin_plots[1]->Fill(theta_diff_rad_vec[i], chi2_vec[i]/ndf_vec[i]);
    chi2_ndf_vs_delta_kin_plots[2]->Fill(phi_diff_rad_vec[i], chi2_vec[i]/ndf_vec[i]);
    
    chi2_vs_rec_kin_plots[0]->Fill(rec_p_vec[i], chi2_vec[i]);
    chi2_vs_rec_kin_plots[1]->Fill(rec_theta_deg_vec[i], chi2_vec[i]);
    chi2_vs_rec_kin_plots[2]->Fill(rec_phi_deg_vec[i], chi2_vec[i]);
    chi2_vs_delta_kin_plots[0]->Fill(p_diff_vec[i], chi2_vec[i]);
    chi2_vs_delta_kin_plots[1]->Fill(theta_diff_rad_vec[i], chi2_vec[i]);
    chi2_vs_delta_kin_plots[2]->Fill(phi_diff_rad_vec[i], chi2_vec[i]);
    
    ndf_vs_rec_kin_plots[0]->Fill(rec_p_vec[i], ndf_vec[i]);
    ndf_vs_rec_kin_plots[1]->Fill(rec_theta_deg_vec[i], ndf_vec[i]);
    ndf_vs_rec_kin_plots[2]->Fill(rec_phi_deg_vec[i], ndf_vec[i]);
    ndf_vs_delta_kin_plots[0]->Fill(p_diff_vec[i], ndf_vec[i]);
    ndf_vs_delta_kin_plots[1]->Fill(theta_diff_rad_vec[i], ndf_vec[i]);
    ndf_vs_delta_kin_plots[2]->Fill(phi_diff_rad_vec[i], ndf_vec[i]);
  }
  //Draw Histograms
  for (int j = 0; j < KINES.size(); j++) {
    chi2_ndf_can->cd(2*j+1);                                                                                                                                                                   
    gPad->SetBottomMargin(0.15);    //So the axis titles fit on canvas                                                                                                                       
    gPad->SetLeftMargin(0.15);                                                                                                                                                               
    chi2_ndf_vs_rec_kin_plots[j]->Draw("Colz");
    chi2_ndf_can->cd(2*j+2);
    gPad->SetBottomMargin(0.15);    //So the axis titles fit on canvas  
    gPad->SetLeftMargin(0.15);
    chi2_ndf_vs_delta_kin_plots[j]->Draw("Colz");
    
    chi2_can->cd(2*j+1);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);
    chi2_vs_rec_kin_plots[j]->Draw("Colz");
    chi2_can->cd(2*j+2);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);
    chi2_vs_delta_kin_plots[j]->Draw("Colz");
    
    ndf_can->cd(2*j+1);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);
    ndf_vs_rec_kin_plots[j]->Draw("Colz");
    ndf_can->cd(2*j+2);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);
    ndf_vs_delta_kin_plots[j]->Draw("Colz");
  }
  //chi2_ndf_can->SaveAs("plots/Chi2_ndf_vs_kins.pdf");
  //chi2_can->SaveAs("plots/Chi2_vs_kins.pdf");
  //ndf_can->SaveAs("plots/NDF_vs_kins.pdf");
  //-----------------------------------------------------------------------------------------------------------------//



  //Particles are generated with an energy range of 1 to 11 GeV, theta from 2 to 38 deg, and phi from 25 to 95 deg   
  double P_bin_min = 1;
  double P_bin_max = 11;
  double P_bin_size = (P_bin_max - P_bin_min)/P_bins;
  double theta_bin_min = 5;
  double theta_bin_max = 35;
  double theta_bin_size = (theta_bin_max - theta_bin_min) / theta_bins;
  double phi_bin_min = 30;
  double phi_bin_max = 90;
  double phi_bin_size = (phi_bin_max - phi_bin_min) / phi_bins;

  //3-dimensional vector for saving bin info. From innermost outward: phi bin, theta bin, and then P bin
  std::vector<std::vector<std::vector<double>>> P_bin_centers_all(P_bins, std::vector<std::vector<double>> (theta_bins, std::vector<double> (phi_bins)));
  std::vector<std::vector<std::vector<double>>> theta_bin_centers_all(P_bins, std::vector<std::vector<double>> (theta_bins, std::vector<double> (phi_bins)));
  std::vector<std::vector<std::vector<double>>> phi_bin_centers_all(P_bins, std::vector<std::vector<double>> (theta_bins, std::vector<double> (phi_bins)));

  //4-dimensional vector for saving (Rec - MC) info for each kinematic variable. From innermost outward: (Rec - MC) for given kinematic, phi bin, theta bin, and then P bin
  std::vector<std::vector<std::vector<std::vector<double>>>> p_diff_binned(P_bins, std::vector<std::vector<std::vector<double>>> (theta_bins, std::vector<std::vector<double>> (phi_bins, std::vector<double> ())));
  std::vector<std::vector<std::vector<std::vector<double>>>> phi_diff_binned(P_bins, std::vector<std::vector<std::vector<double>>> (theta_bins, std::vector<std::vector<double>> (phi_bins, std::vector<double> ())));
  std::vector<std::vector<std::vector<std::vector<double>>>> theta_diff_binned(P_bins, std::vector<std::vector<std::vector<double>>> (theta_bins, std::vector<std::vector<double>> (phi_bins, std::vector<double> ())));

  // Create a ROOT file to store the data
  TFile* C_file = new TFile("covariances.root", "RECREATE");

  // Create a TTree with four branches for the four dimensions
  TTree* C_tree = new TTree("covarianceTree", "Covariance Tree");
  double P, theta, phi, C_P, C_P_err, C_theta, C_theta_err, C_phi, C_phi_err, C_P_phi, C_P_phi_err, C_P_theta, C_P_theta_err, C_theta_phi, C_theta_phi_err, entries;
  C_tree->Branch("P", &P, "P/D");
  C_tree->Branch("theta", &theta, "theta/D");
  C_tree->Branch("phi", &phi, "phi/D");
  C_tree->Branch("C_P", &C_P, "C_P/D");
  C_tree->Branch("C_theta", &C_theta, "C_theta/D");
  C_tree->Branch("C_phi", &C_phi, "C_phi/D");
  C_tree->Branch("C_P_phi", &C_P_phi, "C_P_phi/D");
  C_tree->Branch("C_P_theta", &C_P_theta, "C_P_theta/D");
  C_tree->Branch("C_theta_phi", &C_theta_phi, "C_theta_phi/D");
  C_tree->Branch("event_count", &event_count, "event_count/I");

  //Create a TH3 for each covariance parameter and each covariance error parameter
  TH3D* C_P_hist = new TH3D("C_P_hist", "C_P Histogram", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);
  TH3D* C_phi_hist = new TH3D("C_phi_hist", "C_phi Histogram", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);
  TH3D* C_theta_hist = new TH3D("C_theta_hist", "C theta_Hist", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);
  TH3D* C_P_phi_hist = new TH3D("C_P_phi_hist", "C_P_phi Hist", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);
  TH3D* C_P_theta_hist = new TH3D("C_P_theta_hist", "C_P_theta Hist", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);
  TH3D* C_theta_phi_hist = new TH3D("C_theta_phi_hist", "C_theta_phi Hist", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);
  TH3D* entries_hist = new TH3D("entries_hist", "Entries Histogram", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);


  //Set up TH3 to store the binning info
  TH3* P_theta_and_phi = new TH3D("P_theta_and_phi", "P_theta_and_phi", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);  

  //Loop through events saved from the input hipo file(s)
  for (int i = 0; i < mc_p_vec.size(); i++) {
    Int_t P_bin_num =  P_theta_and_phi->GetXaxis()->FindBin(rec_p_vec[i]);
    Int_t theta_bin_num =  P_theta_and_phi->GetYaxis()->FindBin(rec_theta_deg_vec[i]);
    Int_t phi_bin_num =  P_theta_and_phi->GetZaxis()->FindBin(rec_phi_deg_vec[i]);
    //InfoFile << "P_bin_num = " << P_bin_num << ", theta_bin_num = " << theta_bin_num << ", phi_bin_num = " << phi_bin_num << std::endl;
    
    if (((P_bin_num > 0) && (P_bin_num <= P_bins)) && ((theta_bin_num > 0) && (theta_bin_num <= theta_bins)) && ((phi_bin_num > 0) && (phi_bin_num <= phi_bins))) { 
      //Save the kinematic variable info for this bin to be used in the calculateCovariance method defined at the beginning of this file
      p_diff_binned[P_bin_num-1][theta_bin_num-1][phi_bin_num-1].push_back(p_diff_vec[i]);
      phi_diff_binned[P_bin_num-1][theta_bin_num-1][phi_bin_num-1].push_back(phi_diff_rad_vec[i]);
      theta_diff_binned[P_bin_num-1][theta_bin_num-1][phi_bin_num-1].push_back(theta_diff_rad_vec[i]);
    }
  }

  //Loop through the bins and calculate the covariances
  for (int pBin = 0; pBin < P_bins; pBin++) {
    double P_bin_low = P_bin_min + pBin*P_bin_size;  
    double P_bin_high = P_bin_low + P_bin_size;
    for (int thetaBin = 0; thetaBin < theta_bins; thetaBin++) {
      double theta_bin_low = theta_bin_min + thetaBin*theta_bin_size; 
      double theta_bin_high = theta_bin_low + theta_bin_size;
      for (int phiBin = 0; phiBin < phi_bins; phiBin++) {
	double phi_bin_low = phi_bin_min + phiBin*phi_bin_size;
	double phi_bin_high = phi_bin_low + phi_bin_size;

	//Get variances and covariances of the three kinematic variables: P, theta, and phi 
	std::vector<double> p_diff_this_bin = p_diff_binned[pBin][thetaBin][phiBin];
	std::vector<double> theta_diff_this_bin = theta_diff_binned[pBin][thetaBin][phiBin];
	std::vector<double> phi_diff_this_bin = phi_diff_binned[pBin][thetaBin][phiBin];

	double variance_P = calculateCovariance(p_diff_this_bin, p_diff_this_bin, p_low, p_upp, p_low, p_upp);
	double variance_theta = calculateCovariance(theta_diff_this_bin, theta_diff_this_bin, theta_low, theta_upp, theta_low, theta_upp);
        double variance_phi = calculateCovariance(phi_diff_this_bin, phi_diff_this_bin, phi_low, phi_upp, phi_low, phi_upp);
	double covariance_P_phi = calculateCovariance(p_diff_this_bin, phi_diff_this_bin, p_low, p_upp, phi_low, phi_upp);
	double covariance_P_theta = calculateCovariance(p_diff_this_bin, theta_diff_this_bin, p_low, p_upp, theta_low, theta_upp);
	double covariance_theta_phi = calculateCovariance(phi_diff_this_bin, theta_diff_this_bin, phi_low, phi_upp, theta_low, theta_upp);

	//Get the number of entries in each bin
	int N_entries = p_diff_binned[pBin][thetaBin][phiBin].size();
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Save covariance matrix elements in C_tree and TH3's~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	//Get bin centers
        double current_P_bin_center = (P_bin_high + P_bin_low) / 2;  
	double current_theta_bin_center = (theta_bin_high + theta_bin_low) / 2;
	double current_phi_bin_center = (phi_bin_high + phi_bin_low) / 2;   

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
	event_count = N_entries;
	C_tree->Fill();
	
	//Fill TH3's by setting the bin content of each bin to be equal to the covariance value
	C_P_hist->Fill(current_P_bin_center, current_theta_bin_center, current_phi_bin_center, variance_P);
	C_theta_hist->Fill(current_P_bin_center, current_theta_bin_center, current_phi_bin_center, variance_theta);
	C_phi_hist->Fill(current_P_bin_center, current_theta_bin_center, current_phi_bin_center, variance_phi);
	C_P_phi_hist->Fill(current_P_bin_center, current_theta_bin_center, current_phi_bin_center, covariance_P_phi);
	C_P_theta_hist->Fill(current_P_bin_center, current_theta_bin_center, current_phi_bin_center, covariance_P_theta);
	C_theta_phi_hist->Fill(current_P_bin_center, current_theta_bin_center, current_phi_bin_center, covariance_theta_phi);
	entries_hist->Fill(current_P_bin_center, current_theta_bin_center, current_phi_bin_center, N_entries);
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      }
    }
  }  
  // Save the TH3's to the ROOT file
  C_P_hist->Write();
  C_theta_hist->Write();
  C_phi_hist->Write();
  C_P_phi_hist->Write();
  C_P_theta_hist->Write();
  C_theta_phi_hist->Write();
  entries_hist->Write();

  // Save the tree to the ROOT file
  C_tree->Write();

  // Close the ROOT file
  C_file->Close();
  return 0;
}

void read_Hipo(char inputFile[256], int part_charge, float part_mass, std::vector<float>* mc_p_vec, std::vector<float>* rec_p_vec, std::vector<float>* mc_theta_deg_vec, std::vector<float>* rec_theta_deg_vec, std::vector<float>* mc_phi_deg_vec, std::vector<float>* rec_phi_deg_vec, std::vector<float>* p_diff_vec, std::vector<float>* theta_diff_rad_vec, std::vector<float>* phi_diff_rad_vec, std::vector<float> *chi2_vec, std::vector<int> *ndf_vec, int& event_count)
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
  hipo::bank RECCAL(factory.getSchema("REC::Calorimeter"));
  hipo::bank RECTRACK(factory.getSchema("REC::Track"));
      
  int event_num;
  int good_events = 0;
  int total_events = 0;
  //while (nevents < max_events)
  while(reader.next()==true)
  //while(reader.next()==true && (total_events < 1000))
  {
    total_events++;
    reader.read(event);
    event.getStructure(RECPART);
    event.getStructure(MCPART);
    event.getStructure(RECCAL);
    event.getStructure(RECTRACK);
    
    //Get reconstructed particle info
    std::vector<int> rec_pid;
    std::vector<float> rec_px;
    std::vector<float> rec_py;
    std::vector<float> rec_pz;
    std::vector<int> rec_status;
    std::vector<int> rec_sector;
    std::vector<float> chi2;
    std::vector<int> ndf;
    read_Part_Bank(RECPART, RECCAL, RECTRACK, &rec_pid, &rec_px, &rec_py, &rec_pz, &rec_sector, &chi2, &ndf, event_count);

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
	
	//Theta and Phi
	float rec_theta_rad = rec.Theta();
	float rec_theta_deg = rec_theta_rad*180/TMath::Pi();
	float rec_phi_rad = rec.Phi();
	float rec_phi_deg = rec_phi_rad*180/TMath::Pi();
	/*
	//float abs_dist_deg;
	//float rec_phi_modified_deg = rec_phi_deg;
	//float abs_dist_rad;
        //float rec_phi_modified_rad = rec_phi_rad;
	//If phi>180, it will start from -180.To account for this, any phi>180 is modified to continue counting 
	//up from 180. So if the reconstructed phi is -179, the following lines will convert is to phi=181.
	if (rec_phi_deg < -TMath::Pi()/2) {
	  abs_dist_deg = TMath::Abs(180 - TMath::Abs(rec_phi_deg));
	  rec_phi_modified_deg = 180 + abs_dist_deg;
	  abs_dist_rad = TMath::Abs(TMath::Pi() - TMath::Abs(rec_phi_rad));
          rec_phi_modified_rad = TMath::Pi() + abs_dist_rad;
	}
	*/
	float mc_theta_rad = mc.Theta();
	float mc_theta_deg = mc_theta_rad*180/TMath::Pi();
	float mc_phi_rad = mc.Phi();
	float mc_phi_deg = mc_phi_rad*180/TMath::Pi();
	float theta_diff_rad = rec_theta_rad - mc_theta_rad;
	float theta_diff_deg = rec_theta_deg - mc_theta_deg;
	float phi_diff_rad = rec_phi_rad - mc_phi_rad;
	float phi_diff_deg = rec_phi_deg - mc_phi_deg;

	//if ((rec_theta_deg >= 18) && (rec_theta_deg <= 22) && (rec_phi_deg >= 58) && (rec_phi_deg <= 62)) {
	//if ((rec_p >= 6) && (rec_p <= 7) && (rec_phi_deg >= 58) && (rec_phi_deg <= 62)) {
	//if ((rec_p >= 6) && (rec_p <= 7) && (rec_theta_deg >= 18) && (rec_theta_deg <= 22)) {
	  rec_p_vec->push_back(rec_p);
	  mc_p_vec->push_back(mc_p);
	  p_diff_vec->push_back(p_diff);
	  chi2_vec->push_back(chi2[0]);
	  ndf_vec->push_back(ndf[0]);
	  theta_diff_rad_vec->push_back(theta_diff_rad);
	  phi_diff_rad_vec->push_back(phi_diff_rad);
	  rec_theta_deg_vec->push_back(rec_theta_deg);
	  mc_theta_deg_vec->push_back(mc_theta_deg);
	  rec_phi_deg_vec->push_back(rec_phi_deg);
	  mc_phi_deg_vec->push_back(mc_phi_deg);
	  //}
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
void read_Part_Bank(hipo::bank PartBank, hipo::bank CalBank, hipo::bank RecTrack, std::vector<int>* pid, std::vector<float>* px, std::vector<float>* py, std::vector<float>* pz, std::vector<int> *rec_sector, std::vector<float> *chi2, std::vector<int> *ndf, int &event_count){

  event_count++; 
  int nrows = PartBank.getRows();
  int calBank_rows = CalBank.getRows();
  int rectrack_rows = RecTrack.getRows();
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

  int count = 0;
  for (int j = 0; j < nrows; j++){
    float current_chi2;
    int current_ndf;
    float current_chi2_ndf;
    int current_sec;
    int   current_pid = PartBank.getInt("pid",j);
    float  current_px = PartBank.getFloat("px",j);
    float  current_py = PartBank.getFloat("py",j);
    float  current_pz = PartBank.getFloat("pz",j);
    float  current_chi2pid = PartBank.getFloat("chi2pid",j);
    int current_status = PartBank.getInt("status",j);
    int charge = PartBank.getInt("charge",j);
    for (int k = 0; k < rectrack_rows; k++) { 
      int current_pindex = RecTrack.getInt("pindex",k);
      if (current_pindex == j) { 
	current_chi2 = RecTrack.getFloat("chi2",k);
	current_ndf = RecTrack.getInt("NDF",k);
	current_chi2_ndf = current_chi2/current_ndf;
	current_sec = RecTrack.getInt("sector",k);
      }
    }
    //Rather than requiring the correct pid, just require the correct charge
    if (charge == part_charge) {
      px_vec.push_back(current_px);           
      py_vec.push_back(current_py);    
      pz_vec.push_back(current_pz);  
      pid_vec.push_back(current_pid);
      chi2_vec.push_back(current_chi2);
      ndf_vec.push_back(current_ndf);
      status_vec.push_back(current_status);
      sector_vec.push_back(current_sec);
    }
  }

  //Only save the particle info if there is just one particle of the correct charge per event
  if (pid_vec.size() == 1 && (TMath::Abs(status_vec[0]) >= 2000 && TMath::Abs(status_vec[0]) < 3000) && (sector_vec[0] == 2)) {
    px->push_back(px_vec[0]);      
    py->push_back(py_vec[0]);   
    pz->push_back(pz_vec[0]);   
    pid->push_back(pid_vec[0]);
    chi2->push_back(chi2_vec[0]);
    ndf->push_back(ndf_vec[0]);
    }
}

int main() {
    return covMatrix_extraction();
}

