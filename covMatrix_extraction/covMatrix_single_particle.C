#include <vector>
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "reader.h"
#include <dirent.h>
#include <sys/types.h>
#include <cstring>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"
#include "TMath.h"
#include <iostream>
#include "TPaveStats.h"
#include "TLatex.h"
#include "TTree.h"
#include "TCut.h"


void read_Hipo(char in_data[256], int part_charge, float part_mass, std::vector<float>* mc_p_vec, std::vector<float>* rec_p_vec, std::vector<float>* mc_theta_deg_vec, std::vector<float>* rec_theta_deg_vec, std::vector<float>* mc_phi_deg_vec, std::vector<float>* rec_phi_deg_vec, std::vector<float>* p_diff_vec, std::vector<float>* theta_diff_rad_vec, std::vector<float>* phi_diff_rad_vec, std::vector<float>* mc_MM2_vec, std::vector<float>* rec_MM2_vec, std::vector<float>* mc_tot_E_vec, std::vector<float>* rec_tot_E_vec, std::vector<float>* mc_tot_pz_vec, std::vector<float>* rec_tot_pz_vec, std::vector<float>* chi2_vec, std::vector<int> *ndf_vec);
void read_Part_Bank(hipo::bank PartBank, std::vector<int>* pid, std::vector<float>* px, std::vector<float>* py, std::vector<float>* pz);
void read_Part_Bank(hipo::bank PartBank, hipo::bank CalBank, hipo::bank RecTrack, std::vector<int>* pid, std::vector<float>* px, std::vector<float>* py, std::vector<float>* pz, std::vector<int> *rec_status, std::vector<float> *chi2, std::vector<int> *ndf);

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

std::ofstream covMatrix_elements;
std::ofstream shifts;
std::ofstream InfoFile("event_info.txt");

//Initialize the particle type and binning info
const TString part_type = "#pi^{+}";
const float part_mass = 0.13957;
const int part_charge = 1;
const int P_bins = 1;
const int theta_bins = 1;
const int phi_bins = 1;

const std::vector<TString> KINES = {"P", "#theta", "#phi"};
const std::vector<TString> UNITS = {"[GeV]", "[rad]", "[rad]"};
const std::vector<int> kin_plots_low = {0, 0, 25};
const std::vector<int> kin_plots_high = {12, 40, 95};
//const std::vector<double> kin_delta_plots_low = {-0.08, -0.002, -0.01};
//const std::vector<double> kin_delta_plots_high = {0.08, 0.002, 0.01};
const std::vector<double> kin_delta_plots_low = {-0.5, -0.05, -0.05};
const std::vector<double> kin_delta_plots_high = {0.5, 0.05, 0.05};

int covMatrix_extraction()
{
  //---------------Specify input file or directory path---------------//                                                                                                                                  
  //char in_data[256] = "/work/clas12/reedtg/clas12_kinematic_fitter/covMatrix/outputs/pip_pim_test/pip_pim_binary_field/cooked/";                         
  char in_data[256] = "/volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/pip-sec2/cooked/out_pip-sec2513.rec.hipo";
  //char in_data[256] = "/volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/pip-sec2/cooked/";                                                                                   
  //char in_data[256] = "/volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/test_dir/";
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
  std::vector<float> mc_MM2_vec;
  std::vector<float> mc_tot_E_vec;
  std::vector<float> mc_tot_pz_vec;
  std::vector<float> rec_MM2_vec;
  std::vector<float> rec_tot_E_vec;
  std::vector<float> rec_tot_pz_vec;
  std::vector<float> chi2_vec;
  std::vector<int> ndf_vec;

  DIR *dr;
  struct dirent *en;
  dr = opendir(in_data);                                                                                                                                                
  //If in_data is a directory, loop through all files within it  
  int in_file_count = 0;
  if (dr) {
    while ((en = readdir(dr)) != NULL) {
      if ((strcmp(en->d_name, ".") != 0) && (strcmp(en->d_name, "..") != 0)) {
        in_file_count++;
	std::string dir_path_str(in_data);
        dir_path_str.append(en->d_name);
        //Convert file path string back to char to be used by read_Hipo    
	int string_l = dir_path_str.length();
        char dir_path_char[256];
        strcpy(dir_path_char, dir_path_str.c_str());
	std::cout << "input file = " << dir_path_char << std::endl;
        read_Hipo(dir_path_char, part_charge, part_mass, &mc_p_vec, &rec_p_vec, &mc_theta_deg_vec, &rec_theta_deg_vec, &mc_phi_deg_vec, &rec_phi_deg_vec, &p_diff_vec, &theta_diff_rad_vec, &phi_diff_rad_vec, &mc_MM2_vec, &rec_MM2_vec, &mc_tot_E_vec, &rec_tot_E_vec, &mc_tot_pz_vec, &rec_tot_pz_vec, &chi2_vec, &ndf_vec);
        //Print list of all read hipo files to input_files.txt   
        //InFileList << in_data << en->d_name << std::endl;
      }
    }
    closedir(dr); //close all directory 
  }
  //If in_data is not a directory, then get the single input file    
  else {
    std::cout << "Not a directory. Only input single file: " << in_data << std::endl;
    read_Hipo(in_data, part_charge, part_mass, &mc_p_vec, &rec_p_vec, &mc_theta_deg_vec, &rec_theta_deg_vec, &mc_phi_deg_vec, &rec_phi_deg_vec, &p_diff_vec, &theta_diff_rad_vec, &phi_diff_rad_vec, &mc_MM2_vec, &rec_MM2_vec, &mc_tot_E_vec, &rec_tot_E_vec, &mc_tot_pz_vec, &rec_tot_pz_vec, &chi2_vec, &ndf_vec);
  }



  //---------------------------------------Plots Over All Bins------------------------------------------//
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
    TH2* chi2_ndf_vs_rec_kin_temp = new TH2F(Form("chi2_ndf_%s_plot", KINES[j].Data()), Form("%s;Rec %s %s;#Chi^{2}/NDF",  part_type.Data(), KINES[j].Data(), UNITS[j].Data()), 200, kin_plots_low[j], kin_plots_high[j], 200, 0, 10);
    chi2_ndf_vs_rec_kin_temp->GetXaxis()->SetTitleSize(0.07);
    chi2_ndf_vs_rec_kin_temp->GetYaxis()->SetTitleSize(0.07);
    chi2_ndf_vs_rec_kin_plots.push_back(chi2_ndf_vs_rec_kin_temp);

    TH2* chi2_ndf_vs_delta_kin_temp = new TH2F(Form("chi2_ndf_delta_%s_plot", KINES[j].Data()), Form("%s;#Delta %s (Rec - Gen) %s;#Chi^{2}/NDF", part_type.Data(), KINES[j].Data(), UNITS[j].Data()), 200, kin_delta_plots_low[j], kin_delta_plots_high[j], 200, 0, 10);
    chi2_ndf_vs_delta_kin_temp->GetXaxis()->SetTitleSize(0.07);
    chi2_ndf_vs_delta_kin_temp->GetYaxis()->SetTitleSize(0.07);
    chi2_ndf_vs_delta_kin_plots.push_back(chi2_ndf_vs_delta_kin_temp);

    TH2* chi2_vs_rec_kin_temp = new TH2F(Form("chi2_%s_plot", KINES[j].Data()), Form("%s;Rec %s %s;#Chi^{2}",  part_type.Data(), KINES[j].Data(), UNITS[j].Data()), 200, kin_plots_low[j], kin_plots_high[j], 200, 0, 2000);
    chi2_vs_rec_kin_temp->GetXaxis()->SetTitleSize(0.07);
    chi2_vs_rec_kin_temp->GetYaxis()->SetTitleSize(0.07);
    chi2_vs_rec_kin_plots.push_back(chi2_vs_rec_kin_temp);

    TH2* chi2_vs_delta_kin_temp = new TH2F(Form("chi2_delta_%s_plot", KINES[j].Data()), Form("%s;#Delta %s (Rec - Gen) %s;#Chi^{2}", part_type.Data(), KINES[j].Data(), UNITS[j].Data()), 200, kin_delta_plots_low[j], kin_delta_plots_high[j], 200, 0, 2000);
    chi2_vs_delta_kin_temp->GetXaxis()->SetTitleSize(0.07);
    chi2_vs_delta_kin_temp->GetYaxis()->SetTitleSize(0.07);
    chi2_vs_delta_kin_plots.push_back(chi2_vs_delta_kin_temp);

    TH2* ndf_vs_rec_kin_temp = new TH2F(Form("ndf_%s_plot", KINES[j].Data()), Form("%s;Rec %s %s;NDF",  part_type.Data(), KINES[j].Data(), UNITS[j].Data()), 200, kin_plots_low[j], kin_plots_high[j], 50, 0, 50);
    ndf_vs_rec_kin_temp->GetXaxis()->SetTitleSize(0.07);
    ndf_vs_rec_kin_temp->GetYaxis()->SetTitleSize(0.07);
    ndf_vs_rec_kin_plots.push_back(ndf_vs_rec_kin_temp);

    TH2* ndf_vs_delta_kin_temp = new TH2F(Form("ndf_delta_%s_plot", KINES[j].Data()), Form("%s;#Delta %s (Rec - Gen) %s;NDF", part_type.Data(), KINES[j].Data(), UNITS[j].Data()), 200, kin_delta_plots_low[j], kin_delta_plots_high[j], 50, 0, 50);
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
  chi2_ndf_can->SaveAs("plots/Chi2_ndf_vs_kins.pdf");
  chi2_can->SaveAs("plots/Chi2_vs_kins.pdf");
  ndf_can->SaveAs("plots/NDF_vs_kins.pdf");
  //------------------------------------------------------------------------------------------------------//



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

  //3-dimensional vectors of canvases. From innermost outward: phi bin, theta bin, and then P bin
  std::vector<std::vector<std::vector<TCanvas*>>> kin_can_all_bins(P_bins, std::vector<std::vector<TCanvas*>> (theta_bins, std::vector<TCanvas*> (phi_bins)));
  std::vector<std::vector<std::vector<TCanvas*>>> delta_kin_can_all_bins(P_bins, std::vector<std::vector<TCanvas*>> (theta_bins, std::vector<TCanvas*> (phi_bins)));

  //4-dimensional vectors of histograms. From innermost outward: kinematic variable, phi bin, theta bin, and then P bin 
  std::vector<std::vector<std::vector<std::vector<TH1*>>>> rec_kin_plots_all_bins(P_bins, std::vector<std::vector<std::vector<TH1*>>> (theta_bins, std::vector<std::vector<TH1*>> (phi_bins, std::vector<TH1*> (KINES.size()))));
  std::vector<std::vector<std::vector<std::vector<TH1*>>>> mc_kin_plots_all_bins(P_bins, std::vector<std::vector<std::vector<TH1*>>> (theta_bins, std::vector<std::vector<TH1*>> (phi_bins, std::vector<TH1*> (KINES.size()))));
  std::vector<std::vector<std::vector<std::vector<TH1*>>>> delta_kin_plots_all_bins(P_bins, std::vector<std::vector<std::vector<TH1*>>> (theta_bins, std::vector<std::vector<TH1*>> (phi_bins, std::vector<TH1*> (KINES.size()))));
  std::vector<std::vector<std::vector<std::vector<TH2*>>>> delta_kin_2D_plots_all_bins(P_bins, std::vector<std::vector<std::vector<TH2*>>> (theta_bins, std::vector<std::vector<TH2*>> (phi_bins, std::vector<TH2*> (KINES.size()))));

  //3-dimensional vector for saving bin info. From innermost outward: phi bin, theta bin, and then P bin
  std::vector<std::vector<std::vector<double>>> P_bin_centers_all(P_bins, std::vector<std::vector<double>> (theta_bins, std::vector<double> (phi_bins)));
  std::vector<std::vector<std::vector<double>>> theta_bin_centers_all(P_bins, std::vector<std::vector<double>> (theta_bins, std::vector<double> (phi_bins)));
  std::vector<std::vector<std::vector<double>>> phi_bin_centers_all(P_bins, std::vector<std::vector<double>> (theta_bins, std::vector<double> (phi_bins)));

  //std::vector<int> theta_bin_vec;

  //Loop through the bins, create the histos for each one, and add them to the vectors of histos declared above
  int pushback_count = 0;
  for (int pBin = 0; pBin < P_bins; pBin++) {
    double P_bin_low = P_bin_min + pBin*P_bin_size;
    double P_bin_high = P_bin_low + P_bin_size;
    for (int thetaBin = 0; thetaBin < theta_bins; thetaBin++) {
      double theta_bin_low = theta_bin_min + thetaBin*theta_bin_size;
      double theta_bin_high = theta_bin_low + theta_bin_size;
      for (int phiBin = 0; phiBin < phi_bins; phiBin++) {
	double phi_bin_low = phi_bin_min + phiBin*phi_bin_size;
	double phi_bin_high = phi_bin_low + phi_bin_size;

	//Bin Centers
	double P_bin_center = (P_bin_high + P_bin_low) / 2;
	double theta_bin_center = (theta_bin_high + theta_bin_low) / 2;
	double phi_bin_center = (phi_bin_high + phi_bin_low) / 2;
	P_bin_centers_all[pBin][thetaBin].push_back(P_bin_center);
	theta_bin_centers_all[pBin][thetaBin].push_back(theta_bin_center);
	phi_bin_centers_all[pBin][thetaBin].push_back(phi_bin_center);

	//Initialize Canvases and add to the 3-dimensional vector of canvases, rec_kin_can_all_bins
	TCanvas* kin_canvas = new TCanvas(Form("kin_P%i_th%i_ph%i_canvas", pBin, thetaBin, phiBin), Form("kin_P%i_th%i_ph%i_canvas", pBin, thetaBin, phiBin), 800, 800);
	kin_canvas->Divide(2,3);
	kin_can_all_bins[pBin][thetaBin].push_back(kin_canvas);
	TCanvas* delta_kin_canvas = new TCanvas(Form("delta_kin_P%i_th%i_ph%i_canvas", pBin, thetaBin, phiBin), Form("delta_kin_P%i_th%i_ph%i_canvas", pBin, thetaBin, phiBin), 800, 800);
        delta_kin_canvas->Divide(2,3);
        delta_kin_can_all_bins[pBin][thetaBin].push_back(delta_kin_canvas);
      
	//~~~~~~~~~~~~~~~~~~Set Up Plots~~~~~~~~~~~~~~~~~~~~~~~//
	std::vector<TH1*> mc_kin_plots;
	std::vector<TH1*> rec_kin_plots;
	std::vector<TH1*> kin_delta_plots;
	std::vector<TH2*> kin_delta_2D_plots;
	//TH1* mc_MM2_plot;
	//TH1* rec_MM2_plot;
	for (int j = 0; j < KINES.size(); j++) {
	  TH1* mc_kin_temp = new TH1F(Form("mc_%s_P%i_th%i_ph%i_plot", KINES[j].Data(), pBin, thetaBin, phiBin), Form("MC: %s;%s %s;Counts", part_type.Data(), KINES[j].Data(), UNITS[j].Data()), 200, kin_plots_low[j], kin_plots_high[j]);
	  mc_kin_temp->GetXaxis()->SetTitleSize(0.07);
	  mc_kin_temp->GetYaxis()->SetTitleSize(0.05);
	  mc_kin_plots.push_back(mc_kin_temp);
	  TH1* rec_kin_temp = new TH1F(Form("rec_%s_P%i_th%i_ph%i_plot", KINES[j].Data(), pBin, thetaBin, phiBin), Form("Rec: %s;%s %s;Counts", part_type.Data(), KINES[j].Data(), UNITS[j].Data()), 200, kin_plots_low[j], kin_plots_high[j]);
	  rec_kin_temp->GetXaxis()->SetTitleSize(0.07);
	  rec_kin_temp->GetYaxis()->SetTitleSize(0.05);
	  rec_kin_plots.push_back(rec_kin_temp);
	  TH1* kin_delta_temp = new TH1F(Form("delta_%s_P%i_th%i_ph%i_plot", KINES[j].Data(), pBin, thetaBin, phiBin), Form("%s;#Delta %s (Rec - Gen) %s;Counts", part_type.Data(), KINES[j].Data(), UNITS[j].Data()), 500, kin_delta_plots_low[j], kin_delta_plots_high[j]);
	  kin_delta_temp->GetXaxis()->SetTitleSize(0.07);
	  kin_delta_temp->GetYaxis()->SetTitleSize(0.05);
	  kin_delta_plots.push_back(kin_delta_temp);
	  if (j == KINES.size() - 1) {
	    TH2* kin_delta_2D_temp = new TH2F(Form("delta_%s_delta_%s_P%i_th%i_ph%i_plot", KINES[0].Data(), KINES[j].Data(), pBin, thetaBin, phiBin), Form("%s;#Delta %s (Rec - Gen) %s;#Delta %s (Rec - Gen) %s", part_type.Data(), KINES[j].Data(), UNITS[j].Data(), KINES[0].Data(), UNITS[0].Data()), 500, kin_delta_plots_low[j], kin_delta_plots_high[j], 500, kin_delta_plots_low[0], kin_delta_plots_high[0]);
	    kin_delta_2D_temp->GetXaxis()->SetTitleSize(0.07);
	    kin_delta_2D_temp->GetYaxis()->SetTitleSize(0.07);
	    kin_delta_2D_plots.push_back(kin_delta_2D_temp);
	  }
	  else {
	    TH2* kin_delta_2D_temp = new TH2F(Form("delta_%s_delta_%s_P%i_th%i_ph%i_plot", KINES[j].Data(), KINES[j+1].Data(), pBin, thetaBin, phiBin), Form("%s;#Delta %s (Rec - Gen) %s;#Delta %s (Rec - Gen) %s", part_type.Data(), KINES[j+1].Data(), UNITS[j+1].Data(), KINES[j].Data(), UNITS[j].Data()), 500, kin_delta_plots_low[j+1], kin_delta_plots_high[j+1], 500, kin_delta_plots_low[j], kin_delta_plots_high[j]);
	    kin_delta_2D_temp->GetXaxis()->SetTitleSize(0.07);
	    kin_delta_2D_temp->GetYaxis()->SetTitleSize(0.07);
	    kin_delta_2D_plots.push_back(kin_delta_2D_temp);
	  }
	}
	rec_kin_plots_all_bins[pBin][thetaBin].push_back(rec_kin_plots);
	mc_kin_plots_all_bins[pBin][thetaBin].push_back(mc_kin_plots);
	delta_kin_plots_all_bins[pBin][thetaBin].push_back(kin_delta_plots);
	delta_kin_2D_plots_all_bins[pBin][thetaBin].push_back(kin_delta_2D_plots);
	pushback_count++;
      }
    }
  }

  //Set up TH3 to store the binning info
  TH3* P_theta_and_phi = new TH3F("P_theta_and_phi", "P_theta_and_phi", P_bins, P_bin_min, P_bin_max, theta_bins, theta_bin_min, theta_bin_max, phi_bins, phi_bin_min, phi_bin_max);  

  //Loop through events saved from the input hipo file(s)
  for (int i = 0; i < mc_p_vec.size(); i++) {
    Int_t P_bin_num =  P_theta_and_phi->GetXaxis()->FindBin(rec_p_vec[i]);
    Int_t theta_bin_num =  P_theta_and_phi->GetYaxis()->FindBin(rec_theta_deg_vec[i]);
    Int_t phi_bin_num =  P_theta_and_phi->GetZaxis()->FindBin(rec_phi_deg_vec[i]);
    //InfoFile << "P_bin_num = " << P_bin_num << ", theta_bin_num = " << theta_bin_num << ", phi_bin_num = " << phi_bin_num << std::endl;
    //std::cout << "rec_kin_plots_all_bins size = " << rec_kin_plots_all_bins[0][0].size() << std::endl;
    if (((P_bin_num > 0) && (P_bin_num <= P_bins)) && ((theta_bin_num > 0) && (theta_bin_num <= theta_bins)) && ((phi_bin_num > 0) && (phi_bin_num <= phi_bins))) { 
      rec_kin_plots_all_bins[P_bin_num-1][theta_bin_num-1][2*phi_bin_num-1][0]->Fill(rec_p_vec[i]);
      rec_kin_plots_all_bins[P_bin_num-1][theta_bin_num-1][2*phi_bin_num-1][1]->Fill(rec_theta_deg_vec[i]);
      rec_kin_plots_all_bins[P_bin_num-1][theta_bin_num-1][2*phi_bin_num-1][2]->Fill(rec_phi_deg_vec[i]);
      mc_kin_plots_all_bins[P_bin_num-1][theta_bin_num-1][2*phi_bin_num-1][0]->Fill(mc_p_vec[i]);
      mc_kin_plots_all_bins[P_bin_num-1][theta_bin_num-1][2*phi_bin_num-1][1]->Fill(mc_theta_deg_vec[i]);
      mc_kin_plots_all_bins[P_bin_num-1][theta_bin_num-1][2*phi_bin_num-1][2]->Fill(mc_phi_deg_vec[i]);
      delta_kin_plots_all_bins[P_bin_num-1][theta_bin_num-1][2*phi_bin_num-1][0]->Fill(p_diff_vec[i]);
      delta_kin_plots_all_bins[P_bin_num-1][theta_bin_num-1][2*phi_bin_num-1][1]->Fill(theta_diff_rad_vec[i]);
      delta_kin_plots_all_bins[P_bin_num-1][theta_bin_num-1][2*phi_bin_num-1][2]->Fill(phi_diff_rad_vec[i]);
      delta_kin_2D_plots_all_bins[P_bin_num-1][theta_bin_num-1][2*phi_bin_num-1][0]->Fill(theta_diff_rad_vec[i], p_diff_vec[i]);
      delta_kin_2D_plots_all_bins[P_bin_num-1][theta_bin_num-1][2*phi_bin_num-1][1]->Fill(phi_diff_rad_vec[i], theta_diff_rad_vec[i]);
      delta_kin_2D_plots_all_bins[P_bin_num-1][theta_bin_num-1][2*phi_bin_num-1][2]->Fill(phi_diff_rad_vec[i], p_diff_vec[i]);
    }
  }
  
  //std::cout << "kin_can_all_bins[pBin] size = " << kin_can_all_bins[0][0].size() << std::endl;

  //Loop through the bins and draw the filled histos
  for (int pBin = 0; pBin < P_bins; pBin++) {
    for (int thetaBin = 0; thetaBin < theta_bins; thetaBin++) {
      for (int phiBin = 0; phiBin < phi_bins; phiBin++) {
	std::vector<double> delta_kin_plot_low_vec, delta_kin_plot_upp_vec;
	for (int j = 0; j < KINES.size(); j++) {
	  //Draw the histograms
	  kin_can_all_bins[pBin][thetaBin][2*phiBin+1]->cd(2*j+1);
	  gPad->SetBottomMargin(0.15);    //So the axis titles fit on canvas  
	  gPad->SetLeftMargin(0.15);       
	  rec_kin_plots_all_bins[pBin][thetaBin][2*phiBin+1][j]->Draw();
	  kin_can_all_bins[pBin][thetaBin][2*phiBin+1]->cd(2*j+2);
	  gPad->SetBottomMargin(0.15);    //So the axis titles fit on canvas   
	  gPad->SetLeftMargin(0.15);
	  mc_kin_plots_all_bins[pBin][thetaBin][2*phiBin+1][j]->Draw();
	  gPad->SetBottomMargin(0.15);    //So the axis titles fit on canvas        
	  gPad->SetLeftMargin(0.15);

	  //Set the ranges for these plots according to the mean and RMS
	  delta_kin_can_all_bins[pBin][thetaBin][2*phiBin+1]->cd(2*j+1);
	  delta_kin_plots_all_bins[pBin][thetaBin][2*phiBin+1][j]->Draw();
	  gPad->SetBottomMargin(0.15);    //So the axis titles fit on canvas
	  gPad->SetLeftMargin(0.15);
	  Double_t delta_kin_mean = delta_kin_plots_all_bins[pBin][thetaBin][2*phiBin+1][j]->GetMean();
          Double_t delta_kin_rms = delta_kin_plots_all_bins[pBin][thetaBin][2*phiBin+1][j]->GetRMS();
          double delta_kin_plot_low = delta_kin_mean - delta_kin_rms;
          double delta_kin_plot_upp = delta_kin_mean + delta_kin_rms;
          delta_kin_plot_low_vec.push_back(delta_kin_plot_low);
          delta_kin_plot_upp_vec.push_back(delta_kin_plot_upp);
          delta_kin_plots_all_bins[pBin][thetaBin][2*phiBin+1][j]->GetXaxis()->SetRangeUser(delta_kin_plot_low, delta_kin_plot_upp);
	  //std::cout << "Mean from plot = " << delta_kin_mean << ", RMS from plot = " << delta_kin_rms << std::endl;
	  
	  //Check if rebinning is necessary using Sturge's Rule
	  double original_x_range = kin_delta_plots_high[j] - kin_delta_plots_low[j];
	  double reduced_x_range = delta_kin_plot_upp - delta_kin_plot_low;
	  Int_t entries = delta_kin_plots_all_bins[pBin][thetaBin][2*phiBin+1][j]->GetEntries();
	  int N_bins_orig = delta_kin_plots_all_bins[pBin][thetaBin][2*phiBin+1][j]->GetNbinsX();
	  double N_bins_reduced = reduced_x_range / original_x_range * N_bins_orig;
	  int N_bins_reduced_int = int(N_bins_reduced);
	  int N_bins_Sturges_rule = Sturges_rule(entries);
	  int N_bins_ratio = N_bins_reduced_int / N_bins_Sturges_rule;
	  if (N_bins_ratio >= 2) {
	    delta_kin_plots_all_bins[pBin][thetaBin][2*phiBin+1][j]->Rebin(N_bins_ratio);
	    InfoFile << "Reducing binning by factor of " << N_bins_ratio << std::endl;
	  }
	}
	//Loop again through the kinematic variables to apply x and y axis ranges (as determined above) to the 2D histograms
	for (int j = 0; j < KINES.size(); j++) {
	  delta_kin_can_all_bins[pBin][thetaBin][2*phiBin+1]->cd(2*j+2);
          gStyle->SetTitleFontSize(0.07);
          gPad->SetBottomMargin(0.2);  //So the axis titles fit on canvas      
	  gPad->SetLeftMargin(0.17);
          delta_kin_2D_plots_all_bins[pBin][thetaBin][2*phiBin+1][j]->Draw("Colz");
	  //Change 2D histogram axis ranges according to means and RMS from 1D histograms       
	    if (j == KINES.size() - 1) {
	      delta_kin_2D_plots_all_bins[pBin][thetaBin][2*phiBin+1][j]->GetXaxis()->SetRangeUser(delta_kin_plot_low_vec[j], delta_kin_plot_upp_vec[j]);
	      delta_kin_2D_plots_all_bins[pBin][thetaBin][2*phiBin+1][j]->GetYaxis()->SetRangeUser(delta_kin_plot_low_vec[0], delta_kin_plot_upp_vec[0]);
	      //std::cout << "x low = " << delta_kin_plot_low_vec[j] << ", x high = " << delta_kin_plot_upp_vec[j] << std::endl;
	      //std::cout<< "y low = " << delta_kin_plot_low_vec[0] << ", y high = " << delta_kin_plot_upp_vec[0] << std::endl;
	    }
	    else {
	      delta_kin_2D_plots_all_bins[pBin][thetaBin][2*phiBin+1][j]->GetXaxis()->SetRangeUser(delta_kin_plot_low_vec[j+1], delta_kin_plot_upp_vec[j+1]);
	      delta_kin_2D_plots_all_bins[pBin][thetaBin][2*phiBin+1][j]->GetYaxis()->SetRangeUser(delta_kin_plot_low_vec[j], delta_kin_plot_upp_vec[j]);
	    }
	}
	delta_kin_can_all_bins[pBin][thetaBin][2*phiBin+1]->SaveAs(Form("plots/kinematics_P%i_theta%i_phi%i.pdf(", pBin, thetaBin, phiBin));
	kin_can_all_bins[pBin][thetaBin][2*phiBin+1]->SaveAs(Form("plots/kinematics_P%i_theta%i_phi%i.pdf)", pBin, thetaBin, phiBin));

	//Get variances and covariances of the three kinematic variables: P, theta, and phi   
	Double_t p_var = delta_kin_2D_plots_all_bins[pBin][thetaBin][2*phiBin+1][2]->GetCovariance(2, 2);
        Double_t phi_var = delta_kin_2D_plots_all_bins[pBin][thetaBin][2*phiBin+1][2]->GetCovariance(1, 1);
        Double_t theta_var = delta_kin_2D_plots_all_bins[pBin][thetaBin][2*phiBin+1][1]->GetCovariance(2, 2);
        Double_t p_phi_cov = delta_kin_2D_plots_all_bins[pBin][thetaBin][2*phiBin+1][2]->GetCovariance(1, 2);
        Double_t p_theta_cov = delta_kin_2D_plots_all_bins[pBin][thetaBin][2*phiBin+1][0]->GetCovariance(1, 2);
        Double_t theta_phi_cov = delta_kin_2D_plots_all_bins[pBin][thetaBin][2*phiBin+1][1]->GetCovariance(1, 2);

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Save covariance matrix elements~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	//Get bin centers
	double current_P_bin_center = P_bin_centers_all[pBin][thetaBin][2*phiBin+1];
	double current_theta_bin_center = theta_bin_centers_all[pBin][thetaBin][2*phiBin+1];
	double current_phi_bin_center = phi_bin_centers_all[pBin][thetaBin][2*phiBin+1];

	//Overwrites pre-existing files when on the first bin. For all subsequent bins, the file is appended
	if ((pBin == 0) && (thetaBin == 0) && (phiBin == 0)) {
          covMatrix_elements.open("cov_matrix_txt_files/matrix_elements_pip_sec2.txt");
        }
        else {
          covMatrix_elements.open("cov_matrix_txt_files/matrix_elements_pip_sec2.txt", std::ios::app);
        }
	//Add covariance matrix elements to an output text file
	covMatrix_elements << p_var << " " << theta_var << " " << phi_var << " " << p_phi_cov << " " << p_theta_cov << " " << theta_phi_cov << " " << current_P_bin_center << " " << current_theta_bin_center << " " << current_phi_bin_center << std::endl;
	covMatrix_elements.close();
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      }
    }
  }  
  return 0;
}

  void read_Hipo(char inputFile[256], int part_charge, float part_mass, std::vector<float>* mc_p_vec, std::vector<float>* rec_p_vec, std::vector<float>* mc_theta_deg_vec, std::vector<float>* rec_theta_deg_vec, std::vector<float>* mc_phi_deg_vec, std::vector<float>* rec_phi_deg_vec, std::vector<float>* p_diff_vec, std::vector<float>* theta_diff_rad_vec, std::vector<float>* phi_diff_rad_vec, std::vector<float>* mc_MM2_vec, std::vector<float>* rec_MM2_vec, std::vector<float>* mc_tot_E_vec, std::vector<float>* rec_tot_E_vec, std::vector<float>* mc_tot_pz_vec, std::vector<float>* rec_tot_pz_vec, std::vector<float> *chi2_vec, std::vector<int> *ndf_vec)
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
    read_Part_Bank(RECPART, RECCAL, RECTRACK, &rec_pid, &rec_px, &rec_py, &rec_pz, &rec_sector, &chi2, &ndf);

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
	//float mc_tot_E = mc_pip_E + mc_pim_E;
	//float mc_tot_pz = mc_pz[0] + mc_pz[1];
	float rec_MM2 = (W - rec).M2();
	//float rec_tot_E = rec_pip_E + rec_pim_E;
        //float rec_tot_pz = rec_pz[0] + rec_pz[1];
	
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
	  mc_MM2_vec->push_back(mc_MM2);
	  rec_MM2_vec->push_back(rec_MM2);
	  theta_diff_rad_vec->push_back(theta_diff_rad);
	  phi_diff_rad_vec->push_back(phi_diff_rad);
	  rec_theta_deg_vec->push_back(rec_theta_deg);
	  mc_theta_deg_vec->push_back(mc_theta_deg);
	  rec_phi_deg_vec->push_back(rec_phi_deg);
	  mc_phi_deg_vec->push_back(mc_phi_deg);
	  //mc_tot_E_vec->push_back(mc_tot_E);  
	  //mc_tot_pz_vec->push_back(mc_tot_pz);       
	  //rec_tot_E_vec->push_back(rec_tot_E); 
	  //rec_tot_pz_vec->push_back(rec_tot_pz);
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
 void read_Part_Bank(hipo::bank PartBank, hipo::bank CalBank, hipo::bank RecTrack, std::vector<int>* pid, std::vector<float>* px, std::vector<float>* py, std::vector<float>* pz, std::vector<int> *rec_sector, std::vector<float> *chi2, std::vector<int> *ndf){
  
  int nrows = PartBank.getRows();
  int calBank_rows = CalBank.getRows();
  int rectrack_rows = RecTrack.getRows();

  //Temporary vectors
  std::vector<float> px_vec;
  std::vector<float> py_vec;
  std::vector<float> pz_vec;
  std::vector<int> pid_vec;
  std::vector<float> sector_vec;
  std::vector<float> chi2_vec;
  std::vector<int> ndf_vec;
  std::vector<int> status_vec;

  int count = 0;
  for (int j = 0; j < nrows; j++){
    float current_chi2;
    int current_ndf;
    float current_chi2_ndf;
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
    }
  }

  //Only save the particle info if there is just one particle of the correct charge per event
  if (pid_vec.size() == 1 && (TMath::Abs(status_vec[0]) >= 2000 && TMath::Abs(status_vec[0]) < 3000)) { //y && (chi2_ndf_vec[0] > 1 && chi2_ndf_vec[0] < 10)) {
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

