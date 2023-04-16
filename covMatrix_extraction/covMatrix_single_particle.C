#include <vector>
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "reader.h"
#include <dirent.h>
#include <sys/types.h>
#include <cstring>
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"
#include "TMath.h"
#include <iostream>
#include "TPaveStats.h"
#include "TLatex.h"
#include "TTree.h"
#include "TCut.h"


void read_Hipo(char in_data[256], int part_charge, float part_mass, std::vector<float>* mc_p_vec, std::vector<float>* rec_p_vec, std::vector<float>* mc_theta_deg_vec, std::vector<float>* rec_theta_deg_vec, std::vector<float>* mc_phi_deg_vec, std::vector<float>* rec_phi_deg_vec, std::vector<float>* p_diff_vec, std::vector<float>* theta_diff_rad_vec, std::vector<float>* phi_diff_rad_vec, std::vector<float>* mc_MM2_vec, std::vector<float>* rec_MM2_vec, std::vector<float>* mc_tot_E_vec, std::vector<float>* rec_tot_E_vec, std::vector<float>* mc_tot_pz_vec, std::vector<float>* rec_tot_pz_vec);
void read_Part_Bank(hipo::bank PartBank, std::vector<int>* pid, std::vector<float>* px, std::vector<float>* py, std::vector<float>* pz);
void read_Part_Bank(hipo::bank PartBank, hipo::bank CalBank, std::vector<int>* pid, std::vector<float>* px, std::vector<float>* py, std::vector<float>* pz, std::vector<int> *rec_status);

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

std::ofstream covMatrix_elements;
std::ofstream shifts;
std::ofstream InfoFile("event_info.txt");

//Initialize the particle type and binning info
const TString part_type = "#pi^{+}";
const float part_mass = 0.13957;
const int part_charge = 1;
const int P_bins = 3;
const int theta_bins = 3;
const int phi_bins = 3;

const std::vector<TString> KINES = {"P", "#theta", "#phi"};
const std::vector<TString> UNITS = {"[GeV]", "[rad]", "[rad]"};
const std::vector<int> kin_plots_low = {0, 0, 25};
const std::vector<int> kin_plots_high = {12, 40, 95};
//const std::vector<double> kin_delta_plots_low = {-0.08, -0.002, -0.01};
//const std::vector<double> kin_delta_plots_high = {0.08, 0.002, 0.01};
const std::vector<double> kin_delta_plots_low = {-0.3, -0.01, -0.02};
const std::vector<double> kin_delta_plots_high = {0.3, 0.01, 0.02};

int covMatrix_extraction()
{
  //---------------Specify input file or directory path---------------//                                                                                                                                  
  //char in_data[256] = "/work/clas12/reedtg/clas12_kinematic_fitter/covMatrix/outputs/pip_pim_test/pip_pim_binary_field/cooked/";                         
  //char in_data[256] = "/volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/pip-sec2/cooked/out_pip-sec2513.rec.hipo";
  char in_data[256] = "/volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/pip-sec2/cooked/";                                                                                   
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
        read_Hipo(dir_path_char, part_charge, part_mass, &mc_p_vec, &rec_p_vec, &mc_theta_deg_vec, &rec_theta_deg_vec, &mc_phi_deg_vec, &rec_phi_deg_vec, &p_diff_vec, &theta_diff_rad_vec, &phi_diff_rad_vec, &mc_MM2_vec, &rec_MM2_vec, &mc_tot_E_vec, &rec_tot_E_vec, &mc_tot_pz_vec, &rec_tot_pz_vec);
        //Print list of all read hipo files to input_files.txt                                                                                                                                             
        //InFileList << in_data << en->d_name << std::endl;                                                                                                                                                
      }
    }
    closedir(dr); //close all directory                                                                                                                                                                    
  }
  //If in_data is not a directory, then get the single input file    
  else {
    std::cout << "Not a directory. Only input single file: " << in_data << std::endl;
    read_Hipo(in_data, part_charge, part_mass, &mc_p_vec, &rec_p_vec, &mc_theta_deg_vec, &rec_theta_deg_vec, &mc_phi_deg_vec, &rec_phi_deg_vec, &p_diff_vec, &theta_diff_rad_vec, &phi_diff_rad_vec, &mc_MM2_vec, &rec_MM2_vec, &mc_tot_E_vec, &rec_tot_E_vec, &mc_tot_pz_vec, &rec_tot_pz_vec);
  }

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
  for (int pBin = 0; pBin < P_bins; pBin++) {
    double P_bin_low = P_bin_min + pBin*P_bin_size;
    double P_bin_high = P_bin_low + P_bin_size;
    //std::cout << "P_bin_low = " << P_bin_low << ", P_bin_high = " << P_bin_high << std::endl;
    for (int thetaBin = 0; thetaBin < theta_bins; thetaBin++) {
      double theta_bin_low = theta_bin_min + thetaBin*theta_bin_size;
      double theta_bin_high = theta_bin_low + theta_bin_size;
      //std::cout << "theta_bin_low = "<< theta_bin_low <<", theta_bin_high = " << theta_bin_high<< std::endl;
      for (int phiBin = 0; phiBin < phi_bins; phiBin++) {
	double phi_bin_low = phi_bin_min + phiBin*phi_bin_size;
	double phi_bin_high = phi_bin_low + phi_bin_size;
	//std::cout << "phi_bin_low = "<< phi_bin_low <<", phi_bin_high = " << phi_bin_high<< std::endl;

	//Loop through all saved events and get the info for the current bin
	std::vector<float> mc_p_this_bin, mc_theta_deg_this_bin, mc_phi_deg_this_bin;
	std::vector<float> rec_p_this_bin, rec_theta_deg_this_bin, rec_phi_deg_this_bin;
	std::vector<float> p_diff_this_bin, theta_diff_rad_this_bin, phi_diff_rad_this_bin;
	for (int i = 0; i < mc_p_vec.size(); i++) {
	  if (((rec_p_vec[i] >= P_bin_low) && (rec_p_vec[i] < P_bin_high)) && ((rec_theta_deg_vec[i] >= theta_bin_low) && (rec_theta_deg_vec[i] < theta_bin_high)) && ((rec_phi_deg_vec[i] >= phi_bin_low) && (rec_phi_deg_vec[i] < phi_bin_high))) {
	    mc_p_this_bin.push_back(mc_p_vec[i]);
	    mc_theta_deg_this_bin.push_back(mc_theta_deg_vec[i]);
	    mc_phi_deg_this_bin.push_back(mc_phi_deg_vec[i]);
	    rec_p_this_bin.push_back(rec_p_vec[i]);
	    rec_theta_deg_this_bin.push_back(rec_theta_deg_vec[i]);
	    rec_phi_deg_this_bin.push_back(rec_phi_deg_vec[i]);
	    p_diff_this_bin.push_back(p_diff_vec[i]);
	    theta_diff_rad_this_bin.push_back(theta_diff_rad_vec[i]);
	    phi_diff_rad_this_bin.push_back(phi_diff_rad_vec[i]);
	    InfoFile << "rec_p_vec = " << rec_p_vec[i] << std::endl;
	  }
	}

	auto delta_kin_canvas = new TCanvas("delta_kin_canvas", "delta_kin_canvas", 800, 800);
	delta_kin_canvas->Divide(2, 3);
	auto kin_canvas = new TCanvas("kin_canvas", "kin_canvas", 800, 800);
	kin_canvas->Divide(2, 3);
	auto rec_kin_canvas = new TCanvas("rec_kin_canvas", "rec_kin_canvas", 800, 800);
	auto more_kinematics_can = new TCanvas("more_kinematics_can", "more_kinematics_can", 800, 800);
	more_kinematics_can->Divide(2,3);
  
	//~~~~~~~~~~~~~~~~~~Kinematics plots~~~~~~~~~~~~~~~~~~~~~~~//
	std::vector<TH1*> mc_kin_plots;
	std::vector<TH1*> rec_kin_plots;
	std::vector<TH1*> kin_delta_plots;
	std::vector<TH2*> kin_delta_2D_plots;
	TH1* mc_MM2_plot;
	TH1* rec_MM2_plot;
	for (int j = 0; j < KINES.size(); j++) {
	  TH1* mc_kin_temp = new TH1F(Form("mc_%s_plot", KINES[j].Data()), Form("MC: %s;%s %s;Counts", part_type.Data(), KINES[j].Data(), UNITS[j].Data()), 200, kin_plots_low[j], kin_plots_high[j]);
	  mc_kin_temp->GetXaxis()->SetTitleSize(0.07);
	  mc_kin_temp->GetYaxis()->SetTitleSize(0.05);
	  mc_kin_plots.push_back(mc_kin_temp);
	  TH1* rec_kin_temp = new TH1F(Form("rec_%s_plot", KINES[j].Data()), Form("Rec: %s;%s %s;Counts", part_type.Data(), KINES[j].Data(), UNITS[j].Data()), 200, kin_plots_low[j], kin_plots_high[j]);
	  rec_kin_temp->GetXaxis()->SetTitleSize(0.07);
	  rec_kin_temp->GetYaxis()->SetTitleSize(0.05);
	  rec_kin_plots.push_back(rec_kin_temp);
	  TH1* kin_delta_temp = new TH1F(Form("delta_%s_plot", KINES[j].Data()), Form("%s;#Delta %s (Rec - Gen) %s;Counts", part_type.Data(), KINES[j].Data(), UNITS[j].Data()), 200, kin_delta_plots_low[j], kin_delta_plots_high[j]);
	  kin_delta_temp->GetXaxis()->SetTitleSize(0.07);
	  kin_delta_temp->GetYaxis()->SetTitleSize(0.05);
	  kin_delta_plots.push_back(kin_delta_temp);
	  if (j == KINES.size() - 1) {
	    TH2* kin_delta_2D_temp = new TH2F(Form("delta_%s_delta_%s_plot", KINES[0].Data(), KINES[j].Data()), Form("%s;#Delta %s (Rec - Gen) %s;#Delta %s (Rec - Gen) %s", part_type.Data(), KINES[j].Data(), UNITS[j].Data(), KINES[0].Data(), UNITS[0].Data()), 200, kin_delta_plots_low[j], kin_delta_plots_high[j], 200, kin_delta_plots_low[0], kin_delta_plots_high[0]);
	    kin_delta_2D_temp->GetXaxis()->SetTitleSize(0.07);
	    kin_delta_2D_temp->GetYaxis()->SetTitleSize(0.07);
	    kin_delta_2D_plots.push_back(kin_delta_2D_temp);
	  }
	  else {
	    TH2* kin_delta_2D_temp = new TH2F(Form("delta_%s_delta_%s_plot", KINES[j].Data(), KINES[j+1].Data()), Form("%s;#Delta %s (Rec - Gen) %s;#Delta %s (Rec - Gen) %s", part_type.Data(), KINES[j+1].Data(), UNITS[j+1].Data(), KINES[j].Data(), UNITS[j].Data()), 200, kin_delta_plots_low[j+1], kin_delta_plots_high[j+1], 200, kin_delta_plots_low[j], kin_delta_plots_high[j]);
	    kin_delta_2D_temp->GetXaxis()->SetTitleSize(0.07);
	    kin_delta_2D_temp->GetYaxis()->SetTitleSize(0.07);
	    kin_delta_2D_plots.push_back(kin_delta_2D_temp);
	  }
	}

  
	//--------------------------Fill histograms----------------------------//
	int p_diff_length = p_diff_this_bin.size();
	for (int i = 0; i < p_diff_length; i++) {
	  kin_delta_plots[0]->Fill(p_diff_this_bin[i]);
	  kin_delta_plots[1]->Fill(theta_diff_rad_this_bin[i]);
	  kin_delta_plots[2]->Fill(phi_diff_rad_this_bin[i]);
	}
  
	for (int i = 0; i < mc_p_this_bin.size(); i++) {
	  mc_kin_plots[0]->Fill(mc_p_this_bin[i]);
	  mc_kin_plots[1]->Fill(mc_theta_deg_this_bin[i]);
	  mc_kin_plots[2]->Fill(mc_phi_deg_this_bin[i]);
	}

	for (int i = 0; i < rec_p_this_bin.size(); i++) {
	  rec_kin_plots[0]->Fill(rec_p_this_bin[i]);
	  rec_kin_plots[1]->Fill(rec_theta_deg_this_bin[i]);
	  rec_kin_plots[2]->Fill(rec_phi_deg_this_bin[i]);
	  //rec_phi_vs_theta->Fill(rec_theta_deg_vec[i], rec_phi_deg_vec[i]);
	  //mc_MM2_plot->Fill(mc_MM2_vec[i]);
	  //mc_tot_E_plot->Fill(mc_tot_E_vec[i]);
	  //mc_tot_pz_plot->Fill(mc_tot_pz_vec[i]);
	  //rec_MM2_plot->Fill(rec_MM2_vec[i]);
	  //rec_tot_E_plot->Fill(rec_tot_E_vec[i]);
	  //rec_tot_pz_plot->Fill(rec_tot_pz_vec[i]);
	}  
	//-------------------------------------------------------------------//
  

	//===================================================Draw Histograms===========================================//
	//MC and Rec kinmatics
	for (int jkine = 0; jkine < KINES.size(); jkine++) {
	  kin_canvas->cd(2*jkine+1);
	  gPad->SetBottomMargin(0.15);    //So the axis titles fit on canvas            
	  gPad->SetLeftMargin(0.15);
	  mc_kin_plots[jkine]->Draw();
	  kin_canvas->cd(2*jkine+2);
	  gPad->SetBottomMargin(0.15);    //So the axis titles fit on canvas             
	  gPad->SetLeftMargin(0.15);
	  rec_kin_plots[jkine]->Draw();
	}
  

	//Rec phi vs rec theta
	//rec_kin_canvas->cd();
	//gPad->SetBottomMargin(0.15);
	//gPad->SetLeftMargin(0.15);
	//gPad->SetRightMargin(0.15);
	//rec_phi_vs_theta->Draw("Colz");

  
	int binmax;
	double binmax_bincenter;
	float p_gaus_fit_low = -0.05;
	float p_gaus_fit_high = 0.05;
	float p_fit_low = -0.2;
	float p_fit_high = 0.2;
	float theta_gaus_fit_low = -0.001;
	float theta_gaus_fit_high = 0.001;
	float theta_fit_low = -0.003;
	float theta_fit_high = 0.003;
	float phi_gaus_fit_low = -0.008;
	float phi_gaus_fit_high = 0.008;
	float phi_fit_low = -0.015;
	float phi_fit_high = 0.015;
	std::vector<double> gausfit_means;
	std::vector<double> gaus_fits_low = {p_gaus_fit_low, theta_gaus_fit_low, phi_gaus_fit_low};
	std::vector<double> gaus_fits_high = {p_gaus_fit_high, theta_gaus_fit_high, phi_gaus_fit_high};
	std::vector<double> total_fits_low = {p_fit_low, theta_fit_low, phi_fit_low};
	std::vector<double> total_fits_high = {p_fit_high, theta_fit_high, phi_fit_high};
  
	for (int jkine = 0; jkine < KINES.size(); jkine++) {
	  //////////////1D histograms
	  delta_kin_canvas->cd(2*jkine+1);
	  gStyle->SetTitleFontSize(0.07);
	  gPad->SetBottomMargin(0.2);    //So the axis titles fit on canvas   
	  gPad->SetLeftMargin(0.15);
	  kin_delta_plots[jkine]->Draw();
	  TF1 *gausfunc = new TF1("gausfunc", GaussianPeak, gaus_fits_low[jkine], gaus_fits_high[jkine], 3);
	  gausfunc->SetLineColor(kGreen);
	  gausfunc->SetLineStyle(2);
	  binmax = kin_delta_plots[jkine]->GetMaximumBin();
	  binmax_bincenter = kin_delta_plots[jkine]->GetXaxis()->GetBinCenter(binmax);
	  gausfunc->SetParameters(binmax, binmax_bincenter, 0.001);
	  kin_delta_plots[jkine]->Fit("gausfunc", "q", "", gaus_fits_low[jkine], gaus_fits_high[jkine]);
	  gausfunc->Draw("same");
	  float fit_mean = gausfunc->GetParameter(1);
	  gausfit_means.push_back(fit_mean);
	  
	  //Linear Background
	  TF1 *fitFcn = new TF1("fitFcn", fitFunction_p1_bg, total_fits_low[jkine], total_fits_high[jkine], 5);
	  fitFcn->SetParNames("bg p0", "bg p1", "Amp", "Mean", "Sigma");
	  fitFcn->SetParameter(3, gausfunc->GetParameter(1));
	  fitFcn->SetParameter(4, gausfunc->GetParameter(2));
	  InfoFile << "gausfit mean = " << fit_mean << std::endl;
	  
	  //Clear histogram and refill with same data, but shifting the mean to 0                
	  kin_delta_plots[jkine]->Reset();
	  for (int i = 0; i < p_diff_length; i++) {
	    //InfoFile << KINES[jkine].Data() << std::endl;
	    TString P = "P";
	    TString theta = "#theta";
	    TString phi = "#phi";
	    if (KINES[jkine].Data() == P) {
	      kin_delta_plots[jkine]->Fill(p_diff_this_bin[i] - fit_mean);
	    }
	    else if (KINES[jkine].Data() == theta) {
	      kin_delta_plots[jkine]->Fill(theta_diff_rad_this_bin[i] - fit_mean);
	    }
	    else if (KINES[jkine].Data() == phi) {
	      kin_delta_plots[jkine]->Fill(phi_diff_rad_this_bin[i] - fit_mean);
	    }
	  }
      
	  kin_delta_plots[jkine]->Fit("fitFcn", "q", "", total_fits_low[jkine], total_fits_low[jkine]);
	  Double_t fit_params[5];
	  fitFcn->GetParameters(fit_params);
	  TF1 *pBackFcn = new TF1("pBackFcn", p1Background, total_fits_low[jkine], total_fits_high[jkine], 2);
	  pBackFcn->SetParameters(fit_params);
	  pBackFcn->SetLineStyle(2);
	  pBackFcn->SetLineColor(4);
	  pBackFcn->SetLineWidth(1);
	  pBackFcn->Draw("same");
	  gStyle->SetOptFit(0111);
	  gStyle->SetOptStat("e");
	  
	  ///////////////2D Histograms
	  delta_kin_canvas->cd(2*jkine+2);
	  gStyle->SetTitleFontSize(0.07);
	  gPad->SetBottomMargin(0.2);  //So the axis titles fit on canvas      
	  gPad->SetLeftMargin(0.15);
	  kin_delta_2D_plots[jkine]->Draw("Colz");
	}
  
	/*
	//More Kinematics
	more_kinematics_can->cd(1);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.15);
	mc_MM2_plot->Draw();
	more_kinematics_can->cd(2);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.15);
	rec_MM2_plot->Draw();
	*/
	/*
	more_kinematics_can->cd(3);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.15);
	mc_tot_E_plot->Draw();
	more_kinematics_can->cd(4);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.15);
	rec_tot_E_plot->Draw();
	more_kinematics_can->cd(5);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.15);
	mc_tot_pz_plot->Draw();
	more_kinematics_can->cd(6);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.15);
	rec_tot_pz_plot->Draw();
	*/
	//============================================================================================================================//
	
	//Fill 2D histograms of kinematics deltas with shifted data           
	for (int i = 0; i < phi_diff_rad_this_bin.size(); i++) {
	  kin_delta_2D_plots[2]->Fill(phi_diff_rad_this_bin[i] - gausfit_means[2], p_diff_this_bin[i] - gausfit_means[0]);
	  kin_delta_2D_plots[0]->Fill(theta_diff_rad_this_bin[i] - gausfit_means[1], p_diff_this_bin[i] - gausfit_means[0]);
	  kin_delta_2D_plots[1]->Fill(phi_diff_rad_this_bin[i] - gausfit_means[2], theta_diff_rad_this_bin[i] - gausfit_means[1]);
	}
	//std::cout << "P bin = " << pBin << std::endl;
	delta_kin_canvas->SaveAs(Form("plots/kinematics_P%d_theta%d_phi%d.pdf(", pBin, thetaBin, phiBin));
	kin_canvas->SaveAs(Form("plots/kinematics_P%d_theta%d_phi%d.pdf)", pBin, thetaBin, phiBin));
	//rec_kin_canvas->SaveAs(Form("kinematics.pdf)"));
	//more_kinematics_can->SaveAs(Form("kinematics.pdf)"));

	//Get variances and covariances of the three kinematic variables: P, theta, and phi
	Double_t p_var = kin_delta_2D_plots[2]->GetCovariance(2, 2);
	Double_t phi_var = kin_delta_2D_plots[2]->GetCovariance(1, 1);
	Double_t theta_var = kin_delta_2D_plots[1]->GetCovariance(2, 2);
	Double_t p_phi_cov = kin_delta_2D_plots[2]->GetCovariance(1, 2);
	Double_t p_theta_cov = kin_delta_2D_plots[0]->GetCovariance(1, 2);
	Double_t theta_phi_cov = kin_delta_2D_plots[1]->GetCovariance(1, 2);

	//~~~~~~~~~~~~~~~~~~~Save covarince matrix elements and systematic shifts~~~~~~~~~~~~~~~~~~~~~//
	//Overwrites pre-existing files
	if ((pBin == 0) && (thetaBin == 0) && (phiBin == 0)) {
	  covMatrix_elements.open("cov_matrix_txt_files/matrix_elements_pip_sec2.txt");      
	  shifts.open("systematic_shifts_txt_files/systematic_shifts_pip_sec2.txt");
	}
	else {
	  covMatrix_elements.open("cov_matrix_txt_files/matrix_elements_pip_sec2.txt", std::ios::app);
	  shifts.open("systematic_shifts_txt_files/systematic_shifts_pip_sec2.txt", std::ios::app);
	}
	//Add covariance matrix elements to an output text file
	//First row is for pi+, second row is for pi-
	covMatrix_elements << p_var << " " << theta_var << " " << phi_var << " " << p_phi_cov << " " << p_theta_cov << " " << theta_phi_cov << std::endl;

	//Add the systematic shifts in p, theta, and phi to an output text file
	shifts << gausfit_means[0] << " " << gausfit_means[1] << " " << gausfit_means[2] << std::endl;
  
	covMatrix_elements.close();
	shifts.close();
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      }
    }
  }
  return 0;
}

void read_Hipo(char inputFile[256], int part_charge, float part_mass, std::vector<float>* mc_p_vec, std::vector<float>* rec_p_vec, std::vector<float>* mc_theta_deg_vec, std::vector<float>* rec_theta_deg_vec, std::vector<float>* mc_phi_deg_vec, std::vector<float>* rec_phi_deg_vec, std::vector<float>* p_diff_vec, std::vector<float>* theta_diff_rad_vec, std::vector<float>* phi_diff_rad_vec, std::vector<float>* mc_MM2_vec, std::vector<float>* rec_MM2_vec, std::vector<float>* mc_tot_E_vec, std::vector<float>* rec_tot_E_vec, std::vector<float>* mc_tot_pz_vec, std::vector<float>* rec_tot_pz_vec)
{
  TLorentzVector target(0.0, 0.0, 0.0, 0.0);
  TLorentzVector beam(0.0, 0.0, 10.45, 11.004);
  TLorentzVector W = beam + target;

  //std::cout << "Opening input file " << std::endl;
  hipo::reader  reader;
  reader.open(inputFile);
  hipo::dictionary  factory;
  reader.readDictionary(factory);
  hipo::event      event;
  hipo::bank RECPART(factory.getSchema("REC::Particle"));
  hipo::bank MCPART(factory.getSchema("MC::Particle"));
  hipo::bank RECCAL(factory.getSchema("REC::Calorimeter"));
      
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
    
    //Get reconstructed particle info
    std::vector<int> rec_pid;
    std::vector<float> rec_px;
    std::vector<float> rec_py;
    std::vector<float> rec_pz;
    std::vector<int> rec_status;
    std::vector<int> rec_sector;
    read_Part_Bank(RECPART, RECCAL, &rec_pid, &rec_px, &rec_py, &rec_pz, &rec_sector);

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
	//InfoFile << "rec_p = " << rec_p << ", mc_p = " << mc_p << std::endl;
	float p_diff = rec_p - mc_p;
	float mc_MM2 = (W - mc).M2();
	//float mc_tot_E = mc_pip_E + mc_pim_E;
	//float mc_tot_pz = mc_pz[0] + mc_pz[1];
	float rec_MM2 = (W - rec).M2();
	//float rec_tot_E = rec_pip_E + rec_pim_E;
        //float rec_tot_pz = rec_pz[0] + rec_pz[1];
	mc_MM2_vec->push_back(mc_MM2);
	//mc_tot_E_vec->push_back(mc_tot_E);
	//mc_tot_pz_vec->push_back(mc_tot_pz);
	rec_MM2_vec->push_back(rec_MM2);
        //rec_tot_E_vec->push_back(rec_tot_E);
        //rec_tot_pz_vec->push_back(rec_tot_pz);
	rec_p_vec->push_back(rec_p);
	mc_p_vec->push_back(mc_p);
	p_diff_vec->push_back(p_diff);
	
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
	theta_diff_rad_vec->push_back(theta_diff_rad);
	phi_diff_rad_vec->push_back(phi_diff_rad);
	rec_theta_deg_vec->push_back(rec_theta_deg);
	mc_theta_deg_vec->push_back(mc_theta_deg);
	rec_phi_deg_vec->push_back(rec_phi_deg);
	mc_phi_deg_vec->push_back(mc_phi_deg);
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
void read_Part_Bank(hipo::bank PartBank, hipo::bank CalBank, std::vector<int>* pid, std::vector<float>* px, std::vector<float>* py, std::vector<float>* pz, std::vector<int> *rec_sector){
  
  int nrows = PartBank.getRows();
  int calBank_rows = CalBank.getRows();

  //Temporary vectors
  std::vector<float> px_vec;
  std::vector<float> py_vec;
  std::vector<float> pz_vec;
  std::vector<int> pid_vec;
  std::vector<float> sector_vec;

  int count = 0;
  for (int j = 0; j < nrows; j++){
    int   current_pid = PartBank.getInt("pid",j);
    float  current_px = PartBank.getFloat("px",j);
    float  current_py = PartBank.getFloat("py",j);
    float  current_pz = PartBank.getFloat("pz",j);
    float  current_chi2 = PartBank.getFloat("chi2pid",j);
    int current_status = PartBank.getInt("status",j);
    int charge = PartBank.getInt("charge",j);

    //Rather than requiring the correct pid, just require the correct charge
    if (charge == part_charge) {
      px_vec.push_back(current_px);           
      py_vec.push_back(current_py);    
      pz_vec.push_back(current_pz);  
      pid_vec.push_back(current_pid);
    }
  }

  //Only save the particle info if there is just one particle of the correct charge per event
  if (pid_vec.size() == 1) {
    px->push_back(px_vec[0]);      
    py->push_back(py_vec[0]);   
    pz->push_back(pz_vec[0]);   
    pid->push_back(pid_vec[0]);
  }
}

int main() {
    return covMatrix_extraction();
}

