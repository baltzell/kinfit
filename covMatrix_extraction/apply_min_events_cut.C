#include <TFile.h>
#include <TTree.h>
#include <TH3D.h>
#include <iostream>


void apply_min_events_cut() {
  //TFile *file = TFile::Open("/work/clas12/reedtg/clas12_kinematic_fitter/dev_trevor/pim_covariances_7-31-23.root");
  TFile *file = TFile::Open("/work/clas12/reedtg/clas12_kinematic_fitter/dev_trevor/pip_covariances_8-1-23.root");
  TTree *tree = (TTree*)file->Get("covarianceTree");
  TH3D *C_P_hist = (TH3D*)file->Get("C_P_hist");
  
  //Read the binning and covariance info stored in the ROOT file. The binning below assumes P is x-axis, theta is y-axis, and phi is z-axis.
  Int_t P_Bins = C_P_hist->GetNbinsX();
  Double_t P_Min = C_P_hist->GetXaxis()->GetBinLowEdge(1);
  Double_t P_Max = C_P_hist->GetXaxis()->GetBinUpEdge(P_Bins);
  Int_t theta_Bins = C_P_hist->GetNbinsY();
  Double_t theta_Min = C_P_hist->GetYaxis()->GetBinLowEdge(1);
  Double_t theta_Max = C_P_hist->GetYaxis()->GetBinUpEdge(theta_Bins);
  Int_t phi_Bins = C_P_hist->GetNbinsZ();
  Double_t phi_Min = C_P_hist->GetZaxis()->GetBinLowEdge(1);
  Double_t phi_Max = C_P_hist->GetZaxis()->GetBinUpEdge(phi_Bins);
  
  //Create new TH3's for saving the covariance info with the minimum-events cut applied 
  TH3D *C_P_cut = new TH3D("C_P_cut", "C_P Histogram (cut)", P_Bins, P_Min, P_Max, theta_Bins, theta_Min, theta_Max, phi_Bins, phi_Min, phi_Max);
  TH3D *C_theta_cut = new TH3D("C_theta_cut", "C_theta Histogram (cut)", P_Bins, P_Min, P_Max, theta_Bins, theta_Min, theta_Max, phi_Bins, phi_Min, phi_Max);
  TH3D *C_phi_cut = new TH3D("C_phi_cut", "C_phi Histogram (cut)", P_Bins, P_Min, P_Max, theta_Bins, theta_Min, theta_Max, phi_Bins, phi_Min, phi_Max);
  TH3D *C_P_theta_cut = new TH3D("C_P_theta_cut", "C_P_theta Histogram (cut)", P_Bins, P_Min, P_Max, theta_Bins, theta_Min, theta_Max, phi_Bins, phi_Min, phi_Max);
  TH3D *C_P_phi_cut = new TH3D("C_P_phi_cut", "C_P_phi Histogram (cut)", P_Bins, P_Min, P_Max, theta_Bins, theta_Min, theta_Max, phi_Bins, phi_Min, phi_Max);
  TH3D *C_theta_phi_cut = new TH3D("C_theta_phi_cut", "C_theta_phi Histogram (cut)", P_Bins, P_Min, P_Max, theta_Bins, theta_Min, theta_Max, phi_Bins, phi_Min, phi_Max);
  
  double C_P;
  double C_theta;
  double C_phi;
  double C_P_theta;
  double C_P_phi;
  double C_theta_phi;
  double P_bin_center;
  double theta_bin_center;
  double phi_bin_center;
  Int_t event_count;
  
  // Set the branch addresses                                                                                                                                                                                      
  tree->SetBranchAddress("P", &P_bin_center);
  tree->SetBranchAddress("theta", &theta_bin_center);
  tree->SetBranchAddress("phi", &phi_bin_center);
  tree->SetBranchAddress("C_P", &C_P);
  tree->SetBranchAddress("C_theta", &C_theta);
  tree->SetBranchAddress("C_phi", &C_phi);
  tree->SetBranchAddress("C_P_phi", &C_P_phi);
  tree->SetBranchAddress("C_P_theta", &C_P_theta);
  tree->SetBranchAddress("C_theta_phi", &C_theta_phi);
  tree->SetBranchAddress("event_count", &event_count);
  
  //Loop over TTree and fill the TH3 histogram while applying the minimum-events cut:
  int min_events = 100;
  int bins_with_nan = 0;
  int bins_with_0 = 0;
  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    
    // Apply the minimum cut on the number of events
    if (event_count >= min_events) {

      // Check for unwanted values before filling the histograms
      if (std::isnan(C_P) || std::isnan(C_theta) || std::isnan(C_phi) ||
	  std::isnan(C_P_theta) || std::isnan(C_P_phi) || std::isnan(C_theta_phi)) {
	bins_with_nan++;
	continue;
      }

      if (C_P == 0.0 || C_theta == 0.0 || C_phi == 0.0 ||
	  C_P_theta == 0.0 || C_P_phi == 0.0 || C_theta_phi == 0.0 || 
	  C_P == -0.0 || C_theta == -0.0 || C_phi == -0.0 ||
          C_P_theta == -0.0 || C_P_phi == -0.0 || C_theta_phi == -0.0) {
	bins_with_0++;
	continue;
      }
      // Fill the TH3 histograms
      C_P_cut->Fill(P_bin_center, theta_bin_center, phi_bin_center, C_P);
      C_theta_cut->Fill(P_bin_center, theta_bin_center, phi_bin_center, C_theta);
      C_phi_cut->Fill(P_bin_center, theta_bin_center, phi_bin_center, C_phi);
      C_P_theta_cut->Fill(P_bin_center, theta_bin_center, phi_bin_center, C_P_theta);
      C_P_phi_cut->Fill(P_bin_center, theta_bin_center, phi_bin_center, C_P_phi);
      C_theta_phi_cut->Fill(P_bin_center, theta_bin_center, phi_bin_center, C_theta_phi);
    }
  }
  std::cout << "bins_with_nan = " << bins_with_nan << std::endl;
  std::cout << "bins_with_0 = " << bins_with_0 << std::endl;
  
  //Save the histograms to a ROOT file
  TFile *outputFile = TFile::Open("minEventCut_covariances.root", "RECREATE");
  C_P_cut->Write();
  C_theta_cut->Write();
  C_phi_cut->Write();
  C_P_theta_cut->Write();
  C_P_phi_cut->Write();
  C_theta_phi_cut->Write();
  outputFile->Close();
}

int main() {

  apply_min_events_cut();

  return 0;
}
