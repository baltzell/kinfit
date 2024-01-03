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
#include <TFile.h>
#include "TTree.h"
#include "TCut.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include <set>

/*
TGraph* tgraph_setup(TGraph* graph, double markersize, Int_t markercolor, Int_t markerstyle) {
  graph->SetMarkerColor(markercolor);
  graph->SetMarkerStyle(markerstyle);
  graph->SetMarkerSize(markersize);

  return graph;
}

TGraphErrors* tgraph_setup(TGraphErrors* graph, double markersize, Int_t markercolor, Int_t markerstyle) {
  graph->SetMarkerColor(markercolor);
  graph->SetMarkerStyle(markerstyle);
  graph->SetMarkerSize(markersize);

  return graph;
}
*/
TGraphErrors* tgraph_setup(TGraphErrors* graph, const char* title, double titlesize, double maintitlesize, double markersize, Int_t markercolor, Int_t markerstyle, Int_t maxdigits, TString kin_axis_title,double bottommargin, double topmargin, double leftmargin) {
  graph->SetTitle(title);
  graph->GetYaxis()->SetTitle(title);
  graph->GetYaxis()->SetTitleSize(titlesize);
  graph->GetYaxis()->SetLabelSize(titlesize);
  graph->SetMarkerColor(markercolor); 
  graph->SetMarkerStyle(markerstyle); 
  graph->SetMarkerSize(markersize); 
  graph->GetYaxis()->SetMaxDigits(maxdigits);
  graph->GetXaxis()->SetTitle(kin_axis_title);
  graph->GetXaxis()->SetTitleSize(titlesize);
  graph->GetXaxis()->SetLabelSize(titlesize);
  gPad->SetBottomMargin(bottommargin);
  gPad->SetTopMargin(topmargin);
  gPad->SetLeftMargin(leftmargin);
  gStyle->SetTitleFontSize(maintitlesize);

  return graph;
}

const std::vector<TString> KINES = {"P", "#theta", "#phi"};
const std::vector<TString> UNITS = {"[GeV]", "[rad]", "[rad]"};
int covMatrix_plots()
{
  //Get covariance matrix values from file
  int N_lines = 0;
  //const char* filename = "/work/clas12/reedtg/clas12_kinematic_fitter/covMatrix/pim_covariances.root";
  //const char* filename = "/work/clas12/reedtg/clas12_kinematic_fitter/dev_trevor/pip_covariances_test.root";
  //const char* filename = "/work/clas12/reedtg/clas12_kinematic_fitter/dev_trevor/cov_root_files/getQuantile/covariances_pip_iterative_outlier_cut_10000_in_files_12-2-23.root";
  const char* filename = "/work/clas12/reedtg/clas12_kinematic_fitter/dev_trevor/cov_root_files/covariances_pip_10000_in_files.root";
  TFile *cov_inFile = TFile::Open(filename); 

  // Get the TTree from the file
  TTree* tree = dynamic_cast<TTree*>(cov_inFile->Get("covarianceTree")); 
  if (!tree) {
    std::cerr << "Error retrieving TTree from file: " << cov_inFile << std::endl;
    cov_inFile->Close();
  }

  //Choose whether or not error of covariance values should be included. (Only set to true if these were calculated along with the covariances)
  bool errors = true;

  double C_P;
  double C_P_err;
  double C_theta;
  double C_theta_err;
  double C_phi;
  double C_phi_err;
  double C_P_theta;
  double C_P_theta_err;
  double C_P_phi;
  double C_P_phi_err;
  double C_theta_phi;
  double C_theta_phi_err;
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
  if (errors) {
    tree->SetBranchAddress("C_P_err", &C_P_err);
    tree->SetBranchAddress("C_theta_err", &C_theta_err);
    tree->SetBranchAddress("C_phi_err", &C_phi_err);
    tree->SetBranchAddress("C_P_phi_err", &C_P_phi_err);
    tree->SetBranchAddress("C_P_theta_err", &C_P_theta_err);
    tree->SetBranchAddress("C_theta_phi_err", &C_theta_phi_err);
  }

  //Vectors to be filled with information from covariance tree
  std::vector<double> P_bin_vec;
  std::vector<double> theta_bin_vec;
  std::vector<double> phi_bin_vec;
  std::vector<double> C_P_vec;
  std::vector<double> C_theta_vec;
  std::vector<double> C_phi_vec;
  std::vector<double> C_P_phi_vec;
  std::vector<double> C_P_theta_vec;
  std::vector<double> C_theta_phi_vec;
  std::vector<Int_t> event_count_vec;
  std::vector<double> C_P_err_vec;
  std::vector<double> C_theta_err_vec;
  std::vector<double> C_phi_err_vec;
  std::vector<double> C_P_phi_err_vec;
  std::vector<double> C_P_theta_err_vec;
  std::vector<double> C_theta_phi_err_vec;
  
  // Get the number of entries in the TTree (this is equivalent to the number of bins)
  Long64_t numEntries = tree->GetEntries();
  // Loop over the entries and read the branch values
  for (Long64_t entry = 0; entry < numEntries; ++entry) {
    tree->GetEntry(entry);
    if (entry < 1000) {
      //std::cout << "events = " << event_count << std::endl;
      //std::cout << "P_bin_center = " << P_bin_center << std::endl;
      continue;
    }
    
    P_bin_vec.push_back(P_bin_center);
    //std::cout << "P_bin_center = " << P_bin_center << std::endl;
    theta_bin_vec.push_back(theta_bin_center);
    //std::cout << "theta_bin_center = " << theta_bin_center << std::endl;
    phi_bin_vec.push_back(phi_bin_center);
    //std::cout << "phi_bin_center = " << phi_bin_center << std::endl;
    C_P_vec.push_back(C_P);
    C_theta_vec.push_back(C_theta);
    C_phi_vec.push_back(C_phi);
    C_P_phi_vec.push_back(C_P_phi);
    C_P_theta_vec.push_back(C_P_theta);
    C_theta_phi_vec.push_back(C_theta_phi);
    event_count_vec.push_back(event_count);
    if (errors) {
      //std::cout << "C_P_err: " << C_P_err << std::endl;
      C_P_err_vec.push_back(C_P_err);
      C_theta_err_vec.push_back(C_theta_err);
      C_phi_err_vec.push_back(C_phi_err);
      C_P_phi_err_vec.push_back(C_P_phi_err);
      C_P_theta_err_vec.push_back(C_P_theta_err);
      C_theta_phi_err_vec.push_back(C_theta_phi_err);
    }
  }

  int N_fixed_variables = 2;
  std::string varying_kin = "p";
  //std::string varying_kin = "theta";
  //std::string varying_kin = "phi";

  double fixed_P_bin;
  double fixed_theta_bin;
  double fixed_phi_bin;

  std::set<double> unique_P_bins(P_bin_vec.begin(), P_bin_vec.end());
  std::set<double> unique_theta_bins(theta_bin_vec.begin(), theta_bin_vec.end());
  std::set<double> unique_phi_bins(phi_bin_vec.begin(), phi_bin_vec.end());
  //std::cout << "Number of P bins = " << unique_P_bins.size() << std::endl;

  //Loop through all the bins. The TTree will be cut and scatterplots will be made for each. 
  std::set<double> first_kin_bins;
  std::set<double> second_kin_bins;
  std::string first_kin;
  std::string second_kin;
  
  if (varying_kin == "p") {
    first_kin_bins = unique_theta_bins;
    first_kin = "theta";
    second_kin_bins = unique_phi_bins;
    second_kin = "phi";
  }
  
  else if (varying_kin == "theta") {
    first_kin_bins = unique_P_bins;
    first_kin = "P";
    second_kin_bins = unique_phi_bins;
    second_kin = "phi";
  }
  else if (varying_kin == "phi") {
    first_kin_bins = unique_P_bins;
    first_kin = "P";
    second_kin_bins = unique_theta_bins;
    second_kin = "theta";
  }
  
  TString pdfFileName = Form("plots/getQuantile/C_vs_%s_plots.pdf", TString(varying_kin.c_str()).Data());
  TCanvas *cov_mat_elements_vs_kin_can = new TCanvas("cov_mat_elements_vs_kin_can", "Covariance Matrix Elements as Function of Kinematics", 800, 1000);
  cov_mat_elements_vs_kin_can->Print(pdfFileName + "[");  // Open the PDF file for writing
  cov_mat_elements_vs_kin_can->cd();
  cov_mat_elements_vs_kin_can->Divide(2,3);

  TTree* cutTree = nullptr;
  int kin1_bin_count = 0;
  int kin2_bin_count = 0;
  //for (int i=0; i<3; i++) {
  for (auto kin1_element : first_kin_bins) {
    kin1_bin_count++;
 
    for (auto kin2_element : second_kin_bins) {
      kin2_bin_count++;
        
      // Create a new TTree with the cut data
      char cutExpression[100];  // Buffer to hold the theta cut expression
      sprintf(cutExpression, "(%s > %f && %s < %f) && (%s > %f && %s < %f)", 
	      first_kin.c_str(), kin1_element-0.01, first_kin.c_str(), kin1_element+0.01, 
	      second_kin.c_str(), kin2_element-0.01, second_kin.c_str(), kin2_element+0.01);

      //std::cout << "Cut expression: " << cutExpression << std::endl;
      // Release memory if cutTree is already allocated
      if (cutTree) {
        delete cutTree;
        cutTree = nullptr;
      }
      //cov_mat_elements_vs_kin_can->Clear();

      cutTree = tree->CopyTree(cutExpression);

      
      //Vectors to be filled with information from the cut tree
      std::vector<double> bin_center_vec;
      std::vector<double> P_bin_cut_vec;
      std::vector<double> theta_bin_cut_vec;
      std::vector<double> phi_bin_cut_vec;
      std::vector<double> C_P_cut_vec;
      std::vector<double> C_theta_cut_vec;
      std::vector<double> C_phi_cut_vec;
      std::vector<double> C_P_phi_cut_vec;
      std::vector<double> C_P_theta_cut_vec;
      std::vector<double> C_theta_phi_cut_vec;
      std::vector<int> event_count_cut_vec;
      std::vector<double> C_P_err_cut_vec;
      std::vector<double> C_theta_err_cut_vec;
      std::vector<double> C_phi_err_cut_vec;
      std::vector<double> C_P_phi_err_cut_vec;
      std::vector<double> C_P_theta_err_cut_vec;
      std::vector<double> C_theta_phi_err_cut_vec;
 
      //Loop through the cut tree to get the values for plotting
      TString kin_axis_title;
      int min_events = 1000;
      int passed_min_events = 0;
      Long64_t numCutEntries = cutTree->GetEntries();
      for (Long64_t entry = 0; entry < numCutEntries; ++entry) {
	cutTree->GetEntry(entry);
	//std::cout << "events = " << event_count << std::endl;
	if (event_count >= min_events) {    
	  passed_min_events++;
	  P_bin_cut_vec.push_back(P_bin_center);
	  theta_bin_cut_vec.push_back(theta_bin_center);
	  phi_bin_cut_vec.push_back(phi_bin_center);
	  C_P_cut_vec.push_back(C_P);
	  C_theta_cut_vec.push_back(C_theta);
	  C_phi_cut_vec.push_back(C_phi);
	  C_P_phi_cut_vec.push_back(C_P_phi);
	  C_P_theta_cut_vec.push_back(C_P_theta);
	  C_theta_phi_cut_vec.push_back(C_theta_phi);
	  event_count_cut_vec.push_back(event_count);
	  if (errors) {
	    C_P_err_cut_vec.push_back(C_P_err);
	    C_theta_err_cut_vec.push_back(C_theta_err);
	    C_phi_err_cut_vec.push_back(C_phi_err);
	    C_P_phi_err_cut_vec.push_back(C_P_phi_err);
	    C_P_theta_err_cut_vec.push_back(C_P_theta_err);
	    C_theta_phi_err_cut_vec.push_back(C_theta_phi_err);
	  }
	  
	  //These values below depend on which kinematics are fixed
	  if (varying_kin == "p") {
	    bin_center_vec.push_back(P_bin_center);
	    kin_axis_title = "P";
	  }
	  else if (varying_kin == "theta") {
	    bin_center_vec.push_back(theta_bin_center);
	    kin_axis_title = "#theta";
	  }
	  else if (varying_kin == "phi") {
	    bin_center_vec.push_back(phi_bin_center);
	    kin_axis_title = "#phi";
	  }
	}
	else { 
	  std::cout << "events = " << event_count << std::endl;
	}
      }

      //Get the TH3's from the ROOT file
      TH3D* C_P_hist = dynamic_cast<TH3D*>(cov_inFile->Get("C_P_hist"));
      TH3D* C_phi_hist = dynamic_cast<TH3D*>(cov_inFile->Get("C_phi_hist"));
      TH3D* C_theta_hist = dynamic_cast<TH3D*>(cov_inFile->Get("C_theta_hist"));
      TH3D* C_P_phi_hist = dynamic_cast<TH3D*>(cov_inFile->Get("C_P_phi_hist"));
      TH3D* C_P_theta_hist = dynamic_cast<TH3D*>(cov_inFile->Get("C_P_theta_hist"));
      TH3D* C_theta_phi_hist = dynamic_cast<TH3D*>(cov_inFile->Get("C_theta_phi_hist"));
      //if (!hist) {
      //  cout << "Error: Failed to retrieve the C_P TH3 histogram: " << C_P_hist << " from the file." << endl;
      //  return;
      //}
      
      Int_t C_P_fixed_kin_1_bin;
      Int_t C_P_fixed_kin_2_bin;
      Int_t C_phi_fixed_kin_1_bin;
      Int_t C_phi_fixed_kin_2_bin;
      Int_t C_theta_fixed_kin_1_bin;
      Int_t C_theta_fixed_kin_2_bin;
      Int_t C_P_phi_fixed_kin_1_bin;
      Int_t C_P_phi_fixed_kin_2_bin;
      Int_t C_P_theta_fixed_kin_1_bin;
      Int_t C_P_theta_fixed_kin_2_bin;
      Int_t C_theta_phi_fixed_kin_1_bin;
      Int_t C_theta_phi_fixed_kin_2_bin;

      C_P_fixed_kin_1_bin = C_P_hist->GetYaxis()->FindBin(kin1_element);
      C_P_fixed_kin_2_bin = C_P_hist->GetZaxis()->FindBin(kin2_element);
      C_phi_fixed_kin_1_bin = C_phi_hist->GetYaxis()->FindBin(kin1_element);
      C_phi_fixed_kin_2_bin = C_phi_hist->GetZaxis()->FindBin(kin2_element);
      C_theta_fixed_kin_1_bin = C_theta_hist->GetYaxis()->FindBin(kin1_element);
      C_theta_fixed_kin_2_bin = C_theta_hist->GetZaxis()->FindBin(kin2_element);
      C_P_phi_fixed_kin_1_bin = C_P_phi_hist->GetYaxis()->FindBin(kin1_element);
      C_P_phi_fixed_kin_2_bin = C_P_phi_hist->GetZaxis()->FindBin(kin2_element);
      C_P_theta_fixed_kin_1_bin = C_P_theta_hist->GetYaxis()->FindBin(kin1_element);
      C_P_theta_fixed_kin_2_bin = C_P_theta_hist->GetZaxis()->FindBin(kin2_element);
      C_theta_phi_fixed_kin_1_bin = C_theta_phi_hist->GetYaxis()->FindBin(kin1_element);
      C_theta_phi_fixed_kin_2_bin = C_theta_phi_hist->GetZaxis()->FindBin(kin2_element);

      cov_mat_elements_vs_kin_can->cd();

      TPaveText *titleText = new TPaveText(0.1, 0.95, 0.9, 0.99, "NDC");
      titleText->SetFillColor(0); // Set background color (0 = white)
      titleText->SetTextFont(42); // Set font
      titleText->SetTextSize(0.04); // Set text size
      titleText->SetBorderSize(0);
      TString title;
      if (varying_kin == "p") {
	title.Form("#theta=%.2f, #phi=%.2f", kin1_element, kin2_element);
	//std::cout << "fixed_theta_bin: " << kin1_element << std::endl;
      }
      else if (varying_kin == "theta") {
	title.Form("P=%.2f, #phi=%.2f", kin1_element, kin2_element);
      }
      else if (varying_kin == "phi") {
	title.Form("P=%.2f, #theta=%.2f", kin1_element, kin2_element);
      }

      titleText->AddText(title);
      titleText->Draw();

      //-------------------------Set up the covariance scatterplots------------------------// 
  
    //If using error bars
    //if (errors) {
      TGraphErrors *C_P_plot = new TGraphErrors(passed_min_events, bin_center_vec.data(), C_P_cut_vec.data(), nullptr, C_P_err_cut_vec.data());
      TGraphErrors *C_theta_plot = new TGraphErrors(passed_min_events, bin_center_vec.data(), C_theta_cut_vec.data(), nullptr, C_theta_err_cut_vec.data());
      TGraphErrors *C_phi_plot = new TGraphErrors(passed_min_events, bin_center_vec.data(), C_phi_cut_vec.data(), nullptr, C_phi_err_cut_vec.data());
      TGraphErrors *C_P_phi_plot = new TGraphErrors(passed_min_events, bin_center_vec.data(), C_P_phi_cut_vec.data(), nullptr, C_P_phi_err_cut_vec.data());
      TGraphErrors *C_P_theta_plot = new TGraphErrors(passed_min_events, bin_center_vec.data(), C_P_theta_cut_vec.data(), nullptr, C_P_theta_err_cut_vec.data());
      TGraphErrors *C_theta_phi_plot = new TGraphErrors(passed_min_events, bin_center_vec.data(), C_theta_phi_cut_vec.data(), nullptr, C_theta_phi_err_cut_vec.data());
      //}
      /*
      //If NOT using error bars
      else {
      TGraph *C_P_plot = new TGraph(passed_min_events, bin_center_vec.data(), C_P_cut_vec.data());
      TGraph *C_theta_plot = new TGraph(passed_min_events, bin_center_vec.data(), C_theta_cut_vec.data());
      TGraph *C_phi_plot = new TGraph(passed_min_events, bin_center_vec.data(), C_phi_cut_vec.data());
      TGraph *C_P_phi_plot = new TGraph(passed_min_events, bin_center_vec.data(), C_P_phi_cut_vec.data());
      TGraph *C_P_theta_plot = new TGraph(passed_min_events, bin_center_vec.data(), C_P_theta_cut_vec.data());
      TGraph *C_theta_phi_plot = new TGraph(passed_min_events, bin_center_vec.data(), C_theta_phi_cut_vec.data());
      }
      */
      //-----------------------------------------------------------------------------------//

      //===================================Variances Plots==========================================//
      //Momentum
      cov_mat_elements_vs_kin_can->cd(1);
      C_P_plot = tgraph_setup(C_P_plot, "C_{P}", 0.07, 0.1, 0.7, kBlue, 8, 2, kin_axis_title, 0.15, 0.3, 0.2);
      C_P_plot->Draw("AP");
      
      //Phi
      cov_mat_elements_vs_kin_can->cd(3);
      C_phi_plot = tgraph_setup(C_phi_plot, "C_{#phi}", 0.07, 0.1, 0.7, kBlue, 8, 2, kin_axis_title, 0.15, 0.15, 0.2);
      C_phi_plot->Draw("AP");
      
      //Theta
      cov_mat_elements_vs_kin_can->cd(5);
      C_theta_plot = tgraph_setup(C_theta_plot, "C_{#theta}", 0.07, 0.1, 0.7, kBlue, 8, 2, kin_axis_title, 0.15, 0.15, 0.2);
      C_theta_plot->Draw("AP");
      //=======================================================================================//
      
      //----------------------------------Covariance Plots-------------------------------------//
      //Momentum and Phi
      cov_mat_elements_vs_kin_can->cd(2);
      C_P_phi_plot = tgraph_setup(C_P_phi_plot, "C_{P#phi}", 0.07, 0.1, 0.7, kBlue, 8, 2, kin_axis_title, 0.15, 0.3, 0.2);
      C_P_phi_plot->Draw("AP");
      
      //Momentum and Theta
      cov_mat_elements_vs_kin_can->cd(4);
      C_P_theta_plot = tgraph_setup(C_P_theta_plot, "C_{P#theta}", 0.07, 0.1, 0.7, kBlue, 8, 2, kin_axis_title, 0.15, 0.15, 0.2);
      C_P_theta_plot->Draw("AP");
      
      //Theta and Phi
      cov_mat_elements_vs_kin_can->cd(6);
      C_theta_phi_plot = tgraph_setup(C_theta_phi_plot, "C_{#theta#phi}", 0.07, 0.1, 0.7, kBlue, 8, 2, kin_axis_title, 0.15, 0.15, 0.23);
      C_theta_phi_plot->Draw("AP");
      //-------------------------------------------------------------------------------------//
      
      cov_mat_elements_vs_kin_can->Print(pdfFileName); // Append the current page to the PDF
      
    }
  }
  cov_mat_elements_vs_kin_can->Print(pdfFileName + "]");
  std::cout << "Closing output pdf file." << std::endl;
  delete cov_mat_elements_vs_kin_can;

  if (cutTree) {
    delete cutTree;
    cutTree = nullptr;
  }

  return 0;
}
int main() {
  return covMatrix_plots();
}
