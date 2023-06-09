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


TGraph* tgraph_setup(TGraph* graph, double markersize, Int_t markercolor, Int_t markerstyle) {
  graph->SetMarkerColor(markercolor);
  graph->SetMarkerStyle(markerstyle);
  graph->SetMarkerSize(markersize);

  return graph;
}

TMultiGraph* tmultigraph_setup(TMultiGraph* graph, const char* title, double titlesize, double maintitlesize, Int_t maxdigits, TString kin_axis_title,double bottommargin, double topmargin, double leftmargin) {
  graph->SetTitle(title);
  graph->GetYaxis()->SetTitle(title);
  graph->GetYaxis()->SetTitleSize(titlesize);
  graph->GetYaxis()->SetLabelSize(titlesize);
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
  //FILE *cov_inFile = fopen("/work/clas12/reedtg/clas12_kinematic_fitter/dev_trevor/cov_matrix_txt_files/matrix_elements_pip_sec2.txt", "r");
  //FILE *cov_inFile = fopen("/work/clas12/reedtg/clas12_kinematic_fitter/dev_trevor/cov_matrix_txt_files/matrix_elements_pip_sec2_ver2.txt", "r");
  //FILE *cov_inFile = fopen("/work/clas12/reedtg/clas12_kinematic_fitter/updated_4-16-23/kinfit/covMatrix_extraction/cov_matrix_txt_files/matrix_elements_pip_sec2_ver2.txt", "r");
  //FILE *cov_inFile = fopen("/work/clas12/reedtg/clas12_kinematic_fitter/dev_trevor/cov_matrix_txt_files/matrix_elements_pip_sec2_fixed_P_and_theta_25_phi_bins.txt", "r");
  //FILE *cov_inFile = fopen("/work/clas12/reedtg/clas12_kinematic_fitter/dev_trevor/cov_matrix_txt_files/matrix_elements_pip_sec2_fixed_P_and_phi_25_theta_bins.txt", "r");
  //FILE *cov_inFile = fopen("/work/clas12/reedtg/clas12_kinematic_fitter/dev_trevor/cov_matrix_txt_files/matrix_elements_pip_sec2_fixed_angle_25_P_bins.txt", "r");
  //const char* filename = "/work/clas12/reedtg/clas12_kinematic_fitter/updated_4-16-23/kinfit/covMatrix_extraction/covariances.root";
  const char* filename = "/work/clas12/reedtg/clas12_kinematic_fitter/dev_trevor/covariances.root";
  //const char* filename = "/work/clas12/reedtg/clas12_kinematic_fitter/updated_4-16-23/kinfit/covMatrix_extraction/covariances.root";
  TFile *cov_inFile = TFile::Open(filename); 

  // Get the TTree from the file
  TTree* tree = dynamic_cast<TTree*>(cov_inFile->Get("covarianceTree")); 
  if (!tree) {
    std::cerr << "Error retrieving TTree from file: " << cov_inFile << std::endl;
    cov_inFile->Close();
  }

  //Choose whether or not error of covariance values should be included. (Only set to true if these were calculated along with the covariances)
  bool errors = false;

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
    if (entry < 100) {
      //std::cout << "events = " << event_count << std::endl;
      //std::cout << "P_bin_center = " << P_bin_center << std::endl;
    }
    P_bin_vec.push_back(P_bin_center);
    theta_bin_vec.push_back(theta_bin_center);
    phi_bin_vec.push_back(phi_bin_center);
    C_P_vec.push_back(C_P);
    C_theta_vec.push_back(C_theta);
    C_phi_vec.push_back(C_phi);
    C_P_phi_vec.push_back(C_P_phi);
    C_P_theta_vec.push_back(C_P_theta);
    C_theta_phi_vec.push_back(C_theta_phi);
    event_count_vec.push_back(event_count);
    if (errors) {
      C_P_err_vec.push_back(C_P_err);
      C_theta_err_vec.push_back(C_theta_err);
      C_phi_err_vec.push_back(C_phi_err);
      C_P_phi_err_vec.push_back(C_P_phi_err);
      C_P_theta_err_vec.push_back(C_P_theta_err);
      C_theta_phi_err_vec.push_back(C_theta_phi_err);
    }
  }

  int N_fixed_variables = 2;
  std::string first_fixed_var;

  double fixed_P_bin;
  double fixed_theta_bin;
  double fixed_phi_bin;

  std::set<double> unique_P_bins(P_bin_vec.begin(), P_bin_vec.end());
  std::set<double> unique_P_bins_rounded;
  for (double value : unique_P_bins) {
    double roundedValue = std::round(value * 100.0) / 100.0;;
    unique_P_bins_rounded.insert(roundedValue);
  }
  //  for (int j = 0; j < unique_P_bins.size(); j++) { 
  //  unique_P_bins[j] = std::round(unique_P_bins[j] * 100.0) / 100.0;
  //}
  std::set<double> unique_theta_bins(theta_bin_vec.begin(), theta_bin_vec.end());
  std::set<double> unique_phi_bins(phi_bin_vec.begin(), phi_bin_vec.end());
  //std::cout << "Number of P bins = " << unique_P_bins.size() << std::endl;

  std::string response1;
  std::string response2;
  std::string response3;
  std::cout << "Do you want to fix two of the kinematic variables and plot the covariance values as a function of the third? Enter 'yes' or 'no'." << std::endl;
  std::cin >> response1;

  bool valid_input = false;
  // Process the response
  if (response1 == "yes" || response1 == "y") {
    for (int i = 0; i < N_fixed_variables; i++) {
      if (i == 0) {
	std::cout << "Select a variable to fix. Type 'p', 'theta', or 'phi'." << std::endl;
	std::cin >> response2;
      }
      else if (first_fixed_var == "p" || first_fixed_var == "theta" || first_fixed_var == "phi") {
	std::cout << std::endl;
	std::cout << "Select another variable to fix." << std::endl;
	std::cin >> response2;
	valid_input = true;
      }
      else {
	std::cout << "Not a valid response. Exiting" << std::endl;
	break;
      }
      if (response2 == "p") {
	if (i >= 1 && first_fixed_var == "p") {
	  std::cout << "Momentum was already fixed. Exiting." << std::endl;
	  valid_input = false;
	  break;
	}
	else {
	  int count = 0;
	  std::cout << "These are the " << unique_P_bins.size() << " P bin centers: " << std::endl;
	  for (auto element : unique_P_bins) {
	    count++;
	    std::cout << count << "-" << element << "  ";
	  }
	  std::cout << std::endl;
	  std::cout << "Select theta P value you want to use. Type in the index number." << std::endl;
	  std::cin >> response3;

	  int response3_int = std::stoi(response3);
	  if ((response3_int > 0) && (response3_int <= unique_P_bins.size())) { 
	    int P_index = response3_int - 1;
	    std::set<double>::iterator it = unique_P_bins.begin();
	    std::advance(it, P_index);
	    fixed_P_bin = *it;
	    std::cout << std::endl;
	    std::cout << "Ok. Fixing the momentum at " << fixed_P_bin << " GeV." << std::endl;
            if (i == 0) {
              first_fixed_var = "p";
            }
          }
	  else {
	    std::cout << "Invalid response. Exiting." << std::endl;
	    valid_input = false;
	  }
	}
      }
      else if (response2 == "theta") {
	if (i >= 1 && first_fixed_var == "theta") {
	  std::cout << "Theta was already fixed. Exiting." <<std::endl;
	  break;
	}
	else {
	  int count = 0;
	  std::cout << "These are the " << unique_theta_bins.size() << " theta bin centers: " << std::endl;
	  for (auto element : unique_theta_bins) {
	    count++;
	    std::cout << count << "-" << element << "  ";
	  }
	  std::cout << std::endl;
	  std::cout << "Select the theta value you want to use. Type in the index number." << std::endl;
	  std::cin >> response3;
	  int response3_int = std::stoi(response3);
          if ((response3_int > 0) && (response3_int <= unique_theta_bins.size())) {
	    int theta_index = response3_int - 1;
	    std::set<double>::iterator it = unique_theta_bins.begin();
	    std::advance(it, theta_index);
            fixed_theta_bin = *it;
	    std::cout << std::endl;
	    std::cout << "Ok. Fixing theta at " << fixed_theta_bin << " degrees." << std::endl;
	    if (i == 0) {
	      first_fixed_var = "theta";
	    }
	  }
	  else {
	    std::cout << "Invalid response. Exiting." << std::endl;
	  } 
	}
      }
      else if (response2 == "phi") {
	if (i >= 1 && first_fixed_var == "phi") {
	  std::cout << "Phi was already fixed. Exiting." <<std::endl;
	  valid_input = false;
	  break;
	}
	else {
	  int count = 0;
	  std::cout << "These are the " << unique_phi_bins.size() << " phi bin centers: " << std::endl;
	  for (auto element : unique_phi_bins) {
	    count++;
	    std::cout << count << "-" << element << "  ";
	  }
	  std::cout << std::endl;
	  std::cout << "Select the phi value you want to use. Type in the index number." << std::endl;
	  std::cin >> response3;
	  int response3_int = std::stoi(response3);
          if ((response3_int > 0) && (response3_int <= unique_phi_bins.size())) {  
	    int phi_index = response3_int - 1;
	    std::set<double>::iterator it = unique_phi_bins.begin();
	    std::advance(it, phi_index);
            fixed_phi_bin = *it;
	    std::cout << std::endl;
	    std::cout << "Ok. Fixing phi at " << fixed_phi_bin << " degrees." << std::endl;
	    if (i == 0) {
	      first_fixed_var = "phi";
	    }
	  }
	  else {
	    std::cout << "Invalid response. Exiting." << std::endl;
	    valid_input = false;
	  }
	}
      }
      else {
	std::cout << "Invalid response. Exiting." << std::endl;
	valid_input = false;
      }
    }
  } else if (response1 == "no" || response1 == "n") {
    std::cout << "Ok. Nothing else will be done then." << std::endl;
  } else {
    std::cout << "Invalid response. Exiting" << std::endl;
  }


  //If valid responses were given for fixing two of the kinematic variables, the TTree will be cut and scatterplots will be made. 
  if (valid_input == true) {
    // Create a new TTree with the cut data
    TTree* cutTree;
    TString part = "#pi^{+}";
    TString varying_kin;
    if (fixed_P_bin < 0.000001) {
      std::cout << "Fixed theta bin = " << fixed_theta_bin << " degrees, Fixed phi bin = " << fixed_phi_bin << " degrees" << std::endl;
      varying_kin = "P";
      char cutExpression[100];  // Buffer to hold the theta cut expression
      sprintf(cutExpression, "(theta > %f && theta < %f) && (phi > %f && phi < %f)", fixed_theta_bin-0.01, fixed_theta_bin+0.01, fixed_phi_bin-0.01, fixed_phi_bin+0.01);
      cutTree = tree->CopyTree(cutExpression);
      //cutTree->Print();
    }
    else if (fixed_theta_bin < 0.000001) { 
      std::cout << "Fixed p bin = " << fixed_P_bin << " GeV, Fixed phi bin = " << fixed_phi_bin << " degrees" << std::endl;
      varying_kin = "theta";
      char cutExpression[100];
      sprintf(cutExpression, "(P > %f && P < %f) && (phi > %f && phi < %f)", fixed_P_bin-0.001, fixed_P_bin+0.001, fixed_phi_bin-0.01, fixed_phi_bin+0.01);
      cutTree = tree->CopyTree(cutExpression);
      //cutTree->Print();
    }
    else if (fixed_phi_bin < 0.000001) {
      std::cout << "Fixed p bin = " << fixed_P_bin << " GeV, Fixed theta bin = " << fixed_theta_bin << " degrees" << std::endl;
      varying_kin = "phi";
      char cutExpression[100];
      sprintf(cutExpression, "(P > %f && P < %f) && (theta > %f && theta < %f)", fixed_P_bin-0.001, fixed_P_bin+0.001, fixed_theta_bin-0.01, fixed_theta_bin+0.01);
      cutTree = tree->CopyTree(cutExpression);
      //cutTree->Print();
    }
    
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
	if (varying_kin == "P") {
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
    if (varying_kin == "P") {
      C_P_fixed_kin_1_bin = C_P_hist->GetYaxis()->FindBin(fixed_theta_bin);
      C_P_fixed_kin_2_bin = C_P_hist->GetZaxis()->FindBin(fixed_phi_bin);
      C_phi_fixed_kin_1_bin = C_phi_hist->GetYaxis()->FindBin(fixed_theta_bin);
      C_phi_fixed_kin_2_bin = C_phi_hist->GetZaxis()->FindBin(fixed_phi_bin);
      C_theta_fixed_kin_1_bin = C_theta_hist->GetYaxis()->FindBin(fixed_theta_bin);
      C_theta_fixed_kin_2_bin = C_theta_hist->GetZaxis()->FindBin(fixed_phi_bin);
      C_P_phi_fixed_kin_1_bin = C_P_phi_hist->GetYaxis()->FindBin(fixed_theta_bin);
      C_P_phi_fixed_kin_2_bin = C_P_phi_hist->GetZaxis()->FindBin(fixed_phi_bin);
      C_P_theta_fixed_kin_1_bin = C_P_theta_hist->GetYaxis()->FindBin(fixed_theta_bin);
      C_P_theta_fixed_kin_2_bin = C_P_theta_hist->GetZaxis()->FindBin(fixed_phi_bin);
      C_theta_phi_fixed_kin_1_bin = C_theta_phi_hist->GetYaxis()->FindBin(fixed_theta_bin);
      C_theta_phi_fixed_kin_2_bin = C_theta_phi_hist->GetZaxis()->FindBin(fixed_phi_bin);
    }
    else if (varying_kin == "theta") {
      C_P_fixed_kin_1_bin = C_P_hist->GetXaxis()->FindBin(fixed_P_bin);
      C_P_fixed_kin_2_bin = C_P_hist->GetZaxis()->FindBin(fixed_phi_bin);
      C_phi_fixed_kin_1_bin = C_phi_hist->GetXaxis()->FindBin(fixed_P_bin);
      C_phi_fixed_kin_2_bin = C_phi_hist->GetZaxis()->FindBin(fixed_phi_bin);
      C_theta_fixed_kin_1_bin = C_theta_hist->GetXaxis()->FindBin(fixed_P_bin);
      C_theta_fixed_kin_2_bin = C_theta_hist->GetZaxis()->FindBin(fixed_phi_bin);
      C_P_phi_fixed_kin_1_bin = C_P_phi_hist->GetXaxis()->FindBin(fixed_P_bin);
      C_P_phi_fixed_kin_2_bin = C_P_phi_hist->GetZaxis()->FindBin(fixed_phi_bin);
      C_P_theta_fixed_kin_1_bin = C_P_theta_hist->GetXaxis()->FindBin(fixed_P_bin);
      C_P_theta_fixed_kin_2_bin = C_P_theta_hist->GetZaxis()->FindBin(fixed_phi_bin);
      C_theta_phi_fixed_kin_1_bin = C_theta_phi_hist->GetXaxis()->FindBin(fixed_P_bin);
      C_theta_phi_fixed_kin_2_bin = C_theta_phi_hist->GetZaxis()->FindBin(fixed_phi_bin);
    }
    else if (varying_kin == "phi") {
      C_P_fixed_kin_1_bin = C_P_hist->GetXaxis()->FindBin(fixed_P_bin);
      C_P_fixed_kin_2_bin = C_P_hist->GetYaxis()->FindBin(fixed_theta_bin);
      C_phi_fixed_kin_1_bin = C_phi_hist->GetXaxis()->FindBin(fixed_P_bin);
      C_phi_fixed_kin_2_bin = C_phi_hist->GetYaxis()->FindBin(fixed_theta_bin);
      C_theta_fixed_kin_1_bin = C_theta_hist->GetXaxis()->FindBin(fixed_P_bin);
      C_theta_fixed_kin_2_bin = C_theta_hist->GetYaxis()->FindBin(fixed_theta_bin);
      C_P_phi_fixed_kin_1_bin = C_P_phi_hist->GetXaxis()->FindBin(fixed_P_bin);
      C_P_phi_fixed_kin_2_bin = C_P_phi_hist->GetYaxis()->FindBin(fixed_theta_bin);
      C_P_theta_fixed_kin_1_bin = C_P_theta_hist->GetXaxis()->FindBin(fixed_P_bin);
      C_P_theta_fixed_kin_2_bin = C_P_theta_hist->GetYaxis()->FindBin(fixed_theta_bin);
      C_theta_phi_fixed_kin_1_bin = C_theta_phi_hist->GetXaxis()->FindBin(fixed_P_bin);
      C_theta_phi_fixed_kin_2_bin = C_theta_phi_hist->GetYaxis()->FindBin(fixed_theta_bin);
    }

    TCanvas *cov_mat_elements_vs_kin_can = new TCanvas("cov_mat_elements_vs_kin_can", "Covariance Matrix Elements as Function of Kinematics", 800, 800);
    cov_mat_elements_vs_kin_can->Divide(2,3);
    //const Int_t numInterpolatedPoints = 2;
    //-------------------------Set up the covariance scatterplots------------------------// 
    TGraph *C_P_plot;
    TGraph *C_theta_plot;
    TGraph *C_phi_plot;
    TGraph *C_P_theta_plot;
    TGraph *C_P_phi_plot;
    TGraph *C_theta_phi_plot;
    /*
    TGraphErrors *C_P_plot;
    TGraphErrors *C_theta_plot;
    TGraphErrors *C_phi_plot;
    TGraphErrors *C_P_theta_plot;
    TGraphErrors *C_P_phi_plot;
    TGraphErrors *C_theta_phi_plot;
    */
    //Interpolation TGraphs
    TGraph *C_P_plot_interpolate = new TGraph();
    TGraph *C_phi_plot_interpolate = new TGraph();
    TGraph *C_theta_plot_interpolate = new TGraph();
    TGraph *C_P_phi_plot_interpolate = new TGraph();
    TGraph *C_P_theta_plot_interpolate = new TGraph();
    TGraph *C_theta_phi_plot_interpolate = new TGraph();
    //TMultiGraphs: The tgraphs of the initial data points and the interpolated points will both be added to these
    TMultiGraph* C_P_with_interpolation_plot = new TMultiGraph();
    TMultiGraph* C_phi_with_interpolation_plot = new TMultiGraph();
    TMultiGraph* C_theta_with_interpolation_plot = new TMultiGraph();
    TMultiGraph* C_P_phi_with_interpolation_plot = new TMultiGraph();
    TMultiGraph* C_P_theta_with_interpolation_plot = new TMultiGraph();
    TMultiGraph* C_theta_phi_with_interpolation_plot = new TMultiGraph();
    //If using error bars
    if (errors) {
      C_P_plot = new TGraphErrors(numCutEntries, bin_center_vec.data(), C_P_cut_vec.data(), nullptr, C_P_err_cut_vec.data());
      C_theta_plot = new TGraphErrors(numCutEntries, bin_center_vec.data(), C_theta_cut_vec.data(), nullptr, C_theta_err_cut_vec.data());
      C_phi_plot = new TGraphErrors(numCutEntries, bin_center_vec.data(), C_phi_cut_vec.data(), nullptr, C_phi_err_cut_vec.data());
      C_P_phi_plot = new TGraphErrors(numCutEntries, bin_center_vec.data(), C_P_phi_cut_vec.data(), nullptr, C_P_phi_err_cut_vec.data());
      C_P_theta_plot = new TGraphErrors(numCutEntries, bin_center_vec.data(), C_P_theta_cut_vec.data(), nullptr, C_P_theta_err_cut_vec.data());
      C_theta_phi_plot = new TGraphErrors(numCutEntries, bin_center_vec.data(), C_theta_phi_cut_vec.data(), nullptr, C_theta_phi_err_cut_vec.data());
    }
    //If NOT using error bars
    else {
      C_P_plot = new TGraph(passed_min_events, bin_center_vec.data(), C_P_cut_vec.data());
      C_theta_plot = new TGraph(passed_min_events, bin_center_vec.data(), C_theta_cut_vec.data());
      C_phi_plot = new TGraph(passed_min_events, bin_center_vec.data(), C_phi_cut_vec.data());
      C_P_phi_plot = new TGraph(passed_min_events, bin_center_vec.data(), C_P_phi_cut_vec.data());
      C_P_theta_plot = new TGraph(passed_min_events, bin_center_vec.data(), C_P_theta_cut_vec.data());
      C_theta_phi_plot = new TGraph(passed_min_events, bin_center_vec.data(), C_theta_phi_cut_vec.data());
    }
    //-----------------------------------------------------------------------------------//

    int interpolate_index = 0;
    const Int_t numInterpolatedPoints = 10;  //Number of points to be interpolated in between each data point
    //const Int_t numPoints = C_P_plot->GetN();
    for (Int_t i = 0; i < passed_min_events - 1; ++i) {

      //Get the data points from the covariance plots (only the x values are used in the interpolation method)
      Double_t C_P_x1, C_P_x2, C_P_y1, C_P_y2;
      C_P_plot->GetPoint(i, C_P_x1, C_P_y1);            //Get the x (kin) and y (C) values of a data point 
      C_P_plot->GetPoint(i + 1, C_P_x2, C_P_y2);        //Get the x and y values of the next data point
      Double_t C_phi_x1, C_phi_x2, C_phi_y1, C_phi_y2;
      C_phi_plot->GetPoint(i, C_phi_x1, C_phi_y1);
      C_phi_plot->GetPoint(i + 1, C_phi_x2, C_phi_y2);
      Double_t C_theta_x1, C_theta_x2, C_theta_y1, C_theta_y2;
      C_theta_plot->GetPoint(i, C_theta_x1, C_theta_y1);
      C_theta_plot->GetPoint(i + 1, C_theta_x2, C_theta_y2);
      Double_t C_P_phi_x1, C_P_phi_x2, C_P_phi_y1, C_P_phi_y2;
      C_P_phi_plot->GetPoint(i, C_P_phi_x1, C_P_phi_y1);
      C_P_phi_plot->GetPoint(i + 1, C_P_phi_x2, C_P_phi_y2);
      Double_t C_P_theta_x1, C_P_theta_x2, C_P_theta_y1, C_P_theta_y2;
      C_P_theta_plot->GetPoint(i, C_P_theta_x1, C_P_theta_y1);
      C_P_theta_plot->GetPoint(i + 1, C_P_theta_x2, C_P_theta_y2);
      Double_t C_theta_phi_x1, C_theta_phi_x2, C_theta_phi_y1, C_theta_phi_y2;
      C_theta_phi_plot->GetPoint(i, C_theta_phi_x1, C_theta_phi_y1);
      C_theta_phi_plot->GetPoint(i + 1, C_theta_phi_x2, C_theta_phi_y2);

      //Get the interpolated values in between the data points 
      for (Int_t j = 1; j <= numInterpolatedPoints; ++j) {
	Double_t C_P_xInterpolated = C_P_x1 + (C_P_x2 - C_P_x1) * (j / static_cast<Double_t>(numInterpolatedPoints + 1));
	Double_t C_phi_xInterpolated = C_phi_x1 + (C_phi_x2 - C_phi_x1) * (j / static_cast<Double_t>(numInterpolatedPoints + 1));
	Double_t C_theta_xInterpolated = C_theta_x1 + (C_theta_x2 - C_theta_x1) * (j / static_cast<Double_t>(numInterpolatedPoints + 1));
	Double_t C_P_phi_xInterpolated = C_P_phi_x1 + (C_P_phi_x2 - C_P_phi_x1) * (j / static_cast<Double_t>(numInterpolatedPoints + 1));
	Double_t C_P_theta_xInterpolated = C_P_theta_x1 + (C_P_theta_x2 - C_P_theta_x1) * (j / static_cast<Double_t>(numInterpolatedPoints + 1));
	Double_t C_theta_phi_xInterpolated = C_theta_phi_x1 + (C_theta_phi_x2 - C_theta_phi_x1) * (j / static_cast<Double_t>(numInterpolatedPoints + 1));
	Double_t C_P_yInterpolated;
	Double_t C_phi_yInterpolated;
	Double_t C_theta_yInterpolated;
	Double_t C_P_phi_yInterpolated;
	Double_t C_P_theta_yInterpolated;
	Double_t C_theta_phi_yInterpolated;
	if (varying_kin == "P") {
	  C_P_yInterpolated = C_P_hist->Interpolate(C_P_xInterpolated, fixed_theta_bin, fixed_phi_bin);
	  C_phi_yInterpolated = C_phi_hist->Interpolate(C_phi_xInterpolated, fixed_theta_bin, fixed_phi_bin);
	  C_theta_yInterpolated = C_theta_hist->Interpolate(C_theta_xInterpolated, fixed_theta_bin, fixed_phi_bin);
	  C_P_phi_yInterpolated = C_P_phi_hist->Interpolate(C_P_phi_xInterpolated, fixed_theta_bin, fixed_phi_bin);
	  C_P_theta_yInterpolated = C_P_theta_hist->Interpolate(C_P_theta_xInterpolated, fixed_theta_bin, fixed_phi_bin);
	  C_theta_phi_yInterpolated = C_theta_phi_hist->Interpolate(C_theta_phi_xInterpolated, fixed_theta_bin, fixed_phi_bin);
	}
	else if (varying_kin == "theta") {
          C_P_yInterpolated = C_P_hist->Interpolate(fixed_P_bin, C_P_xInterpolated, fixed_phi_bin);
	  C_phi_yInterpolated = C_phi_hist->Interpolate(fixed_P_bin, C_phi_xInterpolated, fixed_phi_bin);
	  C_theta_yInterpolated = C_theta_hist->Interpolate(fixed_P_bin, C_theta_xInterpolated, fixed_phi_bin);
	  C_P_phi_yInterpolated = C_P_phi_hist->Interpolate(fixed_P_bin, C_P_phi_xInterpolated, fixed_phi_bin);
	  C_P_theta_yInterpolated = C_P_theta_hist->Interpolate(fixed_P_bin, C_P_theta_xInterpolated, fixed_phi_bin);
	  C_theta_phi_yInterpolated = C_theta_phi_hist->Interpolate(fixed_P_bin, C_theta_phi_xInterpolated, fixed_phi_bin);
	}
	else if (varying_kin == "phi") {
          C_P_yInterpolated = C_P_hist->Interpolate(fixed_P_bin, fixed_theta_bin, C_P_xInterpolated);
	  C_phi_yInterpolated = C_phi_hist->Interpolate(fixed_P_bin, fixed_theta_bin, C_phi_xInterpolated);
	  C_theta_yInterpolated = C_theta_hist->Interpolate(fixed_P_bin, fixed_theta_bin, C_theta_xInterpolated);
	  C_P_phi_yInterpolated = C_P_phi_hist->Interpolate(fixed_P_bin, fixed_theta_bin, C_P_phi_xInterpolated);
	  C_P_theta_yInterpolated = C_P_theta_hist->Interpolate(fixed_P_bin, fixed_theta_bin, C_P_theta_xInterpolated);
	  C_theta_phi_yInterpolated = C_theta_phi_hist->Interpolate(fixed_P_bin, fixed_theta_bin, C_theta_phi_xInterpolated);
        }
	//Add interpolated data point to tgraph
	C_P_plot_interpolate->SetPoint(interpolate_index, C_P_xInterpolated, C_P_yInterpolated);
	C_phi_plot_interpolate->SetPoint(interpolate_index, C_phi_xInterpolated, C_phi_yInterpolated);
	C_theta_plot_interpolate->SetPoint(interpolate_index, C_theta_xInterpolated, C_theta_yInterpolated);
	C_P_phi_plot_interpolate->SetPoint(interpolate_index, C_P_phi_xInterpolated, C_P_phi_yInterpolated);
	C_P_theta_plot_interpolate->SetPoint(interpolate_index, C_P_theta_xInterpolated, C_P_theta_yInterpolated);
	C_theta_phi_plot_interpolate->SetPoint(interpolate_index, C_theta_phi_xInterpolated, C_theta_phi_yInterpolated);
	interpolate_index++;
      }
    }

    //===================================Variances Plots==========================================//
    //Momentum
    cov_mat_elements_vs_kin_can->cd(1);
    C_P_plot = tgraph_setup(C_P_plot, 0.5, kBlue, 8);
    C_P_plot_interpolate = tgraph_setup(C_P_plot_interpolate, 0.25, kRed, 8);
    C_P_with_interpolation_plot->Add(C_P_plot);
    C_P_with_interpolation_plot->Add(C_P_plot_interpolate);
    C_P_with_interpolation_plot = tmultigraph_setup(C_P_with_interpolation_plot, "C_{P}", 0.07, 0.1, 2, kin_axis_title, 0.15, 0.15, 0.2);
    C_P_with_interpolation_plot->Draw("AP");

    //Phi
    cov_mat_elements_vs_kin_can->cd(3);
    C_phi_plot = tgraph_setup(C_phi_plot, 0.5, kBlue, 8);
    C_phi_plot_interpolate = tgraph_setup(C_phi_plot_interpolate, 0.25, kRed, 8);
    C_phi_with_interpolation_plot->Add(C_phi_plot);
    C_phi_with_interpolation_plot->Add(C_phi_plot_interpolate);
    C_phi_with_interpolation_plot = tmultigraph_setup(C_phi_with_interpolation_plot, "C_{#phi}", 0.07, 0.1, 2, kin_axis_title, 0.15, 0.15, 0.2);
    C_phi_with_interpolation_plot->Draw("AP");

    //Theta
    cov_mat_elements_vs_kin_can->cd(5);
    C_theta_plot = tgraph_setup(C_theta_plot, 0.5, kBlue, 8);
    C_theta_plot_interpolate = tgraph_setup(C_theta_plot_interpolate, 0.25, kRed, 8);
    C_theta_with_interpolation_plot->Add(C_theta_plot);
    C_theta_with_interpolation_plot->Add(C_theta_plot_interpolate);
    C_theta_with_interpolation_plot = tmultigraph_setup(C_theta_with_interpolation_plot, "C_{#theta}", 0.07, 0.1, 2, kin_axis_title, 0.15, 0.15, 0.2);
    C_theta_with_interpolation_plot->Draw("AP");
    //=======================================================================================//

    //----------------------------------Covariance Plots-------------------------------------//
    //Momentum and Phi
    cov_mat_elements_vs_kin_can->cd(2);
    C_P_phi_plot = tgraph_setup(C_P_phi_plot, 0.5, kBlue, 8);
    C_P_phi_plot_interpolate = tgraph_setup(C_P_phi_plot_interpolate, 0.25, kRed, 8);
    C_P_phi_with_interpolation_plot->Add(C_P_phi_plot);
    C_P_phi_with_interpolation_plot->Add(C_P_phi_plot_interpolate);
    C_P_phi_with_interpolation_plot = tmultigraph_setup(C_P_phi_with_interpolation_plot, "C_{P#phi}", 0.07, 0.1, 2, kin_axis_title, 0.15, 0.15, 0.2);
    C_P_phi_with_interpolation_plot->Draw("AP");

    //Momentum and Theta
    cov_mat_elements_vs_kin_can->cd(4);
    C_P_theta_plot = tgraph_setup(C_P_theta_plot, 0.5, kBlue, 8);
    C_P_theta_plot_interpolate = tgraph_setup(C_P_theta_plot_interpolate, 0.25, kRed, 8);
    C_P_theta_with_interpolation_plot->Add(C_P_theta_plot);
    C_P_theta_with_interpolation_plot->Add(C_P_theta_plot_interpolate);
    C_P_theta_with_interpolation_plot = tmultigraph_setup(C_P_theta_with_interpolation_plot, "C_{P#theta}", 0.07, 0.1, 2, kin_axis_title, 0.15, 0.15, 0.2);
    C_P_theta_with_interpolation_plot->Draw("AP");

    //Theta and Phi
    cov_mat_elements_vs_kin_can->cd(6);
    C_theta_phi_plot = tgraph_setup(C_theta_phi_plot, 0.5, kBlue, 8);
    C_theta_phi_plot_interpolate = tgraph_setup(C_theta_phi_plot_interpolate, 0.25, kRed, 8);
    C_theta_phi_with_interpolation_plot->Add(C_theta_phi_plot);
    C_theta_phi_with_interpolation_plot->Add(C_theta_phi_plot_interpolate);
    C_theta_phi_with_interpolation_plot = tmultigraph_setup(C_theta_phi_with_interpolation_plot, "C_{#theta#phi}", 0.07, 0.1, 2, kin_axis_title, 0.15, 0.15, 0.2);
    C_theta_phi_with_interpolation_plot->Draw("AP");
    //-------------------------------------------------------------------------------------//

    cov_mat_elements_vs_kin_can->SaveAs(Form("plots/C_vs_%s_plots.pdf", varying_kin.Data()));
  }
  return 0;
}
int main() {
  return covMatrix_plots();
}
