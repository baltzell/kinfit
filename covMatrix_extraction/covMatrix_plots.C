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
#include "TGraph.h"


const std::vector<TString> KINES = {"P", "#theta", "#phi"};
const std::vector<TString> UNITS = {"[GeV]", "[rad]", "[rad]"};
int covMatrix_plots()
{
  //Get covariance matrix values from file
  int N_lines = 0;
  FILE *cov_inFile = fopen("/work/clas12/reedtg/clas12_kinematic_fitter/dev_trevor/cov_matrix_txt_files/matrix_elements_pip_sec2.txt", "r");
  int ch;
  
  //Count number of lines in input file
  while (EOF != (ch=getc(cov_inFile))) {
    if ('\n' == ch) {
      ++N_lines;
    }
  }
  
  double C_P_temp;
  double C_theta_temp;
  double C_phi_temp;
  double C_P_theta_temp;
  double C_P_phi_temp;
  double C_theta_phi_temp;
  double P_bin_center_temp;
  double theta_bin_center_temp;
  double phi_bin_center_temp;
  
  double C_P[N_lines];
  double C_theta[N_lines];
  double C_phi[N_lines];
  double C_P_theta[N_lines];
  double C_P_phi[N_lines];
  double C_theta_phi[N_lines];
  double P_bin_center[N_lines];
  double theta_bin_center[N_lines];
  double phi_bin_center[N_lines];
  /*
  TVectorD C_P;
  std::vector<double> C_theta;
  std::vector<double> C_phi;
  std::vector<double> C_P_theta;
  std::vector<double> C_P_phi;
  std::vector<double> C_theta_phi;
  std::vector<double> P_bin_center;
  std::vector<double> theta_bin_center;
  std::vector<double> phi_bin_center;
  */
  int bufferLength = 100;
  char buffer[bufferLength];
  int line_count = 0;
  
  //Fill arrays with info from covariance matrix txt file
  rewind(cov_inFile);      //Move the cursor back to the beginning of the file in case it's already been read
  while(fgets(buffer, bufferLength, cov_inFile)) {
    sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &C_P_temp, &C_theta_temp, &C_phi_temp, &C_P_phi_temp, &C_P_theta_temp, &C_theta_phi_temp, &P_bin_center_temp, &theta_bin_center_temp, &phi_bin_center_temp);
    
    C_P[line_count] = C_P_temp;
    C_theta[line_count] = C_theta_temp;
    C_phi[line_count] = C_phi_temp;
    C_P_phi[line_count] = C_P_phi_temp;
    C_P_theta[line_count] = C_P_theta_temp;
    C_theta_phi[line_count] = C_theta_phi_temp;
    P_bin_center[line_count] = P_bin_center_temp;
    theta_bin_center[line_count] = theta_bin_center_temp;
    phi_bin_center[line_count] = phi_bin_center_temp;
    /*
    C_P.push_back(C_P_temp);   
    C_theta.push_back(C_theta_temp);      
    C_phi.push_back(C_phi_temp);      
    C_P_phi.push_back(C_P_phi_temp);   
    C_P_theta.push_back(C_P_theta_temp);  
    C_theta_phi.push_back(C_theta_phi_temp);   
    P_bin_center.push_back(P_bin_center_temp);   
    theta_bin_center.push_back(theta_bin_center_temp);   
    phi_bin_center.push_back(phi_bin_center_temp);
    */
    line_count++;
  }
  fclose(cov_inFile);
  
  
  //Plots
  TString part = "#pi^{+}";
  TString varying_kin = "P";
  //TString varying_kin = "#theta";
  /*
  std::vector<TString> kin_var = {"P", "#theta", "#phi"};
  std::vector<TString> C_elements = {"C_{P}", "C_{#theta}", "C_{#phi}", "C_{P #phi}", "C_{P #theta}", "C_{#theta #phi}"};
  std::vector<std::vector<double>> C_vec = {C_P, C_theta, C_phi, C_P_phi, C_P_theta, C_theta_phi};
  //C_vec.push_back(C_P);
  
  TCanvas *cov_mat_elements_can = new TCanvas("cov_mat_elements_can", "Covariance Matrix Elements", 800, 800);
  cov_mat_elements_can->Divide(2,3);
  
  std::vector<double> bin_center = P_bin_center;
  for (int i = 0; i < C_elements.size(); i++) {
    TGraph *plot = new TGraph(N_lines, P_bin_center, C_vec[i]);
    plot->GetXaxis()->SetTitleSize(0.07);
    plot->GetYaxis()->SetTitleSize(0.07);
    
    cov_mat_elements_can->cd(i+1);
    gPad->SetBottomMargin(0.15);                                                                                                                              
    gPad->SetLeftMargin(0.15);
    plot->Draw();
  }
  */

  TCanvas *cov_mat_elements_vs_P_can = new TCanvas("cov_mat_elements_vs_P_can", "Covariance Matrix Elements as Function of P", 800, 800);
  cov_mat_elements_vs_P_can->Divide(2,3);

  //map<TString, vector<int> > kin_map;
  double bin_center[N_lines];
  double C_P_div_P[N_lines];
  double C_theta_times_P[N_lines];
  if (varying_kin == "P") { 
    for(int j = 0; j < N_lines; ++j) {
      bin_center[j] = P_bin_center[j];
      C_P_div_P[j] = C_P[j] / P_bin_center[j];
      C_theta_times_P[j] = C_theta[j] * theta_bin_center[j];
    }
  }
  else if (varying_kin == "#theta") {
    for(int j = 0; j < N_lines; ++j) {
      bin_center[j] = theta_bin_center[j];
      C_P_div_P[j] = C_P[j] / P_bin_center[j];
      C_theta_times_P[j] = C_theta[j] * theta_bin_center[j];
    }
  }
  else if (varying_kin == "#phi") {
    for(int j = 0; j < N_lines; ++j) {
      bin_center[j] = phi_bin_center[j];
      C_P_div_P[j] = C_P[j] / P_bin_center[j];
      C_theta_times_P[j] = C_theta[j] * theta_bin_center[j];
    }
  }

  //TGraph *C_P_plot = new TGraph(N_lines, bin_center, C_P_div_P);
  //TGraph *C_theta_plot = new TGraph(N_lines, bin_center, C_theta_times_P);
  TGraph *C_P_plot = new TGraph(N_lines, bin_center, C_P);
  TGraph *C_theta_plot = new TGraph(N_lines, bin_center, C_theta);
  TGraph *C_phi_plot = new TGraph(N_lines, bin_center, C_phi);
  TGraph *C_P_phi_plot = new TGraph(N_lines, bin_center, C_P_phi);
  TGraph *C_P_theta_plot = new TGraph(N_lines, bin_center, C_P_theta);
  TGraph *C_theta_phi_plot = new TGraph(N_lines, bin_center, C_theta_phi);
  //Variances
  cov_mat_elements_vs_P_can->cd(1);
  C_P_plot->SetMarkerColor(4);
  C_P_plot->SetMarkerStyle(8);
  C_P_plot->SetMarkerSize(0.5);
  C_P_plot->SetTitle("C_{P}");
  C_P_plot->GetYaxis()->SetTitle("C_{P}");
  C_P_plot->GetYaxis()->SetTitleSize(0.07);
  C_P_plot->GetYaxis()->SetLabelSize(0.07);
  C_P_plot->GetYaxis()->SetMaxDigits(2);
  C_P_plot->GetXaxis()->SetTitle(varying_kin);
  C_P_plot->GetXaxis()->SetTitleSize(0.07);
  C_P_plot->GetXaxis()->SetLabelSize(0.07);
  //C_P_plot->GetHistogram()->SetMaximum(0.0055);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.15);
  gPad->SetLeftMargin(0.2);
  C_P_plot->Draw("AP");

  cov_mat_elements_vs_P_can->cd(3);
  C_theta_plot->SetMarkerColor(4);
  C_theta_plot->SetMarkerStyle(8);
  C_theta_plot->SetMarkerSize(0.5);
  C_theta_plot->SetTitle("C_{#theta}");
  C_theta_plot->GetYaxis()->SetTitle("C_{#theta}");
  C_theta_plot->GetYaxis()->SetTitleSize(0.07);
  C_theta_plot->GetYaxis()->SetLabelSize(0.07);
  C_theta_plot->GetYaxis()->SetMaxDigits(2);
  C_theta_plot->GetXaxis()->SetTitle(varying_kin);
  C_theta_plot->GetXaxis()->SetTitleSize(0.07);
  C_theta_plot->GetXaxis()->SetLabelSize(0.07);
  //C_theta_plot->GetHistogram()->SetMaximum(0.000007);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.15);
  gPad->SetLeftMargin(0.2);
  C_theta_plot->Draw("AP");

  cov_mat_elements_vs_P_can->cd(5);
  C_phi_plot->SetMarkerColor(4);
  C_phi_plot->SetMarkerStyle(8);
  C_phi_plot->SetMarkerSize(0.5);
  C_phi_plot->SetTitle("C_{#phi}");
  C_phi_plot->GetYaxis()->SetTitle("C_{#phi}");
  C_phi_plot->GetYaxis()->SetTitleSize(0.07);
  C_phi_plot->GetYaxis()->SetLabelSize(0.07);
  C_phi_plot->GetYaxis()->SetMaxDigits(2);
  C_phi_plot->GetXaxis()->SetTitle(varying_kin);
  C_phi_plot->GetXaxis()->SetTitleSize(0.07);
  C_phi_plot->GetXaxis()->SetLabelSize(0.07);
  //C_phi_plot->GetHistogram()->SetMaximum(0.00003);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.15);
  gPad->SetLeftMargin(0.2);
  C_phi_plot->Draw("AP");

  //Covariances
  cov_mat_elements_vs_P_can->cd(2);
  C_P_phi_plot->SetMarkerColor(4);
  C_P_phi_plot->SetMarkerStyle(8);
  C_P_phi_plot->SetMarkerSize(0.5);
  C_P_phi_plot->SetTitle("C_{P#phi}");
  C_P_phi_plot->GetYaxis()->SetTitle("C_{P#phi}");
  C_P_phi_plot->GetYaxis()->SetTitleSize(0.07);
  C_P_phi_plot->GetYaxis()->SetLabelSize(0.07);
  C_P_phi_plot->GetYaxis()->SetMaxDigits(2);
  C_P_phi_plot->GetXaxis()->SetTitle(varying_kin);
  C_P_phi_plot->GetXaxis()->SetTitleSize(0.07);
  C_P_phi_plot->GetXaxis()->SetLabelSize(0.07);
  //C_P_phi_plot->GetHistogram()->SetMaximum(0.00004);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.15);
  gPad->SetLeftMargin(0.2);
  C_P_phi_plot->Draw("AP");

  cov_mat_elements_vs_P_can->cd(4);
  C_P_theta_plot->SetMarkerColor(4);
  C_P_theta_plot->SetMarkerStyle(8);
  C_P_theta_plot->SetMarkerSize(0.5);
  C_P_theta_plot->SetTitle("C_{P#theta}");
  C_P_theta_plot->GetYaxis()->SetTitle("C_{P#theta}");
  C_P_theta_plot->GetYaxis()->SetTitleSize(0.07);
  C_P_theta_plot->GetYaxis()->SetLabelSize(0.07);
  C_P_theta_plot->GetYaxis()->SetMaxDigits(2);
  C_P_theta_plot->GetXaxis()->SetTitle(varying_kin);
  C_P_theta_plot->GetXaxis()->SetTitleSize(0.07);
  C_P_theta_plot->GetXaxis()->SetLabelSize(0.07);
  //C_P_theta_plot->GetHistogram()->SetMaximum(0.000012);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.15);
  gPad->SetLeftMargin(0.2);
  C_P_theta_plot->Draw("AP");

  cov_mat_elements_vs_P_can->cd(6);
  C_theta_phi_plot->SetMarkerColor(4);
  C_theta_phi_plot->SetMarkerStyle(8);
  C_theta_phi_plot->SetMarkerSize(0.5);
  C_theta_phi_plot->SetTitle("C_{#theta#phi}");
  C_theta_phi_plot->GetYaxis()->SetTitle("C_{#theta#phi}");
  C_theta_phi_plot->GetYaxis()->SetTitleSize(0.07);
  C_theta_phi_plot->GetYaxis()->SetLabelSize(0.07);
  C_theta_phi_plot->GetYaxis()->SetMaxDigits(2);
  C_theta_phi_plot->GetXaxis()->SetTitle(varying_kin);
  C_theta_phi_plot->GetXaxis()->SetTitleSize(0.07);
  C_theta_phi_plot->GetXaxis()->SetLabelSize(0.07);
  //C_theta_phi_plot->GetHistogram()->SetMaximum(0.000007);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.15);
  gPad->SetLeftMargin(0.2);
  C_theta_phi_plot->Draw("AP");
  
  gStyle->SetTitleFontSize(0.1);
  cov_mat_elements_vs_P_can->SaveAs(Form("plots/C_vs_%s_plots.pdf", varying_kin.Data()));
  
  return 0;
}
int main() {
  return covMatrix_plots();
}
