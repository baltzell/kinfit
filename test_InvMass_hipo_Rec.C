#include <vector>
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "KinFitter.h"
#include "KinParticle.h"
#include "KinCovariance.h"
#include "KinMomentumCorrections.h"
#include "KinFitTest.h"
#include "reader.h"
#include <dirent.h>
#include "TTree.h"
#include "TH2F.h"
#include <TROOT.h>
#include <sys/types.h>
#include <cstring>
#include <iostream>

// Create Struct for storing kinematic fitter info
struct kinFitInfoStruct
{
	double rec_p1, rec_p2, rec_theta1, rec_theta2, rec_phi1, rec_phi2, mc_p1, mc_p2, mc_theta1, mc_theta2, mc_phi1, mc_phi2;
	double fit_p1, fit_p2, fit_theta1, fit_theta2, fit_phi1, fit_phi2;
	double pull_p1, pull_p2, pull_theta1, pull_theta2, pull_phi1, pull_phi2;
	double confidence_level, cov_p1, cov_p2, cov_theta1, cov_theta2, cov_phi1, cov_phi2, cov_p_theta1, cov_p_theta2, cov_p_phi1, cov_p_phi2, cov_phi_theta1, cov_phi_theta2, convergence_status;
	double cov_err_p1, cov_err_p2, cov_err_theta1, cov_err_theta2, cov_err_phi1, cov_err_phi2, cov_err_p_theta1, cov_err_p_theta2, cov_err_p_phi1, cov_err_p_phi2, cov_err_phi_theta1, cov_err_phi_theta2;
};

// Create Struct for storing all MC events info
struct allGenStruct
{
	double all_mc_p1, all_mc_p2, all_mc_theta1, all_mc_theta2, all_mc_phi1, all_mc_phi2;
};

void read_Hipo(char input_hipo[256], std::vector<int> required_pids, std::vector<double> masses, TLorentzVector target, TLorentzVector beam, KinFitTest test, std::vector<TString> parts, TTree *kinfitTree, kinFitInfoStruct &kinFitData, TTree *allGenTree, allGenStruct &allGenData);
void read_Rec_Part_Bank(hipo::bank PartBank, hipo::bank MCMatchBank, hipo::bank TrackBank, hipo::bank TrajBank, std::vector<int> *pid_list, std::vector<TLorentzVector> *vec_list, std::vector<int> *in_part_index, std::vector<int> *in_part_sector);
void read_MC_Part_Bank(hipo::bank Bank, std::vector<int> required_pids, std::vector<int> *pid_list, std::vector<TLorentzVector> *vec_list, std::vector<int> *in_part_index);
int get_sector(int index, hipo::bank TrackBank);
int get_sector_ver2(int index, hipo::bank TrackBank);
bool pass_fid_cut_DC(int index, hipo::bank TrajBank);
bool pass_limit_cov_matrix(int sector, TLorentzVector vector, double limit_up_P, double limit_down_P, double limit_up_Theta, double limit_down_Theta, double limit_up_Phi, double limit_down_Phi);

// Utility function to generate reverse index map
std::map<int, std::vector<int>> loadMapByIndex(hipo::bank fromBank, const char idxVarName[256])
{
	std::map<int, std::vector<int>> map;
	// if (fromBank != null) {
	for (int iFrom = 0; iFrom < fromBank.getRows(); ++iFrom)
	{
		const int iTo = fromBank.getInt(idxVarName, iFrom);
		map[iTo].push_back(iFrom);
	}
	//}
	return map;
}

TH2F *pass_fid_cut = new TH2F("pass", "pass", 60, -180, 180, 60, 0, 45);
TH2F *fail_fid_cut = new TH2F("pass", "pass", 60, -180, 180, 60, 0, 45);
TH2F *P_pip_P_pim = new TH2F("PpipVSPpim", "PpipVSPpim", 60, 1.0, 11., 60, 1.0, 11.0);

std::ofstream InfoFile("event_info.txt");

int test_InvMass_hipo_Rec()
{

	gROOT->SetBatch(kTRUE);

	// Create a ROOT file to store the kinematic fit info
	TFile *kinfitRoot = new TFile("kinfit.root", "RECREATE");

	// Create a TTree to store kinematic fitting information
	TTree *kinfitTree = new TTree("kinFitInfoTree", "kinFitInfoTree");

	kinFitInfoStruct kinFitData;
	kinfitTree->Branch("rec_p1", &kinFitData.rec_p1, "rec_p1/D");
	kinfitTree->Branch("rec_p2", &kinFitData.rec_p2, "rec_p2/D");
	kinfitTree->Branch("rec_theta1", &kinFitData.rec_theta1, "rec_theta1/D");
	kinfitTree->Branch("rec_theta2", &kinFitData.rec_theta2, "rec_theta2/D");
	kinfitTree->Branch("rec_phi1", &kinFitData.rec_phi1, "rec_phi1/D");
	kinfitTree->Branch("rec_phi2", &kinFitData.rec_phi2, "rec_phi2/D");
	kinfitTree->Branch("mc_p1", &kinFitData.mc_p1, "mc_p1/D");
	kinfitTree->Branch("mc_p2", &kinFitData.mc_p2, "mc_p2/D");
	kinfitTree->Branch("mc_theta1", &kinFitData.mc_theta1, "mc_theta1/D");
	kinfitTree->Branch("mc_theta2", &kinFitData.mc_theta2, "mc_theta2/D");
	kinfitTree->Branch("mc_phi1", &kinFitData.mc_phi1, "mc_phi1/D");
	kinfitTree->Branch("mc_phi2", &kinFitData.mc_phi2, "mc_phi2/D");
	kinfitTree->Branch("fit_p1", &kinFitData.fit_p1, "fit_p1/D");
	kinfitTree->Branch("fit_p2", &kinFitData.fit_p2, "fit_p2/D");
	kinfitTree->Branch("fit_theta1", &kinFitData.fit_theta1, "fit_theta1/D");
	kinfitTree->Branch("fit_theta2", &kinFitData.fit_theta2, "fit_theta2/D");
	kinfitTree->Branch("fit_phi1", &kinFitData.fit_phi1, "fit_phi1/D");
	kinfitTree->Branch("fit_phi2", &kinFitData.fit_phi2, "fit_phi2/D");
	kinfitTree->Branch("pull_p1", &kinFitData.pull_p1, "pull_p1/D");
	kinfitTree->Branch("pull_p2", &kinFitData.pull_p2, "pull_p2/D");
	kinfitTree->Branch("pull_theta1", &kinFitData.pull_theta1, "pull_theta1/D");
	kinfitTree->Branch("pull_theta2", &kinFitData.pull_theta2, "pull_theta2/D");
	kinfitTree->Branch("pull_phi1", &kinFitData.pull_phi1, "pull_phi1/D");
	kinfitTree->Branch("pull_phi2", &kinFitData.pull_phi2, "pull_phi2/D");
	kinfitTree->Branch("confidence_level", &kinFitData.confidence_level, "confidence_level/D");
	kinfitTree->Branch("cov_p1", &kinFitData.cov_p1, "cov_p1/D");
	kinfitTree->Branch("cov_p2", &kinFitData.cov_p2, "cov_p2/D");
	kinfitTree->Branch("cov_theta1", &kinFitData.cov_theta1, "cov_theta1/D");
	kinfitTree->Branch("cov_theta2", &kinFitData.cov_theta2, "cov_theta2/D");
	kinfitTree->Branch("cov_phi1", &kinFitData.cov_phi1, "cov_phi1/D");
	kinfitTree->Branch("cov_phi2", &kinFitData.cov_phi2, "cov_phi2/D");
	kinfitTree->Branch("cov_p_theta1", &kinFitData.cov_p_theta1, "cov_p_theta1/D");
	kinfitTree->Branch("cov_p_theta2", &kinFitData.cov_p_theta2, "cov_p_theta2/D");
	kinfitTree->Branch("cov_p_phi1", &kinFitData.cov_p_phi1, "cov_p_phi1/D");
	kinfitTree->Branch("cov_p_phi2", &kinFitData.cov_p_phi2, "cov_p_phi2/D");
	kinfitTree->Branch("cov_phi_theta1", &kinFitData.cov_phi_theta1, "cov_phi_theta1/D");
	kinfitTree->Branch("cov_phi_theta2", &kinFitData.cov_phi_theta2, "cov_phi_theta2/D");
	kinfitTree->Branch("convergence_status", &kinFitData.convergence_status, "convergence_status/D");
	kinfitTree->Branch("cov_err_p1", &kinFitData.cov_err_p1, "cov_err_p1/D");
	kinfitTree->Branch("cov_err_p2", &kinFitData.cov_err_p2, "cov_err_p2/D");
	kinfitTree->Branch("cov_err_theta1", &kinFitData.cov_err_theta1, "cov_err_theta1/D");
	kinfitTree->Branch("cov_err_theta2", &kinFitData.cov_err_theta2, "cov_err_theta2/D");
	kinfitTree->Branch("cov_err_phi1", &kinFitData.cov_err_phi1, "cov_err_phi1/D");
	kinfitTree->Branch("cov_err_phi2", &kinFitData.cov_err_phi2, "cov_err_phi2/D");
	kinfitTree->Branch("cov_err_p_theta1", &kinFitData.cov_err_p_theta1, "cov_err_p_theta1/D");
	kinfitTree->Branch("cov_err_p_theta2", &kinFitData.cov_err_p_theta2, "cov_err_p_theta2/D");
	kinfitTree->Branch("cov_err_p_phi1", &kinFitData.cov_err_p_phi1, "cov_err_p_phi1/D");
	kinfitTree->Branch("cov_err_p_phi2", &kinFitData.cov_err_p_phi2, "cov_err_p_phi2/D");
	kinfitTree->Branch("cov_err_phi_theta1", &kinFitData.cov_err_phi_theta1, "cov_err_phi_theta1/D");
	kinfitTree->Branch("cov_err_phi_theta2", &kinFitData.cov_err_phi_theta2, "cov_err_phi_theta2/D");

	// Create a TTree to store the generated particle information
	TTree *allGenTree = new TTree("allGenTree", "allGenTree");

	allGenStruct allGenData;
	allGenTree->Branch("all_mc_p1", &allGenData.all_mc_p1, "all_mc_p1/D");
	allGenTree->Branch("all_mc_p2", &allGenData.all_mc_p2, "all_mc_p2/D");
	allGenTree->Branch("all_mc_theta1", &allGenData.all_mc_theta1, "all_mc_theta1/D");
	allGenTree->Branch("all_mc_theta2", &allGenData.all_mc_theta2, "all_mc_theta2/D");
	allGenTree->Branch("all_mc_phi1", &allGenData.all_mc_phi1, "all_mc_phi1/D");
	allGenTree->Branch("all_mc_phi2", &allGenData.all_mc_phi2, "all_mc_phi2/D");

	TLorentzVector target(0.0, 0.0, 0.0, 0.938272);
	TLorentzVector beam(0.0, 0.0, 10.6, 10.6);
	TLorentzVector W = beam + target;
	const std::vector<TString> parts = {"#pi^{+}", "#pi^{-}"};
	const std::vector<double> masses = {0.139, 0.139};
	const std::vector<int> required_pids = {211, -211};
	const std::vector<TString> KINES = {"P", "#theta", "#phi"};

	KinFitTest test("InvMass_hipo", parts, W, 0, 0, 0.77);

	//---------------Specify input file or directory path---------------//
	// char dir_path[256] = "/volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/pim-sec2-100mil_events_6-15-23/cooked/";
	// char dir_path[256] = "/volatile/clas12/reedtg/clas12_kinfitter/rho_to_pippim_gen_events_7-20-23/rho_to_pippim/cooked/out_rho_to_pippim130.rec.hipo";
	// char dir_path[256] = "/volatile/clas12/reedtg/clas12_kinfitter/rho_to_pippim_gen_events_7-20-23/rho_to_pippim/cooked/out_rho_to_pippim130.rec.hipo";
	char dir_path[256] = "/volatile/clas12/reedtg/clas12_kinfitter/rho_to_pippim_gen_events_8-23-23/rho_to_pippim/cooked/";
	// char dir_path[256] = "/volatile/clas12/reedtg/clas12_kinfitter/rho_to_pippim_gen_events/rho_to_pippim_8-26-23/cooked/";
	//------------------------------------------------------------------//

	DIR *dr;
	struct dirent *en;
	dr = opendir(dir_path); // open all or present directory

	// If dir_path is a directory, loop through all files within it, with a max number of file limit
	int in_file_count = 0;
	int max_number_file = 10;

	if (dr)
	{
		while ((en = readdir(dr)) != NULL && in_file_count < max_number_file)
		{
			if ((strcmp(en->d_name, ".") != 0) && (strcmp(en->d_name, "..") != 0))
			{
				in_file_count++;
				std::string dir_path_str(dir_path);
				dir_path_str.append(en->d_name);

				int string_l = dir_path_str.length();
				char dir_path_char[256];
				strcpy(dir_path_char, dir_path_str.c_str());
				std::cout << "File name: " << dir_path_char << std::endl;
				read_Hipo(dir_path_char, required_pids, masses, target, beam, test, parts, kinfitTree, kinFitData, allGenTree, allGenData);
			}
		}
		closedir(dr); // close all directory
	}
	else
	{
		std::cout << "Not a directory. Only single file " << dir_path << std::endl;
		read_Hipo(dir_path, required_pids, masses, target, beam, test, parts, kinfitTree, kinFitData, allGenTree, allGenData);
	}
	test.plot();

	// Save the tree to the ROOT file
	kinfitRoot->Write();

	// Close the ROOT file
	kinfitRoot->Close();

	TCanvas *fid_cut = new TCanvas("fid_cut", "fid_cut", 1000, 1000);
	pass_fid_cut->Draw("col");
	fail_fid_cut->SetMarkerColor(kRed);
	fid_cut->SaveAs("fidCut.pdf");

	TCanvas *fid_cut_fail = new TCanvas("fid_cut_fail", "fid_cut", 1000, 1000);
	fail_fid_cut->Draw("col");
	fid_cut_fail->SaveAs("fidCutfail.pdf");

	TCanvas *can_mom = new TCanvas("can_mom", "can_mpm", 1000, 1000);
	P_pip_P_pim->Draw("col");
	can_mom->SaveAs("can_mom.pdf");

	return 0;
}

void read_Hipo(char inputFile[256], std::vector<int> required_pids, std::vector<double> masses, TLorentzVector target, TLorentzVector beam, KinFitTest test, std::vector<TString> parts, TTree *kinfitTree, kinFitInfoStruct &kinFitData, TTree *allGenTree, allGenStruct &allGenData)
{

	TLorentzVector W = beam + target;

	std::vector<int> constraint_idx = {0, 1};

	// std::cout << "Opening input file " << std::endl;
	hipo::reader reader;
	reader.open(inputFile);

	std::cout << "Number of events in file: " << reader.getEntries() << std::endl;
	if (reader.getEntries() < 1)
		return; // Workaround to prevent corrupted files

	hipo::dictionary factory;
	reader.readDictionary(factory);
	hipo::event event;
	hipo::bank PART(factory.getSchema("REC::Particle"));
	hipo::bank TRACK(factory.getSchema("REC::Track"));
	hipo::bank TRAJ(factory.getSchema("REC::Traj"));
	hipo::bank COV(factory.getSchema("REC::CovMat"));
	hipo::bank MCPART(factory.getSchema("MC::Particle"));
	hipo::bank MC_MATCH(factory.getSchema("MC::GenMatch"));
	hipo::bank RUN(factory.getSchema("RUN::config"));

	int event_num;
	int good_events = 0;
	int total_events = 0;
	int fitted_events = 0;
	int failed_fitting = 0;
	int events_passed_kin_cuts = 0;

	// KinCovariance Covariance_PiPlus("../KinematicFitter/pip_minEventCut_covariances.root");  //("pip_covariances.root");//
	// KinCovariance Covariance_PiMinus("../KinematicFitter/pim_minEventCut_covariances.root"); //("pip_covariances.root");//

	// TString input_pip = "/work/clas12/reedtg/clas12_kinematic_fitter/dev_trevor/pip_minEventCut_covariances.root";
	// TString input_pim = "/work/clas12/reedtg/clas12_kinematic_fitter/dev_trevor/pim_minEventCut_covariances.root";
	TString input_pip = "/work/clas12/reedtg/clas12_kinematic_fitter/dev_trevor/pip_covariances_9-18-23.root";
	TString input_pim = "/work/clas12/reedtg/clas12_kinematic_fitter/dev_trevor/pim_covariances_9-19-23.root";

	KinCovariance Covariance_PiPlus(input_pip);
	KinCovariance Covariance_PiMinus(input_pim);

	KinMomentumCorrections MomCorr_PiPlus(input_pip);
	KinMomentumCorrections MomCorr_PiMinus(input_pim);

	double limit_down_P = 1.2;	  // 2.0;    // 5.0;//
	double limit_up_P = 10.8;	  // 3.0;      // 6.0;//
	double limit_down_Theta = 6.; // 15.0;
	double limit_up_Theta = 34.0; // 17.0;
	double limit_down_Phi = 32.0;
	double limit_up_Phi = 88.0;

	while (reader.next() == true /*&& total_events<10*/)
	{

		total_events++;
		reader.read(event);
		event.getStructure(PART);
		event.getStructure(COV);
		event.getStructure(MCPART);
		event.getStructure(MC_MATCH);
		event.getStructure(RUN);
		event.getStructure(TRACK);
		event.getStructure(TRAJ);

		// Get reconstructed particle info
		std::vector<int> pid_list;
		std::vector<int> input_part_index;
		std::vector<TLorentzVector> vec_list;
		std::vector<int> sector_list;
		read_Rec_Part_Bank(PART, MC_MATCH, TRACK, TRAJ, &pid_list, &vec_list, &input_part_index, &sector_list);

		// Get generated particle info
		std::vector<int> mc_pid_list;
		std::vector<int> mc_input_part_index;
		std::vector<TLorentzVector> mc_vec_list;
		read_MC_Part_Bank(MCPART, required_pids, &mc_pid_list, &mc_vec_list, &mc_input_part_index);

		allGenData.all_mc_p1 = mc_vec_list[0].P();
		allGenData.all_mc_p2 = mc_vec_list[1].P();
		allGenData.all_mc_theta1 = mc_vec_list[0].Theta();
		allGenData.all_mc_theta2 = mc_vec_list[1].Theta();
		allGenData.all_mc_phi1 = mc_vec_list[0].Phi();
		allGenData.all_mc_phi2 = mc_vec_list[1].Phi();
		allGenTree->Fill();

		std::vector<std::vector<double>> pull_vals(parts.size(), std::vector<double>(3));
		double conf_lev;

		bool exitEvent = false;
		// Proceed to kinematic fitting only if reconstructed event meets particle requirements
		if (vec_list.size() != parts.size())
		{
			continue;
		}
		else
		{
			good_events++;
			event_num = RUN.getInt("event", 0);
			std::vector<KinParticle> kin_parts;

			P_pip_P_pim->Fill(vec_list[0].P(), vec_list[1].P());

			for (int ipart = 0; ipart < pid_list.size(); ++ipart)
			{

				// Correct the reconstructed momenta to suppress offsets with Gen momenta
				std::cout << "PID " << pid_list[ipart] << " " << ipart << " " << pid_list.size() << std::endl;
				std::cout << "Before " << vec_list[ipart].P() << " " << vec_list[ipart].Theta() * 180. / 3.141592 << " " << vec_list[ipart].Phi() * 180. / 3.141592 << std::endl;

				

				if (!pass_limit_cov_matrix(sector_list[ipart], vec_list[ipart], limit_up_P, limit_down_P, limit_up_Theta, limit_down_Theta, limit_up_Phi, limit_down_Phi))
				{

					InfoFile << "Particle " << ipart << ", REC: Sector: " << sector_list[ipart] << ",P = " << vec_list[ipart].P() << ", theta: " << vec_list[ipart].Theta() * TMath::RadToDeg() << ", phi: " << vec_list[ipart].Phi() * TMath::RadToDeg() << std::endl;
					InfoFile << "Particle " << ipart << ", MC: P = " << mc_vec_list[ipart].P() << ", theta: " << mc_vec_list[ipart].Theta() * TMath::RadToDeg() << ", phi: " << mc_vec_list[ipart].Phi() * TMath::RadToDeg() << std::endl;

					exitEvent = true;
					break;
				}

				KinCovariance Current_Covmatrix = (pid_list[ipart] == 211) ? Covariance_PiPlus : Covariance_PiMinus; // attributes PiMinus Cov matrix to anything but pi+

				// Correct the reconstructed momenta to suppress offsets with Gen momenta
				if (pid_list[ipart] == 211)
				{
					MomCorr_PiPlus.Correct_4_Vector(sector_list[ipart], &vec_list[ipart]);
				}
				if (pid_list[ipart] == -211)
				{
					MomCorr_PiMinus.Correct_4_Vector(sector_list[ipart], &vec_list[ipart]);
				}

				std::cout << "After " << vec_list[ipart].P() << " " << vec_list[ipart].Theta() * 180. / 3.141592 << " " << vec_list[ipart].Phi() * 180. / 3.141592 << std::endl;

				std::cout << "PID " << pid_list[ipart] << std::endl;
				kin_parts.push_back(KinParticle(vec_list[ipart], 0.139, Current_Covmatrix, sector_list[ipart]));
				// std::cout << "Current_Covmatrix: " << kin_parts[ipart].GetCovMatrix()[1][1] << std::endl;
				// kin_parts[ipart].GetCovMatrix().Print();
			}

			// Skip to the next event if at least one of the particles did not make the kinematic cuts
			if (exitEvent)
			{
				continue;
			}

			std::cout << "here" << std::endl;
			auto kin = new KinFitter({KinParticle(target), KinParticle(beam)}, kin_parts);
			kin->Add_InvMass_Constraint(constraint_idx, 0.77);
			kin->DoFitting(100);

			events_passed_kin_cuts++;
			bool is_background = false;
			double weight = 1.0;

			if (kin->HasConverged())
			{
				std::vector<TLorentzVector> parts_fit = kin->GetFitted4Vectors();
				double conf_lev = kin->GetConfidenceLevel();
				for (int ipart = 0; ipart < parts.size(); ipart++)
				{
					for (int jkine = 0; jkine < KINES.size(); jkine++)
					{
						pull_vals[ipart][jkine] = kin->GetPulls()[ipart * KINES.size() + jkine];
					}
				}
				// std::cout << "parts_fit[0].P(): " << parts_fit[0].P() << std::endl;
				// std::tie(pull_vals, conf_lev) = test.fill_InvariantMass(kin, mc_vec_list, vec_list, constraint_idx, weight, is_background);
				test.fill_InvariantMass(kin, mc_vec_list, vec_list, constraint_idx, weight, is_background);
				fitted_events++;

				kinFitData.fit_p1 = parts_fit[0].P();
				kinFitData.fit_p2 = parts_fit[1].P();
				kinFitData.fit_theta1 = parts_fit[0].Theta();
				kinFitData.fit_theta2 = parts_fit[1].Theta();
				kinFitData.fit_phi1 = parts_fit[0].Phi();
				kinFitData.fit_phi2 = parts_fit[1].Phi();

				kinFitData.pull_p1 = pull_vals[0][0];
				kinFitData.pull_p2 = pull_vals[1][0];
				kinFitData.pull_theta1 = pull_vals[0][1];
				kinFitData.pull_theta2 = pull_vals[1][1];
				kinFitData.pull_phi1 = pull_vals[0][2];
				kinFitData.pull_phi2 = pull_vals[1][2];

				kinFitData.confidence_level = conf_lev;
			}
			else
			{
				failed_fitting++;

				kinFitData.fit_p1 = 0;
				kinFitData.fit_p2 = 0;
				kinFitData.fit_theta1 = 0;
				kinFitData.fit_theta2 = 0;
				kinFitData.fit_phi1 = 0;
				kinFitData.fit_phi2 = 0;

				kinFitData.pull_p1 = 0;
				kinFitData.pull_p2 = 0;
				kinFitData.pull_theta1 = 0;
				kinFitData.pull_theta2 = 0;
				kinFitData.pull_phi1 = 0;
				kinFitData.pull_phi2 = 0;

				kinFitData.confidence_level = -1;
			}

			kinFitData.rec_p1 = vec_list[0].P();
			kinFitData.rec_p2 = vec_list[1].P();
			kinFitData.rec_theta1 = vec_list[0].Theta();
			kinFitData.rec_theta2 = vec_list[1].Theta();
			kinFitData.rec_phi1 = vec_list[0].Phi();
			kinFitData.rec_phi2 = vec_list[1].Phi();

			kinFitData.mc_p1 = mc_vec_list[0].P();
			kinFitData.mc_p2 = mc_vec_list[1].P();
			kinFitData.mc_theta1 = mc_vec_list[0].Theta();
			kinFitData.mc_theta2 = mc_vec_list[1].Theta();
			kinFitData.mc_phi1 = mc_vec_list[0].Phi();
			kinFitData.mc_phi2 = mc_vec_list[1].Phi();

			kinFitData.cov_p1 = kin_parts[0].GetCovMatrix()[0][0];
			kinFitData.cov_p2 = kin_parts[1].GetCovMatrix()[0][0];
			kinFitData.cov_theta1 = kin_parts[0].GetCovMatrix()[1][1];
			kinFitData.cov_theta2 = kin_parts[1].GetCovMatrix()[1][1];
			kinFitData.cov_phi1 = kin_parts[0].GetCovMatrix()[2][2];
			kinFitData.cov_phi2 = kin_parts[1].GetCovMatrix()[2][2];
			kinFitData.cov_p_theta1 = kin_parts[0].GetCovMatrix()[0][1];
			kinFitData.cov_p_theta2 = kin_parts[1].GetCovMatrix()[0][1];
			kinFitData.cov_p_phi1 = kin_parts[0].GetCovMatrix()[0][2];
			kinFitData.cov_p_phi2 = kin_parts[1].GetCovMatrix()[0][2];
			kinFitData.cov_phi_theta1 = kin_parts[0].GetCovMatrix()[1][2];
			kinFitData.cov_phi_theta2 = kin_parts[1].GetCovMatrix()[1][2];

			kinFitData.convergence_status = kin->GetConvergenceStatus();

			kinFitData.cov_err_p1 = kin_parts[0].GetErrCovMatrix()[0][0];
			kinFitData.cov_err_p2 = kin_parts[1].GetErrCovMatrix()[0][0];
			kinFitData.cov_err_theta1 = kin_parts[0].GetErrCovMatrix()[1][1];
			kinFitData.cov_err_theta2 = kin_parts[1].GetErrCovMatrix()[1][1];
			kinFitData.cov_err_phi1 = kin_parts[0].GetErrCovMatrix()[2][2];
			kinFitData.cov_err_phi2 = kin_parts[1].GetErrCovMatrix()[2][2];
			kinFitData.cov_err_p_theta1 = kin_parts[0].GetErrCovMatrix()[0][1];
			kinFitData.cov_err_p_theta2 = kin_parts[1].GetErrCovMatrix()[0][1];
			kinFitData.cov_err_p_phi1 = kin_parts[0].GetErrCovMatrix()[0][2];
			kinFitData.cov_err_p_phi2 = kin_parts[1].GetErrCovMatrix()[0][2];
			kinFitData.cov_err_phi_theta1 = kin_parts[0].GetErrCovMatrix()[1][2];
			kinFitData.cov_err_phi_theta2 = kin_parts[1].GetErrCovMatrix()[1][2];

			kinfitTree->Fill();

			delete kin;
		}
	}

	std::cout << "Number of events that meet particle requirements = " << good_events << ", Number of fitted events = " << fitted_events << std::endl;
	std::cout << "Number of events that passed kinematic cuts = " << events_passed_kin_cuts << ", Number failed fits = " << failed_fitting << std::endl;
}

void read_Rec_Part_Bank(hipo::bank PartBank, hipo::bank MCMatchBank, hipo::bank TrackBank, hipo::bank TrajBank, std::vector<int> *pid_list, std::vector<TLorentzVector> *vec_list, std::vector<int> *in_part_index, std::vector<int> *in_part_sector)
{

	TLorentzVector piplus_vec;
	TLorentzVector piminus_vec;

	for (int index_particle=0; i<PartBank.getRows(); i++)
	{
		int status_pi = PartBank.getInt("status", index_particle);
		int pid_pi = PartBank.getInt("pid", index_particle);
		float px_pi = PartBank.getFloat("px", index_particle);
		float py_pi = PartBank.getFloat("py", index_particle);
		float pz_pi = PartBank.getFloat("pz", index_particle);

		if (((int)(abs(status_pi) / 1000) == 2) && (int)pid_pi == 211 && pass_fid_cut_DC(index_particle, TrajBank) && PartBank.getInt("charge", index_particle) == 1)
		{
			piplus_vec.SetXYZM(px_pi, py_pi, pz_pi, 0.139);
			vec_list->push_back(piplus_vec);
			pid_list->push_back(pid_pi);
			in_part_index->push_back(index_particle);
			in_part_sector->push_back(get_sector(index_particle, TrackBank));
		}

		if (((int)(abs(status_pi) / 1000) == 2) && (int)pid_pi == -211 && pass_fid_cut_DC(index_particle, TrajBank) && PartBank.getInt("charge", index_particle) == -1)
		{
			piminus_vec.SetXYZM(px_pi, py_pi, pz_pi, 0.139);
			vec_list->push_back(piminus_vec);
			pid_list->push_back(pid_pi);
			in_part_index->push_back(index_particle);
			in_part_sector->push_back(get_sector(index_particle, TrackBank));
		}
	}
}

void read_MC_Part_Bank(hipo::bank PartBank, std::vector<int> required_pids, std::vector<int> *pid_list, std::vector<TLorentzVector> *mc_vec_list, std::vector<int> *in_part_index)
{

	TLorentzVector piplus_vec;
	pid_list->push_back(PartBank.getInt("pid", 2));
	in_part_index->push_back(0);
	piplus_vec.SetXYZM(PartBank.getFloat("px", 2), PartBank.getFloat("py", 2), PartBank.getFloat("pz", 2), 0.139);
	mc_vec_list->push_back(piplus_vec);

	TLorentzVector piminus_vec;
	pid_list->push_back(PartBank.getInt("pid", 3));
	in_part_index->push_back(1);
	piminus_vec.SetXYZM(PartBank.getFloat("px", 3), PartBank.getFloat("py", 3), PartBank.getFloat("pz", 3), 0.139);
	mc_vec_list->push_back(piminus_vec);
}

int get_sector(int index, hipo::bank TrackBank)
{
	std::map<int, std::vector<int>> trackBankMap = loadMapByIndex(TrackBank, "pindex");
	auto it = trackBankMap.find(index);
	if (it != trackBankMap.end())
	{
		const std::vector<int> &indices = it->second;
		if (!indices.empty())
		{
			// Assuming sector is the same for all indices in the vector
			int sector = TrackBank.getInt("sector", indices[0]);
			return sector;
		}
	}
	return -1; // Index not found in the map
}

bool pass_fid_cut_DC(int index, hipo::bank TrajBank)
{
	int nrows_traj = TrajBank.getRows();
	for (int i = 0; i < nrows_traj; i++)
	{
		int pindex = TrajBank.getInt("pindex", i);
		int detector = TrajBank.getInt("detector", i);
		int layer = TrajBank.getInt("layer", i);
		float edge = TrajBank.getFloat("edge", i);

		if (pindex == index && detector == 6 && layer == 6 && edge < 5)
			return false;
	}
	return true;
}

bool pass_limit_cov_matrix(int sector, TLorentzVector vector, double limit_up_P, double limit_down_P, double limit_up_Theta, double limit_down_Theta, double limit_up_Phi, double limit_down_Phi)
{
	if (vector.P() < limit_down_P || vector.P() > limit_up_P || vector.Theta() * TMath::RadToDeg() > limit_up_Theta || vector.Theta() * TMath::RadToDeg() < limit_down_Theta)
		return false;

	double Temp_Phi_in_vector = (vector.Phi() * TMath::RadToDeg() < 0. && sector > 1) ? vector.Phi() * TMath::RadToDeg() + 360 : vector.Phi() * TMath::RadToDeg(); // implement that and edge fiducial cuts
	double Phi_in_vector = (sector == 1) ? Temp_Phi_in_vector + 60 : Temp_Phi_in_vector - (sector - 2) * 60.;

	if (Phi_in_vector < limit_down_Phi || Phi_in_vector > limit_up_Phi)
		return false;

	return true;
}

int main()
{
	return test_InvMass_hipo_truth_matching();
}
