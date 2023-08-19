#include <vector>
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "KinFitter.h"
#include "KinParticle.h"
#include "KinCovariance.h"
#include "KinFitTest.h"
#include "reader.h"
#include <dirent.h>
#include "TTree.h"
#include "TH2F.h"
#include <TROOT.h>
#include <sys/types.h>
#include <cstring>

void read_Hipo(char input_hipo[256], std::vector<int> required_pids, std::vector<double> masses, TLorentzVector target, TLorentzVector beam, KinFitTest test, std::vector<TString> parts, TTree* kinfitTree, double& P_part_1, double& P_part_2, double& theta_part_1, double& theta_part_2, double& phi_part_1, double& phi_part_2, double& P_pulls_part_1, double& P_pulls_part_2, double& theta_pulls_part_1, double& theta_pulls_part_2, double& phi_pulls_part_1, double& phi_pulls_part_2, double& confidence_level);
void read_Rec_Part_Bank(hipo::bank PartBank, hipo::bank MCMatchBank, hipo::bank TrackBank, hipo::bank TrajBank, std::vector<int> *pid_list, std::vector<TLorentzVector> *vec_list, std::vector<int> *in_part_index, std::vector<int> *in_part_sector);
void read_MC_Part_Bank(hipo::bank Bank, std::vector<int> required_pids, std::vector<int> *pid_list, std::vector<TLorentzVector> *vec_list, std::vector<int> *in_part_index);
int get_sector(int index, hipo::bank TrackBank);
bool pass_fid_cut_DC(int index, hipo::bank TrajBank);
bool pass_limit_cov_matrix(int sector, TLorentzVector vector, double limit_up_P, double limit_down_P, double limit_up_Theta, double limit_down_Theta, double limit_up_Phi, double limit_down_Phi);

TH2F *pass_fid_cut = new TH2F("pass", "pass", 60, -180, 180, 60, 0, 45);
TH2F *fail_fid_cut = new TH2F("pass", "pass", 60, -180, 180, 60, 0, 45);
TH2F *P_pip_P_pim = new TH2F("PpipVSPpim", "PpipVSPpim", 60, 1.0, 11., 60, 1.0, 11.0);

int test_InvMass_hipo_truth_matching()
{

    gROOT->SetBatch(kTRUE);

    // Create a ROOT file to store the kinematic fit info    
    TFile* kinfitRoot = new TFile("kinfit.root", "RECREATE");

    // Create a TTree to store kinematic fitting information  
    TTree* kinfitTree = new TTree("kinFitInfoTree", "kinFitInfoTree");
    //float P_part_1;
    double P_part_1, P_part_2, theta_part_1, theta_part_2, phi_part_1, phi_part_2, P_pulls_part_1, P_pulls_part_2, theta_pulls_part_1, theta_pulls_part_2, phi_pulls_part_1, phi_pulls_part_2, confidence_level;
    kinfitTree->Branch("P_part_1", &P_part_1, "P_part_1/D");
    kinfitTree->Branch("P_part_2", &P_part_2, "P_part_2/D");
    kinfitTree->Branch("theta_part_1", &theta_part_1, "theta_part_1/D");
    kinfitTree->Branch("theta_part_2", &theta_part_2, "theta_part_2/D");
    kinfitTree->Branch("phi_part_1", &phi_part_1, "phi_part_1/D");
    kinfitTree->Branch("phi_part_2", &phi_part_2, "phi_part_2/D");
    kinfitTree->Branch("P_pulls_part_1", &P_pulls_part_1, "P_pulls_part_1/D");
    kinfitTree->Branch("P_pulls_part_2", &P_pulls_part_2, "P_pulls_part_2/D");
    kinfitTree->Branch("theta_pulls_part_1", &theta_pulls_part_1, "theta_pulls_part_1/D");
    kinfitTree->Branch("theta_pulls_part_2", &theta_pulls_part_2, "theta_pulls_part_2/D");
    kinfitTree->Branch("phi_pulls_part_1", &phi_pulls_part_1, "phi_pulls_part_1/D");
    kinfitTree->Branch("phi_pulls_part_2", &phi_pulls_part_2, "phi_pulls_part_2/D");
    kinfitTree->Branch("confidence_level", &confidence_level, "confidence_level/D");

    TLorentzVector target(0.0, 0.0, 0.0, 0.938272);
    TLorentzVector beam(0.0, 0.0, 10.6, 10.6);
    TLorentzVector W = beam + target;
    const std::vector<TString> parts = {"#pi^{+}", "#pi^{-}"};
    const std::vector<double> masses = {0.139, 0.139};
    const std::vector<int> required_pids = {211, -211};

    KinFitTest test("InvMass_hipo", parts, W, 0, 0, 0.77);

    //---------------Specify input file or directory path---------------//
    //char dir_path[256] = "../KinematicFitter/Input_files_for_test//New_simu_for_kinfit.hipo"; // New_version_dipion_kinfit.hipo";
    char dir_path[256] = "/volatile/clas12/reedtg/clas12_kinfitter/rho_to_pippim_gen_events_7-20-23/rho_to_pippim/cooked/";
    //char dir_path[256] = "/volatile/clas12/reedtg/clas12_kinfitter/rho_to_pippim_gen_events_7-20-23/rho_to_pippim/cooked/out_rho_to_pippim130.rec.hipo";
    //------------------------------------------------------------------//

    DIR *dr;
    struct dirent *en;
    dr = opendir(dir_path); // open all or present directory

    // If dir_path is a directory, loop through all files within it, with a max number of file limit
    int in_file_count = 0;
    int max_number_file = 3;

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
                read_Hipo(dir_path_char, required_pids, masses, target, beam, test, parts, kinfitTree, P_part_1, P_part_2, theta_part_1, theta_part_2, phi_part_1, phi_part_2, P_pulls_part_1, P_pulls_part_2, theta_pulls_part_1, theta_pulls_part_2, phi_pulls_part_1, phi_pulls_part_2, confidence_level);
            }
        }
        closedir(dr); // close all directory
    }
    else
    {
        std::cout << "Not a directory. Only single file " << dir_path << std::endl;
        read_Hipo(dir_path, required_pids, masses, target, beam, test, parts, kinfitTree, P_part_1, P_part_2, theta_part_1, theta_part_2, phi_part_1, phi_part_2, P_pulls_part_1, P_pulls_part_2, theta_pulls_part_1, theta_pulls_part_2, phi_pulls_part_1, phi_pulls_part_2, confidence_level);
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

// std::ofstream InfoFile("event_info.txt");

void read_Hipo(char inputFile[256], std::vector<int> required_pids, std::vector<double> masses, TLorentzVector target, TLorentzVector beam, KinFitTest test, std::vector<TString> parts, TTree* kinfitTree, double& P_part_1, double& P_part_2, double& theta_part_1, double& theta_part_2, double& phi_part_1, double& phi_part_2, double& P_pulls_part_1, double& P_pulls_part_2, double& theta_pulls_part_1, double& theta_pulls_part_2, double& phi_pulls_part_1, double& phi_pulls_part_2, double& confidence_level)
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

    //KinCovariance Covariance_PiPlus("../KinematicFitter/pip_minEventCut_covariances.root");  //("pip_covariances.root");//
    //KinCovariance Covariance_PiMinus("../KinematicFitter/pim_minEventCut_covariances.root"); //("pip_covariances.root");//

    KinCovariance Covariance_PiPlus("/work/clas12/reedtg/clas12_kinematic_fitter/dev_trevor/pip_minEventCut_covariances.root");  //("pip_covariances.root");//
    KinCovariance Covariance_PiMinus("/work/clas12/reedtg/clas12_kinematic_fitter/dev_trevor/pim_minEventCut_covariances.root"); //("pip_covariances.root");//

    double limit_down_P = 1.2;//2.0;    // 5.0;//
    double limit_up_P = 10.8;//3.0;      // 6.0;//
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

	std::vector<double> pull_vals;
	double conf_lev;
        // Proceed to kinematic fitting if reconstructed event meets particle requirements
        if (vec_list.size() == 2)
        {
            good_events++;
            event_num = RUN.getInt("event", 0);
            std::vector<KinParticle> kin_parts;

            P_pip_P_pim->Fill(vec_list[0].P(), vec_list[1].P());

            for (int ipart = 0; ipart < pid_list.size(); ++ipart)
            {
                if (!pass_limit_cov_matrix(sector_list[ipart], vec_list[ipart], limit_up_P, limit_down_P, limit_up_Theta, limit_down_Theta, limit_up_Phi, limit_down_Phi))
                    continue;

                KinCovariance Current_Covmatrix = (pid_list[ipart] == 211) ? Covariance_PiPlus : Covariance_PiMinus; // attributes PiMinus Cov matrix to anything but pi+
                kin_parts.push_back(KinParticle(vec_list[ipart], 0.139, Current_Covmatrix, sector_list[ipart]));
            }

            if (!pass_limit_cov_matrix(sector_list[0], vec_list[0], limit_up_P, limit_down_P, limit_up_Theta, limit_down_Theta, limit_up_Phi, limit_down_Phi))
                    continue;
            if (!pass_limit_cov_matrix(sector_list[1], vec_list[1], limit_up_P, limit_down_P, limit_up_Theta, limit_down_Theta, limit_up_Phi, limit_down_Phi))
                    continue;

            auto kin = new KinFitter({KinParticle(target), KinParticle(beam)}, kin_parts);
            kin->Add_InvMass_Constraint(constraint_idx, 0.77);
            kin->DoFitting(100);

            bool is_background = false;
            double weight = 1.0;
            if (kin->HasConverged())
            {
	        std::tie(pull_vals, conf_lev) = test.fill_InvariantMass(kin, mc_vec_list, vec_list, constraint_idx, weight, is_background);
		fitted_events++;
		//std::cout << conf_lev << std::endl;
		//double P_part_1, P_part_2, theta_part_1, theta_part_2, phi_part_1, phi_part_2, P_pulls_part_1, P_pulls_part_2, theta_pulls_part_1, theta_pulls_part_2, phi_pulls_part_1, phi_pulls_part_2, confidence_level;
		
		
		//Fill the tree with this bin's kinfit info only if the confidence requirement is met
		//Otherwise, the pulls are not saved from KinFitTest.h and trying to retrieve them here will result in segmentation fault
		if (conf_lev > 0.01) {
		  P_part_1 = vec_list[0].P();
		  P_part_2 = vec_list[1].P();
		  theta_part_1 = vec_list[0].Theta();
		  theta_part_2 = vec_list[1].Theta();
		  phi_part_1 = vec_list[0].Phi();
		  phi_part_2 = vec_list[1].Phi();
		  P_pulls_part_1 = pull_vals[0];
		  P_pulls_part_2 = pull_vals[1];
		  theta_pulls_part_1 = pull_vals[3];
		  theta_pulls_part_2 = pull_vals[4];
		  phi_pulls_part_1 = pull_vals[5];
		  phi_pulls_part_2 = pull_vals[6];
		  confidence_level = conf_lev;

		  kinfitTree->Fill();
		}
            }
            delete kin;
        }
    }

    // Save the tree to the ROOT file
    //tree->Write();

    // Close the ROOT file
    //root_file->Close();
    
    std::cout << "Number of fitted events " << fitted_events << std::endl;
}

void read_Rec_Part_Bank(hipo::bank PartBank, hipo::bank MCMatchBank, hipo::bank TrackBank, hipo::bank TrajBank, std::vector<int> *pid_list, std::vector<TLorentzVector> *vec_list, std::vector<int> *in_part_index, std::vector<int> *in_part_sector)
{
    // fail safe for truth matching
    if (MCMatchBank.getRows() < 1)
        return;

    TLorentzVector piplus_vec;
    TLorentzVector piminus_vec;

    int index_pi_plus = MCMatchBank.getInt("pindex", 2);
    int index_pi_minus = MCMatchBank.getInt("pindex", 3);

    // fail safe for truth matching
    if (index_pi_plus == index_pi_minus)
        return;

    if (index_pi_plus < 0 || index_pi_minus < 0)
        return;

    int status_pi_plus = PartBank.getInt("status", index_pi_plus);
    if (((int)(abs(status_pi_plus) / 1000) == 2) && pass_fid_cut_DC(index_pi_plus, TrajBank))
    {

        int pid_pi_plus = PartBank.getInt("pid", index_pi_plus);
        float px_pi_plus = PartBank.getFloat("px", index_pi_plus);
        float py_pi_plus = PartBank.getFloat("py", index_pi_plus);
        float pz_pi_plus = PartBank.getFloat("pz", index_pi_plus);

        piplus_vec.SetXYZM(px_pi_plus, py_pi_plus, pz_pi_plus, 0.139);
        vec_list->push_back(piplus_vec);
        pid_list->push_back(pid_pi_plus);
        in_part_index->push_back(index_pi_plus);
        in_part_sector->push_back(get_sector(index_pi_plus, TrackBank));
    }

    int status_pi_minus = PartBank.getInt("status", index_pi_minus);
    if (((int)(abs(status_pi_minus) / 1000) == 2) && pass_fid_cut_DC(index_pi_minus, TrajBank))
    {
        int pid_pi_minus = PartBank.getInt("pid", index_pi_minus);
        float px_pi_minus = PartBank.getFloat("px", index_pi_minus);
        float py_pi_minus = PartBank.getFloat("py", index_pi_minus);
        float pz_pi_minus = PartBank.getFloat("pz", index_pi_minus);

        piminus_vec.SetXYZM(px_pi_minus, py_pi_minus, pz_pi_minus, 0.139);
        vec_list->push_back(piminus_vec);
        pid_list->push_back(pid_pi_minus);
        in_part_index->push_back(index_pi_minus);
        in_part_sector->push_back(get_sector(index_pi_minus, TrackBank));

        if (piminus_vec.P() > 0.0)
            pass_fid_cut->Fill(piminus_vec.Phi() * TMath::RadToDeg(), piminus_vec.Theta() * TMath::RadToDeg());
    }

    if (((int)(abs(status_pi_minus) / 1000) == 2) && !pass_fid_cut_DC(index_pi_minus, TrajBank))
    {
        float px_pi_minus = PartBank.getFloat("px", index_pi_minus);
        float py_pi_minus = PartBank.getFloat("py", index_pi_minus);
        float pz_pi_minus = PartBank.getFloat("pz", index_pi_minus);

        piminus_vec.SetXYZM(px_pi_minus, py_pi_minus, pz_pi_minus, 0.139);
        fail_fid_cut->Fill(piminus_vec.Phi() * TMath::RadToDeg(), piminus_vec.Theta() * TMath::RadToDeg());
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
    int sector = -1;
    int nrows_track = TrackBank.getRows();
    for (int i = 0; i < nrows_track; i++)
    {
        int pindex = TrackBank.getInt("pindex", i);
        int in_sector = TrackBank.getInt("sector", i);

        if (pindex == index)
            sector = in_sector;
    }

    return sector;
}

bool pass_fid_cut_DC(int index, hipo::bank TrajBank)
{
    int nrows_traj = TrajBank.getRows();
    
    bool pass_cut = true;
    for (int i = 0; i < nrows_traj; i++)
    {
        int pindex = TrajBank.getInt("pindex", i);
        int detector = TrajBank.getInt("detector", i);
        int layer = TrajBank.getInt("layer", i);
        float edge = TrajBank.getFloat("edge", i);

        pass_cut = (pass_cut && !(pindex == index && detector == 6 && layer == 6 && edge < 5.));
    }
    return pass_cut;
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

