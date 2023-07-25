#include <vector>
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "KinFitter.h"
#include "KinParticle.h"
#include "KinCovariance.h"
#include "KinFitTest.h"
#include "reader.h"
#include <dirent.h>
#include <sys/types.h>
#include <cstring>

void read_Hipo(char input_hipo[256], std::vector<int> required_pids, std::vector<double> masses, TLorentzVector target, TLorentzVector beam, KinFitTest test, std::vector<TString> parts);
void read_Rec_Part_Bank(hipo::bank Bank, hipo::bank Match_Bank, hipo::bank Track_Bank, std::vector<int> required_pids, std::vector<int> *pid_list, std::vector<float> *px_list, std::vector<float> *py_list, std::vector<float> *pz_list, std::vector<int> *in_part_index, std::vector<int> *in_part_sector);
void read_MC_Part_Bank(hipo::bank Bank, std::vector<int> required_pids, std::vector<int> *pid_list, std::vector<float> *px_list, std::vector<float> *py_list, std::vector<float> *pz_list, std::vector<int> *in_part_index);

int test_InvMass_hipo_truth_matching()
{

    gROOT->SetBatch(kTRUE);

    TLorentzVector target(0.0, 0.0, 0.0, 0.938272);
    TLorentzVector beam(0.0, 0.0, 10.6, 10.6);
    TLorentzVector W = beam + target;
    const std::vector<TString> parts = {"#pi^{+}", "#pi^{-}"};
    const std::vector<double> masses = {0.139, 0.139};
    const std::vector<int> required_pids = {211, -211};

    KinFitTest test("InvMass_hipo", parts, W, 0, 0, 0.77);

    //---------------Specify input file or directory path---------------//
    char dir_path[256] = "../KinematicFitter/Input_files_for_test/";//New_version_dipion_kinfit.hipo";
    //------------------------------------------------------------------//

    DIR *dr;
    struct dirent *en;
    dr = opendir(dir_path); // open all or present directory

    // If dir_path is a directory, loop through all files within it
    int in_file_count = 0;
    if (dr)
    {
        while ((en = readdir(dr)) != NULL)
        {
            if ((strcmp(en->d_name, ".") != 0) && (strcmp(en->d_name, "..") != 0))
            {
                in_file_count++;
                std::string dir_path_str(dir_path);
                dir_path_str.append(en->d_name);
                // dir_path_str.append("/dst.hipo");
                //  Convert file path string back to char to be used by read_Hipo
                int string_l = dir_path_str.length();
                char dir_path_char[256];
                strcpy(dir_path_char, dir_path_str.c_str());
                read_Hipo(dir_path_char, required_pids, masses, target, beam, test, parts);
                // Print list of all read hipo files to input_files.txt
                // InFileList << dir_path << en->d_name << "/dst.hipo" << std::endl;
            }
        }
        closedir(dr); // close all directory
    }
    else
    {
        std::cout << "Not a directory. Only single file " << dir_path << std::endl;
        read_Hipo(dir_path, required_pids, masses, target, beam, test, parts);
    }
    test.plot();

    return 0;
}

// std::ofstream InfoFile("event_info.txt");

void read_Hipo(char inputFile[256], std::vector<int> required_pids, std::vector<double> masses, TLorentzVector target, TLorentzVector beam, KinFitTest test, std::vector<TString> parts)
{
    TLorentzVector W = beam + target;

    std::vector<int> constraint_idx = {0, 1};

    // std::cout << "Opening input file " << std::endl;
    hipo::reader reader;
    reader.open(inputFile);
    hipo::dictionary factory;
    reader.readDictionary(factory);
    hipo::event event;
    hipo::bank PART(factory.getSchema("REC::Particle"));
    hipo::bank TRACK(factory.getSchema("REC::Track"));
    hipo::bank COV(factory.getSchema("REC::CovMat"));
    hipo::bank MCPART(factory.getSchema("MC::Particle"));
    hipo::bank MC_MATCH(factory.getSchema("MC::GenMatch"));
    hipo::bank RUN(factory.getSchema("RUN::config"));

    int event_num;
    int good_events = 0;
    int total_events = 0;
    int fitted_events = 0;

    KinCovariance Covariance_PiPlus("../KinematicFitter/pip_covariances.root");  //("pip_covariances.root");//
    KinCovariance Covariance_PiMinus("../KinematicFitter/pim_covariances.root"); //("pip_covariances.root");//

    // while (nevents < max_events)
    while (reader.next() == true)
    {

        cout << endl;
        cout << "New event" << endl;

        total_events++;
        reader.read(event);
        event.getStructure(PART);
        event.getStructure(COV);
        event.getStructure(MCPART);
        event.getStructure(MC_MATCH);
        event.getStructure(RUN);
        event.getStructure(TRACK);

        // Get reconstructed particle info
        std::vector<int> pid_list;
        std::vector<int> input_part_index;
        std::vector<float> px_list;
        std::vector<float> py_list;
        std::vector<float> pz_list;
        std::vector<int> sector_list;
        read_Rec_Part_Bank(PART, MC_MATCH, TRACK, required_pids, &pid_list, &px_list, &py_list, &pz_list, &input_part_index, &sector_list);
        std::vector<TLorentzVector> parts_rec;

        // Get generated particle info
        std::vector<int> mc_pid_list;
        std::vector<int> mc_input_part_index;
        std::vector<float> mc_px_list;
        std::vector<float> mc_py_list;
        std::vector<float> mc_pz_list;
        read_MC_Part_Bank(MCPART, required_pids, &mc_pid_list, &mc_px_list, &mc_py_list, &mc_pz_list, &mc_input_part_index);
        std::vector<TLorentzVector> parts_mc;

        // Check if event (reconstructed) contains exactly one of each of required particles
        bool contains_required_parts = false;
        for (int k = 0; k < required_pids.size(); k++)
        {
            int pid_count = count(pid_list.begin(), pid_list.end(), required_pids.at(k));
            // InfoFile << "pid count " << pid_count << std::endl;
            if (pid_count == 1)
            {
                contains_required_parts = true;
                // InfoFile << "PID from 'required pid' list: " << required_pids.at(k) << std::endl;
            }
            else
            {
                contains_required_parts = false;
                // InfoFile << "Break from event for-loop" << std::endl;
                break;
            }
        }

        // Proceed to kinematic fitting if reconstructed event meets particle requirements
        if (contains_required_parts)
        {
            good_events++;
            // std::cout << "Good event" << std::endl;
            event_num = RUN.getInt("event", 0);
            // InfoFile << "Event Number: " << event_num << std::endl;
            std::vector<KinParticle> kin_parts;
            TLorentzVector current_part_vec;
            TLorentzVector current_part_mc_vec;

            for (int ipart = 0; ipart < pid_list.size(); ++ipart)
            {
                double px_part = px_list[ipart], py_part = py_list[ipart], pz_part = pz_list[ipart];
                double current_part_P_mag = sqrt(px_part * px_part + py_part * py_part + pz_part * pz_part);
                double current_part_E = sqrt(current_part_P_mag * current_part_P_mag + 0.139 * 0.139);

                std::cout << "sector_list " << sector_list.size() << endl;
                std::cout << "sector_list " << sector_list[ipart] << endl;
                current_part_vec.SetPxPyPzE(px_part, py_part, pz_part, current_part_E);
                parts_rec.push_back(current_part_vec);

                KinCovariance Current_Covmatrix = (pid_list[ipart] = 211) ? Covariance_PiPlus : Covariance_PiMinus; // attributes PiMinus Cov matrix to anything but pi+
                std::cout << "here" << endl;
                kin_parts.push_back(KinParticle(current_part_vec, 0.139, Current_Covmatrix, sector_list[ipart]));
                std::cout << "there" << endl;
                // InfoFile << "Reconstructed particle (PID " << pid_list[ipart] << ") 3-momentum: " << px_part << " " <<  py_part << " " <<  pz_part << std::endl;
            }

            for (int ipart = 0; ipart < mc_pid_list.size(); ++ipart)
            {
                double mc_px_part = mc_px_list[ipart], mc_py_part = mc_py_list[ipart], mc_pz_part = mc_pz_list[ipart];
                double current_mc_part_P_mag = sqrt(mc_px_part * mc_px_part + mc_py_part * mc_py_part + mc_pz_part * mc_pz_part);
                double current_mc_part_E = sqrt(current_mc_part_P_mag * current_mc_part_P_mag + 0.139 * 0.139);
                current_part_mc_vec.SetPxPyPzE(mc_px_part, mc_py_part, mc_pz_part, current_mc_part_E);
                parts_mc.push_back(current_part_mc_vec);
            }

            auto kin = new KinFitter({KinParticle(target), KinParticle(beam)}, kin_parts);
            kin->Add_InvMass_Constraint(constraint_idx, 0.77);
            kin->DoFitting(100);
            bool is_background = false;
            double weight = 1.0;
            if(kin->HasConverged()){
                test.fill_InvariantMass(kin, parts_mc, parts_rec, constraint_idx, weight, is_background);
                fitted_events++;
            }
            delete kin;
        }
    }

    cout << "Number of fitted events " <<fitted_events<< endl;
}

void read_Rec_Part_Bank(hipo::bank PartBank, hipo::bank MCMatchBank, hipo::bank TrackBank, std::vector<int> required_pids, std::vector<int> *pid_list, std::vector<float> *px_list, std::vector<float> *py_list, std::vector<float> *pz_list, std::vector<int> *in_part_index, std::vector<int> *in_part_sector)
{

    if (MCMatchBank.getRows() < 1)
        return;
    // MCMatchBank.show();

    int index_pi_plus = MCMatchBank.getInt("pindex", 2);
    int index_pi_minus = MCMatchBank.getInt("pindex", 3);

    if (index_pi_plus < 0 || index_pi_minus < 0)
        return;

    int status_pi_plus = PartBank.getInt("status", index_pi_plus);
    int pid_pi_plus = PartBank.getInt("pid", index_pi_plus);
    float px_pi_plus = PartBank.getFloat("px", index_pi_plus);
    float py_pi_plus = PartBank.getFloat("py", index_pi_plus);
    float pz_pi_plus = PartBank.getFloat("pz", index_pi_plus);

    int status_pi_minus = PartBank.getInt("status", index_pi_minus);
    int pid_pi_minus = PartBank.getInt("pid", index_pi_minus);
    float px_pi_minus = PartBank.getFloat("px", index_pi_minus);
    float py_pi_minus = PartBank.getFloat("py", index_pi_minus);
    float pz_pi_minus = PartBank.getFloat("pz", index_pi_minus);

    if (abs(status_pi_plus) < 4000)
    {
        pid_list->push_back(pid_pi_plus);
        in_part_index->push_back(0);
        px_list->push_back(px_pi_plus);
        py_list->push_back(py_pi_plus);
        pz_list->push_back(pz_pi_plus);
    }

    if (abs(status_pi_minus) < 4000)
    {
        pid_list->push_back(pid_pi_minus);
        in_part_index->push_back(1);
        px_list->push_back(px_pi_minus);
        py_list->push_back(py_pi_minus);
        pz_list->push_back(pz_pi_minus);
    }

    int nrows_track = TrackBank.getRows();
    for (int i = 0; i < nrows_track; i++)
    {
        int pindex = TrackBank.getInt("pindex", i);
        int sector = TrackBank.getInt("sector", i);

        if (pindex == index_pi_plus)
        {
            in_part_sector->push_back(sector);
            std::cout << sector << endl;
        }
    }

    for (int i = 0; i < nrows_track; i++)
    {
        int pindex = TrackBank.getInt("pindex", i);
        int sector = TrackBank.getInt("sector", i);
        if (pindex == index_pi_minus)
        {
            in_part_sector->push_back(sector);
        }
    }

    std::cout << "after rec" << endl;
}

void read_MC_Part_Bank(hipo::bank PartBank, std::vector<int> required_pids, std::vector<int> *pid_list, std::vector<float> *px_list, std::vector<float> *py_list, std::vector<float> *pz_list, std::vector<int> *in_part_index)
{
   
    std::cout << "in MC" << endl;

    pid_list->push_back(PartBank.getInt("pid", 2));
    in_part_index->push_back(0);
    px_list->push_back(PartBank.getFloat("px", 2));
    py_list->push_back(PartBank.getFloat("py", 2));
    pz_list->push_back(PartBank.getFloat("pz", 2));

    pid_list->push_back(PartBank.getInt("pid", 3));
    in_part_index->push_back(1);
    px_list->push_back(PartBank.getFloat("px", 3));
    py_list->push_back(PartBank.getFloat("py", 3));
    pz_list->push_back(PartBank.getFloat("pz", 3));

    std::cout << "after MC" << endl;

}

/*int main()
{
    return test_4C_hipo();
}*/